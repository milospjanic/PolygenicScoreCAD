#!/bin/bash

WORKDIR=~/CADcausalitytest
GENENAME="$(cat $1)"
GENE=$(pwd)/$1
CHR=$2
#DISTANCE=$2
TRESHOLD=$3
NELSON=UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt
CC4D=~/CADcausalitytest/CARDIOGRAMC4D
EXPR=$WORKDIR/HCASMC_expr
GENOTYPES=$WORKDIR/HCASMC_genotypes
VCF=$WORKDIR/HCASMC_genotypes/vcf
REV=$EXPR/reverse

if [ ! -d $WORKDIR ]
then
mkdir $WORKDIR
fi

if [ ! -d $CC4D ]
then
mkdir $CC4D
fi

cd $CC4D

if [ ! -f $NELSON ]
then
wget https://www.dropbox.com/s/hjmkk5rjj0tfswy/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt?dl=0
echo "Unpacking CARDIOGRAM plus C4D data, Nelson et al."
unzip ...
fi

echo "Selecting p-value threshold" $TRESHOLD

cut -f2,3,5,10 UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt > UKBB.GWAS_SNP_effect.alele_pval.txt
awk -v var="$TRESHOLD" '{if ($4<=var) print $0}' UKBB.GWAS_SNP_effect.alele_pval.txt > SNP_effect.alele_pval.threshold.txt

#create genotype files

if [ ! -d $GENOTYPES ]
then
mkdir $GENOTYPES
fi


cd $GENOTYPES


if [ ! -d $VCF ]
then 
mkdir $VCF
cd $VCF
wget https://www.dropbox.com/s/nnytxlbx1v0gh8y/phased_and_imputed.tar
echo "Unpacking genome vcf files..."

tar -xvf phased_and_imputed.tar
gunzip phased_and_imputed*
fi

cd $VCF

grep CHROM phased_and_imputed.chr$CHR.vcf > HEADER.txt
awk '{$1=$2=$3=$4=$5=$6=$7=$8=$9=""; print $0}' HEADER.txt > HEADER.txt.cut
#cat HEADER.txt.cut


while read -r a b c d; do
	grep "$a	" phased_and_imputed.chr$b.vcf > SNP.txt
	#cat SNP.txt
        #clean genotype files
        sed -i -E 's/:[0-9.,:]*//g' SNP.txt
	#sed -i -E 's/:[0-9.]*:[0-9.]*,[0-9.]*,[0-9.]*//g' SNP.txt 
	#cat SNP.txt
       
#grab reference and alternative aleles
        REF="$(awk '{printf $4}' SNP.txt)"
        ALT="$(awk '{printf $5}' SNP.txt)"
	#echo $REF
	#echo $ALT

        EFFECT=$c
        #echo $EFFECT
	
	sed -i "s/0|0/$REF$REF/g" SNP.txt
        sed -i -E "s/0\|[1|2]/$REF$ALT/g" SNP.txt
        sed -i -E "s/[1|2]\|0/$REF$ALT/g" SNP.txt
        sed -i "s/1|1/$ALT$ALT/g" SNP.txt
	#cat SNP.txt
        
        awk '{$1=$2=$3=$4=$5=$6=$7=$8=$9=""; print $0}' SNP.txt > SNP.txt.cut 
        #cat SNP.txt.cut

	grep -v -P "[ATCG,]{3,}" SNP.txt.cut > SNP.txt.cut2
	sed -i "s/"$EFFECT"/1/g" SNP.txt.cut2
        sed -i -E "s/[ATGC]/0/g" SNP.txt.cut2
        
#cat SNP.txt.cut2

        cat HEADER.txt.cut SNP.txt.cut2 > GENOTYPES.$a.txt
	#cat GENOTYPES.$a.txt

done < $CC4D/SNP_effect.alele_pval.threshold.txt

cat GENOTYPES.rs* > GENOTYPES.combined
cat GENOTYPES.*:* >> GENOTYPES.combined

#rm GENOTYPES.rs*

grep -v "1020301 102901" GENOTYPES.combined > GENOTYPES.combined.even
#awk 'NR%2==0' GENOTYPES.combined > GENOTYPES.combined.even

sed -i 's/11/2/g' GENOTYPES.combined.even
sed -i -E 's/(10|01)/1/g' GENOTYPES.combined.even
sed -i 's/00/0/g' GENOTYPES.combined.even

cat HEADER.txt.cut GENOTYPES.combined.even > GENOTYPES.combined.even.HEADER

####awk 'BEGIN{print "count", "lineNum"}{print gsub(/t/,"") "\t" NR}' file
####awk -F'|' -v fld=2 'BEGIN{print "count", "lineNum"}{print gsub(/t/,"",$fld) "\t" NR}' file



if [ ! -d $EXPR ]
then
mkdir $EXPR
fi

cd $EXPR

#check if reverse folder with expression levels counted on reverse strand exists, if not, download from Dropbox link

if [ ! -d $REV ]
then
wget https://www.dropbox.com/s/edm0ykexjmue5yf/reverse.zip
echo "Unpacking expression files..."
unzip reverse.zip
fi

#write R script to get ENSEMBL id, needs biomaRt in R

echo "#!/usr/bin/Rscript
library(biomaRt)
listMarts(host=\"grch37.ensembl.org\")
ensembl = useMart(\"ENSEMBL_MART_ENSEMBL\",dataset=\"hsapiens_gene_ensembl\", host=\"grch37.ensembl.org\")
id_merge = getBM(attributes=c(\"ensembl_gene_id\",\"external_gene_name\"),mart=ensembl)
write.table(id_merge, file=\"id_merge.txt\", sep = \"\t\", quote =F, col.names=F, row.names=F)
" > script.r

#run R script

chmod 775 script.r
./script.r

#Use awk to append gene names

awk 'NR==FNR {h[$2] = $1; h2[$2] = $2; next} {print h[$1]}' id_merge.txt $GENE > genename

#remove temporary files

rm id_merge.txt
rm script.r

#get gene counts for gene of interest

while read line; do
                set $line
                find . -name *gene.count | xargs grep $line > COUNTS.txt
done < genename

sed -E "s/^.\/reverse\///g" COUNTS.txt | sed -E "s/\/.*\.[0-9]*//g" > COUNTS.txt.cut

#get total gene counts per sample

find . -name *gene.count | xargs -I % awk 'BEGIN {FS = " "} ; {sum+=$2} END {print sum}' % > TOTAL.txt
find . -name *gene.count | xargs -I % echo % > SAMPLES.txt 

sed -E "s/^.\/reverse\///g" SAMPLES.txt | sed -E "s/\/.*\.[0-9a-zA-Z]*//g" > SAMPLES.txt.cut

paste SAMPLES.txt.cut TOTAL.txt > TOTALCOUNTS.txt

awk 'NR==FNR {h[$1] = $0; next} {if(h[$1]) print h[$1]"\t"$0}' COUNTS.txt.cut TOTALCOUNTS.txt > TABLE.txt

awk '{print $1 "\t" $2/$4*1000000}' TABLE.txt > TABLE.RPM.txt

cp TABLE.RPM.txt $VCF

cd $VCF

#create ggplot2 graph 

echo "#!/usr/bin/Rscript
library(\"ggplot2\")
data<-read.table (file=\"GENOTYPES.combined.even.HEADER\", head=T, check.names=F)
data<-rbind(data, Total = colSums(data))

data.tr<-t(data)
#data.sum<-rowSums(data.tr)
expr<-read.table(\"TABLE.RPM.txt\", row.names=1, check.names=F)

data.merge<-merge(x = data.tr, y = expr, by = \"row.names\", all = F)

write.table(file=\"data.merge.table.txt\", data.merge)
pdf(\"output.pdf\")
ggplot(data.merge, aes(Total, V2, color = Total)) +
  geom_point(shape = 16, size = 5, show.legend = FALSE, alpha = .4) +
  theme_minimal() +
  scale_color_gradient(low = \"#0091ff\", high = \"#f0650e\") +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=26))+ labs(title = \"$GENENAME Causality Test\", x=\"CAD Risk SNPs Total\", y=\"$GENENAME expression\") + theme(plot.title = element_text(size = rel(2)))+
  geom_smooth(method=lm)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(data.merge\$V2)+max(data.merge\$V2)/5))+
  annotate(x=min(data.merge\$Total)+(max(data.merge\$Total)-min(data.merge\$Total))/5, y=max(data.merge\$V2)+(max(data.merge\$V2)-min(data.merge\$V2))/6,label=paste(\"R = \", round(cor (data.merge\$V2,data.merge\$Total),2)),geom=\"text\", size=8, col=\"darkred\")

dev.off()
"> script.R

chmod 775 script.R
./script.R


