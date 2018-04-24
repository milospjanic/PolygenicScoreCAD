# PolygenicScoreCAD

PolygenicScoreCAD for CAD (coronary artery disease) is a combined bash/awk/R script for testing causality of a gene for a given trait, in this case CAD, given the directionality of change of expression level with the increasing number of global risk GWAS allels. Script uses CAD GWAS data from the most recent meta analysis (Nelson et al. Association analyses based on false discovery rate implicate new loci for coronary artery disease. Nat Genet 2017. doi: 10.1038/ng.3913) to select risk SNPs and define risk alleles and HCASMC eQTL data to perform regression analysis.

User provides as arguments:

1. gene of interest
2. chromosome where the gene is located
3. defines p-value threshold to select genome-wide SNPs from Nelson et al. meta analysis for CAD.

# Usage

To run the script download the .sh file
<pre>
wget https://github.com/milospjanic/PolygenicScoreCAD/blob/master/PolygenicScoreCAD.sh
chmod 755 PolygenicScoreCAD.sh
</pre>

Place the script in your home or any other folder. The script will create ~/CADcausalitytest as its working folder, and three subfolders: CARDIOGRAMC4D, HCASMC_expr and HCASMC_genotypes. HCASMC_expr will contain per-gene RNAseq read counts for each HCASMC sample, while HCASMC_genotypes/vcf/ contains whole genome sequencing vcf files of HCASMC samples. CARDIOGRAMC4D folder will contain summary data from Nelson et al. Script will check if all three folders and data sets are present and if not download.

In the folder CARDIOGRAMC4D, file **chrN.GENE.region.txt** will contain all the SNPs from Nelson et al. in the selected region from the tested gene. File **SNP_effect.alele_pval.threshold.txt** will contain SNPs selected with a p-value threshold:

<pre>
rs7483886	G	2.13e-06
rs11532052	T	4.91e-06
rs11226029	G	3.52e-10
rs11226031	C	2.02e-06
...
</pre>

In the HCASMC_expr folder, file **TABLE.RPM.txt** contains per-gene RNAseq read counts for each HCASMC sample. In the HCASMC_genotypes folder, file GENOTYPES.combined.even.HEADER contains risk SNP matrix with counts representing 0, 1 or 2 risk SNPs:

<pre>

1020301 102901 1042702 1051601 1060602 10705 112201 1278 1346 1347 1369 1386 1448 1483 1497 1522 1559 1576 1587 1596 177089 1795 1923 200212 2030801 2040401 20805 2102 2109 2115 2135 2139 2161 2228 2282 2305 2356 24156 2435 2463 2477 2510 289727 2989 3003 3100203 3101801 313605 59386145 59885590 7103002 8072501 8100901 9052004 9070202 9071501 9090701 CA1401
         0 1 1 0 2 1 0 0 0 2 1 1 0 0 0 0 1 0 1 0 2 0 1 1 1 1 0 1 0 1 1 0 2 2 1 0 1 1 2 1 2 0 1 1 0 1 2 1 2 0 1 0 0 1 0 1 1 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
         0 1 1 0 2 1 0 0 0 2 1 1 0 0 0 0 1 0 1 0 2 0 1 1 1 1 0 1 0 1 1 0 2 2 1 0 1 1 2 1 2 0 1 1 0 1 2 1 2 0 1 0 0 1 0 1 1 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
</pre>

Running the script. Place gene name in a file gene.txt, and provide gene.txt, chromosome, distance and p-value threshold as arguments:

<pre>
./GeneCausalityTestCAD.sh gene.txt 11 100000 0.00001
</pre>

# Dependences

UniqueHaplotypeTestCAD requires biomaRt installed in R. In case the biomaRt repository is down run this command to substitute the repository URL with the archive version of biomaRt:

<pre>
sed -i 's/grch37.ensembl.org/jul2016.archive.ensembl.org/g' UniqueHaplotypeTestCAD.sh
</pre>

# Examples

If gene is not causal for the trait (CAD) it will not show correlation with the increased number of risk SNPs:

![alt text](https://github.com/milospjanic/GeneCausalityTestCAD/blob/master/test24.png)

If a gene is causal for the trait and it's increased expression is positively corelated with the trait there should be a positive correlation on the graph:

![alt text](https://github.com/milospjanic/GeneCausalityTestCAD/blob/master/test28.png)

On the other hand, gene might be causal and showing negative correlation with the trait:

![alt text](https://github.com/milospjanic/GeneCausalityTestCAD/blob/master/test29.png)


