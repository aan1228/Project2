# Project2 BASH Script

## Housekeeping

First, the directory structure was constructed.

```bash
mkdir /home/micb405/Group4/project2
mkdir /home/micb405/Group4/project2/Fq
mkdir /home/micb405/Group4/project2/EMIRGE
mkdir /home/micb405/Group4/project2/output
mkdir /home/micb405/Group4/project2/output/Megahit
mkdir /home/micb405/Group4/project2/output/Maxbin
mkdir /home/micb405/Group4/project2/output/MASH
mkdir /home/micb405/Group4/project2/output/PROKKA
```
```Fq``` contains FastQC files of the raw FASTQ data.

```EMIRGE``` contains outputs from Emirge, multiple sequence alignment inputs/outputs, and phylogeny tree files.

```Megahit``` contains outputs from Megahit.

```Maxbin``` contains all the bins generated by Maxbin.

```MASH``` contains taxonomical assignments made by MASH or LAST.

```PROKKA``` contains the annotations of the 20 bins that passed our thresholds (a directory for each bin), a directory called ```RPKM``` which contains outputs used to generate RPKM data, and ```Prokka_All``` which contains the annotations for all of the metagenome. 

![alt text](https://i.imgur.com/qvpdJPC.png)

We then used FastQC to check quality of FASTQ files.

```bash
cd /home/micb405/data/project_2
for i in $(ls *165*); do
fastqc $i -o /home/micb405/Group4/project2/Fq
done;
```

## Metagenome Assembly and Binning

We used Megahit to create contigs, with settings suggested by Connor.
```bash
cd /home/micb405/data/project_2
megahit -1 SI072_LV_165m_DNA_R1.fastq.gz -2 SI072_LV_165m_DNA_R2.fastq.gz --k-min 27 --k-max 147 --k-step 20 --min-contig-len 1000 -m 0.09 -t 2 --out-dir /home/micb405/Group4/project2/output/Megahit
```

Then, we used Maxbin to bin the assemblies by taxonomy.
```bash
cd /home/micb405/Group4/project2/output
perl5.26.0 /home/micb405/resources/project_2/MaxBin-2.2.4/run_MaxBin.pl -contig Megahit/final.contigs.fa -thread 2 -out Maxbin/SI072_LV_165m_DNA_binned -reads /home/micb405/data/project_2/SI072_LV_165m_DNA_R1.fastq.gz -reads2 /home/micb405/data/project_2/SI072_LV_165m_DNA_R2.fastq.gz
```

checkM was run by Connor to determine completeness and contamination of assemblies. We filtered only assemblies greater than 10% completion and less than 5% contamintation from checkM output. Using a script provided by Connor, the list of these bins passing threshold was generated, called ```GT10Complete_LT5Contam_MAGs_checkM.tsv```
```bash
cd /home/micb405/Group4/project2/output/MASH
awk -F"\t" '{ if ($12>10 && $13<5) print $0 }' Group4_checkM_stdout.tsv >GT10Complete_LT5Contam_MAGs_checkM.tsv
```

## Taxonomy Assignment

We ran MASH on the filtered checkM output, based on either RefSeq or Saanich databases. Script was provided by Connor.
For RefSeq:
```bash
cd /home/micb405/Group4/project2/output/MASH
while read line
do
bin=$( echo $line | awk '{ print $1 }')
sid=$( echo $bin | awk -F. '{ print $1 }')
if [ -f ../Maxbin/$bin.fasta ]
    then
    mash dist -v 1E-8 /home/micb405/resources/project_2/refseq.genomes.k21s1000.msh ../Maxbin/$bin.fasta
fi
done< GT10Complete_LT5Contam_MAGs_checkM.tsv > RefSeq_Mash_output.tsv
```
For Saanich:
```bash
while read line 
do  
bin=$( echo $line | awk '{ print $1 }')
sid=$( echo $bin | awk -F. '{ print $1 }')
if [ -f ../Maxbin/$bin.fasta ]
     then
     mash dist -v 1E-8 /home/micb405/resources/project_2/Saanich_QCd_SAGs_k21s1000.sig.msh ../Maxbin/$bin.fasta
fi
done< GT10Complete_LT5Contam_MAGs_checkM.tsv > Saanich_Mash_output.tsv
```
The Saanich and RefSeq MASH outputs were then merged.
```bash
cat RefSeq_Mash_output.tsv Saanich_Mash_output.tsv | sort -t$'\t' -k2,2 | awk '{ if(!x[$2]++) {print $0; dist=($3-1)} else { if($3<dist) print $0} }' > Mash_classifications.BEST.tsv 
```
This was then manually annotated based on GenBank IDs, and called ```Mash_classifications.BEST.annotated.tsv```

LAST was also run using the SILVA databse on the filtered checkM output to complete our taxonomic assignments.
```bash
cd /home/micb405/Group4/project2/output
while read line; do bin=$( echo $line | awk '{ print $1 }'); sid=$( echo $bin | awk -F. '{ print $1 }'); if [ -f Maxbin/$bin.fasta ]; then best_hit=$(lastal -f TAB -P 4 /home/micb405/resources/project_2/db_SILVA_128_SSURef_tax_silva Maxbin/$bin.fasta | grep -v "^#" | head -1); echo $bin,$sid,$best_hit | sed 's/,\| /\t/g'; fi; done<MASH/GT10Complete_LT5Contam_MAGs_checkM.tsv > /MASH/LAST_SILVA_alignments.BEST.tsv
```

## Nitrogen Cycle Gene Annotation

PROKKA was used to annotate genes found in the bins passing our threshold.
```bash
cd /home/micb405/Group4/project2/output
while read line
do
bin=$( echo $line | awk '{ print $1 }')
sid=$( echo $bin | awk -F. '{ print $1 }')
if [ -f /Maxbin/$bin.fasta ]
    then
    prokka /Maxbin/$bin.fasta --outdir PROKKA/$bin --prefix $bin
fi
done< MASH/GT10Complete_LT5Contam_MAGs_checkM.tsv
```
This was also done on the ```final.contigs.fa``` file generated by Megahit to determine if there were nitrogen genes not in our bins passing threshold. Output is in the directory ```Prokka_All```.



## RPKM Abundance of Nitrogen Cycle Genes

The reference set was created for alignment. Script provided by Connor. Requires a ```nitrogen_cyclers.txt``` list of genes of interest.

```bash
cd /home/micb405/Group4/project2/output/PROKKA
while read line; do grep $line */*tsv >>bin_nitrogen_cycler_genes.txt; done<nitrogen_cyclers.txt
while read line
do ffn=$( echo $line | awk -F':' '{ print $1 }' | sed 's/.tsv/.ffn/g' )
prefix=$( echo $line | awk '{ print $1 }' | awk -F':' '{ print $2 }' )
grep "$prefix" $ffn; done<bin_nitrogen_cycler_genes.txt >bin_nitrogen_cycler_headers.txt
cat */*ffn >tmp_All_bins.ffn

/home/micb405/resources/project_2/FastaSubsetter.py -l bin_nitrogen_cycler_headers.txt \
-i tmp_All_bins.ffn -o bin_nitrogen_cycler_genes.ffn -m 1 -v
```

Then, use ```bwa``` to align.

```bash
bwa index bin_nitrogen_cycler_genes.ffn
nohup bwa mem -t 12 bin_nitrogen_cycler_genes.ffn \
/home/micb405/data/project_2/SI072_LV_165m_DNA_R1.fastq.gz \
/home/micb405/data/project_2/SI072_LV_165m_DNA_R2.fastq.gz \
1>bin_nitrogen_cycler_genes_165m.sam 2>bin_nitrogen_cycler_genes.bwa.stderr &
/home/micb405/resources/project_2/rpkm -c bin_nitrogen_cycler_genes.ffn \
-a bin_nitrogen_cycler_genes_165m.sam \
-o bin_nitrogen_cycler_genes_165m_RPKM.csv \
--verbose
```

Make a "matching file" to put taxonomic information alongside PROKKA genes.
```bash
cd /home/micb405/Group4/project2/output/PROKKA
paste <(cat bin_nitrogen_cycler_genes.txt | cut -f 1 | cut -d'/' -f9 | cut -d':' -f2) <(cat bin_nitrogen_cycler_genes.txt | cut -f 1 | cut -d'/' -f8) <(cat bin_nitrogen_cycler_genes.txt | cut -f 3) <(cat bin_nitrogen_cycler_genes.txt | cut -f 1 | cut -d'/' -f8 | cut -d'.' -f2) | sort -k1 > nitrogen_mapping_file.txt
paste <(cat ../MASH/Mash_classifications.BEST.annotated.tsv | cut -f -1 | cut -d';' -f3,4 | sed 's/\s//g') <(cat ../MASH/Mash_classifications.BEST.annotated.tsv | cut -f -2 | cut -d'/' -f8 | sed 's/.fasta//g' | sort -k2) > taxonomy_mapper.tsv 
join -j 2 taxonomy_mapper.tsv nitrogen_mapping_file.txt 
awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,h[$1]}' bin_nitrogen_cycler_genes_165m_RPKM.tsv nitrogen_mapping_file.txt | sort | uniq | tr -s " " | tr ' ' '\t' >> gene_abund_to_bin_rpkm.tsv
cat bin_nitrogen_cycler_genes_165m_RPKM.csv | sed 's/,/\t/g' > bin_nitrogen_cycler_genes_165m_RPKM.tsv
echo -e 'gene\tbin_ID\tprokka_gene_id\tbin_num\trkpm_abundance' > gene_abund_to_bin_rpkm.tsv
awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,h[$1]}' bin_nitrogen_cycler_genes_165m_RPKM.tsv nitrogen_mapping_file.txt | sort | uniq | tr -s " " | tr ' ' '\t' >> gene_abund_to_bin_rpkm.tsv
```

## Constructing Phylogeny
EMIRGE output was provided by Connor. MUSCLE and FastTree were used to perform multiple sequence alignment and tree construction, respectively. The ```165M_fasttree.tree``` file was uploated to [iTOL](http://itol.embl.de/upload.cgi) for visualization.
```bash
cd /home/micb405/Group4/project2/EMIRGE
muscle -in SI072_LV_165m_DNA.emirge_taxa.fasta -out MUSCLE_165m.mfa
FastTree -nt MUSCLE_165M_new.mfa 1>165m_fasttree.tree
```
