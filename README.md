# Deroceras laeve genome
Scripts used in the paper in Genes, Genomes and Genetics by Miranda-Rodriguez et al. 2025

- [Lastz alignments](#lastz-alignments-of-deroceras-assembly)
## Assembly.  Figure 1 and supplementary Figure 1
The PacBio sequence was assembled with Verkko. The final contig assembly is available at zenodo 
10.5281/zenodo.15594897

This file is doubly compressed to save space, first unzip to get a 2bit file:
`gunzip derlae1_hic.2bit.gz`

Then you can get the UCSC executables https://hpc.nih.gov/apps/Genome_Browser.html to get a fasta file and a file with the contig sizes:

`twoBitInfo derlae1_hic.2bit derLae1_hic.chrom.sizes`

Should give something like:
```
ptg000001l      2741356
ptg000002l      1482693
ptg000003l      2021067
ptg000004l      349208
ptg000005l      15547907
```
And `twoBitToFa derlae1_hic.2bit derLae1_hic.fasta` will get you the fasta file used for HiC assembly.

### HiC genome assembly
We followed the pipeline outlined by the Aiden Lab in https://www.dnazoo.org/methods
First, the HiC data, available from the SRA at must be aligned to the fasta contig file with Juicer https://github.com/aidenlab/juicer

However, you must first adapt the section of the auxiliary script `generate_site_positions.py` to your specific experiments. In our case, the restriction enzymes used for the HiC protocol were:
DpnII — GATC
DdeI — CTNAG
MseI — TTAA
HinfI — GANTC

```
  patterns = {
    'HindIII'     : 'AAGCTT',
    'DpnII'       : 'GATC',
    'MboI'        : 'GATC',
    'Sau3AI'      : 'GATC',
    'PhaseGen'    : ['GATC', 'CTNAG', 'TTAA', 'GANTC' ]
  }
```

When ready, just run the preparatory script in your favorite environment. We use Sun Grid Engine to manage our HPC, so you might have to adapt the script to the specific resources. From now one, we just input the general command lines and modules used.

`qsub dpnindex.sge`

Next, customize the juicer.sh script if necessary for your specific protocol. To get the cut sites, you have to think about the HiC protocol. Check the cut sites for every restriction enzyme, and see how it would cut and see how both ends would be 3' end filled. For example, HinfI cuts at G^ANTC. If you take both strands and fill them to the 3' you will see that you have an end 1 "5'...GANT3'" and end 2 "5'GANT3'". Then, if you take all the combinations (substituting N by A,T,G,C) and taking into account all four restriciton enzymes, you get the possible sites:

```
## Set ligation junction based on restriction enzyme
if [ -z "$ligation" ]
then
    case $site in
        HindIII) ligation="AAGCTAGCTT";;
        DpnII) ligation="GATCGATC";;
        MboI) ligation="GATCGATC";;
        NcoI) ligation="CCATGCATGG";;
        none) ligation="XXXX";;
        PhaseGen) ligation="(CTAAGATC|CTTAGATC|CTCAGATC|CTGAGATC|GATCGATC|TTAGATC|GAATGATC|GATTGATC|GACTGATC|GAGTGATC|CTAATAAG|CTTATAAG|CTCATAAG|CTGATAAG|GATCTAAG|TTATAAG|GAATTAAG|GATTTAAG|GACTTAAG|GAGTTAAG|CTAATTAG|CTTATTAG|CTCATTAG|CTGATTAG|GATCTTAG|TTATTAG|GAATTTAG|GATTTTAG|GACTTTAG|GAGTTTAG|CTAATCAG|CTTATCAG|CTCATCAG|CTGATCAG|GATCTCAG|TTATCAG|GAATTCAG|GATTTCAG|GACTTCAG|GAGTTCAG|CTAATGAG|CTTATGAG|CTCATGAG|CTGATGAG|GATCTGAG|TTATGAG|GAATTGAG|GATTTGAG|GACTTGAG|GAGTTGAG|CTAATAA|CTTATAA|CTCATAA|CTGATAA|GATCTAA|TTATAA|GAATTAA|GATTTAA|GACTTAA|GAGTTAA|CTAAAATC|CTTAAATC|CTCAAATC|CTGAAATC|GATCAATC|TTAAATC|GAATAATC|GATTAATC|GACTAATC|GAGTAATC|CTAAATTC|CTTAATTC|CTCAATTC|CTGAATTC|GATCATTC|TTAATTC|GAATATTC|GATTATTC|GACTATTC|GAGTATTC|CTAAACTC|CTTAACTC|CTCAACTC|CTGAACTC|GATCACTC|TTAACTC|GAATACTC|GATTACTC|GACTACTC|GAGTACTC|CTAAAGTC|CTTAAGTC|CTCAAGTC|CTGAAGTC|GATCAGTC|TTAAGTC|GAATAGTC|GATTAGTC|GACTAGTC|GAGTAGTC)";;
        *)  ligation="XXXX"
            echo "$site not listed as recognized enzyme. Using $site_file as site file"
            echo "Ligation junction is undefined"
    esac
fi

```

Use this modified juicer.sh script as above and run with. The genome assembly cookbook mentions to finish at the early stage to get only the `merged_nodups.txt` file, which is the only thing we need to start the chromosome assembly:

```
module load bwa/0.7.17 samtools/1.20
juserdir=/path/juicer-1.6/CPU/
export PATH=$PATH:$juserdir
genoma=/pathtogenome/0_assembly/derLae1_hic.fasta
gsizes=/pathtogenome/0_assembly/derLae1_hic.chrom.sizes
restriction_sites=/results_of_generate_site_script/derLae1_hic_PhaseGen.txt
juicer.sh -g derLae1_hic -s PhaseGen -z $genoma -p $gsizes -y $restriction_sites -D $juserdir -S early -t 60
```

### 3d-dna chromosome assembly
The elegant scripts by (Dudchenko et al., Science, 2017) only need awk installed and java to run juicebox for the manipulation of *.hic files. The use of `parallel` also speeds things up.
The 3D-DNA pipeline was used to error-correct, anchor, order and orient the pieces in the draft assembly. The 3D-DNA visualization module, in conjunction with Juicer Tools, was used to create contact maps for the draft and the final genome assemblies. 
Already the "round 0" of 3D-DNA results in an assembly with 31 chromosomes. The resolution parameters just give a balance between how much repeat or exogenous sequence is included in the scaffolds. This is is one reason why we took care to figure out what is in the unplaced scaffolds after the assembly. And why we aligned the chromosomes to closely related genomes for quality control.
```
module load parallel/20180122
java -version
awk -V | head -n 1
sort --version | head -n 1
export PATH=$PATH:/path_to_3ddna/3d-dna
draft=/pathtogenome/0_assembly/derLae1_hic.fasta
run-asm-pipeline.sh -i 13000 -r 2 --sort-output --fast-start --editor-coarse-resolution 50000 --editor-coarse-stringency 15 --polisher-coarse-stringency 15 --polisher-coarse-resolution 100000 --splitter-coarse-resolution 500000 --build-gapped-map --editor-repeat-coverage 12 $draft /path_to_juicer_output/merged_nodups.txt
```

After reviewing the resulting *.hic file, export the new assembly file, that can be used to get the final maps and sequences.
```
draft=/path_to_genome/0_assembly/derLae1_hic.fasta
export PATH=$PATH:/path_to_3ddna/3d-dna
run-asm-pipeline-post-review.sh --sort-output --build-gapped-map -r derLae1_hic.rawchrom.review.assembly $draft /path_to_juicer_output/merged_nodups.txt
```

Finally, to get the data on Hi-C coverage for Figure 1A and B, we use juicebox_tools to get a wig file with the coverage values at 1kb resolution.

```
bash /path_to_3ddna/3d-dna/visualize/juicebox_tools.sh dump norm SCALE $hic_file assembly BP 1000 "coverage_1kb.wig"
```

### LastZ alignments of Deroceras assembly

The genome of Achatina fulica is available at [this link](https://gigadb.org/dataset/100647).
For this Figure S1, we need [the fasta assembly](https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/100001_101000/100647/Achatina_rebuild.fasta) and the [repeat element annotation gff file](https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/100001_101000/100647/repeat.gff)
Sometimes, the chromosomes names are changed after the article submission, and you have to change them in either the GFF file, or the fasta file, whichever is more convenient with a custom script or bash command.

We use a HARD masked version of the Achatina fulica fasta file, that can be obtained by processing with bedtools or UCSC's `maskOutFa Achatina_rebuild.fasta processed_repeat_gff.bed -clip -soft Achatina_rebuild.fasta.soft.mask` and then `maskOutFa Achatina_rebuild.fastasoft.mask hard Achatina_rebuild.fasta.HARD.mask`.

Since the repeats were masked in the Achatina genome, we proceeded to do repeatmasking of the Deroceras laeve assembly. This was done in three steps, to get a hard masked file, but it was later analyzed and expanded and is decribed in a different section of the paper, so if you are interested in the details, go to the [repeats section](#repeat-identification-with-repeatmasker-and-tandem-repeat-finder).

In short, after RepeatModeler and RepeatMasker, we run a tandem repeat finder identification round and use this masked genome for the whole genome alignments. We again use the strategy of splitting the genome by chromosome or scaffold to parallelize and speed things up, so in the chromosomes/ directory are fasta files like HiC_scaffold_1 ... HiC_scaffold_3497. Another time saving strategy is to keep in `Scaffolds_NOT100p_repeats.txt` the names of the HiC scaffolds that are not completely covered by interspersed or tandem repeats. This list is very easy to obtain with UCSCs `faSize -veryDetailed chromosomes/*` and simple filters.

```
module load lastz/1.04.15 parallel/20180122
fulica=/path_to_afulica_genome/Achatina_rebuild.fasta.HARD.mask
mkdir toFulica
parallel --colsep=" " --will-cite --jobs 80% -a Scaffolds_NOT100p_repeats.txt "lastz ${fulica}[multiple] chromosomes/{1}.fa --output=toFulica/{1}.axt --step=20 --gfextend --chain --gapped --format=axt --notransition --allocate:traceback=1.99G --rdotplot=toFulica/{1}.csv"
```

The axt format is cool to check out the alignment sequences, although they occupy a lot of space. What is used to plot the alignments are the individual csv files. There is a csv file for each HiC scaffold listed in  `Scaffolds_NOT100p_repeats.txt` with the form

`head toFulica/HIC_SCAFFOLD_1.csv`
```
seq1    HIC_SCAFFOLD_1
2200    11203684
2413    11203897
NA      NA
2414    11203905
2419    11203910
NA      NA
```
Where every two rows represent an alignment segment, the numbers in the first column represent the start and end coordinates of the alignment in the A fulica assembly and the second column are the start and end of each segment in each Deroceras laeve scaffold. These csv files, and two tables with the cromosome and scaffold sizes for D laeve and A fulica are all that is needed to reproduce the dotplot of whole genome alignment with this [R script](https://github.com/jerolon/Dlaeve_genome/blob/main/Alignments/plotAlignments.R).

### LastZ D laeve whole genome self-alignment

Generally the same strategy. Use the hardmasked genome from before and split the query by scaffolds to be able to parallelize and speed things up. We use the same modules
`
module load gcc/5.1.0 bedtools/2.27.1 lastz/1.04.15 parallel/20180122
`

We further use a feature of LastZ to mask all sequences that appear in an alignment more than N times. This is called [self masking](https://www.bx.psu.edu/~rsharris/lastz/README.lastz-1.04.15.html#adv_selfmasking). Using N=3. Also `--allocate:traceback=1.99G` or it runs out of memory because there are too many possible seeds in a self-alignment.
```
parallel --colsep=" " --will-cite --jobs 80% -a Scaffolds_NOT100p_repeats.txt "lastz ${assembly}[multiple]   chromosomes/{1}.fa --format=none --gfextend --nochain --gapped --masking=3 --progress+masking=10k --outputmasking+:soft=lastz_identified_repeats/lastz_repeats_{1}.dat --step=20 --allocate:traceback=1.99G"
#Chromosomes 1,9,14,15 were masked with notransition nogaped because they took more than 25 hours
```

With the info contained in the *.dat files, we soft mask genome sequence (which was already hard masked with repeat masker). This is only to avoid too many seeds and saving memory and processing time, these soft masked sequences can still align as mentioned in the lastz documentation.

```
###The sed line turns the .dat file into a bed like file to soft mask the chrmosomes wiith bedtools #Save the soft masked parts to to the directory chromosomes_lastzmask
echo "finished lastz soft mask. Masking with bedtools"
parallel --colsep=" " --will-cite --jobs 80% -a Scaffolds_NOT100p_repeats.txt "bedtools maskfasta -soft -fi chromosomes/{1}.fa -fo chromosomes_lastzmask/{1}.fa -bed <(sed 's/ /\t/g' lastz_identified_repeats/lastz_repeats_{1}.dat)
##use the lastz info to mask the file with the complete genome
echo "finished masking pieces. Masking big"
cat lastz_identified_repeats/lastz_repeats_*.dat | sed 's/ /\t/g' > big_mask.bed
bedtools maskfasta -soft -fi $assembly -fo genoma_lastz_softmasked.fa -bed big_mask.bed
```
And use lastz for self alignment
```
#Use parallel to independently align each scaffold to the whole genome. Save as axt (optional) and as csv (easy to plot in R)
#Here you can play with the parameters --chain --nochain and others, depending on your purposes.
parallel --colsep=" " --will-cite --jobs 80% -a Scaffolds_NOT100p_repeats.txt "lastz genoma_lastz_softmasked.fa[multiple] chromosomes_lastzmask/{1}.fa --format=axt --output=toSelfgap/{1}.axt --rdotplot=toSelfgap/self_{1}.csv --gfextend --nochain --gapped --step=20 --allocate:traceback=1.99G"
```

### Repeat identification with RepeatMasker and Tandem Repeat Finder

First, we used a [repeat library from DFAM](https://www.dfam.org/browse?clade=2697495&clade_ancestors=true&clade_descendants=true) that included the curated repeats of the ancestors and descendants of the spiralia taxon. We installed the DFAM tools in the HPC cluster using singularity. The `-s` option uses a slow a bit more sensitive search.

```
module load singularity/3.7.0
bindings=/path_to_hic_assembly:/genome,/working_directory:/data,/path_to_dfam_library:/Libraries

singularity exec --bind $bindings /cm/shared/apps/singularity/images/3.7.0/dfam-tetools-latest.sif bash -c 'cd /data && RepeatMasker /genome/derLae1_hic.FINAL.fasta -pa 5 -s -gff -no_is -lib /Librerias/Spiralia_Dfam.lib'
```

This step only masks 10% of the genome. The logic is to hard mask these curated previously annotated repeats, and feed the resulting sequence to repeat modeller. In order to save time and avoid identifying de novo repeats that are already known. The genome hard masked with the Spiralian DFMA repeats is then fed to repeat masker, including the LTR modelling pipeline. The next commands were also run with singularity, but this is ommited for clarity.

```
BuildDatabase -name dlaeve_repmod -engine ncbi derLae1_hic.spiraliaMasked.fasta
RepeatModeler -pa 60 -database dlaeve_repmod -LTRStruct
```

This outputs the sequences of de novo identified repeat elements. RepeatModeller runs a simple repeat classifier that reports the repeat family as >rnd-N_family-MMMM#LTR/ERV or similar. However, for non-model organisms this classifier is not good enough and 1110 out of 1705 are classified as "Unknown". These sequences were extracted and classified with [TEclass](https://www.bioinformatics.uni-muenster.de/tools/teclass//index.hbi?). The results were combined with this 
[script](https://github.com/jerolon/Dlaeve_genome/blob/main/Repeats/TEclass_changefa.pl) to add the TEclass prediction with the highest probability to each Unknown familiy. These are identified in the final dlaeve-families.fa library because they contain a string like \[Classified by TEclass with prob=0.476\].

https://zenodo.org/records/15603231

Subsequently, the RepeatModeler library is used to run RepeatMasker on the genome previously masked of spiralia repeats.

```
RepeatMasker spiralia/derLae1_hic.spiraliaMasked.fasta -pa 16 -s -gff -no_is -dir repModeller -lib repModeller/dlaeve-families.fa'
```

The resulting *.cat files from each stage (DFAM-spiralia, RepeatModeller) are concatenated in that order, so that the curated DFAM repeats are annotated with piority as suggested by [Darren Card](https://darencard.net/blog/2017-05-16-maker-genome-annotation/). The concatenated round2Repeatmasker.cat.gz file is then processed by the script ProcessRepeats. This output is the result shown in the paper for the RepeatModeller pipeline. In total, the combination masks 44.34% of total bases, but as stated in the text, a simple look into the masked sequence shows that tandem, simple and low complexity repeats are not masked by this method.

```
ProcessRepeats -a round2RepeatMasker.cat.gz
```

Therefore, we run a third round, on the fasta genome `derLae1_hic.2ndroundRMasker.fasta` masked by RepeatMasker as just described in the previous section, looking only for tandem repeats. The fasta file is split by chromosomes or scaffolds in order to parallelize an run faster.

```
module load UCSC-executables/12may2022/ parallel/20180122 gcc/5.1.0 bedtools/2.27.1

mkdir trftempBorrar
echo "split in SplitChromosomes"
faSplit byname derLae1_hic.2ndroundRMasker.fasta SplitChromosomes/
mkdir masked
parallel --progress --eta "trfBig {} masked/{/.}.masked.fasta -bedAt=masked/{/.}.bed -tempDir=trftempBorrar -trf=trf409.linux64 -l=6" ::: SplitChromosomes/*.fa
cat masked/*.masked.fasta > derLae.3rdRoundTRF.masked.fasta
cat masked/*.bed > derLae.tandemRepeats.bed
bedtools sort -chrThenSizeD -i derLae.tandemRepeats.bed > derLae.TRFsofrted.bed
```

