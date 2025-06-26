# Deroceras laeve genome
Scripts used in the paper in Genes, Genomes and Genetics by Miranda-Rodriguez et al. 2025

- [Lastz alignments](#lastz-alignments-of-deroceras-assembly)
- [Mitochondrial Genome](#mitochondrial-genome)
- [Repeat identification](#repeat-identification)
- [Genome annotation and gene functions](#genome-annotation-and-gene-function)
- [small RNA annotation](#small-rna-annotation)
  
## Assembly.  Figure 1 and supplementary Figure 1
The PacBio sequence was assembled with Verkko. 

To resolve gaps in the assembly, we aligned ONT reads to the Verkko-generated assembly graph using
[GraphAligner](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02157-2). Since
the Verkko output is homopolymer compressed, we first converted the ONT fastq file to fasta then applied
homopolymer compression to the ONT reads using `seqtk hpc` command. The resulting alignment is
recorded in the `graphaligner.gaf` file, we manually inspected to identify the correct ONT reads traversals
across the gap regions. These validated paths were then extracted and provided to Verkko and an input for
gap patching. The following command was used to conduct the alignment:

```
GraphAligner -t 24 -g assembly.homopolymer-compressed.gfa -f ont_hpc.fasta -a graphaliner.gaf --seeds-
mxm-windowsize 5000 --seeds-mxm-length 30 --seeds-mem-count 10000 --bandwidth 15 --multimap-
score-fraction 0.99 --precise-clipping 0.85 --min-alignment-score 5000 --hpc-collapse-reads --discard-
cigar --clip-ambiguous-ends 100 --overlap-incompatible-cutoff 0.15 --max-trace-count 5 --mem-index-no-
wavelet-tree
```

An example of a read spanning a gap in the assembly appears as follows:

```
1880c0f9-b3d2-4e41-a606-b3cf4de00428 133770 59 130460 + <utig4-577<gapont1-1-len-46542-cov-3>utig4-1953 98221777
98091240 98221777 0 0 60
```

In this example, the read `1880c0f9-b3d2-4e41-a606-b3cf4de00428` aligns across two unigs `utig4-577`
and `utig4-1953`, indicating a correct traversal that is used to patch the gap. The read spanning the gap is
46kb (46542) long. You then prepare a tab delimited `paths.gaf` file as follows: `Chr1 <utig4-577<gapont1-1-len-46542-cov-3>utig4-1953 Hap`. Run Verkko again using `--paths path.gaf`.

The resulting assembly will have all the gaps patched.

The final contig assembly is available at zenodo 
10.5281/zenodo.15594897

This file is doubly compressed to save space, first unzip to get a 2bit file:
`gunzip derlae1_hic.2bit.gz`

Then you can get the [UCSC executables](https://hpc.nih.gov/apps/Genome_Browser.html) to get a fasta file and a file with the contig sizes:

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

### Assembly statistics

N50 statistics were calculated using quast version 5.2.0.
Additionally, BUSCO was run in genome mode in two steps, as adviced in the developer documentation:

```
#modules
module load anaconda3/2021.05
source activate busco542

#copy the augustus_config_path into a location with write permissions
export AUGUSTUS_CONFIG_PATH=/custom/augustus_config_path

busco -i /path_to_genome/0_assembly/derLae1_hic.fasta \
      --auto-lineage-euk --augustus --long \
      --augustus_parameters='--progress=true,--AUGUSTUS_CONFIG_PATH=/custom/augustus_config_path/' \
      -l metazoa_odb10 -m genome \
      -o genome_metazoa_BUSCO -c 30

#optimize and train augustus with the results from the first run
module load augustus/3.3.2
export PATH=/pathto/augustus/3.3.2/scripts:$PATH

optimize_augustus.pl --kfold=30 --cpus=30 --species=BUSCO_genome_metazoa_BUSCO --metapars=/custom/augustus_config_path/species/BUSCO_genome_metazoa_BUSCO/BUSCO_genome_metazoa_BUSCO_metapars.cfg ./genome_metazoa_BUSCO/run_metazoa_odb10/augustus_output/training_set.db

etraining --species=BUSCO_genome_metazoa_BUSCO ./genome_metazoa_BUSCO/run_metazoa_odb10/augustus_output/training_set.db
```

Next, run a second round of BUSCO with either the metazoan or mollusca database:

```
export AUGUSTUS_CONFIG_PATH=/custom/augustus_config_path

busco -i /path_to_genome/0_assembly/derLae1_hic.fasta \
      -l [mollusca_odb10,metazoa_odb_10] -m genome \
      --out busco_[mollusca,metazoa]_round2 -c 30 --augustus_species BUSCO_genome_metazoa_BUSCO --augustus_parameters='--progress=true,--AUGUSTUS_CONFIG_PATH=/custom/augustus_config_path/'

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

The axt format is good to check out the alignment sequences, although they occupy a lot of space. What is used to plot the alignments are the individual csv files. There is a csv file for each HiC scaffold listed in  `Scaffolds_NOT100p_repeats.txt` with the form

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
And use lastz for self alignment. The use of the `--inner=2000` option does not seem to have made that much difference. What is important is `--nochain` because chaining would just give you the trivial alignment and not find much outside the diagonal. It is also possible to use the options `--self --notrivial` but we have trouble with lastz accepting a multifasta file in that mode.

```
#Use parallel to independently align each scaffold to the whole genome. Save as axt (optional) and as csv (easy to plot in R)
#Here you can play with the parameters --chain --nochain and others, depending on your purposes.
parallel --colsep=" " --will-cite --jobs 80% -a Scaffolds_NOT100p_repeats.txt "lastz genoma_lastz_softmasked.fa[multiple] chromosomes_lastzmask/{1}.fa --format=axt --output=toSelfgap/{1}.axt --rdotplot=toSelfgap/self_{1}.csv --gfextend --inner=200 --nochain --gapped --step=20 --allocate:traceback=1.99G"
```

the directory toSelfgap has a self_HiC_scaffold_[NN].csv file for each scaffold with the segments aligned to the whole genome. To convert this into the dotplot in Figure 1C, use the R script [plotSelfAlignment.R](Alignments/plotSelfAlignment.R). It is pretty straight forward, you can play with plotting just segments of different sizes or playing with the color and alpha (which in this plot, depend on aligned segment size) in this line:

```
minalsize <- 200
p <- filter(alignments, width >= minalsize) %>% ggplot() + geom_point(aes(x= laeve_ass_0, y = laevew_0, alpha = width, col = log(width)), size = 0.1) + scale_color_gradient2(low = "lightgray", high = "black", mid = "gray", space = "Lab", midpoint = log(10), limit = c(0, log(500))) + theme_bw()
```
Or try a 2D density plot, which we did not include in the manuscript:
```
p <- filter(alignments, width >= minalsize) %>% ggplot() + geom_bin2d(aes(x= laeve_ass_0, y = laevew_0), bins = 3000) + theme_bw() + scale_fill_gradient2(low = "lightgray", high = "red", mid = "blue", space = "Lab", midpoint = 2000)
```

### Synteny comparison to other limacidae whole genome assemblies

The chromosome-sized Hi-C scaffolds of *Deroceras laeve* are numbered according to size, with HiCscaffold1 the longest and HiCscaffold31 the shortest chromosomes. To aid future evolutionary comparisons with the closely related *Deroceras lasithionense*, we use [Ragtag](https://github.com/malonge/RagTag?tab=readme-ov-file) to find homologous chromosomes between the two species. We use `ragtag-2.1.0` and `minimap2/2.24` and exclude unplaced scaffolds in both species:

```
ragtag.py scaffold ../Dlasithionense/GCA_964271515.2_xgDerLasi3.1_genomic.fna derLae1_hic.FINAL.fasta -e dlasi_unanchored.txt  -j dlaeve_unanchored.txt -t 2

```

The resulting `ragtag.scaffold.fasta` contains the *Deroceras laeve* sequences in the ordering and orientation defined by the *Deroceras lasithionense* genome. There is no change in the sequence, but Chr1,2,3..etc, are not ordered by length anymore, they are just the homologs of the respective *Deroceras lasithionense* chromosomes. The `ragtag.scaffold.agp` file is provided as supplementary material and can be used to transform GFF files from the Hi-C coordinates to the Chromosome coordinates as:

```
ragtag.py updategff [-c] <genes.gff> <ragtag.agp>
```

## Mitochondrial genome

To identify the mitochondrial genome, we used a reference that is available from NCBI [NC_072953.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_072953.1/) and used blastn to get matches from our genome. This gave several Kb length hits to HiC_scaffold_1563. Next, we used lastz with default parameters to align NC_072953.1 to scaffold 1563 to get the plot in Figure 2B. From this, the concatamerization was obvious. Next, we loaded NC_072953.1 in SnapGene along with all its annotated features from NCBI directly and used the blast and lastz results as a guide to align Scaffold 1563 along its length. From figure 2C, and observing the sequence, it was obvious that scaffold 1563 covered the mitochondrial genome almost twice, and even though it had some mutations relative to the reference, scaffold 1563 constituted a single sequence. As stated in the text, we used snapgene `Replace Original with Aligned` function to get a single, circularized, annotated mitochondrial sequence for the INB-UNAM strain of *Deroceras laeve*.

### Limacid Mitochondrial phylogeny and dN/dS calculations

Apart from the reference *D. laeve* mitochondial genome, we use the organlle genomes of [*Deroceras reticulatum* NC035495](https://www.ncbi.nlm.nih.gov/nuccore/NC035495) and [*Ambigolimax valentianus* NC072954]. Since the three sequences come from the [NCBI Organelle RefSeq Project](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA927338) they are standardized to start at the first base of the Cox1 gene. Therefore, in SnapGene, where you can choose the start base of a circular chromosome, we also fix the start of the sequence of our strain in this way and export as fasta. Finally, the mitochondrial genome of [*Deroceras lasithionense*](https://www.ncbi.nlm.nih.gov/nuccore/OZ186211), has a similar problem to ours, in that it is a raw Hi-C assembly, however, this genome is 32Kb, compared to the 14Kb of the rest, and moreover, it does not seem to be concatamerized in the coding sequences, the expansion corresponds to the control sequences. In any case, we proceded to find the start of Cox1 and make that the start of the sequence. We also found the trnK(ttt) gene right before Cox1 and make that the end of the sequence. With this alignment, we produce the tree in Figure 2F. Finally, we align each pair of protein coding sequences from both strains of Deroceras laeve (NC_072953.1 & this study) and we calculate the dN/dS rate of non-synonymous to synonymous substitutions. 

Alignment is performed for each protein and output in fasta format using mafft with default parameters, e.g. `mafft ATP6.fasta > ATP6_aligned.fasta`. Then, each alignment is processed with `codeml`, part of the paml suite using this [bash script](Mitochondria/runcodeml.sh). 

## small RNA-seq mapping and annotation

The annotation of micro-RNAs and piwi-interacting RNAs is available at
[https://github.com/pepap/DL-genome-sRNAs/tree/main](https://github.com/pepap/DL-genome-sRNAs/tree/main)

## Repeat identification

### RepeatMasker and Tandem Repeat Finder

First, we used a [repeat library from DFAM](https://www.dfam.org/browse?clade=2697495&clade_ancestors=true&clade_descendants=true) that included the curated repeats of the ancestors and descendants of the spiralia taxon. We installed the DFAM tools in the HPC cluster using singularity. The `-s` option uses a slow a bit more sensitive search.

```
module load singularity/3.7.0
bindings=/path_to_hic_assembly:/genome,/working_directory:/data,/path_to_dfam_library:/Libraries

singularity exec --bind $bindings /cm/shared/apps/singularity/images/3.7.0/dfam-tetools-latest.sif bash -c 'cd /data && RepeatMasker /genome/derLae1_hic.FINAL.fasta -pa 5 -s -gff -no_is -lib /Librerias/Spiralia_Dfam.lib'
```

This step only masks 10% of the genome. The logic is to hard mask these curated previously annotated repeats, and feed the resulting sequence to repeat modeller. In order to save time and avoid identifying de novo repeats that are already known. The genome hard masked with the spiralian DFMA repeats is then fed to repeat masker, including the LTR modelling pipeline. The next commands were also run with singularity, but this is ommited for clarity.

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

### HiTE repeat identification

We used [HiTE: a fast and accurate dynamic boundary adjustment approach for full-length transposable element detection and annotation](https://pmc.ncbi.nlm.nih.gov/articles/PMC11219922/) as an alternative pipeline for interspersed repeat annotation:

```
bindings=/path_to_hic_assembly/:/data

singularity run -B $bindings hite_3.2.0.sif bash -c ' cd /data && python /HiTE/main.py \
--genome /data/derLae1_hic.FINAL.fasta \
--outdir /data/repeats/HiTE \
--thread 80 \
--domain 1 \
--plant 0 \
--annotate 1 \
--intact_anno 1'
```

HiTE automatically runs RepeatMasker with the `--annotate 1` option. According to the logs, the version that came bundled with HiTE 3.2 was RepeatMasker 4.1.1. It uses its own output file `confident_TE.cons.fa` as a custom repeat library. It outputs a HiTE.out and HiTE.gff file showing the location of the annotated Transposable elements. The option `--intact_anno 1` outputs a HiTE_intact.sorted.gff3 of elements that are judged to be intact copies. We use all the confident repeats to get the divergence landscapes shown in Figure 3 A-D, because if we only used the intact, we could be throwing out the most divergent ones. For the genomic repeat density in Figure 3 E-J, and Figure S2, we use the intact transponsable elements because they are fewer in number and possibly more interesting. 

The repeat divergence landscapes can be calculated with RepeatMasker's own auxiliary scripts. Unfortunately, for repeat divergence we need the *.align file and HiTe does not run RepeatMasker with the `-a` option. Therefore, we run RepeatMasker again.

```
bindings=$(pwd):/data,/path_to_hic_assembly:/genoma
echo "Getting the .align file"
singularity exec --bind $bindings /singularity/images/3.7.0/dfam-tetools-latest.sif bash -c 'cd /data && RepeatMasker -a -nolow -no_is -lib confident_TE.cons.fa /genoma/derLae1_hic.FINAL.fasta'

echo "Finished RepMask. Calculating divergence"

singularity exec --bind $bindings /singularity/images/3.7.0/dfam-tetools-latest.sif bash -c 'cd /data && perl /opt/RepeatMasker/util/calcDivergenceFromAlign.pl -s derLae1_hic.FINAL.fasta.divsum derLae1_hic.FINAL.fasta.align'

echo "calculando repeatLandscape"

singularity exec --bind $bindings /cm/shared/apps/singularity/images/3.7.0/dfam-tetools-latest.sif bash -c 'cd /data && perl /opt/RepeatMasker/util/createRepeatLandscape.pl -div derLae1_hic.FINAL.fasta.divsum -twoBit derLa1_hic.FIINAL.2bit > derLaeveHite.html'
```

The  they output an interactive html file that can be checked [here](https://jerolon.github.io/derLaeve_rm.html) for the [Repeats identified in the previous pipeline](#repeatmasker-and-tandem-repeat-finder) and [here] for the repeats identified by HiTE. Comparing them side by side shows the higher sensitivity of the RepeatModeller pipeline, that can detect highly divergent LINEs. For the paper, we use the information contained in the *.divsum file to make our own repeat landscapes with custom colors and different sub-families using the R script [repeatLandscape_hite.R](Repeats/repeatLandscape_hite.R). To know the line where the relevant info starts, use `grep "Coverage for each repeat class and divergence (Kimura)" *.divsum -n`.

### Plot of intact transposable elements: Figure 3 and Figure S2

No we have all the files needed to reproduce the figures 3 and S2 that show the genomic density of different kind of reapeats along the chromosomes and unplaced scaffolds. We also use the Table S6, see [small RNA-seq mapping and annotation](#small-rna-seq-mapping-and-annotation) to plot the location of piRNA clusters. The script to make this plots is [plotRepeatKaryo.R](Repeats/plotRepeatKaryo.R). The input needed, apart from the table of chromosome sizes is:

- `HiTE_intact.sorted.gff3` file with Intact TEs we get from HiTE.
- `derLae.TRFsorted.bed` from Tandem Repeat Finder
- `derLae1_hic.FINAL.fasta` The genome assembly we use to see the GC % of TRFs and which TRFs contain telomeric motifs.

## Genome annotation and gene function

We ran BRAKER3 with the following settings:

```
braker.pl --species=derLae1 --genome=/path_to_hic_assembly/derLae1_hic.FINAL.fasta --rnaseq_sets_dirs=/Transcripts
 --rnaseq_sets_ids=Dc2 Dc3, Dc4, Dc6, Di1, Di4, Di5, Di6, Di7 --prot_seq=Metazoa.fa –gff
```

`Metazoa.fa` is a fasta file downloaded from OrthoDB, and the directory Transcripts contains RNA-seq data from body wall and tail experiments.

### PASA assemblies

As described in the text, RNA-seq data comes from two different labs, one is Single-end and the other is Paired-end, are different in coverage and cover different tissues. What they have in common is that they retain strand information. Therefore, we used TRINITY to assemble separately *de novo* transcriptomes for each batch. We also make a Genome guided transcriptome using the head, juvenile and ovotestis dataset.

We used [PASA](https://github.com/PASApipeline/PASApipeline/wiki) to align the transcripts to the genome and combine them into assemblies.
The transcriptomes are concatenated into a single `transcripts.fasta`. The IDs of *de novo* assemblies are listed in `tdn.accs`. The script to align and build the database is [assembly_pasa.sge](Annotation/assembly_pasa.sge).

This gets us a database of all transcripts, with splicing info and transdecoder rating of their protein-coding potential.

### Evidence modeller

We combine the Braker3 annotations with the pasa Assemblies using Evidence modeller: 
```
module load evidence-modeler/2.1.0

EVidenceModeler \
    --sample_id derLae \
    --weights weightsfile.txt \
    --genome derLae1_fullsoftmask.fasta \
    --gene_predictions gene_predictions.gff \
    --transcript_alignments ../pasa_db.db.pasa_assemblies.gff3 \
    --repeat Full_mask/Dlaeve.full_mask.out.gff \
    --segmentSize 1000000 \
    --overlapSize 30000 \
```

As described in the [pasa documentation](https://github.com/PASApipeline/PASApipeline/wiki/PASA_abinitio_training_sets), we extract the predicted ORFs into a GFF3 file, and we concatenate it with the GFF3 prediction from braker into the file `gene_predictions.gff`. The braker predictions have "AGAT, gmst, or GeneMark.hmm3" in the third column of the GFF3, while the pasa ORFs have "transdecoder". We also pass the GFF3 of pasa assemblies as transcript alignments to Evidence Modeller. Therefore, the weights file looks like this:

```
ABINITIO_PREDICTION     AGAT    1
ABINITIO_PREDICTION     gmst    1
ABINITIO_PREDICTION     GeneMark.hmm3   1
TRANSCRIPT      assembler-pasa_db.db    10
OTHER_PREDICTION        transdecoder    5
```

The resulting `derLae.EVM.pep` file with peptides is benchmarked against busco databases metazoa_odb10 and mollusca_odb10 and is reported in figure 4A as just **EVM**.

### PASA enrichment of EVM annotation.

Next, we used the PASA database generated before, to enrich and correct the EVM annotations, adding UTRs, mRNA variants, fuse split genes, etc.

```
singularity exec --bind $bindings --env TMPDIR=/data/temp/ /singularity/images/3.7.0/pasapipeline_2.5.3.sif bash -c 'cd /data && /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c annotCompare.config -A -g derLae1_fullsoftmask.fasta -t transcripts.fasta.clean -L --annots reference.gff3 --CPU 34'
```

Where `reference.gff3` is just a soft link to the gff3 output of EVM. This results in the final annotation file we use that is reported as EVM post PASA. It is marginally better in terms of Busco scores, but it contains more transcript evidence, transcript isoforms and UTR information.

We also use AGAT's [agat_sp_statistcs](https://agat.readthedocs.io/en/latest/tools/agat_sp_statistics.html) perl script to get the statistics reported in the text.

### Gene functions

Gene annotations were generated using eggNOG-mapper. We used the coding protein FASTA file generated after the annotation and aligned it to the eggNOG5 database. Best hit gene names were extracted from the table and added to the GFF/GTF files.

EggNog also gives the orthologs that the query matches to its database. Since it specializes in assigning fine grained orthologs to queries, we use this information to compare within our own proteome, which genes share orthologs with the script [jaccard_orthologs.R](Annotation/jaccard_orthologs.R).

The construction of the matrix might seem confusing, but refer to the documentation of sparseMatrix. Basically, we construct a vector "rows" and a vector "columns", of the same length, and in each position, rows gives the gene id and columns give the EGGNOG ortholog Id. Then sparseMatrix transforms that into matrix.

A function `jaccard_sparse` computes the jaccard intersection between gene i and gene j, by counting how many orthologs the two genes have in common. 

```
jaccard_sparse <- function(i, j, matrix) {
  intersect_count <- sum(matrix[i, ] & matrix[j, ])
  union_count <- sum(matrix[i, ] | matrix[j, ])
  return(intersect_count / union_count)
}
```

Try all 100,000,000+ combinations of genes, and report only the pairs of genes that substantially share orthologs (i.e. jaccard > 0.5). This script took around two weeks in the cluster. 

### Circos plot of ohnologs, paralogs and clusters. Figure 5.

We make a soft link to the annotation gff3 file called `Dlaeve_annotated.gff3`.
With the script [create_circos_links_orthologs.R](Annotation/create_circos_links_orthologs.R), we import the gff3 to make a TxDb file that we later use for other analyses, like the hox-cluster plot and some intron-exon statistics mentioned in the text. We also import the tsv file with the homologous pairs from *D. laeve* and their jaccard score, to convert to a file that is readable by circos, `links.txt`.

Another file we prepare for the Circos plot, and for Table S7, is the clusters of genes, that in less than 5Mb have 5 or more members that are connected between them by having common orthologs. In the script [ortholog_density_clustering.R](Annotation/ortholog_density_clustering.R) follows these steps:

- Load the ortholog pairs with their jaccard simmilarity. Load the TxDb object with the gene coordinates.
- Get only pairs genes that share the same chromosome.
- Make weights for each connected pair as jaccard_simmilarity/log2(genomic distance)
- Make a graph (network) and put the weights to the connections
- Using those weights, cluster the graph with louvain algorithm
- Filter the clusters with more than 5 members, that span less than 5Mb, and get their probable gene function from the functions.tsv file from Eggnog

Finally, when manually looking for Hox genes, it was obvious that most were in either chromosome 6 or chromosome 25, so we manually prepared a file with their coordinates with grep and awk called `hox_cluster.txt`.

The file [minimal.conf](Annotation/minimal.conf) has the instructions for circos to make the plot in figure 5:

```
module load circos/0.69-6
circos -m minimal.conf
```



