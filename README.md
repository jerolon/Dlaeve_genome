# Dlaeve_genome
Scripts used in the paper in Genes, Genomes and Genetics by Miranda-Rodriguez et al. 2025

##Assembly.  Figure 1 and supplementary Figure 1
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

To get the cut sites, you have to think about the HiC protocol. Check the cut sites for every restriction enzyme, and see how it would cut and see how both ends would be 3' end filled. For example, HinfI cuts at G^ANTC. If you take both strands and fill them to the 3' you will see that you have an end 1 "5'...GANT3'" and end 2 "5'GANT3'". Then, if you take all the combinations (substituting N by A,T,G,C) and taking into account all four restriciton enzymes, you get the possible sites:

(CTAAGATC|CTTAGATC|CTCAGATC|CTGAGATC|GATCGATC|TTAGATC|GAATGATC|GATTGATC|GACTGATC|GAGTGATC|CTAATAAG|CTTATAAG|CTCATAAG|CTGATAAG|GATCTAAG|TTATAAG|GAATTAAG|GATTTAAG|GACTTAAG|GAGTTAAG|CTAATTAG|CTTATTAG|CTCATTAG|CTGATTAG|GATCTTAG|TTATTAG|GAATTTAG|GATTTTAG|GACTTTAG|GAGTTTAG|CTAATCAG|CTTATCAG|CTCATCAG|CTGATCAG|GATCTCAG|TTATCAG|GAATTCAG|GATTTCAG|GACTTCAG|GAGTTCAG|CTAATGAG|CTTATGAG|CTCATGAG|CTGATGAG|GATCTGAG|TTATGAG|GAATTGAG|GATTTGAG|GACTTGAG|GAGTTGAG|CTAATAA|CTTATAA|CTCATAA|CTGATAA|GATCTAA|TTATAA|GAATTAA|GATTTAA|GACTTAA|GAGTTAA|CTAAAATC|CTTAAATC|CTCAAATC|CTGAAATC|GATCAATC|TTAAATC|GAATAATC|GATTAATC|GACTAATC|GAGTAATC|CTAAATTC|CTTAATTC|CTCAATTC|CTGAATTC|GATCATTC|TTAATTC|GAATATTC|GATTATTC|GACTATTC|GAGTATTC|CTAAACTC|CTTAACTC|CTCAACTC|CTGAACTC|GATCACTC|TTAACTC|GAATACTC|GATTACTC|GACTACTC|GAGTACTC|CTAAAGTC|CTTAAGTC|CTCAAGTC|CTGAAGTC|GATCAGTC|TTAAGTC|GAATAGTC|GATTAGTC|GACTAGTC|GAGTAGTC).

When ready, just run the preparatory script in your favorite environment. We use Sun Grid Engine to manage our HPC, so you might have to adapt the script to the specific resources. From now one, we just input the general command lines and modules used.

`qsub dpnindex.sge`








