#!/bin/bash

#-------------prepare-----------------------#

ref_pre_miRNA='/data1/linxing/workspace/DMY_sci_data/mirdeep2/quantifier/hairpin.mekada.fa'
ref_mature_miRNA='/data1/linxing/workspace/DMY_sci_data/mirdeep2/quantifier/mature.mekada.fa'
refgenome='/data1/linxing/workspace/DMY_sci_data/mirdeep2/Ola_genome.fa'
index_refgenome='/data1/linxing/workspace/DMY_sci_data/mirdeep2/Ola'
species='ola'

#-----------1. fastq2fasta----------------------------#
ls /data1/linxing/workspace/DMY_sci_data/sRNA_data/*fq |while read id;do sed -n '1~4s/^@/>/p;2~4p' ${id} > `basename ${id} .fq`.fa;done



#-----------2. align to refgenome by bowtie---------#
cd /data1/linxing/workspace/DMY_sci_data/mirdeep2/

nohup bowtie-build Ola_genome.fa Ola &

ls /data1/linxing/workspace/DMY_sci_data/sRNA_data/*.fa |while read id;do mapper.pl ${id} -c -j -l 18 -m -p Ola -s `basename ${id} .fa`.collapsed.fa -t `basename ${id} .fa`_collapsed_vs_genome.arf -o 20; done

:<<EOF
paremeters:
-e: input fq
-i: convert rna to dna alphabet
-k: clip 3' adatper
-l: min seq len
-m: collapse reads
-p: map to genome
-q: allow 1 mismatch
-r: defualt 5,a read is allowed to map up to this number of positions in the genome
-t: mappping outfile name
-s: processed reads outfile name,  collapsed reads
-v: progress report
-o: number of threads to use for bowtie
EOF


#----3. microRNA quantification----#
cd /data1/linxing/workspace/DMY_sci_data/mirdeep2/
mkdir quantifier
cd /data1/linxing/workspace/DMY_sci_data/mirdeep2/quantifier

ls /data1/linxing/workspace/DMY_sci_data/mirdeep2/*collapsed.fa |while read id;do quantifier.pl -p hairpin.mekada.fa -m mature.mekada.fa -r ${id} -y `basename ${id} .clean.collapsed.fa`;done

:<<EOF
paremeters of quantifier.pl
input file：  A fasta file with deep sequencing reads,
                          a fasta file of the corresponding genome,
                          a file of mapped reads to the genome in miRDeep2 arf format,
                          an optional fasta file with known miRNAs of the analysing species，
              and an option fasta file of known miRNAs of related species.

-a int        minimum read stack height that triggers analysis.
              Using this option disables automatic estimation
              of the optimal value.

-t species    species being analyzed - this is used to link to
              the appropriate UCSC browser
-v            remove directory with temporary files
-s file       File with known miRBase star sequences
-P            use this switch if mature_ref_miRNAs contain miRBase
              v18 identifiers (5p and 3p) instead of previous
              ids from v17
EOF
