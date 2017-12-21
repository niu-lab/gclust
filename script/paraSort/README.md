paraSort
============
This sub-program can be used when the genome.fna file is  too large, which fulfill the same functional role with 'sortgenome.pl'. The program sort the genomes in decreasing order of length in the distributed environment. Please make sure that the total memory is available for the genome.fna file.

Install
---------
```bash
cd paraSort
make
```
Example 
----------
1. extract index and length 
```bash
perl extract_index_length.pl --genomes-file data/viral.genomic.fna --index-file data/viral.index.length
```
2. sort index
```bash
perl sort_genome_index.pl --index-file data/viral.index.length --sortedindex-file data/viral.index.length.sorted
```
3. sort genomes
```bash
mpirun -np 4 ./para_sortgenome  --input data/viral.genomic.fna  --index data/viral.index.length.sorted --output data/viral.genomic.parasort.fna
```
