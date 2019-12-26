Gclust: A Parallel Clustering Tool for Microbial Genomic Data
===========
Gclust is a parallel program for clustering complete or draft genomic sequences, where clustering is ac-celerated with a novel parallelization strategy and a fast sequence comparison algorithm using sparse suffix arrays (SSAs). Moreover, genome identity measures between two sequences are calculated based on their maximal exact matches (MEMs). Gclust is freely available for non-commercial use at https://github.com/niu-lab/gclust. We also introduce a web server for clustering user-uploaded genomes at http://niulab.scgrid.cn/gclust.

Usage
-----

        Version 1.0
        Usage: ./gclust [options] <genomes-file>

Options:

       -minlen   <int>       Set the minimum length of an exact match, if not set, default = 20
       -both     <no-args>   Compute forward and reverse complement matches, default = forward
       -nuc      <no-args>   Match only the characters a, c, g, or t
       -sparse   <int>       Step of sparse suffix array, default = 1
       -threads  <int>       Number of threads to use, default = 1
       -chunk    <int>       Chunk size for one time clustering, the unit is MB, default = 100
       -nchunk   <int>       Chunk number loaded one time for remaining genomes alignment, default = 2
       -loadall  <int>       Loading total genomes one time
       -rebuild  <int>       Rebuild suffix array after clustering into one chunk, default = 1

Clustering cutoff:

       -memiden  <int>       Set extended MEM idendity(eMEMi) or MEM idendity for clustering, default = 90, i.e., the genomes were clustered at 90% eMEMi under the condition of '-ext = 1'

MEM extension options:

       -ext       <int>      Extension options, 0: no extension, 1: gapped extension, 2: un-gapped extension, default = 1
       -mas       <int>      Reward for a nucleotide match, default = 1
       -umas      <int>      Penalty for a nucleotide mismatch, default = -1
       -gapo      <int>      Cost to open a gap, default = -1
       -gape      <int>      Cost to extend a gap, default = -1
       -drops     <int>      X dropoff value for extension, default = 1

Installation
-------

Clone the gclust repos, and build the `gclust` binary:

    git clone https://github.com/niu-lab/gclust
    cd gclust
    make

Now you can put the resulting binary where your `$PATH` can find it. If you have su permissions, then
I recommend dumping it in the system directory for locally compiled packages:
    
    sudo mv gclust /usr/local/bin/

If you have su permissions, then you may just add an environment variable to ~/.bashrc:
    
    export PATH=/your_install_path/gclust:$PATH

Example
-------
1. Sort the input genomes in decreasing order of length:    
```bash    
perl script/sortgenome.pl --genomes-file data/viral.1.1.genomic.fna --sortedgenomes-file data/viral.1.1.genomic.sort.fna
```
2. Run gclust:  
```bash    
./gclust -minlen 20 -both -nuc -threads 8 -ext 1 -sparse 2 data/viral.1.1.genomic.sort.fna > data/viral.1.1.genomic.sort.fna.clustering.out
```
Description:
Find all MEMs on forward and reverse strands of length 20 or greater, matching only a, c, t, or g, with step of sparse suffix array 32 using 8 threads parallel computing gapped extension for MEMs seeds.

Output
-------
The output file will be written to ./data directory 
     viral.1.1.genomic.fna.clustering.out

Contact
-------
Please contact Dr. Beifang Niu by niubf@cnic.cn if you have any questions.
`
