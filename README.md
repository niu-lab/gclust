Gclust
===========
Gclust (A Parallel Clustering Tool for Microbial Genomic Data), a program for clustering the rapid growth of complete or draft genome sequences. Using a sparse suffix array algorithm and a genomic-distance identity taking into account the criteria of diversity, which is based on extension DNA maximal exact matches (MEMs), Gclust creates clusters under the given set of genome sequences and extension MEM identity. It takes less than 10 minutes for the clustering of the 9578 complete microbial genome sequences with average 27Kbp length on 24-core Intel(R) Xeon(R) CPU E5-2680 v3 @2.50GHz with 16 threads parallel computing. It offers the possibility of clustering the rapid growth of complete or draft microbial genomes in the future. 

Gclust is specially designed for genome sized sequences clustering and introduced one kind of genomic-distance identity taking into account the criteria of diversity. The fast sparse suffix array construction algorithm was used in finding MEMs between query genome sequenceand representative genome sequences. The dynamic programming extension of MEMs is also supported for genome sequence identity computing. Our implementation supports multithreads parallel computing. 

Gclust was written in C++ and uses the SeqAn library. It is currently maintained by Dr. Beifang Niu (niubf@cnic.cn).

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

Install
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
