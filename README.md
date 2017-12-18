Gclust
===========
Gclust (Genome sequence clustering), a program for clustering the rapid growth of complete or draft genome sequences. Using a sparse suffix array algorithm and a genomic-distance identity taking into account the criteria of diversity, which is based on extension DNA maximal exact matches (MEM), Gclust creates clusters under the given set of genome sequences and extension MEM identity. It takes less than 7 hours for the clustering of the 1560 complete microbial genome sequences with average 3.4MB length on Intel(R) Xeon(R) CPU 2.27GHz with 8 threads parallel computing. It offers the possibility of clustering the rapid growth of complete or draft microbial genomes in the future. 

Gclust is special designed for genome sized sequences clustering and introduced one kind of genomic-distance identity taking into account the criteria of diversity. The fast sparse suffix array construction algorithm was used in finding MEMs between query genome sequence and representative genome sequences. The dynamic programming extension of MEMs is also supported for genome sequence identity computing. Our implementation supports multithreads parallel computing. 

Gclust was written in C++ and uses the SeqAn library and the libdivsufsort library. It is currently maintained by Dr. Beifang Niu (bniu@sccas.cn).

Usage
-----

        Version 1.0
        Usage: ./gclust [options] <genomes-file>

Options:

       -minlen   <int>       set the minimum length of a exact match, if not set, default=20
       -both     <no-args>   compute forward and reverse complement matches, the default value is forward
       -nuc      <no-args>   match only the characters a, c, g, or t
       -sparse   <int>       step of sparse suffix array, default=1
       -threads  <int>       number of threads to use, default=1
       -chunk    <int>       chunk size for one time clustering, the unit is MB, default=100
       -nchunk   <int>       chunk number loaded one time for remaining genomes alignment, default=2
       -loadall  <int>       loading total genomes one time
       -rebuild  <int>       rebuild suffix array after clustering into one chunk, default=1

Clustering cutoff:

       -memiden  <int>       MEMs identity for clustering (default=90, 90% MEMs identity)

MEM extension options:

       -ext       <int>      Extension options, 0: No extension, 1: Gapped extension, 2: Ungapped extension, the default is 1
       -mas       <int>      Reward for a nucleotide match, the default value is 1
       -umas      <int>      Penalty for a nucleotide mismatch, the default value is -1
       -gapo      <int>      Cost to open a gap, the default value is -1
       -gape      <int>      Cost to extend a gap, the default value is -1
       -drops     <int>      X dropoff value for extension, the default value is 1

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
 
	perl script/sortgenome.pl --genomes-file data/viral.1.1.genomic.fna --sortedgenomes-file data/viral.1.1.genomic.sort.fna
 
2. Run gclust:

	./gclust -minlen 20 -both -nuc -threads 8 -ext 1 -sparse 32 data/viral.1.1.genomic.sort.fna > data/viral.1.1.genomic.sort.fna.clustering.out

Description:
Find all MEMs on forward and reverse strands of length 20 or greater, matching only a, c, t, or g, with step of sparse suffix array 32 using 8 threads parallel computing gapped extension for MEMs seeds.

Output
-------
The output file will be written to ./data directory 
     viral.1.1.genomic.fna.clustering.out

Contact
-------
Please contact Dr. Beifang Niu by bniu@sccas.cn if you have any questions.

