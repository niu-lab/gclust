Gclust: A Parallel Clustering Tool for Microbial Genomic Data
-----
Gclust is a parallel program for clustering complete or draft genomic sequences, where clustering is accelerated with a novel parallelization strategy and a fast sequence comparison algorithm using sparse suffix arrays (SSAs). Moreover, genome identity measures between two sequences are calculated based on their maximal exact matches (MEMs). Gclust is freely available for non-commercial use at https://github.com/niu-lab/gclust. We also introduce a web server for clustering user-uploaded genomes at http://niulab.scgrid.cn/gclust.

Usage
-----

        Version 1.0
        Usage: ./gclust [options] <genomes-file>

Options:

       -minlen   <int>       Set the minimum length for exact match, if not set, default = 20
       -both     <no-args>   Compute forward and reverse complement matches, default = forward
       -nuc      <no-args>   Match only the characters a, c, g, or t
       -sparse   <int>       Set the step of sparse suffix array, default = 1
       -threads  <int>       Set the number of threads to use, default = 1
       -chunk    <int>       Set the chunk size for one time clustering, default = 100, where the unit is million base pairs (Mbp)
       -nchunk   <int>       Set the chunk number loaded one time for remaining genomes alignment, default = 2
       -loadall  <int>       Load the total genomes one time
       -rebuild  <int>       Rebuild suffix array after clustering into one chunk, default = 1

Clustering cutoff:

       -memiden  <int>       Set the value of extended maximal exact match (MEM) idendity or non-extended MEM idendity for clustering, default = 90

Extension options of MEM:

       -ext       <int>      Set the extension type of MEM, where '0' means no extension, '1' means gapped extension and '2' means un-gapped extension, default = 1
       -mas       <int>      Set the reward value for a nucleotide match, default = 1
       -umas      <int>      Set the penalty value for a nucleotide mismatch, default = -1
       -gapo      <int>      Set the cost value to open a gap, default = -1
       -gape      <int>      Set the cost value to extend a gap, default = -1
       -drops     <int>      Set the X dropoff value for extension, default = 1

Installation
-------

Clone the gclust repos, and build the `gclust` binary:

    git clone https://github.com/niu-lab/gclust
    cd gclust
    make

Now you can put the resulting binary where your `$PATH` can find it. If you have root permissions, then
I recommend dumping it in the system directory for locally compiled packages:
    
    sudo mv gclust /usr/local/bin/

If you have root permissions, then you may just add an environment variable to ~/.bashrc:
    
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

Output
-------
If the user does not specify an output file name, the clustering results will be placed in the './data' directory by default, and the default file name is 'viral.1.1.genomic.sort.fna'.

Contact
-------
Please contact Dr. Beifang Niu by niubf@cnic.cn if you have any questions.
