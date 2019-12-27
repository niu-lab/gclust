
/****************************
 
 Gclust: For genomes clustering based on all Maximal Exact Matches or MEM extension between genome sequences.
 By Beifang Niu 2011.10.12

****************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <cctype>
#include "fasta.hpp"
#include "paraSA.hpp"

using namespace std;

void usage(string prog);

int K = 1; // Note: Using sparse suffix array for larger chunk size.
int Nchunk = 2; // load one time for remaining genomes clustering.
int chunk = 100; // block size for clustering part by part.
int min_len = 20; // Default minimum exact match length.
int MEMiden = 90; // Default identity cutoff.
int total_threads = 1; // Threads number.

// MEM extension parameters
int ext = 1; // noextension, gap or ungap extension
int mas = 1; // Match score
int umas = -1; // Mismatch cost
int gapo = -1; // Gap open penalty
int gape = -1; // Gap extension penalty
int drops = 1; // Maximum score drop

bool rev_comp = false;
bool nucleotides_only = false;
bool rebuild = false; // Rebuild suffix array into one part.
bool loadall = false; // load all genomes one time, need more memory.

paraSA *sa, *saa; // Suffix array.

vector<Genome> refseqs, allrefseqs; // Part genomes and total part genomes.
vector<GenomeClustInfo> totalgenomes; // Total genomes.
vector<vector<match_t> > matchlist; // Parallel buffer for match_t.
vector<vector<mumi_unit> > mumilist; // parallel buffer for mumi_unit.

struct threads_arg // Multithreads parallel parameters passing.
{ 
	int skip;
	int skip0; // Zebra-style distribution.
	long chunk;
	long begin; // Begin alignment position.
	bool part; // If internal part clustering.

};

// Note: just test distances between genomes.
void testDistanceBgenomes(vector<GenomeClustInfo> &totalgenomes)
{
	for (long i=0; i<(long)totalgenomes.size(); i++)
	{
		if (totalgenomes[i].clusters.size()>0)
		{	
			for (long j=0; j<(long)totalgenomes[i].clusters.size(); j++)
			{
				hit thit;
				thit = totalgenomes[i].clusters[j];
				cerr<<thit.id<<"\t"<<thit.identity<<endl;
			}
		}
	}

}

// Note: output clustering information as cd-hit format.
void outputClusteringInfoSimple(vector<GenomeClustInfo> &totalgenomes)
{
	long clusters = 0;
	for (long i=0; i<(long)totalgenomes.size(); i++)
	{	
		if (totalgenomes[i].rep)
		{
			long clusterunit=0;
			// Output representive.
			cout<<">Cluster "<<clusters<<endl;
			cout<<clusterunit<<"\t"<<totalgenomes[i].size<< \
				"nt, >"<<totalgenomes[i].descript<< \
				"... *"<<endl;
			if (totalgenomes[i].clusterunits.size()>0) 
			{
				for (long j=0;j<(long)totalgenomes[i].clusterunits.size();j++) 
				{
					clusterunit++;
					cout<<clusterunit<<"\t"<< \
						totalgenomes[totalgenomes[i].clusterunits[j].id].size<<"nt, >"<< \
						totalgenomes[totalgenomes[i].clusterunits[j].id].descript \
						<<"... at "<<totalgenomes[i].clusterunits[j].strand<<"/"<< \
						totalgenomes[i].clusterunits[j].identity<<endl;
				}
			}
			clusters++;
		}
	}
	cerr<<"Total clusters: "<<clusters<<endl;

}

// Note: collect clustering information.
void getClusteringInfoOnepart(vector<GenomeClustInfo> &totalgenomes,
															long begin,
															long chunk,
															bool inpart,
															bool &clusterhit)
{	
	long b,e;
	clusterhit=false;
	b=begin;
	e=begin+chunk;
	for (long i=b; i<e; i++)
	{	
		if (!totalgenomes[i].rep){ continue; }	
		if (totalgenomes[i].clusters.size()>0)
		{	
			for (long j=0; j<(long)totalgenomes[i].clusters.size(); j++)
			{
				hit thit;
				thit=totalgenomes[i].clusters[j];
				if (totalgenomes[thit.id].rep)
				{
					hit thit0;
					thit0.id=i;
					thit0.identity=thit.identity;
					thit0.strand=thit.strand;
					totalgenomes[thit.id].clusterunits.push_back(thit0);
					totalgenomes[i].rep=false;
					clusterhit=true;
					break;
				}
			}
			// Note: important clear
			totalgenomes[i].clusters.clear();
		}
	}
}

// Note: one genome as reference (internal part).
void *single_thread(void *arg_)
{
	Genome tg;
	threads_arg *arg = (threads_arg *)arg_;

	// Match information container.
	vector<match_t> &matches=matchlist[arg->skip0];
	// Mem index container.
	vector<mumi_unit> &mumis=mumilist[arg->skip0];

	long seq_cnt = 0;
	long beginclust = arg->begin;
	long chunk = arg->chunk;
	long sizeadd = 0;
	long edge = long(refseqs.size()-1);
	bool ifhit = false;
	bool ispart = arg->part;
	double cutoff=(double)MEMiden/100;
	string *P=new string; 
	edge=long(refseqs.size()-1); 

	while(1) 
	{
		if ( seq_cnt > edge ){ break; }
		if ( arg->skip0 == 0 ){
			if (seq_cnt % 100 ==0){
				cerr<<"...... "<<seq_cnt<<" done"<<endl;
			}
		}
		// paralle part.
		if( seq_cnt % arg->skip == arg->skip0 ) 
		{
			ifhit=false;
			tg=refseqs[seq_cnt];
			if ( totalgenomes[tg.id].rep )
			{
				*P=tg.cont;
				// Filter 'n'.
				if (nucleotides_only)
				{ 
					filter_n(*P);
				}
				// 100% ?
				if (MEMiden==100)
				{
					saa->MEMperfect(*P, matches, tg.size, tg.id);
				}else{
					saa->MEM(*P, matches, min_len, tg.id);
				}
				sizeadd += saa->load_match_info(tg.id, matches, mumis, true, tg.size);
				matches.clear();
				if ((double)sizeadd/tg.size >= cutoff)
				{
					ifhit=ComputeMemIdentity(totalgenomes, 
																	 allrefseqs, 
																	 mumis, 
																	 beginclust, 
																	 tg.id, 
																	 MEMiden, 
																	 ispart, 
																	 chunk, 
																	 '+', 
																	 ext, 
																	 mas,
																	 umas, 
						                       gapo, 
						                       gape,
						                       drops);
				}
				mumis.clear();
        sizeadd=0;

				if ((ispart)||(!ifhit))
				{
					if(rev_comp) {
						reverse_complement(*P, nucleotides_only);
						// 100% ?
						if (MEMiden==100)
						{
							saa->MEMperfect(*P, matches, tg.size, tg.id);
						}else{
							saa->MEM(*P, matches, min_len, tg.id);
						}
						// Loading match information - strand.
					  sizeadd += saa->load_match_info(tg.id, matches, mumis, true, tg.size);
						matches.clear();

						if ((double)sizeadd/tg.size >= cutoff)
						{
							ComputeMemIdentity(totalgenomes, 
								                 allrefseqs, 
								                 mumis,
																 beginclust, 
								                 tg.id, 
								                 MEMiden,
																 ispart, 
								                 chunk, 
								                 '-', 
								                 ext, 
								                 mas,
																 umas, 
								                 gapo, 
								                 gape, 
								                 drops);
						}
						sizeadd=0;
						mumis.clear();
					}
				}
			}
		}
		seq_cnt++;
		delete P; 
		P = new string;

	}
	delete P;
	pthread_exit(NULL);

}

int main(int argc, char* argv[]) 
{
	time_t start, end;
        start=time(NULL);
	// Version notice.
	cerr<<"\nGclust version 1.0\n"<<endl;
	// Collect parameters from the command line.
	while (1) 
	{
		static struct option long_options[] =
		{ 
			{"minlen", 1, 0, 0}, // 0
			{"both", 0, 0, 0}, // 1
			{"nuc", 0, 0, 0}, // 2
			{"threads", 1, 0, 0}, // 3
			{"chunk", 1, 0, 0}, //4
			{"memiden", 1, 0, 0}, //5
			{"nchunk", 1, 0, 0}, //6
			{"loadall", 0, 0, 0}, //7
			{"rebuild", 0, 0, 0}, //8

			// Seed extension part
			{"mas", 1, 0, 0}, //9
			{"umas", 1, 0, 0}, //10
			{"gapo", 1, 0, 0}, //11
			{"gape", 1, 0, 0}, //12
			{"drops", 1, 0, 0}, //13
			{"ext", 1, 0, 0}, //14

			// Sparse step of suffix array
			{"sparse", 1, 0, 0,}, //15

			{0, 0, 0, 0}

		};

		int longindex = -1;
		int c = getopt_long_only(argc, argv, "", long_options, &longindex);
		if(c == -1) break; // Done parsing flags.
		else if(c == '?'){ // If the user entered junk, let him know. 
			cerr << "Invalid parameters." << endl;
			usage(argv[0]);
		}else {
			// Branch on long options.
			switch(longindex) 
			{ 
				case 0: min_len = atol(optarg); break;
				case 1: rev_comp = true;	break;
				case 2: nucleotides_only = true; break;
				case 3: total_threads = atoi(optarg) ; break;
				case 4: chunk = atoi(optarg) ; break;
				case 5: MEMiden = atoi(optarg) ; break;
				case 6: Nchunk = atoi(optarg) ; break;
				case 7: loadall = true ; break;
				case 8: rebuild = true ; break;

				// Seed extension part
				case 9: mas = atoi(optarg) ; break;
				case 10: umas = atoi(optarg) ; break;
				case 11: gapo = atoi(optarg) ; break;
				case 12: gape = atoi(optarg) ; break;
				case 13: drops = atoi(optarg) ; break;
				case 14: ext = atoi(optarg) ; break;
					
				// Sparse step of suffix array
				case 15: K = atoi(optarg) ; break;

				default: break; 
			}
		}
	}

	// Only using all maximal matches for clustering. 
	if (argc - optind != 1) usage(argv[0]);
	if(total_threads <= 0) 
	{ 
		cerr << "invalid number of threads specified" << endl; 
		exit(1); 
	}
	// no extension when 100% match
	if (MEMiden == 100){ ext = 0; }
	// Allocate memory for multithreads.
	for (int i=0;i<total_threads; i++)
	{
		vector<match_t> matches;
		vector<mumi_unit> mumis;
		matchlist.push_back(matches);
		mumilist.push_back(mumis);
		matchlist[i].reserve(MAX_THREADCONTAINER);
		mumilist[i].reserve(MAX_THREADCONTAINER);
	}

	// Genome file.
	string ref_fasta = argv[optind]; 
	// Load total genomes part information.
	load_total_genomes(ref_fasta, totalgenomes);
	// Load total part genomes one time.
	if (loadall) load_part_genomes_all(ref_fasta, allrefseqs); 

	vector<long> refdescr;
	vector<long> startpos;
	long chunksize;
	long dchunk;
	//long genomes; //warning: variable 'genomes' set but not used
	long begin = 0;
	long cbegin = 0;
	bool ifend=false;
	bool clusterhit=false;
	Genome tg; 
	string ref;

	// set chunk size for clustering chunk by chunk.
	chunksize=(long)chunk*PART_BASE;

	// Note: Take a fixed value of total genomes number.
	// Main clustering loop.
	while (1) 
	{
		if (loadall){ 
			load_part_genomes_internal_mem(allrefseqs, 
				                             refseqs, 
																		 totalgenomes, 
				                             cbegin, 
																		 chunk, 
				                             chunksize, 
				                             ifend,
																		 MEMiden);
		}else{
			load_part_genomes_internal(ref_fasta, 
				                         refseqs, 
																 totalgenomes, 
				                         cbegin, 
																 chunk, 
				                         chunksize, 
				                         ifend,
			                           MEMiden);
		}
		cerr<<"\nLoad genomes: "<<refseqs.size()<<endl;
		cerr<<"\nchunk: "<<chunk<<endl;
		refdescr.clear(); startpos.clear();
		// Clear container.
		ref=""; 
		// Make part suffix array.
		make_block_ref(refseqs, ref, totalgenomes, refdescr, startpos);
		cerr<<"Creating suffix array ......\n"<<endl;
		saa = new paraSA(ref, refdescr, startpos, true, K);
		cerr<<"\nFinished creating suffix array ......\n"<<endl;
		//genomes=refseqs.size();

		// Part internal clustering || parallel part.
		pthread_attr_t attr;  pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		vector<threads_arg> args(total_threads);
		vector<pthread_t> thread_ids(total_threads);

		// Initialize additional thread data.
		for(int i=0; i<total_threads; i++)
		{
			args[i].skip = total_threads;
			args[i].skip0 = i;
			args[i].begin = cbegin;
			args[i].part = true;
			args[i].chunk = chunk; 
		}

		// Create joinable threads to find MEMs.
		for(int i=0; i<total_threads; i++) 
			pthread_create(&thread_ids[i], &attr, single_thread, (void *)&args[i]);
		// Wait for all threads to terminate.
		for(int i=0; i<total_threads; i++) pthread_join(thread_ids[i], NULL);
		// Collect clustering information into one chunk.
		getClusteringInfoOnepart(totalgenomes, cbegin, chunk, true, clusterhit);
		if (ifend) { break;} // Finished.

		if ( (rebuild)&&(clusterhit) )
		{
			delete saa; ref="";
			refdescr.clear(); 
			startpos.clear();	
			// Make part suffix array.
			make_block_ref(refseqs, ref, totalgenomes, refdescr, startpos);
			cerr<<"Creating suffix array ......\n"<<endl;
			saa = new paraSA(ref, refdescr, startpos, true, K);
			cerr<<"\nFinished creating suffix array ......\n"<<endl;
		}

		refseqs.clear();
		begin = cbegin+chunk;
		if (loadall) { dchunk=(long)allrefseqs.size();
		}else{ dchunk=(long)(chunk*Nchunk);}

		// Make alignment for last genomes.
		while (1)
		{
			cerr<<"\n=================="<<endl;
			cerr<<"begin alignment "<<begin<<"\n"<<endl;

			if (loadall) { 
				load_part_genomes_mem(allrefseqs, 
					                    refseqs, 
					                    totalgenomes, 
					                    begin, 
					                    dchunk);
			}else{
			  load_part_genomes(ref_fasta, 
					                refseqs, 
					                totalgenomes, 
					                begin, 
					                dchunk);
			}
			//Parallel part.
			pthread_attr_t attr;  pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			vector<threads_arg> args(total_threads);
			vector<pthread_t> thread_ids(total_threads);
			// Initialize additional thread data.
			for(int i = 0; i < total_threads; i++) 
			{
				args[i].skip = total_threads;
				args[i].skip0 = i;
				args[i].begin = cbegin;
				args[i].part = false;
				args[i].chunk = chunk;
			}
			//Create joinable threads to find MEMs.
			for(int i = 0; i < total_threads; i++) 
				pthread_create(&thread_ids[i], &attr, single_thread, (void *)&args[i]);
			//Wait for all threads to terminate.
			for(int i = 0; i < total_threads; i++) pthread_join(thread_ids[i], NULL);
			if ((long)refseqs.size()<dchunk){
				getClusteringInfoOnepart(totalgenomes, 
					                       begin, 
					                       (long)refseqs.size(), 
					                       false, 
					                       clusterhit);
			}else{
				getClusteringInfoOnepart(totalgenomes, 
					                       begin, 
					                       dchunk, 
					                       false, 
					                       clusterhit);
			}
			if ((long)refseqs.size()<dchunk) { break; }
			begin+=dchunk;
			refseqs.clear();
		}
		delete saa;
		refseqs.clear();
		cbegin=cbegin+chunk;

	}//end while(1)

	//testDistanceBgenomes(totalgenomes); 
	// Collect clustering information.
	cerr<<"\n==========================="<<endl; 
	cerr<<"Output clustering information ......\n"<<endl;
	// Output with CD-HIT format.
	outputClusteringInfoSimple(totalgenomes); 
	cerr<<"The finish.\n"<<endl;
        end=time(NULL);
        cerr<<"time:"<<end-start<<" s"<<endl;
	return 0; // The end.
}

void usage(string prog) 
{
	cerr << "Gclust is a clustering program for genome, draft assembly contigs,";
	cerr << " which algorithm is based on all Maximal Exact Matches(MEMs) between genome sequences." << endl;
	cerr << endl;
  cerr << "Usage: " << prog << " [options] <genomes-file> " << endl;
	cerr << endl;
  cerr << "Options:" << endl;
	cerr << endl;
  cerr << "-minlen        Set the minimum length for exact match, if not set, default = 20" << endl;
  cerr << "-both          Compute forward and reverse complement matches, default = forward" << endl;
  cerr << "-nuc           Match only the characters a, c, g, or t" << endl;
	cerr << "-sparse        Set the step of sparse suffix array, default = 1" <<endl;
  cerr << "-threads       Set the number of threads to use, default = 1" << endl;
	cerr << "-chunk         Set the chunk size for one time clustering, default = 100, where the unit is million base pairs (Mbp)" << endl;
	cerr << "-nchunk        Set the chunk number loaded one time for remaining genomes alignment, default = 2" << endl;
	cerr << "-loadall       Load the total genomes one time" << endl;
	cerr << "-rebuild       Rebuild suffix array after clustering into one chunk, default = 1" << endl;
	cerr << endl;
  cerr << "Clustering cutoff:" << endl;
	cerr << endl;
	cerr << "-memiden       Set the value of extended maximal exact match (MEM) idendity or non-extended MEM idendity for clustering, default = 90" << endl;
  cerr << endl;
	cerr << "Extension options of MEM:" << endl;
	cerr << endl;
	cerr << "-ext       Set the extension type of MEM, where '0' means no extension, '1' means gapped extension and '2' means un-gapped extension, default = 1" << endl;
	cerr << "-mas       Set the reward value for a nucleotide match, default = 1" << endl;
	cerr << "-umas      Set the penalty value for a nucleotide mismatch, default = -1" << endl;
	cerr << "-gapo      Set the cost value to open a gap, default = -1" << endl;
	cerr << "-gape      Set the cost value to extend a gap, default = -1" << endl;
	cerr << "-drops     Set the X dropoff value for extension, default = 1" << endl;
	cerr << endl;
  cerr << "Example usage:" << endl;
  cerr << endl;
  cerr << "./gclust -minlen 20 -both -nuc -threads 8 -ext 1 -sparse 2 data/viral.1.1.genomic.sort.fna > data/viral.1.1.genomic.sort.fna.clustering.out" << endl;
	cerr <<endl;
  cerr << "Find all MEMs on forward and reverse strands" << endl;
  cerr << "of length 20 or greater, matching only a, c, t, or g" << endl;
	cerr << "using 8 threads parallel computing" << endl;
	cerr << "gapped extension for MEMs seeds." <<endl;
  cerr << endl;
  exit(1);

}
