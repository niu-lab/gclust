#ifndef __paraSA_hpp__
#define __paraSA_hpp__

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <limits>

using namespace std;

#include "fasta.hpp"

// MUMI index cutoff unit.
struct mumi_unit
{
	long ref;
	long id;
	bool dgdelete;
	long g0sens; //rep num
	long g0init;
	long g0long;
	long g0fin;
	long g1sens; //rep num
	long g1init;
	long g1long;
	long g1fin;

};

// Order G0.
bool compareg0(const mumi_unit lhs, const mumi_unit rhs);
// Order G1.
bool compareg1(const mumi_unit lhs, const mumi_unit rhs);

// Compute genome identity.
bool ComputeMemIdentity(vector<GenomeClustInfo> &totalgenomes,
												vector<Genome> &allpartgenomes,
												vector<mumi_unit> &mumiunits,
												long beginclust,
												long id,
												int MEMiden,
												bool part,
												long chunk,
												char strand,
												int ext, //extension options
												int mas,
												int umas,
												int gapo,
												int gape,
												int drops);

void collectg0(vector<mumi_unit> &singleunits);
void collectg1(vector<mumi_unit> &singleunits);

void Remove_symetrically(vector<mumi_unit> &singleunits,
												 vector<mumi_unit> &removes);

void doublecollectg0(vector<mumi_unit> &removes);
void doublecollectg1(vector<mumi_unit> &removes);
void trimendg0(vector<mumi_unit> &dremoves);
void trimendg1(vector<mumi_unit> &dremoves);
void merging(vector<mumi_unit> &ddremoves, vector<mumi_unit> &mergeunits);
double tell_me(vector<mumi_unit> &ddremoves,long size);

// seed extension part
// Seed extension with two directions
void seedextensions(vector<mumi_unit> &ddremoves,
									 vector<Genome> &allpartgenomes,
									 long id, 
									 long iref, 
									 char strand,
									 int ext, //extension options
									 int mas,
									 int umas,
									 int gapo,
									 int gape,
									 int drops);

struct vec_uchar 
{
  struct item_t
	{
    item_t(size_t i, int v) { idx = i; val = v; }
    size_t idx; int val;
    bool operator < (item_t t) const { return idx < t.idx; }
  };
  vector<unsigned char> vec;  // LCP values from 0-65534
  vector<item_t> M;
  void resize(size_t N) { vec.resize(N); }
  
	// Vector X[i] notation to get LCP values.
  int operator[] (size_t idx) 
	{
    if(vec[idx] == numeric_limits<unsigned char>::max()) 
      return lower_bound(M.begin(), M.end(), item_t(idx,0))->val;
    else 
      return vec[idx]; 
  }
  // Actually set LCP values, distingushes large and small LCP
  // values.
  void set(size_t idx, int v) 
	{
		if(v >= numeric_limits<unsigned char>::max()) 
		{
      vec.at(idx) = numeric_limits<unsigned char>::max();
      M.push_back(item_t(idx, v));
    }else{ 
			vec.at(idx) = (unsigned char)v; 
		}
  }
	// Once all the values are set, call init. This will assure the
  // values >= 255 are sorted by index for fast retrieval.
  void init() 
	{ 
		sort(M.begin(), M.end()); 
		cerr << "M.size()=" << M.size() << endl;
	}

};

// Match find by findMEM. 
struct match_t 
{
  match_t() { query = 0, len = 0; refseq=0; refpos=0;}
  match_t(long q, long l, long m, long n) { 
		query = q; len = l; refseq=m; refpos=n; 
	}
	long query; // position in query
  long len; // length of match
	long refseq;
	long refpos;

};

// Match find by findMEM. 
struct postmatch_t 
{
	long id; // reference id
  long ref; // position in reference sequence
  long query; // position in query
  long len; // length of match
	bool rc;  // direction

};

// depth : [start...end] 
struct interval_t 
{
	interval_t() { start = 1; end = 0; depth = -1; }
  interval_t(long s, long e, long d) { start = s; end = e; depth = d; }
  void reset(long e) { start = 0; end = e; depth = 0; }
  long depth, start, end;
  long size() { return end - start + 1; }

};

// Suffix array.
struct paraSA 
{
	vector<long> &descr; // Descriptions of concatenated sequences.
  vector<long> &startpos; // Lengths of concatenated sequences.
  long maxdescrlen; // Maximum length of the sequence description, used for formatting.
  bool _4column; // Use 4 column output format.
  long N; //!< Length of the sequence.
  long logN; // ceil(log(N)) 
  long NKm1; // N/K - 1
  string &S; //!< Reference to sequence data.
  vector<unsigned int> SA;  // Suffix array.
  vector<int> ISA;  // Inverse suffix array.
	vec_uchar LCP; // Simulates a vector<int> LCP.
  long K; // suffix sampling, stable K = 1.

	// Suffix sort part.
	// New qsort in class
	int *I;   /* group array, ultimately suffix array.*/
	int *V;   /* inverse array, ultimately inverse of I.*/
	int r;    /* number of symbols aggregated by transform.*/
	int h;    /* length of already-sorted prefixes.*/
	// Change them to inline function
	int KEY(int *p)
	{
		return(V[*(p)+(h)]);
	}
	void SWAP(int *p, int *q)
	{
		int tmp;
		tmp=*(p);
		*(p)=*(q);
		*(q)=tmp;
	}
	int* MED3(int *a, int *b, int *c)
	{
		return(KEY(a)<KEY(b) ?                        \
        (KEY(b)<KEY(c) ? (b) : KEY(a)<KEY(c) ? (c) : (a))       \
        : (KEY(b)>KEY(c) ? (b) : KEY(a)>KEY(c) ? (c) : (a)));
	}
	// End of suffix sort part.
  // Maps a hit in the concatenated sequence set to a position in that sequence.
  void from_set(long hit, long &seq, long &seqpos) 
	{
		// Use binary search to locate index of sequence and position
    // within sequence.
    vector<long>::iterator it = upper_bound(startpos.begin(), startpos.end(), hit);
    seq = distance(startpos.begin(), it) - 1;
    it--;
    seqpos = hit - *it;
	}
  // Constructor builds sparse suffix array. 
  paraSA(string &S_, 
				 vector<long> &descr_,
				 vector<long> &startpos_,
				 bool __4column, 
				 long K_);

	// suffix sort part
	void suffixsort(int *x, int *p, int n, int k, int l);
	void update_group(int *pl, int *pm);
	void select_sort_split(int *p, int n);
	int choose_pivot(int *p, int n);
	void sort_split(int *p, int n);
	void bucketsort(int *x, int *p, int n, int k);
	int transform(int *x, int *p, int n, int k, int l, int q);

  // Modified Kasai et all for LCP computation.
  void computeLCP();

  // Radix sort required to construct transformed text for sparse SA construction.
  void radixStep(int *t_new, 
								 int *SA, 
								 long &bucketNr,
								 long *BucketBegin,
								 long l,
								 long r,
								 long h);
  
	// Binary search for left boundry of interval.
  inline long bsearch_left(char c, long i, long s, long e);
  // Binary search for right boundry of interval.
  inline long bsearch_right(char c, long i, long s, long e);
  // Simple suffix array search.
  inline bool search(string &P, long &start, long &end);

	// Simple top down traversal of a suffix array.
  inline bool top_down(char c, long i, long &start, long &end);
  inline bool top_down_faster(char c, long i, long &start, long &end);

  // Traverse pattern P starting from a given prefix and interval
  // until mismatch or min_len characters reached.
  inline void traverse(string &P, long prefix, interval_t &cur, int min_len);
  
	// Simulate a suffix link.
  inline bool suffixlink(interval_t &m);
  // Expand ISA/LCP interval. Used to simulate suffix links.
  inline bool expand_link(interval_t &link) 
	{
    // Threshold link expansion.
		long thresh = 2 * link.depth * logN, exp = 0; 
    long start = link.start;
    long end = link.end;
    while(LCP[start] >= link.depth) 
		{ 
      exp++; 
      if(exp >= thresh) return false; 
      start--; 
    }
    while(end < NKm1 && LCP[end+1] >= link.depth) 
		{ 
      exp++; 
      if(exp >= thresh) return false; 
      end++; 
    }
    link.start = start; link.end = end;
    return true;
  }

  // Given a position i in S, finds a left maximal match of minimum
  // length within stable 1 step.
	inline void find_Lmaximal(string &P, 
														long prefix,
														long i,
														long len,
														vector<match_t> &matches,
														int min_len,
														long id);

  // Given an interval where the given prefix is matched up to a
  // mismatch, find all MEMs up to a minimum match depth.
  void collectMEMs(string &P,
									 long prefix,
									 interval_t mli,
									 interval_t xmi,
									 vector<match_t> &matches,
									 int min_len,
									 long id);
	// 100% match
  void collectMEMsperfect(string &P, 
													long prefix, 
													interval_t mli, 
													interval_t xmi, 
													vector<match_t> &matches, 
													int min_len, 
													long id);

  // Find all MEMs given a prefix pattern offset k.
  void findMEM(long k, 
							 string &P,
							 vector<match_t> &matches,
							 int min_len,
							 long id);

	void findMEMperfect(long k, 
											string &P, 
											vector<match_t> &matches, 
											int min_len, 
											long id);
  
	//void findMAM(string &P, vector<match_t> &matches, int min_len, bool print);
  inline bool is_leftmaximal(string &P, long p1, long p2);
  // Find Maximal Exact Matches (MEMs) 
  void MEM(string &P, vector<match_t> &matches, int min_len, long id);
	void MEMperfect(string &P, vector<match_t> &matches, int min_len, long id);

  // Maximal Unique Match (MUM) 
  // void MUM(string &P, vector<match_t> &unique, int min_len, bool print);
	long load_match_info(long id, 
											 vector<match_t> &buf,
											 vector<mumi_unit> &mumiunits,
											 bool rc,
											 long qlen);

};

#endif // __paraSA_hpp__
