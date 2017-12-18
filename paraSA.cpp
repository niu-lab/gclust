
#include<seqan/seeds.h>
//#include<seqan/include/seqan/seeds.h>
//#include<seqan/include/seqan/sequence.h>
#include<cstdlib>
#include<stdio.h>
#include <math.h>
#include <pthread.h>
#include <limits.h>

#include "fasta.hpp"
#include "paraSA.hpp"

using namespace seqan;

// seed extension part | seed content test
template<typename TSeed, typename TSeq>
void writeSeed(TSeed & seed, TSeq const & seq0, TSeq const & seq1)
{
	cout << "Ref seed from position "<<leftPosition(seed, 0);
	cout << " to " << rightPosition(seed, 0)<<" : ";
	cout << infix(seq0, leftPosition(seed, 0), rightPosition(seed, 0)+1)<<endl;
	cout << "query seed from position " <<leftPosition(seed, 1);
	cout << " to " << rightPosition(seed, 1)<<" : ";
	cout << infix(seq1, leftPosition(seed, 1), rightPosition(seed, 1)+1)<<endl;
	cout <<endl;
}

// Order G0 | based on ref
bool compareg0(const mumi_unit lhs, const mumi_unit rhs)
{
	if( lhs.g0init != rhs.g0init ) return lhs.g0init < rhs.g0init;
	if( lhs.g0long != rhs.g0long ) return rhs.g0long < lhs.g0long;
	return lhs.g1init < rhs.g1init;

};

// Order G1 | based on query
bool compareg1(const mumi_unit lhs, const mumi_unit rhs)
{
	if( lhs.g1init != rhs.g1init ) return lhs.g1init < rhs.g1init;
	if( lhs.g1long != rhs.g1long ) return rhs.g1long < lhs.g1long;
	return lhs.g0init < rhs.g0init;

};

// Compute Mems identity.
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
												int drops)
{

	double miniden;
	double distance;
	bool hitted = false;
	long addsize, sid;

	vector<mumi_unit> singleunits;
	miniden=(double)MEMiden/100;

	if (!part){
		sid=beginclust+chunk;
	}else{
		sid=id;
	}

	for (long i=beginclust;i<sid;i++)
	{

		addsize=0;
		singleunits.clear();

		// Collect single ref genome matches.
		for (long j=0;j<(long)mumiunits.size();j++){
			if (mumiunits[j].ref==i){
				singleunits.push_back(mumiunits[j]);
				addsize+=mumiunits[j].g0long;
			}
		}

		//if ((double)addsize/totalgenomes[id].size < miniden ){ continue; }

		// Sort matches base on ref genome.
		sort(singleunits.begin(), singleunits.end(),compareg0);
		collectg0(singleunits);

		// Sort matches base on query genome.
		sort(singleunits.begin(), singleunits.end(),compareg1);

		collectg1(singleunits);
		vector<mumi_unit> removes;

		// Remove symetrical matches.
		Remove_symetrically(singleunits,removes);
		sort(removes.begin(), removes.end(), compareg0);

		doublecollectg0(removes);
		sort(removes.begin(), removes.end(), compareg1);

		doublecollectg1(removes);
		vector<mumi_unit> dremoves;

		Remove_symetrically(removes,dremoves);
		sort(dremoves.begin(), dremoves.end(), compareg0);
		trimendg0(dremoves);

		//preG1_postG0
		sort(dremoves.begin(), dremoves.end(),compareg1);
		vector<mumi_unit> ddremoves;
		Remove_symetrically(dremoves,ddremoves);

		//treat_chG1
		sort(ddremoves.begin(), ddremoves.end(),compareg1);
		trimendg1(ddremoves);

		//continue if no unit
		if ((long)ddremoves.size()==0){ continue; }

		//merge near units for extension
		vector<mumi_unit> mergeunits;
	  merging(ddremoves, mergeunits);

		if (ext!=0){
			//seed extension part | checking and extending
			seedextensions(mergeunits, allpartgenomes, id, i, strand, 
										 ext, mas, umas, gapo, gape, drops);
		}

		distance=tell_me(mergeunits,totalgenomes[id].size);

		//loading clustering info
		if (distance >= miniden){
			hit thit;
			thit.id=i;
			thit.strand=strand;
			thit.identity=distance;
			totalgenomes[id].clusters.push_back(thit);
			hitted=true;
			if (!part){ break; } // Hitted and exit
		} // Do MuMi computing.

	}

	return hitted;

}

// Collect ref matches.
void collectg0(vector<mumi_unit> &singleunits)
{

	long j, i=0;
	while (i<(long)(singleunits.size()-1))
	{
		if (singleunits[i].g0fin >= singleunits[i+1].g0fin)
		{
			singleunits[i+1].dgdelete=true;
			j=i;
			j++;
			while ((j<(long)(singleunits.size()-1))&&\
				(singleunits[i].g0fin >= singleunits[j+1].g0fin))
			{
				singleunits[j+1].dgdelete=true;
				j++;
			}
			i=j+1;
		}else{
			i++;
		}
	}

}


// Collect query matches.
void collectg1(vector<mumi_unit> &singleunits)
{
	long j, i=0;
	while (i<(long)(singleunits.size()-1))
	{
		if (singleunits[i].g1fin >= singleunits[i+1].g1fin)
		{
			singleunits[i+1].dgdelete=true;
			j=i;
			j++;
			while ((j<(long)(singleunits.size()-1))&&\
				(singleunits[i].g1fin >= singleunits[j+1].g1fin))
			{
				singleunits[j+1].dgdelete=true;
				j++;
			}
			i=j+1;
		}else{
			i++;
		}
	}
}


// Remove symetrical matches.
void Remove_symetrically(vector<mumi_unit> &singleunits,
												 vector<mumi_unit> &removes)
{
	for (long i=0;i<(long)singleunits.size();i++)
	{
		if (!singleunits[i].dgdelete){
			removes.push_back(singleunits[i]);
		}
	}
}

// Double collect ref matches.
void doublecollectg0(vector<mumi_unit> &removes)
{
	for (long i=removes.size()-2;i>0;i--)
	{
		if ( (removes[i-1].g0fin >= removes[i+1].g0init-1)&&\
			(removes[i].g0fin <= removes[i+1].g0fin))
		{
			removes[i].dgdelete=true;
		}
	}
}

// Double collecting ...
void doublecollectg1(vector<mumi_unit> &removes)
{
	for (long i=removes.size()-2;i>0;i--)
	{
		if ( (removes[i-1].g1fin >= removes[i+1].g1init-1)&&\
			(removes[i].g1fin <= removes[i+1].g1fin))
		{
			removes[i].dgdelete=true;
		}

	}
}

// Trim match part.
void trimendg0(vector<mumi_unit> &dremoves)
{
		long cptrim=0, tmplen=0;
		for (long i=(dremoves.size()-2);i>=0;i--)
		{
			if (dremoves[i].g0fin >= dremoves[i+1].g0init)
			{
				 cptrim++;
				 tmplen=dremoves[i].g0fin-dremoves[i+1].g0init+1;
				 dremoves[i].g0fin=dremoves[i+1].g0init-1;
				 dremoves[i].g0long=dremoves[i].g0fin-dremoves[i].g0init+1;
				 if (dremoves[i].g1sens==1){
					dremoves[i].g1fin=dremoves[i].g1fin-tmplen;
					dremoves[i].g1long=dremoves[i].g1fin-dremoves[i].g1init+1;
				 }else{
					dremoves[i].g1init=dremoves[i].g1init+tmplen;
					dremoves[i].g1long=dremoves[i].g1fin-dremoves[i].g1init+1;
				 }
			}
		}

}

// Trim ...
void trimendg1(vector<mumi_unit> &dremoves)
{
		long cptrim=0, tmplen=0;
		for (long i=(dremoves.size()-2);i>=0;i--)
		{
			if (dremoves[i].g1fin >= dremoves[i+1].g1init)
			{
				 cptrim++;
				 tmplen=dremoves[i].g1fin-dremoves[i+1].g1init+1;
				 dremoves[i].g1fin=dremoves[i+1].g1init-1;
				 dremoves[i].g1long=dremoves[i].g1fin-dremoves[i].g1init+1;
				 if (dremoves[i].g1sens==1){
					dremoves[i].g0fin=dremoves[i].g0fin-tmplen;
					dremoves[i].g0long=dremoves[i].g0fin-dremoves[i].g0init+1;
				 }else{
					dremoves[i].g0init=dremoves[i].g0init+tmplen;
					dremoves[i].g0long=dremoves[i].g0fin-dremoves[i].g0init+1;
				 }
			}
		}

}


//merge near units.
void merging(vector<mumi_unit> &ddremoves,
						 vector<mumi_unit> &mergeunits)
{

	//merge units gap <= 3bps
	mumi_unit tmu, tmu0, tmu1;
	tmu = ddremoves[0];

	//only single unit
	if (ddremoves.size()==1)
	{
		mergeunits.push_back(tmu);
		return;
	}

	//merge part
	for (long i=1; i<(long)ddremoves.size();i++)
	{
		tmu0=ddremoves[i];
		if ( (((tmu0.g0init-tmu.g0fin)>0)&&((tmu0.g0init-tmu.g0fin)<=4))&&\
			((tmu0.g1init-tmu.g1fin)<=4))
		{
			tmu1=tmu;
			tmu1.g0long=tmu.g0long+tmu0.g0long;
			tmu1.g0fin=tmu0.g0fin;
			tmu1.g1long=tmu1.g0long;
			tmu1.g1fin=tmu0.g1fin;
			tmu=tmu1;
		}else{
			mergeunits.push_back(tmu);
			tmu=tmu0;
		}
	}
	//last one unit
	mergeunits.push_back(tmu);

	return;

}

// Get distance.
double tell_me(vector<mumi_unit> &ddremoves,long size)
{
	long add=0;
	for (long i=0;i<(long)ddremoves.size();i++)
	{
		add+=ddremoves[i].g1long;
	}
	return (double)add/size;
}

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
									 int drops)

{
	// Total Units 
	long units=(long)ddremoves.size();
	long g1leftband, g1rightband;
	long g1leftbegin=0;
	long g0infixbegin, g1infixbegin;
	long g0infixend, g1infixend;
	long seedfrom, seedto;
	long extensionlen;
	long g1rightend;
  mumi_unit tmu0;
	
	// For right extension
	//g1rightend = allpartgenomes[id].size-1;

	// Get DNA content
	DnaString seq0, seq1;
	if(allpartgenomes.size() != 0) 
	{//add judgment statement for variable allpartgenomes
		seq0 = allpartgenomes[iref].cont;
		seq1 = allpartgenomes[id].cont;
	}
	else
	{
		return; 
	}


	// Get reverse complement DNA sequence 
	if (strand == '-')
	{
		reverseComplement(seq1);
	}

	// left seed extension part
	for (long i=0; i<units;i++)
	{

		/*
		cout<<ddremoves[i].ref<<"\t"<<
					id<<"\t"<<
					ddremoves[i].g0sens<<"\t"<<
					ddremoves[i].g0init<<"\t"<<
					ddremoves[i].g0long<<"\t"<<
					ddremoves[i].g0fin<<"\t"<<
					ddremoves[i].g1sens<<"\t"<<
					ddremoves[i].g1init<<"\t"<<
					ddremoves[i].g1long<<"\t"<<
					ddremoves[i].g1fin<<"\t"<<" | ";
		//cout<<endl;
		*/

		tmu0 = ddremoves[i];
		//left extension
		if ( (tmu0.g1init-2) > g1leftbegin )
		{
			g1infixbegin = g1leftbegin;
			//length of candidate extension part
			g1leftband=tmu0.g1init-g1leftbegin-1;
			if (tmu0.g0init > (g1leftband+10))
			{
				g0infixbegin=tmu0.g0init-g1leftband-10-1;
			}else{
				g0infixbegin=0;
			}

			g1infixend=tmu0.g1init;
			g0infixend=tmu0.g0init;
			//get infix part instead whole sequence for extension
			typedef Infix<DnaString>::Type TInfix;
			TInfix infix0 = infix(seq0, g0infixbegin, g0infixend);
			TInfix infix1 = infix(seq1, g1infixbegin, g1infixend);
			
			long seedbegin0;
			if (g0infixbegin==0)
			{
				seedbegin0=g0infixend-1;
			}else{
				seedbegin0=g1leftband+10;
			}

			typedef Seed<> TSeed;
			TSeed seed(seedbegin0, g1leftband, 1);
			seedfrom=leftPosition(seed, 1);
			//extension parameters 
			typedef int TScore;
			
			//Score<TScore, Simple> scoreMatrix(1, -1, -1);
			//TScore scoreDropOff = 1;
			Score<TScore, Simple> scoreMatrix(mas, umas, gape, gapo);
			TScore scoreDropOff = drops;

			//Gapped or ungapped extension
			if (ext==1){
				extendSeed(seed, scoreDropOff, scoreMatrix, infix0, infix1, 0, GappedXDrop());
			}else{
				extendSeed(seed, scoreDropOff, scoreMatrix, infix0, infix1, 0, UngappedXDrop());
			}

			seedto=leftPosition(seed, 1);
			extensionlen=seedfrom-seedto;
			ddremoves[i].g1init-=extensionlen;
			ddremoves[i].g1long+=extensionlen;

			//global alignment of extension part | if need
			/*
			Align<Infix<DnaString>::Type > align;
			resize(rows(align), 2);
			assignSource(row(align, 0), infix(infix0, leftPosition(seed, 0), \
				rightPosition(seed, 0)+1));
			assignSource(row(align, 1), infix(infix1, leftPosition(seed, 1), \
				rightPosition(seed, 1)+1));
			std::cout << std::endl << "Banded Alignment:" << std::endl;
			std::cout << "Score: " << globalAlignment(align, stringSet(align), scoreMatrix, \
				-leftDiagonal(seed) - 2, -rightDiagonal(seed) + 2, \
				BandedNeedlemanWunsch()) << std::endl;
			std::cout << align;
			*/
		}

		g1leftbegin = tmu0.g1fin;

	}
	
	//right extension part 
	g1rightend = allpartgenomes[id].size-1;
	for (long i=(units-1); i>=0; i--)
	{
		tmu0 = ddremoves[i];
		//right extension 
		//do extension only when gap >=3
		if ( (tmu0.g1fin+2) < g1rightend ) 
		{
			g1infixbegin = tmu0.g1fin-1;
			g0infixbegin = tmu0.g0fin-1;
			g1rightband=g1rightend-tmu0.g1fin+2;

			g0infixend = g0infixbegin+g1rightband+10;
			g1infixend = g1infixbegin+g1rightband;
			
			typedef Infix<DnaString>::Type TInfix;
			TInfix infix0 = infix(seq0, g0infixbegin, g0infixend);
			TInfix infix1 = infix(seq1, g1infixbegin, g1infixend);
		
			typedef Seed<> TSeed;
			TSeed seed(0, 0, 1);

			//extension parameters 
			typedef int TScore;
			//Score<TScore, Simple> scoreMatrix(1, -1, -1);
			//TScore scoreDropOff = 1;
			Score<TScore, Simple> scoreMatrix(mas, umas, gape, gapo);
			TScore scoreDropOff = drops;
			//Gapped or ungapped extension
			if (ext==1){
				extendSeed(seed, scoreDropOff, scoreMatrix, infix0, infix1, 1, GappedXDrop());
			}else{
				extendSeed(seed, scoreDropOff, scoreMatrix, infix0, infix1, 1, UngappedXDrop());
			}

			extensionlen=rightPosition(seed, 1);
			//update unit info
			ddremoves[i].g1fin+=extensionlen;
			ddremoves[i].g1long+=extensionlen;
		}
		g1rightend = tmu0.g1init-2;
	}

}

paraSA::paraSA(string &S_,
							 vector<long> &descr_,
							 vector<long> &startpos_,
							 bool __4column, 
							 long K_) : descr(descr_), startpos(startpos_), S(S_) 
{
  _4column = __4column;
  K = K_;

	if(S.length() % K != 0) {
    long appendK = K - S.length() % K ;
    for(long i = 0; i < appendK; i++) S += '$';
  }
  // Make sure last K-sampled characeter is this special character as well!!
	// Append "special" end character. Note: It must be lexicographically less.
  for(long i = 0; i < K; i++) S += '$'; 
  N = S.length();

		// Sparse suffix array construction part
	if(K > 1) 
	{
    long bucketNr = 1;
    int *intSA = new int[N/K+1];  for(int i = 0; i < N/K; i++) intSA[i] = i; // Init SA.
    int* t_new = new int[N/K+1];
    long* BucketBegin = new long[256]; // array to save current bucket beginnings
    radixStep(t_new, intSA, bucketNr, BucketBegin, 0, N/K-1, 0); // start radix sort
    t_new[N/K] = 0; // Terminate new integer string.
    delete[] BucketBegin;
    // Suffix sort integer text.
    cerr<<"suffix sorting ...."<<endl;
    suffixsort(t_new, intSA, N/K, bucketNr, 0);
    cerr <<"suffix sorting done ...." << endl;
    delete[] t_new;
    // Translate suffix array. 
    SA.resize(N/K);
    for (long i=0; i<N/K; i++) SA[i] = (unsigned int)intSA[i+1] * K;
    delete[] intSA;
    // Build ISA using sparse SA. 
    ISA.resize(N/K);             
    for(long i = 0; i < N/K; i++) { ISA[SA[i]/K] = i; }

  }else{

		SA.resize(N);
		ISA.resize(N);
		int char2int[UCHAR_MAX+1]; // Map from char to integer alphabet.
		// Zero char2int mapping.
		for (int i=0; i<=UCHAR_MAX; i++) char2int[i]=0;
		// Determine which characters are used in the string S.
		for (long i = 0; i < N; i++) {	
			char2int[(int)S[i]]=1;
		}
		// Count the size of the alphabet. 
		int alphasz = 0;
		for(int i=0; i <= UCHAR_MAX; i++) {
			if (char2int[i]) char2int[i]=alphasz++;
			else char2int[i] = -1;
		}
		// Remap the alphabet. 
		for(long i = 0; i < N; i++) ISA[i] = (int)S[i]; 
		for (long i = 0; i < N; i++) 	{
			ISA[i]=char2int[ISA[i]] + 1; 
		}
		// First "character" equals 1 because of above plus one, l=1 in suffixsort(). 
		int alphalast = alphasz + 1;
		// Use LS algorithm to construct the suffix array.
		int *SAint = (int*)(&SA[0]);
		cerr<<"suffix sorting ...."<<endl;
		suffixsort(&ISA[0], SAint , N-1, alphalast, 1);
		cerr <<"suffix sorting done ...." << endl;

	}

  cerr << "N=" << N << endl;
  // Adjust to "sampled" size. 
  logN = (long)ceil(log(N/K) / log(2.0));
	cerr<<"logN \t"<<logN<<endl;
  LCP.resize(N/K);
  cerr << "N/K=" << N/K << endl;
  // Use algorithm by Kasai et al to construct LCP array.
  computeLCP();  // SA + ISA -> LCP
  LCP.init();
  NKm1 = N/K-1;
	cerr<< "NKm1= "<<NKm1<<endl;

}

// Suffix sort part.
void paraSA::update_group(int *pl, int *pm)
{
   int g;
   g=pm-I;                      /* group number.*/
   V[*pl]=g;                    /* update group number of first position.*/
   if (pl==pm)
      *pl=-1;                   /* one element, sorted group.*/
   else
      do                        /* more than one element, unsorted group.*/
         V[*++pl]=g;            /* update group numbers.*/
      while (pl<pm);
}

void paraSA::select_sort_split(int *p, int n) 
{
   int *pa, *pb, *pi, *pn;
	 int f, v;
   pa=p;                        /* pa is start of group being picked out.*/
   pn=p+n-1;                    /* pn is last position of subarray.*/
   while (pa<pn) {
      for (pi=pb=pa+1, f=KEY(pa); pi<=pn; ++pi)
         if ((v=KEY(pi))<f) {
            f=v;                /* f is smallest key found.*/
            SWAP(pi, pa);       /* place smallest element at beginning.*/
            pb=pa+1;            /* pb is position for elements equal to f.*/
         } else if (v==f) {     /* if equal to smallest key.*/
            SWAP(pi, pb);       /* place next to other smallest elements.*/
            ++pb;
         }
      update_group(pa, pb-1);   /* update group values for new group.*/
      pa=pb;                    /* continue sorting rest of the subarray.*/
   }
   if (pa==pn) {                /* check if last part is single element.*/
      V[*pa]=pa-I;
      *pa=-1;                   /* sorted group.*/
   }
}

int paraSA::choose_pivot(int *p, int n) 
{
   int *pl, *pm, *pn;
   int s;
   
   pm=p+(n>>1);                 /* small arrays, middle element.*/
   if (n>7) {
      pl=p;
      pn=p+n-1;
      if (n>40) {               /* big arrays, pseudomedian of 9.*/
         s=n>>3;
         pl=MED3(pl, pl+s, pl+s+s);
         pm=MED3(pm-s, pm, pm+s);
         pn=MED3(pn-s-s, pn-s, pn);
      }
      pm=MED3(pl, pm, pn);      /* midsize arrays, median of 3.*/
   }
   return KEY(pm);
}

void paraSA::sort_split(int *p, int n)
{
   int *pa, *pb, *pc, *pd, *pl, *pm, *pn;
	 int f, v, s, t;
   if (n<7) {                   /* multi-selection sort smallest arrays.*/
      select_sort_split(p, n);
      return;
   }
   v=choose_pivot(p, n);
   pa=pb=p;
   pc=pd=p+n-1;
   while (1) {                  /* split-end partition.*/
      while (pb<=pc && (f=KEY(pb))<=v) {
         if (f==v) {
            SWAP(pa, pb);
            ++pa;
         }
         ++pb;
      }
      while (pc>=pb && (f=KEY(pc))>=v) {
         if (f==v) {
            SWAP(pc, pd);
            --pd;
         }
         --pc;
      }
      if (pb>pc)
         break;
      SWAP(pb, pc);
      ++pb;
      --pc;
   }
   pn=p+n;
   if ((s=pa-p)>(t=pb-pa))
      s=t;
   for (pl=p, pm=pb-s; s; --s, ++pl, ++pm)
      SWAP(pl, pm);
   if ((s=pd-pc)>(t=pn-pd-1))
      s=t;
   for (pl=pb, pm=pn-s; s; --s, ++pl, ++pm)
      SWAP(pl, pm);

   s=pb-pa;
   t=pd-pc;
   if (s>0)
      sort_split(p, s);
   update_group(p+s, p+n-t-1);
   if (t>0)
      sort_split(p+n-t, t);
}

void paraSA::bucketsort(int *x, int *p, int n, int k)
{
   int *pi, i, c, d, g;

   for (pi=p; pi<p+k; ++pi)
      *pi=-1;                   /* mark linked lists empty.*/
   for (i=0; i<=n; ++i) {
      x[i]=p[c=x[i]];           /* insert in linked list.*/
      p[c]=i;
   }
   for (pi=p+k-1, i=n; pi>=p; --pi) {
      d=x[c=*pi];               /* c is position, d is next in list.*/
      x[c]=g=i;                 /* last position equals group number.*/
      if (d>=0) {               /* if more than one element in group.*/
         p[i--]=c;              /* p is permutation for the sorted x.*/
         do {
            d=x[c=d];           /* next in linked list.*/
            x[c]=g;             /* group number in x.*/
            p[i--]=c;           /* permutation in p.*/
         } while (d>=0);
      } else
         p[i--]=-1;             /* one element, sorted group.*/
   }
}

int paraSA::transform(int *x, int *p, int n, int k, int l, int q)
{
   int b, c, d, e, i, j, m, s;
   int *pi, *pj;
   
   for (s=0, i=k-l; i; i>>=1)
      ++s;                      /* s is number of bits in old symbol.*/
   e=INT_MAX>>s;                /* e is for overflow checking.*/
   for (b=d=r=0; r<n && d<=e && (c=d<<s|(k-l))<=q; ++r) {
      b=b<<s|(x[r]-l+1);        /* b is start of x in chunk alphabet.*/
      d=c;                      /* d is max symbol in chunk alphabet.*/
   }
   m=(1<<(r-1)*s)-1;            /* m masks off top old symbol from chunk.*/
   x[n]=l-1;                    /* emulate zero terminator.*/
   if (d<=n) {                  /* if bucketing possible, compact alphabet.*/
      for (pi=p; pi<=p+d; ++pi)
         *pi=0;                 /* zero transformation table.*/
      for (pi=x+r, c=b; pi<=x+n; ++pi) {
         p[c]=1;                /* mark used chunk symbol.*/
         c=(c&m)<<s|(*pi-l+1);  /* shift in next old symbol in chunk.*/
      }
      for (i=1; i<r; ++i) {     /* handle last r-1 positions.*/
         p[c]=1;                /* mark used chunk symbol.*/
         c=(c&m)<<s;            /* shift in next old symbol in chunk.*/
      }
      for (pi=p, j=1; pi<=p+d; ++pi)
         if (*pi)
            *pi=j++;            /* j is new alphabet size.*/
      for (pi=x, pj=x+r, c=b; pj<=x+n; ++pi, ++pj) {
         *pi=p[c];              /* transform to new alphabet.*/
         c=(c&m)<<s|(*pj-l+1);  /* shift in next old symbol in chunk.*/
      }
      while (pi<x+n) {          /* handle last r-1 positions.*/
         *pi++=p[c];            /* transform to new alphabet.*/
         c=(c&m)<<s;            /* shift right-end zero in chunk.*/
      }
   } else {                     /* bucketing not possible, don't compact.*/
      for (pi=x, pj=x+r, c=b; pj<=x+n; ++pi, ++pj) {
         *pi=c;                 /* transform to new alphabet.*/
         c=(c&m)<<s|(*pj-l+1);  /* shift in next old symbol in chunk.*/
      }
      while (pi<x+n) {          /* handle last r-1 positions.*/
         *pi++=c;               /* transform to new alphabet.*/
         c=(c&m)<<s;            /* shift right-end zero in chunk.*/
      }
      j=d+1;                    /* new alphabet size.*/
   }
   x[n]=0;                      /* end-of-string symbol is zero.*/
   return j;                    /* return new alphabet size.*/
}

// LS suffix sorter (integer alphabet). 
void paraSA::suffixsort(int *x, int *p, int n, int k, int l)
{

	int *pi, *pk;
	int i, j, s, sl;
	V=x;                         /* set global values.*/
	I=p;
	if (n>=k-l) 
	{                /* if bucketing possible,*/
		j=transform(V, I, n, k, l, n);
		bucketsort(V, I, n, j);   /* bucketsort on first r positions.*/
	} else 
	{
		transform(V, I, n, k, l, INT_MAX);
		for (i=0; i<=n; ++i) 
			I[i]=i;                /* initialize I with suffix numbers.*/
		h=0;
		sort_split(I, n+1);       /* quicksort on first r positions.*/
	}
	h=r;                         /* number of symbols aggregated by transform.*/
	while (*I>=-n) 
	{
		pi=I;                     /* pi is first position of group.*/
		sl=0;                     /* sl is negated length of sorted groups.*/
		do {
			if ((s=*pi)<0) 
			{
				pi-=s;              /* skip over sorted group.*/
				sl+=s;              /* add negated length to sl.*/
			} else {
				if (sl) 
				{
					*(pi+sl)=sl;     /* combine sorted groups before pi.*/
					sl=0;
				}
				pk=I+V[s]+1;        /* pk-1 is last position of unsorted group.*/
				sort_split(pi, pk-pi);
				pi=pk;              /* next group.*/
			}
		} while (pi<=I+n);
		if (sl)                   /* if the array ends with a sorted group.*/
			*(pi+sl)=sl;           /* combine sorted groups at end of I.*/
		h=2*h;                    /* double sorted-depth.*/
	}

	for (i=0; i<=n; ++i)         /* reconstruct suffix array from inverse.*/
		I[V[i]]=i;

}
// End of suffix sort code.

// Uses the algorithm of Kasai et al 2001 which was described in
// Manzini 2004 to compute the LCP array. Modified to handle sparse
// suffix arrays and inverse sparse suffix arrays.
void paraSA::computeLCP() 
{
  long h=0;
  for(long i = 0; i < N; i+=K) 
	{ 
    long m = ISA[i/K]; 
    if(m==0) LCP.set(m, 0); // LCP[m]=0;
    else{
      long j = SA[m-1];
      while(i+h < N && j+h < N && S[i+h] == S[j+h])  h++;
      LCP.set(m, h); //LCP[m] = h;
    }
    h = max(0L, h - K);
  }
}

// Implements a variant of American flag sort (McIlroy radix sort).
// Recurse until big-K size prefixes are sorted. Adapted from the C++
// source code for the wordSA implementation from the following paper:
// Ferragina and Fischer. Suffix Arrays on Words. CPM 2007.
void paraSA::radixStep(int *t_new, 
											 int *SA,
											 long &bucketNr,
											 long *BucketBegin,
											 long l,
											 long r,
											 long h) 
{
  if(h >= K) return;
  // first pass: count
  vector<long> Sigma(256, 0); // Sigma counts occurring characters in bucket
	// count characters
  for (long i = l; i <= r; i++) Sigma[ S[ SA[i]*K + h ] ]++; 
  BucketBegin[0] = l; 
	for (long i = 1; i < 256; i++) 
	{ 
		BucketBegin[i] = Sigma[i-1] + BucketBegin[i-1]; 
	} // accumulate
  // second pass: move (this variant does *not* need an additional array!)
  unsigned char currentKey = 0;    // character of current bucket
  long end = l-1+Sigma[currentKey]; // end of current bucket
  long pos = l;                     // 'pos' is current position in bucket
  while (1) 
	{
    if (pos > end) 
		{ // Reached the end of the bucket.
      if (currentKey == 255) break; // Last character?
      currentKey++; // Advance to next characer.
      pos = BucketBegin[currentKey]; // Next bucket start.
      end += Sigma[currentKey]; // Next bucket end.
    }else{
      // American flag sort of McIlroy et al. 1993. BucketBegin keeps
      // track of current position where to add to bucket set.
      int tmp = SA[ BucketBegin[ S[ SA[pos]*K + h ] ] ]; 
			// Move bucket beginning to the right, and replace 
      SA[ BucketBegin[ S[ SA[pos]*K + h] ]++ ] = SA[pos];  
      SA[ pos ] = tmp; // Save value at bucket beginning.
			// Advance to next position if the right character.
      if (S[ SA[pos]*K + h ] == currentKey) pos++; 
    }
  }
  // recursively refine buckets and calculate new text:
  long beg = l; end = l-1;
  for (long i = 1; i < 256; i++) 
	{ // step through Sigma to find bucket borders
		end += Sigma[i];
		if (beg <= end) 
		{
			if(h == K-1) 
			{
				for (long j = beg; j <= end; j++) 
				{
					t_new[ SA[j] ] = bucketNr; // set new text
				}
				bucketNr++;
			}else {
				// recursive refinement
				radixStep(t_new, SA, bucketNr, BucketBegin, beg, end, h+1); 
			}
			beg = end + 1; // advance to next bucket
		}
  }
}

// Binary search for left boundry of interval.
long paraSA::bsearch_left(char c, long i, long s, long e) 
{
  if(c == S[SA[s]+i]) return s;
  long l = s, r = e;
  while (r - l > 1) 
	{
    long m = (l+r) / 2;
    if (c <= S[SA[m] + i]) r = m;
    else l = m;
  }
  return r;
}

// Binary search for right boundry of interval.
long paraSA::bsearch_right(char c, long i, long s, long e) 
{
  if(c == S[SA[e]+i]) return e;
  long l = s, r = e;
  while (r - l > 1) 
	{
    long m = (l+r) / 2;
    if (c < S[SA[m] + i]) r = m;
    else l = m;
  }
  return l;
}

// Simple top down traversal of a suffix array.
bool paraSA::top_down(char c, long i, long &start, long &end) 
{
  if(c < S[SA[start]+i]) return false;
  if(c > S[SA[end]+i]) return false;
  long l = bsearch_left(c, i, start, end);
  long l2 = bsearch_right(c, i, start, end);
  start = l; end = l2;
  return l <= l2;
}

// Top down traversal of the suffix array to match a pattern.  NOTE:
// NO childtab as in the enhanced suffix array (ESA).
bool paraSA::search(string &P, long &start, long &end) 
{
  start = 0; end = N - 1;
  long i = 0;
  while(i < (long)P.length()) 
	{
    if(top_down(P[i], i, start, end) == false) 
		{
      return false;
    }
    i++;
  }
  return true;

}

// Traverse pattern P starting from a given prefix and interval
// until mismatch or min_len characters reached.
void paraSA::traverse(string &P, long prefix, interval_t &cur, int min_len) 
{
  if(cur.depth >= min_len) return;
  while(prefix+cur.depth < (long)P.length()) 
	{
    long start = cur.start; 
		long end = cur.end;
    // If we reach a mismatch, stop.
    if(top_down_faster(P[prefix+cur.depth], cur.depth, start, end) == false) return;
    // Advance to next interval.
    cur.depth += 1; 
		cur.start = start; 
		cur.end = end;
    // If we reach min_len, stop.
    if(cur.depth == min_len) return;
  }

}

// Given SA interval apply binary search to match character c at
// position i in the search string. Adapted from the C++ source code
// for the wordSA implementation from the following paper: Ferragina
// and Fischer. Suffix Arrays on Words. CPM 2007.
bool paraSA::top_down_faster(char c, long i, long &start, long &end) 
{
  long l, r, m, r2=end, l2=start, vgl;
  bool found = false;
  long cmp_with_first = (long)c - (long)S[SA[start]+i];
  long cmp_with_last = (long)c - (long)S[SA[end]+i];
  
	if(cmp_with_first < 0) 
	{ 
    l = start+1; l2 = start; // pattern doesn't occur!  
  }else if(cmp_with_last > 0) 
	{ 
    l = end+1; l2 = end; 
    // pattern doesn't occur!  
  }else {

		// search for left border:
		l = start; r = end;
		if (cmp_with_first == 0) 
		{
			found = true; 
			r2 = r;
		}else{
			while (r - l > 1) 
			{
				m = (l+r) / 2;
				vgl = (long)c - (long)S[SA[m] + i];
				if (vgl <= 0) 
				{
					if (!found && vgl == 0) 
					{
						found = true;
						l2 = m; 
						r2 = r; // search interval for right border
					}
					r = m;
				}
				else l = m;
			}
			l = r;
		}
		
		// search for right border (in the range [l2:r2])
		if (!found) 
		{ 
			l2 = l - 1; // pattern not found => right border to the left of 'l' 
		}
		
		if (cmp_with_last == 0) 
		{ 
			l2 = end; // right border is the end of the array 
		}else 
		{
			while (r2 - l2 > 1) 
			{
				m = (l2 + r2) / 2;
				vgl = (long)c - (long)S[SA[m] + i];
				if (vgl < 0) r2 = m;
				else l2 = m;
			}
		}
  }
  start = l;
  end = l2;
  return l <= l2;
}

// Suffix link simulation using ISA/LCP heuristic.
bool paraSA::suffixlink(interval_t &m) 
{
  m.depth -= K;
  if( m.depth <= 0) return false;
  m.start = ISA[SA[m.start] / K + 1];  
  m.end = ISA[SA[m.end] / K + 1]; 
  return expand_link(m);

}

// For a given offset in the prefix k, find all MEMs.
void paraSA::findMEM(long k, 
										 string &P,
										 vector<match_t> &matches,
										 int min_len,
										 long id) 
{

	if(k < 0 || k >= K) { 
		cerr << "Invalid k." << endl; 
		return; 
	}

	long prefix = k; // Offset all intervals at different start points.
	interval_t mli(0,N/K-1,0); // min length interval
  interval_t xmi(0,N/K-1,0); // max match interval

  // Right-most match used to terminate search.
  int min_lenK = min_len - (K-1);

  while( prefix <= (long)P.length() - (K-k) )
	{
    traverse(P, prefix, mli, min_lenK); // Traverse until minimum length matched.
    if(mli.depth > xmi.depth) xmi = mli;
    if(mli.depth <= 1) { mli.reset(N/K-1); xmi.reset(N/K-1); prefix+=K; continue; }

		if(mli.depth >= min_lenK)
		{ 
			// Traverse until mismatch.
      traverse(P, prefix, xmi, P.length()); 
			// Using LCP info to find MEM length.
      collectMEMs(P, prefix, mli, xmi, matches, min_len, id); 
      // When using ISA/LCP trick, depth = depth - K. prefix += K. 
      prefix+=K;	
      if( suffixlink(mli) == false ) { mli.reset(N/K-1); xmi.reset(N/K-1); continue; }
      suffixlink(xmi);

		}else {
      prefix+=K;
      if( suffixlink(mli) == false ) { mli.reset(N/K-1); xmi.reset(N/K-1); continue; }
      xmi = mli;
    }
  }

}

// For a given offset in the prefix k, find all MEMs. (100%)
void paraSA::findMEMperfect(long k, 
														string &P,
														vector<match_t> &matches,
														int min_len,
														long id) 
{

	if(k < 0 || k >= K) { cerr << "Invalid k." << endl; return; }
  long prefix = k; // Offset all intervals at different start points.
	interval_t mli(0,N/K-1,0); // min length interval
  interval_t xmi(0,N/K-1,0); // max match interval
  // Right-most match used to terminate search.
  int min_lenK = min_len - (K-1);

	//100% part
	traverse(P, prefix, mli, min_lenK); // Traverse until minimum length matched.
	if(mli.depth > xmi.depth) xmi = mli;
	if(mli.depth >= min_lenK){ 
			// Traverse until mismatch.
      traverse(P, prefix, xmi, P.length()); 
			// Using LCP info to find MEM length.
      collectMEMsperfect(P, prefix, mli, xmi, matches, min_len, id); 
      // When using ISA/LCP trick, depth = depth - K. prefix += K. 
	}

}

// Use LCP information to locate right maximal matches. Test each for
// left maximality.
void paraSA::collectMEMs(string &P,
												 long prefix,
												 interval_t mli,
												 interval_t xmi,
												 vector<match_t> &matches,
												 int min_len,
												 long id)
{
	// All of the suffixes in xmi's interval are right maximal.
  for(long i = xmi.start; i <= xmi.end; i++) 
		find_Lmaximal(P, prefix, SA[i], xmi.depth, matches, min_len, id);
	if(mli.start == xmi.start && mli.end == xmi.end) return;
  while(xmi.depth >= mli.depth) 
	{
		// Attempt to "unmatch" xmi using LCP information.
		if(xmi.end+1 < N/K) xmi.depth = max(LCP[xmi.start], LCP[xmi.end+1]);
		else xmi.depth = LCP[xmi.start];

		// If unmatched XMI is > matched depth from mli, then examine rmems.
		if(xmi.depth >= mli.depth) 
		{
			// Scan RMEMs to the left, check their left maximality..
			while(LCP[xmi.start] >= xmi.depth) 
			{ 
				xmi.start--; 
				find_Lmaximal(P, prefix, SA[xmi.start], xmi.depth, matches, min_len, id);
			}
			// Find RMEMs to the right, check their left maximality.
			while(xmi.end+1 < N/K && LCP[xmi.end+1] >= xmi.depth) 
			{ 
				xmi.end++;
				find_Lmaximal(P, prefix, SA[xmi.end], xmi.depth, matches, min_len, id);
			}
		}
  }
}


// Use LCP information to locate right maximal matches. Test each for
// left maximality.
void paraSA::collectMEMsperfect(string &P, long prefix, interval_t mli, \
	interval_t xmi, vector<match_t> &matches, int min_len, long id) 
{
	//All of the suffixes in xmi's interval are right maximal.
  for(long i = xmi.start; i <= xmi.end; i++) 
		find_Lmaximal(P, prefix, SA[i], xmi.depth, matches, min_len, id);
}


// Finds left maximal matches given a right maximal match at position i.
void paraSA::find_Lmaximal(string &P, 
													 long prefix,
													 long i,
													 long len,
													 vector<match_t> &matches,
													 int min_len,
													 long id) 
{
		// Advance to the left up to K steps.
		for(long k = 0; k < K; k++) 
		{
			// If we reach the end and the match is long enough, print.
			if(prefix == 0 || i == 0) 
			{
				if(len >= min_len)
				{
					long refseq=0, refpos=0;
					from_set(i, refseq, refpos);
					if (descr[refseq]<id){
						matches.push_back(match_t(prefix, len, refseq, refpos));
					}
				}
				return; // Reached mismatch, done.
			}else if(P[prefix-1] != S[i-1]){
				// If we reached a mismatch, print the match if it is long enough.
				if(len >= min_len) 
				{
					long refseq=0, refpos=0;
					from_set(i, refseq, refpos);
					if (descr[refseq]<id){
						matches.push_back(match_t(prefix, len, refseq, refpos));
					}
				}
				return; // Reached mismatch, done.
			}
			prefix--; i--; len++; // Continue matching.

		}

}

// Load matching information.
long paraSA::load_match_info(long id, 
														 vector<match_t> &buf,
														 vector<mumi_unit> &mumiunits,
														 bool rc,
														 long qlen) 
{
	match_t m;
	mumi_unit mu;
	long refseq, refpos, bufs, addsize=0;
	bufs=(long)buf.size();

	// from_set is slow!!!
	for(long i=0; i<bufs; i++)
	{
		m=buf[i];
		refseq=m.refseq;
		refpos=m.refpos;
		addsize+=m.len;

		mu.ref=descr[refseq];

		mu.g0sens=1;
		mu.g0init=refpos + 1L;
		mu.g0long=m.len;
		mu.g0fin=refpos+m.len;

		if (rc){
			mu.g1sens=1;
			mu.g1init=m.query + 1L;
			mu.g1long=m.len;
			mu.g1fin=m.query+m.len;
		}else{
			mu.g1sens=2;
			mu.g1fin=qlen-m.query;
			mu.g1long=m.len;
			mu.g1init=mu.g1fin-m.len+1L;
		}

		mu.dgdelete=false; // G0 G1 filter
		mumiunits.push_back(mu); // Load one.
	}
	return(addsize);

}

// Finds maximal almost-unique matches (MAMs) 
bool paraSA::is_leftmaximal(string &P, long p1, long p2) 
{
  if(p1 == 0 || p2 == 0) return true;
  else return P[p1-1] != S[p2-1];
}

// Multithreads for large memory -- doing next
struct thread_data
{
  vector<long> Kvalues; // Values of K this thread should process.
  paraSA *sa; // Suffix array + aux informaton
  int min_len; // Minimum length of match.
  string *P; // Query string.
};

// Maximal Exact Matches (MEMs) 
void paraSA::MEM(string &P,
								 vector<match_t> &matches,
								 int min_len,
								 long id) 
{		
		for(int k = 0; k < K; k++) 
		{
			findMEM(k, P, matches, min_len, id); 
		}
}

// Maximal Exact Matches 100% (MEMs) 
void paraSA::MEMperfect(string &P, 
												vector<match_t> &matches, 
												int min_len, 
												long id) 
{		
		for(int k = 0; k < K; k++) 
		{
			findMEMperfect(k, P, matches, min_len, id);
		}

}
