#include "para_sortgenome.h"

ParallelSort::ParallelSort(){
        MPI_Comm_rank(MPI_COMM_WORLD, &myid_);
        MPI_Comm_size(MPI_COMM_WORLD, &pnum_);
}

void ParallelSort::Read(const string& filename){
	ifstream fin1(filename.c_str()); //.fna
	string line;
	int i = -1;
	string tit="";
	string seq="";
	Genome gmo;
	while (getline(fin1, line)) {  // Each line is a genome
		if(line[0]=='>'){ // title line
			i++;	
			if((i > 0) && ((i-1) % pnum_ == myid_)){
				gmo.title = tit;
				gmo.sequence = seq;
				genomes_.push_back(gmo);
				tit="";
				seq="";
				
			}
		}
		if(i % pnum_ == myid_){
			if(line[0]=='>'){
				tit=line;
			}else{
				seq=seq+'\n';
				seq+=line;	
			}
		}
        }
	if(i % pnum_ == myid_){
		gmo.title = tit;
		gmo.sequence = seq;
		genomes_.push_back(gmo);
	}
  num_total_docs_ = i+1;
  if(myid_==0){ cerr << "Total genomes: " << num_total_docs_ << endl; } 
//if(myid_ ==0 ){ for(int i=0; i< 2; i++){ cerr<<"*****"<< genomes_[i].title<<genomes_[i].sequence << endl;  } } //test line
}
void ParallelSort::Write(const string& filename1, const string& filename2){ // filename1:output file, filename2: index 
        ifstream fin2(filename2.c_str()); //sorted index value
	string line;
        vector<int>ids;
        while (getline(fin2, line)){
                istringstream ss(line);
                int index;
                int value;
                ss >> index >> value;
                ids.push_back(index);
        }
	//if(myid_ == 0 ){ cerr << ids[0] << endl;} // test line
	for(int i = 0; i < ids.size(); ++i)
	{
		std::ofstream fout;
		if(ids[i] % pnum_ == myid_){ 
			if (i == 0) {
                            fout.open(filename1.c_str());
                        } else {
                            fout.open(filename1.c_str(), std::ios::app);
	                }
			int ii=ids[i]/pnum_;
			if(ii<genomes_.size()){ 
				fout << genomes_[ii].title; 
				fout << genomes_[ii].sequence << endl;
			}else{
				cerr << "boundary error in Write()" << endl;
			}
			//if(myid_ ==0 && i<1){ cerr <<"ids:" << i << " " << ids[i] << ", ids[i]/pnum_:" <<  ids[i]/pnum_ <<endl << genomes_[ids[i]/pnum_].title<< endl;} // test line
      			fout.close();
		}
	MPI_Barrier(MPI_COMM_WORLD);
	}
}

using namespace std;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  string FLAGS_input = "";
  string FLAGS_index = "";
  string FLAGS_output = "";
  for (int i = 1; i < argc; ++i) {
    if (0 == strcmp(argv[i], "--input")) {
      FLAGS_input = argv[i+1];
      ++i;
    } else if (0 == strcmp(argv[i], "--index")){
      FLAGS_index = argv[i+1];
       ++i;
    } else if (0 == strcmp(argv[i], "--output")) {
      FLAGS_output = argv[i+1];
      ++i;
    }
  }
  if (FLAGS_input == "") {
    cerr << "--input must not be empty" << endl;
    MPI_Finalize();
    return 1;
  }
  if (FLAGS_index == "") {
    cerr << "--index must not be empty" << endl;
    MPI_Finalize();
    return 1;
  }
  if (FLAGS_output == "") {
    cerr << "--output must not be empty" << endl;
    MPI_Finalize();
    return 1;
  }
//cerr << FLAGS_input << endl << FLAGS_index << endl << FLAGS_output << endl;
  ParallelSort parasort;
  parasort.Read(FLAGS_input);
  parasort.Write(FLAGS_output, FLAGS_index);
  MPI_Finalize();
  return 0;
}


