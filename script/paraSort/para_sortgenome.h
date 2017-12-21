#include<iostream>
#include<fstream>
#include<sstream>
#include<list>
#include<map>
#include<string>
#include<vector>
#include<cmath>
#include<ctime>
#include<stdlib.h>
#include<mpi.h>
#include<list>
#include<algorithm>
#include<functional>
#include<iterator>
#include<string.h>
using namespace std;

struct Genome{
//        int index;
	string title;
        string sequence;
//	int len;
};

/*struct IndexValue {
  int index;
  //int local_index;
  double value;
  //IndexValue() {}
  //IndexValue(int i, double v) : index(i), value(v) {
 // }
 // 
  //bool operator > (const IndexValue& Piv) const; //decresing sort
};*/

/*bool IndexValue::opeartor>(const IndexValue& Piv) const
{
	return value>Piv.value;
}

struct CmpIndexValueByValueDsc {
  bool operator()(const IndexValue& p1,
                  const IndexValue& p2) {
    return p1.value > p2.value;
  }
};
*/
class ParallelSort{
public:
	ParallelSort();
	void Read(const string& filename);
	void Write(const string& filename1, const string& filename2);
private:
	vector<Genome> genomes_;
	int num_total_docs_;
	int myid_;
	int pnum_;
	//vector<IndexValue>iv_;
};

