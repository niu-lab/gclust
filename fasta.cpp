#include <fstream>
#include <iostream>
#include <algorithm>
#include "fasta.hpp"

// Filter 'n' in genome.
void filter_n(string &seq_rc) 
{
	// Bit by bit.
  for(long i = 0; i < (long)seq_rc.length(); i++) {
    // Adapted from Kurtz code in MUMmer v3. 
    switch(seq_rc[i]) 
		{
			case 'a': case 't': case 'g': case 'c': break;
			default:
				seq_rc[i] = '~';
		}
  }
}

// Return the reverse complement of sequence. This allows searching
// the plus strand of instances on the minus strand.
void reverse_complement(string &seq_rc, bool nucleotides_only) 
{
  // Reverse in-place.
  reverse(seq_rc.begin(), seq_rc.end());
  for(long i = 0; i < (long)seq_rc.length(); i++) 
	{
    // Adapted from Kurtz code in MUMmer v3. 
    switch(seq_rc[i]) 
		{
    case 'a': seq_rc[i] = 't'; break;
    case 'c': seq_rc[i] = 'g'; break;
    case 'g': seq_rc[i] = 'c'; break;
    case 't': seq_rc[i] = 'a'; break;
    case 'r': seq_rc[i] = 'y'; break; /* a or g */
    case 'y': seq_rc[i] = 'r'; break; /* c or t */
    case 's': seq_rc[i] = 's'; break; /* c or g */
    case 'w': seq_rc[i] = 'w'; break; /* a or t */
    case 'm': seq_rc[i] = 'k'; break; /* a or c */
    case 'k': seq_rc[i] = 'm'; break; /* g or t */
    case 'b': seq_rc[i] = 'v'; break; /* c, g or t */
    case 'd': seq_rc[i] = 'h'; break; /* a, g or t */
    case 'h': seq_rc[i] = 'd'; break; /* a, c or t */
    case 'v': seq_rc[i] = 'b'; break; /* a, c or g */
    default:  
      if(!nucleotides_only) seq_rc[i] = 'n'; 
      break; /* anything */
    }
  }

}

// trim a string.
void trim(string &line, long &start, long &end) 
{
  // Trim leading spaces. 
  for(long i = start; i < (int)line.length(); i++) 
	{ 
    if(line[i] != ' ') 
			{ 
				start = i; 
				break;
			} 
  }
  // Trim trailing spaces.
  for(long i = line.length() - 1; i >= 0; i--) { 
    if(line[i] != ' ') { 
			end = i;
			break;
		} 
    if(i == 0) break;
  }
}

// load data.
void load_fasta(string filename, 
								string &S, 
								vector<string> &descr,
								vector<long> &startpos) 
{
  string meta, line;
  long length = 0;

  // Everything starts at zero.
  startpos.push_back(0);
  ifstream data(filename.c_str());
  if(!data.is_open()) { 
		cerr << "unable to open " << filename << endl; 
		exit(1); 
	} 
  while(!data.eof()) 
	{
    getline(data, line); // Load one line at a time.
    if(line.length() == 0) continue;
    long start = 0, end = line.length() - 1;
    // Meta tag line and start of a new sequence.
    if(line[0] == '>') {
      // Save previous sequence and meta data.
      if(length > 0) {
				descr.push_back(meta);
				//cout<<meta<<endl;
				S += '`'; // ` character used to separate strings
				startpos.push_back(S.length());
				//lengths.push_back(length+1);
      }
      // Reset parser state.
      start = 1; meta = ""; length = 0;
    }
    trim(line, start, end);
    // Collect meta data.
    if(line[0] == '>') {
      for(long i = start; i <= end; i++) { 
				if(line[i] == ' ') break; 
				meta += line[i]; 
			}
    }else { // Collect sequence data.
      length += end - start + 1;
      for(long i = start; i <= end; i++) { 
				S += std::tolower(line[i]);
      }
    }
  }
  if(length > 0) {
    descr.push_back(meta);
		//cout<<meta<<endl;
  }  
  cerr << "# S.length=" << S.length() << endl;
  for(long i = 0; i < (long)descr.size(); i++) {
    cerr << "# " << descr[i] << " " << startpos[i] << endl;
  }

}

// Load part genomes of total for clustering.
// Previous genomes have been processed need to skip.
void load_part_genomes(string filename, 
											 vector<Genome> &partgenomes,
											 vector<GenomeClustInfo> &totalgenomes,
											 long previous, 
											 long number)
{
	string meta, line, S;
  long length = 0;
	Genome tg;
  ifstream data(filename.c_str());
  if(!data.is_open()) { 
		cerr << "unable to open " << filename << endl; 
		exit(1); 
	}
	// Skip previous genomes.
	long i=0;
	while(!data.eof()) {
		getline(data, line);
		if(line[0] == '>') {
			i++;
			if (i>previous){
				break;
			}
		}
	}
	// Genome id from previous.
	long id=i-1; 
	long loadnum=0;
	// First one sequence.
	meta = "";
	long start = 1, end = line.length() - 1;
	for(i = start; i <= end; i++) { 
		if(line[i] == ' ') break; 
		meta += line[i]; 
	}
	S="";
  while(!data.eof()) {
    getline(data, line); // Load one line at a time.
    if(line.length() == 0) continue;
    long start = 0, end = line.length() - 1;
    // Meta tag line and start of a new sequence.
    if(line[0] == '>') {
			// Save previous sequence and meta data.
			if(length > 0) {
				tg.descript=meta;
				tg.id=id++;
				tg.size=S.length();
				tg.cont=S;
				partgenomes.push_back(tg);
			  loadnum++;
				S="";
				// Check genomes loaded one time.
				if (loadnum>=number){
					break;
				}
      }
      // Reset parser state.
      start = 1; meta = ""; length = 0;
    }
    trim(line, start, end);
    // Collect meta data.
    if(line[0] == '>') {
      for(long i = start; i <= end; i++) { 
				if(line[i] == ' ') break; 
				meta += line[i]; 
			}
    }else { // Collect sequence data.
      length += end - start + 1;
      for(long i = start; i <= end; i++) { 
				S += std::tolower(line[i]);
      }
    }
  }
	// last sequence
	if (loadnum<number){
		if(length > 0) {
			tg.descript=meta;
			tg.id=id++;
			tg.size=S.length();
			tg.cont=S;
			partgenomes.push_back(tg);
			loadnum++;
		}  
	}

}

// Load genomes from memory.
void load_part_genomes_mem(vector<Genome> &allpartgenomes,
													 vector<Genome> &partgenomes,
													 vector<GenomeClustInfo> &totalgenomes,
													 long previous, 
													 long number)
{	
	Genome tg;
	long loadnumber=0;
	long genomes=allpartgenomes.size();
	for (long i=previous; i<genomes; i++){
		tg=allpartgenomes[i];
		partgenomes.push_back(tg);
		loadnumber++;
		if (loadnumber>=number){
			break;
		}
	}

}

// Load genomes from file.
void load_part_genomes_internal(string filename,
																vector<Genome> &partgenomes,
																vector<GenomeClustInfo> &totalgenomes,
																long previous,
																int &number, 
																long totalsize, 
																bool &ifend,
																int memiden)
{
	string meta, line, S;
	Genome tg;
	long length, loadgenomes, id, loadnumber, start, end, sizeadd, stablenumber;
	length=loadgenomes=loadnumber=sizeadd=0;
	if (memiden==100)
	{
		stablenumber=MAX_PARTNUMBERFORPERFECT;
	}else{
		stablenumber=MAX_PARTNUMBER;
	}
	bool beforeend=false;
	ifend=false;
  ifstream data(filename.c_str());
  if(!data.is_open()) { 
		cerr << "unable to open " << filename << endl; 
		exit(1); 
	}
	// Skip previous genomes.
	long i=0;
	while(!data.eof()) {
		getline(data, line);
		if(line[0] == '>') {
			i++;
			if (i>previous){ break;}
		}
	}
	id=i-1; // Genome id from previous.

	// First sequence.
	meta = "";
	start = 1; end = line.length() - 1;
	for(i = start; i <= end; i++) { 
		if(line[i] == ' ') break; 
		meta += line[i]; 
	}
	S="";
	while(!data.eof()) 
	{
    getline(data, line); // Load one line at a time.
    if(line.length() == 0) continue;
    start = 0; end = line.length() - 1;
    // Meta tag line and start of a new sequence.
    if(line[0] == '>') {
			// Save previous sequence and meta data.
			if(length > 0) {
				tg.descript=meta;
				tg.id=id++;
				tg.size=S.length();
				tg.cont=S;
				if (totalgenomes[tg.id].rep){
					sizeadd+=tg.size;
					if ((sizeadd<=totalsize)&&(loadgenomes<stablenumber)){
						partgenomes.push_back(tg);
						loadnumber++;
						loadgenomes++;
					}else{
						beforeend=true;
						break;
					}
				}else{
					partgenomes.push_back(tg);
					loadnumber++;
				}
				S="";
			}
      // Reset parser state.
      start = 1; meta = ""; length = 0;
    }
    trim(line, start, end);
    if(line[0] == '>') {
      for(long i = start; i <= end; i++) { 
				if(line[i] == ' ') break; 
				meta += line[i]; 
			}
    }else { // Collect sequence data.
      length += end - start + 1;
      for(long i = start; i <= end; i++) { 
				S += std::tolower(line[i]);
      }
    }
  }
	// The last one
	if(length > 0) {
			tg.descript=meta;
			tg.id=id++;
			tg.size=S.length();
			tg.cont=S;
	}
	if (!beforeend){
		
		partgenomes.push_back(tg);
		loadnumber++;
		if (totalgenomes[tg.id].rep){
			loadgenomes++;
		}
		ifend=true;
	}
	number=loadnumber;

}

// Load genomes from memory.
void load_part_genomes_internal_mem(vector<Genome> &allpartgenomes,
																		vector<Genome> &partgenomes,
																		vector<GenomeClustInfo> &totalgenomes,
																		long previous,
																		int &number,
																		long totalsize,
																		bool &ifend,
																		int memiden)
{
	Genome tg;
	long loadnumber, loadgenomes, genomes, sizeadd, stablenumber, i;
	loadnumber=loadgenomes=sizeadd=0;
	genomes=allpartgenomes.size();
	if (memiden==100)
	{
		stablenumber=MAX_PARTNUMBERFORPERFECT;
	}else{
		stablenumber=MAX_PARTNUMBER;
	}
	ifend=false;
	for (i=previous; i<genomes; i++){
		tg=allpartgenomes[i];
		if (totalgenomes[tg.id].rep){
			sizeadd+=tg.size;
			if ((sizeadd<=totalsize)&&(loadgenomes<stablenumber)){
				partgenomes.push_back(tg);
				loadnumber++;
				loadgenomes++;
			}else{
				break;
			}
		}else{
			partgenomes.push_back(tg);
			loadnumber++;
		}
	}
	number=loadnumber;
	if (i==genomes){
		ifend=true;
	}

}

// Load all genomes same as loading part genomes
void load_part_genomes_all(string filename, vector<Genome> &partgenomes)
{
	string meta, line, S;
  long length = 0;
	Genome tg;
  ifstream data(filename.c_str());
  if(!data.is_open()) { 
		cerr << "unable to open " << filename << endl; 
		exit(1); 
	}
	// Skip previous genomes.
	long i=0;
	while(!data.eof()) {
		getline(data, line);
		if(line[0] == '>') {
			i++;
			if (i>0){
				break;
			}
		}
	}
	// Genome id from previous.
	long id=0;
	long loadnum=0;
	// First one sequence.
	meta = "";
	long start = 1, end = line.length() - 1;
	for(i = start; i <= end; i++) { 
		if(line[i] == ' ') break; 
		meta += line[i]; 
	}
	S="";
  while(!data.eof()) {
    getline(data, line); // Load one line at a time.
    if(line.length() == 0) continue;
    long start = 0, end = line.length() - 1;
    // Meta tag line and start of a new sequence.
    if(line[0] == '>') {
			// Save previous sequence and meta data.
			if(length > 0) {
				tg.descript=meta;
				tg.id=id++;
				tg.size=S.length();
				tg.cont=S;
				partgenomes.push_back(tg);
			  loadnum++;
				S="";
      }
      // Reset parser state.
      start = 1; meta = ""; length = 0;
    }
    trim(line, start, end);
    // Collect meta data.
    if(line[0] == '>') {
      for(long i = start; i <= end; i++) { 
				if(line[i] == ' ') break; 
				meta += line[i]; 
			}
    }else { // Collect sequence data.
      length += end - start + 1;
      for(long i = start; i <= end; i++) { 
				S += std::tolower(line[i]);
      }
    }
  }
	// last sequence
	if(length > 0) 
	{
		tg.descript=meta;
		tg.id=id++;
		tg.size=S.length();
		tg.cont=S;
		partgenomes.push_back(tg);
		loadnum++;
	}

}

//Test part genomes loading vector 
void test_part(vector<Genome> &partgenomes)
{
	Genome tg;
	int s=partgenomes.size();
	for (int i=0;i<s;i++){	
		tg=partgenomes[i];
		cout<<"===========\n";
		cout<<tg.descript<<endl;
		cout<<tg.id<<endl;
		cout<<tg.size<<endl;
		cout<<tg.cont<<endl;
		cout<<"===========\n";
	}
}

// Make one suffix array from one block.
void make_block_ref(vector<Genome> &partgenomes, string &S,
										vector<GenomeClustInfo> &totalgenomes,
										vector<long> &descr, 
										vector<long> &startpos)
{
	long pos = 0;
	Genome tg;
	long s=partgenomes.size();
	S="";
	startpos.push_back(0);
	for (long i=0;i<s;i++){	
		tg=partgenomes[i];
		if (!totalgenomes[tg.id].rep){
			continue;
		}
		S += tg.cont;
		S += '`';
		pos = pos+tg.size+1;
		startpos.push_back(pos);
		descr.push_back(tg.id);
	}
	startpos.pop_back();
	int k = S.length(); 
	S = S.substr(0,k-1);
	cerr<<"\n===="<<endl;
	cerr<<"S "<<S.length()<<endl;
	cerr<<"startpos "<<startpos.size()<<endl;
	cerr<<"descr "<<descr.size()<<endl;
	cerr<<"====\n"<<endl;

}

// Load total genomes one time.
void load_total_genomes(string filename, 
												vector<GenomeClustInfo> &totalgenomes)
{
	long length, maxlen, totallen, minlen;
	length=maxlen=totallen=0;
	minlen=MAX_GENOME;
	int id, loadnum; // Genome id from previous.
	id=loadnum=0;
	string meta, line;
	GenomeClustInfo tg;
  // Everything starts at zero.
  ifstream data(filename.c_str());
  if(!data.is_open()) { 
		cerr << "unable to open " << filename << endl; 
		exit(1); 
	}
  while(!data.eof()) {
    getline(data, line); // Load one line at a time.
    if(line.length() == 0) continue;
    long start = 0, end = line.length() - 1;
    if(line[0] == '>') { // Meta tag line and start of a new sequence.
			if(length > 0) { // Save previous sequence and meta data.
				tg.descript=meta;
				tg.id=id;
				//tg.pid=id;
				tg.rep=true;
				tg.size=length;
				totallen+=length;
				totalgenomes.push_back(tg);
				if (length>maxlen){	maxlen=length; }
				if (length<minlen){ minlen=length; }
				id++;
				loadnum++;
      }
      start = 1; meta = ""; length = 0; // Reset parser state.
    }
    trim(line, start, end);
    if(line[0] == '>') { // Collect meta data.
      for(long i = start; i <= end; i++) { 
				if(line[i] == ' ') break; 
				meta += line[i];
			}
    }else { // Collect sequence data.
      length += end - start + 1;
    }
  }
	
	if(length > 0) { // Last one.
		tg.descript=meta;
		tg.id=id;
		tg.rep=true;
		tg.size=length;
		totalgenomes.push_back(tg);
		totallen+=length;
		if (length>maxlen){	maxlen=length; }
		if (length<minlen){ minlen=length; }
		loadnum++;
		id++;
	}

	cerr<<"=====================\n\n";
	cerr<<"Genomes: "<<loadnum<<endl;
	cerr<<"Maximum Length: "<<maxlen<<endl;
	cerr<<"Minimum Length: "<<minlen<<endl;
	cerr<<"Average Length: "<<totallen/loadnum<<endl;
	cerr<<"\n";
	cerr<<"=====================\n\n";

}

