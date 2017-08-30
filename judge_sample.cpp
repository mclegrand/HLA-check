
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <string>
#include <map>
#include <set>
#include <limits>
#include <unistd.h>
#include "omp.h"
#include <getopt.h>
using namespace std;

#define A 0
#define B 1
#define C 2
#define DRB1 3
#define DPA1 4
#define DPB1 5
#define DQA1 6
#define DQB1 7


string prefix;  
string alnfile = "info.txt";



// hla for ppl
vector<int> hla1_;
vector<int> hla2_;

//
map<string, string> nuc;
vector<std::string> hlalist;
vector<std::string> hladef;
string reftruc;
string refnuc;
string nucfilename;


vector<string>data_rs;//list of rs as in data file
vector<int>data_rs_pos;//list of rs pos
vector<vector<float> >data; //postprobs
vector<char> data_allele1; //list of allele1
vector<char> data_allele2; //list of allele2
vector<char> data_allele3; //list of allele3

int numberofnames;










void main_loop(float* , string& , string& , int ,const vector<string> &h, const vector<string> &, const vector<string> &, const vector<string> &, bool do_all = true);

inline char invert(char in) {
    if(in=='A')return 'T';
    if(in=='T')return 'A';
    if(in=='G')return 'C';
    if(in=='C')return 'G';
    return in;
}


int main(int argc, char** argv) {
    if(argc<4) {
        printf("%s A impute2_file hla_file [files/]\n", argv[0]);
        return 1;
    }
    if(argc<5) prefix="files/"; else prefix = argv[4];

    int chosen_hla = -1;
    if(!strcmp(argv[1],"A"))chosen_hla=A;
    if(!strcmp(argv[1],"B"))chosen_hla=B;
    if(!strcmp(argv[1],"C"))chosen_hla=C;
    if(!strcmp(argv[1],"DRB1"))chosen_hla=DRB1;
    if(!strcmp(argv[1],"DPA1"))chosen_hla=DPA1;
    if(!strcmp(argv[1],"DPB1"))chosen_hla=DPB1;
    if(!strcmp(argv[1],"DQA1"))chosen_hla=DQA1;
    if(!strcmp(argv[1],"DQB1"))chosen_hla=DQB1;
    if(chosen_hla == -1) {
        printf("HLA must be A,B,C,DRB1,DPA1,DPB1,DQA1, or DQB1\n");
        return 1;
    }

    ifstream aln_file;
    aln_file.open(prefix + alnfile);
    if (!aln_file.is_open()) {
        printf("Could not open alignment file at %s\n", (prefix+alnfile).c_str());
        return 1;
    }
    string hla, nm_name, coding_sequence, nm_sequence;
    while(hla != argv[1])
        aln_file >> hla >> nm_name >> reftruc >> nucfilename >> coding_sequence >> nm_sequence;
    ifstream nucfile;
    nucfile.open(prefix + nucfilename);
    if (!nucfile.is_open()) {
        printf("Could not open hla sequences definition at %s\n", (prefix + nucfilename).c_str());
        return 1;
    }

    ifstream rslist;
    rslist.open(prefix + "rs.txt");
    if (!rslist.is_open()) {
        printf("Could not open the rs list at %s\n", (prefix + "rs.txt").c_str());
        return 1;
    }

    ifstream datafile;
    datafile.open(argv[2]);
    if (!datafile.is_open()) {
        printf("Could not open data file %s\n", argv[2]);
        return 1;
    }

    ifstream hlafile;
    hlafile.open(argv[3]);
    if (!hlafile.is_open()) {
        printf("Could not open hla file %s\n",argv[3]);
        return 1;
    }


    {
        string line, a;
        int x,y;
        while(getline(hlafile,line)) {
            stringstream ss(line);
            string a,b,c,d,e,f;
            ss >> a >> b >> c >> d >> e >> f;
            int order[] = {0,2,1,3,6,7,4,5};//order in files is A, C, B, DR, DQ, DP (chromosomal order)
            for(int i=-1; i<order[chosen_hla]; i++)
                ss >> x >> y;
            if(x<100)x*=100;
            if(y<100)y*=100;
            hla1_.push_back(x);
            hla2_.push_back(y);
        }
    }


    string line;
    string a,b;
    while ( nucfile >> a >> b ) {
        int pos=a.find_first_of('*');
        a=a.erase(0,pos+1);
        hlalist.push_back(a);
        hladef.push_back(b);
        nuc[a]=b;
    }
    refnuc = nuc[reftruc];
    string x[13];
    map<string, int> rslist_nuc;
    while ( 1 ) {
        for(int i=0; i<=12; i++) {
            if(!(rslist >> a))goto endlecture;
            x[i]=a;
        }
        if(nm_name!=x[4])continue;
        int rspos_nm = stoi(x[7]);
        int rspos_f =0;
        for(int i=0;; i++) {
            if(nm_sequence[i]=='A' || nm_sequence[i]=='T' || nm_sequence[i]=='G' || nm_sequence[i]=='C') {
                rspos_nm--;
                if(rspos_nm<0) {
                    rspos_f=i-1;
                    break;
                }
            }
        }
        int rspos_nuc =0;
        for(int i=0; i<rspos_f; i++) {
            if(coding_sequence[i]=='A' || coding_sequence[i]=='T' || coding_sequence[i]=='G' || coding_sequence[i]=='C') {
                rspos_nuc++;
            }
        }
        int refpos_nuc2=0;
        for(int i=0;; i++) {
            if(refnuc[i]=='A' || refnuc[i]=='T' || refnuc[i]=='G' || refnuc[i]=='C' ) {
                rspos_nuc--;
                if(rspos_nuc<0) {
                    refpos_nuc2=i;
                    break;
                }
            }
        }
        rslist_nuc[x[0]]=refpos_nuc2;
    }
endlecture:



/////reading data file
    vector<string> oldlines;

    string cur_rs="";
    bool end=false;
    while(!end) {
        string rs,allele1,allele2;
        if(!getline(datafile,line)) {
            end=true;
            rs = "rs0";
        }
        else {
            stringstream x(line);
            x>>a;
            x>>rs;
            x>>a;
            x>>allele1;
            x>>allele2;
        }
        if(rs != cur_rs) {
            stringstream x,y;
            char a1,a2,a3,a4;
            vector<float> tmp;
            float aa,bb,cc,dd,ee,ff,f;
            //do stuff with oldlines
            switch(oldlines.size()) {
            case 1://biallelic
                x.str(oldlines[0]);
                x>>a;
                x>>cur_rs;
                x>>a;
                x>>a1;
                x>>a2;
                data_allele1.push_back(a1);
                data_allele2.push_back(a2);
                data_allele3.push_back('-');
                data_rs.push_back(cur_rs);
                data_rs_pos.push_back(rslist_nuc[cur_rs]);
                tmp.clear();
                while(x>>f)tmp.push_back(f);
                data.push_back(tmp);
                break;
            case 2://triallelic
                x.str(oldlines[0]);
                y.str(oldlines[1]);
                x>>a;
                x>>cur_rs;
                //std::cout << oldlines[0] << std::endl << oldlines[1] << std::endl;
                x>>a;
                x>>a1;
                x>>a2;
                y>>a;
                y>>cur_rs;
                y>>a;
                y>>a3;
                y>>a4;
                if(a1 != a3 || a2==a4) {
                    std::cerr << "wrong triallelic encoding at " << cur_rs.c_str() << std::endl;
                    oldlines.clear();
                    continue;
                }
                data_allele1.push_back(a1);
                data_allele2.push_back(a2);
                data_allele3.push_back(a4);
                data_rs.push_back(cur_rs);
                data_rs_pos.push_back(rslist_nuc[cur_rs]);
                tmp.clear();
                while(x>>aa>>bb>>cc) {
                    y>>dd>>ee>>ff;
                    f = aa*dd+bb*dd+dd*cc+aa*ee+bb*ee+ff*aa;
                    //AA AB BB AC BC CC
                    tmp.push_back(aa*dd/f);
                    tmp.push_back(bb*dd/f);
                    tmp.push_back(cc*dd/f);
                    tmp.push_back(aa*ee/f);
                    tmp.push_back(bb*ee/f);
                    tmp.push_back(ff*aa/f);
                }
                data.push_back(tmp);

            }

            oldlines.clear();
        }

        if(allele1.size()!=1 || allele2.size()!=1)continue;
        if(rs == "6")continue;
        if(rslist_nuc.end()==rslist_nuc.find(rs))continue;
        cur_rs = rs;
        oldlines.push_back(line);
    }
    std::cerr << "#SNPs used:" << data.size() << std::endl;



///////end reading data file

//for hla B and C, invert snps because of strand
    if(chosen_hla==B || chosen_hla==C || chosen_hla==DRB1 || chosen_hla==DQB1 || chosen_hla==DPA1 ) {
        for(int i=0; i<data_allele1.size(); i++) {
            data_allele1[i]=invert(data_allele1[i]);
            data_allele2[i]=invert(data_allele2[i]);
            data_allele3[i]=invert(data_allele3[i]);
        }
    }

    if(data.empty()) {
        printf("No SNP found in your data for selected HLA gene: data format failure ?\n");
        return 1;
    }
    numberofnames = data[0].size()/3;

//openmp calls
    omp_set_dynamic(0);
    omp_set_num_threads(8);

    #pragma omp parallel for 
    for(int i=0; i<numberofnames; i++) {
        float scoremin = -numeric_limits<float>::infinity() ;
        string hla1min="";
        string hla2min = "";

        main_loop(&scoremin, hla1min, hla2min, i,hlalist ,hlalist ,hladef,hladef,false);


        char pr[100];
        sprintf(pr, "%d %f %s %s ",i+1, scoremin,hla1min.c_str(),hla2min.c_str());


        vector<string> compat1;
        vector<string> compat2;
        vector<string> hladef1;
        vector<string> hladef2;

        for(int ii=0; ii<hlalist.size(); ii++) {
            //if(hlalist[ii][5]>=0 && hlalist[ii][5]<=9)continue;
            istringstream iss(hlalist[ii]);
            int h2,h4;
            char xxxxx;
            iss >> h2 >> xxxxx >> h4; //xxxxx is ':', but whatever.
            if(h4>100)continue;//FIXME someday
            int h = 100*h2+h4;
            //int h = stoi(hlalist[ii].substr(0,2))*100+stoi(hlalist[ii].substr(3,2));
            if( h == hla1_[i] ||  100*(h/100) == hla1_[i] ) {
                compat1.push_back(hlalist[ii]);
                hladef1.push_back(hladef[ii]);
            }
            if( h == hla2_[i] ||  100*(h/100) == hla2_[i] ) {
                compat2.push_back(hlalist[ii]);
                hladef2.push_back(hladef[ii]);
            }
        }

        float realmin = scoremin;
        scoremin = -numeric_limits<float>::infinity();
        main_loop(&scoremin, hla1min, hla2min, i,compat1 ,compat2 ,hladef1,hladef2);

        #pragma omp critical
        printf("%s %d %d %f\n",pr, hla1_[i], hla2_[i],scoremin-realmin);

    }

    return 0;
}//main





void main_loop(float* scoremin, string& hla1min, string& hla2min, int i,const vector<string> &hlax1, const vector<string> & hlax2, const vector<string> &defx1, const vector<string> &defx2, bool do_all ) {


    for(int ii=0; ii<hlax1.size(); ii++) {
        string hla1 = hlax1[ii];
        string seq1 = defx1[ii];
        int jj=do_all?0:ii;
        for(; jj<hlax2.size(); jj++) {
            string hla2 = hlax2[jj];
            string seq2 = defx2[jj];
            float scoretot = 0;

            for(int k=0; k<data.size(); k++) {
                if(scoretot<*scoremin)break;
                float score = 0;
                string rs = data_rs[k];
                char allele1 = data_allele1[k];
                char allele2 = data_allele2[k];
                char allele3 = data_allele3[k];
                int pos = data_rs_pos[k];
                char h1 = seq1[pos];
                if(h1 == '-') {
                    h1=refnuc[pos];
                }
                char h2 = seq2[pos];
                if(h2 == '-') {
                    h2=refnuc[pos];
                }
                if(h2 == '*' && h1=='*') {
                    continue;
                }
                if(allele3 == '-') {
                    if((h1!= allele1 && h1 != allele2 && h1 != '*')||(h2!= allele1 && h2 != allele2 && h2 != '*') ) {
                        score=1;
                        goto scoring;
                    }
                    if(h1 == '*') {
                        if(h2 == allele1) score += (data[k][3*i+2]);//*A
                        else if(h2 == allele2) score += (data[k][3*i  ]);//*B
                        else printf("wut1 %s\n", rs.c_str());
                    } else if(h2 == '*') {
                        if(h1 == allele1) score += (data[k][3*i+2]);//A*
                        else if(h1 == allele2) score += (data[k][3*i  ]);//B*
                        else printf("wut2 %s\n",rs.c_str());
                    } else { // pas de *
                        if((h1==allele1 && h2 == allele2) || (h1==allele2 && h2 == allele1) ) score += (data[k][3*i  ]+data[k][3*i+2]);//AB ou BA
                        else if(h1 == allele1 && h2 == allele1)                               score += (data[k][3*i+1]+data[k][3*i+2]);//AA
                        else if(h1 == allele2 && h2 == allele2)                               score += (data[k][3*i  ]+data[k][3*i+1]);//BB
                        else printf("wut3 %c %c %c %c %s \n",h1,h2,allele1,allele2,rs.c_str());
                    }
                } else {
                    if((h1!= allele1 && h1 != allele2 && h1 != allele3 && h1 != '*')||(h2!= allele1 && h2 != allele2 && h2 != allele3 && h2 != '*') ) {
                        score=1;
                        goto scoring;
                    }
                    if(h1 == '*') {
                        if(h2 == allele1) score +=      (data[k][6*i+2]+data[k][6*i+4]+data[k][6*i+5]);//*A
                        else if(h2 == allele2) score += (data[k][6*i  ]+data[k][6*i+3]+data[k][6*i+5]);//*B
                        else if(h2 == allele3) score += (data[k][6*i  ]+data[k][6*i+1]+data[k][6*i+2]);//*C
                        else continue;//printf("wut4\n");
                    } else if ( h2 == '*') {
                        if(h1 == allele1) score +=      (data[k][6*i+2]+data[k][6*i+4]+data[k][6*i+5]);//A*
                        else if(h1 == allele2) score += (data[k][6*i  ]+data[k][6*i+3]+data[k][6*i+5]);//B*
                        else if(h1 == allele3) score += (data[k][6*i  ]+data[k][6*i+1]+data[k][6*i+2]);//C*
                        else continue;//printf("wut5\n");
                    } else {
                        if((h1==allele1 && h2 == allele2) || (h1==allele2 && h2 == allele1) )      score+= (1-data[k][6*i+1]); //AB ou BA
                        else if((h1==allele3 && h2 == allele2) || (h1==allele2 && h2 == allele3) ) score+= (1-data[k][6*i+4]); //CB ou BC
                        else if((h1==allele1 && h2 == allele3) || (h1==allele3 && h2 == allele1) ) score+= (1-data[k][6*i+3]); //AC ou CA
                        else if(h1 == allele1 && h2 == allele1) score += (1-data[k][6*i  ]);//AA
                        else if(h1 == allele2 && h2 == allele2) score += (1-data[k][6*i+2]);//BB
                        else if(h1 == allele3 && h2 == allele3) score += (1-data[k][6*i+5]);//CC
                        else continue;// printf("wut6\n");
                    }
                }
scoring:
                if(score>1.001 || score<0) printf("wat %f %s\n",score,data_rs[k].c_str());

                //scoretot += (0.99-(1/(score+0.01)));nope
                scoretot+=(-score);
                //scoretot+=log((1-score)*0.9+.1);
                //scoretot+=log(1.01-score);

                //la fonction doit être NÉGATIVE et DÉCROISSANTE SUR [0,1] :
                //1 = "pas bon du tout"; 0 = "perfect match"
            }
            if(scoretot > *scoremin) {
                (*scoremin)=scoretot;
                hla1min=hla1;
                hla2min=hla2;
            }
        }
    }
}//main_loop


