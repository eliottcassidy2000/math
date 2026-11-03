#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using Mat=std::array<long long,9>;
std::vector<Mat> bank;
std::vector<std::string> names;
std::vector<int> signs;
uint64_t tested=0, trace_candidates=0;
int hits=0;
Mat mul(const Mat&a,const Mat&b){
    Mat c{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)c[3*i+j]+=a[3*i+k]*b[3*k+j];
    return c;
}
void visit(const Mat&m,int depth,int det,std::string word){
    if(depth){
        ++tested;
        long long tr=m[0]+m[4]+m[8];
        int q=0,want_second=0,want_det=0;
        if(depth<=2&&tr==5){q=3;want_second=3;want_det=-1;}
        if(depth<=4&&tr==19){q=6;want_second=19;want_det=1;}
        if(depth<=7&&tr==77){q=9;want_second=75;want_det=-1;}
        if(q){
            ++trace_candidates;
            __int128 second=(__int128)m[0]*m[4]+(__int128)m[0]*m[8]+(__int128)m[4]*m[8]
                            -(__int128)m[1]*m[3]-(__int128)m[2]*m[6]-(__int128)m[5]*m[7];
            if(second==want_second&&det==want_det){++hits;std::cout<<"HIT q="<<q<<" word="<<word<<"\n";}
        }
    }
    if(depth==7)return;
    for(size_t i=0;i<bank.size();++i)visit(mul(m,bank[i]),depth+1,det*signs[i],word+" "+names[i]);
}
int main(int argc,char**argv){
    if(argc!=2)return 2;
    std::ifstream in(argv[1]);std::string name;int sign;
    while(in>>name>>sign){Mat m{};for(auto &x:m)in>>x;names.push_back(name);signs.push_back(sign);bank.push_back(m);}
    if(bank.size()!=7)return 3;
    visit(Mat{1,0,0,0,1,0,0,0,1},0,1,"");
    std::cout<<"{\"tested_words\":"<<tested<<",\"trace_candidates\":"<<trace_candidates<<",\"matching_cycles\":"<<hits<<"}\n";
    return hits?1:0;
}
