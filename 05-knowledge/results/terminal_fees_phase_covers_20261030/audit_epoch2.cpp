// Independent F_2 certificate checker. No search, no assertions.
// Transition arithmetic uses successive modular halving; density uses a sum.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using I = int64_t;
void require(bool b,const char*s){if(!b)throw std::runtime_error(s);}
int main(int argc,char**argv){try{
 require(argc==3,"usage: audit E certificate");
 int E=std::stoi(argv[1]),S=2;
 require(E>=3&&E<=9&&S>=1&&S<=4,"range");
 int q=1,r=1;for(int i=0;i<E;++i)q*=5;for(int i=0;i<S;++i)r*=5;
 I D=200LL*q,N=125LL*q-5,neg=-(1LL<<60);
 std::vector<int> residue,ordinary(r,-1),pending(r,-1);
 for(int m=1;m<r;++m)if(m%5){ordinary[m]=residue.size();residue.push_back(m);}
 int ordinary_count=residue.size();
 std::vector<int> orbit;int u=1;do{orbit.push_back(u);u=u*14%25;}while(u!=1);
 std::sort(orbit.begin(),orbit.end());
 for(int m:orbit)if(m!=1){pending[m]=residue.size();residue.push_back(m);}
 require(residue.size()==29,"mode count");
 size_t count=residue.size()*q;
 std::ifstream f(argv[2],std::ios::binary|std::ios::ate);require(bool(f),"open");
 require(f.tellg()==std::streamoff(8*count),"file length");f.seekg(0);
 std::vector<I> V(count);char bytes[8];
 for(size_t i=0;i<count;++i){f.read(bytes,8);uint64_t u=0;for(int j=0;j<8;++j)u|=uint64_t((unsigned char)bytes[j])<<(8*j);V[i]=(u>>63)?-1-I(~u):I(u);}
 require(bool(f),"read");
 std::vector<int>height(q);height[0]=E;
 for(int z=1;z<q;++z){int t=z;while(t%5==0){++height[z];t/=5;}}
 std::vector<I>density(E+1);I term=D/50;
 for(int h=3;h<=E;++h){density[h]=density[h-1]+term;term/=5;}
 I lo=0,hi=neg;uint64_t valid=0,excluded=0,edges=0,entries=0;
 for(int mode=0;mode<(int)residue.size();++mode)for(int z=0;z<q;++z){I x=V[(size_t)mode*q+z];
   if(mode>=ordinary_count&&z%125){require(x==neg,"excluded slot");++excluded;}
   else{require(x>neg/2&&x<=0,"valid finite nonpositive potential");lo=std::min(lo,x);hi=std::max(hi,x);++valid;}}
 auto half=[](I x,int modulus){x%=modulus;if(x&1)x+=modulus;return int(x/2);};
 auto edge=[&](int mode,int z,int a,int m,int zz,I source){
   int H=height[zz];
   bool active=std::find(orbit.begin(),orbit.end(),m)!=orbit.end();
   bool pend=H>=3&&active&&m!=1&&(a>1||mode>=ordinary_count);
   int target_mode=pend?pending[m]:ordinary[m];require(target_mode>=0,"target mode");
   require(!pend||zz%125==0,"domain closure");
   I gross=0;
   if(H>=1&&m%5==1)gross+=D;
   if(H>=2&&m==1)gross+=D;
   if(H>=3&&active){
     if(a>1||mode>=ordinary_count){if(m==1)gross+=D*(H-2);}
     else gross+=density[H];
   }
   I target=V[(size_t)target_mode*q+zz];
   require(source+gross-N*a<=target,"edge or entry inequality");
 };
 // x=z-1 evolves by x -> (3x+1)/2^a, independently of producer formula.
 for(int mode=0;mode<(int)residue.size();++mode){
   for(int z=0;z<q;++z){if(mode>=ordinary_count&&z%125)continue;
     I source=V[(size_t)mode*q+z];int m=3*residue[mode]%r,zz=(3LL*z+q-2)%q;
     for(int a=1;a<=64;++a){m=half(m,r);zz=half(zz,q);edge(mode,z,a,m,(zz+1)%q,source);++edges;}
   }
 }
 // Initial all-one orbit: enumerate by actual multiplication by 3 and halving.
 int initial=1;do{
   int mode=ordinary[initial],m=3*initial%r,zz=q-2;
   for(int a=1;a<=64;++a){m=half(m,r);zz=half(zz,q);if(a>=2){edge(mode,0,a,m,(zz+1)%q,0);++entries;}}
   initial=half(3LL*initial,r);
 }while(initial!=1);
 require(hi+D*E-65*N<=lo,"transition tail");require(D*E-65*N<=lo,"entry tail");
 I terminal=V[(size_t)ordinary[1]*q];
 std::cout<<"{\"status\":\"PASS\",\"E\":"<<E<<",\"multiplier_precision\":"<<S
 <<",\"valid_states\":"<<valid<<",\"excluded_slots\":"<<excluded<<",\"edges\":"<<edges
 <<",\"entries\":"<<entries<<",\"terminal_numerator\":"<<terminal<<",\"denominator\":"<<D
 <<",\"minimum\":"<<lo<<",\"maximum\":"<<hi<<",\"terminal_fee_two\":"<<(terminal<=-2*D?"true":"false")<<"}\n";
 return 0;
 }catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}
