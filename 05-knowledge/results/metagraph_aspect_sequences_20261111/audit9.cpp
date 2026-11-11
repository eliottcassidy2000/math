// Independent all-ordinary-source seed replay: push path DP, scalar subset
// convolution for U, and Gray-code selected-set updates for W.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
using namespace std;using I=int64_t;
void need(bool q,const char*s){if(!q)throw runtime_error(s);}
struct Witness{uint64_t code;int dual,mask;I H,U;int exponent;};
int main(int argc,char**argv){need(argc==3,"audit9 input witnesses.tsv");ifstream f(argv[1]);uint64_t code,conv;int aut,scH;array<int,9> original;map<I,Witness> cert;
 uint64_t classes=0,sources=0,units=0,orbits=0;constexpr int n=9,N=512,full=511;
 while(f>>code>>conv>>aut>>scH){for(auto&r:original)f>>r;
  for(int isdual=0;isdual<(code==conv?1:2);isdual++){
   classes++;array<int,9> rows=original;for(int v=0;v<n;v++)if(isdual)rows[v]=full^rows[v]^(1<<v);
   array<array<I,9>,N> end{},begin{};array<I,N> h{};h[0]=1;
   for(int v=0;v<n;v++)end[1<<v][v]=begin[1<<v][v]=1;
   for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(S>>v&1){
    h[S]+=end[S][v];int unused=full^S;
    for(int u=0;u<n;u++)if(unused>>u&1){
     if(rows[v]>>u&1)end[S|(1<<u)][u]+=end[S][v];
     else begin[S|(1<<u)][u]+=begin[S][v];
    }
   }
   I H=h[full];array<I,9> u{};array<array<I,9>,9> Q{};
   for(int S=1;S<N;S++)for(int j=0;j<n;j++)u[j]+=begin[S][j]*h[full^S];
   for(int S=1;S<full;S++)for(int i=0;i<n;i++)if(end[S][i])for(int j=0;j<n;j++)if(begin[full^S][j]){
    I t=end[S][i]*begin[full^S][j];Q[i][j]+=t;Q[j][i]+=t;
   }
   int previous=0,count=0;I W=0,U=0;
   for(int k=1;k<N;k++){
    int O=k^(k>>1),j=__builtin_ctz((unsigned)(O^previous));int sign=(O>>j&1)?1:-1;
    for(int i=0;i<n;i++)if(previous>>i&1)W+=sign*Q[i][j];U+=sign*u[j];previous=O;
    if(W==H-1&&U>0){count++;I c=U;int e=1;while(c%2==0)c/=2,e++;
     Witness w{code,isdual,O,H,U,e};auto it=cert.find(c);if(it==cert.end()||e<it->second.exponent)cert[c]=w;
    }
   }
   need(count%aut==0,"free seed action");sources+=count>0;units+=count;orbits+=count/aut;
  }
 }
 ofstream out(argv[2]);out<<"c\tbase_order\tcanonical_code\tconverse\tselected_mask\tH\tU\tfirst_exponent\n";
 for(auto[c,w]:cert)out<<c<<"\t9\t"<<w.code<<'\t'<<w.dual<<'\t'<<w.mask<<'\t'<<w.H<<'\t'<<w.U<<'\t'<<w.exponent<<'\n';
 cout<<"{\"ordinary_sources_checked\":"<<classes<<",\"unit_sources\":"<<sources<<",\"unit_selections_sum\":"<<units<<",\"unit_orbits_sum\":"<<orbits<<",\"unit_odd_coefficients\":"<<cert.size()<<",\"missing_through_4097\":[";
 bool first=true;for(int c=1;c<=4097;c+=2)if(!cert.count(c)){if(!first)cout<<',';first=false;cout<<c;}cout<<"]}\n";
}
