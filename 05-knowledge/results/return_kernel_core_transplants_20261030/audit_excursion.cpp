#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <cstdlib>
using namespace std;
int main(int argc,char**argv){
 if(argc!=5)return 2;int E=stoi(argv[1]),cut=stoi(argv[4]);if(E<1||E>4)return 2;
 int z=1;for(int k=0;k<E;k++)z*=5;int q=5*z,S=q*q;long long D=8*q,N=5*q-5,INF=numeric_limits<long long>::max()/4;
 vector<long long> h(2*S,INF),c;string cell;stringstream ss(argv[3]);while(getline(ss,cell,','))c.push_back(stoll(cell));if(c.size()!=25)return 3;
 ifstream input(argv[2]);int f,m,b;long long d,lo=0,hi=*max_element(c.begin(),c.end());long long rows=0;
 while(input>>f>>m>>b>>d){if(f<0||f>1||m<=0||m>=q||m%5==0||b<0||b>=q)return 4;if(d!=INF && (d < -1000000000000000LL || d > 1000000000000000LL))return 9;h[f*S+m*q+b]=d;rows++;if(d!=INF){lo=min(lo,d);hi=max(hi,d);}}
 if(h[q]!=0){cerr<<"root potential changed\n";return 5;}
 vector<int> v(q);for(int k=0;k<q;k++){int x=k,t=0;if(x==0)t=E;else while(x%5==0&&t<E){x/=5;t++;}v[k]=t;}
 long long tested=0,reachable=0;
 for(int flag=0;flag<2;flag++)for(int u=1;u<q;u++)if(u%5)for(int t=0;t<q;t++){
  int idx=flag*S+u*q+t;bool terminal=u%z==1&&t%z==0;
  if(terminal && idx!=q)continue;if(h[idx]==INF)continue;reachable++;
  int im=3*u%q,ib=(3*t+1)%q;
  for(int a=1;a<=cut;a++){
   im=(im+(im%2)*q)/2;ib=(ib+(ib%2)*q)/2; // exact modular halving, separate from producer inverses
   int flag2=flag || a!=1;long long edge=N*a-D*min(v[(im+q-1)%q],v[ib]),want;
   bool final=im%z==1&&ib%z==0;
   if(final){if(!flag2)continue;want=c[5*((im-1)/z)+ib/z];}
   else want=h[flag2*S+im*q+ib];
   if(want==INF || want>h[idx]+edge){cerr<<"Bellman failure "<<flag<<" "<<u<<" "<<t<<" "<<a<<"\n";return 6;}
   tested++;
  }
 }
 // Every omitted a is >1, hence targets the fully reachable non-one channel.
 // All residue representatives reduce by period4z; no-reduction tail test is stronger.
 long long tail=lo+N*(cut+1)-D*E-hi;
 if(tail<0){cerr<<"all-letter tail failed\n";return 7;}
 for(int u=1;u<q;u++)if(u%5)for(int t=0;t<q;t++)if(!(u%z==1&&t%z==0)&&h[S+u*q+t]==INF)return 8;
 cout<<"{\"E\":"<<E<<",\"bellman_edges\":"<<tested<<",\"reachable_states\":"<<reachable<<",\"input_rows\":"<<rows<<",\"tail_margin\":"<<tail<<"}\n";
}
