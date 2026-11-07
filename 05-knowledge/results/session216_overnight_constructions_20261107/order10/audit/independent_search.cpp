#include <fstream>
#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <set>
using namespace std;using I=long long;using P=array<I,9>;
int main(int argc,char**argv){ifstream bank(argv[1]),cand(argv[2]);map<pair<int,int>,P>expected;int ix,mask;while(cand>>ix>>mask){P p;for(auto&x:p)cand>>x;if(!expected.emplace(make_pair(ix,mask),p).second)throw runtime_error("duplicate input");}int total;bank>>total;I retained=0,cutchecks=0;set<int>starts{2,4,8,10,14,16,20,22};
for(ix=0;ix<total;ix++){int n;bank>>n;vector<unsigned>a(n),rev(n);for(auto&x:a)bank>>x;for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[j]>>i&1)rev[i]|=1u<<j;int N=1<<n,F=N-1;vector<I>E(N*n),S(N*n),H(N),C(N);H[0]=1;for(int j=0;j<n;j++)E[(1<<j)*n+j]=S[(1<<j)*n+j]=1;
// Independent forward path extension in the graph and its converse.
for(int A=1;A<N;A++)for(int j=0;j<n;j++){I x=E[A*n+j],y=S[A*n+j];H[A]+=x;unsigned todo=a[j]&~unsigned(A);if(x)while(todo){unsigned z=todo&-todo;todo-=z;E[(A|z)*n+__builtin_ctz(z)]+=x;}todo=rev[j]&~unsigned(A);if(y)while(todo){unsigned z=todo&-todo;todo-=z;S[(A|z)*n+__builtin_ctz(z)]+=y;}}
I d=H[F];if(d%3)continue;I b=0;for(int A=0;A<N;A++)b+=H[A]*H[F^A];if(b%12)continue;
// Deliberately omit the producer's minimum-positive-start early pruning.
for(int A=0;A<N;A++){int B=A;while(true){C[A]+=H[B]*H[A^B];if(!B)break;B=(B-1)&A;}}
vector<I>lcoef(n),rcoef(n),tcoef(n);for(int j=0;j<n;j++)for(int A=1;A<N;A++){lcoef[j]+=E[A*n+j]*H[F^A];rcoef[j]+=S[A*n+j]*H[F^A];tcoef[j]+=S[A*n+j]*C[F^A];}
for(mask=0;mask<N;mask++){I st=0,en=0,L=d,R=d,T=b;for(int j=0;j<n;j++)if(mask>>j&1){st+=S[F*n+j];R+=rcoef[j];T+=tcoef[j];}else{en+=E[F*n+j];L+=lcoef[j];}cutchecks++;if(!starts.count(st)||en%3||L%12!=6||R%3||T%6)continue;
vector<I>incoming(N),outgoing(N),suffix(N);for(int A=1;A<N;A++)for(int j=0;j<n;j++)if(mask>>j&1)outgoing[A]+=S[A*n+j];else incoming[A]+=E[A*n+j];I h=st+en;for(int A=1;A<N;A++)h+=incoming[A]*outgoing[F^A];if(h%3)continue;
// Convolve outgoing-root suffix with ordinary leftover bin before joining prefix.
for(int A=1;A<N;A++){int B=A;while(B){suffix[A]+=outgoing[B]*H[A^B];B=(B-1)&A;}}
I cover=L+R-d;for(int A=1;A<N;A++)cover+=incoming[A]*suffix[F^A];cover*=2;if(cover%12)continue;P p{h,st,en,L,R,d,cover,b,T};auto it=expected.find({ix,mask});if(it==expected.end()||it->second!=p){cerr<<"mismatch "<<ix<<" "<<mask<<endl;return 2;}expected.erase(it);retained++;}
}if(!expected.empty())return 3;cout<<"{\"deletion_cores\":"<<total<<",\"cuts_without_minstart_prune\":"<<cutchecks<<",\"exact_surviving_profiles\":"<<retained<<"}"<<endl;}
