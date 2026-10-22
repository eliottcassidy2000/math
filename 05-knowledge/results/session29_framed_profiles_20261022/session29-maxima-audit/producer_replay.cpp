#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <vector>
using namespace std;
struct Term {vector<unsigned> wedge; int v,e;};
void run(int m,unsigned q){int E=m*(m-1)/2,id[8][8],ctrl[28],e=0;vector<vector<int>>stars(m);for(int i=0;i<m;i++)for(int j=i+1;j<m;j++){id[i][j]=id[j][i]=e;ctrl[e]=(q>>e&1)?j:i;stars[ctrl[e]].push_back(e++);}vector<int>free;for(auto s:stars)for(int j=1;j<(int)s.size();j++)free.push_back(s[j]);vector<Term>ts;vector<int>p(m);iota(p.begin(),p.end(),0);do{Term t;for(int j=1;j<m-1;j++){int a=id[p[j-1]][p[j]],b=id[p[j]][p[j+1]];if(ctrl[a]==p[j]&&ctrl[b]==p[j])t.wedge.push_back((1u<<a)|(1u<<b));}t.e=id[p[m-2]][p[m-1]];t.v=ctrl[t.e]==p[m-1]?p[m-1]:-1;ts.push_back(t);}while(next_permutation(p.begin(),p.end()));long long best=-1,bA=0,bB[8]={};unsigned bsign=0;for(unsigned mask=0;mask<(1u<<free.size());mask++){unsigned sign=0;for(unsigned j=0;j<free.size();j++)if(mask>>j&1)sign|=1u<<free[j];long long A=0,B[8]={};for(auto &t:ts){bool ok=true;for(auto w:t.wedge){auto x=sign&w;if(!x||x==w){ok=false;break;}}if(!ok)continue;int wt=1<<t.wedge.size();A+=wt;if(t.v>=0)B[t.v]+=(sign>>t.e&1)?-wt:wt;}long long f=A;for(int v=0;v<m;v++)f+=abs(B[v]);if(f>best){best=f;bA=A;bsign=sign;copy(B,B+m,bB);}}unsigned finalsign=bsign;for(int v=0;v<m;v++)if(bB[v]<0)for(int e:stars[v])finalsign^=1u<<e;cout<<m<<" "<<q<<" "<<best<<" "<<bA<<" "<<finalsign<<" "<<free.size();for(int v=0;v<m;v++)cout<<" "<<bB[v];cout<<endl;}
int oldmain(int argc,char**argv){int m=atoi(argv[1]),N=atoi(argv[2]);unsigned seed=290063;mt19937 gen(seed);int E=m*(m-1)/2;set<unsigned>qs={0,2}; // transitive and initial C3, then ordered tail
// consecutive triangles, cycles, converse, and seeded random labelled controllers
unsigned tri=0;int e=0;for(int i=0;i<m;i++)for(int j=i+1;j<m;j++,e++)if(i/3==j/3&&j==i+2)tri|=1u<<e;qs.insert(tri);unsigned cyc=0;e=0;for(int i=0;i<m;i++)for(int j=i+1;j<m;j++,e++)if(j-i>m/2)cyc|=1u<<e;qs.insert(cyc);qs.insert(((1u<<E)-1)^cyc);while((int)qs.size()<N)qs.insert(gen()&((1u<<E)-1));for(auto q:qs)run(m,q);return 0;}

int main(int argc,char**argv){int m=atoi(argv[1]);if(argc>2)return oldmain(argc,argv);int E=m*(m-1)/2;vector<char>seen(1<<E);vector<pair<int,int>>es;int ids[8][8];for(int i=0;i<m;i++)for(int j=i+1;j<m;j++){ids[i][j]=ids[j][i]=es.size();es.push_back({i,j});}for(unsigned q=0;q<(1u<<E);q++)if(!seen[q]){run(m,q);vector<int>p(m);iota(p.begin(),p.end(),0);do{unsigned r=0;for(int k=0;k<E;k++){auto [i,j]=es[k];int a=p[i],b=p[j];bool rev=(q>>k&1)^(a>b);if(rev)r|=1u<<ids[a][b];}seen[r]=1;}while(next_permutation(p.begin(),p.end()));}}
