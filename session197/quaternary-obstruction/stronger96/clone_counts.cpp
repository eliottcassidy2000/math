#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
using namespace std;
long long H(const vector<int>&a,long long cap){int n=a.size(),full=(1<<n)-1;vector<long long>d((1<<n)*n);for(int v=0;v<n;v++)d[(1<<v)*n+v]=1;for(int s=1;s<=full;s++)for(int u=0;u<n;u++)if(d[s*n+u])for(int z=a[u]&~s;z;z&=z-1){int b=z&-z,v=__builtin_ctz((unsigned)b);auto &t=d[(s|b)*n+v];t=min(cap,t+d[s*n+u]);}long long h=0;for(int v=0;v<n;v++)h=min(cap,h+d[full*n+v]);return h;}
vector<int> clone(const vector<int>&a,int r,int m){vector<int>lab;for(int i=0;i<(int)a.size();i++){if(i==r)for(int j=0;j<m;j++)lab.push_back(i);else lab.push_back(i);}vector<int>b(lab.size());for(int i=0;i<(int)lab.size();i++)for(int j=i+1;j<(int)lab.size();j++){bool f=lab[i]==lab[j]||a[lab[i]]>>lab[j]&1;b[f?i:j]|=1<<(f?j:i);}return b;}
int main(int argc,char**argv){if(argc!=2)return 1;ifstream f(argv[1]);string line;long long graphs=0,eligible=0,roots=0,first=0,second=0,third=0;while(getline(f,line)){replace(line.begin(),line.end(),',',' ');istringstream in(line);int n,h;in>>n>>h;vector<int>a(n);for(int &x:a)in>>x;graphs++;if(H(a,362)!=h)throw runtime_error("bank count");for(int beta:{-1,1}){if((h-beta)%4)continue;int A=(h-beta)/4;if(A<=0||A>90||A%6)continue;eligible++;for(int r=0;r<n;r++){roots++;if(H(clone(a,r,2),16*A+beta+1)!=16*A+beta)continue;first++;if(H(clone(a,r,3),64*A+beta+1)!=64*A+beta)continue;second++;if(H(clone(a,r,4),256*A+beta+1)!=256*A+beta)continue;third++;cout<<"HIT "<<n<<" "<<h<<" "<<A<<" "<<beta<<" "<<r<<"\n";}}}cout<<"PASS graphs "<<graphs<<" eligible_graph_signs "<<eligible<<" roots "<<roots<<" m2_survivors "<<first<<" m3_survivors "<<second<<" m4_survivors "<<third<<"\n";if(third)throw runtime_error("candidate exists");}
