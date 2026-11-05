#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
int rankat(const vector<int>&a,int r){vector<int>O,I;int n=a.size();for(int v=0;v<n;v++)if(v!=r){if(a[r]>>v&1)O.push_back(v);else I.push_back(v);}vector<int>match(I.size(),-1);int rank=0;for(int u:O){vector<int>vis(I.size());auto aug=[&](auto&&self,int v)->bool{for(int j=0;j<(int)I.size();j++)if(!vis[j]&&(a[v]>>I[j]&1)){vis[j]=1;if(match[j]<0||self(self,match[j])){match[j]=v;return true;}}return false;};rank+=aug(aug,u);}return rank;}
uint64_t addcap(uint64_t a,uint64_t b){const uint64_t L=1ULL|(1ULL<<12)|(1ULL<<24)|(1ULL<<36)|(1ULL<<48),mask=(1ULL<<60)-1;uint64_t s=a+b,f=((s>>10)|(s>>11))&L;return((s&~(4095*f))|(f<<10))&mask;}
array<int,5>runs(const vector<int>&a,int r){int N=a.size()-1,full=(1<<N)-1;vector<int>E;for(int v=0;v<(int)a.size();v++)if(v!=r)E.push_back(v);int O=0;vector<int>b(N);for(int u=0;u<N;u++){if(a[r]>>E[u]&1)O|=1<<u;for(int v=0;v<N;v++)if(a[E[u]]>>E[v]&1)b[u]|=1<<v;}int I=full^O;vector<uint64_t>d((1<<N)*N);for(int v=0;v<N;v++)d[(1<<v)*N+v]=1+((O>>v&1)?1ULL<<12:0);for(int S=1;S<=full;S++)for(int v=0;v<N;v++){uint64_t z=d[S*N+v];if(!z)continue;for(int D=b[v]&~S;D;D&=D-1){int bit=D&-D;auto &w=d[(S|bit)*N+__builtin_ctz((unsigned)bit)];w=addcap(w,z);}if(I>>v&1)for(int D=O&~S;D;D&=D-1){int bit=D&-D;auto&w=d[(S|bit)*N+__builtin_ctz((unsigned)bit)];w=addcap(w,z<<12);}}uint64_t z=0;for(int v=0;v<N;v++){z=addcap(z,d[full*N+v]);if(I>>v&1)z=addcap(z,d[full*N+v]<<12);}array<int,5>x{};for(int k=0;k<5;k++)x[k]=(z>>(12*k))&4095;return x;}
int main(int argc,char**argv){if(argc!=2)return 1;ifstream in(argv[1]);long long tested=0,hits=0;string line;while(getline(in,line)){replace(line.begin(),line.end(),',',' ');vector<int>x;size_t at=0;while(at<line.size()){while(at<line.size()&&line[at]==' ')at++;if(at==line.size())break;size_t z;x.push_back(stoi(line.substr(at),&z));at+=z;}vector<int>a(x.begin()+2,x.end());for(int r=0;r<x[0];r++)if(rankat(a,r)==3){auto v=runs(a,r);tested++;int A=v[4],beta=v[1]-4*A;if(A>0&&A<=125&&v[3]==4*A&&v[2]==6*A&&(beta==1||beta==-1)){hits++;cout<<"HIT order "<<x[0]<<" C "<<x[1]<<" A "<<A<<" beta "<<beta<<" root "<<r<<" rows";for(int row:a)cout<<" "<<row;cout<<endl;}}}cout<<"DONE rank3_roots "<<tested<<" hits "<<hits<<"\n";return 0;}
