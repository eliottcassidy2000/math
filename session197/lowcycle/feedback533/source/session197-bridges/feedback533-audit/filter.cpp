#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;
long long H(const vector<int>&a){int n=a.size(),M=1<<n;vector<vector<long long>>d(M,vector<long long>(n));for(int S=1;S<M;S++)for(int v=0;v<n;v++)if(S>>v&1){int R=S^(1<<v);if(!R){d[S][v]=1;continue;}for(int U=R&a[v];U;U&=U-1)d[S][v]+=d[R][__builtin_ctz((unsigned)(U&-U))];}long long z=0;for(auto x:d.back())z+=x;return z;}
int main(int argc,char**argv){if(argc!=2)return 2;ifstream in(argv[1]);long long count=0,hits=0;string line;while(getline(in,line)){replace(line.begin(),line.end(),',',' ');istringstream s(line);int n,c;s>>n>>c;vector<int>a(n);for(auto&x:a)s>>x;if(c>2)continue;count++;if(H(a)!=533)continue;hits++;cout<<n<<','<<c;for(int x:a)cout<<','<<x;cout<<'\n';}cerr<<"graphs "<<count<<" hits "<<hits<<'\n';if(count!=12077||hits!=6)throw runtime_error("unexpected complete inventory");}
