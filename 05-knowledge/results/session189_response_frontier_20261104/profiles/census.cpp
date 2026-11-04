#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
using namespace std;
bool prime(const array<unsigned,7>&r){for(unsigned m=1;m<127;m++){if(__builtin_popcount(m)<2)continue;bool yes=true;for(int i=0;i<7;i++)if(!(m>>i&1)){unsigned x=r[i]&m;if(x!=0&&x!=m){yes=false;break;}}if(yes)return false;}return true;}
int code(const array<unsigned,7>&r,const array<int,7>&p){int c=0,b=0;for(int i=0;i<7;i++)for(int j=i+1;j<7;j++,b++)if(r[p[i]]>>p[j]&1)c|=1<<b;return c;}
int main(int argc,char**argv){if(argc!=2)return 2;int n;cin>>n;vector<uint16_t> owner(1<<21);vector<uint64_t>orbit(n);for(int k=0;k<n;k++){array<unsigned,7>r;for(auto&x:r)cin>>x;if(!prime(r))return 3;array<int,7>p={0,1,2,3,4,5,6};array<uint64_t,7>S{},E{};do{int c=code(r,p);if(owner[c]&&owner[c]!=k+1)return 4;if(!owner[c]){owner[c]=k+1;orbit[k]++;}bool ok=true;for(int j=0;j<6;j++)if(!(r[p[j]]>>p[j+1]&1)){ok=false;break;}if(ok){S[p[0]]++;E[p[6]]++;}}while(next_permutation(p.begin(),p.end()));cout<<orbit[k];for(auto x:S)cout<<" "<<x;for(auto x:E)cout<<" "<<x;cout<<"\n";}uint64_t pc=0;for(int c=0;c<(1<<21);c++){array<unsigned,7>r{};int b=0;for(int i=0;i<7;i++)for(int j=i+1;j<7;j++,b++)if(c>>b&1)r[i]|=1<<j;else r[j]|=1<<i;bool pr=prime(r);pc+=pr;if(pr!=(owner[c]!=0)){cerr<<"coverage mismatch "<<c<<"\n";return 5;}}ofstream f(argv[1],ios::binary);for(auto x:owner){f.put(x&255);f.put(x>>8);}f.close();cerr<<"complete labeled prime count "<<pc<<"\n";}
