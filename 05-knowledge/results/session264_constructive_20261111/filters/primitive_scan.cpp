#include <iostream>
#include <vector>
#include <cstdint>
#include <string>
using namespace std;
int main(){int n;while(cin>>n){vector<long long>r(8);vector<vector<long long>>v(n,vector<long long>(8));for(auto &x:r)cin>>x;for(int d=1;d<n;d++)for(auto &x:v[d])cin>>x; uint64_t mask=0,last=0,count=0;bool found=false;for(uint64_t t=0;t<(1ULL<<n);t++){if(t){uint64_t now=t^(t>>1);int i=__builtin_ctzll(now^last);bool add=(now>>i)&1;for(int j=0;j<n;j++)if(j!=i){int delta= j<i ? ((((now>>j)&1)&&!add)-(((last>>j)&1)&&add)) : ((add&&!((now>>j)&1))-(!add&&!((last>>j)&1)));if(delta)for(int k=0;k<8;k++)r[k]+=delta*v[abs(j-i)][k];}last=now;}bool zero=true;for(auto x:r)if(x){zero=false;break;}if(zero){cout<<"FOUND "<<n<<" ";for(int i=0;i<n;i++)cout<<((last>>i)&1);cout<<"\n";found=true;count++;if(count>=12)break;}}cout<<"FINISHED "<<n<<" hits "<<count<<endl;}}
