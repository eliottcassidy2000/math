#include <iostream>
#include <vector>
#include <cstdint>
using namespace std;using U=uint64_t;
vector<U> dp(const vector<unsigned>&r){int n=r.size(),N=1<<n;vector<U>out(N*n);for(int v=0;v<n;v++)out[(1<<v)*n+v]=1;for(int mask=1;mask<N;mask++)for(int v=0;v<n;v++){U x=out[mask*n+v];if(!x)continue;unsigned todo=r[v]&~mask;while(todo){unsigned b=todo&-todo;todo-=b;int u=__builtin_ctz(b);out[(mask|b)*n+u]+=x;}}return out;}
bool reaches(const vector<unsigned>&r){unsigned seen=1,last=0,full=(1u<<r.size())-1;while(last!=seen){last=seen;for(int i=0;i<(int)r.size();i++)if(seen>>i&1)seen|=r[i];}return seen==full;}
int main(){int n,t;while(cin>>n>>t){vector<unsigned>r(n),rev(n);for(auto&x:r)cin>>x;for(int i=0;i<n;i++)rev[i]=r[i]^(((1u<<n)-1)^(1u<<i));auto E=dp(r),S=dp(rev);int full=(1<<n)-1;vector<U>H(1<<n);H[0]=1;for(int mask=1;mask<=full;mask++)for(int v=0;v<n;v++)H[mask]+=E[mask*n+v];U B=0;int rest=full^(1<<t);for(int a=rest;;a=(a-1)&rest){B+=H[a]*S[(full^a)*n+t];if(!a)break;}cout<<H[full]<<' '<<E[full*n+t]<<' '<<B<<' '<<(reaches(r)&&reaches(rev))<<'\n';}}
