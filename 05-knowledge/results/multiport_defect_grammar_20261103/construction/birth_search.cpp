#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;using I=int64_t;using T=vector<unsigned>;
pair<I,vector<I>>response(const T&t){int n=t.size(),full=(1<<n)-1;vector<I>d((full+1)*n),h(full+1);h[0]=1;
#ifdef FORWARD_AUDIT
for(int v=0;v<n;v++)d[(1<<v)*n+v]=1;
for(int s=1;s<=full;s++)for(unsigned q=s;q;q&=q-1){int v=countr_zero(q);I k=d[s*n+v];h[s]+=k;if(k)for(unsigned out=t[v]&~s;out;out&=out-1)d[(s|(out&-out))*n+countr_zero(out)]+=k;}
#else
for(int s=1;s<=full;s++)for(unsigned q=s;q;q&=q-1){int v=countr_zero(q);int previous=s^(1<<v);I k=previous?0:1;for(unsigned z=previous;z;z&=z-1){int u=countr_zero(z);if(t[u]>>v&1)k+=d[previous*n+u];}d[s*n+v]=k;h[s]+=k;}
#endif
vector<I>f(n);for(int s=1;s<=full;s++)for(unsigned q=s;q;q&=q-1){int w=countr_zero(q);f[w]+=d[s*n+w]*h[full^s];}return{h[full],f};}
T birth(T t,unsigned cut){int n=t.size();for(int v=0;v<n;v++)if(cut>>v&1)t[v]|=1<<n;t.push_back(((1<<n)-1)^cut);return t;}
bool strong(const T&t){unsigned full=(1u<<t.size())-1;for(bool reverse:{false,true}){unsigned seen=1;while(true){unsigned next=seen;for(int v=0;v<(int)t.size();v++)if(seen>>v&1)next|=reverse?(full^t[v]^(1u<<v)):t[v];if(next==seen)break;seen=next;}if(seen!=full)return false;}return true;}
int main(){int n,w,moving;cin>>n>>w>>moving;if(n<1||n>8)return 2;T source(n);for(auto&v:source)cin>>v;auto[h,fs]=response(source);I f=fs[w];int tested=0,hits=0;cout<<"{\"source_H\":"<<h<<",\"source_F\":"<<f<<",\"moving_mark\":"<<(moving?"true":"false")<<",\"hits\":[";
for(unsigned x=0;x<(1u<<n);x++){auto X=birth(source,x);for(unsigned y=0;y<(1u<<(n+1));y++){auto Y=birth(X,y);tested++;if(!strong(Y))continue;auto[hy,fy]=response(Y);if(hy!=3*h+2&&hy!=3*h+4)continue;for(int mark=0;mark<n+2;mark++){if((!moving&&mark!=w)||fy[mark]!=3*f)continue;if(hits++)cout<<',';cout<<"{\"offset\":"<<hy-3*h<<",\"mark\":"<<mark<<",\"cut_x\":"<<x<<",\"cut_y\":"<<y<<",\"intermediate_strong\":"<<(strong(X)?"true":"false")<<",\"rows\":[";for(int i=0;i<n+2;i++){if(i)cout<<',';cout<<Y[i];}cout<<"]}";}}}
cout<<"],\"tested\":"<<tested<<"}\n";}
