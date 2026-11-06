#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std; using ll=long long; using U=uint32_t;
void need(bool p,const string&s){if(!p)throw runtime_error(s);}
struct Family{string name;vector<U>P,T;U mask;int extra,i,j;};
const vector<Family> F={
 {"A",{2,24,3,37,13,23},{2,8,67,101,15,87,19},15,1,1,4},
 {"B",{32,5,1,39,47,70,31},{0,5,129,39,47,199,159,27},5,0,0,5},
 {"C",{2,4,41,3,15,91,31},{0,5,41,131,15,219,159,23},27,1,0,1}};
void valid(const vector<U>&a){int n=a.size();for(int i=0;i<n;i++){need(a[i]<(1u<<n)&&!(a[i]>>i&1),"row");for(int j=0;j<i;j++)need(((a[i]>>j&1)+(a[j]>>i&1))==1,"pair");}}
vector<U> module(const vector<U>&core,U mask,const vector<U>&u){int p=core.size(),n=u.size();auto a=core;a.resize(p+n);for(int j=0;j<n;j++){a[p+j]=u[j]<<p;for(int i=0;i<p;i++)if(mask>>i&1)a[p+j]|=1u<<i;else a[i]|=1u<<(p+j);}valid(a);return a;}
vector<U> del(const vector<U>&a,int r){vector<int>ix;for(int i=0;i<(int)a.size();i++)if(i!=r)ix.push_back(i);vector<U>b(ix.size());for(int i=0;i<(int)ix.size();i++)for(int j=0;j<(int)ix.size();j++)if(a[ix[i]]>>ix[j]&1)b[i]|=1u<<j;valid(b);return b;}
vector<ll> all_h(const vector<U>&a){int n=a.size(),N=1<<n;vector<ll>d(N*n),h(N);h[0]=1;for(int v=0;v<n;v++)d[(1<<v)*n+v]=1;for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(S>>v&1){int R=S^(1<<v);for(int w=0;w<n;w++)if((R>>w&1)&&(a[w]>>v&1))d[S*n+v]+=d[R*n+w];h[S]+=d[S*n+v];}return h;}
pair<ll,ll> counts(const vector<U>&a){int n=a.size(),N=1<<n;vector<ll>d(N*n),c(N*n);for(int v=0;v<n;v++)d[(1<<v)*n+v]=1;c[n]=1;for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(S>>v&1){U todo=a[v]&~S&((1u<<n)-1);while(todo){U b=todo&-todo;todo-=b;int w=__builtin_ctz(b);d[(S|b)*n+w]+=d[S*n+v];c[(S|b)*n+w]+=c[S*n+v];}}ll h=0,k=0;for(int v=0;v<n;v++){h+=d[(N-1)*n+v];if(a[v]&1)k+=c[(N-1)*n+v];}return {h,k};}
void edit_check(const Family&f,const vector<U>&p,const vector<U>&t,int un){int a=f.P.size(),b=f.T.size();auto map=[&](int v){return v<a?v:v+1;};for(int i=0;i<a+un;i++)for(int j=i+1;j<a+un;j++){bool flip=(i==min(f.i,f.j)&&j==max(f.i,f.j));need(((p[i]>>j&1)^(t[map(i)]>>map(j)&1))==flip,"old edit");}need(b==a+1,"birth");}
int main(){int modules=0,graphs=0,rewrites=0;ll largest=0;vector<string> samples;
 for(int n=0;n<=4;n++){int e=n*(n-1)/2;for(int code=0;code<(1<<e);code++){vector<U>u(n);int bit=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){if(code>>bit++&1)u[i]|=1u<<j;else u[j]|=1u<<i;}valid(u);auto h=all_h(u);int full=(1<<n)-1;vector<ll>b2(1<<n);for(int S=0;S<=full;S++)for(int A=S;;A=(A-1)&S){b2[S]+=h[A]*h[S^A];if(!A)break;}ll b3=0;for(int A=full;;A=(A-1)&full){b3+=h[A]*b2[full^A];if(!A)break;}ll H=h[full],B2=b2[full],B3=b3,B0=n==0;modules++;
 vector<U>up=u;for(auto &row:up)row|=1u<<n;up.push_back(0);valid(up);auto hpadded=all_h(up);int fp=(1<<(n+1))-1;ll paddedB2=0;for(int S=0;S<=fp;S++)paddedB2+=hpadded[S]*hpadded[fp^S];need(hpadded[fp]==H&&paddedB2==2*B2,"padding bin law");if(n){need(B2>=(n+1)*H && B2%2==0,"coupling bounds");auto rebased=module(F[2].P,F[2].mask,up);auto rec=counts(rebased);graphs++;need(rec.first==54*B2-H&&rec.second==H,"C count rebase");}else{auto rebased=module(F[2].P,F[2].mask,up);need(counts(rebased).first==53,"empty rebase hostile");graphs++;}
 for(const auto&f:F){auto p=module(f.P,f.mask,u),t=module(f.T,f.mask|(U(f.extra)<<f.P.size()),u);edit_check(f,p,t,n);auto hp=counts(p),ht=counts(t);graphs+=2;ll ep,et,ecp,ect;
 if(f.name=="A"){ep=H+17*B2+8*B3-B0;et=3*H+34*B2+16*B3-2*B0;ecp=H+B2;ect=H+2*B2;}
 else if(f.name=="B"){ep=5*H+20*B2;et=11*H+40*B2;ecp=H;ect=3*H;}
 else{ep=27*B2-H-B0;et=54*B2-H-2*B0;ecp=H;ect=3*H;}
 need(hp.first==ep&&ht.first==et,"path formula");need(ht.first==2*hp.first+H,"correction");need(hp.second==ecp,"parent cycles");
 if(f.name!="A"){need(t[0]==0,"sink");auto r=del(t,0);auto hr=counts(r);graphs++;need(hr.first==et&&ht.second==0&&hr.second==ect,"reduced child");need(del(p,0)==del(del(t,0),int(f.T.size())-2),"retained graph");}else need(ht.second==ect,"A cycles");
 largest=max({largest,hp.first,ht.first});rewrites++;
 if(n==0||(n==3&&code==2))samples.push_back("{\"family\":\""+f.name+"\",\"module_order\":"+to_string(n)+",\"module_code\":"+to_string(code)+",\"h\":"+to_string(H)+",\"B2\":"+to_string(B2)+",\"B3\":"+to_string(B3)+",\"parent_H\":"+to_string(hp.first)+",\"child_H\":"+to_string(ht.first)+",\"parent_C\":"+to_string(ecp)+",\"child_C\":"+to_string(ect)+"}");
 }}}
 cout<<"{\n  \"status\":\"PASS\",\n  \"module_universe\":\"Every labelled tournament of orders 0 through 4\",\n  \"modules\":"<<modules<<",\n  \"graphs\":"<<graphs<<",\n  \"rewrite_pairs\":"<<rewrites<<",\n  \"largest_path_count\":"<<largest<<",\n  \"samples\":[\n";for(int i=0;i<(int)samples.size();i++)cout<<"    "<<samples[i]<<(i+1<(int)samples.size()?",":"")<<"\n";cout<<"  ]\n}\n";
}
