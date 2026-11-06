#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;using I=long long;using G=vector<unsigned>;
void need(bool ok,const char*s){if(!ok)throw runtime_error(s);}
void valid(const G&a){int n=a.size();need(n>0&&n<31,"order");for(int i=0;i<n;i++){need(a[i]<(1u<<n)&&!(a[i]>>i&1),"row");for(int j=0;j<i;j++)need(((a[i]>>j&1)+(a[j]>>i&1))==1,"pair");}}
struct Data{vector<I>h;vector<vector<I>>s,e;};
Data profile(const G&a){valid(a);int n=a.size(),N=1<<n;Data z{vector<I>(N),vector<vector<I>>(N,vector<I>(n)),vector<vector<I>>(N,vector<I>(n))};z.h[0]=1;for(int i=0;i<n;i++)z.s[1<<i][i]=z.e[1<<i][i]=1;for(int S=1;S<N;S++)for(int i=0;i<n;i++)if(S>>i&1){int T=S^(1<<i);if(T)for(int j=0;j<n;j++)if(T>>j&1){if(a[i]>>j&1)z.s[S][i]+=z.s[T][j];else z.e[S][i]+=z.e[T][j];}z.h[S]+=z.s[S][i];}return z;}
G child(const G&a,int r){vector<int>order{r};for(int i=0;i<(int)a.size();i++)if(i!=r)order.push_back(i);G u(a.size());for(int i=0;i<(int)a.size();i++)for(int j=0;j<(int)a.size();j++)if(a[order[i]]>>order[j]&1)u[i]|=1u<<j;G out{2,8,3,5};out.resize(a.size()+3);for(int j=1;j<(int)u.size();j++){int x=j+3;for(int i=0;i<4;i++)if(i||u[0]>>j&1)out[i]|=1u<<x;else out[x]|=1u<<i;for(int k=1;k<j;k++){int y=k+3;if(u[j]>>k&1)out[x]|=1u<<y;else out[y]|=1u<<x;}}valid(out);return out;}

struct Row{I B,e,V,d,s;bool operator==(const Row&o)const{return B==o.B&&e==o.e&&V==o.V&&d==o.d;}};
Row row(const G&a,int r){auto z=profile(a);int full=(1<<a.size())-1,T=full^(1<<r);I L=0,R=0,beta=0;for(int S=1;S<=full;S++){L+=z.e[S][r]*z.h[full^S];R+=z.s[S][r]*z.h[full^S];}for(int S=T;;S=(S-1)&T){beta+=z.h[S]*z.h[T^S];if(!S)break;}return {z.h.back()+L-z.s.back()[r],z.e.back()[r],R+beta,z.h[T],z.s.back()[r]};}
G pad(G a,int r,int bit){int n=a.size();for(int i=0;i<n;i++)if(!bit||i!=r)a[i]|=1u<<n;a.push_back(bit?(1u<<r):0);valid(a);return a;}
Row step(Row q,int bit){if(bit)return {q.B+q.V-q.e,q.d,2*q.V-q.d,q.d,q.s};return {q.B,0,2*q.V,q.d,q.s};}
G amalgam(const G&a,const G&u){int n=a.size();G out=a;out.resize(n+u.size()-1);for(int j=1;j<(int)u.size();j++){int x=n+j-1;for(int i=0;i<n;i++)if(i||u[0]>>j&1)out[i]|=1u<<x;else out[x]|=1u<<i;for(int k=1;k<j;k++){int y=n+k-1;if(u[j]>>k&1)out[x]|=1u<<y;else out[y]|=1u<<x;}}valid(out);return out;}
int main(){G seed[2]={{10,20,33,38,13,19},{10,12,17,52,35,7}};int words=0,graphs=0;for(int length=0;length<=4;length++)for(int mask=0;mask<(1<<length);mask++){Row states[2];I counts[2];for(int i=0;i<2;i++){G a=seed[i];Row q=row(a,0);need(q.B==111&&q.e==9&&q.V==161&&q.d==13&&q.s==6+i,"seed row");for(int k=0;k<length;k++){int bit=mask>>k&1;a=pad(a,0,bit);q=step(q,bit);}Row real=row(a,0);need(real==q&&real.s==q.s,"physical padding state");states[i]=real;G out=child(a,0);counts[i]=profile(out).h.back();need(counts[i]==3*real.B+2*real.s,"compiled readout");graphs+=2;}need(states[0]==states[1]&&counts[1]-counts[0]==2,"common base adjacent pair");words++;}
 int right=0;for(int n=1;n<=4;n++)for(int code=0;code<(1<<(n*(n-1)/2));code++){G u(n);int bit=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){if(code>>bit++&1)u[i]|=1u<<j;else u[j]|=1u<<i;}auto z=profile(u);int full=(1<<n)-1,T=full^1;I h=z.h.back(),ss=z.s.back()[0],e=z.e.back()[0],d=z.h[T],L=0,R=0,beta=0;for(int S=1;S<=full;S++){L+=z.e[S][0]*z.h[full^S];R+=z.s[S][0]*z.h[full^S];}for(int S=T;;S=(S-1)&T){beta+=z.h[S]*z.h[T^S];if(!S)break;}Row expected{13*h+148*L-4*ss-46*d,n==1?9:13*e,161*beta+13*(R-beta),13*d,0};for(int i=0;i<2;i++){G a=amalgam(seed[i],u);Row q=row(a,0);need(q==expected&&q.s==(6+i)*d,"arbitrary right row");G c=child(a,0);need(profile(c).h.back()==3*q.B+2*q.s,"arbitrary right child");graphs+=2;}right++;}
 cout<<"{\"status\":\"PASS\",\"all_binary_suffix_words_through_length\":4,\"words\":"<<words<<",\"literal_graph_counts\":"<<graphs<<",\"initial_common_row\":[111,9,161,13],\"preserved_endpoint_choices\":[6,7],\"all_labelled_right_inputs_orders_1_through_4\":75,\"coverage_claim\":false}\n";
}
