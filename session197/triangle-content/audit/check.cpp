#include <iostream>
#include <vector>
#include <numeric>
#include <functional>
#include <stdexcept>
using namespace std;using I=long long;using Graph=vector<int>;using Poly=vector<I>;
void need(bool b){if(!b)throw runtime_error("independent check failed");}
Graph tournament(int n,int code){Graph a(n);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(code>>k&1)a[i]|=1<<j;else a[j]|=1<<i;return a;}
Graph tr(int n){Graph a(n);for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)a[i]|=1<<j;return a;}
Graph cyclic(Graph a,Graph b,Graph c){int n=a.size(),m=b.size(),k=c.size();Graph r;for(int x:a)r.push_back(x|(((1<<m)-1)<<n));for(int x:b)r.push_back((x<<n)|(((1<<k)-1)<<(n+m)));for(int x:c)r.push_back((x<<(n+m))|((1<<n)-1));return r;}
Graph W(Graph a){return cyclic(a,{0},{0});}Graph E(Graph a){return cyclic(a,tr(2),{0});}
bool strong(Graph a){int n=a.size();for(int i=0;i<n;i++){int s=1<<i;while(true){int t=s;for(int j=0;j<n;j++)if(s>>j&1)t|=a[j];if(t==s)break;s=t;}if(s!=(1<<n)-1)return false;}return true;}
I walks(Graph a){int n=a.size();if(!n)return 1;__int128 out=0;for(int S=1;S<(1<<n);S++){vector<__int128>v(n);for(int i=0;i<n;i++)if(S>>i&1)v[i]=1;for(int len=1;len<n;len++){vector<__int128>w(n);for(int i=0;i<n;i++)if(v[i])for(int U=a[i]&S;U;U&=U-1)w[__builtin_ctz((unsigned)U)]+=v[i];v.swap(w);}__int128 z=accumulate(v.begin(),v.end(),(__int128)0);out+=((n-__builtin_popcount((unsigned)S))%2?-z:z);}need(out>=0&&out<=9223372036854775807LL);return (I)out;}
Poly oddcycle_profile(Graph a){int n=a.size(),N=1<<n;vector<I>cy(N);for(int root=0;root<n;root++){function<void(int,int,int)>go=[&](int last,int mask,int len){if(len>=3&&len%2&&(a[last]>>root&1))cy[mask]++;int U=a[last]&~mask&~((1<<(root+1))-1);for(;U;U&=U-1){int v=__builtin_ctz((unsigned)U);go(v,mask|(1<<v),len+1);}};go(root,1<<root,1);}vector<Poly>dp(N,Poly(n+1));dp[0][0]=1;for(int S=1;S<N;S++){int bit=S&-S;for(int j=0;j<n;j++)dp[S][j+1]+=dp[S^bit][j];for(int C=S;C;C=(C-1)&S)if((C&bit)&&cy[C])for(int j=0;j<n;j++)dp[S][j+1]+=2*cy[C]*dp[S^C][j];}return dp.back();}
I eval(Poly p,I x){I z=0;for(int j=(int)p.size()-1;j>=0;j--)z=z*x+p[j];return z;}
vector<I>newton(Poly p){vector<I>v;for(int x=0;x<(int)p.size();x++)v.push_back(eval(p,x));vector<I>f;while(v.size()){f.push_back(v[0]);vector<I>w;for(int j=1;j<(int)v.size();j++)w.push_back(v[j]-v[j-1]);v.swap(w);}return f;}
I G(Poly p){auto f=newton(p);I ans=0,fac=1;for(int r=1;r<(int)f.size();r++){if(r>1)fac*=r-1;ans=gcd(ans,fac*f[r]);}return ans;}
int main(){int profiles=0,outputs=0;for(int n=1;n<=5;n++)for(int code=0;code<(1<<(n*(n-1)/2));code++){Graph a=tournament(n,code);auto p=oddcycle_profile(a);int d=(n+1)/2;I g=G(p),bank=0,sbank=0;Graph t=a;for(int j=0;j<=max(1,d-1);j++){I h=walks(t);if(j<d)bank=gcd(bank,h);if(j){need(strong(t));sbank=gcd(sbank,h);outputs++;}t=W(t);}auto e=E(a);I he=walks(e);need(strong(e));sbank=gcd(sbank,he);outputs++;need(bank==g&&sbank==g);need(3*walks(W(a))-he==4*walks(a));profiles++;}
Graph v={382,420,88,706,810,972,18,597,204,327};Poly p=oddcycle_profile(v);need(p==Poly({0,0,2772,0,4020,0,1000,0,64,0,1}));Graph w=W(v),z=cyclic(tr(2),tr(2),w);I hw=walks(w),hz=walks(z);need(hw==164673&&hz==196004745&&hz%81==54&&strong(w)&&strong(z));need(walks(W({0}))==3&&walks(W(W({0})))==15&&walks(E({0}))==5);cout<<"{\"status\":\"PASS\",\"odd_cycle_profiles\":"<<profiles<<",\"strong_walk_counts\":"<<outputs<<",\"large_wrapper\":"<<hw<<",\"double_cyclic\":"<<hz<<",\"double_cyclic_mod81\":"<<hz%81<<"}\n";}
