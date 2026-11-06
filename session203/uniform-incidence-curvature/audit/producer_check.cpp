#include<bits/stdc++.h>
using namespace std;using L=long long;using U=unsigned;
L H(vector<U>a){int n=a.size(),N=1<<n;vector<L>d(N*n);for(int v=0;v<n;v++)d[(1<<v)*n+v]=1;for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(d[S*n+v]){U z=a[v]&~U(S);while(z){U b=z&-z;z-=b;d[(S|b)*n+__builtin_ctz(b)]+=d[S*n+v];}}L h=0;for(int v=0;v<n;v++)h+=d[(N-1)*n+v];return n?h:1;}
vector<U>make(vector<U>a,int k,U polarity,U uniform,U word,int m){int n=a.size();a.resize(n+m);for(int j=0;j<m;j++){int v=n+j;for(int i=0;i<n;i++){bool out=i<k?(((word>>j)&1)^((polarity>>i)&1)):((uniform>>(i-k))&1);if(out)a[i]|=1u<<v;else a[v]|=1u<<i;}for(int i=0;i<j;i++)a[n+i]|=1u<<v;}return a;}
vector<L>skeleton(vector<U>a,int k,U uniform){int n=a.size()-k;vector<U>w(n);for(int i=0;i<n;i++)w[i]=a[k+i]>>k;vector<int>p(n);iota(p.begin(),p.end(),0);vector<L>total(n+2);do{vector<L>q={1};for(int g=0;g<=n;g++){int l=g?p[g-1]:-1,r=g<n?p[g]:-1;bool e=l<0||r<0||(w[l]>>r&1);bool f=(l<0||(uniform>>l&1))&&(r<0||!(uniform>>r&1));vector<L>z(q.size()+1);for(int i=0;i<(int)q.size();i++){z[i]+=e*q[i];z[i+1]+=f*q[i];}q=z;}for(int i=0;i<(int)q.size();i++)total[i]+=q[i];}while(next_permutation(p.begin(),p.end()));return total;}
L fact(int k){L z=1;while(k)z*=k--;return z;}
L binom(int n,int r){if(r<0||r>n)return 0;L z=1;for(int i=1;i<=r;i++)z=z*(n-i+1)/i;return z;}
int valuation(L x){if(x<0)x=-x;int r=0;while(x&&x%2==0){x/=2;r++;}return r;}
int main(){mt19937 gen(2032026);int cases=0,graphs=0;for(int k=1;k<=4;k++)for(int extra=0;extra<=(k==4?1:3);extra++)for(int rep=0;rep<(k==4?1:3);rep++){
 int n=k+extra,m=2*k;vector<U>a(n);for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(gen()&1)a[i]|=1u<<j;else a[j]|=1u<<i;U pol=gen()&((1u<<k)-1),uniform=gen()&((1u<<extra)-1);L coeff=0;
 for(U bits=0;bits<(1u<<m);bits++){L h=H(make(a,k,pol,uniform,bits,m));coeff+=__builtin_popcount(bits)%2?-h:h;graphs++;}
 auto ar=skeleton(a,k,uniform);L odd=0;for(int r=1;r<=k;r++)if(r<(int)ar.size())odd+=ar[r]*binom(k-1,r-1)*(1LL<<(r-1));L constant=fact(k)*fact(2*k)/(1LL<<(k-1)),expected=(k%2?-1:1)*constant*odd;
 if(coeff!=expected||odd%2!=1||valuation(coeff)!=2*valuation(fact(k))+1){cerr<<"FAIL "<<k<<" "<<extra<<" "<<coeff<<" "<<expected<<" "<<odd<<endl;return 1;}
 cout<<"CASE k="<<k<<" uniform="<<extra<<" polarity="<<pol<<" uniform_mask="<<uniform<<" coefficient="<<coeff<<" constant="<<constant<<" odd="<<odd<<" valuation="<<valuation(coeff)<<" rows=";for(U v:a)cout<<v<<",";cout<<" a_r=";for(L v:ar)cout<<v<<",";cout<<endl;cases++;
 }
 cout<<"PASS cases="<<cases<<" literal_graphs="<<graphs<<endl;
}
