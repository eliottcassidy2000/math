#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
using U=uint64_t;
struct Canon{U lo,hi,merged;int aut;bool sc;};
Canon canonical(const vector<int>&r){
 int n=r.size(); vector<int>p(n);iota(p.begin(),p.end(),0);
 U lo=~U(0),hi=0;int aut=0;
 do{U v=0;for(int q=1;q<n;q++)for(int a=0;a<q;a++)v=(v<<1)|((r[p[a]]>>p[q])&1);
    if(v<lo){lo=v;aut=1;}else if(v==lo)aut++;hi=max(hi,v);
 }while(next_permutation(p.begin(),p.end()));
 U op=((U(1)<<(n*(n-1)/2))-1)^hi;
 return {lo,hi,min(lo,op),aut,lo==op};
}
vector<U> ends(const vector<int>&r){int n=r.size(),N=1<<n;vector<U>d(N*n);
 for(int m=1;m<N;m++)for(int v=0;v<n;v++)if(m>>v&1){int q=m^(1<<v);U z=q?0:1;for(int u=0;u<n;u++)if(q>>u&1&&r[u]>>v&1)z+=d[q*n+u];d[m*n+v]=z;}
 return vector<U>(d.end()-n,d.end());}
int main(){int count;cin>>count;cout<<"[";
 for(int g=0;g<count;g++){int n;cin>>n;vector<int>r(n),op(n);for(int&x:r)cin>>x;
 for(int v=0;v<n;v++)op[v]=((1<<n)-1)^r[v]^(1<<v);
 auto e=ends(r),s=ends(op);auto c=canonical(r);if(g)cout<<",";
 cout<<"{\"canonical\":"<<c.lo<<",\"merged\":"<<c.merged<<",\"aut\":"<<c.aut<<",\"self_converse\":"<<c.sc<<",\"H\":"<<accumulate(e.begin(),e.end(),U(0))<<",\"cards\":[";
 for(int x=0;x<n;x++){vector<int>sub;for(int v=0;v<n;v++)if(v!=x){int row=0,k=0;for(int w=0;w<n;w++)if(w!=x){if(r[v]>>w&1)row|=1<<k;k++;}sub.push_back(row);}auto q=canonical(sub);
 if(x)cout<<",";cout<<"{\"deleted\":"<<x<<",\"merged\":"<<q.merged<<",\"starts\":"<<s[x]<<",\"ends\":"<<e[x]<<"}";
 }cout<<"]}";
 }cout<<"]\n";
}
