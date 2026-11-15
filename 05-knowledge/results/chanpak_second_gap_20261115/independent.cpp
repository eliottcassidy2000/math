#include <bits/stdc++.h>
using namespace std;
struct Chain { int h; vector<int>a,c; };
vector<Chain> listc[16];
int n,all,ncap,h,edge[32];
bool dfs(int u,int remain){
 if(!remain)return ++h==ncap;
 for(int x=edge[u]&remain; x; x&=x-1) {
  if(dfs(__builtin_ctz(x),remain^(x&-x)))return true;
 }
 return false;
}
void add(const vector<int>& gaps,int n) {
 vector<int>a(n),junction{0},middle; int l=0,h=1;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)a[i]|=1<<j;
 for(int s:gaps){int r=l+s-1;for(int i=l+1;i<r;i++)middle.push_back(i);
  a[l]^=1<<r;a[r]^=1<<l;h*=1+(1<<(s-2));junction.push_back(r);l=r;
 }
 reverse(junction.begin(),junction.end());for(int i:middle)junction.push_back(i);
 for(int i=0;i<n;i++)assert(a[junction[i]]&(1<<junction[(i+1)%n]));
 listc[n].push_back({h,a,junction});
}
void gen(int d,vector<int>&v,int n){if(!d){add(v,n);return;}
 for(int j=2;j<=4;j++)if(d>=j){v.push_back(j+1);gen(d-j,v,n);v.pop_back();}}
int f[16]={0,1,1,3,5,9};
int main(){
 vector<int>g;for(int i=6;i<=14;i++)f[i]=f[i-3]*5;
 for(int i=3;i<=14;i++)gen(i-1,g,i);
 for(n=9;n<=14;n++)for(int m:{4,5}){
  int b=n+1-m;
  if(b<6)continue;
  int p=m-1,q=b-1;all=(1<<n)-1;ncap=f[n]+5;int total=0,hit=0,mn=ncap;
  vector<int>words;for(int z=0;z<(1<<(p+q));z++)if(__builtin_popcount((unsigned)z)==p)words.push_back(z);
  for(const Chain&A:listc[m])for(const Chain&B:listc[b])if(A.h==f[m]&&B.h==f[b]&&A.h*B.h==f[n]){
   for(int v=0;v<m;v++)for(int w=0;w<b;w++){
    vector<int>ap(m),bp(b),base(n);
    for(int i=0;i<m;i++)ap[i]=A.c[(v+i)%m];
    for(int i=0;i<b;i++)bp[i]=B.c[(w+i)%b];
    for(int i=0;i<m;i++)for(int j=0;j<m;j++)if(A.a[ap[i]]&(1<<ap[j]))base[i]|=1<<j;
    for(int i=0;i<b;i++)for(int j=0;j<b;j++)if(B.a[bp[i]]&(1<<bp[j]))base[i?p+i:0]|=1<<(j?p+j:0);
    for(int z:words){
     copy(base.begin(),base.end(),edge);
     int i=0,j=0;for(int k=0;k<p+q;k++){
      if(z&(1<<k)){++i;for(int r=j+1;r<=q;r++)edge[i]|=1<<(p+r);}
      else {++j;for(int r=i+1;r<=p;r++)edge[p+j]|=1<<r;}
     }
     h=0;for(int t=0;t<n;t++)if(dfs(t,all^(1<<t)))break;
     hit+=(h==f[n]+4);mn=min(mn,h);total++;
    }
   }
  }
  cout<<"n="<<n<<" a="<<m<<" matrices="<<total<<" hits="<<hit<<"\n"<<flush;
  assert(!hit);
 }
}
