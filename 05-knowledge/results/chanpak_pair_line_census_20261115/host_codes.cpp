#include <bits/stdc++.h>
using namespace std; using ul=uint64_t;
ul canon(const array<int,10>& a,int n) {
  array<pair<ul,int>,10> v;
  for(int i=0;i<n;i++){
    int d=__builtin_popcount(a[i]);ul z=ul(d)<<40;
    for(int m=a[i];m;m&=m-1) z+=1ULL<<(4*__builtin_popcount(a[__builtin_ctz(m)]));
    v[i]={z,i};
  }
  sort(v.begin(),v.begin()+n);
  ul best=ULLONG_MAX;
  function<void(int)> calc=[&](int l) {
    if(l==n) {ul x=0;int p=0;
      for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) x|=ul((a[v[i].second]>>v[j].second)&1)<<p++;
      best=min(best,x);return;
    }
    int r=l+1;while(r<n && v[r].first==v[l].first)++r;
    do {calc(r);} while(next_permutation(v.begin()+l,v.begin()+r));
  };
  calc(0);return best;
}
array<int,10> decode(ul x,int n){
  array<int,10>a={};int p=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if((x>>p++)&1)a[i]|=1<<j;else a[j]|=1<<i;
  return a;
}
pair<ul,ul> profile(const array<int,10>& a,int n){
  int full=(1<<n)-1;ul dp[1024][10]={},h[1024]={};h[0]=1;
  for(int i=0;i<n;i++)dp[1<<i][i]=1;
  for(int s=1;s<=full;s++)for(int w=s;w;w&=w-1){
    int i=__builtin_ctz(w); ul t=dp[s][i]; h[s]+=t;
    if(t)for(int m=a[i]&~s;m;m&=m-1){int j=__builtin_ctz(m);dp[s|1<<j][j]+=t;}
  }
  ul p=0;for(int s=1;s<full;s++)p+=h[s]*h[full^s];
  return {h[full],p};
}
bool strong(const array<int,10>& a,int n){
 int full=(1<<n)-1;
 for(int v=0;v<n;v++){
   int reached=1<<v,pending=reached;
   while(pending){int u=__builtin_ctz(pending);pending&=pending-1;
   int fresh=a[u]&~reached;pending|=fresh;reached|=fresh;}
   if(reached!=full)return false;
 }
 return true;
}
int main(int argc,char**argv){
 int limit=argc>1?atoi(argv[1]):9;
 if(limit<1||limit>10)return 2;
 unordered_set<ul> olds={0};
 int expected[]={0,1,1,2,4,12,56,456,6880,191536,9733056};
 for(int n=1;n<=limit;n++){
   if(n>1){unordered_set<ul> next;
     vector<ul> vv(olds.begin(),olds.end());
     #pragma omp parallel
     {
     unordered_set<ul> local;
     #pragma omp for schedule(dynamic,512)
     for(size_t ix=0;ix<vv.size();++ix){ul x=vv[ix];auto base=decode(x,n-1);
       for(int m=0;m<1<<(n-1);m++){
         auto a=base;a[n-1]=m;
         for(int i=0;i<n-1;i++)if(!(m>>i&1))a[i]|=1<<(n-1);
         if(n>=9) {
           // Extension is covered by deleting ANY minimum outdegree vertex.
           int d=__builtin_popcount((unsigned)m);
           bool skip=false;
           for(int i=0;i<n-1;++i) if(__builtin_popcount((unsigned)a[i])<d){skip=true;break;}
           if(skip) continue;
         }
         local.insert(canon(a,n));
       }
     }
     #pragma omp critical
     next.merge(local);
     }
     olds.swap(next);
   }
   assert(olds.size()==expected[n]); cerr<<"Generated "<<n<<" "<<olds.size()<<endl;
   for (ul x:olds) cout << n <<" " << x << "\n";
 }
}
