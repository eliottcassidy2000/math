#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
// Partition by the exact set of repeated children. All words have unweighted
// multiplicity one. Endpoint counts are stored as bitplanes. Complete cover
// factors are applied only to completed words, at precision 8, 4, or 2.
int side(const vector<int>& rows) {
 int answer=0;
 auto run=[&](vector<int> repeat){
  int nr=repeat.size(), ns=13-nr, precision=3-nr;
  vector<int> order, incoming(13), steps(13);
  for(int i=0;i<13;i++){bool used=false;for(int j:repeat)used|=i==j;if(!used)order.push_back(i);}
  for(int i:repeat)order.push_back(i);
  for(int i=0;i<13;i++){
   steps[i]=1<<(i<ns?i:ns+2*(i-ns));
   for(int j=0;j<13;j++)if(rows[order[j]]>>order[i]&1)incoming[i]|=1<<j;
  }
  int states=1<<(ns+2*nr); vector<uint64_t> dp(states);
  for(int s=1;s<states;s++){
   uint64_t out=0;
   for(int i=0;i<13;i++){
    if(i<ns?!(s&steps[i]):!((s/steps[i])&3))continue;
    int prev=s-steps[i];uint64_t d=dp[prev];int value=prev==0;
    for(int b=0;b<precision;b++)value+=__builtin_popcount((d>>(13*b))&incoming[i])<<b;
    for(int b=0;b<precision;b++)if(value>>b&1)out|=1ull<<(i+13*b);
   }
   dp[s]=out;
  }
  int possibilities=1<<nr;
  for(int c=0;c<possibilities;c++){
   int state=(1<<ns)-1;vector<int> visits(9,1);
   for(int j=0;j<nr;j++){int k=2+((c>>j)&1);state+=k*steps[ns+j];visits[repeat[j]]=k;}
   uint64_t d=dp[state];int total=0;
   for(int b=0;b<precision;b++)total+=__builtin_popcount((d>>(13*b))&8191u)<<b;
   answer=(answer+(1<<nr)*total)&7;
  }
 };
 run({});
 for(int i=0;i<9;i++)run({i});
 for(int i=0;i<9;i++)for(int j=i+1;j<9;j++)run({i,j});
 return answer;
}
int main(){int n;cin>>n;for(int z=0;z<n;z++){
 int node,a,c;cin>>node>>a>>c;vector<int> q(9),s(13),t(13);for(int &r:q)cin>>r;
 for(int i=0;i<9;i++){s[i]=q[i]|(((a>>i)&1)<<9)|(((c>>i)&1)<<10)|(1<<11);t[i]=(511^q[i]^(1<<i))|(((a>>i)&1)<<9)|(((c>>i)&1)<<10)|(1<<11);}
 s[9]=t[9]=(511^a)|(1<<10)|(1<<11);s[10]=t[10]=(511^c)|(1<<11);
 s[11]=t[11]=1<<12;s[12]=t[12]=511|(1<<9)|(1<<10);
 auto start=chrono::steady_clock::now();int left=side(s)*3%8,right=side(t)*3%8;
 cout<<node<<'\t'<<a<<'\t'<<c<<'\t'<<left<<'\t'<<right<<'\n';cout.flush();
 if(z%100==0)cerr<<z+1<<'/'<<n<<" pair_seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
}}
