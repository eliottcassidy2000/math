#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
#include <fstream>
#include <string>
using namespace std;
bool prime(int p){if(p<2)return false;for(int d=2;d*d<=p;d++)if(p%d==0)return false;return true;}
int main(int argc,char**argv){
 const int X=4620;
 map<pair<int,int>,vector<int>> bank;
 for(int n=1;n<=7;n++){
 vector<pair<int,int>> E;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)E.push_back({i,j});
 int L=1<<n;vector<char> seen(1<<E.size());int classes=0;
 for(int code=0;code<(int)seen.size();code++)if(!seen[code]){
  vector<int>A(n);
  for(int k=0;k<(int)E.size();k++){auto [i,j]=E[k];if(code>>k&1)A[i]|=1<<j;else A[j]|=1<<i;}
  vector<int> q(n);iota(q.begin(),q.end(),0);
  do{int z=0;for(int k=0;k<(int)E.size();k++){auto [i,j]=E[k];z|=((A[q[i]]>>q[j])&1)<<k;}seen[z]=1;}while(next_permutation(q.begin(),q.end()));
  classes++;
  int D[128][7]{};for(int i=0;i<n;i++)D[1<<i][i]=1;
  for(int s=1;s<L;s++)for(int i=0;i<n;i++)if(D[s][i])for(int w=A[i]&~s;w;w&=w-1){int b=w&-w;D[s|b][__builtin_ctz((unsigned)b)]+=D[s][i];}
  vector<int> h(L);h[0]=1;for(int s=1;s<L;s++)h[s]=accumulate(D[s],D[s]+n,0);
  int P=0;for(int s=0;s<L;s++)P+=h[s]*h[(L-1)^s];
  const int a[]={0,13,25,9,27,5}, b[]={0,18,284,79,127,34};
 for(int rail=1;rail<=5;rail++){
  int N=a[rail]*h.back()+b[rail]*P;
  if(N<=X&&prime(N)){vector<int> witness={n,h.back(),P};witness.insert(witness.end(),A.begin(),A.end());bank.emplace(make_pair(rail,N),witness);}
  }
 }
 cout<<"classes "<<n<<" "<<classes<<"\n";
 }
 int pcount=0;for(int i=3;i<=X;i+=2)pcount+=prime(i);
 cout<<"odd primes "<<pcount<<" hit "<<bank.size()<<"\n";
 for(auto& [p,v]:bank){cout<<p.first<<" "<<p.second;for(int a:v)cout<<" "<<a;cout<<"\n";}
 if(argc>1){
 ofstream out(argv[1]);out<<"prime,template_ids\n";
 for(int p=5;p<=X;p+=2)if(p!=7&&prime(p)){
  string hits;for(int id=1;id<=5;id++)if(bank.count({id,p}))hits+=char('0'+id);
  out<<p<<","<<(hits.empty()?"none_in_this_bank":hits)<<"\n";
 }
 }

}
