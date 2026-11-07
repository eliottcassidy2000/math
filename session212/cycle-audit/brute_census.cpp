// Independent permutation oracle. No subset DP or graph canonicalization here.
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 if(argc!=2)throw runtime_error("usage: brute_census classes_directory");
 for(int n=3;n<=8;n++){
  ifstream f(string(argv[1])+"/classes"+to_string(n)+".rows");
  if(!f)throw runtime_error("missing input");
  map<int,long long>hist;vector<int>a(n),p(n-1);long long seen=0;
  while(f>>a[0]){
   for(int i=1;i<n;i++)if(!(f>>a[i]))throw runtime_error("short row");
   for(int i=0;i<n;i++)for(int j=0;j<n;j++)
    if((i==j && (a[i]>>j&1)) || (i!=j && ((a[i]>>j&1)==(a[j]>>i&1))))throw runtime_error("not tournament");
   for(int i=1;i<n;i++)p[i-1]=i;
   int c=0;
   do{
    if(!(a[0]>>p[0]&1) || !(a[p.back()]&1))continue;
    bool ok=true;
    for(int j=1;j<n-1;j++)if(!(a[p[j-1]]>>p[j]&1)){ok=false;break;}
    c+=ok;
   }while(next_permutation(p.begin(),p.end()));
   hist[c]++;seen++;
  }
  cout<<n<<' '<<seen;
  for(auto [c,m]:hist)cout<<' '<<c<<':'<<m;
  cout<<'\n';
 }
}
