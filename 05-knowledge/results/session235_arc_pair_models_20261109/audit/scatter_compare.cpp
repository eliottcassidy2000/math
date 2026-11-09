// Independent direct scatter: visit all oriented supergraphs of each permutation.
// No zeta transform or endpoint DP. Compare producer binary tables in census/.
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <functional>
#include <stdexcept>
using namespace std;
int main(int argc,char**argv){
 string dir=argc>1?argv[1]:"../census";
 for(int n=1;n<=6;n++){
  int m=n*(n-1)/2,ix[6][6]; vector<int> pw(m+1,1);
  for(int k=0;k<m;k++)pw[k+1]=3*pw[k];
  int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)ix[i][j]=ix[j][i]=k++;
  vector<int> h(pw[m]),p(n);iota(p.begin(),p.end(),0); long long contributions=0;
  do{
   vector<bool> fixed(m); int seed=pw[m]-1;
   for(int t=1;t<n;t++){int u=p[t-1],v=p[t];int q=ix[u][v];fixed[q]=true;seed-=(u<v?2:1)*pw[q];}
   vector<int> free;for(int q=0;q<m;q++)if(!fixed[q])free.push_back(q);
   function<void(int,int)> scatter=[&](int t,int code){if(t==(int)free.size()){h[code]++;contributions++;return;}int z=pw[free[t]];scatter(t+1,code);scatter(t+1,code-z);scatter(t+1,code-2*z);};
   scatter(0,seed);
  }while(next_permutation(p.begin(),p.end()));
  vector<int> expected(h.size());ifstream f(dir+"/counts"+to_string(n)+".bin",ios::binary);
  if(!f.read((char*)expected.data(),expected.size()*sizeof(int)))throw runtime_error("Cannot read binary table");
  for(int c=0;c<pw[m];c++)if(h[c]!=expected[c])throw runtime_error("Count mismatch");
  cout<<"PASS n="<<n<<" all_oriented="<<pw[m]<<" direct_path_contributions="<<contributions<<" H7="<<count(h.begin(),h.end(),7)<<" H21="<<count(h.begin(),h.end(),21)<<"\n";
 }
}
