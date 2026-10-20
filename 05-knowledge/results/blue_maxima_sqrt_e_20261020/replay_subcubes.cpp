#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 if(argc!=4)return 2;int k=stoi(argv[1]);if(k<1||k>6)return 2;int m=k*(k-1)/2,sz=1<<m;
 vector<int> reps;ifstream inrep(argv[2]);int x;while(inrep>>x)reps.push_back(x);
 vector<uint32_t> expected(sz*reps.size());ifstream in(argv[3],ios::binary);in.read((char*)expected.data(),expected.size()*4);if(in.gcount()!=expected.size()*4)return 3;
 int id[7][7],z=0;for(int i=0;i<k;i++)for(int j=i+1;j<k;j++)id[i][j]=id[j][i]=z++;
 uint64_t updates=0,candidates=0;auto begin=chrono::steady_clock::now();
 for(size_t r=0;r<reps.size();r++){
  int a=reps[r];vector<uint32_t> found(sz,0);vector<int> p(k);iota(p.begin(),p.end(),0);
  do{for(int signs=0;signs<(1<<(k-1));signs++){
   bool ok=true;int required=0,selected=0;
   for(int t=0;t<k-1;t++){
    int i=p[t],j=p[t+1],s=(signs>>t)&1,u=(signs>>(t+1))&1,q=id[i][j];
    if(s==u){int arc=((a>>q)&1)^(i>j);if(arc==s){ok=false;break;}}
    else{selected|=1<<q;if(s==0)required|=1<<q;}
   }
   if(!ok)continue;candidates++;int free=(sz-1)^selected,sub=free;
   while(true){found[required|sub]++;updates++;if(sub==0)break;sub=(sub-1)&free;}
  }}while(next_permutation(p.begin(),p.end()));
  for(int b=0;b<sz;b++)if(found[b]!=expected[r*sz+b]){cerr<<"Mismatch "<<k<<" "<<a<<" "<<b<<" "<<found[b]<<" "<<expected[r*sz+b]<<endl;return 4;}
 }
 cout<<"{\"k\":"<<k<<",\"representatives\":"<<reps.size()<<",\"pointwise_checks\":"<<expected.size()<<",\"all_equal\":true,\"compatible_prefixes\":"<<candidates<<",\"subcube_updates\":"<<updates<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-begin).count()<<"}\n";
}
