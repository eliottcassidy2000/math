#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <set>
#include <string>
using namespace std;
using U=unsigned;
using W=long long;
constexpr int N=11, S=1<<N;
W endsDP[S][N], starts[S][N], into[S][N], outof[S][N];
set<int> missing={13,41,43,55,67,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647};
bool triangle_only=false;
long long countByOrder[N+1]{};
void test(const string& word){
 int n=word.size()+1,root=0; U out[N]{},in[N]{},tree[N]{};
 for(int i=1;i<n;i++){
  char c=word[i-1];
  for(int j=0;j<i;j++)if(j!=root||c=='P')out[j]|=1u<<i;else out[i]|=1u<<j;
  tree[root]|=1u<<i;tree[i]|=1u<<root;
  if(c=='M')root=i;
 }
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
  if(out[i]>>j&1)in[j]|=1u<<i;
  else in[i]|=1u<<j;
 }
 int all=(1<<n)-1;
 for(int x=0;x<n;x++)into[0][x]=outof[0][x]=1;
 for(U A=1;A<=(U)all;A++){
  unsigned b=A;
  for(int j=0;j<n;j++)endsDP[A][j]=starts[A][j]=0;
  while(b){
   int j=__builtin_ctz(b);b&=b-1;
   U Q=A^(1u<<j);
   if(!Q){endsDP[A][j]=starts[A][j]=1;continue;}
   for(U z=Q&in[j];z;z&=z-1)endsDP[A][j]+=endsDP[Q][__builtin_ctz(z)];
   for(U z=Q&out[j];z;z&=z-1)starts[A][j]+=starts[Q][__builtin_ctz(z)];
  }
  for(int j=0;j<n;j++){
   into[A][j]=outof[A][j]=0;
   for(U z=A&in[j];z;z&=z-1)into[A][j]+=endsDP[A][__builtin_ctz(z)];
   for(U z=A&out[j];z;z&=z-1)outof[A][j]+=starts[A][__builtin_ctz(z)];
  }
 }
 W h=0;for(int j=0;j<n;j++)h+=endsDP[all][j];
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(!(tree[i]>>j&1)&&(!triangle_only||(tree[i]&tree[j]))){
  assert((out[i]>>j&1)==1);
  U left=all^((1u<<i)|(1u<<j));
  W hh=h;
  for(U A=left;;A=(A-1)&left){
   hh+=into[A][j]*outof[left^A][i]-into[A][i]*outof[left^A][j];
   if(!A)break;
  }
  countByOrder[n]++;
  if(missing.erase(hh))cout<<hh<<" "<<word<<" "<<i<<" "<<j<<" "<<h<<endl;
 }
}
int main(int argc,char**argv){
 assert(argc==1||(argc==2&&string(argv[1])=="--triangle"));
 triangle_only=argc==2;
 for(int n=3;n<=N&&!missing.empty();n++){
  function<void(string)> gen=[&](string w){
   if(n==11&&missing.empty())return;
   if((int)w.size()==n-1){test(w);return;}
   for(char c:string("PKM")){
    if(!w.empty()&&w.back()=='M'&&c=='K')continue;
    gen(w+c);
   }
  };
  gen("");
  cerr<<n<<" : "<<countByOrder[n]<<" raw chord tests; "<<missing.size()<<" remaining"<<endl;
 }
 for(int h:missing)cerr<<"Missing "<<h<<endl;
 return !missing.empty();
}
