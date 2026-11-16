#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;using ll=long long;
int main(int argc,char**argv){ if(argc!=2)return 2; const int gg=stoi(argv[1]);
 string line;
 while(getline(cin,line)){
  istringstream z(line);string tag;if(!(z>>tag)||tag!="triple")continue;
  int s,m;ll wantH,wantE,wantO;z>>s>>wantH>>wantE>>wantO>>tag>>m>>tag;
  vector<int> a;int k;while(z>>k)a.push_back(k);
  int n=a.size(),r=n-1,F=(1<<n)-1;
  vector<int> P((F+1)*n),H(F+1);H[0]=1;
  for(int i=0;i<n;i++){
   if(a[i]>>i&1)throw runtime_error("loop");
   for(int j=0;j<i;j++)if(((a[i]>>j&1)+(a[j]>>i&1))!=1)throw runtime_error("not tournament");
  }
  ll E=0,O=0;
  for(int S=1;S<=F;S++)for(int p=S;p;p&=p-1){int i=__builtin_ctz(p),R=S^(1<<i),d=!R;
   for(int b=R&a[i];b;b&=b-1)d+=P[R*n+__builtin_ctz(b)];
   P[S*n+i]=d;H[S]+=d;
  }
  for(int S=1<<r;S<=F;S++){
   ll v=ll(P[S*n+r])*H[F^S];if(__builtin_popcount(S)%2)O+=v;else E+=v;
  }
  if(H[F]*m!=wantH||m*E!=wantE||m*O!=wantO||m*H[F^(1<<r)]!=gg)throw runtime_error("profile mismatch");
  cout<<s<<" "<<wantH<<" "<<wantE<<" "<<wantO<<" rows "<<n<<"\n";
 }
}
