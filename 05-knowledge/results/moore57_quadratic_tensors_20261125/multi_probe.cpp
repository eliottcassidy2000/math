#define main supplied_main
#include "slot_dlx.cpp"
#undef main
#include <set>
using Word=array<int,52>;
int main(int argc,char**argv){
 int n;cin>>n; if(n!=53||argc<4)return 2;
 int C[52][56][3];bool ok[52][56]{};
 int old,num,marker[3],y;
 cin>>old>>num>>y>>marker[0]>>marker[1]>>marker[2];
 for(int i=0;i<52;i++){
  cin>>old>>num;
  while(num--){
   cin>>y>>C[i][y][0]>>C[i][y][1]>>C[i][y][2];
   ok[i][y]=true;
   for(int g=0;g<3;g++)if(C[i][y][g]==marker[g])ok[i][y]=false;
  }
 }
 int idx[3][56]{};
 for(int g=0;g<3;g++){set<int>S;for(int i=0;i<52;i++)for(int j=0;j<56;j++)if(ok[i][j])S.insert(C[i][j][g]);if(S.size()!=52)return 4;
  int k=(g+1)*52;for(int j:S)idx[g][j]=k++;
 }
 array<uint64_t,52> base{};
 ifstream f(argv[1]);int g,k;
 while(f>>g>>k)for(auto&t:base){f>>y;t|=1ULL<<y;}
 ifstream p(argv[2]);vector<Word>W;
 while(p>>g>>k>>y){Word w;for(int&z:w)p>>z;W.push_back(w);}
 const int how=stoi(argv[3]);int lim=argc>4?stoi(argv[4]):10000,reps=argc>5?stoi(argv[5]):(int)W.size();
 mt19937 rng(83172);int cap=0;
 for(int t=0;t<reps;t++){
  auto B=base;
  for(int i=0;i<how;i++){
   const auto&w=W[how==1?t:rng()%W.size()];for(int h=0;h<52;h++)B[h]|=1ULL<<w[h];
  }
  DLX d(208,lim);
  for(int h=0;h<52;h++)for(uint64_t v=B[h];v;v&=v-1){
   int y=__builtin_ctzll(v);
   if(!ok[h][y])return 5;
   d.add(56*h+y,{h,idx[0][C[h][y][0]],idx[1][C[h][y][1]],idx[2][C[h][y][2]]});
  }
  if(d.go()){
   Word w;for(int x:d.path)w[x/56]=x%56;
   cout<<"FULL iteration "<<t<<endl<<"7";for(int y:w)cout<<" "<<y;
   cout<<endl;return 0;
  }
  if(d.visits>(uint64_t)lim)++cap;
 }
 cout<<"words per union "<<how<<" trials "<<reps<<" capped "<<cap<<endl;
}
