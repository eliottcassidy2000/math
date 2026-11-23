// Exact finite Q36 decoder test. Birth status through 36 is inherited from
// COLOR_WORD_BOUNDARY.md and COLOR_FIRST_FAILURE.md; this test adds the sign bit.
#include <bits/stdc++.h>
using namespace std;
using U=uint64_t;
const int N=36;
U A[N]; int C[N],p[N]; 
map<string,vector<array<int,N>>> fibres;
bool born(int n){return n==1||n==15||n==16||n==17||n==23||n>=25;}
void dfs(U R,int i,int v,int first){
 if(!R){
  if(v<=first)return;
  string k; k+=char(first); k+=char(v);
  for(int t=1;t<N;t++)k+=char(C[max(p[t],p[t-1])]);
  array<int,N> x; copy(p,p+N,x.begin());
  fibres[k].push_back(x); return;
 }
 auto S=R|(1ULL<<v);int ones=0;
 for(U Y=R;Y;Y&=Y-1){
  int j=__builtin_ctzll(Y);
  U Z=A[j]&S;
  if(!Z)return;
  if(!(Z&(Z-1))){if(++ones>1 || (Z==(1ULL<<v) && (R&(R-1))))return;}
 }
 for(U Y=A[v]&R;Y;Y&=Y-1){
   int j=__builtin_ctzll(Y);p[i]=j;
   dfs(R^(1ULL<<j),i+1,j,first);
 }
}
int main(){
 for(int i=0;i<N;i++){
  int d=0;for(int j=0;j<N;j++) if(i!=j && int(sqrt(i+j+2))*int(sqrt(i+j+2))==i+j+2){
    A[i]|=1ULL<<j; d+=j<i;
  }
  C[i]=born(i+1)?0:d==2?1:2;
 }
 for(int i=0;i<N;i++){p[0]=i;dfs(((1ULL<<N)-1)^(1ULL<<i),1,i,i);}
 int n=0,collisions=0;
 U sign_all=(1ULL<<N)-1, edge_all=(1ULL<<(N-1))-1;
 vector<pair<int,int>> endpoints;
 for(auto &[k,v]:fibres){
  n+=v.size();if(v.size()<2)continue;
  if(v.size()!=2) return 2;
  ++collisions;
  endpoints.emplace_back(int(k[0])+1,int(k[1])+1);
  U signs=0,edges=0;
  array<int,2> totals{};
  for(int j=0;j<N;j++){
    signs |= U(((v[0][j]+1)^(v[1][j]+1))>>1&1)<<j;
    if(j+1<N){
      int x=(v[0][j]+1)&(v[0][j+1]+1)&1;
      int y=(v[1][j]+1)&(v[1][j+1]+1)&1;
      totals[0]+=x; totals[1]+=y; edges|=U(x^y)<<j;
    }
  }
  if(totals[0]!=totals[1]) return 3;
  sign_all &= signs; edge_all &= edges;
  cout<<"position 19 vertices "<<v[0][18]+1<<' '<<v[1][18]+1<<'\n';
  cout<<"endpoints "<<int(k[0])+1<<' '<<int(k[1])+1<<" fibre "<<v.size()<<'\n';
  for(auto p:v){for(int i=0;i<N;i++)cout<<p[i]+1<<(i+1<N?',':'\n');}
 }
 cout<<n<<" paths "<<fibres.size()<<" labels\n";
 if(n!=12960||fibres.size()!=12957||collisions!=3) return 4;
 if(endpoints!=vector<pair<int,int>>{{6,32},{12,21},{21,32}}) return 5;
 if(sign_all!=(1ULL<<18)||edge_all!=((1ULL<<16)|(1ULL<<22))) return 6;
 cout<<"All-pair sign position: 19; all-pair odd--odd edge positions: 17, 23\n";
}
