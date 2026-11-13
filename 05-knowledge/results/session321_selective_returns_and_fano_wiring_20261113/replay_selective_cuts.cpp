#include<bits/stdc++.h>
using namespace std;using i128=__int128_t;
string dec(i128 x){
 if(x<0)return "-"+dec(-x);if(!x)return "0";string s;while(x){s+=char('0'+x%10);x/=10;}reverse(s.begin(),s.end());return s;
}
vector<int> inflate(vector<int> Q,vector<int> U){
 int n=U.size(),m=Q.size(),N=n+m-1; vector<int> A(N);
 for(int i=0;i<N;i++)for(int j=0;j<N;j++)if(i!=j){
  bool w;
  if(i<n&&j<n)w=U[i]>>j&1;
  else{int a=i<n?0:i-n+1,b=j<n?0:j-n+1;w=Q[a]>>b&1;}
  if(w)A[i]|=1<<j;
 }
 return A;
}
vector<i128> fs(vector<int> A){
 int n=A.size(),N=1<<n;vector<i128> h(N),x(N),y(N),F;
 vector<array<i128,16>> D(N);
 for(int i=0;i<n;i++)D[1<<i][i]=1;
 for(int M=1;M<N;M++)for(int i=0;i<n;i++)if(D[M][i]){
  h[M]+=D[M][i];int b=A[i]&~M;
  while(b){int e=b&-b;b-=e;D[M|e][__builtin_ctz(e)]+=D[M][i];}
 }
 x[0]=1;
 for(int j=1;j<=n;j++){
  fill(y.begin(),y.end(),0);
  for(int M=1;M<N;M++) if(__builtin_popcount(M)>=j) for(int S=M;S;S=(S-1)&M) y[M]+=h[S]*x[M^S];
  x.swap(y); F.push_back(x[N-1]);
 }
 return F;
}

vector<array<i128,16>> ends(const vector<int>& A){
 int n=A.size(),N=1<<n;vector<array<i128,16>> D(N);
 for(int i=0;i<n;i++)D[1<<i][i]=1;
 for(int M=1;M<N;M++)for(int i=0;i<n;i++)if(D[M][i]){
  int b=A[i]&~M;
  while(b){int e=b&-b;b-=e;D[M|e][__builtin_ctz(e)]+=D[M][i];}
 }
 return D;
}
vector<i128> cuts(vector<int> A){
 int n=A.size(),N=1<<n,M=N-1;
 auto B=A;for(int i=0;i<n;i++)B[i]=M^A[i]^(1<<i);
 auto E=ends(A),S=ends(B);
 i128 C[16][16]{},h=0,best=-1;int num=0;vector<int> examples;vector<i128> out(N);
 for(auto a:E[M])h+=a;
 for(int X=1;X<M;X++){
  int Y=M^X;
  for(int i=0;i<n;i++)if(E[X][i])
   for(int j=0;j<n;j++)C[i][j]+=E[X][i]*S[Y][j];
 }
 for(int U=1;U<M;U++){
  i128 k=0;for(int i=0;i<n;i++)if(U>>i&1)k+=S[M][i];else {
   k+=E[M][i];for(int j=0;j<n;j++)if(U>>j&1)k+=C[i][j];
  }
  out[U]=k;k-=h;
  if(best<0||k<best){best=k;num=0;examples.clear();}
  if(best==k){num++;if(examples.size()<5)examples.push_back(U);}
 }
 cout<<"h="<<dec(h)<<" min "<<dec(best)<<" cuts "<<num;
 for(int a:examples)cout<<" "<<a;cout<<" start11 "<<dec(S[M][11])<<endl; if(num!=1 || examples!=vector<int>{2048}) throw runtime_error("minimum cut"); return out;
}

void cmp(vector<int> R,vector<int> X,vector<int> Y){
 auto a=cuts(inflate(R,X)),b=cuts(inflate(R,Y));vector<pair<int,i128>> matched;
 for(int U=1;U<(int)a.size()-1;U++)if(a[U]==b[U])matched.push_back({U,a[U]});
 vector<pair<int,i128>> want=R.size()==5?vector<pair<int,i128>>{{7167,1068387}}:vector<pair<int,i128>>{{2048,7655505},{5120,13775313},{6144,9189287}};
 if(matched!=want)throw runtime_error("equalizing cuts");
 for(auto [u,v]:matched)cout<<u<<" "<<dec(v)<<endl;
 cout<<"number same cuts "<<matched.size()<<endl;
}
int main(){
 cmp({24,9,3,4,14},{388,449,202,195,335,351,129,48,204},{284,33,506,418,138,17,59,355,114});
 cmp({44,33,2,6,15,28},{766,248,1018,832,264,280,944,312,3,434},{32,909,1,821,295,454,415,541,133,373});
}
