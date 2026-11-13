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
void run(vector<int> Q,vector<int> X,vector<int> Y){
 auto A=fs(inflate(Q,X)),B=fs(inflate(Q,Y));
 for(int j=0;j<(int)A.size();j++)cout<<j+1<<" "<<dec(A[j])<<" "<<dec(B[j])<<" "<<dec(A[j]-B[j])<<"\n";
}

void req(bool v){if(!v)throw runtime_error("replay mismatch");}
vector<i128> packing(const vector<int>& A){
 int n=A.size(),N=1<<n;
 vector<array<i128,16>> D(N),P(N);
 vector<i128> C(N);
 for(int M=1;M<N;M++){
  int root=__builtin_ctz(M);
  if(!(M&(M-1)))D[M][root]=1;
  for(int i=0;i<n;i++)if(D[M][i]){
   if(A[i]>>root&1) C[M]+=D[M][i];
   int b=A[i]&~M&~((1<<root)-1);
   while(b){int e=b&-b;b-=e;D[M|e][__builtin_ctz(e)]+=D[M][i];}
  }
 }
 P[0][0]=1;
 for(int M=1;M<N;M++){
  int root=M&-M, rest=M^root;
  for(int k=0;k<n;k++)P[M][k+1]+=P[rest][k];
  for(int S=rest;S;S=(S-1)&rest)if(!(__builtin_popcount(S)&1)&&C[S|root])
   for(int k=0;k<n;k++)P[M][k+1]+=2*C[S|root]*P[rest^S][k];
 }
 return vector<i128>(P[N-1].begin(),P[N-1].begin()+n+1);
}
i128 powev(vector<i128> p,int z){
 i128 ans=0;for(auto i=p.rbegin();i!=p.rend();++i)ans=ans*z+*i;
 return ans;
}
i128 count_slots(vector<i128> f,int z){
 i128 v=0,choose=1;
 for(int j=1;j<=z&&j<=(int)f.size();j++){
  choose=choose*(z-j+1)/j;v+=choose*f[j-1];
 }
 return v;
}
void audit(vector<int> Q,vector<int> X,vector<int> Y,vector<int> r,i128 h,i128 d4){
 auto A=inflate(Q,X),B=inflate(Q,Y);auto a=fs(A),b=fs(B);
 auto p=packing(A),q=packing(B);
 req(a[0]==h&&b[0]==h&&b[3]-a[3]==d4);
 int n=a.size();for(int j=1;j<n;j++)req(j<n-2?b[j]>a[j]:b[j]==a[j]);
 for(int z=0;z<=n;z++){
  i128 R=0;for(auto i=r.rbegin();i!=r.rend();++i)R=R*z*z+*i;
  req(powev(p,z)==count_slots(a,z)&&powev(q,z)==count_slots(b,z));
  req(powev(q,z)-powev(p,z)==2*z*(z*z-1)*R);
 }
 cout<<"order "<<n<<", H="<<dec(h)<<", F4 gap "<<dec(d4)<<"; two independent profiles agree\n";
}
int main(){
 audit({24,9,3,4,14},{388,449,202,195,335,351,129,48,204},{284,33,506,418,138,17,59,355,114},{240,708,340,38,1},532623,26045760);
 audit({44,33,2,6,15,28},{766,248,1018,832,264,280,944,312,3,434},{32,909,1,821,295,454,415,541,133,373},{720,7016,4752,828,48,1},4591771,808292160);
}
