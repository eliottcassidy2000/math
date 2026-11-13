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

vector<int> add(vector<int> R,int u){
 int n=R.size();for(int i=0;i<n;i++)if(!(u>>i&1))R[i]|=1<<n;
 R.push_back(u);return R;
}
void req(bool b){if(!b)throw runtime_error("profile mismatch");}
void check(vector<int> Q, vector<int> X, vector<int> Y, int cut, vector<int> coeff, int H){
 auto a=fs(inflate(add(Q,cut),X)),b=fs(inflate(add(Q,cut),Y));
 int n=a.size();req(a[0]==H&&b[0]==H);
 for(int j=1;j<n;j++)req(j<n-2?b[j]>a[j]:b[j]==a[j]);
 for(int z=0;z<=n;z++){
  i128 comb=1, v=0,R=0;
  for(int j=1;j<=z;j++){comb=comb*(z-j+1)/j;v+=comb*(b[j-1]-a[j-1]);}
  for(auto p=coeff.rbegin();p!=coeff.rend();++p)R=R*z*z+*p;
  req(v==2*z*z*(z*z-1)*R);
 }
 cout<<"n="<<n<<" quotient cut="<<cut<<" H="<<H<<": strict F2 through F"<<n-2<<endl;
}
int main(){
 vector<int> Q{24,9,3,4,14},X{388,449,202,195,335,351,129,48,204},Y{284,33,506,418,138,17,59,355,114};
 check(Q,X,Y,27,{916,1292,472,42,1},1068387);
 Q={44,33,2,6,15,28};X={766,248,1018,832,264,280,944,312,3,434};Y={32,909,1,821,295,454,415,541,133,373};
 check(Q,X,Y,4,{3648,11906,6054,940,50,1},7655505);
 check(Q,X,Y,10,{8784,21466,9818,1398,68,1},13775313);
 check(Q,X,Y,12,{4368,14722,7012,1034,52,1},9189287);
}
