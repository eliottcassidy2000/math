#include <bits/stdc++.h>
using namespace std;using L=long long;
vector<int> inflate(vector<int> R,vector<int> U) {
 int m=U.size(), n=m+R.size()-1;vector<int>A(n);
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j){
 int a=i<m?0:i-m+1,b=j<m?0:j-m+1;
 if(a==b?(U[i]>>j&1):(R[a]>>b&1))A[i]|=1<<j;
 }return A;
}
vector<pair<L,L>> calc(vector<int>A){
 int n=A.size(),N=(1<<n)-1; vector<array<L,16>> F(N+1),G(N+1),K(N+1);vector<L> H(N+1);
 H[0]=1;for(int i=0;i<n;i++)F[1<<i][i]=G[1<<i][i]=1;
 for(int s=1;s<=N;s++)for(int i=0;i<n;i++){
  if(auto d=F[s][i]){H[s]+=d;for(int R=A[i]&~s;R;R&=R-1)F[s|(R&-R)][__builtin_ctz(R)]+=d;}
  if(auto d=G[s][i])for(int R=(N^A[i]^(1<<i))&~s;R;R&=R-1)G[s|(R&-R)][__builtin_ctz(R)]+=d;
 }
 // K[V][j] counts two optional slots over V with first nonempty start j.
 for(int s=1;s<=N;s++)for(int t=s;t;t=(t-1)&s) {
  L w=H[s^t];for(int r=t;r;r&=r-1){int j=__builtin_ctz(r);K[s][j]+=w*G[t][j];}
 }
 L x[2][16]={}, y[2][16]={}, b[2][16][16]={};L base=0;
 for(int s=0;s<=N;s++){
  base=H[N];
  if(s==0)continue;
  for(int r=s;r;r&=r-1){int i=__builtin_ctz(r);x[1][i]+=H[N^s]*G[s][i];y[1][i]+=H[N^s]*F[s][i];
   for(int z=N^s;z;z&=z-1){int j=__builtin_ctz(z);b[0][i][j]+=F[s][i]*G[N^s][j];b[1][i][j]+=F[s][i]*K[N^s][j];}
  }
 }
 for(int i=0;i<n;i++){x[0][i]=G[N][i];y[0][i]=F[N][i];}
 vector<pair<L,L>> out(N+1);
 for(int s=0;s<=N;s++){
  L q[2]={};
  for(int t=0;t<2;t++)for(int i=0;i<n;i++){
   if(s>>i&1){q[t]+=x[t][i];continue;}
   q[t]+=y[t][i];for(int r=s;r;r&=r-1)q[t]+=b[t][i][__builtin_ctz(r)];
  }
  out[s]={q[0],2*(q[1]+base)};
 }
 return out;
}
int main(int argc,char**argv){
 vector<int>C,X,Y;
 if(argc>1){C={44,33,2,6,15,28};X={766,248,1018,832,264,280,944,312,3,434};Y={32,909,1,821,295,454,415,541,133,373};}
 else{C={24,9,3,4,14};X={388,449,202,195,335,351,129,48,204};Y={284,33,506,418,138,17,59,355,114};}
 auto A=inflate(C,X),B=inflate(C,Y);auto a=calc(A),b=calc(B);
 map<pair<L,L>,int> D;map<L,vector<int>> S; for(int i=1;i<(int)b.size()-1;i++){D[b[i]]=i;S[b[i].first].push_back(i);}
 cout<<"base "<<a[0].first<<" "<<a[0].second<<" "<<b[0].first<<" "<<b[0].second<<endl;
 L best=LLONG_MAX;int x=0,y=0,same=0,neutralA=0,neutralB=0,matched=0,exact=0,positive=0,negative=0; L dmin=LLONG_MAX,dmax=LLONG_MIN; int mincut=-1,maxcut=-1;
 for(int i=1;i<(int)a.size()-1;i++){
   L delta=a[i].first-b[i].first;if(delta>0)positive++;if(delta<0)negative++; if(delta<dmin){dmin=delta;mincut=i;} if(delta>dmax){dmax=delta;maxcut=i;}
   if(a[i].first==a[0].first)neutralA++;
   if(b[i].first==b[0].first)neutralB++;
   if(a[i].first==b[i].first){same++;cout<<"same cut "<<i<<" "<<a[i].first<<" "<<a[i].second-b[i].second<<endl;}
   if(auto it=D.find(a[i]);it!=D.end()){exact++;cout<<"both "<<i<<" "<<it->second<<" "<<a[i].first<<" "<<a[i].second<<endl;}
   if(auto it=S.find(a[i].first);it!=S.end())for(int j:it->second){
    matched++;L d=llabs(a[i].second-b[j].second);if(d<best){best=d;x=i;y=j;}
   }
 }
 cout<<"same mask "<<same<<" H match "<<matched<<" exact "<<exact<<" neutral "<<neutralA<<" "<<neutralB<<endl;
 cout<<"signs "<<positive<<" "<<negative<<" range "<<mincut<<" "<<dmin<<" "<<maxcut<<" "<<dmax<<endl;
 if(x)cout<<"closest "<<x<<" "<<y<<" "<<a[x].first<<" "<<a[x].second<<" "<<b[y].second<<endl;
if(argc>1){
  if(same!=3||matched!=25||exact||neutralA||neutralB||positive!=15942||negative!=16821||best!=2735736) return 3;
 }else{
  if(same!=1||matched!=51||exact||neutralA||neutralB||positive!=3960||negative!=4229||best!=16080) return 4;
 }
}
