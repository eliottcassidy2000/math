#include <array>
#include <algorithm>
#include <bit>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;
constexpr int P=13,Z=(1<<P)-1,F=18;
int rot(int m,int k){return ((m<<k)|(m>>(P-k)))&Z;}
int rev(int m){int ret=m&1;for(int i=1;i<P;i++)if(m>>i&1)ret|=1<<(P-i);return ret;}
int c(int a,int b,int k){return P-2*popcount((unsigned)(a^rot(b,k)));}
int v(int m,int i){return m>>i&1?1:-1;}
int mod(int q){return (q%3+3)%3;}
int sig(int m){int a=0;for(int i=1;i<=6;i++)a=3*a+mod(c(m,m,i));return a;}
int canon(int m){int q=m;for(int k=0;k<P;k++){int f=rot(m,k);q=min({q,f,f^Z});}return q;}
bool eq(int a,int b,int r,int s){
 const int first=mod(c(a,rev(b),0)-c(r,rev(s),0));
 for(int k=1;k<P;k++)if(mod(c(a,rev(b),k)-c(r,rev(s),k))!=first)return false;
 return true;
}
bool pass(int w,int v){int d=w+(v<0?1:-2);return d>=10&&d<=46&&d%3==1;}
int profilesP[4][3]={{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1}};
map<array<int,10>,bool> cache;
bool rowbound(const array<int,4>& n,int Mat[3][3]){
 array<int,10> key{};
 for(int i=0;i<4;i++)key[i]=n[i];
 for(int i=0,u=4;i<3;i++)for(int j=i;j<3;j++)key[u++]=Mat[i][j];
 auto it=cache.find(key);if(it!=cache.end())return it->second;
 int tot[16]={},rhs[16]={};
 for(int u=1;u<16;u++)for(int j=0;j<4;j++)if(u>>j&1)rhs[u]+=9*n[j]*(n[j]-1);
 for(int r=0;r<4;r++)if(n[r]){
   int best[16];for(int i=0;i<16;i++)best[i]=-10000;
   for(int a=0;a<=n[0]-(r==0);a++)for(int b=0;b<=n[1]-(r==1);b++)
   for(int c=0;c<=n[2]-(r==2);c++)for(int d=0;d<=n[3]-(r==3);d++){
    int k[4]={a,b,c,d};bool ok=true;
    for(int i=0;i<3;i++){
      int t=0;for(int j=0;j<4;j++)t+=(2*k[j]-n[j])*profilesP[j][i];
      for(int j=0;j<3;j++)t+=profilesP[r][j]*Mat[j][i];
      if(mod(t)!=0){ok=false;break;}
    }
    if(!ok)continue;
    int sc[4];for(int j=0;j<4;j++)sc[j]=2*k[j]*(n[j]-k[j])-(r==j?3*k[j]:0);
    int sum[16]={};for(int u=1;u<16;u++){
     sum[u]=sum[u&(u-1)]+sc[countr_zero((unsigned)u)];
     best[u]=max(best[u],sum[u]);
    }
   }
   if(best[15]==-10000)return cache[key]=false;
   for(int u=1;u<16;u++)tot[u]+=n[r]*best[u];
 }
 for(int u=1;u<16;u++)if(tot[u]<rhs[u])return cache[key]=false;
 return cache[key]=true;
}
void test_three(){
 vector<int> all[729],reps[729];
 for(int m=0;m<=Z;m++){all[sig(m)].push_back(m);if(m==canon(m))reps[sig(m)].push_back(m);}
 long long checked=0,rankok=0,within=0,deltas=0,profiles=0,caps=0,complete=0;
 map<array<int,4>,long long> shape;set<array<int,10>> survivors;
 for(int mask=0;mask<64;mask++){
   int g=0;for(int k=1;k<=6;k++)if(mask>>(k-1)&1)g|=(1<<k)|(1<<(P-k));
   const int key=sig(g);
   for(int f:reps[key])for(int h:reps[key])if(h>=f)for(int e:all[key]){
      checked++;
      if(!eq(rev(f),h,g,e)||!eq(f,e,g,h)||!eq(h,rev(e),g,f))continue;
      rankok++;
      int C[3][3]={{g,f,h},{rev(f),g,e},{rev(h),rev(e),g}};
      bool ok=true;
      for(int i=0;i<3&&ok;i++)for(int k=1;k<P;k++){
         int dot=0;for(int j=0;j<3;j++)dot+=c(C[i][j],C[i][j],k);
         if(!pass((3*P-dot)/2,v(g,k))){ok=false;break;}
      }
      if(!ok)continue;
      within++;
      int opts[3][3]={};
      for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)for(int d=0;d<=F;d++){
         bool yes=true;
         for(int k=0;k<P;k++){
            int dot=0;for(int x=0;x<3;x++)dot+=c(C[i][x],C[j][x],k);
            if(!pass((3*P-dot)/2+d,v(C[i][j],k))){yes=false;break;}
         }
         if(yes)opts[i][j]|=1<<d;
      }
      if(!opts[0][1]||!opts[0][2]||!opts[1][2])continue;
      deltas++;
      for(int a=0;a<=18;a++)for(int b=0;b<=18-a;b++)for(int c=0;c<=18-a-b;c++){
       int d=18-a-b-c;
       if((opts[0][1]>>(c+d)&1)&&(opts[0][2]>>(b+d)&1)&&(opts[1][2]>>(b+c)&1)){
         profiles++;
         if(max({a,b,c,d})>6)continue;
         caps++;
         int M[3][3];for(int i=0;i<3;i++)for(int j=0;j<3;j++)M[i][j]=mod(2*popcount((unsigned)C[i][j])-13);
         if(!rowbound({a,b,c,d},M))continue;
         complete++;
         array<int,10> K{a,b,c,d};for(int m=0,t=4;m<3;m++)for(int n=m;n<3;n++)K[t++]=M[m][n];survivors.insert(K);
         array<int,4> z{a,b,c,d};sort(z.begin(),z.end());shape[z]++;
       }
      }
   }
 }
 cout<<"triples "<<checked<<" nonconstant square-zero "<<rankok<<" within "<<within<<" across "<<deltas<<" profile choices "<<profiles<<" capped "<<caps<<" row-bound "<<complete<<" distinct constraints "<<cache.size()<<endl;
 const set<array<int,10>> expected{{3,5,5,5,2,1,1,2,1,2},{5,3,5,5,2,1,2,2,2,2},{5,5,3,5,2,2,1,2,2,2},{5,5,5,3,2,2,2,2,1,2},{6,6,6,0,0,0,0,0,0,0}};
 assert(survivors==expected&&complete==17&&caps==552);
 for(auto A:survivors){for(int x:A)cout<<x<<' ';cout<<'\n';}
 for(auto [a,n]:shape){for(int j:a)cout<<j<<",";cout<<" n "<<n<<endl;}
}
struct Type {int r,p;array<int,4> k;};
void parity_test(){
 array<int,4> n{3,5,5,5};
 const int M[3][3]={{2,1,1},{1,2,1},{1,1,2}};
 vector<Type>T[4][2];
 for(int r=0;r<4;r++)
 for(int a=0;a<=n[0]-(r==0);a++)for(int b=0;b<=n[1]-(r==1);b++)
 for(int c=0;c<=n[2]-(r==2);c++)for(int d=0;d<=n[3]-(r==3);d++){
   array<int,4> k{a,b,c,d};bool ok=true;
   for(int i=0;i<3;i++){
    int t=0;for(int j=0;j<4;j++)t+=(2*k[j]-n[j])*profilesP[j][i];
    for(int j=0;j<3;j++)t+=profilesP[r][j]*M[j][i];
    if(mod(t)){ok=false;break;}
   }
   int p=(a+b+c+d)%2;
   if(ok)T[r][p].push_back({r,p,k});
 }
 int cpt=0,reached=0;
 for(int a=0;a<=3;a++)for(int b=1;b<=4;b++)for(int c=1;c<=4;c++)for(int d=1;d<=4;d++){
  array<int,4> np{a,b,c,d};++cpt;
  int rhs[16]={};for(int s=1;s<16;s++)for(int j=0;j<4;j++)if(s>>j&1){
   const int u=np[j],v=n[j]-u;
   rhs[s]+=6*(u*(u-1)+v*(v-1))+9*u*v;
  }
  int tot[16]={};
  for(int r=0;r<4;r++)for(int p=0;p<2;p++){
   const int size=(p?np[r]:n[r]-np[r]);if(!size)continue;
   int best[16];for(int i=0;i<16;i++)best[i]=-10000;
   for(auto x:T[r][p]){
    int score[4];
    for(int j=0;j<4;j++)score[j]=x.k[j]*(n[j]-x.k[j])-(j==r?3*max(0,x.k[j]+1-size):0);
    int z[16]={};for(int u=1;u<16;u++){
      z[u]=z[u&(u-1)]+score[countr_zero((unsigned)u)];
      best[u]=max(best[u],z[u]);
    }
   }
   for(int u=1;u<16;u++)tot[u]+=size*best[u];
  }
  bool ok=true;
  for(int u=1;u<16;u++)if(tot[u]<rhs[u]){ok=false;break;}
  if(ok){reached++;cout<<a<<" "<<b<<" "<<c<<" "<<d<<endl;}
 }
 assert(cpt==256&&reached==0);
 cout<<"parity splits "<<cpt<<", remaining "<<reached<<endl;
}
void test_four(){
 int tested=0,out=0;
 for(int m=0;m<256;m++)if(popcount((unsigned)m)==5){
  int p[5],n=0;
  for(int i=0;i<8;i++)if(m>>i&1)p[n++]=i;
  for(int s=0;s<1024;s++){
   tested++;
   int R[5]={},k=0;
   for(int i=0;i<5;i++)for(int j=i+1;j<5;j++,k++)if(s>>k&1){R[i]|=1<<j;R[j]|=1<<i;}
   bool ok=true;
   for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)
    if((13*popcount((unsigned)(p[i]^p[j]))+popcount((unsigned)(R[i]^R[j])))%3)ok=false;
   if(ok)out++;
  }
 }
 assert(tested==57344&&out==0);
 cout<<"four-cycle fixed profile/sign fixtures "<<tested<<" surviving "<<out<<endl;
}
int main(){test_four();test_three();parity_test();}
