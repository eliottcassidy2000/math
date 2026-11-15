#include<bits/stdc++.h>
using namespace std; using LL=long long;
int main(int argc,char**argv){
const LL alpha=argc>1?atoll(argv[1]):6;
  int n;uint64_t bits;
  int a[12],b[12];
  static LL en[1024][10],st[1024][10],h[1024],ca[10],cu[10],cf[10],Q2[1024],cw[10][10],cw3[10][10];
  map<tuple<LL,LL>,map<LL,tuple<int,uint64_t,int,LL,LL,LL,LL>>> same;
  while(cin>>n>>bits){
    int f=(1<<n)-1;
    fill(a,a+n,0);fill(b,b+n,0);
    int t=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
      int x=i,y=j; if(!(bits>>t++&1))swap(x,y); a[x]|=1<<y;b[y]|=1<<x;
    }
    memset(en,0,sizeof(en));memset(st,0,sizeof(st));
    for(int i=0;i<n;i++)en[1<<i][i]=st[1<<i][i]=1;
    h[0]=1;
    for(int s=1;s<=f;s++){
      h[s]=0;
      for(int m=s;m;m&=m-1){
        int i=__builtin_ctz(m);
        h[s]+=en[s][i];
        if(en[s][i])for(int u=a[i]&~s;u;u&=u-1) en[s|(u&-u)][__builtin_ctz(u)]+=en[s][i];
        if(st[s][i])for(int u=b[i]&~s;u;u&=u-1) st[s|(u&-u)][__builtin_ctz(u)]+=st[s][i];
      }
    }
    for(int s=0;s<=f;s++){
      Q2[s]=0;
      for(int x=s;;x=(x-1)&s){Q2[s]+=h[x]*h[s^x];if(!x)break;}
    }
    for(int i=0;i<n;i++){
      ca[i]=cu[i]=cf[i]=0;
      for(int j=0;j<n;j++)cw[i][j]=cw3[i][j]=0;
    }
    for(int z=1;z<=f;z++){
      int r=f^z;
      for(int i=0;i<n;i++)cf[i]+=h[r]*en[z][i];
      for(int o=z;o;o&=o-1){
        int j=__builtin_ctz(o);LL v=st[z][j];if(!v)continue;
        ca[j]+=h[r]*v;
        cu[j]+=Q2[r]*v;
        for(int x=r;x;x=(x-1)&r){
          for(int q=x;q;q&=q-1){
            int i=__builtin_ctz(q);LL w=v*en[x][i];if(!w)continue;
            cw3[i][j]+=w*h[r^x];if(x==r)cw[i][j]+=w;
          }
        }
      }
    }
    for(int mask=1;mask<=f;mask++){
      LL M=0,R=Q2[f],U=0,r=h[f],F=0;
      for(int o=mask;o;o&=o-1){
        int i=__builtin_ctz(o);M+=ca[i];U+=cu[i];F+=cf[i];
        for(int p=mask;p;p&=p-1){
          int j=__builtin_ctz(p);R-=cw3[i][j];r-=cw[i][j];
        }
      }
      same[{U,2*alpha*R-M}][r+2*alpha*F]={n,bits,mask,M,R,r,F};
    }
  }
  for(auto &row:same){
    if(row.second.size()<2)continue;
    auto [U,A]=row.first;
    cout<<U<<' '<<A;
    for(auto &p:row.second){auto [n,bits,mask,M,R,r,F]=p.second;
      cout<<" / "<<n<<' '<<p.first<<' '<<bits<<' '<<mask<<' '<<M<<' '<<R<<' '<<r<<' '<<F;
    }
    cout<<'\n';
  }
}
