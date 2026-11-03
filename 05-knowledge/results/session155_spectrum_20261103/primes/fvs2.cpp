#include <bits/stdc++.h>
using namespace std; const int B=100000; int base[B+1], mn[B+1]; uint32_t wit[B+1];
void rec(int n,int h,int a,int p){ if(p==0&&n>=2&&h<=B)base[h]=1; if(n==21)return; if(a<=B&&h+a+p<=B)rec(n+1,h+a+p,2*a,0);if(a<=(B-1)/2)rec(n+1,h,2*a+1,1);}
int main(){ rec(1,1,1,1); for(int m=0;m<=8;m++){int n=m+2; int M=1<<n,lim=1<<(2*m+1); long long ns=0,unique=0;vector<int>d(M*n);for(int mask=0;mask<lim;mask++){int adj[10]={};for(int i=0;i<m;i++)for(int j=i+1;j<m;j++)adj[i]|=1<<j;for(int i=0;i<m;i++)for(int a=0;a<2;a++) {int v=m+a; if((mask>>(2*i+a))&1)adj[v]|=1<<i; else adj[i]|=1<<v;} if(mask>>(2*m)&1)adj[m]|=1<<(m+1); else adj[m+1]|=1<<m;
 int reach=1;for(int t=0;t<n;t++){int r=reach;for(int i=0;i<n;i++)if(r>>i&1)reach|=adj[i];} if(reach!=M-1)continue; int rev=1; for(int t=0;t<n;t++){int r=rev;for(int i=0;i<n;i++)if(adj[i]&r)rev|=1<<i;}if(rev!=M-1)continue;ns++;
 fill(d.begin(),d.end(),0);for(int i=0;i<n;i++)d[(1<<i)*n+i]=1;
 for(int s=1;s<M;s++)for(int i=0;i<n;i++)if((s>>i&1)&&d[s*n+i]){int avail=adj[i]&~s; while(avail){int j=__builtin_ctz((unsigned)avail); avail&=avail-1;d[(s|1<<j)*n+j]+=d[s*n+i];}}
 int h=0;for(int i=0;i<n;i++)h+=d[(M-1)*n+i];if(h<=B&&!mn[h]){mn[h]=n;wit[h]=mask;unique++;}}
 fprintf(stderr,"m=%d strong=%lld new=%lld\n",m,ns,unique);}
 vector<bool>p(B+1,true);p[0]=p[1]=0;for(int i=2;i*i<=B;i++)if(p[i])for(int j=i*i;j<=B;j+=i)p[j]=false;int cntp=0,cnt7p=0;for(int v=3;v<=B;v+=2)if(mn[v]&&!base[v]){if(p[v])cntp++; if(v%7==0&&p[v/7])cnt7p++;}printf("new prime %d; new 7p %d\n",cntp,cnt7p);for(int v=3;v<=3000;v+=2)if(mn[v]&&!base[v]&&(p[v]||(v%7==0&&p[v/7])))printf("%d %s n=%d mask=%u\n",v,p[v]?"prime":"7p",mn[v],wit[v]);}
