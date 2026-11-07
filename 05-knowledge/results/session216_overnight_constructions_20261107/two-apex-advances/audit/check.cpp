#include <bits/stdc++.h>
using namespace std;
using V=array<uint64_t,36>;
int nx[4][36][3]; bool fin[36];
void init(){for(int p=0;p<2;p++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int c=0;c<2;c++){
 int id=p*18+a*6+b*2+c; fin[id]=(a!=1 && b!=1 && (b!=0||p==0));
 for(int m=0;m<4;m++)for(int g=0;g<3;g++){
 int q[3]={a,b,c}; int prev=g==1?p:1-p, next=g==0?p:1-p;
 nx[m][id][g]=-1;
 if(q[g]==0 && g && ((m>>prev)&1)==0)continue;
 q[g]=g==2?1:(((m>>next)&1)?1:2);
 nx[m][id][g]=p*18+q[0]*6+q[1]*2+q[2];
 }
}}
uint64_t H(vector<int> w){V v={};v[0]=v[18]=1;for(int m:w){V z={};for(int i=0;i<36;i++)for(int g=0;g<3;g++)if(nx[m][i][g]>=0)z[nx[m][i][g]]+=v[i];v=z;}uint64_t r=0;for(int i=0;i<36;i++)if(fin[i])r+=v[i];return r;}
uint64_t subset(vector<int>w){int n=w.size()+2;vector<int> out(n);out[0]|=2;for(int i=0;i<(int)w.size();i++){for(int a=0;a<2;a++)if((w[i]>>a)&1)out[a]|=1<<(i+2);else out[i+2]|=1<<a;for(int j=i+1;j<(int)w.size();j++)out[i+2]|=1<<(j+2);}vector<uint64_t> d((1<<n)*n);for(int i=0;i<n;i++)d[(1<<i)*n+i]=1;for(int s=1;s<(1<<n);s++)for(int i=0;i<n;i++)if(d[s*n+i]){int t=out[i]&~s;while(t){int j=__builtin_ctz(t);t&=t-1;d[(s|(1<<j))*n+j]+=d[s*n+i];}}uint64_t r=0;for(int i=0;i<n;i++)r+=d[((1<<n)-1)*n+i];return r;}
int main(int argc,char**argv){init();ifstream f(argv[1],ios::binary);uint32_t x;size_t h=0,checked=0,dp=0;mt19937 rng(216);while(f.read((char*)&x,4)){if(x!=UINT32_MAX){int n=x>>28;if(n>14||(x&0xfffffff)>>(2*n)){cerr<<"bad format";return 1;}vector<int>w;for(int j=0;j<n;j++)w.push_back((x>>(2*j))&3);if(H(w)!=h){cerr<<"bad H "<<h;return 2;}if((h%2&&h<201)||(h%2&&h>2327473)||(rng()%12000==0)){if(subset(w)!=h){cerr<<"bad DP "<<h;return 3;}dp++;}checked++;}else if(h%2&&h!=7&&h!=21){cerr<<"missing "<<h;return 4;}h++;}vector<int>w={1,2,0,0,2,1,0,0,2,0,3,3,0,1,0};cout<<"records="<<h<<" checked="<<checked<<" subset="<<dp<<" boundary_transfer="<<H(w)<<" boundary_DP="<<subset(w)<<"\n";return 0;}
