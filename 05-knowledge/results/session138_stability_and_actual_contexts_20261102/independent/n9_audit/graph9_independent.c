/* Independent low-memory n9 audit: S7 BFS, equitable-refinement canonical8,
   vertex-addition matching/wedge recurrences, full-DP snapshots. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static int e7[21][2],e8[28][2],e9[36][2],ix7[9][9],ix8[9][9],ix9[9][9];
static unsigned char seen[1u<<21];static unsigned queue[5040],aug8[133632];
static unsigned canadj[9];static int cell[8][8],clen[8],ncells,order[8];static unsigned canbest;
static int dp8[512][5],dp9[512][5],mm[256][4],pp[256],ds[256],sz[256];
static void need(int b,const char*s){if(!b){fprintf(stderr,"FAIL %s\n",s);exit(1);}}
static int pop(unsigned b){return __builtin_popcount(b);}
static int gcd(int a,int b){while(b){int t=a%b;a=b;b=t;}return a;}
static int cmppu(const void*a,const void*b){unsigned x=*(const unsigned*)a,y=*(const unsigned*)b;return (x>y)-(x<y);}
static int nextperm(int*a,int n){int i=n-2,j,t;while(i>=0&&a[i]>=a[i+1])i--;if(i<0){for(j=0;j<n/2;j++){t=a[j];a[j]=a[n-1-j];a[n-1-j]=t;}return 0;}j=n-1;while(a[j]<=a[i])j--;t=a[i];a[i]=a[j];a[j]=t;for(j=i+1;j<(n+i+1)/2;j++){t=a[j];a[j]=a[n+i-j];a[n+i-j]=t;}return 1;}
static void canonrec(int k,int off){
 if(k==ncells){unsigned g=0;for(int e=0;e<28;e++)g|=((canadj[order[e8[e][0]]]>>order[e8[e][1]])&1u)<<e;if(g<canbest)canbest=g;return;}
 do{for(int j=0;j<clen[k];j++)order[off+j]=cell[k][j];canonrec(k+1,off+clen[k]);}while(nextperm(cell[k],clen[k]));
}
static int sigless(const int*a,const int*b){for(int i=0;i<10;i++)if(a[i]!=b[i])return a[i]<b[i];return 0;}
static unsigned canonical8(const unsigned*a){
 int col[8],next[8],sig[8][10],ids[8];for(int i=0;i<8;i++)col[i]=pop(a[i]);
 for(;;){memset(sig,0,sizeof sig);for(int i=0;i<8;i++){sig[i][0]=col[i];ids[i]=i;for(int j=0;j<8;j++)if(a[i]>>j&1)sig[i][1+col[j]]++;}
  for(int i=1;i<8;i++){int x=ids[i],j=i;while(j&&sigless(sig[x],sig[ids[j-1]])){ids[j]=ids[j-1];j--;}ids[j]=x;}
  int c=0;for(int i=0;i<8;i++){if(i&&memcmp(sig[ids[i]],sig[ids[i-1]],sizeof sig[0]))c++;next[ids[i]]=c;}
  if(!memcmp(next,col,sizeof col))break;memcpy(col,next,sizeof col);
 }
 memset(clen,0,sizeof clen);ncells=0;for(int i=0;i<8;i++){cell[col[i]][clen[col[i]]++]=i;if(col[i]>=ncells)ncells=col[i]+1;}memcpy(canadj,a,9*sizeof(unsigned));canbest=~0u;canonrec(0,0);return canbest;
}
static void decode(uint64_t g,int n,unsigned*a){memset(a,0,9*sizeof(unsigned));int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(g>>k&1){a[i]|=1u<<j;a[j]|=1u<<i;}}
static void matching(const unsigned*a,int n,int dp[512][5]){
 memset(dp,0,512*5*sizeof(int));dp[0][0]=1;for(unsigned s=1;s<(1u<<n);s++){unsigned vb=s&-s,t=s^vb;int v=__builtin_ctz(vb);memcpy(dp[s],dp[t],5*sizeof(int));for(unsigned nn=a[v]&t;nn;nn&=nn-1){unsigned ub=nn&-nn;for(int j=1;j<=4;j++)dp[s][j]+=dp[t^ub][j-1];}}
}
static void fullcounts(const unsigned*a,int n,int dp[512][5],int*out){
 matching(a,n,dp);unsigned all=(1u<<n)-1;for(int j=0;j<4;j++)out[j]=dp[all][j+1];out[4]=out[5]=0;
 for(int v=0;v<n;v++){int d=pop(a[v]);out[4]+=d*(d-1)/2;for(unsigned nn=a[v];nn;nn&=nn-1){unsigned ub=nn&-nn;for(unsigned kk=nn^ub;kk;kk&=kk-1){unsigned wb=kk&-kk;out[5]+=dp[all^((1u<<v)|ub|wb)][1];}}}
}
int main(void){
 int z=0;for(int i=0;i<7;i++)for(int j=i+1;j<7;j++){ix7[i][j]=ix7[j][i]=z;e7[z][0]=i;e7[z++][1]=j;}
 z=0;for(int i=0;i<8;i++)for(int j=i+1;j<8;j++){ix8[i][j]=ix8[j][i]=z;e8[z][0]=i;e8[z++][1]=j;}
 z=0;for(int i=0;i<9;i++)for(int j=i+1;j<9;j++){ix9[i][j]=ix9[j][i]=z;e9[z][0]=i;e9[z++][1]=j;}
 unsigned swap[6][5][2];for(int v=0;v<6;v++){int t=0;for(int x=0;x<7;x++)if(x!=v&&x!=v+1){swap[v][t][0]=1u<<ix7[v][x];swap[v][t++][1]=1u<<ix7[v+1][x];}}
 unsigned covered=0,reps7=0,augcount=0;
 for(unsigned g=0;g<(1u<<21);g++)if(!seen[g]){
  unsigned count=1;queue[0]=g;seen[g]=1;for(unsigned h=0;h<count;h++)for(int v=0;v<6;v++){unsigned image=queue[h];for(int t=0;t<5;t++){unsigned p=swap[v][t][0],q=swap[v][t][1];if(!!(image&p)!=!!(image&q))image^=p|q;}if(!seen[image]){need(count<5040,"orbit queue overflow");seen[image]=1;queue[count++]=image;}}
  covered+=count;reps7++;unsigned base[9];decode(g,7,base);
  for(unsigned attach=0;attach<128;attach++){unsigned a[9];memcpy(a,base,sizeof a);a[7]=attach;for(int i=0;i<7;i++)if(attach>>i&1)a[i]|=128u;need(augcount<133632,"augmentation overflow");aug8[augcount++]=canonical8(a);}
 }
 need(covered==2097152&&reps7==1044&&augcount==133632,"n7 coverage");qsort(aug8,augcount,sizeof(unsigned),cmppu);unsigned reps8=0;for(unsigned i=0;i<augcount;i++)if(!i||aug8[i]!=aug8[i-1])aug8[reps8++]=aug8[i];need(reps8==12346,"n8 orbit count");fprintf(stderr,"independent n8 representatives %u\n",reps8);
 uint64_t cases=0,digest=1469598103934665603ULL,zero_masks[32];int minimum=1000000,maximum=-1,zeros=0,snapshots=0,maxnum=0,maxden=1,zero_degrees[32][9];
 for(unsigned gg=0;gg<reps8;gg++){
  unsigned g=aug8[gg],base[9];decode(g,8,base);int old[6],deg[8],contains[8]={0},leaf[8]={0},lin[8],pair[8][8]={{0}};fullcounts(base,8,dp8,old);
  for(int v=0;v<8;v++){deg[v]=pop(base[v]);for(unsigned nn=base[v];nn;nn&=nn-1){unsigned ub=nn&-nn;for(unsigned kk=nn^ub;kk;kk&=kk-1){unsigned wb=kk&-kk;contains[v]++;contains[__builtin_ctz(ub)]++;contains[__builtin_ctz(wb)]++;}}}
  for(int v=0;v<8;v++)for(int u=0;u<8;u++)if(v!=u){pair[v][u]=dp8[255^((1<<v)|(1<<u))][1];if(base[v]>>u&1)leaf[v]+=pair[v][u];}
  for(int v=0;v<8;v++)lin[v]=old[4]-contains[v]+leaf[v];for(int j=0;j<4;j++)mm[0][j]=old[j];pp[0]=old[5];ds[0]=sz[0]=0;
  uint64_t base9=0;for(int e=0;e<28;e++)if(g>>e&1)base9|=1ull<<ix9[e8[e][0]][e8[e][1]];
  for(unsigned attach=0;attach<256;attach++){
   if(attach){unsigned vb=attach&-attach,t=attach^vb;int v=__builtin_ctz(vb);sz[attach]=sz[t]+1;ds[attach]=ds[t]+deg[v];for(int j=0;j<4;j++)mm[attach][j]=mm[t][j]+dp8[255^vb][j];pp[attach]=pp[t]+lin[v];for(unsigned nn=t;nn;nn&=nn-1)pp[attach]+=pair[v][__builtin_ctz(nn)];}
   int w=old[4]+ds[attach]+sz[attach]*(sz[attach]-1)/2,*m=mm[attach];need(m[1]+w==m[0]*(m[0]-1)/2,"edge-pair partition");int b=81*m[0]-64*m[1]+16*m[2]-4*m[3]+16*pp[attach],value=b+19*w;
   need(value>=0,"negative n9 residual");if(value<minimum)minimum=value;if(value>maximum)maximum=value;if(w&&(int64_t)(-b)*maxden>(int64_t)maxnum*64*w){maxnum=-b;maxden=64*w;int d=gcd(maxnum,maxden);maxnum/=d;maxden/=d;}
   if(cases<16||cases%10007==0||!value){unsigned a[9];memcpy(a,base,sizeof a);a[8]=attach;for(int i=0;i<8;i++)if(attach>>i&1)a[i]|=256u;int cc[6];fullcounts(a,9,dp9,cc);for(int j=0;j<4;j++)need(cc[j]==m[j],"full matching snapshot");need(cc[4]==w&&cc[5]==pp[attach],"full wedge snapshot");snapshots++;}
   if(!value){need(zeros<32,"zero list overflow");uint64_t h=base9;zero_degrees[zeros][0]=sz[attach];for(int i=0;i<8;i++){zero_degrees[zeros][i+1]=deg[i]+((attach>>i)&1);if(attach>>i&1)h|=1ull<<ix9[i][8];}for(int i=1;i<9;i++){int x=zero_degrees[zeros][i],j=i;while(j&&zero_degrees[zeros][j-1]>x){zero_degrees[zeros][j]=zero_degrees[zeros][j-1];j--;}zero_degrees[zeros][j]=x;}zero_masks[zeros++]=h;}
   digest^=(uint64_t)value;digest*=1099511628211ULL;cases++;
  }
 }
 need(cases==3160576&&minimum==0&&zeros==14&&maxnum==19&&maxden==64,"n9 summary mismatch");
 printf("{\"status\":\"PASS\",\"method\":\"S7 adjacent-transposition BFS; equitable-refinement canonical8; vertex-addition matching/wedge recurrences\",\"reps7\":%u,\"covered7\":%u,\"augmentations8\":%u,\"reps8\":%u,\"augmentations9\":%llu,\"minimum\":%d,\"maximum\":%d,\"zero_augmentations\":%d,\"maximum_required_alpha\":\"%d/%d\",\"full_dp_snapshots\":%d,\"residual_digest\":%llu,\"zero_witnesses\":[",reps7,covered,augcount,reps8,(unsigned long long)cases,minimum,maximum,zeros,maxnum,maxden,snapshots,(unsigned long long)digest);
 for(int z=0;z<zeros;z++){if(z)putchar(',');printf("{\"mask\":%llu,\"degrees\":[",(unsigned long long)zero_masks[z]);for(int i=0;i<9;i++){if(i)putchar(',');printf("%d",zero_degrees[z][i]);}printf("]}");}printf("]}\n");return 0;
}
