/* Exact every-N source-weighted clock correlation ledger. No orbit assumptions:
 * capped or overflowing routes abort; every positive result has a finite scope.
 * cc -O3 -Wall -Wextra root_wiener_scan_20261018o.c -o /tmp/root-wiener */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <limits.h>

static int64_t norm(int64_t a,int64_t b){return a*a-a*b+b*b;}
static int64_t bil(int64_t a,int64_t b,int64_t c,int64_t d){return 2*a*c-a*d-b*c+2*b*d;}
static void fail(const char *s){fprintf(stderr,"%s\n",s);exit(2);}
static uint64_t seed(uint64_t x){while(x%8==5)x=(x-1)/4;return x;}
static uint64_t edge(uint64_t x,int *dt){
 if(x>(UINT64_MAX-1)/3)fail("overflow guard: unresolved orbit");
 int a=x%8==1?2:1; x=(3*x+1)>>a;*dt=a;
 while(x%8==5){x=(x-1)/4;*dt+=2;}return x;
}

int main(int argc,char **argv){
 if(argc!=3)fail("usage: root-wiener FRONTIER MAX_POWER (12..22)");
 int power=atoi(argv[2]);if(power<12||power>22)fail("unsupported power");
 int bound=1<<power,limit=4*bound,nv=0,maxv=0;
 FILE *f=fopen(argv[1],"r");if(!f)fail("frontier missing");
 int vs[1000];char line[256];
 while(fgets(line,sizeof(line),f)){int x;if(line[0]!='#'&&sscanf(line,"%d",&x)==1){
  if(nv==1000)fail("frontier overflow");vs[nv++]=x;if(x>maxv)maxv=x;
 }}fclose(f);if(nv!=627)fail("wrong frontier size");
 int *vi=malloc((maxv+1)*sizeof(int));for(int i=0;i<=maxv;i++)vi[i]=-1;
 for(int i=0;i<nv;i++)vi[vs[i]]=i;
 int *port=malloc(limit*sizeof(int)),*clock=calloc(limit,sizeof(int));
 if(!vi||!port||!clock)fail("allocation");for(int i=0;i<limit;i++)port[i]=-1;
 int maxclock=0;uint64_t path[10000];int dts[10000];
 for(int s=1;s<limit;s+=2){if(s%8==5)continue;
  uint64_t x=s;int k=0,u=-1,t=0;
  while(1){
   if(x<=(uint64_t)maxv&&vi[x]>=0){u=vi[x];break;}
   if(x<(uint64_t)limit&&port[x]>=0){u=port[x];t=clock[x];break;}
   if(k==10000)fail("capped unresolved route");path[k]=x;x=edge(x,&dts[k]);k++;
  }
  while(k){--k;t+=dts[k];if(path[k]<(uint64_t)limit){port[path[k]]=u;clock[path[k]]=t;}}
  if(port[s]<0){port[s]=u;clock[s]=t;}if(clock[s]>maxclock)maxclock=clock[s];
 }
 int width=maxclock+1;
 int64_t *qa=calloc(nv*width,sizeof(int64_t)),*qb=calloc(nv*width,sizeof(int64_t));
 int64_t *sa=calloc(limit,sizeof(int64_t)),*sb=calloc(limit,sizeof(int64_t));
 int64_t *wa=calloc(nv,sizeof(int64_t)),*wb=calloc(nv,sizeof(int64_t));
 int64_t *lag=calloc(width,sizeof(int64_t));
 if(!qa||!qb||!sa||!sb||!wa||!wb||!lag)fail("allocation");
 const int units[3][2]={{1,0},{0,1},{-1,-1}};
 int64_t S=0,D=0,Pi=0,maxnum=0,maxden=1;int peakN=0,bad=0,firstbad=0;
 printf("{\"status\":\"FINITE-EXACT\",\"bound_exclusive\":%d,\"max_clock\":%d,\"samples\":[",bound,maxclock);
 int first=1;
 for(int n=1;n<bound;n++){
  int ms[3]={4*n-3,4*n-1,2*n-1};
  for(int e=0;e<3;e++){
   int m=ms[e],s=(int)seed(m),u=port[s],t=clock[s],r=(m+(e==2?2:1))%3;
   int a=units[r][0],b=units[r][1];if(e==2){a=-a;b=-b;}
   S+=bil(sa[s],sb[s],a,b)+norm(a,b);sa[s]+=a;sb[s]+=b;
   int offset=u*width;
   for(int v=0;v<width;v++)if(v!=t){
    lag[abs(v-t)]+=bil(a,b,qa[offset+v],qb[offset+v]);
   }
   D+=bil(qa[offset+t],qb[offset+t],a,b)+norm(a,b);
   qa[offset+t]+=a;qb[offset+t]+=b;
   if(t%2){a=-a;b=-b;}
   Pi+=bil(wa[u],wb[u],a,b)+norm(a,b);wa[u]+=a;wb[u]+=b;
  }
  int64_t W=D,ledger=D;for(int l=1;l<width;l++){W+=llabs(lag[l]);ledger+=(l%2?-1:1)*lag[l];}
  if(Pi!=ledger)fail("lag ledger failed");
  if(n>=4096){
   if((__int128)W*maxden>(__int128)maxnum*S){maxnum=W;maxden=S;peakN=n;}
   if(2*W>3*S){bad++;if(!firstbad)firstbad=n;}
  }
  if((n&(n-1))==0||n==1204||n==5997||n==7910){
   printf("%s{\"N\":%d,\"source_norm\":%"PRId64",\"clock_diagonal\":%"PRId64",\"pi_norm\":%"PRId64",\"wiener\":%"PRId64"}",first?"":",",n,S,D,Pi,W);first=0;
  }
 }
 printf("],\"checked_N\":%d,\"violations_3_over_2\":%d,\"first_violation\":%d,\"peak_N\":%d,\"peak_numerator\":%"PRId64",\"peak_denominator\":%"PRId64",\"unresolved_routes\":0}\n",bound-4096,bad,firstbad,peakN,maxnum,maxden);
 return 0;
}
