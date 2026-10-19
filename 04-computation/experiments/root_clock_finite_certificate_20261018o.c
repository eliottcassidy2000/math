/* Verify all clock caps and native-depth caps through T, at every scale below
 * the PROVED all-N tail cutoff supplied by the complete inverse-tree receipt.
 * Input routes are independently replayed by the Python producer. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
static void fail(const char*s){fprintf(stderr,"%s\n",s);exit(2);}
static int64_t norm(int64_t a,int64_t b){return a*a-a*b+b*b;}
static int64_t bil(int64_t a,int64_t b,int64_t c,int64_t d){return 2*a*c-a*d-b*c+2*b*d;}
static int seed(int x){while(x%8==5)x=(x-1)/4;return x;}
int main(int argc,char**argv){
 if(argc!=2)fail("usage: clock-cert ROUTES");FILE*f=fopen(argv[1],"r");if(!f)fail("routes missing");
 int nv,T,bound;int64_t k2;if(fscanf(f,"%d %d %d %"SCNd64,&nv,&T,&bound,&k2)!=4)fail("header");
 if(nv!=627||T<0||T>16||bound<4096||5LL*bound-32<8*k2)fail("invalid bound");
 int width=T+1,states=width*width,limit=4*bound;
 int *rp=malloc(limit*sizeof(int)),*rt=calloc(limit,sizeof(int)),*rd=calloc(limit,sizeof(int));
 if(!rp||!rt||!rd)fail("allocation");for(int i=0;i<limit;i++)rp[i]=-1;
 int s,u,t,d,p,count=0;int64_t*countu=calloc(nv,sizeof(int64_t));if(!countu)fail("allocation");
 while(fscanf(f,"%d %d %d %d %d",&s,&u,&t,&d,&p)==5){
  if(s<1||!(s%2)||s%8==5||u<0||u>=nv||t<0||t>T||d<0||d>t)fail("invalid route");
  if(s<limit){if(rp[s]>=0)fail("duplicate route");rp[s]=u;rt[s]=t;rd[s]=d;}countu[u]++;count++;
 }fclose(f);int64_t kcheck=0;for(u=0;u<nv;u++)kcheck+=countu[u]*countu[u];if(kcheck!=k2)fail("fibre bound mismatch");
 int64_t *sa=calloc(limit,sizeof(int64_t)),*sb=calloc(limit,sizeof(int64_t));
 size_t size=(size_t)states*nv*width;
 int64_t *qa=calloc(size,sizeof(int64_t)),*qb=calloc(size,sizeof(int64_t));
 int64_t *lag=calloc(states*width,sizeof(int64_t)),*D=calloc(states,sizeof(int64_t));
 int64_t *peak=calloc(states,sizeof(int64_t)),*den=calloc(states,sizeof(int64_t));int*peakn=calloc(states,sizeof(int));
 if(!sa||!sb||!qa||!qb||!lag||!D||!peak||!den||!peakn)fail("allocation");
 for(int k=0;k<states;k++)den[k]=1;
 const int units[3][2]={{1,0},{0,1},{-1,-1}};int64_t S=0;int64_t bad=0;int firstbad=0;
 for(int n=1;n<bound;n++){
  int ms[3]={4*n-3,4*n-1,2*n-1};
  for(int e=0;e<3;e++){
   int m=ms[e],s0=seed(m),r=(m+(e==2?2:1))%3;
   int a=units[r][0],b=units[r][1];if(e==2){a=-a;b=-b;}
   S+=bil(sa[s0],sb[s0],a,b)+norm(a,b);sa[s0]+=a;sb[s0]+=b;
   if(rp[s0]<0)continue;
   u=rp[s0];t=rt[s0];d=rd[s0];
   for(int tc=t;tc<=T;tc++)for(int dc=d;dc<=T;dc++){
    int k=tc*width+dc;size_t offset=((size_t)k*nv+u)*width;
    for(int v=0;v<width;v++)if(v!=t)lag[k*width+abs(v-t)]+=bil(a,b,qa[offset+v],qb[offset+v]);
    D[k]+=bil(qa[offset+t],qb[offset+t],a,b)+norm(a,b);qa[offset+t]+=a;qb[offset+t]+=b;
   }
  }
  if(n>=4096)for(int k=0;k<states;k++){
   int64_t W=D[k];for(int l=1;l<width;l++)W+=llabs(lag[k*width+l]);
   if((__int128)W*den[k]>(__int128)peak[k]*S){peak[k]=W;den[k]=S;peakn[k]=n;}
   if(W>S){bad++;if(!firstbad)firstbad=n;}
  }
 }
 printf("{\"status\":\"finite part of all-N bounded-clock certificate\",\"target_factor\":1,\"clock_cap\":%d,\"native_depth_caps\":%d,\"sources\":%d,\"N_start\":4096,\"N_end_exclusive\":%d,\"tail_K2\":%"PRId64",\"violations\":%"PRId64",\"first_violation\":%d,\"peaks\":[",T,width,count,bound,k2,bad,firstbad);
 for(int k=0;k<states;k++)printf("%s[%d,%d,%d,%"PRId64",%"PRId64"]",k?",":"",k/width,k%width,peakn[k],peak[k],den[k]);
 printf("]}\n");return bad?1:0;
}
