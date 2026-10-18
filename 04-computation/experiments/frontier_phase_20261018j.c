/* Exact every-N three-phase first-arrival Gram scan. Independent incremental
   source construction from frontier_phase_20261018j.py. Routes have a cap and
   checked integer arithmetic. No unconditional convergence assumption. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#define CLK 4096
#define CAP 10000
static int64_t seed(int64_t x){while((x&7)==5)x=(x-1)/4;return x;}
static int64_t nxt(int64_t x,int *t){
 assert(x>0 && x<=(INT64_MAX-1)/3);
 int a=(x&7)==1?2:1; x=(3*x+1)>>a; *t=a;
 while((x&7)==5){x=(x-1)/4;*t+=2;}return x;
}
static void upd(int64_t *g,int64_t *v,int *d){
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)g[3*i+j]+=v[i]*d[j]+v[j]*d[i]+(int64_t)d[i]*d[j];
 for(int i=0;i<3;i++)v[i]+=d[i];
 assert(v[0]+v[1]+v[2]==0);
}
static int64_t tr(int64_t *g){return g[0]+g[4]+g[8];}
static void arr(int64_t *g){printf("[");for(int i=0;i<9;i++)printf("%s%lld",i?",":"",(long long)g[i]);printf("]");}
int main(int argc,char**argv){
 assert(argc==3); int power=atoi(argv[2]); assert(power>=3&&power<=26);
 int64_t lim=1LL<<power,M=4*lim; int Vs[4096],nv=0,maxV=0;
 FILE *f=fopen(argv[1],"r");assert(f);char line[256];
 while(fgets(line,sizeof(line),f)){int x;if(line[0]!='#'&&sscanf(line,"%d",&x)==1){assert(nv<4096);Vs[nv++]=x;if(x>maxV)maxV=x;}}
 fclose(f);assert(nv==627);
 int *vi=malloc((maxV+1)*sizeof(int));assert(vi);for(int i=0;i<=maxV;i++)vi[i]=-1;
 for(int i=0;i<nv;i++)vi[Vs[i]]=i;
 for(int i=0;i<nv;i++){int t;int64_t y=nxt(Vs[i],&t);assert(y<=maxV&&vi[y]>=0);}
 int16_t *port=malloc((M/2)*sizeof(int16_t)),*clock=malloc((M/2)*sizeof(int16_t));
 int8_t *coeff=calloc(M/2*3,sizeof(int8_t));assert(port&&clock&&coeff);
 int maxclock=0;
 for(int64_t i=0;i<M/2;i++){
   int64_t x=2*i+1;port[i]=-1;clock[i]=0;if((x&7)==5)continue;
   int total=0,steps=0;
   while(!(x<=maxV&&vi[x]>=0)){
     assert(++steps<CAP);int t;x=nxt(x,&t);total+=t;assert(total<CLK);
   }
   port[i]=vi[x];clock[i]=total;if(total>maxclock)maxclock=total;
 }
 fprintf(stderr,"routes complete: %lld odd source slots; max clock %d\n",(long long)(M/2),maxclock);
 int64_t *arrival=calloc(nv*3,sizeof(int64_t)),*hist=calloc((size_t)nv*CLK*3,sizeof(int64_t));assert(arrival&&hist);
 int64_t g[9]={0},s[9]={0},h[9]={0};
 int64_t count[32]={0},nratio[32]={0},nte[32]={0},ntc[32]={0};
 int64_t viol[3]={0},first[3]={0}; double maxratio[32]={0},maxt[32]={0},maxc[32]={0};
 printf("{\"status\":\"FINITE-EXACT\",\"limit_exclusive\":%lld,\"dyadic_checks\":[",(long long)lim);int sep=0;
 for(int64_t N=1;N<lim;N++){
   for(int event=0;event<3;event++){
     int64_t m=event==0?4*N-3:event==1?4*N-1:2*N-1;
     int64_t x=seed(m),i=(x-1)/2;int u=port[i],t=clock[i];assert(u>=0);
     int d[3];for(int r=0;r<3;r++){int res=(m+r)%3;d[r]=event==2?(res==1?-2:1):(res==2?2:-1);}
     int64_t old[3];for(int r=0;r<3;r++)old[r]=coeff[3*i+r];
     upd(s,old,d);for(int r=0;r<3;r++){assert(old[r]>=-6&&old[r]<=6);coeff[3*i+r]=old[r];if(t&1)d[r]=-d[r];}
     upd(g,arrival+3*u,d);upd(h,hist+((size_t)t*nv+u)*3,d);
   }
   int k=63-__builtin_clzll(N);count[k]++;
   int64_t st=tr(s),gt=tr(g),ht=tr(h);assert(st>=0&&gt>=0&&ht>=0);
   double ratio=(double)gt/st,te=(double)gt/N,tc=(double)ht/N;
   if(ratio>maxratio[k]){maxratio[k]=ratio;nratio[k]=N;}
   if(te>maxt[k]){maxt[k]=te;nte[k]=N;}
   if(tc>maxc[k]){maxc[k]=tc;ntc[k]=N;}
   if(N>=4096){int fail[3]={2*gt>3*st,gt>st,ht>st};for(int j=0;j<3;j++)if(fail[j]){viol[j]++;if(!first[j])first[j]=N;}}
   if((N&(N-1))==0){printf("%s{\"N\":%lld,\"source\":",sep++?",":"",(long long)N);arr(s);printf(",\"arrival\":");arr(g);printf(",\"clock\":");arr(h);printf("}");}
 }
 printf("],\"octaves\":[");
 for(int k=0;k<power;k++)printf("%s{\"k\":%d,\"max_trace_ratio\":%.12g,\"N_ratio\":%lld,\"max_trace_per_N\":%.12g,\"N_trace\":%lld,\"max_clock_trace_per_N\":%.12g,\"N_clock\":%lld}",k?",":"",k,maxratio[k],(long long)nratio[k],maxt[k],(long long)nte[k],maxc[k],(long long)ntc[k]);
 printf("],\"violations\":[%lld,%lld,%lld],\"first_violations\":[%lld,%lld,%lld],\"max_clock\":%d}\n",(long long)viol[0],(long long)viol[1],(long long)viol[2],(long long)first[0],(long long)first[1],(long long)first[2],maxclock);
 free(arrival);free(hist);free(port);free(clock);free(coeff);free(vi);return 0;
}
