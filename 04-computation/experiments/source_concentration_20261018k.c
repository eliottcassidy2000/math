/* Source-specific cone and sibling packet concentration scan (20261018k).
   Forked from frontier_phase_20261018j.c. Independent incremental
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
 assert(argc==3||argc==4); int prefix=argc==4; int power=atoi(argv[2]); assert(power>=3&&power<=26);
 int64_t lim=1LL<<power,M=4*lim; int Vs[4096],nv=0,maxV=0;
 FILE *f=fopen(argv[1],"r");assert(f);char line[256];
 while(fgets(line,sizeof(line),f)){int x;if(line[0]!='#'&&sscanf(line,"%d",&x)==1){assert(nv<4096);Vs[nv++]=x;if(x>maxV)maxV=x;}}
 fclose(f);assert(nv==627);
 int *vi=malloc((maxV+1)*sizeof(int));assert(vi);for(int i=0;i<=maxV;i++)vi[i]=-1;
 for(int i=0;i<nv;i++)vi[Vs[i]]=i;
 for(int i=0;i<nv;i++){int t;int64_t y=nxt(Vs[i],&t);assert(y<=maxV&&vi[y]>=0);}
 int16_t *port=malloc((M/2)*sizeof(int16_t)),*clock=malloc((M/2)*sizeof(int16_t));
 int16_t *depth=malloc((M/2)*sizeof(int16_t));assert(depth);
 int8_t *coeff=calloc(M/2*3,sizeof(int8_t));assert(port&&clock&&coeff);
 int maxclock=0,maxdepth=0;
 for(int64_t i=0;i<M/2;i++){
   int64_t x=2*i+1;port[i]=-1;clock[i]=0;depth[i]=0;if((x&7)==5)continue;
   int total=0,steps=0;
   while(!(x<=maxV&&vi[x]>=0)){
     assert(++steps<CAP);int t;x=nxt(x,&t);total+=t;assert(total<CLK);
   }
   port[i]=vi[x];clock[i]=total;if(total>maxclock)maxclock=total;
   depth[i]=steps;if(steps>maxdepth)maxdepth=steps;
 }
 fprintf(stderr,"routes complete: %lld odd source slots; max clock %d\n",(long long)(M/2),maxclock);
 int64_t *arrival=calloc(nv*3,sizeof(int64_t)),*hist=calloc((size_t)nv*CLK*3,sizeof(int64_t));assert(arrival&&hist);
 int64_t *pv=NULL,*pt=NULL,*pe=NULL; int64_t sib=0;
 int64_t coneviol=0,conefirst=0,conelast=0,packetviol=0,packetfirst=0,packetlast=0;
 int64_t conemax=LLONG_MIN,conemaxN=0,packetmax=LLONG_MIN,packetmaxN=0; int conemaxD=0,packetmaxD=0;
 int64_t pviol[2]={0},pfirst[2]={0};double pratio=0;int64_t pn=0;int pd=0;
 if(prefix){pe=calloc(maxdepth+1,sizeof(int64_t));assert(pe);pv=calloc((size_t)nv*(maxdepth+1)*3,sizeof(int64_t));pt=calloc(maxdepth+1,sizeof(int64_t));assert(pv&&pt);}
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
     int64_t partner=(x%8==1)?2*x+1:((x%16==3)?(x-1)/2:0);
     if(partner>0 && partner<M && !(x<=maxV&&vi[x]>=0) && !(partner<=maxV&&vi[partner]>=0))
       sib+=(int64_t)d[0]*coeff[3*((partner-1)/2)];
     upd(s,old,d);for(int r=0;r<3;r++){assert(old[r]>=-6&&old[r]<=6);coeff[3*i+r]=old[r];if(t&1)d[r]=-d[r];}
     upd(g,arrival+3*u,d);upd(h,hist+((size_t)t*nv+u)*3,d);
     if(prefix)for(int z=depth[i];z<=maxdepth;z++){
       int64_t *v=pv+((size_t)z*nv+u)*3;
       pe[z]+=2*v[0]*d[0]+(int64_t)d[0]*d[0];
       for(int r=0;r<3;r++){pt[z]+=2*v[r]*d[r]+(int64_t)d[r]*d[r];v[r]+=d[r];}
     }
   }
   int k=63-__builtin_clzll(N);count[k]++;
   int64_t st=tr(s),gt=tr(g),ht=tr(h);assert(st>=0&&gt>=0&&ht>=0);
   double ratio=(double)gt/st,te=(double)gt/N,tc=(double)ht/N;
   if(ratio>maxratio[k]){maxratio[k]=ratio;nratio[k]=N;}
   if(te>maxt[k]){maxt[k]=te;nte[k]=N;}
   if(tc>maxc[k]){maxc[k]=tc;ntc[k]=N;}
   if(N>=4096){int fail[3]={2*gt>3*st,gt>st,ht>st};for(int j=0;j<3;j++)if(fail[j]){viol[j]++;if(!first[j])first[j]=N;}}
   if(prefix&&N>=4096){
     assert(pt[maxdepth]==gt && pe[maxdepth]==g[0]);
     int64_t base=3*pe[0]-pt[0],packet=s[0]-2*sib;
     int cv=0,pf=0;
     for(int z=0;z<=maxdepth;z++){
       int64_t q=3*pe[z]-pt[z]-base,dd=pe[z]-packet;
       if(q>0)cv=1;if(dd>0)pf=1;
       if(q>conemax){conemax=q;conemaxN=N;conemaxD=z;}
       if(dd>packetmax){packetmax=dd;packetmaxN=N;packetmaxD=z;}
     }
     if(cv){coneviol++;if(!conefirst)conefirst=N;conelast=N;}
     if(pf){packetviol++;if(!packetfirst)packetfirst=N;packetlast=N;}
     int64_t mx=0;int md=0;for(int z=0;z<=maxdepth;z++)if(pt[z]>mx){mx=pt[z];md=z;}
     double r=(double)mx/st;if(r>pratio){pratio=r;pn=N;pd=md;}
     int fail[2]={2*mx>3*st,mx>st};for(int j=0;j<2;j++)if(fail[j]){pviol[j]++;if(!pfirst[j])pfirst[j]=N;}
   }
   if((N&(N-1))==0){printf("%s{\"N\":%lld,\"source\":",sep++?",":"",(long long)N);arr(s);printf(",\"arrival\":");arr(g);printf(",\"clock\":");arr(h);if(prefix){printf(",\"prefix_actual\":[");for(int z=0;z<=maxdepth;z++)printf("%s%lld",z?",":"",(long long)pe[z]);printf("],\"prefix_trace\":[");for(int z=0;z<=maxdepth;z++)printf("%s%lld",z?",":"",(long long)pt[z]);printf("],\"packet_source_energy\":%lld",(long long)(s[0]-2*sib));}printf("}");}
 }
 printf("],\"octaves\":[");
 for(int k=0;k<power;k++)printf("%s{\"k\":%d,\"max_trace_ratio\":%.12g,\"N_ratio\":%lld,\"max_trace_per_N\":%.12g,\"N_trace\":%lld,\"max_clock_trace_per_N\":%.12g,\"N_clock\":%lld}",k?",":"",k,maxratio[k],(long long)nratio[k],maxt[k],(long long)nte[k],maxc[k],(long long)ntc[k]);
 printf("],\"violations\":[%lld,%lld,%lld],\"first_violations\":[%lld,%lld,%lld],\"max_clock\":%d",(long long)viol[0],(long long)viol[1],(long long)viol[2],(long long)first[0],(long long)first[1],(long long)first[2],maxclock);
 printf(",\"prefix_test\":%s,\"max_depth\":%d,\"prefix_violations\":[%lld,%lld],\"prefix_first_violations\":[%lld,%lld],\"max_prefix_ratio\":%.12g,\"prefix_ratio_N\":%lld,\"prefix_ratio_depth\":%d,",prefix?"true":"false",maxdepth,(long long)pviol[0],(long long)pviol[1],(long long)pfirst[0],(long long)pfirst[1],pratio,(long long)pn,pd);
 printf("\"orientation\":{\"violations\":%lld,\"first\":%lld,\"last\":%lld,\"max_increment\":%lld,\"N\":%lld,\"depth\":%d},\"packet_bound\":{\"violations\":%lld,\"first\":%lld,\"last\":%lld,\"max_deficit\":%lld,\"N\":%lld,\"depth\":%d}}\n",(long long)coneviol,(long long)conefirst,(long long)conelast,(long long)conemax,(long long)conemaxN,conemaxD,(long long)packetviol,(long long)packetfirst,(long long)packetlast,(long long)packetmax,(long long)packetmaxN,packetmaxD);
 fprintf(stderr,"cone violations %lld first %lld last %lld; max %lld N %lld d %d; packet violations %lld first %lld last %lld; max %lld N %lld d %d\n",(long long)coneviol,(long long)conefirst,(long long)conelast,(long long)conemax,(long long)conemaxN,conemaxD,(long long)packetviol,(long long)packetfirst,(long long)packetlast,(long long)packetmax,(long long)packetmaxN,packetmaxD);
 free(pe);free(arrival);free(hist);free(port);free(clock);free(depth);free(coeff);free(vi);free(pv);free(pt);return 0;
}
