/* Exact bulk scan for disconnected high edges with gcd<=3 and raw p<264.

   Every supplied context is body-safe, so all nominal first-clause teeth are
   full and the boundary correction in THM-3352 is identically zero.  This is
   a literal integer port of the general floor-moment/triangle-sum engine.
*/
#include <errno.h>
#include <limits.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;
typedef __int128 i128;
typedef struct { i128 m0,m1,m2; } Moments;
typedef struct { int L,j,e,f; } Context;
typedef struct { int g,P,Q; } Task;
typedef struct { i64 num,den; int ci; } Result;

static const i64 TARGET_NUM=186636088362LL;
static const i64 TARGET_DEN=58865718786875LL;
static const char EXPECTED_CONTEXT_SHA256[]="efea6bd97522fe1c31a5a88ca9f3223f9e7a8c08e3be85c493e9f62fdfaf06e4";
static const char EXPECTED_TASK_SHA256[]="2771f2f901f2f052952343fd77412114ae1d1d99543f42539bc28d0f0f1948af";
static const char REFERENCE_ENGINE_SHA256[]="afd417297131401254769e1ef172d89c109ad2f9a843ea55e2badc3e7891435b";
static const char CONTEXT_COMPILER_SHA256[]="aebfe98ab72f7eb0dc1718dfb17529e5b3f9288c6ed97d57f048bf3b29281f19";
static Context *contexts; static int ncontexts;
static Task *tasks; static Result *results; static int ntasks;
static volatile int next_task=0;

static void die(const char *s){ perror(s); exit(2); }
static i64 gcd64(i64 a,i64 b){ if(a<0)a=-a; while(b){i64 t=a%b;a=b;b=t;}return a; }
static i128 abs128(i128 x){return x<0?-x:x;}
static i128 gcd128(i128 a,i128 b){a=abs128(a);b=abs128(b);while(b){i128 t=a%b;a=b;b=t;}return a;}
static void print128(FILE *f,i128 x){
 char b[64];int n=0;if(x==0){fputc('0',f);return;}if(x<0){fputc('-',f);x=-x;}
 while(x){b[n++]=(char)('0'+x%10);x/=10;}while(n)fputc(b[--n],f);
}

static Moments floor_moments(i64 n,i64 m,i64 a,i64 b){
 if(n==0)return (Moments){0,0,0};
 i128 s1=(i128)n*(n-1)/2, s2=(i128)n*(n-1)*(2*n-1)/6;
 i64 qa=a/m,a0=a%m,qb=b/m,b0=b%m;
 i128 base0=(i128)qa*s1+(i128)qb*n;
 i128 base1=(i128)qa*s2+(i128)qb*s1;
 i128 base2=(i128)qa*qa*s2+(i128)2*qa*qb*s1+(i128)qb*qb*n;
 if(a0==0)return (Moments){base0,base1,base2};
 i64 height=(a0*(n-1)+b0)/m;
 if(height==0)return (Moments){base0,base1,base2};
 Moments u=floor_moments(height,a0,m,m-b0+a0-1);
 i128 r0=(i128)n*height-u.m0;
 i128 r1=(i128)height*s1-(u.m2-u.m0)/2;
 i128 r2=(i128)n*height*height-2*u.m1-u.m0;
 return (Moments){base0+r0,base1+r1,base2+2*qa*r1+2*qb*r0+r2};
}

static void residue_prefix(i64 n,i64 m,i64 a,i64 b,i64 threshold,
                           Moments base,i128 total,i128 *count,i128 *sum){
 if(threshold<=0){*count=0;*sum=0;return;}
 if(threshold>=m){*count=n;*sum=total;return;}
 Moments sh=floor_moments(n,m,a,b+m-threshold);
 i128 d0=sh.m0-base.m0,d1=sh.m1-base.m1;
 i128 y0d=(sh.m2-base.m2-d0)/2;
 i128 high_sum=(i128)a*d1+(i128)b*d0-(i128)m*y0d;
 *count=n-d0;*sum=total-high_sum;
}

static i128 triangle_sum(i64 n,i64 m,i64 a,i64 b,i64 peak,i64 L,
                         Moments base,i128 total){
 if(peak<=0)return 0;
 i64 radius=(peak-1)/L,turns=radius/m,tail=radius%m;
 i128 answer=(i128)n*((i128)2*turns*peak-(i128)L*m*turns*turns);
 i128 lc,ls,bc,bs;
 residue_prefix(n,m,a,b,tail+1,base,total,&lc,&ls);
 answer+=((i128)peak-(i128)L*turns*m)*lc-(i128)L*ls;
 residue_prefix(n,m,a,b,m-tail,base,total,&bc,&bs);
 i128 hc=n-bc,hs=total-bs;
 answer+=((i128)peak-(i128)L*(turns+1)*m)*hc+(i128)L*hs;
 return answer;
}

static Result mass_one(Context c,int p,int q){
 int e=c.e,f=c.f;if((i64)c.L*p-e>(i64)c.L*q-f){int t=p;p=q;q=t;t=e;e=f;f=t;}
 i64 L=c.L,z=L*p-e,w=L*q-f;
 i64 r=(i64)e*c.j%L,s=(i64)f*c.j%L;
 i64 det=r*w-s*z;
 if(det%L){fprintf(stderr,"nonintegral phase\n");exit(3);}
 i64 b=(det/L)%z;if(b<0)b+=z;
 i64 a=w%z;
 Moments base=floor_moments(p,z,a,b);
 i128 total=(i128)a*p*(p-1)/2+(i128)b*p-(i128)z*base.m0;
 i64 unit=L/14,outer=unit*(z+w),inner=unit*(w-z);
 i128 num=triangle_sum(p,z,a,b,outer,L,base,total)
          -triangle_sum(p,z,a,b,inner,L,base,total);
 i128 den=(i128)z*w;
 if(num<0||num>den||num>LLONG_MAX||den<=0||den>LLONG_MAX){fprintf(stderr,"mass range failure\n");exit(3);}
 return (Result){(i64)num,(i64)den,0};
}

static int query_mode(const char *inpath,const char *outpath){
 FILE *in=fopen(inpath,"r");if(!in)die(inpath);
 FILE *out=fopen(outpath,"w");if(!out)die(outpath);
 int L,j,e,p,f,q;long long rows=0;
 while(fscanf(in,"%d %d %d %d %d %d",&L,&j,&e,&p,&f,&q)==6){
  int u=L/14,re=e*j%L,rf=f*j%L;
  if(L%14||L>=4592||e<1||e>14||f<1||f>14||p<1||p>2103||q<1||q>2103||
     !(u<=re&&re+e<=L-u&&u<=rf&&rf+f<=L-u)){
   fprintf(stderr,"unsafe query row %lld\n",rows);return 3;
  }
  Result x=mass_one((Context){L,j,e,f},p,q);i64 gg=gcd64(x.num,x.den);
  fprintf(out,"%lld %lld\n",x.num/gg,x.den/gg);rows++;
 }
 if(!feof(in)){fprintf(stderr,"malformed query input\n");return 3;}
 fclose(in);fclose(out);fprintf(stderr,"query_rows %lld\n",rows);return 0;
}

static void *worker(void *unused){(void)unused;
 for(;;){
  int ti=__sync_fetch_and_add(&next_task,1);if(ti>=ntasks)break;
  Task t=tasks[ti];Result best={-1,1,-1};
  int p=t.g*t.P,q=t.g*t.Q;
  for(int ci=0;ci<ncontexts;ci++){
   Result x=mass_one(contexts[ci],p,q);x.ci=ci;
   if(best.num<0 || (i128)x.num*best.den<(i128)best.num*x.den){best=x;}
  }
  results[ti]=best;
 }
 return NULL;
}

static void read_contexts(const char *path){
 FILE *f=fopen(path,"r");if(!f)die(path);int cap=4096;ncontexts=0;
 contexts=malloc((size_t)cap*sizeof(*contexts));if(!contexts)die("malloc contexts");
 Context c;while(fscanf(f,"%d %d %d %d",&c.L,&c.j,&c.e,&c.f)==4){
  if(c.L%14 || c.L>=4592){fprintf(stderr,"bad context\n");exit(3);}
  int u=c.L/14,re=c.e*c.j%c.L,rf=c.f*c.j%c.L;
  if(!(u<=re&&re+c.e<=c.L-u&&u<=rf&&rf+c.f<=c.L-u)){fprintf(stderr,"unsafe context\n");exit(3);}
  contexts[ncontexts++]=c;
 }fclose(f);if(ncontexts!=2530){fprintf(stderr,"context count %d\n",ncontexts);exit(3);}
}
static int coprime(int a,int b){while(b){int t=a%b;a=b;b=t;}return a==1;}
static void make_tasks(void){
 int cap=220000;tasks=malloc((size_t)cap*sizeof(*tasks));if(!tasks)die("malloc tasks");ntasks=0;
 for(int g=1;g<=3;g++)for(int P=1;g*P<264;P++)for(int Q=P+1;Q<8*P;Q++)
  if(P+Q>=8&&coprime(P,Q))tasks[ntasks++]=(Task){g,P,Q};
 if(ntasks!=201377){fprintf(stderr,"task count %d\n",ntasks);exit(3);}
 results=calloc((size_t)ntasks,sizeof(*results));if(!results)die("calloc results");
}

int main(int argc,char **argv){
 if(argc==4&&!strcmp(argv[1],"--query"))return query_mode(argv[2],argv[3]);
 if(argc!=4){fprintf(stderr,"usage: %s CONTEXTS LEDGER THREADS\n",argv[0]);return 2;}
 int nth=atoi(argv[3]);if(nth<1||nth>64)return 2;
 read_contexts(argv[1]);make_tasks();
 pthread_t *th=malloc((size_t)nth*sizeof(*th));if(!th)die("malloc threads");
 for(int i=0;i<nth;i++)if(pthread_create(&th[i],NULL,worker,NULL))die("pthread_create");
 for(int i=0;i<nth;i++)pthread_join(th[i],NULL);
 FILE *out=fopen(argv[2],"w");if(!out)die(argv[2]);
 Result global={-1,1,-1},byg[4]={{0}};int bygi[4]={-1,-1,-1,-1};int failures=0;
 for(int ti=0;ti<ntasks;ti++){
  Task t=tasks[ti];Result x=results[ti];Context c=contexts[x.ci];
  i128 clear=(i128)x.num*TARGET_DEN-(i128)TARGET_NUM*x.den;
  if(clear<=0)failures++;
  if(global.num<0||(i128)x.num*global.den<(i128)global.num*x.den){global=x;global.ci=ti;}
  if(bygi[t.g]<0||(i128)x.num*byg[t.g].den<(i128)byg[t.g].num*x.den){byg[t.g]=x;bygi[t.g]=ti;}
  i64 gg=gcd64(x.num,x.den);
  fprintf(out,"%d %d %d %lld %lld %d %d %d %d\n",t.g,t.P,t.Q,x.num/gg,x.den/gg,c.L,c.j,c.e,c.f);
 }
 fclose(out);
 printf("DISCONNECTED RAW HEAD EXACT BULK SCAN\ncontexts %d tasks %d threads %d failures %d\n",ncontexts,ntasks,nth,failures);
 printf("context_expected_sha256 %s\n",EXPECTED_CONTEXT_SHA256);
 printf("task_semantic_sha256 %s\n",EXPECTED_TASK_SHA256);
 printf("reference_engine_sha256 %s\n",REFERENCE_ENGINE_SHA256);
 printf("context_compiler_sha256 %s\n",CONTEXT_COMPILER_SHA256);
 printf("universe g=1..3;gP<264;gcd(P,Q)=1;P<Q<8P;P+Q>=8;all 2530 L<4592 contexts\n");
 for(int g=1;g<=3;g++){
  int ti=bygi[g];Task t=tasks[ti];Result x=results[ti];Context c=contexts[x.ci];i64 gg=gcd64(x.num,x.den);
  i128 mn=(i128)x.num*TARGET_DEN-(i128)TARGET_NUM*x.den,md=(i128)x.den*TARGET_DEN,mg=gcd128(mn,md);
  printf("g %d weakest_channel %d %d mass %lld/%lld context %d %d %d %d margin ",g,t.P,t.Q,x.num/gg,x.den/gg,c.L,c.j,c.e,c.f);print128(stdout,mn/mg);putchar('/');print128(stdout,md/mg);putchar('\n');
 }
 int ti=global.ci;Task t=tasks[ti];Result x=results[ti];Context c=contexts[x.ci];i64 gg=gcd64(x.num,x.den);
 printf("global_weakest %d %d %d mass %lld/%lld context %d %d %d %d\n",t.g,t.P,t.Q,x.num/gg,x.den/gg,c.L,c.j,c.e,c.f);
 return failures?1:0;
}
