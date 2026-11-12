/* checks that the listed order-8 tournaments are tournaments and pairwise non-isomorphic (canonical form = min over all 8! relabellings) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned int u32;
static int perms[40320][8]; static int np=0;
static void gen(int*p,int k,int used){ if(k==8){ memcpy(perms[np++],p,sizeof(int)*8); return;} for(int i=0;i<8;i++) if(!((used>>i)&1)){ p[k]=i; gen(p,k+1,used|1<<i);} }
static int cmp(const void*a,const void*b){ u32 x=*(u32*)a,y=*(u32*)b; return x<y?-1:x>y; }
int main(void){ int p[8]; gen(p,0,0); static u32 can[7000]; int cnt=0; int mk[8];
  while(scanf("%d %d %d %d %d %d %d %d",&mk[0],&mk[1],&mk[2],&mk[3],&mk[4],&mk[5],&mk[6],&mk[7])==8){
    for(int a=0;a<8;a++) for(int b=a+1;b<8;b++) if(((mk[a]>>b)&1)+((mk[b]>>a)&1)!=1){ printf("bad tournament line %d\n",cnt); return 1;}
    u32 best=0xffffffff;
    for(int q=0;q<np;q++){ int*P=perms[q]; u32 code=0; int k=0; for(int a=0;a<8;a++) for(int b=a+1;b<8;b++){ if((mk[P[a]]>>P[b])&1) code|=1u<<k; k++; } if(code<best) best=code; }
    can[cnt++]=best; }
  qsort(can,cnt,sizeof(u32),cmp); int dist=cnt>0; for(int i=1;i<cnt;i++) if(can[i]!=can[i-1]) dist++;
  printf("lines %d distinct classes %d\n",cnt,dist); return 0; }
