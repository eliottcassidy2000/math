#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#define N 1729
#define LIMIT (2u*N*N)
static uint8_t visited[(LIMIT+7)/8];
static uint32_t queue[LIMIT];
static inline int test(uint32_t i){return (visited[i>>3]>>(i&7))&1;}
static inline void mark(uint32_t i){visited[i>>3]|=1u<<(i&7);}
int main(void){
 uint32_t head=0,tail=0, s=25*N+34;
 mark(s);queue[tail++]=s;
 while(head<tail){
  uint32_t v=queue[head++], e=v/(N*N), x=(v/N)%N,y=v%N;
  uint32_t xq=(11*x+4*y)%N,yq=(18*x+7*y)%N;
  uint32_t xb=(5*x+2*y)%N,yb=(8*x+3*y)%N;
  uint32_t t[2]={e*N*N+xq*N+yq,(1-e)*N*N+xb*N+yb};
  for(int j=0;j<2;j++)if(!test(t[j])){mark(t[j]);queue[tail++]=t[j];}
 }
 size_t bad=0,missing=0;
 for(int e=0;e<2;e++)for(int x=0;x<N;x++)for(int y=0;y<N;y++){
  int prim=(x%7||y%7)&&(x%13||y%13)&&(x%19||y%19);
  int hit=test((uint32_t)e*N*N+x*N+y); if(hit&&!prim)bad++;if(!hit&&prim)missing++;
 }
 printf("Forward base BFS: reached=%u, nonprimitive=%zu, missing primitive=%zu, expected=%d\n",tail,bad,missing,2*48*168*360);
 return bad || missing || tail!=2*48*168*360;
}
