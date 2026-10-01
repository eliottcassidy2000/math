// procgen_edim_20261001_canon6.c -- canonical form of subsets of V(Q6) under Aut(Q6) (order 46080):
// minimum 64-bit image mask; also prints the setwise stabilizer order. Reads hex masks (or 'mask=HEX' lines).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static uint64_t T[720][8][256];
static const uint64_t MJ[6]={0x5555555555555555ULL,0x3333333333333333ULL,0x0F0F0F0F0F0F0F0FULL,0x00FF00FF00FF00FFULL,0x0000FFFF0000FFFFULL,0x00000000FFFFFFFFULL};
static inline uint64_t transl(uint64_t m,int t){ for(int j=0;j<6;j++) if(t>>j&1){int s=1<<j; m=((m&MJ[j])<<s)|((m>>s)&MJ[j]);} return m; }
int main(){ int perms[720][6],np=0; int q[6];
  for(q[0]=0;q[0]<6;q[0]++)for(q[1]=0;q[1]<6;q[1]++)for(q[2]=0;q[2]<6;q[2]++)for(q[3]=0;q[3]<6;q[3]++)for(q[4]=0;q[4]<6;q[4]++)for(q[5]=0;q[5]<6;q[5]++){
    int u=0,ok=1;for(int i=0;i<6;i++){if(u>>q[i]&1)ok=0;u|=1<<q[i];} if(ok) memcpy(perms[np++],q,sizeof q);}
  for(int p=0;p<720;p++){ int img[64]; for(int v=0;v<64;v++){int u=0;for(int i=0;i<6;i++)if(v>>i&1)u|=1<<perms[p][i];img[v]=u;}
    for(int b=0;b<8;b++)for(int val=0;val<256;val++){uint64_t r=0;for(int j=0;j<8;j++)if(val>>j&1)r|=1ULL<<img[b*8+j];T[p][b][val]=r;} }
  char line[512];
  while(fgets(line,sizeof line,stdin)){ char*p=strstr(line,"mask="); p=p?p+5:line; uint64_t m=strtoull(p,0,16);
    uint64_t best=~0ULL; int stab=0;
    for(int t=0;t<64;t++){ uint64_t mt=transl(m,t); for(int g=0;g<720;g++){ uint64_t x=0; for(int b=0;b<8;b++) x|=T[g][b][(mt>>(8*b))&255]; if(x<best)best=x; if(x==m)stab++; } }
    printf("%016llx %d %d\n",(unsigned long long)best,stab,__builtin_popcountll(m)); }
  return 0; }
