#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
int main(int argc, char **argv){ FILE *fp=fopen(argv[1],"rb"); static unsigned char inP[100000];
  FILE *s=fopen(argv[2],"r"); long v; while (fscanf(s,"%ld",&v)==1) inP[v]=1; fclose(s);
  long long cnt[6]={0}, fc[6]={0}, minH[6], maxNF[6]; for(int i=0;i<6;i++){minH[i]=1LL<<40; maxNF[i]=-1;}
  int n=10; unsigned char rec[34];
  while (fread(rec,34,1,fp)==1){ unsigned long long code; uint32_t H; memcpy(&code,rec,8); memcpy(&H,rec+8,4);
    int out[10]={0}, b=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(code>>b&1) out[i]++; else out[j]++; b++; }
    int md=9; for(int i=0;i<n;i++){ int a=out[i], c=n-1-out[i]; if(a<md) md=a; if(c<md) md=c; }
    cnt[md]++; if(!inP[H]) fc[md]++; if(H<minH[md]) minH[md]=H; if(inP[H] && H>maxNF[md]) maxNF[md]=H; }
  printf("order 10: min degree -> classes, first creators, least H, largest non-FC H\n");
  for(int i=0;i<6;i++) if(cnt[i]) printf("   %d: %lld %lld %lld %lld\n", i, cnt[i], fc[i], minH[i], maxNF[i]);
  return 0; }
