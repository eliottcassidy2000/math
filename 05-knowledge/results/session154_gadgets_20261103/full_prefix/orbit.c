#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#define REQUIRE(c) do { if (!(c)) { fprintf(stderr, "verification failed at %d\n", __LINE__); exit(1); } } while (0)
#define M 1183
#define N (M*M)
int main(void){
 uint8_t *seen=calloc(2*N,1); int32_t *queue=malloc((size_t)2*N*sizeof(int32_t)); int32_t *parent=malloc((size_t)2*N*sizeof(int32_t)); uint8_t *letter=malloc((size_t)2*N); int32_t head=0,tail=0; int start=25*M+34; seen[start]=1;parent[start]=-1;queue[tail++]=start;
 while(head<tail){int k=queue[head++], p=k/N,v=k%N,Y=v/M,X=v%M;
  int Yq=(11*Y+4*X)%M,Xq=(18*Y+7*X)%M;
  int Yb=(5*Y+2*X)%M,Xb=(8*Y+3*X)%M;
  int a=p*N+Yq*M+Xq,b=(1-p)*N+Yb*M+Xb;
  if(!seen[a])seen[a]=1,parent[a]=k,letter[a]='q',queue[tail++]=a;
  if(!seen[b])seen[b]=1,parent[b]=k,letter[b]='b',queue[tail++]=b;
 }
 int valid=0; for(int p=0;p<2;p++)for(int y=0;y<M;y++)for(int x=0;x<M;x++){
  int expected=((y%7!=0 || x%7!=0) && (y%13!=0 || x%13!=0));
  REQUIRE(seen[p*N+y*M+x]==expected); valid+=expected;
 }
 REQUIRE(valid==2725632 && tail==valid);
 int witness=-1;for(int i=0;i<tail;i++){int k=queue[i],p=k/N,y=(k%N)/M;if(p==0 && (3*(y+6)+2)%M==7){witness=k;break;}}
 REQUIRE(witness>=0); char path[10000]; int plen=0;for(int k=witness;parent[k]!=-1;k=parent[k])path[plen++]=letter[k];
 printf("witness offset2 residue7 length=%d forward-transition-word=",plen);for(int i=plen-1;i>=0;i--)putchar(path[i]);putchar('\n');
 int count[2]={0,0}; int hits[2][M];memset(hits,0,sizeof hits);
 for(int k=0;k<tail;k++){int id=queue[k],p=id/N,Y=(id%N)/M;count[p]++; if(!p){int h=(Y+6)%M;hits[0][(3*h+2)%M]=1;hits[1][(3*h+4)%M]=1;}}
 printf("orbit=%d parity_even=%d odd=%d\n",tail,count[0],count[1]);
 for(int i=0;i<2;i++){int z=0,a=0;for(int j=0;j<M;j++)z+=hits[i][j];printf("offset %d distinct=%d missing:",i?4:2,z); for(int j=0;j<M;j++)if(!hits[i][j])printf(" %d",j),a++;printf(" (total %d)\n",a);}
 int distinct=0;for(int j=0;j<M;j++)distinct+=hits[0][j]||hits[1][j];printf("union=%d\n",distinct);
 for(int i=0;i<2;i++){int c=0; printf("lane7 offset %d missing quotients mod169:",i?4:2);for(int j=0;j<169;j++)if(!hits[i][7*j])printf(" %d",j),c++;printf(" (total %d)\n",c);}
 int h=25,x=34;long per=0; int only[2][M];memset(only,0,sizeof only);do{int a=(3*(h+6)+2)%M,b=(3*(h+6)+4)%M;only[0][a]=1;only[1][b]=1;int ny=(11*h+4*x)%M,nx=(18*h+7*x)%M;h=ny;x=nx;per++;}while(h!=25||x!=34); printf("q2 orbit period %ld\n",per);for(int i=0;i<2;i++){int z=0,c=0; for(int j=0;j<M;j++)z+=only[i][j];printf("q2 offset %d distinct %d missing lane7 quotients:",i?4:2,z);for(int j=0;j<169;j++)if(!only[i][7*j])printf(" %d",j),c++; printf(" (total %d)\n",c);}
 free(seen);free(queue);free(parent);free(letter);return 0;
}
