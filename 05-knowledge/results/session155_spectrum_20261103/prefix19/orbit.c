/* Exact forward and reverse BFS certificate for modulus 7*13*19. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define REQUIRE(c) do { if (!(c)) { fprintf(stderr,"verification failed at %d\n",__LINE__); exit(1); } } while (0)
static int imod(int a,int m){a%=m;return a<0?a+m:a;}
static int inv(int a,int m){for(int x=1;x<m;x++)if(a*x%m==1)return x;abort();}
int main(void){
 const int m=1729; size_t cap=(size_t)2*m*m;
 for(int reverse=0;reverse<2;reverse++){
  unsigned char *vis=calloc(cap,1); uint32_t *que=malloc(cap*sizeof(uint32_t)); REQUIRE(vis&&que);
  int t=inv(5,m); int A[2][4]={{11,4,18,7},{5,2,8,3}};
  if(reverse){int J[2][4]={{7*t%m,imod(-4*t,m),imod(-18*t,m),11*t%m},{imod(-3,m),2,8,imod(-5,m)}};for(int i=0;i<2;i++)for(int j=0;j<4;j++)A[i][j]=J[i][j];}
  size_t read=0,write=0; uint32_t initial=25*m+34;que[write++]=initial;vis[initial]=1;
  while(read<write){uint32_t h=que[read++];int y=h%m,x=(h/m)%m,p=h/(m*m);
   for(int j=0;j<2;j++){int u=(A[j][0]*x+A[j][1]*y)%m,v=(A[j][2]*x+A[j][3]*y)%m;uint32_t zz=((p^j)*m+u)*m+v;if(!vis[zz]){REQUIRE(write<cap);vis[zz]=1;que[write++]=zz;}}
  }
  for(int p=0;p<2;p++)for(int x=0;x<m;x++)for(int y=0;y<m;y++){
   int primitive=(x%7||y%7)&&(x%13||y%13)&&(x%19||y%19);REQUIRE(vis[((p*m+x)*m+y)]==primitive);
  }
  REQUIRE(write==5806080);printf("%s BFS modulus %d: %zu states; exact primitive-pair predicate in both parities\n",reverse?"reverse":"forward",m,write);free(vis);free(que);
 }
}
