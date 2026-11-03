#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define M 1183
#define N (M*M)
#define Q (2*N)
int main(){ uint8_t *seen=calloc(Q,1); uint32_t *q=malloc(Q*sizeof *q); int t=0; q[t++]=25*M+34; seen[q[0]]=1; int iv=0; while((5*iv)%M!=1)iv++;
 for(int i=0;i<t;i++){int p=q[i]/N,u=q[i]%N,y=u/M,x=u%M; int a=((7*y-4*x)%M+M)%M*iv%M,b=((-18*y+11*x)%M+M)%M*iv%M;int d=((-3*y+2*x)%M+M)%M,e=((8*y-5*x)%M+M)%M;int ids[]={p*N+a*M+b,(1-p)*N+d*M+e};for(int j=0;j<2;j++)if(!seen[ids[j]]){seen[ids[j]]=1;q[t++]=ids[j];}}
 int c[2]={0};int hit[M]={0},outside=0;
 for(int i=0;i<t;i++){int p=q[i]/N,u=q[i]%N,y=u/M,x=u%M;c[p]++; if(!p)hit[(y+6)%M]=1; if((y%7==0&&x%7==0)||(y%13==0&&x%13==0))outside++;}
 int h=0;for(int i=0;i<M;i++)h+=hit[i]; printf("reverse orbit=%d, even=%d, odd=%d, covered=%d, nonprimitive=%d, expected each=%d\n",t,c[0],c[1],h,outside,48*(169*169-13*13));free(q);free(seen);
}
