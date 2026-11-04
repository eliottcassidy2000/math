// deckhash.c: stream allpairs rows (h e B strong), group every n rows (one tournament, roots in vertex order),
// sort the n triples, hash the sorted deck with two independent 64-bit hashes, and count distinct decks.
// Output: tournaments, distinct decks (lower bound on the number of classes the root-response deck separates;
// equality with the class count = the deck is a complete invariant at that order, up to a 2^-128 hash collision risk).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef struct { u64 h, e, B; } tri;
static int cmp(const void*a,const void*b){ const tri*x=a,*y=b; if(x->h!=y->h) return x->h<y->h?-1:1; if(x->e!=y->e) return x->e<y->e?-1:1; if(x->B!=y->B) return x->B<y->B?-1:1; return 0; }
static u64 fnv(const void*p,size_t len,u64 seed){ const unsigned char*s=p; u64 h=seed^0xcbf29ce484222325ULL; for(size_t i=0;i<len;i++){ h^=s[i]; h*=0x100000001b3ULL; h^=h>>29; } return h; }
typedef struct { u64 a,b; } key;
static key *tab; static size_t cap, cnt;
static int insert(key k){ size_t i=(k.a^(k.b*0x9E3779B97F4A7C15ULL))&(cap-1); while(tab[i].a||tab[i].b){ if(tab[i].a==k.a&&tab[i].b==k.b) return 0; i=(i+1)&(cap-1);} tab[i]=k; cnt++; return 1; }
int main(int argc,char**argv){ int n=atoi(argv[1]); cap=1ULL<<25; tab=calloc(cap,sizeof(key)); tri d[16]; int k=0; u64 tours=0; int st;
  while(scanf("%llu %llu %llu %d",&d[k].h,&d[k].e,&d[k].B,&st)==4){ k++; if(k==n){ qsort(d,n,sizeof(tri),cmp); key K={fnv(d,sizeof(tri)*n,1),fnv(d,sizeof(tri)*n,0x5bd1e995ULL)}; if(!K.a&&!K.b) K.a=1; insert(K); tours++; k=0; } }
  printf("order %d: %llu tournaments, %zu distinct root-response decks\n",n,tours,cnt); return 0; }
