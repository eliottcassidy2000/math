/* Independent audit of THM-4756 (sierp session 2026-10-22).
 * For each tournament T read from stdin (gentourng upper-triangle format, one per line)
 * or enumerated (all labelled, mode "all n"), compute:
 *   F(T,x) by BRUTE FORCE over all n! orderings (if brute=1) or by subset DP;
 *   P[X] = H(T[X]) for all X (subset DP), ordered path-cover numbers c_m (set-partition DP);
 *   hc[C] = directed Hamiltonian cycles of T[C];
 *   W4 = sum over disjoint odd-cycle sets R of 4^|R| 2^(n-|V(R)|),  OCF  sum 2^|R|;
 * and checks:
 *   (1) F(T,x) == sum_m c_m (x-1)^(n-m)   exactly (integers)
 *   (2) sum_X H(T[X]) H(T[V\X]) == W4     exactly
 *   (3) OCF: H == sum_R 2^|R|
 *   (4) c_2 == ((-1)^n + 1) H  mod 3
 *   (5) F_k == H * A(n,k) mod 3 for all k     (A = Eulerian numbers, computed by recurrence)
 *   (6) F_k == C(n-1,k) mod 2
 *   (7) u_3 = c_3/6 == H + hc(V) (n odd), 0 (n even) mod 3   [the uncommitted 'second level']
 *   (8) also: fraction of T with 3 | H (THM-086 Theorem B "30-35%")
 * Usage:  ./a  brute|dp  all  n          (all labelled tournaments on n vertices)
 *         ./a  brute|dp  stdin n         (gentourng lines on stdin)
 *         ./a  dp  random n count seed   (random labelled tournaments)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef __int128 i128;
typedef long long ll;
#define MAXN 16
static int n;
static int adj[MAXN][MAXN]; /* adj[i][j]=1 iff i->j */
static ll A_eul[MAXN+1][MAXN+1];
static ll binom[MAXN+1][MAXN+1];
static ll Fk[MAXN];
static ll *P;      /* P[X] = number of directed Ham paths of T[X], P[0]=1 */
static ll *HC;     /* HC[X] = directed Ham cycles of T[X] (|X|>=3) */
static ll *pathend; /* pathend[X*n+v] = # paths with vertex set X ending at v */
static long long nfail[9];
static long long ntot = 0, n3H = 0;

static void eulerian(void){
  memset(A_eul,0,sizeof A_eul);
  A_eul[1][0]=1;
  for(int m=2;m<=MAXN;m++) for(int k=0;k<m;k++){
    ll a = (k<m-1? (ll)(k+1)*A_eul[m-1][k]:0) + (k>0? (ll)(m-k)*A_eul[m-1][k-1]:0);
    A_eul[m][k]=a; }
  for(int a=0;a<=MAXN;a++){binom[a][0]=1; for(int b=1;b<=a;b++) binom[a][b]=binom[a-1][b-1]+(b<=a-1?binom[a-1][b]:0);}
}

static void brute_F(void){
  int perm[MAXN]; for(int i=0;i<n;i++) perm[i]=i;
  memset(Fk,0,sizeof Fk);
  /* iterate permutations in lexicographic order */
  while(1){
    int f=0; for(int i=0;i+1<n;i++) f+=adj[perm[i]][perm[i+1]];
    Fk[f]++;
    int i=n-2; while(i>=0 && perm[i]>perm[i+1]) i--; if(i<0) break;
    int j=n-1; while(perm[j]<perm[i]) j--; int t=perm[i];perm[i]=perm[j];perm[j]=t;
    for(int a=i+1,b=n-1;a<b;a++,b--){t=perm[a];perm[a]=perm[b];perm[b]=t;}
  }
}

static void dp_F(void){
  /* f[S][v][k] */
  int N=1<<n;
  ll *f = calloc((size_t)N*n*n, sizeof(ll));
  for(int v=0;v<n;v++) f[((size_t)(1<<v)*n+v)*n+0]=1;
  for(int S=1;S<N;S++) for(int v=0;v<n;v++) if(S>>v&1){
    ll *fv=&f[((size_t)S*n+v)*n];
    for(int w=0;w<n;w++) if(!(S>>w&1)){
      int S2=S|(1<<w); int e=adj[v][w]; ll *fw=&f[((size_t)S2*n+w)*n];
      for(int k=0;k+e<n;k++) fw[k+e]+=fv[k];
    }}
  memset(Fk,0,sizeof Fk);
  for(int v=0;v<n;v++) for(int k=0;k<n;k++) Fk[k]+=f[((size_t)(N-1)*n+v)*n+k];
  free(f);
}

static void paths_cycles(void){
  int N=1<<n;
  memset(pathend,0,sizeof(ll)*(size_t)N*n);
  for(int v=0;v<n;v++) pathend[(size_t)(1<<v)*n+v]=1;
  for(int S=1;S<N;S++) for(int v=0;v<n;v++) if((S>>v&1) && pathend[(size_t)S*n+v]){
    ll c=pathend[(size_t)S*n+v];
    for(int w=0;w<n;w++) if(!(S>>w&1) && adj[v][w]) pathend[(size_t)(S|1<<w)*n+w]+=c; }
  P[0]=1; for(int S=1;S<N;S++){ll s=0; for(int v=0;v<n;v++) s+=pathend[(size_t)S*n+v]; P[S]=s;}
  /* Ham cycles of T[C]: paths starting at min(C) ending at u with u->min(C); computed by restricted DP */
  for(int S=0;S<N;S++) HC[S]=0;
  ll *g = calloc((size_t)N*n,sizeof(ll));
  for(int s=0;s<n;s++){
    /* paths starting at s using only vertices > s */
    for(int S=0;S<N;S++) for(int v=0;v<n;v++) g[(size_t)S*n+v]=0;
    g[(size_t)(1<<s)*n+s]=1;
    for(int S=1<<s;S<N;S++){ if(!(S>>s&1)) continue; if(S & ((1<<s)-1)) continue;
      for(int v=0;v<n;v++) if((S>>v&1) && g[(size_t)S*n+v]){
        ll c=g[(size_t)S*n+v];
        for(int w=s+1;w<n;w++) if(!(S>>w&1) && adj[v][w]) g[(size_t)(S|1<<w)*n+w]+=c; }
      if(__builtin_popcount(S)>=3){ ll h=0; for(int u=0;u<n;u++) if((S>>u&1) && u!=s && adj[u][s]) h+=g[(size_t)S*n+u]; HC[S]=h; }
    }
  }
  free(g);
}

static void process(void){
  int N=1<<n;
  paths_cycles();
  ll H=P[N-1];
  /* ordered path covers c_m via set-partition DP: U[m][S] unordered */
  int M=n;
  ll *U = calloc((size_t)(M+1)*N,sizeof(ll));
  U[0]=1; /* U[0][0]=1 */
  for(int m=1;m<=M;m++) for(int S=1;S<N;S++){
    int low=S&(-S); int rest=S^low; ll s=0;
    /* B contains low, B subset S */
    for(int sub=rest;;sub=(sub-1)&rest){ int B=sub|low; s+=P[B]*U[(size_t)(m-1)*N+(S^B)]; if(sub==0)break; }
    U[(size_t)m*N+S]=s; }
  ll c[MAXN+2]; ll fact=1;
  for(int m=1;m<=n;m++){ fact*=m; c[m]=fact*U[(size_t)m*N+N-1]; }
  /* (1) F == sum c_m (x-1)^(n-m): coefficient of x^k: sum_m c_m C(n-m,k)(-1)^(n-m-k) */
  int bad1=0;
  for(int k=0;k<n;k++){ i128 s=0; for(int m=1;m<=n;m++){ int j=n-m; if(k<=j){ i128 t=(i128)c[m]*binom[j][k]; if((j-k)&1) s-=t; else s+=t; } }
    if(s!=(i128)Fk[k]) bad1=1; }
  if(Fk[n-1]!=H) bad1=1;
  if(bad1) nfail[1]++;
  /* (2) sum_X P[X]P[V\X] vs W4 ; (3) OCF */
  i128 lhs=0; for(int X=0;X<N;X++) lhs+=(i128)P[X]*P[(N-1)^X];
  i128 *W4=calloc(N,sizeof(i128)); i128 *W2=calloc(N,sizeof(i128));
  W4[0]=1; W2[0]=1;
  for(int S=1;S<N;S++){ int low=S&(-S); int rest=S^low; i128 s4=0,s2=0;
    for(int sub=rest;;sub=(sub-1)&rest){ int C=sub|low; int pc=__builtin_popcount(C);
      if(pc>=3 && (pc&1) && HC[C]){ s4+=(i128)4*HC[C]*W4[S^C]; s2+=(i128)2*HC[C]*W2[S^C]; }
      if(sub==0)break; }
    W4[S]=s4; W2[S]=s2; }
  i128 rhs=0, ocf=0; for(int S=0;S<N;S++){ rhs+=W4[S]*((i128)1<<(n-__builtin_popcount(S))); ocf+=W2[S]; }
  if(lhs!=rhs) nfail[2]++;
  if(ocf!=(i128)H) nfail[3]++;
  /* (4) c_2 mod 3 */
  if(n>=2){ ll want = ((n%2==0)?2:0)*H; if(((c[2]-want)%3+3)%3!=0) nfail[4]++; }
  /* (5),(6) */
  int bad5=0,bad6=0;
  for(int k=0;k<n;k++){ if(((Fk[k]-H*(A_eul[n][k]%3))%3+3)%3!=0) bad5=1; if(((Fk[k]-binom[n-1][k])%2+2)%2!=0) bad6=1; }
  if(bad5) nfail[5]++; if(bad6) nfail[6]++;
  /* (7) u3 */
  if(n>=3){ ll u3=c[3]/6; ll hc=HC[N-1]; ll want=(n%2)?(H+hc):0; if(((u3-want)%3+3)%3!=0) nfail[7]++; }
  ntot++; if(H%3==0) n3H++;
  free(U); free(W4); free(W2);
}

static unsigned long long rs;
static unsigned long long rnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }

int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage\n");return 1;}
  int brute=!strcmp(argv[1],"brute");
  n=atoi(argv[3]);
  eulerian();
  P=malloc(sizeof(ll)<<n); HC=malloc(sizeof(ll)<<n); pathend=malloc(sizeof(ll)*((size_t)n<<n));
  if(!strcmp(argv[2],"all")){
    int E=n*(n-1)/2; for(long long code=0; code<(1LL<<E); code++){
      int b=0; for(int i=0;i<n;i++){adj[i][i]=0; for(int j=i+1;j<n;j++){ int bit=(code>>b)&1; b++; adj[i][j]=bit; adj[j][i]=1-bit; }}
      if(brute) brute_F(); else dp_F();
      process(); }
  } else if(!strcmp(argv[2],"stdin")){
    char line[1024];
    while(fgets(line,sizeof line,stdin)){ if(line[0]=='>'||line[0]=='\n') continue;
      int b=0; for(int i=0;i<n;i++){adj[i][i]=0; for(int j=i+1;j<n;j++){ int bit=line[b]-'0'; b++; adj[i][j]=bit; adj[j][i]=1-bit; }}
      if(brute) brute_F(); else dp_F();
      process(); }
  } else if(!strcmp(argv[2],"random")){
    long cnt=atol(argv[4]); rs=strtoull(argv[5],0,10)*2654435761ULL+12345;
    for(long t=0;t<cnt;t++){
      for(int i=0;i<n;i++){adj[i][i]=0; for(int j=i+1;j<n;j++){ int bit=rnd()&1; adj[i][j]=bit; adj[j][i]=1-bit; }}
      if(brute) brute_F(); else dp_F();
      process(); }
  }
  printf("n=%d mode=%s/%s tournaments=%lld | fails: (1)F=sum c_m(x-1)^(n-m):%lld (2)sumHH=W4:%lld (3)OCF:%lld (4)c2:%lld (5)F=H*A_n mod3:%lld (6)F=(1+x)^(n-1) mod2:%lld (7)u3:%lld | frac 3|H = %lld/%lld = %.4f\n",
    n,argv[1],argv[2],ntot,nfail[1],nfail[2],nfail[3],nfail[4],nfail[5],nfail[6],nfail[7],n3H,ntot,(double)n3H/ntot);
  return 0;
}
