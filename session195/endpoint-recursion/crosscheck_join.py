from itertools import combinations

def data(rows):
 n=len(rows);N=1<<n
 dp=[[0]*n for _ in range(N)];sp=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for v in range(n):dp[1<<v][v]=sp[1<<v][v]=1
 for S in range(1,N):
  for v in range(n):
   if S>>v&1 and S!=(1<<v):
    M=S^(1<<v)
    dp[S][v]=sum(dp[M][u] for u in range(n) if M>>u&1 and rows[u]>>v&1)
    sp[S][v]=sum(sp[M][u] for u in range(n) if M>>u&1 and rows[v]>>u&1)
  H[S]=sum(dp[S])
 L=[sum(dp[S][v]*H[(N-1)^S] for S in range(N)) for v in range(n)]
 R=[sum(sp[S][v]*H[(N-1)^S] for S in range(N)) for v in range(n)]
 return H,dp[-1],sp[-1],L,R

def tournaments(n):
 pairs=list(combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  rows=[0]*n
  for k,(i,j) in enumerate(pairs):
   if mask>>k&1:rows[i]|=1<<j
   else:rows[j]|=1<<i
  yield rows

def check():
 bank=[]
 for n in range(1,4):
  for T in tournaments(n):bank.append((T,data(T)))
 checks=0
 for A,da in bank:
  for G,dg in bank:
   na=len(A);ng=len(G);HA,ea,sa,la,ra=da;HG,eg,sg,lg,rg=dg
   for a in range(na):
    for b in range(ng):
     rows=[row|(((1<<ng)-1)<<na) for row in A]+[row<<na for row in G]
     rows[a]^=1<<(na+b);rows[na+b]|=1<<a
     D=data(rows);S=((1<<ng)-1)^(1<<b)
     B=sum(HG[V]*HG[S^V] for V in range(1<<ng) if V&S==V)
     x=(HA[-1],ea[a],ra[a]);h,e,r=x
     pred=(HG[-1]*h-sg[b]*e+lg[b]*r,eg[b]*h-(ng==1)*e,rg[b]*h+B*r)
     actual=(D[0][-1],D[1][na+b],D[4][na+b])
     if pred!=actual:raise RuntimeError((A,G,a,b,pred,actual))
     checks+=1
 print({'all_labelled_A_and_G_orders':'1..3','rooted_pair_checks':checks,'status':'PASS'})
if __name__=='__main__':check()
