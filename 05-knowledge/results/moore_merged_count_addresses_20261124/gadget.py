from itertools import combinations

def profiles(n,mask):
 bits=list(combinations(range(n),2))
 A=[[0]*n for i in range(n)]
 for k,(i,j) in enumerate(bits):
  A[i][j]=1-((mask>>k)&1);A[j][i]=1-A[i][j]
 M=1<<n
 en=[[0]*n for s in range(M)];st=[[0]*n for s in range(M)]
 H=[0]*M;H[0]=1
 for s in range(1,M):
  for i in range(n):
   if s>>i&1:
    t=s^(1<<i)
    en[s][i]=sum(en[t][j]*A[j][i] for j in range(n)) if t else 1
    st[s][i]=sum(st[t][j]*A[i][j] for j in range(n)) if t else 1
  H[s]=sum(en[s])
 f=M-1
 Q=sum(H[s]*H[f^s] for s in range(M))
 left=[sum(en[s][i]*H[f^s] for s in range(M)) for i in range(n)]
 right=[sum(st[s][i]*H[f^s] for s in range(M)) for i in range(n)]
 cross=[[sum(en[s][i]*st[f^s][j] for s in range(M)) for j in range(n)] for i in range(n)]
 for b in range(M):
  I=[i for i in range(n) if not (b>>i)&1];O=[i for i in range(n) if b>>i&1]
  K=sum(left[i] for i in I)
  R=H[f]+sum(right[i] for i in O)
  ha=sum(en[f][i] for i in I)+sum(st[f][j] for j in O)+sum(cross[i][j] for i in I for j in O)
  yield H[f],Q,K,R,ha
