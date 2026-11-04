from itertools import product

def hp(rows):
 n=len(rows);dp=[{} for _ in range(1<<n)]
 if not n:return 1
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for u,c in dp[mask].items():
   nxt=rows[u]&~mask
   while nxt:
    b=nxt&-nxt;nxt-=b;v=b.bit_length()-1;dp[mask|b][v]=dp[mask|b].get(v,0)+c
 return sum(dp[-1].values())
def ext(R,a,c):
 n=len(R);S=list(R)+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:S[u]|=1<<n
  else:S[n]|=1<<u
  if c>>u&1:S[u]|=1<<(n+1)
  else:S[n+1]|=1<<u
 return S

def calc(m,n):
 A=[sum(1<<j for j in range(i+1,m+n)) for i in range(m+n)]
 B=[sum(1<<j for j in range(i+1,m)) for i in range(m)]
 B += [sum(1<<j for j in range(i+1,m+n))|((1<<m)-1) for i in range(m,m+n)]
 U=(1<<m)-1;V=((1<<n)-1)<<m
 vals=[]
 for x,y in product(range(4),repeat=2):
  a=(U if x&1 else 0)|(V if x&2 else 0);c=(U if y&1 else 0)|(V if y&2 else 0)
  d=hp(ext(A,a,c))-hp(ext(B,a,c));vals.append(d)
 return vals
if __name__=='__main__':
 for m,n in product(range(1,5),repeat=2):
  vals=calc(m,n);print(m,n,sorted(set(vals)),sorted(set(v%4 for v in vals)))
