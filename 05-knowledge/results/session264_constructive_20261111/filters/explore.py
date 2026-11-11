from math import comb

def cycl(n,j,c):
 p=[0]*(n+1)
 for k in range(n-j+1):p[k]+=c*comb(n-j,k);p[k+j]+=c*comb(n-j,k)
 return p

def wrap(n,cs):
 N=n+3;p=[comb(N,i) for i in range(N+1)]
 def add(v):
  for i,x in enumerate(v):p[i]+=x
 for j,c in cs.items():add(cycl(N,j+1,2*c))
 for s in range(1,n+1,2):add(cycl(N,s+2,2*comb(n,s)))
 for k in range(0,n+1,2):
  d=comb(n,k)+2*sum(c*comb(n-j,k-j) for j,c in cs.items() if j<=k)
  add(cycl(N,k+3,2*d))
 for j,c in cs.items():
  for s in range(1,n-j+1,2):
   coeff=4*c*comb(n-j,s);left=N-j-1-s-2
   for t in range(left+1):
    for shift in [0,j+1,s+2,j+s+3]:p[t+shift]+=coeff*comb(left,t)
 return p

def cycles(w):
 n=len(w)
 return {j:sum(comb(b-a-1,j-2) for a in range(n) for b in range(a+1,n) if w[a]=='1' and w[b]=='0' and b-a-1>=j-2) for j in range(2,n+1,2)}

if __name__=='__main__':
 assert wrap(5,cycles('10010'))==[231]+[462]*7+[231]
 print('Affine one-wrap constructor: PASS')
