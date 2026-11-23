from itertools import product

def mat(n,m,f):return tuple(tuple(f(i,j) for j in range(m)) for i in range(n))
def mul(a,b,p):
 return mat(len(a),len(b[0]),lambda i,j:sum(a[i][k]*b[k][j] for k in range(len(b)))%p)
def inv(a,p):
 n=len(a);R=[list(x)+[int(i==j) for j in range(n)] for i,x in enumerate(a)]
 for i in range(n):
  j=next(j for j in range(i,n) if R[j][i]%p)
  R[i],R[j]=R[j],R[i];u=pow(R[i][i]%p,-1,p)
  R[i]=[u*x%p for x in R[i]]
  for j in range(n):
   if j!=i:
    u=R[j][i];R[j]=[(a-u*b)%p for a,b in zip(R[j],R[i])]
 return tuple(tuple(x[n:]) for x in R)
def matrices(p):
 d=p//2
 def w(s):return mat(d,d,lambda i,j:((i+1)**2*(i==j)+(i+1)*(int(j==i+1)-int(j==i-1)+s*(i==j==d-1)))%p)
 P=mat(d,d,lambda i,j:(i+1)**2*(i==j)%p)
 BP=mat(2*d,2*d,lambda i,j:(i%d+1)**2*(i==j)%p)
 we,wo=w(1),w(-1)
 BW=mat(2*d,2*d,lambda i,j: (we[i][j] if i<d and j<d else wo[i-d][j-d] if i>=d and j>=d else -d*(i==d-1 and j==2*d-1) if i<d else 0)%p)
 C=mat(2*d,d,lambda i,j:((((pow(i+1,2*(p-1),p*p)-1)//p)%p) if i<d else -2*pow(i-d+1,-1,p))*(i%d==j)%p)
 return ((P,BP),(we,BW)),C
def check(p):
 (P,W),C=matrices(p);d=p//2
 aa=[P[0],W[0]]; bi=[inv(P[1],p),inv(W[1],p)]
 piv={};words=[]
 for w in ['PPPPWWWPPPWPPPP', 'WWPPWWWPPWPWWPWWPWW', 'WPWPWPWPPWWPWPWWW', 'WPWPPPWWW', 'PWWPPPPWWPWPPWPP', 'WPPWPWWPWPPPPPPPWPWPPWW', 'WPP', 'PWW', 'WWPWWPPPWWPPWWPPWWPWWPWWPPWP', 'PPP', 'WWPWWWWWPWPWWPPPWWWWWPPWWWPPP', 'WPPWPPWPWW', 'WPWWPPWPPPWPPPWP', 'PPPWWPWWWPPPWWP', 'PPPWPWPWWWPPPPPPPPWPPPPPWPPP', 'WWWPWWWPWPWPWWPP', 'PPPW', 'WPPWPPPWWPPWPWP', 'PWPWWPPPWWWPWPWWWPWPWWWP', 'PWPPWPWPPWPWWPWPP', 'PPPWWWWPWWPPWP', 'PWWWPWPPWPPPPWWPPP', 'PPPPWWWWWP', 'WPPPPWW', 'PPPPPWPWPWPPPPWP', 'PPWWWPWWPPWWWWPPWP', 'PWPPW', 'WWPWPWPWPWWWWW']:
  F=C
  for s in reversed(w):
   i='PW'.index(s)
   F=mul(mul(bi[i],F,p),aa[i],p)
  x=[row[0] for row in F]
  for i in range(2*d):
   if not x[i]:continue
   if i in piv:
    x=[(a-x[i]*b)%p for a,b in zip(x,piv[i])]
   else:
    u=pow(x[i],-1,p);piv[i]=[a*u%p for a in x];words.append(w);break
  if len(piv)==2*d:break
 assert len(piv)==2*d
 print(p,len(piv),words)
if __name__=='__main__':
 for p in (5,11,13,17,19,29):check(p)
