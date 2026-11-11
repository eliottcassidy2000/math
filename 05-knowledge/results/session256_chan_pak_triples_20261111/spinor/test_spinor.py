from itertools import product
import json
O=((1,0,0,0),(0,2,0,1),(0,0,0,1),(0,0,0,1))
I=((1,1,1,0),(0,2,0,0),(0,0,0,0),(0,0,0,1))
def mul(A,B,q=None):
 R=tuple(tuple(sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)) for i in range(4))
 return tuple(tuple(v%q for v in row) for row in R) if q else R
def maps(X,Y,a,q=None):
 if a=='A': return mul(Y,X,q),Y
 if a=='C': return X,mul(mul(Y,X,q),X,q)
 return mul(Y,X,q),mul(mul(X,Y,q),X,q)
def matrices(addr,q=None):
 X,Y=O,I
 for a in reversed(addr): X,Y=maps(X,Y,a,q)
 return mul(Y,X,q)
def output(M): return M[0][0]+M[0][3]
subs={'A':{'O':'OI','I':'I'},'B':{'O':'OI','I':'OIO'},'C':{'O':'O','I':'OOI'}}
def word(addr):
 w='OI'
 for a in addr:w=''.join(subs[a][c] for c in w)
 return w
def spin(addr):
 m,n=2,1
 for a in addr:
  m,n={'A':(2*m-n,m),'B':(2*m+n,m),'C':(m+2*n,n)}[a]
 return m,n
def count(w):
 h,p,b=1,0,0
 for c in w:
  if c=='O':p,b=2*p+1,1
  else:h,p,b=h+p+b,2*p,0
 return h
# Direct permutation subset DP, independent of the cut formula.
def tournament_h(w):
 n=len(w)+1; apex=n-1; rows=[0]*n
 for i in range(apex):
  for j in range(i+1,apex):rows[i]|=1<<j
  if w[i]=='O':rows[apex]|=1<<i
  else:rows[i]|=1<<apex
 D=[{} for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for mask in range(1,1<<n):
  for i,v in D[mask].items():
   poss=rows[i]&~mask
   while poss:
    b=poss&-poss;poss-=b;j=b.bit_length()-1
    d=D[mask|b];d[j]=d.get(j,0)+v
 return sum(D[-1].values())
num=direct=0;examples=[]
for k in range(7):
 for t in product('ABC',repeat=k):
  a=''.join(t);w=word(a);m,n=spin(a);h=count(w)
  assert (len(w),w.count('I'))==(m,n)
  assert output(matrices(a))==h
  for q in (3,5,8,15):assert output(matrices(a,q))%q==h%q
  if len(w)<=13:assert tournament_h(w)==h;direct+=1
  if k<=2:examples.append({'address':a,'m':m,'n':n,'word':w,'H':h})
  num+=1
# All cuts, not only morphic words, exercise the newly proved formula.
allcuts=0
for k in range(9):
 for w in map(''.join,product('OI',repeat=k)):
  assert tournament_h(w)==count(w);allcuts+=1
print(json.dumps({'addresses_through_depth_6':num,'direct_morphic_graphs':direct,'all_cuts_through_length_8':allcuts,'examples':examples},indent=2))
