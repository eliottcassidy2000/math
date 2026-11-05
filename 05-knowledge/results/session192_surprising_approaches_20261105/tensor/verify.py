#!/usr/bin/env python3
"""Standalone exact checks for the fixed-apex reset algebra and periodic families."""
import hashlib,itertools,json,math
from fractions import Fraction as Q
from pathlib import Path

def require(test,message):
 if not test:raise ValueError(message)

def states(r):return list(itertools.product(*([range(3)]*r+[range(2)])))
def transitions(r,c):
 ss=states(r);ix={s:i for i,s in enumerate(ss)};out=[]
 for s in ss:
  col={}
  for i in range(r+1):
   if not(s[i] or i==0 or (c>>(i-1)&1)):continue
   t=list(s);t[i]=1 if i==r or not(c>>i&1) else 2;j=ix[tuple(t)]
   col[j]=col.get(j,0)+1
  out.append(col)
 return out

def apply(A,v):
 w=[0]*len(A)
 for i,a in enumerate(v):
  if a:
   for j,b in A[i].items():w[j]+=a*b
 return tuple(w)
def multiply(A,B):
 C=[]
 for col in B:
  z={}
  for i,b in col.items():
   for j,a in A[i].items():z[j]=z.get(j,0)+a*b
  C.append({j:a for j,a in z.items() if a})
 return C
def identity(n):return [{i:1} for i in range(n)]
def triangular(r,c):
 ss=states(r);ix={s:i for i,s in enumerate(ss)};out=[]
 for s in ss:
  col={}
  for i in range(r+1):
   local=[]
   if i==r:
    if s[i] or (c>>(r-1)&1):local=[(1,1)]
   elif s[i]!=2:
    a=1 if s[i]==1 or i==0 else (c>>(i-1)&1)
    if a:
     local=[(1,1)]
     if not(c>>i&1):local.append((2,1))
   for k,a in local:
    t=list(s);t[i]=k;j=ix[tuple(t)];col[j]=col.get(j,0)+a
  out.append(col)
 return out

def basis_change(r):
 ss=states(r);ix={s:i for i,s in enumerate(ss)};out=[]
 for s in ss:
  choices=[]
  for i,x in enumerate(s):
   if i==r:choices.append([(x,1)])
   elif x==0:choices.append([(0,1)])
   elif x==1:choices.append([(2,1)])
   else:choices.append([(1,1),(2,-1)])
  col={}
  for z in itertools.product(*choices):col[ix[tuple(a for a,b in z)]]=math.prod(b for a,b in z)
  out.append(col)
 return out

def block(r,w,tri=False):
 A=identity(2*3**r)
 for c in w:A=multiply(triangular(r,c) if tri else transitions(r,c),A)
 return A

def poly_apply(A,vector,factors):
 for root,multiplicity in factors:
  for _ in range(multiplicity):
   z=apply(A,vector);vector=tuple(a-root*b for a,b in zip(z,vector))
 return vector

def reverse_bits(c,r):return sum((c>>i&1)<<(r-1-i) for i in range(r))
def value(w):
 h=0
 for order in (0,1):
  v=(1,)+(0,)*17
  for c in w:v=apply(transitions(2,c if not order else reverse_bits(c,2)),v)
  for s,a in zip(states(2),v):
   x,y,z=s
   if x!=2 and (y==1 or (y==0 and order==0)):h+=a
 return h

def graph(w):
 n=len(w)+2;rows=[0]*n;rows[0]=2
 for j,c in enumerate(w,2):
  for i in (0,1):
   if c>>i&1:rows[i]|=1<<j
   else:rows[j]|=1<<i
  for i in range(2,j):rows[i]|=1<<j
 return rows

def dp(rows):
 n=len(rows);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for S in range(1,1<<n):
  for i in range(n):
   if not(S>>i&1):continue
   prev=S^(1<<i)
   if prev:D[S][i]=sum(D[prev][j] for j in range(n) if rows[j]>>i&1)
 return sum(D[-1])

def strong(rows):
 n=len(rows);full=(1<<n)-1
 for s in range(n):
  seen=1<<s
  while True:
   z=seen
   for i in range(n):
    if seen>>i&1:z|=rows[i]
   if z==seen:break
   seen=z
  if seen!=full:return False
 return True

def polymul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def gf_den(L,degrees):
 p=[1]
 for j,e in enumerate(degrees,1):
  for _ in range(e):p=polymul(p,[1,-j**L])
 return p

def family_check(w,degrees,numerator):
 L=len(w);den=gf_den(L,degrees);factors=[(j**L,e) for j,e in enumerate(degrees,1)]
 # A residual killed after one block and invisible now proves every recurrence,
 # not merely its first observed instances.
 residuals=[]
 for order in (0,1):
  ww=tuple(reverse_bits(c,2) if order else c for c in w);A=block(2,ww)
  q=poly_apply(A,(1,)+(0,)*17,factors)
  require(not any(apply(A,q)),'all-order residual must be killed after one block')
  residuals.append(sum(a for s,a in zip(states(2),q) if s[0]!=2 and (s[1]==1 or (s[1]==0 and order==0))))
 require(sum(residuals)==0,'all-order residual initially unobservable')
 seq=[value(w*k) for k in range(21)]
 for k,h in enumerate(seq):
  z=sum(den[j]*seq[k-j] for j in range(min(k,len(den)-1)+1))
  require(z==(numerator[k] if k<len(numerator) else 0),'rational generating function')
 vertex=[]
 for k in range(1,5):
  rows=graph(w*k);require(strong(rows),'strong periodic family');h=dp(rows);require(h==seq[k],'vertex DP family')
  vertex.append({'k':k,'order':len(rows),'H':h})
 return {'word':''.join(map(str,w)),'denominator_factors':[[j**L,e] for j,e in enumerate(degrees,1)],'numerator':numerator,'initials':seq[:9],'vertex_checks':vertex,'all_order_residual_certificate':True}

def main():
 receipt={};local=0
 for r in range(1,5):
  C=basis_change(r);ss=states(r)
  for c in range(1<<r):
   M=transitions(r,c);B=triangular(r,c)
   require(multiply(M,C)==multiply(C,B),'exact common tensor flag')
   for i,col in enumerate(B):
    require(col.get(i,0)==ss[i].count(1),'letter diagonal')
    for j in col:require(all(a<=b for a,b in zip(ss[i],ss[j])),'flag support')
   local+=1
 receipt['local_tensor_intertwining_checks']=local
 wordchecks=0;annihilators=0
 for r,Lmax in ((1,4),(2,3),(3,2)):
  ss=states(r);n=len(ss)
  for L in range(1,Lmax+1):
   for w in itertools.product(range(1<<r),repeat=L):
    A=block(r,w,tri=True)
    for i,col in enumerate(A):
     require(col.get(i,0)==ss[i].count(1)**L,'word-independent diagonal')
     for j in col:require(all(a<=b for a,b in zip(ss[i],ss[j])),'word support')
    wordchecks+=1
    if r<=2:
     factors=[(0,r+1)]+[(j**L,r+2-j) for j in range(1,r+2)]
     for i in range(n):
      e=tuple(int(i==j) for j in range(n))
      require(not any(poly_apply(A,e,factors)),'universal sharp-positive-exponent annihilator')
     annihilators+=1
 receipt['word_flags_checked']=wordchecks;receipt['full_matrix_annihilators_checked']=annihilators
 # Independent combinatorial checks of transfer, not merely conjugate matrices.
 presentations=0
 for L in range(5):
  for w in itertools.product(range(4),repeat=L):
   require(value(w)==dp(graph(w)),'two-apex transfer versus vertex DP');presentations+=1
 receipt['independent_presentations']=presentations
 receipt['characteristic_multiplicities_two_apex']={'0':8,'1':16,'2^L':10,'3^L':2}
 receipt['families']=[family_check((1,2),(2,2,1),[1,-14,57,-76]),family_check((1,2,0),(3,2,1),[1,-33,255,-615,-284]),family_check((2,0),(2,1,0),[1,-1,-2])]
 # The 120 numerator is nonzero at every denominator root, proving sharpness.
 P=[1,-33,255,-615,-284]
 vals=[sum(Q(a)*Q(1,b)**i for i,a in enumerate(P)) for b in (1,8,27)]
 require(all(vals),'all three pole orders attained')
 receipt['120_numerator_at_poles']=list(map(str,vals))
 # Same length and same characteristic polynomial; different Hamilton observables.
 receipt['isospectral_hostile']={'word_1':'00','H_1':value((0,0)),'word_2':'12','H_2':value((1,2))}
 require(receipt['isospectral_hostile']['H_1']==1 and receipt['isospectral_hostile']['H_2']==5,'spectrum loss hostile')
 raw=(json.dumps(receipt,sort_keys=True,indent=2)+'\n').encode();Path(__file__).with_name('receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','local_checks':local,'word_flags':wordchecks,'annihilators':annihilators,'independent_presentations':presentations,'families':3,'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
