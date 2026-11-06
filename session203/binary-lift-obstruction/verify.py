#!/usr/bin/env python3
"""Finite controls for a proved format obstruction, not F coverage."""
import itertools,json,math
from counting import labelled,nine
from compiler import move,base

def require(ok,why):
 if not ok:raise RuntimeError(why)

def strong(a):
 n=len(a)
 for v in range(n):
  seen=1<<v
  while True:
   more=seen
   for i in range(n):
    if seen>>i&1:more|=a[i]
   if more==seen:break
   seen=more
  if seen!=(1<<n)-1:return False
 return True

def apply(p,w):
 for op in w:p=move(p,op)
 return p

def main():
 counts=dict(inputs=0,positive_inputs=0,factorial_checks=0,gate_checks=0,word_tests=0,large_lifts=0,dead_terminal_controls=0)
 for n in range(1,6):
  for a in labelled(n):
   counts['inputs']+=1;p=nine(a);h,s,e,L,R,d,C,b,T=p;B=base(p)
   require(h<=math.factorial(n) and L<=math.factorial(n),('permutation bound',a,p))
   if not s:continue
   counts['positive_inputs']+=1
   for c in range(1,11):
    if h<=c:require(B<=2*math.factorial(c),('h factorial',a,p,c));counts['factorial_checks']+=1
   if n==1:continue
   if strong(a):require(h-e>=n-1,('cycle cuts',a,p))
   g=C-3*L-h
   if g>=0:
    DJ=base(move(p,'J'))-2*B
    require(DJ>=h-e+s,('J bound',a,p,DJ))
    for c in range(1,11):
     if DJ<=c:require(B<=2*math.factorial(c),('J factorial',a,p,c));counts['factorial_checks']+=1
   for op in ('K1','K2'):
    D=base(move(p,op))-2*B
    require(D>=h,('gate bound',a,op,p,D));counts['gate_checks']+=1
    for c in range(1,11):
     if D<=c:require(B<=2*math.factorial(c),('gate factorial',a,op,p,c));counts['factorial_checks']+=1
 # Every state below is physical by construction, not an arbitrary cone vector.
 seed=(1,1,0,1,2,1,4,2,3);states={seed}
 front={seed}
 for depth in range(5):
  front={move(p,op)for p in front for op in ('O','Z','J')};states|=front
 states|={move(p,op)for p in list(states)for op in ('K1','K2')}
 words=[w for k in range(4)for w in itertools.product(('O','Z','J','K1','K2'),repeat=k)]
 for p in states:
  B=base(p);require(p[6]-3*p[3]-p[0]>=0,('grammar cone',p))
  for a in range(1,4):
   q=apply(p,('O',)*a)
   require(base(move(q,'Z'))==B+2**a*(p[4]+p[7]) and base(move(q,'Z'))>3*B,('OZ amplification',p,a))
  for w in words:
   q=apply(p,w);Bout=base(q);counts['word_tests']+=1
   for c in range(1,7):
    if B<=2*math.factorial(c) or abs(Bout-2*B)>c:continue
    counts['large_lifts']+=1
    require(len(w)>=1 and w[0]=='Z'and all(op=='O'for op in w[1:]),('large lift classification',p,w,c,q))
    if len(w)>1:
     require(base(move(q,'Z'))>2*Bout+c,('terminal O not dead',p,w,c,q));counts['dead_terminal_controls']+=1
    else:
     h,s,e,L,R,d,C,b,T=q;delta=s+(R-h+e-d)+(b-L-d)
     require(d-e==0 and 1<=delta<=c,('ray initialization',p,w,c,q,delta))
     x=q
     for k in range(5):
      require(base(x)==2**k*(Bout+delta)-delta,('ray',q,k,delta));x=move(x,'Z')
 require(counts['large_lifts']>0 and counts['dead_terminal_controls']>0,('missing controls',counts))
 return dict(status='PASS',scope='Finite controls of bounded-offset binary-lift obstruction only; no F counterexample or coverage theorem.',counts=counts,physical_states=len(states),word_length_limit=3)

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
