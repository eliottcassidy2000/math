"""Signed source cylinders, clocks and first-hit guards; no cycle completeness claim."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as Q
from itertools import product
import json
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(v,m):
 global CHECKS
 CHECKS+=1
 if not v:raise RuntimeError(m)
def val(n):
 need(n!=0,'nonzero numerator')
 n=abs(n);a=0
 while n%2==0:n//=2;a+=1
 return a
def step(n,b):
 a=val(3*n+b)
 return (3*n+b)//2**a,a
def cylinder(word,b):
 A=0;B=0;r=0
 for a in word:
  B=3*B+b*2**A;A+=a;r+=1
 # 2^A U^r(n)=3^r n+B. Endpoint must be odd: one class modulo 2^(A+1).
 modulus=2**(A+1);residue=(2**A-B)*pow(3**r,-1,modulus)%modulus
 return residue,modulus,A,B
def cycle_from_word(word,b):
 A=sum(word);B=0;prefix=0
 for a in word:B=3*B+b*2**prefix;prefix+=a
 den=2**A-3**len(word)
 return Q(B,den),A,den
def main(write=True):
 global CHECKS
 CHECKS=0
 rows=[]
 for K in range(2,15):
  counts=Counter((min(val(3*n+1),K),min(val(3*n-1),K)) for n in range(1,2**K,2))
  size=2**(K-1)
  for k in range(2,K):need(counts[(1,k)]==counts[(k,1)]==2**(K-1-k),'exact joint signed cylinder')
  need(counts[(1,K)]==counts[(K,1)]==1,'two capped tail cylinders')
  mean=Q(sum(a*c for (a,b),c in counts.items()),size);joint=Q(sum(a*b*c for (a,b),c in counts.items()),size)
  need(mean==2-Q(1,2**(K-1)),'capped geometric mean')
  need(joint-mean**2==-(1-Q(1,2**(K-1)))**2,'shared-source negative covariance')
  rows.append({'precision':K,'mean':str(mean),'covariance':str(joint-mean**2),'classes':len(counts)})
 cylinders=[]
 for b in (-1,1):
  for length in range(1,5):
   for word in product(range(1,4),repeat=length):
    residue,modulus,A,B=cylinder(word,b)
    need(residue%2==1,'odd cylinder source')
    for lift in range(3):
     n=residue+lift*modulus;x=n;actual=[]
     for _ in word:x,a=step(x,b);actual.append(a)
     need(tuple(actual)==word and 2**A*x==3**length*n+B,'source-preserving word and clock')
    cylinders.append({'sheet':b,'word':word,'source_modulus':modulus,'source_residue':residue,'valuation_time':A,'affine_intercept':B})
 known=[(1,(2,)),(-1,(1,)),(-1,(1,2)),(-1,(1,1,1,2,1,1,4))]
 cycles=[]
 for b,word in known:
  n,A,den=cycle_from_word(word,b);need(n.denominator==1 and n>0,'supplied cycle integer')
  x=int(n);states=[]
  for a in word:
   states.append(x);x,aa=step(x,b);need(aa==a,'cycle valuation word')
  need(x==n and len(set(states))==len(states),'supplied cycle closes simply')
  need(step(-int(n),-b)[0]==-step(int(n),b)[0],'signed forward conjugacy')
  cycles.append({'sheet':b,'states':states,'odd_steps':len(word),'shortcut_steps':A,'denominator':den})
 need(step(3,1)[0]==step(7,-1)[0]==5,'same endpoint control')
 need(step(5,1)[0]==1 and step(5,-1)[0]==7,'sheet needed after shared endpoint')
 offspring=[((0,0),Q(1,3)),((0,3),Q(1,3)),((3,0),Q(1,3))]
 need(Q(1,4)/(1-Q(1,4))==Q(1,3),'even-valuation probability')
 need(sum(a*p for (a,b),p in offspring)==sum(b*p for (a,b),p in offspring)==1,'critical signed marginal offspring')
 need(sum(a*b*p for (a,b),p in offspring)-1==-1,'shared-source offspring covariance')
 out={'status':'FINITE-EXACT; Haar cylinder law, not a positive-integer orbit law',
      'joint_first_step':rows,'word_cylinders':cylinders,'supplied_cycles':cycles,
      'joint_signed_offspring':{'states':[s for s,p in offspring],'probabilities':[str(p) for s,p in offspring],
                                'covariance':-1,'scope':'same-source first step; not a coupled native tree theorem'},'checks':CHECKS}
 if write:(HERE/'arithmetic.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({'checks':CHECKS,'cycles':cycles,'last_covariance':rows[-1]}))
 return out
if __name__=='__main__':main()
