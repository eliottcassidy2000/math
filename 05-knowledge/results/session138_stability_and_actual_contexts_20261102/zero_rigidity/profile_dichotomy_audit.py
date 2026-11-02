"""Bounded whole-theorem controls, separate from the 15-core certificate."""
import sys,itertools,json,math,hashlib
from pathcover_probe import profile,induce
from reference_paths import c3_blowup,extend,H

def coeff_formula(k,kind,r):
 def p(m):return -2*m*(m**3+2*m)**(k-1) if kind=='A' else -4*m*m*(m**3+2*m)**(k-2)
 return sum((-1)**(r-j)*math.comb(r,j)*p(j) for j in range(r+1))
def outer(U):
 n=len(U);a=tuple(v|(1<<n) for v in U)+(1<<(n+1),(1<<n)-1)
 return a
rows=[]
for k in range(1,4):
 q=tuple(sum(1<<j for j in range(i+1,k)) for i in range(k));T=c3_blowup(q);f=profile(T);counts={'equal':0,'low':0,'high_H':0,'high_equal_H':0}
 for root in range(k):
  V=[i for i in range(3*k) if i!=3*root];B=induce(T,V)
  for t in itertools.product(range(4),repeat=k):
   S=0
   for i,n in enumerate(t):
    if i==root:
     for j in range(2):
      if n>>j&1:S|=1<<V.index(3*i+j+1)
    else:
     for j in range(n):S|=1<<V.index(3*i+j)
   U=extend(B,S);g=profile(U);d=[b-a for a,b in zip(f,g)];outside=[t[i] for i in range(k) if i!=root];mono=not any(t[i]>0 and t[j]<3 for i in range(k) for j in range(i+1,k) if i!=root and j!=root)
   canonical=t[root]==1 and all(t[i]==(0 if i<root else 3) for i in range(k) if i!=root)
   low=(t[root]==2 and all(t[i]==(0 if i<root else 3) for i in range(k) if i!=root)) or (t[root]==0 and all(t[i]==0 for i in range(root)) and mono) or (t[root]==3 and all(t[i]==3 for i in range(root+1,k)) and mono)
   if canonical:
    if any(d):raise RuntimeError('canonical full profile')
    label='equal'
   elif low:
    if max(d)>0 or d[0]>=0:raise RuntimeError('low comparison')
    kind='B' if any(v in (1,2) for v in outside) else 'A'
    if d!=[coeff_formula(k,kind,r) for r in range(1,3*k+1)]:raise RuntimeError('exact low polynomial')
    label='low'
   else:
    if min(d)<0 or not any(d) or not(d[0]>0 or(d[0]==0 and d[1]>0)):raise RuntimeError('high comparison')
    label='high_H' if d[0]>0 else 'high_equal_H'
   counts[label]+=1
   actual=H(outer(U));predicted=3*g[0]+g[1]
   if actual!=predicted:raise RuntimeError('literal outer context formula')
   if (actual==3*f[0]+f[1])!=canonical:raise RuntimeError('restricted strong zero rigidity')
   rows.append((k,root,t,label,actual))
 print(json.dumps({'k':k,'rows':k*4**k,'counts':counts,'all_literal_outer_contexts_match':True},sort_keys=True))
print('PASS '+json.dumps({'rows':len(rows),'sha256':hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest()},sort_keys=True))
