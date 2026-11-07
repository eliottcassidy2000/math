from fractions import Fraction as Q
import json

def g(n):
 if not n:return 1
 u=n//2
 return g(u)+(u+(u&1) if n%2==0 else 0)

def block(w,u):
 v=u; inc=0
 for bit in w:
  if bit=='0':inc+=v+(v&1)
  v=2*v+int(bit)
 return v,inc
checks=0; drift_min=None
for k in range(2,11):
 B=2**k
 for r in range(1,B-1):
  w=f'{r:0{k}b}'; p=r&1
  _,K=block(w,p); _,other=block(w,p+2); C=(other-K)//2; K-=C*p
  drift=Q(K)-Q(C*r,B-1)
  assert drift>0
  # independent virtual cyclic point, sum zero-run contributions
  x=-Q(r,B-1); parity=p; ds=Q(0)
  for bit in w:
   if bit=='0':ds+=x+parity
   x=2*x+int(bit);parity=int(bit)
  assert x==-Q(r,B-1) and ds==drift
  drift_min=drift if drift_min is None else min(drift_min,drift)
  for u in (0,1,2,7):
   vals=[];v=u
   for m in range(6):
    vals.append(g(v)); v,_=block(w,v)
   # after first iteration parity stable: second differences satisfy
   # d_{m+1}-B*d_m=(1-B)*drift
   for m in range(1,4):
    assert vals[m+2]-(B+1)*vals[m+1]+B*vals[m]==(1-B)*drift
    checks+=1
for m in range(1,101):
 v=2*(4**m-1)//3
 assert 9*g(v)==4**(m+1)+6*m+5
print(json.dumps({'status':'PASS','mixed_words_lengths':'2..10','mixed_word_count':sum(2**k-2 for k in range(2,11)),'exact_recurrence_checks':checks,'smallest_drift':str(drift_min),'alternating_checks':100},indent=2))
