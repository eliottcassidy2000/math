import json
from fractions import Fraction

def f(n):
 if not n:return 0
 q,r=divmod(n,2)
 return f(q)+r*(q+q%2)
def g(n):
 if not n:return 1
 q,r=divmod(n,2)
 return g(q)+(1-r)*(q+q%2)
def bitpairs(n):
 bits=[i for i in range(n.bit_length()) if n>>i&1]
 return sum(2**(j-i-1)+(j==i+1) for i in bits for j in bits if j>i)
counts={'carrier_identity':0,'generic_envelopes':0,'alternating_formula':0,'pump_numeric_identity':0}
for n in range(16384):
 assert f(n)==bitpairs(n)
 assert f(n)+g(n)==n+1-n%2
 counts['carrier_identity']+=1
for m in range(1,12):
 can=[0]+[2**m-2**j for j in range(m)]
 for K in (1,2,3,7,29):
  for h in (1,2,3,7,29):
   for S in (0,1,2,9,100):
    score=lambda p:K*p-h*f(p)-S*(p%2)
    M=max(map(score,can))
    assert max(map(score,range(2**m)))==M
    assert all(p in can for p in range(2**m) if score(p)==M)
    counts['generic_envelopes']+=1
for m in range(1,100):
 v=2*(4**m-1)//3
 for a,b,c in ((1,1,1),(3,2,5),(99,29,263)):
  F=a*g(v)-b+c*(v+1)
  assert 9*F==(4*a+6*c)*4**m+6*a*m+5*a+3*c-9*b
  assert 9*F-(6*a+9*c)*v==6*a*m+9*(a+c-b)
  counts['alternating_formula']+=1
# Every pump U V^j W has numerical value C*2^(d*j)+D, including zero V.
for U in ('','0','1','101'):
 for V in ('0','1','01','10','1010'):
  for W in ('','0','1','101'):
   val=lambda s:int(s,2) if s else 0
   d=len(V);e=len(W)
   C=(val(U)+Fraction(val(V),2**d-1))*2**e
   D=val(W)-Fraction(val(V)*2**e,2**d-1)
   for j in range(8):
    assert val(U+V*j+W)==C*2**(d*j)+D
    counts['pump_numeric_identity']+=1
hostile_scores={p:p-10*f(p)+100*(p%2) for p in range(8)}
assert max(hostile_scores,key=hostile_scores.get)==1
print(json.dumps({'status':'PASS','checks':counts,'negative_S_counterexample':{'m':3,'K':1,'h':10,'S':-100,'unique_maximizer':1},'scope':'Exact formula checks and independent all-length arguments in AUDIT.md; no all-target coverage conclusion.'},indent=2))
