from check_automaton import *
from fractions import Fraction
A=[[False,True],[False,False]]
# A sequence u M^m f of dimension D is uniquely determined by D consecutive
# zero residuals for any linear recurrence, by Cayley-Hamilton.
def certify(prefix,block,coeffs,start=0):
 v,acc,tr=compile_automaton(2,A);D=len(v)
 for c in prefix:v=step(v,tr[c])
 h=[]
 for m in range(start+D+len(coeffs)):
  h.append(sum(x for x,b in zip(v,acc) if b))
  for c in block:v=step(v,tr[c])
 for m in range(start,start+D):
  assert sum(c*h[m+j] for j,c in enumerate(coeffs))==0,(prefix,block,m)
 return h
h=certify([1,2,0],[2],[2,-3,1],1)
assert h[1:3]==[29,59]
print('PROVED via dimension36 recurrence: H(1,2,0,2^m)=15*2^m-1 for m>=1')
h=certify([3,1,2],[0],[-6,11,-6,1],1)
assert h[1:4]==[43,131,379]
print('PROVED via dimension36 recurrence: H(3,1,2,0^m)=12*3^m+8*2^m-9 for m>=1')
for p,b in [([1,2],[0,2]),([1,2,0,2],[1,0]),([1],[2,0])]:
 h=certify(p,b,[-144,376,-337,123,-19,1])
 # Characteristic (x-9)(x-4)^2(x-1)^2.
 print('PROVED periodic-block recurrence',p,b,'initial',h[:5])
