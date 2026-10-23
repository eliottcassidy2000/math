import sympy as s
from collections import defaultdict
x=s.symbols('x'); n=57
f=sum(x**a-x**(57-a) for a in [1,7,49])
assert s.gcd(f,x**57-1)==x-1
assert s.rem(f.subs(x,x**7)-f,x**57-1)==0
print('EXACT Frobenius invariant alternating polynomial gcd:',s.gcd(f,x**57-1))
P=s.zeros(n)
for j in range(n):P[(j+1)%n,j]=1
F=s.zeros(n)
for j in range(n):F[(7*j)%n,j]=1
B=s.zeros(57,56)
for j in range(56):B[j,j]=1;B[56,j]=-1
O=sum((P**a-P**(57-a) for a in [1,7,49]),s.zeros(n))
assert O.T==-O and P.T*O*P==O and F.T*O*F==O
assert (B.T*O*B).det()!=0
T=P**19
assert s.trace(T)==0
# fixed W dimension = ambient fixed dimension minus one.
assert n-(T-s.eye(n)).rank()-1==18
print('EXACT W dimension56, restricted form nondegenerate, C57 and Frobenius invariance, order3 fixed18 trace-1 PASS')
d=defaultdict(list)
for a in range(9):
 for b in range(9-a):
  c=8-a-b
  if (b+2*c)%3==0:d[a*(a-1)+2*b*c].append((a,b,c))
assert set(d)=={2,14,20,32,56}
traces={(3*k-56)//2 for k in d}
assert traces=={-25,-7,2,20,56}
assert 18 not in d and -1 not in traces
print('EXACT SL8 order3 fixed dimensions:',sorted(d))
print('EXACT E7 minuscule order3 possible traces including identity:',sorted(traces))
print('ALL SINGER ADDENDUM TESTS PASS')
