from fractions import Fraction as F
from math import factorial

s_num,s_den=116,125
z=F(181,500)
assert z**s_den * 3**s_num>1
A=z+z**3+z**4/(1-2*z)
assert A+4*z*z<1
assert 1/(1-z)<F(8,5)
x=F(8,5)
# Sum n <= 12, and dominate all higher term ratios by x/14.
e=sum(x**n/factorial(n) for n in range(13))+(x**13/factorial(13))/(1-x/14)
assert e<5

states=('', '0','1','00','01','10','11')
b=(10000,10000,9015,10000,9015,8869,8274)
t=F(6865,10000)
lam=F(1902,1000)
assert t**s_den>F(2,3)**s_num
assert lam**s_den<2**s_num
for i,a in enumerate(states):
    m=sum((t*b[0] if len(a)==2 and a[0]==x=='1' else b[states.index((a+x)[-2:])])
          for x in '01')
    assert m<=lam*b[i]
print(f'grammar bound {float(A+4*z*z):.9f}<1; scanner bound {float(lam)}<2^({s_num}/{s_den})')
