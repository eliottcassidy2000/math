"""Numerical consistency only; the report contains the all-order proof."""
import mpmath as m,json
m.mp.dps=50;b=m.mpf(10);a=m.log(b)
def w(x):
 if x==0:return m.mpf(0)
 if b*x==1:return 1/(b-1)
 return (b*x-1)/((b-1)*m.log(b*x))
rows=[]
for n in range(9):
 direct=m.quad(lambda u:1/((n+1)*a+m.log(u)),[1,b])/(b-1)
 moment=m.quad(lambda x:x**n*w(x),[0,1/b,1])
 assert abs(direct-moment)<m.mpf('1e-45')
 rows.append(dict(n=n,direct=str(direct),moment=str(moment)))
q=lambda x:x*(x-1)**2*(2*x-1)
positive=m.quad(lambda x:q(x)**2*w(x),[0,1/b,1])
assert positive>0
print(json.dumps({'precision':50,'moments':rows,'li_square_norm':str(positive),'exact_prime_same_square':'-1152073/2250000000'},indent=2))
