"""Exception-safe exact check of the canonical triangle half, no assertions."""
import sympy as S
import json

a,b,h=S.symbols('a b h',nonzero=True)
N=a*b/2
X=N*(a+h)/b
Y=2*N*N*(a+h)/b**2
slope=(3*X**2-N*N)/(2*Y)
xx=slope*slope-2*X
yy=slope*(X-xx)-Y
results={}
for label,expr in [('direct_half_x',xx-h*h/4),('direct_half_y',yy-h*(a*a-b*b)/8)]:
 numerator=S.cancel(expr).as_numer_denom()[0]
 remainder=S.rem(S.Poly(numerator,h),S.Poly(h*h-a*a-b*b,h))
 if not remainder.is_zero:
  raise RuntimeError(label)
 results[label]='PASS modulo h²=a²+b²'
print(json.dumps(results,indent=2))
