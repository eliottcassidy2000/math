"""Exact identities for the recovered fourfold DG unit-response bridge.

This does not import a prior producer and has no repository side effects.
The all-degree nonvanishing argument is in incoming_synthesis.md.
"""
import sys
import sympy as s

sys.stdout.reconfigure(newline='\n')

x, t, a, h, b0, z, w = s.symbols('x t a h b0 z w', nonzero=True)
gates = 0

def check(value):
    global gates
    if s.cancel(value) != 0:
        raise RuntimeError(str(value))
    gates += 1

F = a*(x-h)**4*t + a*(x*x-4*h*x) + b0
f0 = b0-3*a*h*h
g = F-f0
G0 = 1/(3*a*(x-h)**3)
K = z**3*t+z-2*h

def jac(A,B):
    return s.diff(A,x)*s.diff(B,t)-s.diff(A,t)*s.diff(B,x)

check(g.subs(x,z+h)-a*z*K)
check(jac(F,G0)-1)
P3 = a*a*((x-h)**3*t+x-3*h)**3/3
check(g**3*G0-P3)
check(jac(F,P3)-g**3)
check(s.diff(F,x).subs(x,h)+2*a*h)
check(K.subs(z,0)+2*h)

# Entire scalar principal part in the actual parameter g.
c3 = -8*a*a*h**3/3
c2 = -2*a*h
R = s.cancel((G0-c3/g**3-c2/g**2).subs(x,z+h))
num, den = s.fraction(R)
if s.cancel(den.subs(z,0)) == 0:
    raise RuntimeError('Remainder is not regular at E1')
gates += 1
check(s.cancel((K**3+6*h*z*K+8*h**3)/z**3)
      -(1+3*t*(z-2*h)**2+6*h*t*z+3*t*t*z**3*(z-2*h)+t**3*z**6))

# E2 is the other reduced component; G0 there is w^3/(3a).
check(((x-h)**3*t+x-3*h).subs({x:h+1/w,t:2*h*w**3-w*w}))
check(G0.subs(x,h+1/w)-w**3/(3*a))
print('Fourfold DG unit response: exact primary order 3 by the accompanying component proof.')
print('Scalar principal part: -8*a^2*h^3/(3*g^3) - 2*a*h/g^2; no g^-1 term.')
print(f'Always-active exact gates: {gates}')
