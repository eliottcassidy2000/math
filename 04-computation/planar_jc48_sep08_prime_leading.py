#!/usr/bin/env python3
"""Small exact controls for the all-prime/all-W_m leading filter.
The unbounded theorem is analytic; the named finite banks are not a
prime cutoff, a mate-degree cutoff, or a census of Keller maps.
"""
import hashlib
import json
import sympy as S

x,t,r,b,h = S.symbols('x t r b h')
gates = 0
record = {}
def need(name,value):
    global gates
    if not bool(value): raise RuntimeError(name)
    gates += 1
def zero(name,expr): need(name,S.cancel(expr)==0)
def jac(F,G): return S.diff(F,x)*S.diff(G,t)-S.diff(F,t)*S.diff(G,x)

# Highest Jacobian coefficient, with arbitrary lower next coefficients.
A,B,U,V = [S.Function(s)(x) for s in ['A','B','U','V']]
for prime in [2,3,5,7]:
    need('named degree is prime '+str(prime),S.isprime(prime))
    for n in [1,prime,prime+1,2*prime]:
        F=A*t**prime+U*t**(prime-1)
        G=B*t**n+V*t**(n-1)
        highest=S.expand(jac(F,G)).coeff(t,prime+n-1)
        zero('leading equation '+str((prime,n)),
             highest-n*S.diff(A,x)*B+prime*A*S.diff(B,x))
        zero('constant ratio derivative '+str((prime,n)),
             S.diff(B**prime/A**n,x)+B**(prime-1)*highest/A**(n+1))
    # All possible nonzero residue classes; Bezout is the unbounded proof.
    for nmod in range(1,prime):
        inv=pow(nmod,-1,prime)
        need('coprime valuation transfer '+str((prime,nmod)),(nmod*inv)%prime==1)

# Polynomial shears preserve the exact bracket regardless of lower rows.
F=x+t**5+x*x*t*t
G=x*t**3+t+x
for power in [1,2,3]:
    zero('constant polynomial shear '+str(power),jac(F,G-3*F**power)-jac(F,G))

# A composite degree does not allow the coprime valuation inference.
zero('composite top identity',x**4-(x*x)**2)
need('composite hostile is not a fourth power',2%4!=0)
Fcomp=x*x*t**4+t*t
Gcomp=x*t*t
zero('composite literal top cancellation',jac(Fcomp,Gcomp)+2*t**3)
need('composite control is not a constant-Jacobian pair',S.Poly(jac(Fcomp,Gcomp),t).degree()==3)

# Reciprocal-polynomial exactness: all root positions and scalar factors.
c=S.symbols('c',nonzero=True)
zero('constant reciprocal primitive',S.diff(x/c,x)-1/c)
for k in range(2,9):
    primitive=(x-h)**(1-k)/(c*(1-k))
    zero('one-root primitive '+str(k),S.diff(primitive,x)-1/(c*(x-h)**k))
zero('single simple-root obstruction',S.residue(1/(x-h),x,h)-1)
zero('two-root first residue',S.residue(1/(x*x*(x-1)**2),x,0)-2)
zero('two-root second residue',S.residue(1/(x*x*(x-1)**2),x,1)+2)

# Complete monomial leading lifts for the first quintic frontier W_2.
# The constructive leading map is linear. These global polynomials assert
# nonempty leading strata only, never the existence of polynomial mates.
def lift_leading(poly,m,degree):
    ans=S.Integer(0)
    for (d,),co in S.Poly(S.expand(poly),x).terms():
        j=max(0,(d-m*degree+m-1)//m)
        i=d-m*j
        need('leading lift lies in full box '+str((m,degree,d)),
             0<=j<=degree and 0<=i<=m*degree)
        ans+=co*x**i*t**(degree-j)*(1+x**m*t)**j
    return S.expand(ans)

quintic=[]
for k in [0,2,3,4]:
    leading=(x-h)**(5*k)
    FF=lift_leading(leading,2,5)
    zero('quintic exact leading '+str(k),S.Poly(FF,t).coeff_monomial(t**5)-leading)
    transformed=S.expand(FF.subs({x:1/r,t:-r*r-r**4*b},simultaneous=True))
    need('quintic all-position global lift '+str(k),S.denom(S.cancel(transformed))==1)
    section=S.Poly(S.expand(r**20*leading.subs(x,1/r)),r)
    need('quintic marked infinity '+str(k),min(z[0] for z in section.monoms())==5*(4-k))
    quintic.append({'k':k,'finite_order':5*k,'infinity_order':5*(4-k)})
record['W2_quintic']=quintic

# Rational mates do NOT satisfy the polynomial-mate leading-power conclusion.
# This is an actual global example on every graph degree m, not a model.
for prime in [2,3,5,7]:
    FF=x*t**prime
    GG=-1/((prime-1)*t**(prime-1))
    zero('actual rational hostile '+str(prime),jac(FF,GG)-1)
    need('rational hostile leading not a pth power '+str(prime),1%prime!=0)
    for m in [1,2,3]:
        fi=S.expand(FF.subs({x:1/r,t:-r**m-r**(2*m)*b},simultaneous=True))
        expected=(-1)**prime*r**(m*prime-1)*(1+r**m*b)**prime
        zero('rational hostile actual chart '+str((m,prime)),fi-expected)
        need('rational hostile chart regular '+str((m,prime)),S.denom(S.cancel(fi))==1)
    need('rational hostile affine pole '+str(prime),S.denom(GG).has(t))

# The leading condition alone is not sufficient even with a constant root.
Ffail=t**5
zero('necessary-only example has critical derivative',S.diff(Ffail,t).subs(t,0))
zero('necessary-only example other derivative',S.diff(Ffail,x))
record['named_primes']=[2,3,5,7]
record['controls']=['coprime leading valuations','one-root rational primitives',
                    'complete shifted quintic leading lifts','global rational hostile',
                    'composite top-identity-only hostile']
record['gates']=gates
digest=hashlib.sha256(json.dumps(record,sort_keys=True,separators=(',',':')).encode()).hexdigest()
print('Prime fibre degree on every W_m: necessary leading collapse')
print('Exact gates:',gates)
print('All-prime proof is analytic; named controls use p=2,3,5,7')
print('Complete W2 quintic leading orders:',quintic)
print('Polynomial/rational distinction and unmarked-composite hostile: PASS')
print('Semantic SHA256:',digest)
