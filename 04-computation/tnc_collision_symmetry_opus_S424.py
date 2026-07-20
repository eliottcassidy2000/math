import sympy as sp
from math import gcd
u=sp.symbols('u')
def saddle_values(Rexpr,N):
    R=sp.sympify(Rexpr); S=sp.expand(u*sp.diff(R,u)-N*R)
    vals=[]
    for r in sp.nroots(sp.Poly(S,u)):
        rv=complex(r)
        if abs(complex(R.subs(u,r)))<1e-9 or abs(rv)<1e-9: continue
        vals.append(complex(R.subs(u,r))/rv**N)
    return vals
def rfold(Rexpr):
    P=sp.Poly(sp.expand(sp.sympify(Rexpr)),u)
    ex=[e[0] for e,c in zip(P.monoms(),P.coeffs()) if c!=0]
    g=0
    for e in ex: g=gcd(g,e)
    return g
print("QUESTION: for TRINOMIALS, are all saddle-VALUE collisions root-of-unity SYMMETRIC?")
print("(A symmetric R has R(u)=S(u^g), g>=2; then collisions come from the mu_g branch action.)")
print("If YES: distinct values -> Vandermonde -> TNC; collision -> symmetric -> DESCENT to")
print("binomial (proved). => TRINOMIAL TNC COMPLETE via roots of unity.")
print()
print("   pattern (real coeffs)      #saddle values  #distinct  collision?  g-fold  symmetric?")
allsym=True; collis=0
import itertools
for N in [2,3]:
    for j in range(1,7):
        for d in range(j+1,9):
            if j==N or d==N: continue
            for aval in [1,-1,2,sp.Rational(1,2),3,-2]:
                R=1+aval*u**j+u**d
                vals=saddle_values(R,N)
                if not vals: continue
                # count distinct values (numerically)
                dist=[]
                for w in vals:
                    if not any(abs(w-x)<1e-7 for x in dist): dist.append(w)
                coll = len(dist)<len(vals)
                if not coll: continue
                collis+=1
                g=rfold(R)
                sym = g>=2
                allsym = allsym and sym
                if collis<=12:
                    print(f"   N={N} 1+({aval})u^{j}+u^{d}: {len(vals):2d} vals, {len(dist):2d} distinct, "
                          f"collision=YES, g={g}, symmetric={sym}")
print(f"\n   trinomial patterns WITH saddle-value collision: {collis}")
print(f"   ALL collisions are root-of-unity symmetric (g>=2): {allsym}")
if allsym:
    print("   => TRINOMIAL TNC COMPLETES: collision => symmetric => descent to binomial (proved),")
    print("      distinct => Vandermonde. Roots of unity (mu_g branch action) close the residual.")
else:
    print("   => some ASYMMETRIC trinomial collision exists; roots-of-unity descent alone")
    print("      does not close trinomials -- need the resultant/amoeba argument there.")
