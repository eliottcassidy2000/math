# Test klein's domination: E[P^m] = sum_a gamma_a a!  (gamma_a=[z^a w^a]P^m). Does the TOP term
# dominate, or does mass shift to a saddle -- especially with NON-CONSTANT leading coefficients?
from fractions import Fraction as Fr
from math import factorial
def L(poly):
    s=Fr(0)
    for (p,q),c in poly.items():
        if p==q: s+=c*factorial(p)
    return s
def gammas(poly):  # {a: [z^a w^a] coeff}
    d={}
    for (p,q),c in poly.items():
        if p==q: d[p]=d.get(p,Fr(0))+c
    return d
def pmul(p,q):
    r={}
    for k1,u in p.items():
        for k2,v in q.items():
            k=(k1[0]+k2[0],k1[1]+k2[1]); r[k]=r.get(k,Fr(0))+u*v
    return {k:v for k,v in r.items() if v!=0}
def padd(*ps):
    r={}
    for p in ps:
        for k,v in p.items(): r[k]=r.get(k,Fr(0))+v
    return {k:v for k,v in r.items() if v!=0}
def scal(p,c): return {k:v*Fr(c) for k,v in p.items() if v!=0}
Z={(1,0):Fr(1)}; W={(0,1):Fr(1)}; one={(0,0):Fr(1)}; U={(1,1):Fr(1)}
def Um(k):
    r=one
    for _ in range(k): r=pmul(r,U)
    return r

# CASE 1: two-sided LEADING symbol (top degree has both +1 and -1): P = z + w + b(u), b const.
# CASE 2: leading symbol ONE-sided, lower straddles, NON-CONSTANT coeff: charges {-2,-1,0,1,2}
#   P = z^2 (top, one-sided pos) + w*(1+u) (non-constant, charge -1) + ...  build a genuine two-sided P.
def report(name, P, mmax=8):
    print(f"\n=== {name} ===")
    for m in range(2,mmax+1):
        Pm=one
        for _ in range(m): Pm=pmul(Pm,P)
        g=gammas(Pm)
        if not g: print(f"  m={m}: no charge-0 (E=0)"); continue
        terms={a: g[a]*factorial(a) for a in g}
        E=sum(terms.values())
        amax=max(g); 
        # dominant a by |term|
        adom=max(terms, key=lambda a: abs(terms[a]))
        top=abs(terms[amax]); tot=sum(abs(v) for v in terms.values())
        print(f"  m={m}: E[P^m]={'0' if E==0 else '!=0'}({str(E)[:12]}), a_max={amax}, DOMINANT a={adom}"
              f" {'(=top)' if adom==amax else '(SADDLE below top!)'}, |top|/|total|={float(top/tot):.3f}")

# two-sided leading symbol (should: top dominates, E!=0)
report("CASE1: z + w (two-sided top)", padd(Z,W))
# one-sided top (z^2), two-sided lower, NON-CONSTANT coeff on charge -1: P = z^2 + (1+u)w + w^2? 
# charges: z^2=+2, (1+u)w = -1 (deg 1 and 3), w^2=-2.  top degree 2 = {z^2, w^2} = two-sided actually.
# Make top ONE-sided: P = z^2 + z + (1+u) w  (charges +2,+1,-1; top deg 2 = z^2 only, one-sided pos)
report("CASE2: z^2 + z + (1+u)w  [one-sided top z^2, non-const charge-1]",
       padd(pmul(Z,Z), Z, pmul(padd(one,U),W)))
# genuinely two-sided, one-sided top, NON-CONSTANT leading: P = z^2*(1) + (u)*w  won't be two-sided enough
# P = z^2 + w*(1+u) : charges +2 (top, one-sided) and -1 (deg1,3). two-sided. non-const charge-1 coeff.
report("CASE3: z^2 + (1+u)w  [one-sided top, non-const charge-1, genuinely two-sided]",
       padd(pmul(Z,Z), pmul(padd(one,U),W)))
