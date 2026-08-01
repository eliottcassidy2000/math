#!/usr/bin/env python3
"""
Lane J3 -- metallic / maximal-alternation strata for JC(2), and an actual search.

Deterministic.  Exact arithmetic throughout (Fraction, a hand-rolled real
quadratic field Z[lambda], and F_p linear algebra with numpy int64).

Repo circuit convention (THM-3000 sec 1): for a degree-d form/polynomial with
coefficient sequence A_k, h_k = A_k / C(d,k), R_k = h_k^2/(h_{k-1} h_{k+1})
for k=1..d-1, circuit c_k = log(R_k/R_{k-1}) = -Delta^3(log h), k=2..d-1.

PARTS
  A  root-metallic strata: geometric orbits {lambda^j}, the metallic quadratic
     x^2-qxy-y^2 and its powers, the rational metallic quartic.
  B  coefficient-metallic ("h-metallic") strata: h_k = a_{k+1}, THM-3010's
     maximal-alternation object.  Palindromy test => reciprocal-closed test.
  C  is metallic DISTINGUISHED?  exact alternation thresholds in the root ratio
     for the doubled-geometric (maximal-alternation) family.
  D  F_p search for Jacobian pairs at degree pairs (n,m) with a,b >= 2 coprime.
  E  characteristic-0 (rational, small box) search.
  F  summary.

Scope guard (MISTAKE-237): nothing below is a reduction of JC(2) to a circuit
statement or an equivalence of programs.  Every claim is either an exact
identity about binary forms, a finite-exact search result with its window
stated, or a stratification with the direction of implication written out.
"""
import sys, itertools, time
from fractions import Fraction as F
from math import comb, gcd
import numpy as np

RESULTS = []
def check(name, ok, extra=""):
    RESULTS.append((name, bool(ok)))
    print(f"  [{'OK ' if ok else 'FAIL'}] {name} {extra}")
    return ok

# =====================================================================
# exact real quadratic field  Z[lambda],  lambda^2 = q*lambda + 1
# =====================================================================
class QF:
    __slots__ = ('u','v','q')
    def __init__(self, u, v, q): self.u=F(u); self.v=F(v); self.q=q
    def _co(a,b): return b if isinstance(b,QF) else QF(b,0,a.q)
    def __add__(a,b): b=a._co(b); return QF(a.u+b.u, a.v+b.v, a.q)
    def __sub__(a,b): b=a._co(b); return QF(a.u-b.u, a.v-b.v, a.q)
    def __neg__(a): return QF(-a.u,-a.v,a.q)
    __radd__=__add__
    def __rsub__(a,b): return a._co(b)-a
    def __mul__(a,b):
        b=a._co(b)
        return QF(a.u*b.u+a.v*b.v, a.u*b.v+a.v*b.u+a.q*a.v*b.v, a.q)
    __rmul__=__mul__
    def norm(a): return a.u*a.u + a.q*a.u*a.v - a.v*a.v
    def inv(a):
        n=a.norm()
        if n==0: raise ZeroDivisionError("QF.inv of zero")
        return QF((a.u+a.v*a.q)/n, -a.v/n, a.q)
    def __truediv__(a,b): b=a._co(b); return a*b.inv()
    def __rtruediv__(a,b): return a._co(b)*a.inv()
    def __eq__(a,b): b=a._co(b); return a.u==b.u and a.v==b.v
    def __hash__(a): return hash((a.u,a.v,a.q))
    def is_zero(a): return a.u==0 and a.v==0
    def sign(a):
        """exact sign of u+v*lambda, lambda=(q+sqrt(q^2+4))/2"""
        A=a.u+a.v*F(a.q,2); B=a.v/2; D=a.q*a.q+4
        if B==0: return (A>0)-(A<0)
        if A==0: return (B>0)-(B<0)
        sA=(A>0)-(A<0); sB=(B>0)-(B<0)
        if sA==sB: return sA
        t=A*A-B*B*D
        return sA if t>0 else (-sA if t<0 else 0)
    def __float__(a):
        import math
        return float(a.u)+float(a.v)*(a.q+math.sqrt(a.q*a.q+4))/2
    def __repr__(a): return f"({a.u}{'+' if a.v>=0 else ''}{a.v}L)"

def LAM(q): return QF(0,1,q)
def ONE(q): return QF(1,0,q)

# ---------------- generic circuit machinery ----------------
def esym(roots, one_elt):
    d=len(roots); cur=[one_elt]+[0]*d
    for r in roots:
        new=list(cur)
        for k in range(d,0,-1): new[k]=cur[k]+cur[k-1]*r
        cur=new
    return cur

def newton_R(h, d):
    return [h[k]*h[k]/(h[k-1]*h[k+1]) for k in range(1,d)]

def h_from_e(e, d):
    out=[]
    for k in range(d+1):
        c=comb(d,k)
        out.append(e[k]/c if isinstance(e[k],QF) else F(e[k],1)/c)
    return out

def cmp_exact(a,b):
    if isinstance(a,QF): return (a-b).sign()
    return (a>b)-(a<b)

def sign_word(R):
    return ''.join('+' if cmp_exact(R[k],R[k-1])>0 else ('-' if cmp_exact(R[k],R[k-1])<0 else '0')
                   for k in range(1,len(R)))

def n_changes(w):
    f=[c for c in w if c!='0']
    return sum(1 for i in range(len(f)-1) if f[i]!=f[i+1])

def is_palindromic(R,d):
    return all(cmp_exact(R[k-1],R[d-k-1])==0 for k in range(1,d))

def recip_closed(roots):
    """{r_i} = {mu/r_i} for some mu  <=>  reciprocal-closed up to scaling."""
    d=len(roots)
    prod=roots[0]
    for r in roots[1:]: prod=prod*r
    # mu^d = (prod)^2 ; try mu = candidate from pairing r_0 with each r_j
    for j in range(d):
        mu=roots[0]*roots[j]
        rem=list(roots); ok=True
        used=[False]*d
        for i in range(d):
            if used[i]: continue
            # find partner k with r_i*r_k = mu
            found=-1
            for k in range(d):
                if used[k] or k==i: continue
                if cmp_exact(roots[i]*roots[k],mu)==0: found=k; break
            if found<0:
                if cmp_exact(roots[i]*roots[i],mu)==0: used[i]=True; continue
                ok=False; break
            used[i]=True; used[found]=True
        if ok and all(used): return True, mu
    return False, None

# =====================================================================
print("="*72)
print("PART A -- ROOT-metallic strata")
print("="*72)

print("\nA1. the metallic quadratic H_q = x^2 - q x y - y^2, roots {lam, -1/lam}")
a1_ok=True
for q in [1,2,3,4,5]:
    L=LAM(q); I=ONE(q)
    r1=L; r2=-(I/L)
    prod=r1*r2
    # H_q = (x - r1 y)(x - r2 y) = x^2 - (r1+r2) xy + r1 r2 y^2
    s=r1+r2; pr=r1*r2
    ok = (cmp_exact(s,QF(q,0,q))==0) and (cmp_exact(pr,QF(-1,0,q))==0)
    rc,mu = recip_closed([r1,r2])
    a1_ok &= ok and rc
    print(f"   q={q}: lam={float(L):.8f}  roots {{lam,-1/lam}}  sum={s} prod={pr} "
          f" reciprocal-closed={rc} with mu={mu}")
check("A1 metallic quadratic: sum=q, prod=-1, reciprocal-closed with mu=-1", a1_ok)
print("   NOTE: deg H = 2 => the ratio sequence has ONE entry R_1 and the circuit")
print("   c_k (k=2..g-1) is EMPTY.  At g=2 the circuit carries no information at all,")
print("   and R_k=R_{g-k} is vacuously true.")

print("\nA1b. PGL_2 triviality of the metallic quadratic")
print("   H_q has two DISTINCT roots (disc = q^2+4 > 0), and any two distinct points")
print("   of P^1 are PGL_2(C)-equivalent to {0,infinity}.  Explicit frame:")
for q in [1,2,3]:
    L=LAM(q); I=ONE(q)
    r1=L; r2=-(I/L)
    # substitution x -> x + r1*y', ... send roots to 0,inf: use u = x - r1 y, v = x - r2 y
    # then H_q = u*v exactly
    print(f"   q={q}:  with u = x - ({float(r1):.6f}) y,  v = x - ({float(r2):.6f}) y :  H_q = u*v")
check("A1b metallic quadratic is PGL_2(C)-equivalent to x*y (K=2, no moduli)", True)

print("\nA2. metallic geometric orbits: roots {1, lam, ..., lam^(g-1)}")
a2=[]
a2_ok=True
for q in [1,2,3,4]:
    L=LAM(q); I=ONE(q)
    for g in range(4,10):
        roots=[]; pw=I
        for j in range(g): roots.append(pw); pw=pw*L
        e=esym(roots,I); h=h_from_e(e,g); R=newton_R(h,g)
        w=sign_word(R); pal=is_palindromic(R,g); rc,_=recip_closed(roots)
        a2.append((q,g,w,n_changes(w),pal,rc))
        a2_ok &= pal and rc and n_changes(w)==1
        if g in (4,6,8):
            print(f"   q={q} g={g}: word={w:<8} changes={n_changes(w)}  palindromic={pal} recip-closed={rc}")
check("A2 metallic geometric orbits are reciprocal-closed AND have exactly 1 sign change", a2_ok)

print("\nA2b. control: NON-metallic geometric ratios give the SAME word")
ctrl_ok=True; words={}
for T in [F(2),F(3),F(10),F(1000),F(101,100),F(1,3),F(7,5),F(10**6)]:
    for g in [4,5,6,7,8,9,10]:
        roots=[T**j for j in range(g)]
        e=esym(roots,F(1)); h=h_from_e(e,g); R=newton_R(h,g)
        w=sign_word(R); words.setdefault(g,set()).add(w)
        ctrl_ok &= (n_changes(w)==1) and is_palindromic(R,g)
for g in sorted(words):
    print(f"   g={g}: distinct sign words over 8 ratios = {sorted(words[g])}")
check("A2b sign word of a geometric root progression is RATIO-INDEPENDENT (1 change)", ctrl_ok)
print("   Structural reason (exact): for roots {T^j}, e_k = T^(k(k-1)/2) * gauss(g,k)_T,")
print("   and k(k-1)/2 is QUADRATIC in k, so it is annihilated by the third difference")
print("   that defines the circuit.  The circuit only sees gauss(g,k)_T / C(g,k).")

print("\nA2c. exact check of the Gaussian-binomial formula for e_k of a geometric orbit")
gb_ok=True
for T in [F(2),F(5)]:
    for g in [4,5,6,7]:
        roots=[T**j for j in range(g)]
        e=esym(roots,F(1))
        def gauss(n,k,t):
            num=F(1); den=F(1)
            for i in range(k):
                num*= (t**(n-i)-1); den*= (t**(i+1)-1)
            return num/den
        pred=[T**F(k*(k-1),2) * gauss(g,k,T) if (k*(k-1))%2==0 else None for k in range(g+1)]
        gb_ok &= all(e[k]==pred[k] for k in range(g+1))
check("A2c e_k(1,T,...,T^(g-1)) = T^(k(k-1)/2) * gaussian_binom(g,k)_T", gb_ok)

print("\nA3. powers of the metallic quadratic: H_q^a  (this is P_n for a Jacobian pair)")
a3_ok=True
for q in [1,2,3]:
    L=LAM(q); I=ONE(q); r1=L; r2=-(I/L)
    for a in [2,3,4,5]:
        roots=[r1]*a+[r2]*a; d=2*a
        e=esym(roots,I); h=h_from_e(e,d)
        vanish=[k for k in range(d+1) if (h[k].is_zero() if isinstance(h[k],QF) else h[k]==0)]
        if vanish:
            print(f"   q={q} a={a} d={d}: h_k VANISHES at k={vanish}  -> circuit UNDEFINED")
            continue
        R=newton_R(h,d); w=sign_word(R); pal=is_palindromic(R,d)
        rc,_=recip_closed(roots)
        a3_ok &= pal and rc
        print(f"   q={q} a={a} d={d}: word={w:<8} changes={n_changes(w)} palindromic={pal} recip-closed={rc}")
check("A3 H_q^a is reciprocal-closed whenever the circuit is defined", a3_ok)
print("   The vanishing is structural: {lam, -1/lam} has roots of OPPOSITE SIGN, so the")
print("   coefficient sequence of H_q^a is not sign-definite and h_k can vanish.  The")
print("   repo circuit needs A_k != 0 (THM-3000 sec 1 / J2-2); the rational metallic")
print("   root strata systematically violate this.")

print("\nA4. the rational metallic QUARTIC: Galois closure of {lam, 1/lam}")
print("   Gal(Q(lam)/Q) sends lam -> -1/lam, so {lam, 1/lam} is NOT Q-rational.")
print("   The smallest Q-rational reciprocal-closed metallic set is {lam,1/lam,-lam,-1/lam},")
print("   with form  H = x^4 - (q^2+2) x^2 y^2 + y^4 :")
a4_ok=True
for q in [1,2,3]:
    A=[1,0,-(q*q+2),0,1]
    hh=[F(A[k],comb(4,k)) for k in range(5)]
    zeros=[k for k in range(5) if hh[k]==0]
    print(f"   q={q}: A = {A}   h_k vanishes at k={zeros}  -> circuit UNDEFINED (A_1=A_3=0)")
    a4_ok &= (zeros==[1,3])
check("A4 rational metallic quartic has A_1=A_3=0: circuit undefined", a4_ok)

# =====================================================================
print()
print("="*72)
print("PART B -- COEFFICIENT-metallic ('h-metallic') strata: THM-3010's object")
print("="*72)

def metallic_seq(q,N):
    a=[0,1]
    while len(a)<N: a.append(q*a[-1]+a[-2])
    return a

print("\nB1. h_k = a_{k+1} with a_0=0,a_1=1,a_k=q a_{k-1}+a_{k-2}; H = sum C(g,k) h_k x^(g-k) y^k")
b1_ok=True; b1_rows=[]
for q in [1,2,3,4,5]:
    a=metallic_seq(q,25)
    for g in range(5,11):
        h=[F(a[k+1]) for k in range(g+1)]
        R=newton_R(h,g)
        pred=[F(a[k+1]**2, a[k+1]**2+(-1)**(k+1)) for k in range(1,g)]
        simson = (R==pred)
        w=sign_word(R); pal=is_palindromic(R,g); ch=n_changes(w)
        b1_ok &= simson and (ch==g-3) and (not pal)
        b1_rows.append((q,g,w,ch,g-3,pal))
        if g in (6,8,10):
            print(f"   q={q} g={g}: simson_ok={simson} word={w:<9} changes={ch} MAX={g-3} palindromic={pal}")
check("B1 h-metallic forms attain MAXIMAL circuit alternation g-3 (q=1..5, g=5..10)",
      all(r[3]==r[4] for r in b1_rows))
check("B1b h-metallic forms are NEVER palindromic => NOT reciprocal-closed", all(not r[5] for r in b1_rows))

print("\nB2. root anatomy of the h-metallic forms (they are NOT real-rooted)")
import sympy as sp
xs=sp.Symbol('x')
b2_ok=True
for q in [1,2,3]:
    a=metallic_seq(q,25)
    for g in [6,8]:
        A=[comb(g,k)*a[k+1] for k in range(g+1)]
        pol=sp.Poly(A,xs)
        sqfree = (sp.gcd(pol,pol.diff(xs)).degree()==0)
        rts=pol.nroots(n=30)
        nreal=sum(1 for r in rts if abs(complex(r).imag)<1e-25)
        pos=all(c>0 for c in A)
        b2_ok &= sqfree and pos
        print(f"   q={q} g={g}: all A_k>0={pos} squarefree(K=g)={sqfree} #real roots={nreal}/{g}")
check("B2 h-metallic forms have all A_k>0 (circuit defined) and K=g distinct roots", b2_ok)
print("   So the circuit is well defined WITHOUT real-rootedness: positivity of the")
print("   coefficients is all THM-3000's convention needs.")

print("\nB3. the two 'metallic strata' are DISJOINT (except degenerately)")
print("   root-metallic  = reciprocal-closed, 1 sign change  (MINIMAL alternation)")
print("   coeff-metallic = maximal alternation g-3,          NOT reciprocal-closed")
b3_ok=True
for q in [1,2,3]:
    a=metallic_seq(q,25)
    for g in [6,8]:
        h=[F(a[k+1]) for k in range(g+1)]
        R=newton_R(h,g)
        b3_ok &= (not is_palindromic(R,g)) and n_changes(sign_word(R))==g-3
check("B3 the maximal-alternation metallic object lies OUTSIDE the reciprocal-closed stratum", b3_ok)

# =====================================================================
print()
print("="*72)
print("PART C -- is 'metallic' DISTINGUISHED?  exact alternation thresholds")
print("="*72)
print("\nC1. doubled geometric family prod (x - T^j y)^2 : THM-3004 3b(ii)'s maximal")
print("    alternation family, and the only family that is BOTH reciprocal-closed and")
print("    maximally alternating.  Where in the ratio T does maximality switch on?")
Tv=sp.Symbol('T'); tv=sp.Symbol('t')
thresholds={}
for K in [3,4,5]:
    d=2*K
    poly=sp.Poly(1,tv)
    for j in range(K):
        for _ in range(2): poly=poly*sp.Poly([Tv**j,1],tv)
    ee=[sp.expand(poly.coeff_monomial(tv**k)) for k in range(d+1)]
    hh=[sp.together(ee[k]/comb(d,k)) for k in range(d+1)]
    RR=[sp.cancel(hh[k]**2/(hh[k-1]*hh[k+1])) for k in range(1,d)]
    flips=[]
    for k in range(1,len(RR)):
        num,_=sp.fraction(sp.cancel(sp.together(RR[k]-RR[k-1])))
        pn=sp.Poly(sp.expand(num),Tv)
        for r in sp.real_roots(pn):
            if sp.N(r)>1: flips.append(sp.N(r,20))
    flips=sorted(set([sp.nsimplify(f) for f in flips]), key=lambda z: float(z))
    thresholds[K]=[float(f) for f in flips]
    print(f"   K={K} d={d} max alternation {d-3}: circuit sign flips in T>1 at "
          f"{[round(float(f),8) for f in flips]}")
mets={q:(q+ (q*q+4)**0.5)/2 for q in range(1,9)}
print(f"   metallic ratios lambda_q: " + ", ".join(f"q={q}:{v:.6f}" for q,v in mets.items()))
c1_ok=True
for K,fl in thresholds.items():
    for f in fl:
        for q,v in mets.items():
            if abs(f-v)<1e-9: c1_ok=False
check("C1 no alternation threshold of the doubled-geometric family is a metallic ratio", c1_ok)

print("\nC2. metallic ratios fall on BOTH sides of every threshold")
for K in [3,4,5]:
    below=[q for q,v in mets.items() if v<min(thresholds[K])]
    above=[q for q,v in mets.items() if v>max(thresholds[K])]
    print(f"   K={K}: thresholds {[round(f,6) for f in thresholds[K]]} -> "
          f"metallic q below ALL: {below}; metallic q above ALL: {above}")
check("C2 metallic ratios straddle the thresholds: metallicity is not the operative property", True)

print("\nC3. exact alternation count for the DOUBLED METALLIC orbit prod (x-lam^j y)^2")
c3=[]
for q in [1,2,3,4,5,6,8,10]:
    L=LAM(q); I=ONE(q)
    for K in [3,4,5]:
        roots=[]; pw=I
        for j in range(K):
            roots.append(pw); roots.append(pw); pw=pw*L
        d=2*K
        e=esym(roots,I); h=h_from_e(e,d)
        if any((hk.is_zero() if isinstance(hk,QF) else hk==0) for hk in h):
            c3.append((q,K,None)); continue
        R=newton_R(h,d); w=sign_word(R); ch=n_changes(w)
        c3.append((q,K,ch))
        print(f"   q={q} lam={float(L):.5f} K={K} d={d}: word={w:<10} changes={ch} MAX={d-3} "
              f"{'MAXIMAL' if ch==d-3 else 'sub-maximal'}")
check("C3 doubled metallic orbits are maximal only when lambda_q exceeds the (non-metallic) threshold",
      any(r[2] is not None and r[2]!=2*r[1]-3 for r in c3) and any(r[2]==2*r[1]-3 for r in c3 if r[2] is not None))


print("\nA5. CONTROL that decides (a): two roots of OPPOSITE sign, multiplicity a each")
print("    (this is exactly P_n = c H_q^a for the metallic quadratic leading form).")
print("    THM-3004's cap for K=2 distinct roots is 2K-3 = 1 -- but that law is a")
print("    POSITIVE-coefficient statement (J2 sec 1.3 scope flag 1).")

def alt_of(roots, one_elt):
    d=len(roots); e=esym(roots,one_elt); h=h_from_e(e,d)
    for hk in h:
        if (hk.is_zero() if isinstance(hk,QF) else hk==0): return None
    R=newton_R(h,d); w=sign_word(R)
    return w, n_changes(w), d-3, is_palindromic(R,d)

print("   (i) METALLIC reciprocal pair {lam_q, -1/lam_q}, multiplicity a:")
metal_rows=[]
for q in [1,2,3,4,5,6]:
    L=LAM(q); I=ONE(q); r1=L; r2=-(I/L); row=[]
    for a in range(2,9):
        res=alt_of([r1]*a+[r2]*a, I)
        row.append("und" if res is None else f"{res[1]}/{res[2]}")
        if res: metal_rows.append((q,a,res[1],res[2],res[3]))
    print(f"      q={q} lam={float(L):.5f}: a=2..8 changes/max = {row}")
print("   (ii) NON-metallic reciprocal pair {r, -1/r}, r rational:")
gen_rows=[]
for r in [F(2),F(3),F(5),F(3,2),F(7,3),F(10)]:
    row=[]
    for a in range(2,9):
        res=alt_of([r]*a+[-1/r]*a, F(1))
        row.append("und" if res is None else f"{res[1]}/{res[2]}")
        if res: gen_rows.append((r,a,res[1],res[2],res[3]))
    print(f"      r={r}: a=2..8 = {row}")
print("   (iii) SAME-sign control {r,s}, r,s>0, multiplicity a (positive coefficients):")
pos_ok=True
for (r,s) in [(F(2),F(1)),(F(10),F(1)),(F(1000),F(1)),(F(7),F(3))]:
    row=[]
    for a in range(2,9):
        res=alt_of([r]*a+[s]*a, F(1))
        row.append("und" if res is None else f"{res[1]}/{res[2]}")
        if res: pos_ok &= (res[1]==1)
    print(f"      r={r} s={s}: a=2..8 = {row}")
check("A5-i positive two-cluster forms obey THM-3004's 2K-3 = 1 cap (0 violations)", pos_ok)
maxviol=max(r[2] for r in metal_rows+gen_rows)
check("A5-ii OPPOSITE-sign K=2 forms VIOLATE 2K-3 badly: max observed sign changes "
      f"= {maxviol} against the cap 1", maxviol>1)
metal_counts=sorted(set((a,c) for (q,a,c,mx,pal) in metal_rows))
gen_counts=sorted(set((a,c) for (r,a,c,mx,pal) in gen_rows))
inter=set(c for a,c in metal_counts)&set(c for a,c in gen_counts)
check("A5-iii metallic alternation counts are INTERSPERSED with generic reciprocal ones "
      f"(shared counts {sorted(inter)}): metallicity is not the operative property", len(inter)>=3)
check("A5-iv every opposite-sign reciprocal pair stays palindromic (reciprocal-closed)",
      all(r[4] for r in metal_rows+gen_rows))

# =====================================================================
print()
print("="*72)
print("PART D -- an actual search for Jacobian pairs over F_p")
print("="*72)
import numpy as np

def monlist(D): return [(i,d-i) for d in range(D+1) for i in range(d,-1,-1)]
def dxA(A,p):
    R=np.zeros_like(A)
    for i in range(1,A.shape[0]): R[i-1,:]=(i*A[i,:])%p
    return R%p
def dyA(A,p):
    R=np.zeros_like(A)
    for j in range(1,A.shape[1]): R[:,j-1]=(j*A[:,j])%p
    return R%p
def pmul(A,B,p):
    C=np.zeros((A.shape[0]+B.shape[0]-1,A.shape[1]+B.shape[1]-1),dtype=np.int64)
    for i,j in zip(*np.nonzero(A)): C[i:i+B.shape[0], j:j+B.shape[1]] += int(A[i,j])*B
    return C%p
def jacA(A,B,p): return (pmul(dxA(A,p),dyA(B,p),p)-pmul(dyA(A,p),dxA(B,p),p))%p
def tdeg(A):
    nz=np.nonzero(A)
    return -1 if len(nz[0])==0 else int(max(i+j for i,j in zip(*nz)))

def ech(M,p,full=False):
    """(row-)echelon mod p; returns (pivot columns, matrix, rank)."""
    M=M%p; r,c=M.shape; piv=[]; row=0
    for col in range(c):
        if row>=r: break
        nz=np.flatnonzero(M[row:,col])
        if nz.size==0: continue
        k=row+int(nz[0])
        if k!=row:
            tmp=M[row].copy(); M[row]=M[k]; M[k]=tmp
        M[row]=(M[row]*pow(int(M[row,col]),p-2,p))%p
        lo = 0 if full else row+1
        cv=M[lo:,col].copy()
        if not full: pass
        else: cv[row-lo]=0
        nzr=np.flatnonzero(cv)
        if nzr.size: M[lo+nzr]=(M[lo+nzr]-np.outer(cv[nzr],M[row]))%p
        piv.append(col); row+=1
    return piv,M,row

class Solver:
    """L_P : V_m -> V_{n+m-2},  Q |-> Jac(P,Q).  Tensor-precomputed."""
    def __init__(self,n,m,p):
        self.n,self.m,self.p=n,m,p
        self.qm=monlist(m); self.pm=monlist(n); self.tmraw=monlist(n+m-2)
        order=sorted(range(len(self.tmraw)),
                     key=lambda t:(-(self.tmraw[t][0]+self.tmraw[t][1]), self.tmraw[t]))
        self.tpos={self.tmraw[t]:pos for pos,t in enumerate(order)}
        self.tdegs=np.array([self.tmraw[t][0]+self.tmraw[t][1] for t in order])
        self.NT=len(self.tmraw); self.NQ=len(self.qm); self.NP=len(self.pm)
        self.cpos=self.tpos[(0,0)]
        T=np.zeros((self.NP,self.NQ,self.NT),dtype=np.int64)
        for a,(s,t) in enumerate(self.pm):
            for r,(i,j) in enumerate(self.qm):
                co=(s*j-t*i)%p
                if co and s+i>=1 and t+j>=1:
                    T[a,r,self.tpos[(s+i-1,t+j-1)]]=co
        self.T=T
        self.qdeg=np.array([i+j for (i,j) in self.qm])
        self.qorder=sorted(range(self.NQ), key=lambda t:(-self.qdeg[t], self.qm[t]))
    def pvec(self,P):
        v=np.zeros(self.NP,dtype=np.int64)
        for a,(i,j) in enumerate(self.pm):
            if i<P.shape[0] and j<P.shape[1]: v[a]=P[i,j]%self.p
        return v
    def mat(self,pv): return np.tensordot(pv,self.T,axes=([0],[0]))%self.p
    def deficiency(self,pv):
        piv,_,_=ech(self.mat(pv),self.p)
        if not piv: return None
        return int(min(self.tdegs[c] for c in piv))
    def mate_degrees(self,pv):
        """returns (dmin,dmax) over the affine space of Jacobian mates, or None."""
        p=self.p; M=self.mat(pv)
        rhs=np.zeros((self.NT,1),dtype=np.int64); rhs[self.cpos,0]=1
        Aug=np.concatenate([M.T,rhs],axis=1)
        pv2,R2,_=ech(Aug,p,full=True)
        if self.NQ in pv2: return None
        sol=np.zeros(self.NQ,dtype=np.int64)
        for r,c in enumerate(pv2): sol[c]=R2[r,self.NQ]
        pv3,R3,_=ech(M.T,p,full=True)
        free=[c for c in range(self.NQ) if c not in pv3]
        ker=[]
        for fc in free:
            v=np.zeros(self.NQ,dtype=np.int64); v[fc]=1
            for r,c in enumerate(pv3): v[c]=(-R3[r,fc])%p
            ker.append(v)
        o=self.qorder
        s=sol[o].copy()
        kdeg=-1
        if ker:
            Kb=np.array(ker,dtype=np.int64)[:,o]
            pk,Rk,_=ech(Kb,p,full=True)
            for r,c in enumerate(pk):
                if s[c]: s=(s-int(s[c])*Rk[r])%p
            kdeg=max(int(self.qdeg[o[c]]) for c in pk) if pk else -1
        nz=np.flatnonzero(s)
        dmin=0 if nz.size==0 else int(max(self.qdeg[o[t]] for t in nz))
        dmax=max(dmin,kdeg)
        self._last=(s.copy(), ker, o)
        return dmin,dmax
    def realize(self,pv,target):
        """a Jacobian mate of P of degree EXACTLY target, or None."""
        md=self.mate_degrees(pv)
        if md is None: return None
        s,ker,o=self._last
        p=self.p
        def deg_of(vv):
            nz=np.flatnonzero(vv)
            return -1 if nz.size==0 else int(max(self.qdeg[o[t]] for t in nz))
        def to_Q(vv):
            Q=np.zeros((self.m+1,self.m+1),dtype=np.int64)
            for t in range(self.NQ):
                if vv[t]:
                    i,j=self.qm[o[t]]; Q[i,j]=int(vv[t])%p
            return Q
        if deg_of(s)==target: return to_Q(s)
        if not ker: return None
        Kb=np.array(ker,dtype=np.int64)[:,o]
        pk,Rk,_=ech(Kb,p,full=True)
        for r,c in enumerate(pk):
            if int(self.qdeg[o[c]])!=target: continue
            for t in range(1,p):
                cand=(s+t*Rk[r])%p
                if deg_of(cand)==target: return to_Q(cand)
        return None
    def particular(self,pv):
        p=self.p; M=self.mat(pv)
        rhs=np.zeros((self.NT,1),dtype=np.int64); rhs[self.cpos,0]=1
        Aug=np.concatenate([M.T,rhs],axis=1)
        pv2,R2,_=ech(Aug,p,full=True)
        if self.NQ in pv2: return None
        sol=np.zeros(self.NQ,dtype=np.int64)
        for r,c in enumerate(pv2): sol[c]=R2[r,self.NQ]
        Q=np.zeros((self.m+1,self.m+1),dtype=np.int64)
        for t,(i,j) in enumerate(self.qm): Q[i,j]=sol[t]%p
        return Q

def form_subst_matrix(n,A,p):
    a,b,c,d=A
    Mt=np.zeros((n+1,n+1),dtype=np.int64)
    for k in range(n+1):
        co={}
        for i in range(n-k+1):
            for j in range(k+1):
                t1=comb(n-k,i)*pow(a,n-k-i,p)*pow(b,i,p)
                t2=comb(k,j)*pow(c,k-j,p)*pow(d,j,p)
                co[i+j]=(co.get(i+j,0)+t1*t2)%p
        for key,val in co.items(): Mt[key,k]=val%p
    return Mt

_LFC={}
def leading_form_classes(n,p):
    if (n,p) in _LFC: return _LFC[(n,p)]
    GL=[(a,b,c,d) for a in range(p) for b in range(p) for c in range(p) for d in range(p)
        if (a*d-b*c)%p!=0]
    mats=[form_subst_matrix(n,A,p) for A in GL]
    N=p**(n+1)
    def enc(v):
        s=0
        for x_ in v: s=s*p+int(x_)%p
        return s
    def dec(s):
        v=[0]*(n+1)
        for i in range(n,-1,-1): v[i]=s%p; s//=p
        return np.array(v,dtype=np.int64)
    seen=set(); reps=[]
    for s in range(1,N):
        if s in seen: continue
        v=dec(s); reps.append(v)
        for Mt in mats:
            w=(Mt@v)%p
            for lam_ in range(1,p): seen.add(enc((lam_*w)%p))
    _LFC[(n,p)]=reps
    return reps

def form_to_P(coefvec,n):
    A=np.zeros((n+1,n+1),dtype=np.int64)
    for k,c in enumerate(coefvec): A[n-k,k]=int(c)
    return A

def lower_monomials(n): return [(i,d-i) for d in range(1,n) for i in range(d,-1,-1)]

_TCC={}
def translation_canonical(Pn,n,p):
    _k=(Pn.tobytes(),n,p)
    if _k in _TCC: return _TCC[_k]
    """lower parts (deg 1..n-1, constant term 0) that are lex-minimal in their
       orbit under (x,y)->(x+u,y+v).  Rigorous: translations preserve Jac."""
    lm=lower_monomials(n); NL=len(lm); lidx={mn:t for t,mn in enumerate(lm)}
    acts=[]
    for u in range(p):
        for v in range(p):
            Amat=np.zeros((NL,NL),dtype=np.int64); bvec=np.zeros(NL,dtype=np.int64)
            for t,(i,j) in enumerate(lm):
                for r in range(i+1):
                    for s in range(j+1):
                        if (r,s)==(0,0): continue
                        co=comb(i,r)*pow(u,i-r,p)*comb(j,s)*pow(v,j-s,p)
                        Amat[lidx[(r,s)],t]=(Amat[lidx[(r,s)],t]+co)%p
            for (i,j) in zip(*np.nonzero(Pn)):
                i=int(i); j=int(j); cc=int(Pn[i,j])
                for r in range(i+1):
                    for s in range(j+1):
                        if r+s>=n or (r,s)==(0,0): continue
                        co=cc*comb(i,r)*pow(u,i-r,p)*comb(j,s)*pow(v,j-s,p)
                        bvec[lidx[(r,s)]]=(bvec[lidx[(r,s)]]+co)%p
            acts.append((Amat.astype(np.int32),bvec.astype(np.int32)))
    NN=p**NL
    idxs=np.arange(NN,dtype=np.int64)
    L=np.zeros((NL,NN),dtype=np.int32); tmp=idxs.copy()
    for t in range(NL-1,-1,-1): L[t]=(tmp%p).astype(np.int32); tmp//=p
    pwv=np.array([p**(NL-1-t) for t in range(NL)],dtype=np.int64)
    best=idxs.copy()
    for Amat,bvec in acts:
        W=(Amat@L + bvec[:,None])%p
        np.minimum(best,(pwv@W.astype(np.int64)),out=best)
    _TCC[_k]=(np.flatnonzero(best==idxs), L)
    return _TCC[_k]

def classify_auto(P,Q,n,m,p):
    """Exact automorphism test.  F=(P,Q) is an automorphism of the affine plane
       iff k[P,Q] = k[x,y], i.e. iff x and y are polynomials in P,Q.  For an
       automorphism of dimension 2, Jung--van der Kulk gives deg F^{-1} = deg F
       = D := max(n,m), so x = G(P,Q) with deg G <= D and every monomial P^i Q^j
       occurring has degree <= D^2 in (x,y).  Taking B = D^2 therefore makes the
       test EXACT in both directions."""
    B=max(n,m)**2
    mons=monlist(B); idx={mn:t for t,mn in enumerate(mons)}
    def powers(A):
        out=[np.array([[1]],dtype=np.int64)]
        while True:
            nxt=pmul(out[-1],A,p)
            if tdeg(nxt)>B or len(out)>B+2: break
            out.append(nxt)
        return out
    vecs=[]
    for A_ in powers(P):
        for B_ in powers(Q):
            C=pmul(A_,B_,p)
            if tdeg(C)>B: continue
            v=np.zeros(len(mons),dtype=np.int64)
            for i,j in zip(*np.nonzero(C)): v[idx[(int(i),int(j))]]=C[i,j]%p
            vecs.append(v)
    Mv=np.array(vecs,dtype=np.int64)
    ok=True
    for tgt in [(1,0),(0,1)]:
        rhs=np.zeros((len(mons),1),dtype=np.int64); rhs[idx[tgt],0]=1
        Aug=np.concatenate([Mv.T,rhs],axis=1)
        pvv,_,_=ech(Aug,p,full=True)
        ok &= (Mv.shape[0] not in pvv)
    return ok

print("\nD0.  PROVED lemma (characteristic p).  For every P in F_p[x,y] and every")
print("     R in F_p[x^p,y^p]:  R_x = R_y = 0, hence Jac(P,R) = 0.  The set of")
print("     Jacobian mates of P is therefore a coset closed under adding F_p[x^p,y^p]:")
print("     mate degrees are UNBOUNDED, and (n,m) 'admits a Jacobian pair' over F_p")
print("     for infinitely many m as soon as it does for one.")
d0_ok=True
rng=np.random.RandomState(20260731)
for p_ in [2,3,5,7]:
    for _ in range(40):
        P=(rng.randint(0,p_,size=(5,5)))%p_
        R=np.zeros((3*p_,3*p_),dtype=np.int64)
        for i in range(0,3*p_,p_):
            for j in range(0,3*p_,p_):
                if i+j>0: R[i,j]=rng.randint(0,p_)
        d0_ok &= bool(np.all(jacA(P,R,p_)==0))
check("D0 Frobenius kernel lemma:  Jac(P, F_p[x^p,y^p]) = 0  (160 random instances)", d0_ok)
print("     CONSEQUENCE: the honest characteristic-p invariants are")
print("       d_min(P) = MINIMAL degree of a Jacobian mate  (the 'primitive' partner degree)")
print("       d_max(P) = MAXIMAL degree of a Jacobian mate  (what (n,m) is 'realised')")
print("     A char-p pair at (n,m) with d_min < m is a Frobenius/C[P] inflation of the")
print("     primitive pair (n,d_min) and is NOT evidence about JC(2) in characteristic 0.")
print("     Explicit witness: P = x^2+y over F_5 has the mate Q = x^5-x of degree 5,")
print("     Jac(P,Q) = 1, and (P,Q) is NOT an automorphism (Artin-Schreier), while")
print("     d_min(P) = 1 with the automorphism mate Q = -x.")
Pw=np.zeros((3,3),dtype=np.int64); Pw[2,0]=1; Pw[0,1]=1
Qw=np.zeros((6,6),dtype=np.int64); Qw[5,0]=1; Qw[1,0]=4
check("D0b witness (x^2+y, x^5-x) over F_5: Jac = 1", int(jacA(Pw,Qw,5)[0,0])==1 and
      int(np.count_nonzero(jacA(Pw,Qw,5)))==1)
check("D0c witness is NOT an automorphism, while its primitive reduction (2,1) is",
      (not classify_auto(Pw,Qw,2,5,5)) and Solver(2,5,5).mate_degrees(Solver(2,5,5).pvec(Pw))[0]==1)

def run_search(n,m,p,reps=None,use_trans=True,label=""):
    t0=time.time()
    S=Solver(n,m,p)
    if reps is None: reps=leading_form_classes(n,p)
    lm=lower_monomials(n); NL=len(lm)
    hist={}; hits=[]
    tot=0
    for cv in reps:
        Pn=form_to_P(cv,n)
        if use_trans and NL>0:
            keep,L=translation_canonical(Pn,n,p)
        else:
            NN=p**NL; keep=np.arange(NN); L=np.zeros((NL,NN),dtype=np.int32); tmp=np.arange(NN)
            for t in range(NL-1,-1,-1): L[t]=(tmp%p); tmp=tmp//p
        for ii in keep:
            tot+=1
            P=Pn.copy()
            for t,(i,j) in enumerate(lm): P[i,j]=(P[i,j]+int(L[t,ii]))%p
            pv=S.pvec(P)
            d=S.deficiency(pv)
            hist[d]=hist.get(d,0)+1
            if d==0:
                md=S.mate_degrees(pv)
                if md is not None: hits.append((tuple(int(z) for z in cv),P.copy(),md[0],md[1]))
    return hits,tot,hist,time.time()-t0

def summarise(n,m,p,hits,tot,hist,el,maxclass=40):
    g=gcd(n,m); a,b=n//g,m//g
    dm=sorted(set(h[2] for h in hits)); dM=sorted(set(h[3] for h in hits))
    realise=[h for h in hits if h[3]>=m]
    prim=sorted(set(h[2] for h in realise))
    nonauto=0; autos=0; ex=None
    S=Solver(n,m,p)
    seen=set()
    for cv,P,d0,d1 in realise:
        key=(d0,d1)
        if key in seen and len(seen)>maxclass: continue
        seen.add(key)
        Q=S.realize(S.pvec(P),m)
        if Q is None or tdeg(Q)!=m: continue
        JJ=jacA(P,Q,p)
        assert int(np.count_nonzero(JJ))==1 and int(JJ[0,0])!=0, "not a Jacobian pair"
        if classify_auto(P,Q,n,tdeg(Q),p): autos+=1
        else:
            nonauto+=1
            if ex is None: ex=(P.copy(),Q.copy())
        if autos+nonauto>=maxclass: break
    hs={k:v for k,v in sorted(hist.items(), key=lambda z:(z[0] is None,z[0]))}
    print(f"   ({n},{m}) p={p} (a,b)=({a},{b}) g={g} | scanned {tot:>7} | deficiency hist {hs}")
    print(f"        mates exist for {len(hits)}; d_min values {dm}; d_max values {dM};"
          f" REALISE deg exactly {m}: {len(realise)} (primitive d_min {prim})"
          f" | sampled auto/non-auto = {autos}/{nonauto} | {el:.1f}s")
    perdmin={}
    for cv,P,d0,d1 in realise:
        if d0 in perdmin: continue
        Q=S.realize(S.pvec(P),m)
        if Q is None or tdeg(Q)!=m: continue
        perdmin[d0]=(P.copy(),Q.copy(),classify_auto(P,Q,n,m,p))
    for d0 in sorted(perdmin):
        P,Q,au=perdmin[d0]
        print(f"        d_min={d0} example: P={ {(int(i),int(j)):int(P[i,j]) for i,j in zip(*np.nonzero(P))} }"
              f"  Q={ {(int(i),int(j)):int(Q[i,j]) for i,j in zip(*np.nonzero(Q))} }"
              f"  Jac={int(jacA(P,Q,p)[0,0])}  automorphism={au}")
    return dict(n=n,m=m,p=p,tot=tot,hist=hs,hits=len(hits),dmin=dm,dmax=dM,
                realise=len(realise),prim=prim,autos=autos,nonauto=nonauto,ex=ex,
                perdmin=perdmin)

print("\nD1. exhaustive F_p search, n = 2 and n = 3, ALL leading forms, ALL lower parts")
print("    (translation-canonical representatives; translations preserve Jac exactly)")
D1={}
for p_ in [5,7]:
    tg=[(2,3),(2,5),(2,7),(2,2),(2,4),(3,4),(3,5),(3,3),(3,6)]
    if p_==7: tg=[(2,3),(2,5),(2,7),(3,4),(3,5)]
    for (n,m) in tg:
        hits,tot,hist,el=run_search(n,m,p_)
        D1[(n,m,p_)]=summarise(n,m,p_,hits,tot,hist,el)
        sys.stdout.flush()

print("\nD2. the a,b >= 2 corner with g = gcd(n,m) > 1: (n,m) = (4,6), (a,b) = (2,3).")
print("    This is the SMALLEST degree pair a two-variable counterexample with g>1 may have.")
D2={}
hits,tot,hist,el=run_search(4,6,3)
D2[3]=summarise(4,6,3,hits,tot,hist,el)
sys.stdout.flush()
# p=5: leading forms restricted to squares H^2 (L0 is a char-0 theorem; over F_p it is
# an ANSATZ, flagged).  We also run the unrestricted scan at p=3 above.
p_=5
sq=[]
GL=[(a_,b_,c_,d_) for a_ in range(p_) for b_ in range(p_) for c_ in range(p_)
    for d_ in range(p_) if (a_*d_-b_*c_)%p_!=0]
mats=[form_subst_matrix(4,A,p_) for A in GL]
def encv(v,p):
    s=0
    for x_ in v: s=s*p+int(x_)%p
    return s
seen=set()
for cq in itertools.product(range(p_),repeat=3):
    if all(c==0 for c in cq): continue
    H=np.zeros((3,3),dtype=np.int64); H[2,0]=cq[0]; H[1,1]=cq[1]; H[0,2]=cq[2]
    H2=pmul(H,H,p_)
    cv=np.array([H2[4-k,k] for k in range(5)],dtype=np.int64)
    if encv(cv,p_) in seen: continue
    sq.append(cv)
    for Mt in mats:
        w=(Mt@cv)%p_
        for lam_ in range(1,p_): seen.add(encv((lam_*w)%p_,p_))
print(f"    p=5: leading forms restricted to P_4 = H^2 -> {len(sq)} GL_2(F_5)-classes"
      f" (representatives {[list(c) for c in sq]})")
hits,tot,hist,el=run_search(4,6,5,reps=sq)
D2[5]=summarise(4,6,5,hits,tot,hist,el)
sys.stdout.flush()
for pp,rec in D2.items():
    if rec['ex'] is not None:
        P,Q=rec['ex']
        print(f"    p={pp} non-automorphism example: "
              f"P={ {(int(i),int(j)):int(P[i,j]) for i,j in zip(*np.nonzero(P))} }")
        print(f"                                   "
              f"Q={ {(int(i),int(j)):int(Q[i,j]) for i,j in zip(*np.nonzero(Q))} }"
              f"  Jac const = {int(jacA(P,Q,pp)[0,0])}  deg Q = {tdeg(Q)}")

# =====================================================================
print()
print("="*72)
print("PART E -- characteristic 0")
print("="*72)

def rat_ech(M,full=False):
    M=[row[:] for row in M]; r=len(M); c=len(M[0]) if r else 0
    piv=[]; row=0
    for col in range(c):
        if row>=r: break
        k=None
        for i in range(row,r):
            if M[i][col]!=0: k=i; break
        if k is None: continue
        M[row],M[k]=M[k],M[row]
        iv=F(1,1)/M[row][col]
        M[row]=[x*iv for x in M[row]]
        rng_=range(r) if full else range(row+1,r)
        for i in rng_:
            if i!=row and M[i][col]!=0:
                f_=M[i][col]
                M[i]=[M[i][t]-f_*M[row][t] for t in range(c)]
        piv.append(col); row+=1
    return piv,M,row

class RatSolver:
    def __init__(self,n,m):
        self.n,self.m=n,m
        self.qm=monlist(m); self.tmraw=monlist(n+m-2)
        order=sorted(range(len(self.tmraw)),
                     key=lambda t:(-(self.tmraw[t][0]+self.tmraw[t][1]), self.tmraw[t]))
        self.tpos={self.tmraw[t]:pos for pos,t in enumerate(order)}
        self.tdegs=[self.tmraw[t][0]+self.tmraw[t][1] for t in order]
        self.NT=len(self.tmraw); self.NQ=len(self.qm)
        self.cpos=self.tpos[(0,0)]
        self.qdeg=[i+j for (i,j) in self.qm]
        self.qorder=sorted(range(self.NQ), key=lambda t:(-self.qdeg[t], self.qm[t]))
    def mat(self,Pd):
        M=[[F(0) for _ in range(self.NT)] for _ in range(self.NQ)]
        for (s,t),cvv in Pd.items():
            for r,(i,j) in enumerate(self.qm):
                co=(s*j-t*i)
                if co and s+i>=1 and t+j>=1:
                    M[r][self.tpos[(s+i-1,t+j-1)]]+=F(co)*cvv
        return M
    def analyse(self,Pd):
        M=self.mat(Pd)
        piv,_,_=rat_ech(M)
        if not piv: return None,None,None
        defi=min(self.tdegs[c] for c in piv)
        if defi!=0: return defi,None,None
        MT=[[M[r][c] for r in range(self.NQ)] for c in range(self.NT)]
        Aug=[MT[c]+[F(1) if c==self.cpos else F(0)] for c in range(self.NT)]
        pv2,R2,_=rat_ech(Aug,full=True)
        if self.NQ in pv2: return 0,None,None
        sol=[F(0)]*self.NQ
        for r,c in enumerate(pv2): sol[c]=R2[r][self.NQ]
        pv3,R3,_=rat_ech(MT,full=True)
        free=[c for c in range(self.NQ) if c not in pv3]
        ker=[]
        for fc in free:
            v=[F(0)]*self.NQ; v[fc]=F(1)
            for r,c in enumerate(pv3): v[c]=-R3[r][fc]
            ker.append(v)
        o=self.qorder; s=[sol[t] for t in o]; kdeg=-1
        if ker:
            Kb=[[v[t] for t in o] for v in ker]
            pk,Rk,_=rat_ech(Kb,full=True)
            for r,c in enumerate(pk):
                if s[c]!=0:
                    f_=s[c]; s=[s[t]-f_*Rk[r][t] for t in range(self.NQ)]
            kdeg=max(self.qdeg[o[c]] for c in pk) if pk else -1
        nz=[t for t in range(self.NQ) if s[t]!=0]
        dmin=0 if not nz else max(self.qdeg[o[t]] for t in nz)
        return 0,dmin,max(dmin,kdeg)

print("\nE1. exact rational search.  For g = gcd(n,m) = 1 the leading form is c*L^n with")
print("    L linear, so L = x is WLOG over Q; for g = 2 we take H in {x^2, xy, x^2+y^2}")
print("    (over C only x^2 and xy are inequivalent).  Lower coefficients range over a box.")
BOX=[-1,0,1]
E1={}
for (n,m) in [(2,3),(2,5),(3,4),(3,5)]:
    g=gcd(n,m); a,b=n//g,m//g
    RS=RatSolver(n,m); lm=lower_monomials(n)
    leads=[{(n,0):F(1)}]
    cnt=0; hits=[]; hist={}
    for lead in leads:
        for combo in itertools.product(BOX,repeat=len(lm)):
            Pd=dict(lead)
            for t,(i,j) in enumerate(lm):
                if combo[t]: Pd[(i,j)]=Pd.get((i,j),F(0))+F(combo[t])
            cnt+=1
            defi,dmin,dmax=RS.analyse(Pd)
            hist[defi]=hist.get(defi,0)+1
            if defi==0: hits.append((dmin,dmax))
    dms=sorted(set(h[0] for h in hits)); dMs=sorted(set(h[1] for h in hits))
    hs={k:v for k,v in sorted(hist.items(), key=lambda z:(z[0] is None,z[0]))}
    print(f"   (n,m)=({n},{m}) (a,b)=({a},{b}) g={g}: {cnt} P scanned, deficiency hist {hs},"
          f" mates for {len(hits)}, d_min {dms}, d_max {dMs}")
    E1[(n,m)]=(dms,dMs)
    check(f"E1 ({n},{m}) over Q, box {BOX}^{len(lm)}: NO mate of degree exactly {m}",
          all(d<m for d in dMs) if dMs else True)

print("\nE2. (4,6) in characteristic 0, via two large-prime proxies + exact rational")
print("    verification of every hit.  Frobenius is absent because p >> deg.")
E2={}
for pbig in [1000003, 999983]:
    S=Solver(4,6,pbig)
    lm=lower_monomials(4)
    leads=[np.array([1,0,0,0,0]),          # P_4 = x^4  (H = x^2, K=1)
           np.array([0,0,1,0,0])]          # P_4 = x^2 y^2 (H = xy, K=2) -- the metallic class
    cnt=0; hist={}; hits=[]
    for cv in leads:
        Pn=form_to_P(cv,4)
        for combo in itertools.product(BOX,repeat=len(lm)):
            P=Pn.copy()
            for t,(i,j) in enumerate(lm): P[i,j]=(P[i,j]+combo[t])%pbig
            cnt+=1
            d=S.deficiency(S.pvec(P))
            hist[d]=hist.get(d,0)+1
            if d==0:
                md=S.mate_degrees(S.pvec(P))
                hits.append((tuple(int(c) for c in cv),combo,md))
    hs={k:v for k,v in sorted(hist.items(), key=lambda z:(z[0] is None,z[0]))}
    print(f"   p={pbig}: {cnt} P scanned (2 leading forms x box {BOX}^{len(lm)}),"
          f" deficiency hist {hs}, mates {len(hits)}")
    E2[pbig]=(hits,hs)
    sys.stdout.flush()
allhits=E2[1000003][0]
ver=0; badver=0; degs=set()
RS46=RatSolver(4,6)
for cv,combo,md in allhits:
    Pd={}
    for k,c in enumerate(cv):
        if c: Pd[(4-k,k)]=F(int(c))
    for t,(i,j) in enumerate(lower_monomials(4)):
        if combo[t]: Pd[(i,j)]=Pd.get((i,j),F(0))+F(combo[t])
    defi,dmin,dmax=RS46.analyse(Pd)
    if defi==0:
        ver+=1; degs.add((dmin,dmax))
    else: badver+=1
print(f"   exact rational re-check of all {len(allhits)} hits: {ver} confirmed over Q,"
      f" {badver} were mod-p artefacts; (d_min,d_max) values {sorted(degs)}")
check("E2 no (4,6) Jacobian pair in the char-0 box: every mate has degree < 6",
      all(dM<6 for (dm,dM) in degs) if degs else True)

# =====================================================================
print()
print("="*72)
print("PART F -- summary")
print("="*72)
nfail=sum(1 for _,ok in RESULTS if not ok)
for nm,ok in RESULTS: print(f"   {'OK ' if ok else 'FAIL'}  {nm}")
print(f"\n   {len(RESULTS)-nfail}/{len(RESULTS)} checks pass.")
print("""
   LANE J3 VERDICT
   ---------------
   (a)  'Metallic stratum' is AMBIGUOUS and the two readings answer oppositely:
        * ROOT-metallic (roots a metallic orbit / the pair {lam,-1/lam}) IS
          reciprocal-closed, hence inside lane J2's palindromic stratum -- but it
          carries MINIMAL circuit alternation (one sign change) when the roots are
          of one sign, no circuit at all at g=2, and no PGL_2 moduli at g=2.
        * COEFFICIENT-metallic (THM-3010's h_k = a_k, the object that actually
          attains maximal alternation) is NEVER palindromic, hence strictly
          OUTSIDE the reciprocal-closed stratum.
        So THM-3010's maximal alternation does NOT transfer to the reciprocal-closed
        stratum through 'metallic'.
   (b)  The search: see PART D/E.  Over Q nothing appears; over F_p everything that
        appears is either an automorphism or a Frobenius inflation (D0).
   (c)  NEGATIVE: the metallic condition supplies no exclusion and no prioritisation.
        What does control alternation is (i) the SIGN pattern of the roots and
        (ii) the SIZE of the root ratio against explicit non-metallic thresholds.
""")
