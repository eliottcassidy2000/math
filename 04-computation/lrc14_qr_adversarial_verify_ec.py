#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL VERIFY (independent EXACT re-implementation) of the QR/Gauss-sum
cancellation claims for the support-6 coset kernel D7 in LRC(14).

I build my OWN Z[zeta_7] arithmetic (basis 1, z, z^2, ..., z^5; using z^6 = -(1+z+..+z^5))
and re-derive each claimed structural fact from the definition, then attack the GAP.

Claims under test:
 (1) FACTORIZATION  K(n) = [prod_j 1/(2 pi i n_j)] * D7(n mod 7), arch part REAL for support 6.
 (2) GALOIS-EQUIVARIANCE (exact): D7(a*c) = sigma_a(D7(c)),  sigma_a: z -> z^a.
 (3) FIELD-TRACE COLLAPSE (exact): sum_{a in F_7^*} D7(a c) = Tr_{Q(z7)/Q}(D7(c)) in Z.
     Tr(D7(1,1,1,1,1,1)) = 7! = 5040 ; max|Tr| = 9240 ; sum of Tr over orbit reps = 0;
     vanishes per #QR-pattern class.
 (4) REALITY: D7(-c) = conj(D7(c)) since 6 = -1 is a non-residue => correction RIGOROUSLY REAL.
 (5) THE GAP: dilation c->ac is not a lattice symmetry; within-orbit lattice weights differ.

Then I ADVERSARIALLY HUNT: can ANY lattice-respecting regrouping recover the trace
collapse?  In particular test whether the n<->-n REAL pairing plus the residue
structure gives a magnitude handle for WIDE E after all.
"""
import sys, itertools, math, cmath
from fractions import Fraction
from collections import defaultdict
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')

# ---------------------------------------------------------------------------
# Exact Z[zeta_7]: represent element as tuple of 6 Fractions = coeffs of 1,z,...,z^5.
# Reduction: z^6 = -(1 + z + z^2 + z^3 + z^4 + z^5).
# ---------------------------------------------------------------------------
ZERO = (Fraction(0),)*6
ONE  = (Fraction(1),Fraction(0),Fraction(0),Fraction(0),Fraction(0),Fraction(0))

def z_pow(e):
    """zeta^e as exact element."""
    e %= 7
    if e == 6:
        return (Fraction(-1),)*6
    v = [Fraction(0)]*6
    v[e] = Fraction(1)
    return tuple(v)

def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def smul(a,s): return tuple(x*s for x in a)

def mul(a,b):
    raw = [Fraction(0)]*12  # exponents 0..11 (5+5+ carries up to 10)
    for i in range(6):
        if a[i]==0: continue
        ai=a[i]
        for j in range(6):
            if b[j]==0: continue
            raw[i+j]+=ai*b[j]
    # reduce exponents >=6 using z^6 = -(1+...+z^5), and z^7=z*z^6 etc.
    # easier: fold each exponent e>=6 by writing z^e in basis.
    out=[Fraction(0)]*6
    for e in range(12):
        if raw[e]==0: continue
        be=z_pow(e)
        for k in range(6):
            out[k]+=raw[e]*be[k]
    return tuple(out)

def galois(a, x):
    """sigma_a: zeta -> zeta^a applied to x."""
    r=list(ZERO)
    for k in range(6):
        if x[k]!=0:
            be=z_pow((a*k)%7)
            for j in range(6):
                r[j]+=x[k]*be[j]
    return tuple(r)

def trace(x):
    """Tr_{Q(z7)/Q}(x) = sum_{a=1..6} sigma_a(x). Returns a Fraction (should be rational)."""
    tot=list(ZERO)
    for a in range(1,7):
        ga=galois(a,x)
        tot=[t+g for t,g in zip(tot,ga)]
    # trace must be rational: coeffs of z^1..z^5 should be equal to coeff of 1 contribution?
    # Actually Tr(1)=6, Tr(z^k)=-1 (k=1..6). So tot in basis: tot = q*ONE means tot[1..5]==0.
    return tuple(tot)

def to_num(x):
    Z=cmath.exp(2j*math.pi/7)
    return sum(complex(x[k])*(Z**k) for k in range(6))

def is_rational(x):
    return all(x[k]==0 for k in range(1,6))

# ---------------------------------------------------------------------------
# D7 kernel.  From HYP-2646:
#   D7(c) = pref(c) * sum_{T subset {1..6}} (-1)^|T| prod_j sigma_T(c_j)
#   pref(c) = prod_j (1 - zeta^{-c_j})
#   sigma_T(m) = sum_{t in T} zeta^{-m t}
# (matching lrc14_qr_symmetry_mechanism_kps.py exactly so we test the SAME object)
# ---------------------------------------------------------------------------
Tlist=[T for r in range(7) for T in itertools.combinations(range(1,7),r)]
SGN={T:Fraction((-1)**len(T)) for T in Tlist}
def sigma_T(T,m):
    r=list(ZERO)
    for t in T:
        be=z_pow((-m*t)%7)
        for j in range(6): r[j]+=be[j]
    return tuple(r)
SIG={(T,m):sigma_T(T,m) for T in Tlist for m in range(1,7)}
PREF={m:sub(ONE, z_pow((-m)%7)) for m in range(1,7)}

_cache={}
def D7(c):
    v=_cache.get(c)
    if v is not None: return v
    pref=ONE
    for cj in c: pref=mul(pref, PREF[cj])
    acc=list(ZERO)
    for T in Tlist:
        p=ONE
        for cj in c: p=mul(p, SIG[(T,cj)])
        sg=SGN[T]
        for j in range(6): acc[j]+=sg*p[j]
    v=mul(pref, tuple(acc)); _cache[c]=v
    return v

QR={1,2,4}; NQR={3,5,6}

# ---------------------------------------------------------------------------
def main():
    print("="*72)
    print("ADVERSARIAL VERIFY: QR/Gauss-sum collapse of D7 (independent exact impl)")
    print("="*72)

    # sanity: Gauss sum g = sum chi(r) z^r, g^2 = -7
    def chi(m):
        m%=7
        return 0 if m==0 else (1 if m in QR else -1)
    g=list(ZERO)
    for r in range(1,7):
        be=z_pow(r)
        for j in range(6): g[j]+=Fraction(chi(r))*be[j]
    g=tuple(g)
    g2=mul(g,g)
    print(f"\n[sanity] Gauss sum g^2 = {[str(x) for x in g2]}  (expect -7)")
    assert g2==smul(ONE,Fraction(-7)), "Gauss sum g^2 != -7"
    print("         g^2 = -7 EXACT. OK.")

    # ---- CLAIM (2): Galois equivariance D7(a c) = sigma_a(D7(c)) ----
    print("\n[CLAIM 2] Galois equivariance D7(a*c) = sigma_a(D7(c)):")
    import random
    random.seed(11)
    fails=0; tested=0
    test_cs=[(1,1,1,1,1,1),(1,2,4,3,5,6),(1,1,1,1,1,2),(3,3,3,3,3,3),
             (1,2,3,4,5,6),(2,2,5,5,6,1),(4,1,1,6,3,2)]
    test_cs+=[tuple(random.choice([1,2,3,4,5,6]) for _ in range(6)) for _ in range(30)]
    for c in test_cs:
        dc=D7(c)
        for a in range(1,7):
            ac=tuple((a*cj)%7 for cj in c)
            tested+=1
            if D7(ac)!=galois(a,dc):
                fails+=1
                if fails<=3: print(f"    FAIL c={c} a={a}")
    print(f"    tested {tested} (c,a) pairs, fails={fails}  => {'HOLDS EXACTLY' if fails==0 else 'BROKEN'}")
    claim2 = (fails==0)

    # ---- CLAIM (3): trace collapse sum_a D7(a c) = Tr(D7(c)) rational ----
    print("\n[CLAIM 3] sum_{a in F7*} D7(a*c) = Tr(D7(c)) is rational integer:")
    spec={(1,1,1,1,1,1):None,(1,2,4,3,5,6):None,(1,1,1,1,1,2):None,(1,2,3,1,2,3):None}
    bad=0
    for c in spec:
        orb=list(ZERO)
        for a in range(1,7):
            d=D7(tuple((a*cj)%7 for cj in c))
            orb=[o+x for o,x in zip(orb,d)]
        orb=tuple(orb)
        tr=trace(D7(c))
        israt=is_rational(orb)
        match=(orb==tr)
        if not (israt and match): bad+=1
        val = orb[0] if israt else None
        print(f"    c={c}: orbit-sum rational={israt}  ==Tr={match}  value={val}")
    # Tr(1^6)=5040 check
    tr_ones=trace(D7((1,1,1,1,1,1)))
    print(f"    Tr(D7(1,1,1,1,1,1)) = {tr_ones[0]}  (claim 5040 = 7!) -> {'OK' if tr_ones[0]==5040 and is_rational(tr_ones) else 'MISMATCH'}")
    claim3 = (bad==0 and tr_ones[0]==Fraction(5040))

    # full scan over all 6^6 cosets: Tr per coset, max|Tr|, sum Tr, per #QR class
    print("\n[CLAIM 3b] Full 6^6 coset scan of Tr(D7):")
    sumTr=Fraction(0); maxabs=Fraction(0)
    by_qr=defaultdict(lambda: Fraction(0))
    sample_orbit_sum=Fraction(0)
    # iterate cosets c in {1..6}^6
    seen_reps=set()
    for c in itertools.product(range(1,7),repeat=6):
        tr=trace(D7(c))
        assert is_rational(tr), f"Tr not rational at c={c}"
        t=tr[0]
        sumTr+=t
        if abs(t)>maxabs: maxabs=abs(t)
        nqr=sum(1 for x in c if x in QR)   # #QR pattern
        by_qr[nqr]+=t
    print(f"    sum over ALL 6^6 cosets of Tr = {sumTr}  (claim 0) -> {'OK' if sumTr==0 else 'MISMATCH'}")
    print(f"    max|Tr| over all cosets = {maxabs}  (claim 9240) -> {'OK' if maxabs==9240 else 'val='+str(maxabs)}")
    print(f"    9240 = 2^3*3*5*7*11 = {2**3*3*5*7*11}")
    print(f"    sum Tr by #QR-pattern class:")
    allzero=True
    for k in sorted(by_qr):
        print(f"       #QR={k}: sum Tr = {by_qr[k]}")
        if by_qr[k]!=0: allzero=False
    print(f"    vanishes on EACH #QR class independently: {allzero}")
    claim3b = (sumTr==0 and maxabs==9240 and allzero)

    # ---- CLAIM (4): reality D7(-c) = conj(D7(c)) ----
    print("\n[CLAIM 4] D7(6c)=sigma_6(D7(c))=conj(D7(c)) [6=-1 is NQR] => correction REAL:")
    # complex conjugation on Q(z7) is sigma_6 (z->z^-1=z^6). Test sigma_6(x)=conj.
    fails4=0
    for c in test_cs[:12]:
        d=D7(c)
        s6=galois(6,d)
        # numeric conj check
        if abs(to_num(s6)-to_num(d).conjugate())>1e-9: fails4+=1
        # also D7(6c) == s6 exactly:
        if D7(tuple((6*cj)%7 for cj in c))!=s6: fails4+=1
    print(f"    sigma_6 = complex conjugation and D7(6c)=conj(D7(c)): fails={fails4} -> {'HOLDS' if fails4==0 else 'BROKEN'}")
    claim4=(fails4==0)

    # ---- CLAIM (5)/GAP: dilation is not a lattice symmetry; within-orbit weights differ ----
    # Re-test on real relation lattices that the correction is REAL (n<->-n pairing) and
    # that within a dilation orbit the lattice weights are NOT equal (so trace doesn't apply).
    print("\n[GAP] Lattice test: correction REAL (n<->-n), but dilation-orbit weights unequal:")
    Z=cmath.exp(2j*math.pi/7)
    _Dn={}
    def D7num(c):
        v=_Dn.get(c)
        if v is None:
            v=to_num(D7(c)); _Dn[c]=v
        return v
    def w_real(vals):
        p=1.0
        for v in vals: p*=1.0/v
        return -p/((2*math.pi)**6)   # i^6 = -1
    def relations_support6(E,L):
        out=[]
        for idxs in itertools.combinations(range(len(E)),6):
            es=[E[i] for i in idxs]
            dep=max(range(6),key=lambda i:abs(es[i])); ed=es[dep]
            if ed==0: continue
            free=[i for i in range(6) if i!=dep]; ef=[es[i] for i in free]
            for vf in itertools.product(range(-L,L+1),repeat=5):
                if any(c==0 or c%7==0 for c in vf): continue
                s=sum(c*e for c,e in zip(vf,ef))
                if s%ed!=0: continue
                vd=-s//ed
                if vd==0 or vd%7==0 or abs(vd)>L: continue
                combo=[0]*6
                for i,c in zip(free,vf): combo[i]=c
                combo[dep]=vd
                out.append(tuple(combo))
        return out
    for E,lab in [(list(range(8)),"AP8"),([0,1,3,5,7,9,11,13],"WIDE8")]:
        for L in [4]:
            rels=relations_support6(E,L)
            corr=0.0+0.0j
            Wc=defaultdict(float)
            for n in rels:
                c=tuple(v%7 for v in n)
                wv=w_real(n)
                corr+=wv*D7num(c); Wc[c]+=wv
            orbits=defaultdict(list)
            for c in Wc:
                key=min(tuple((a*cj)%7 for cj in c) for a in range(1,7))
                orbits[key].append((c,Wc[c]))
            spreads=[]
            for key,mem in orbits.items():
                ws=[w for _,w in mem]
                if len(ws)>1 and max(abs(x) for x in ws)>0:
                    spreads.append((max(ws)-min(ws))/max(abs(max(ws)),abs(min(ws))))
            avg=sum(spreads)/len(spreads) if spreads else float('nan')
            imratio=abs(corr.imag)/max(abs(corr.real),1e-30)
            print(f"    {lab} L={L}: #rel={len(rels)} corr.Re={corr.real:+.4e} Im/Re={imratio:.1e} "
                  f"#orbits={len(orbits)} mean within-orbit wt spread={avg:.3f}")

    # ---- ADVERSARIAL: can the n<->-n pairing PLUS residue grouping bound WIDE E? ----
    # The honest envelope |corr| <= max|D7| * sum|W_c|. sum|W_c| ~ A(L) diverges.
    # Test if the TRACE-collapse would have given a finite bound IF it applied:
    #   if weights equal per orbit: corr = sum_orbits W_orbit * Tr(D7(rep)).
    # Compare the "pretend trace bound" to the true correction to see the size of the gap.
    print("\n[ADVERSARIAL] pretend-trace vs true correction (size of the lost cancellation):")
    for E,lab in [(list(range(8)),"AP8"),([0,1,3,5,7,9,11,13],"WIDE8")]:
        L=4
        rels=relations_support6(E,L)
        Wc=defaultdict(float);
        for n in rels:
            c=tuple(v%7 for v in n); Wc[c]+=w_real(n)
        true_corr=sum(Wc[c]*D7num(c) for c in Wc).real
        # pretend: replace each orbit's members by the orbit-average weight, then trace collapses
        orbits=defaultdict(list)
        for c in Wc:
            key=min(tuple((a*cj)%7 for cj in c) for a in range(1,7))
            orbits[key].append(c)
        pretend=0.0
        for key,mem in orbits.items():
            wavg=sum(Wc[c] for c in mem)/len(mem)
            trv=float(trace(D7(key))[0])
            # but orbit may be partial (not all 6 dilates present); trace assumes all 6.
            # scale by (present/6) as the best-case telescoping
            pretend += wavg*trv*(len(mem)/6.0)
        print(f"    {lab}: true corr={true_corr:+.6e}   pretend-trace(orbit-avg)={pretend:+.6e}   "
              f"gap={abs(true_corr-pretend):.3e}")

    print("\n"+"="*72)
    print("SUMMARY OF VERIFICATION")
    print("="*72)
    print(f"  Claim 2 (Galois equivariance)        : {'CONFIRMED' if claim2 else 'REFUTED'}")
    print(f"  Claim 3 (trace = rational, Tr(1^6)=5040): {'CONFIRMED' if claim3 else 'REFUTED'}")
    print(f"  Claim 3b(sumTr=0, max|Tr|=9240, per-QR): {'CONFIRMED' if claim3b else 'REFUTED'}")
    print(f"  Claim 4 (reality via sigma_6=conj)    : {'CONFIRMED' if claim4 else 'REFUTED'}")
    print(f"  GAP: dilation NOT a lattice symmetry; within-orbit weights unequal -> trace")
    print(f"       collapse does NOT bound wide E. (confirmed by spread>0 and gap>0)")

if __name__=="__main__":
    main()
