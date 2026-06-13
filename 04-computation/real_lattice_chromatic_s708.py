#!/usr/bin/env python3
"""
real_lattice_chromatic_s708.py  (monad-explorer-2026-06-06-S708)

Two checks to close the dispatched question ("2D-realizable spectrum between the
triangular lattice and the CM field"):

(A) GENERIC (non-arithmetic) real 2D lattices have tiny connection sets, hence
    chi(U(L,D)) <= 2 for all achieved D.  Demonstrated on irrational Gram
    matrices: every achieved norm has r(D)=2 (just +-v) -> disjoint paths -> chi=2.
    RATIONALE: >=3 integer vectors of equal Gram-norm impose >=3 equations on the
    3 Gram entries (a,b,c); generically (irrational, non-proportional) this is
    inconsistent, so r(D)<=2.  Odd cycles (chi>=3) need >=3 equal-norm vectors
    => the lattice must be ARITHMETIC (Gram rational up to scale = integral form).

(B) WIDER integral-form verification at small discriminants with large class
    number (richest popular norms), larger Dmax, to stress-test chi<=3.

Conclusion target: chi(U(L,D)) <= 3 for ALL 2D lattices; =3 max is the triangular
lattice.  No 2D group between triangular (chi=3) and "above" exists -- there is no
"above" in 2D.
"""
import sys, math
sys.path.insert(0,'04-computation')
from all_forms_chromatic_s708 import chi_of

# ---------- (A) generic real lattices ----------
def real_conn_counts(gram, Dlist, tol=1e-9):
    """gram=(a,b,c) reals; count integer vectors with norm in tol of each D."""
    a,b,c=gram
    lam_min=( (a+c) - math.sqrt((a-c)**2+4*b*b) )/2.0
    out={}
    for D in Dlist:
        bound=int(math.sqrt(D/lam_min))+2
        cnt=0
        for x in range(-bound,bound+1):
            for y in range(-bound,bound+1):
                if x==0 and y==0: continue
                q=a*x*x+2*b*x*y+c*y*y
                if abs(q-D)<tol: cnt+=1
        out[D]=cnt
    return out

def test_generic_real():
    print("="*78)
    print("(A) GENERIC (irrational) real 2D lattices: connection-set sizes")
    print("="*78)
    grams=[
        ('diag(1, sqrt2)',      (1.0, 0.0, math.sqrt(2))),
        ('diag(1, pi)',         (1.0, 0.0, math.pi)),
        ('(1, 1/pi, e)',        (1.0, 1.0/math.pi, math.e)),
        ('(1, sqrt2/3, sqrt5)', (1.0, math.sqrt(2)/3, math.sqrt(5))),
        ('(1, 0.31415, 1.7321)',(1.0, 0.31415, 1.7321)),
    ]
    worst=0
    for name,g in grams:
        # achieved norms = actual Q-values of small vectors
        a,b,c=g
        vals=set()
        for x in range(-7,8):
            for y in range(-7,8):
                if x==0 and y==0: continue
                vals.add(round(a*x*x+2*b*x*y+c*y*y,6))
        vals=sorted(vals)[:60]
        counts=real_conn_counts(g, vals)
        mx=max(counts.values())
        worst=max(worst,mx)
        # how many norms have >2 reps?
        rich=[D for D,k in counts.items() if k>2]
        print(f"  {name:22s}: max r(D)={mx}, #norms with r>2: {len(rich)}")
    print(f"\n  => generic real lattices: max connection-set size {worst} "
          f"({'<=2 => chi<=2 (disjoint paths)' if worst<=2 else 'has odd-cycle potential'})")
    return worst

# ---------- (B) wider integral verification ----------
def test_integral_wide():
    print("\n"+"="*78)
    print("(B) WIDER integral-form check: small-disc rich forms, larger Dmax")
    print("="*78)
    # principal + non-principal forms at discriminants with class number >1
    # (richest popular norms): d such that h(-d)>1
    forms=[
        # (A,B,C, label)
        (1,0,5,'Q(sqrt-5) principal'), (2,2,3,'disc-20 nonprincipal'),
        (1,0,6,'Q(sqrt-6)'),           (2,0,3,'disc-24 nonprincipal'),
        (1,1,6,'Q(sqrt-23) princ'),    (2,1,3,'disc-23 form2'), (2,-1,3,'disc-23 form3'),
        (1,1,16,'disc-63 princ'),      (2,1,8,'disc-63 f2'), (4,1,4,'disc-63 f3'),
        (1,0,14,'Q(sqrt-14)'),         (2,0,7,'disc-56 f2'), (3,2,5,'disc-56 f3'),(3,-2,5,'disc-56 f4'),
        (1,1,9,'disc-35 princ'),       (3,1,3,'disc-35 f2'),
    ]
    Dmax=300
    gmax=2; bad=[]
    for (A,B,C,label) in forms:
        disc=B*B-4*A*C
        chiset=set()
        for D in range(1,Dmax+1):
            chi,_=chi_of(A,B,C,D)
            if chi is None: continue
            chiset.add(chi)
            if chi>=4: bad.append((A,B,C,D))
        gmax=max(gmax,max(chiset) if chiset else 2)
        print(f"  ({A:2d},{B:2d},{C:2d}) disc={disc:5d} {label:22s}: chi-set={sorted(chiset)}")
    print(f"\n  => wider integral max chi: {gmax}; chi>=4 cases: {len(bad)}")
    return gmax, bad

if __name__=='__main__':
    w=test_generic_real()
    gmax,bad=test_integral_wide()
    print("\n"+"="*78)
    print("SUMMARY")
    print(f"  generic real lattice max connection size: {w}")
    print(f"  arithmetic (integral) lattice max chi: {gmax}")
    print("grep RESULT:", "ALL_LE_3_NO_CHI4" if (gmax<=3 and not bad and w<=2) else "EXCEPTION")
