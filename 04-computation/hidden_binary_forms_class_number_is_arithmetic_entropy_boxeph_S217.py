#!/usr/bin/env python3
"""hidden_binary_forms_class_number_is_arithmetic_entropy_boxeph_S217.py -- boxeph-2026-07-21-S217

HISTORICAL / REFUTED SYNTHESIS (MISTAKE-230, MISTAKE-229).  The reduced-form
tables remain useful finite checks.  The rational-prime ``log2(h)`` entropy
claim is invalid because inverse proper classes represent the same integers,
and the transfer to the THM-2053 LRC gate is invalid because that gate has no
discriminant-minus-seven form.  This script intentionally retains the original
calculation so the failure is reproducible; its interpretive conclusions are
not current mathematics.

Hidden binary forms + information theory. The binary form of disc -p (S215/S216 anisotropic gate) does not
live alone: its DISCRIMINANT has a whole CLASS GROUP of forms (Gauss composition). The hidden forms are the
NON-PRINCIPAL classes. Information-theoretic reading:

  CLASS NUMBER h(D) = the ARITHMETIC ENTROPY of the form = the bits BEYOND local (Legendre) data needed to
  decide which form represents a prime p. For h=1 (Heegner: -3,-7,-11 = Paley primes 3,7,11) representation
  is decided by ONE bit -- the Legendre symbol (D/p) -- so ZERO hidden entropy: local determines global,
  the residual is RIGID (S216). For h>1 the primes with (D/p)=1 SPLIT among the h forms; which one costs
  log2(h) extra bits (the class group Cl(D) = Gal(Hilbert class field), invisible to any local test).

LRC(14)=2*7 -> disc -7 -> h=1 -> ZERO arithmetic entropy: the S216 anisotropic gate's residual is fully
determined by local S215 Legendre data. That is WHY 7 is rigid; the difficulty is 'no hidden bits to hide in'.
"""
from math import gcd, log2

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)

sep("STATUS: HISTORICAL / REFUTED SYNTHESIS -- READ MISTAKE-230 AND MISTAKE-229")
print("  Retained: reduced-form tables and finite representation checks.")
print("  Retracted: rational primes split among h oriented classes; log2(h) is their entropy; any LRC transfer.")

def legendre(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1
def is_reduced(a,b,c): return (-a< b<=a<=c) and not (a==c and b<0)
def reduced_forms(D):   # D<0: list primitive reduced (a,b,c), b^2-4ac=D  (the class group elements)
    forms=[]; a=1
    while a*a <= -D//3 + 1:
        for b in range(-a, a+1):
            if (b*b-D)%(4*a)==0:
                c=(b*b-D)//(4*a)
                if a<=c and gcd(gcd(a,abs(b)),c)==1 and is_reduced(a,b,c):
                    forms.append((a,b,c))
        a+=1
    return forms
def represents(form,p,B=60):
    a,b,c=form
    for x in range(-B,B+1):
        for y in range(-B,B+1):
            if a*x*x+b*x*y+c*y*y==p: return True
    return False

# ==========================================================================
sep("A  the CLASS GROUP = hidden binary forms; class number h(D) = arithmetic ENTROPY log2 h")
for D in (-3,-4,-7,-8,-11,-15,-23,-31,-47):
    fs=reduced_forms(D); h=len(fs)
    ent=log2(h)
    tag=" HEEGNER h=1: ZERO entropy (rigid)" if h==1 else f" {h-1} HIDDEN non-principal form(s), entropy log2({h})={ent:.2f} bits"
    print(f"  disc {D:4d}: h={h}  forms={fs} ;{tag}")
print("  => 3,7,11 (Paley/Heegner, disc -3,-7,-11) have ONLY the principal form: zero arithmetic entropy.")
print("     -23 (h=3), -31 (h=3), -47 (h=5): hidden non-principal forms = the class-group information.")

# ==========================================================================
sep("B  representation-LOCALITY: h=1 decided by Legendre (1 bit); h>1 splits (D/p)=1 primes among the forms")
def form_rep_table(D, pmax=120):
    fs=reduced_forms(D)
    primes=[p for p in range(2,pmax) if all(p%k for k in range(2,int(p**0.5)+1))]
    rows=[]
    for p in primes:
        if p== -D or D%p==0: continue
        leg=legendre(D,p)   # (D/p); =1 => p represented by SOME form of disc D
        reps=[i for i,f in enumerate(fs) if represents(f,p)]
        rows.append((p,leg,reps))
    return fs,rows
# disc -7 (h=1): represented <=> (D/p)=1, always by the unique form
fs7,rows7=form_rep_table(-7)
ok7=all((len(r[2])>0) == (r[1]==1) for r in rows7)
print(f"  disc -7 (h=1): 'represented by the unique form (1,1,2)'  <=>  (D/p)=1 (ONE Legendre bit)? {ok7}")
print("     sample:", [(p,'(D/p)=%+d'%leg,'rep' if reps else '---') for p,leg,reps in rows7[:8]])
# disc -23 (h=3): (D/p)=1 primes SPLIT among the 3 forms
fs23,rows23=form_rep_table(-23)
split=[0,0,0]
splitprimes={0:[],1:[],2:[]}
for p,leg,reps in rows23:
    if leg==1 and reps:
        # which class (principal index 0 = (1,1,6)); group by first representing form
        cls=reps[0]; split[cls]+=1; splitprimes[cls].append(p)
allleg1_rep = all((r[1]==1)==(len(r[2])>0) for r in rows23)
print(f"  disc -23 (h=3): forms {fs23}; (D/p)=1 <=> represented? {allleg1_rep}")
print(f"     but (D/p)=1 primes SPLIT among the {len(fs23)} classes: counts {split}")
print(f"     principal (1,1,6) primes: {splitprimes[0][:8]}")
print(f"     non-principal (2,+-1,3) primes: {sorted(splitprimes[1]+splitprimes[2])[:8]}")
print("  => for h=1 the Legendre bit DETERMINES representation; for h>1, WHICH form costs log2(h) extra bits")
print("     (the class group = Gal(Hilbert class field), invisible to any local/Legendre test). = arithmetic entropy.")

# ==========================================================================
sep("C  the info-theoretic principle, and LRC(14)=2*7")
print("""  CLASS NUMBER = ARITHMETIC ENTROPY (bits beyond local Legendre to decide representation):
    h=1 (Heegner -3,-7,-11 = Paley 3,7,11): 0 bits -- local (S215 Legendre) DETERMINES global. RIGID.
    h>1 (-23:log2 3=1.58, -47:log2 5=2.32 bits): hidden non-principal forms = class-group entropy.
  LRC(14)=2*7 -> disc -7 -> h=1 -> ZERO arithmetic entropy: codex's THM-2053 anisotropic gate residual
  (S216) is fully pinned by local Legendre data (S215) -- there are NO hidden bits for a counterexample to
  hide in. The difficulty of 7 is exactly its rigidity: a class-number-1 gate has no slack.""")

# ==========================================================================
sep("D  bonus hidden entropy: the tournament SCORE distribution -- transitive (max) vs Paley (0)")
def score_entropy(scores):
    from collections import Counter
    n=len(scores); cnt=Counter(scores)
    return -sum((c/n)*log2(c/n) for c in cnt.values())
for n in (7,11,13):
    transitive=list(range(n))                       # scores 0..n-1, all distinct
    paley=[(n-1)//2]*n                               # Paley/regular: all equal (n-1)/2
    print(f"  n={n}: transitive score-dist entropy = {score_entropy(transitive):.3f} bits (=log2 {n}, MAX spread) ;"
          f"  Paley/regular = {score_entropy(paley):.3f} bits (delta, MIN)")
print("  => the reify-ladder poles are entropy extremes: transitive/AP = max score-entropy (the nullcone")
print("     vertex, the rank-11 gate); Paley = zero score-entropy (the symmetric pole, disc -p, i*sqrt p).")

sep("SUMMARY")
print("""  Hidden binary forms: the disc -p gate (S216) sits in a CLASS GROUP of forms (Gauss composition); the
  non-principal classes are the hidden forms. INFORMATION THEORY: class number h(D) = the ARITHMETIC ENTROPY
  = bits beyond local Legendre needed to decide representation. Heegner (h=1: -3,-7,-11 = Paley 3,7,11) = 0
  bits = local determines global = the S216 residual is RIGID. LRC(14)=2*7 -> -7 -> h=1 -> no hidden bits.
  Bonus: the tournament score distribution's entropy separates the reify-ladder poles (transitive max /
  Paley zero). Both readings: the anisotropic gate is hard BECAUSE it has zero hidden entropy -- nowhere for
  a counterexample to hide, and no class-group slack to exploit; only the local Legendre/Paley data (S215).""")
