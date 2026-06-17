#!/usr/bin/env python3
"""
ANGLE E -- Cross-modulus REGIONS/SECTIONS unify with the covering reduction (THM-523).

The "regions/sections" reframe of LRC(14): at grid time tau=a/q, the circle splits into q
equal SECTIONS; runner i sits in section r_i = v_i*a mod q. The observer (section 0) is
LONELY at a/q  <=>  no r_i == 0  <=>  no v_i == 0 (mod q)  (a "q-witness", q<=14).

So "regions at modulus q" IS exactly the THM-523 q-witness lemma.
The COVERING condition (S contains a multiple of EVERY q in 2..14) is precisely
"no single grid-modulus q<=14 has a clear observer section" -- every coarse grid is blocked.

Three sub-questions:
  (1) For covering sets, confirm every modulus q in 2..14 has a runner in section 0
      (i.e. the q-grid observer is always hit). The lonely time must live off ALL these grids.
  (2) The actual lonely tau* has denominator = some v_a +- v_b (a FINE section count). Relate
      this fine modulus to the speeds: it is the # of fine sections in which a BINDING PAIR
      becomes equidistant from the observer.
  (3) MIXED-RADIX / CRT: can a lonely point be assembled by choosing residues mod different
      q simultaneously? When the moduli are pairwise coprime, CRT gives a single residue mod
      their product -> a grid point of that fine modulus. We test whether this constructs a
      genuine lonely point (M >= 1/14 witness) for sub-families.

Everything EXACT (Fraction). stdlib only.
"""
from fractions import Fraction as F
from math import gcd, comb
from functools import reduce
from itertools import combinations

# ---------- exact M tool (verbatim from prompt) ----------
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r; return r if r<=F(1,2) else 1-r
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def M(S):
    b=F(0); at=None
    for t in cand(S):
        v=g(S,t)
        if v>b: b=v; at=t
    return b,at

def lcm(a,b): return a*b//gcd(a,b)
def lcml(xs): return reduce(lcm, xs, 1)

# ---------- region/section helpers ----------
def sections(S,q,a):
    """section index of each runner on the q-grid at time a/q."""
    return [ (v*a) % q for v in S ]

def q_witness_clear(S,q):
    """Is there a unit a (gcd(a,q)=1) at which section 0 is clear? Equivalently no v==0 mod q.
       (THM-523: if S omits a multiple of q then tau=1/q is lonely with gap>=1/q.)
       Returns list of units a in (Z/q)* for which the observer section is clear."""
    units=[a for a in range(1,q) if gcd(a,q)==1]
    clear=[a for a in units if all((v*a)%q!=0 for v in S)]
    return clear

def is_covering(S, qs=range(2,15)):
    """S is covering iff it contains a multiple of every q in qs (2..14)."""
    miss=[q for q in qs if not any(v%q==0 for v in S)]
    return (len(miss)==0), miss

# ============================================================================
print("="*78)
print("ANGLE E -- CROSS-MODULUS REGIONS == THE q-WITNESS (THM-523 reframed)")
print("="*78)
print("""
IDENTITY: grid tau=a/q -> q sections; observer lonely at a/q  <=>  no v_i==0 (mod q).
  -> "regions at modulus q" IS the q-witness.
  -> COVERING (mult of every q in 2..14)  <=>  every coarse grid q<=14 has section 0 BLOCKED.
""")

# Test families ----------------------------------------------------------------
tight   = list(range(1,14))                       # perfect SDR, M=1/14 (boundary, LRC holds tight)
hc84    = [1,2,3,4,5,6,7,8,9,10,11,13,84]         # covering hard core, 84=6*14
drop6u98= [1,2,3,4,5,7,8,9,10,11,12,13,98]        # interior drop j=6, u=98=7*14
# another covering hard core: must hit every q in 2..14. {1..13} already hits 2..13; need mult of 14.
hc70    = [1,2,3,4,5,6,7,8,9,10,11,13,70]         # 70=5*14, also covering (drop 12, add 70 keeps 12? no 12 dropped)
# fix: keep 12, drop nothing of 2..13 coverage; {1..11,13}=12 elts misses only 12 from {1..13}; is 12 covered? 12|? need mult of 12. none. so add a mult of lcm? simplest: {1..13} minus one + a big 14-multiple
hc_a    = [1,2,3,4,5,6,7,8,9,10,11,12,14]         # {1..12,14}: covers 2..12 and 14 and 7 and... 13? NO mult of 13 -> not covering
# Build a clean covering 13-set: take {1..11,13} (covers 2..11,13) then need 12 and 14 -> add 84=lcm(12,14)? 84=12*7=14*6, single number covers both 12&14. -> {1..11,13,84} = hc84. good (already have).

FAMS = [
    ("tight AP {1..13}", tight),
    ("covering hardcore 84 {1..11,13,84}", hc84),
    ("interior-drop6 u98 {1..5,7..13,98}", drop6u98),
]

print("-"*78)
print("(1) PER-MODULUS OBSERVER-SECTION TABLE: is section 0 clear at some unit a?")
print("-"*78)
for name,S in FAMS:
    cov, miss = is_covering(S)
    print(f"\n[{name}]  covering={cov}" + (f"  (misses q={miss})" if miss else ""))
    print(f"  {'q':>3} | clear-units a (gcd(a,q)=1, observer section 0 empty) | witness?")
    blocked_all=True
    for q in range(2,15):
        clear=q_witness_clear(S,q)
        wit = "YES tau=1/q lonely" if clear else "NO (section0 always hit)"
        if clear: blocked_all=False
        # show whether S has a multiple of q
        hasmult = any(v%q==0 for v in S)
        print(f"  {q:>3} | {str(clear):<48} | {wit}   [mult of q in S: {hasmult}]")
    print(f"  => ALL coarse grids q<=14 BLOCK the observer? {blocked_all}")
    if blocked_all:
        Mv, t = M(S)
        print(f"     (so any lonely time must be OFF-GRID; exact M={Mv}={float(Mv):.5f}, tau*={t})")

# ============================================================================
print("\n"+"="*78)
print("(2) THE FINE MODULUS: lonely tau* denominator = v_a +- v_b (a BINDING PAIR)")
print("="*78)
print("""
M(S) is attained at tau* where TWO runners are equidistant from the observer.
At a crossing tau=k/(v_a+-v_b), runners a,b both sit at distance ||k v_a/(v_a+-v_b)||.
The denominator d = v_a+-v_b is the # of FINE sections; the binding pair carves the gap.
""")
def binding_pair(S, tstar):
    """find pair (and sign) whose v_a+-v_b equals the denom of tstar, and that are the
       two closest runners to observer (the equidistant binders)."""
    den = tstar.denominator
    Ssort=sorted(set(S))
    # the two runners realizing the min distance:
    dists=sorted(((nrm(v*tstar),v) for v in Ssort))
    binders=[v for d,v in dists if d==dists[0][0]]
    pairs=[]
    for i in range(len(Ssort)):
        for j in range(i+1,len(Ssort)):
            a,b=Ssort[i],Ssort[j]
            for s,d in (("+",a+b),("-",b-a)):
                if d==den:
                    pairs.append((a,b,s,d))
    return den, binders, pairs

for name,S in FAMS:
    Mv,t=M(S)
    den,binders,pairs=binding_pair(S,t)
    print(f"\n[{name}]  M={Mv}={float(Mv):.5f}  tau*={t}")
    print(f"  fine modulus (denominator of tau*) = {den}")
    print(f"  closest runners (binders, equidistant from observer) = {binders}  at dist {nrm(binders[0]*t)}")
    print(f"  speed pairs (a,b,sign) with v_a+-v_b == {den}: {pairs if pairs else 'NONE matching exactly'}")
    # also report the gap width vs 1/14
    print(f"  gap M={Mv}  vs threshold 1/14={F(1,14)}:  {'LONELY (>1/14, LRC FAILS here)' if Mv>F(1,14) else ('TIGHT (=1/14)' if Mv==F(1,14) else 'safe (<1/14)')}")

# ============================================================================
print("\n"+"="*78)
print("(3) MIXED-RADIX / CRT SECTION SCHEME: can residues mod different q build a lonely pt?")
print("="*78)
print("""
The covering-system DODGE: choose, for several coprime moduli q1,q2,..., a target residue
class for tau (i.e. pick which FINE section the observer's gap lands in). CRT glues a residue
mod Q=prod(qi). The question: does there exist tau with ||v tau||>=1/14 for all v, found by
specifying residues mod a few small moduli simultaneously?

We make this CONCRETE two ways:
  (A) Per-modulus "clear-section" sets. For modulus q, the set of grid points a/q that keep
      the OBSERVER clear is U_q={a: no v==0 mod q at a} (the q-witness units). For a SINGLE
      coarse q this only guarantees gap>=1/q at the lattice point a/q. CRT across coprime q's
      does NOT improve this -- a/q is the SAME point regardless of how many residues we name,
      because all coarse grids share the rational a/q. So CRT over COARSE moduli q<=14 gives
      no new lonely point beyond what one q already gives. (We verify the gap is exactly the
      single-q witness, not boosted.)
  (B) The REAL lonely points are at FINE moduli d=v_a+-v_b (>14). Test: pick the binding-pair
      fine modulus d and a residue k with gcd(k,d) handled; does k/d beat 1/14? This is just
      evaluating M at the fine grid -- CRT is the bookkeeping of WHICH fine section.
""")

def clear_sections_fine(S, d):
    """On the FINE d-grid, which residues k/d (0<k<d) have observer gap >= 1/14?
       returns sorted list of (k, gap) with gap>=1/14."""
    out=[]
    for k in range(1,d):
        t=F(k,d)
        gp=g(S,t)
        if gp>=F(1,14):
            out.append((k,gp))
    return out

print("\n(3A) CRT over coarse moduli adds nothing -- the witness gap is single-q:")
for name,S in FAMS:
    cov,_=is_covering(S)
    # find coprime coarse moduli that S omits (clear for at least one); for covering sets none exist
    omitted=[q for q in range(2,15) if not any(v%q==0 for v in S)]
    print(f"  [{name}] omitted coarse moduli (q with clear observer) = {omitted}", end="")
    if not omitted:
        print("  -> NONE (covering): no coarse-grid lonely point, CRT has no coarse ingredients.")
    else:
        # take any one, show gap = 1/q witness
        q=omitted[0]; a=q_witness_clear(S,q)[0]
        print(f"  -> e.g. q={q}, a={a}: gap at {a}/{q} = {g(S,F(a,q))} (= single-q witness {F(1,q)})")

print("\n(3B) FINE-grid lonely sections (the genuine lonely points), by binding-pair modulus:")
for name,S in FAMS:
    Mv,t=M(S)
    den=t.denominator
    if Mv<=F(1,14):
        print(f"  [{name}] M={Mv}<=1/14 -> NO lonely point on any grid (LRC holds). fine modulus {den} has no clearing section >1/14.")
        continue
    cs=clear_sections_fine(S,den)
    print(f"  [{name}] fine modulus d={den}: lonely sections k/d with gap>=1/14:")
    for k,gp in cs:
        # CRT description: residue of k modulo factors of d
        fac=[]
        dd=den; p=2
        primes=[]
        while p*p<=dd:
            if dd%p==0:
                primes.append(p)
                while dd%p==0: dd//=p
            p+=1
        if dd>1: primes.append(dd)
        crt_desc={pp:(k%pp) for pp in primes}
        print(f"     k={k}: tau={F(k,den)}  gap={gp}={float(gp):.5f}  | CRT residues k mod prime-factors {crt_desc}")
    print(f"     => these ARE the constructive lonely points; CRT just names which fine section k lands in.")

# ============================================================================
print("\n"+"="*78)
print("(3C) DOES MIXED-RADIX BUILD A LONELY POINT FROM SCRATCH? -- honest test")
print("="*78)
print("""
Claim to test: 'pick residues mod several small q simultaneously to dodge the covering system'.
Concretely: choose a target tau by fixing tau mod (1/q) sections for q=2,3,5,7 (pairwise coprime,
product 210). For EACH q, a 'safe' section means the q-grid points nearest tau keep all runners off
the observer band of width 1/14. Then ask: does some CRT-glued tau in (0,1/2] beat 1/14?

This is equivalent to scanning the d=210 grid (and its refinements) -- which the M tool already
does exactly via its candidate set. So the honest statement: CRT/mixed-radix is a SEARCH
RE-PARAMETERIZATION, not a new construction. It locates lonely points iff M(S)>1/14, which holds
exactly for the covering hard cores (off-grid) and fails for tight/SDR cases.
""")
# demonstrate: scan Q=210 grid for each family, compare best gap to true M
for name,S in FAMS:
    Mv,t=M(S)
    Q=210
    best=F(0); bk=None
    for k in range(1,Q):
        tt=F(k,Q)
        if tt>F(1,2): continue
        gp=g(S,tt)
        if gp>best: best=gp; bk=k
    print(f"  [{name}] best gap on coarse CRT grid Q=210: {best}={float(best):.5f} at {bk}/210 ; true M={Mv}={float(Mv):.5f}")
    print(f"     -> coarse-CRT grid {'REACHES' if best==Mv else 'MISSES'} the true optimum (true optimum lives at fine modulus {t.denominator}).")

# ============================================================================
print("\n"+"="*78)
print("(4) SWEEP: across MANY covering hard cores, is the fine modulus always v_a+-v_b? "
      "and is it ALWAYS off every coarse grid q<=14?")
print("="*78)
def gen_covering_cores():
    """13-sets = {1..13} with one element j (2..13) replaced by a multiple of 14 that restores
       coverage. Replacing j by 14*m: still covering iff coverage of j's role is kept.
       {1..13}\{j} covers every q in 2..14 EXCEPT possibly the q's that ONLY j provided.
       j provides multiples-of-q for q | j (and q where j is the only multiple <=13).
       Adding 14*m always covers 2,7,14. Simplest robust family: drop j in {a multiple of 14 fixes}.
       We just generate candidates and filter is_covering."""
    out=[]
    base=set(range(1,14))
    for j in range(2,14):
        for m in range(1,13):  # 14m up to 168
            w=14*m
            if w in base: continue
            S=sorted((base-{j})|{w})
            if len(S)!=13: continue
            cov,_=is_covering(S)
            if cov: out.append((j,w,S))
    return out

cores=gen_covering_cores()
print(f"  generated {len(cores)} covering hard cores of form ({{1..13}}\\{{j}})∪{{14m}}")
fails=0; total=0; modulus_is_pair=0; offgrid_ok=0
worst=(F(2),None); best_loose=(F(0),None)
examples=[]
for j,w,S in cores:
    total+=1
    Mv,t=M(S)
    den=t.denominator
    # is fine modulus a v_a+-v_b?
    Ss=sorted(S); ispair=False
    for i in range(len(Ss)):
        for jx in range(i+1,len(Ss)):
            if Ss[i]+Ss[jx]==den or Ss[jx]-Ss[i]==den: ispair=True
    if ispair: modulus_is_pair+=1
    # off every coarse grid?
    blocked=all(not q_witness_clear(S,q) for q in range(2,15))
    if blocked: offgrid_ok+=1
    if Mv<F(1,14): fails+=1   # would be an LRC counterexample (none expected)
    if Mv<worst[0]: worst=(Mv,(j,w,t,den))
    if Mv>best_loose[0]: best_loose=(Mv,(j,w,t,den))
    if len(examples)<6: examples.append((j,w,Mv,t,den,ispair,blocked))
print(f"  LRC counterexamples (M<1/14): {fails}/{total}")
print(f"  fine modulus equals some v_a+-v_b: {modulus_is_pair}/{total}")
print(f"  every coarse grid q<=14 blocked (covering -> off-grid only): {offgrid_ok}/{total}")
print(f"  tightest (closest to 1/14) covering core: M={worst[0]}={float(worst[0]):.5f}  (j,w,tau*,d)={worst[1]}")
print(f"  loosest covering core: M={best_loose[0]}={float(best_loose[0]):.5f}  (j,w,tau*,d)={best_loose[1]}")
print(f"  sample (j, w=14m, M, tau*, fineMod d, d=v_a+-v_b?, allCoarseBlocked?):")
for e in examples:
    print(f"    j={e[0]:>2} w={e[1]:>3}  M={str(e[2]):>7}={float(e[2]):.5f}  tau*={str(e[3]):>7}  d={e[4]:>4}  pair={e[5]}  blocked={e[6]}")

print("\n"+"="*78)
print("DONE.")
print("="*78)
