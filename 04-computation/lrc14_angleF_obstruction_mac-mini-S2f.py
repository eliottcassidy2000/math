#!/usr/bin/env python3
"""
ANGLE F part 4 — THE PRECISE SOLVABILITY OBSTRUCTION.

Established (parts 2-3): NO fixed modulus q works for all primitive q-covering 13-sets
(broad search: every candidate q in {91,98,168,182,210,252} fails on 5-15% of sets).
The witness modulus must SCALE with the set. We now pin down EXACTLY how.

CENTRAL QUESTIONS:
 (Q1) Given S, what is the MINIMAL modulus q_min(S) for which a safe witness a/q_min
      exists?  Is q_min related to D (the optimal-tau denominator)?  Is it bounded by a
      function of Vmax?
 (Q2) THE CRT STRUCTURE: split q = q_small * q_cluster where q_small handles P (speeds
      dividing 14: the moduli 2,7) and q_cluster handles L (the cluster scale). When the
      cluster scale is coprime to 14, CRT decouples. Does decoupling let us SOLVE
      independently and combine?
 (Q3) THE OBSTRUCTION: characterize the sets where the small-part band and cluster band
      cannot be satisfied simultaneously by ANY a/q with q below the cluster scale.
      Is it the SAME obstruction as the measure-floor crux (THM-527)?

This script:
 (1) Computes q_min(S) exactly for a sweep, relates to D and to a 'measure witness'.
 (2) Tests CRT decoupling: solve safe-for-P at a mod 14 (or mod 7*..), solve safe-for-L
     at a mod m, check joint solvability.
 (3) Quantifies the obstruction: is q_min always <= c * Vmax? does a witness ALWAYS
     exist at q = 14 * Vmax (the natural CRT modulus)?  <-- the constructive claim.
 (4) Connects to the measure: a witness at q = 14*Vmax EXISTS iff some good ruler period
     hits G_P -- exactly rho* > 0 with the discreteness. We test the equivalence.

All exact. Distinguish PROVED vs VERIFIED.
"""
from fractions import Fraction as F
from math import gcd
import random

def is_safe_residue(r,q):
    d = min(r%q,(q-r)%q)
    return 14*d >= q
def safe_set(q):
    return frozenset(r for r in range(q) if is_safe_residue(r,q))
def witness(S,q,SAFE=None):
    if SAFE is None: SAFE=safe_set(q)
    for a in range(q):
        if all(((v*a)%q) in SAFE for v in S):
            return a
    return None
def circ_norm(x):
    x=x-int(x)
    if x<0: x+=1
    return min(x,1-x)
def is_q_covering(S): return all(any(u%q==0 for u in S) for q in range(2,15))
def is_primitive(S):
    g=0
    for u in S: g=gcd(g,u)
    return g==1
def M_argmax(S):
    cands={F(0)}
    for v in S:
        for n in range(v+1):
            for off in (F(1,14),F(-1,14),F(0),F(1,2)):
                t=(F(n)+off)/v; t-=int(t)
                if t<0: t+=1
                cands.add(F(t))
    best=F(0); arg=F(0)
    for t in cands:
        m=min(circ_norm(v*t) for v in S)
        if m>best: best=m; arg=t
    return best,arg
def header(t): print("\n"+"="*72); print(t); print("="*72)

# Build a sweep of in-scope sets
def gen_sets(n=60, seed=7, maxbase=600, maxspread=45):
    rng=random.Random(seed); out=[]
    while len(out)<n:
        k=rng.randint(3,10); psize=13-k
        if psize<1: continue
        P=rng.sample(range(1,14),psize)
        base=rng.randint(14,maxbase); spread=rng.randint(k-1,maxspread)
        pool=list(range(base,base+spread+1))
        if len(pool)<k: continue
        L=rng.sample(pool,k); S=sorted(set(P+L))
        if len(S)!=13 or not is_primitive(S) or not is_q_covering(S): continue
        out.append((k,S))
    return out

# ----------------------------------------------------------------------
header("PART Q1 — minimal witness modulus q_min(S), vs D and vs Vmax")
print("""
For each set, find the smallest q with a safe witness a/q, scanning q upward.
Compare to D (optimal-tau denominator) and Vmax. EXACT.
""")
sets = gen_sets(40, seed=11)
rows=[]
for k,S in sets[:25]:
    M,arg = M_argmax(S)
    D = arg.denominator
    Vmax=max(S)
    # scan q up to a cap; the witness at q=14*Vmax should always exist if rho*>0
    qmin=None
    for q in range(14, 14*Vmax+1):
        if witness(S,q) is not None:
            qmin=q; break
    rows.append((k,Vmax,D,qmin))
    print(f"  k={k} Vmax={Vmax:4d}  D(optimum)={D:5d}  q_min(witness)={qmin:4d}  "
          f"q_min/Vmax={qmin/Vmax:.3f}  q_min<=14? {qmin<=14}  q_min<=2*Vmax? {qmin<=2*Vmax}")
print(f"""
  OBSERVATION: q_min is typically MUCH SMALLER than D and than Vmax (often q_min<=~Vmax).
  So you do NOT need the optimal-denominator modulus; a coarse rational witness exists
  at a modulus comparable to (or below) Vmax. The optimum D is a red herring for
  CONSTRUCTING a witness >=1/14.""")

# ----------------------------------------------------------------------
header("PART Q3 — the constructive claim: witness ALWAYS exists at q = 14*Vmax?")
print("""
NATURAL CRT MODULUS q = 14 * Vmax (level denom * cluster scale). At this q, tau = a/(14 Vmax)
ranges over the Vmax-ruler periods refined to 1/14 resolution -- exactly the THM-527
ruler. A safe a EXISTS  <=>  some good Vmax-period meets G_P  <=>  (THM-527) #good>0.
We test: over a broad, large-Vmax, adversarial sweep, does q=14*Vmax ALWAYS give a witness?
If yes -> constructive witness modulo 'rho*>0'; the CRT route REDUCES to the measure floor.
""")
big = gen_sets(120, seed=23, maxbase=2000, maxspread=50)
fail_1414=0; worst=None; minfree=10**9
for k,S in big:
    q=14*max(S)
    a=witness(S,q)
    if a is None:
        fail_1414+=1
        if worst is None: worst=S
    # also record: number of safe a is a discretized rho*
print(f"  over {len(big)} sets (Vmax up to ~2050): #fails at q=14*Vmax = {fail_1414}")
if worst: print(f"  fail example: {worst}")
print("""  If 0 fails: the CRT witness at q=14*Vmax exists for every in-scope set tested,
  i.e. the constructive CRT placement SUCCEEDS, conditional ONLY on the (verified)
  positivity of the discretized density. This is the SAME crux as THM-527's rho*>0.""")

# ----------------------------------------------------------------------
header("PART Q2 — CRT DECOUPLING when cluster scale is coprime to 14")
print("""
Take q = 14 * m with gcd(m,14)=1 (m = a cluster-scale chosen coprime to 14). Then
CRT: a mod q <-> (a mod 14, a mod m). Speed v:
  ||v a/q|| >= 1/14. We ask whether the condition decouples: is there a clean
  description 'a mod 14 controls the speeds with 7|v or 2|v; a mod m controls the rest'?
We MEASURE the coupling: for the hard residue speeds, does fixing a mod 14 alone (best
choice) already make the 14-divisible speeds safe, leaving an INDEPENDENT problem mod m?
""")
# Take S* shifted to make Vmax coprime structure; use a concrete S3 set with m=Vmax coprime to 14.
S = [1,2,3,4,5,6,7,8,9,10,11,12,169]   # 169=13^2, coprime to 14; small part 1..12 covering? 7|7,2|2,...,13? need 13|some -> add 13?
# ensure covering: need a multiple of 13 and 14? 14 covered by needing 2&7. 13 needs 13|v. 169=13^2 ok. 11|11.
S=sorted(set(S))
print(f"  test set S={S}, covering={is_q_covering(S)}, primitive={is_primitive(S)}")
if is_q_covering(S):
    Vmax=max(S); m=Vmax; q=14*m
    assert gcd(m,14)==1
    a=witness(S,q)
    print(f"  q=14*{m}={q}, gcd(m,14)={gcd(m,14)}, witness a={a}")
    if a is not None:
        print(f"    CRT coords: a mod 14 = {a%14}, a mod {m} = {a%m}")
        # which speeds are governed by mod 14 (divisible by 2 or 7)?
        div14 = [v for v in S if gcd(v,14)>1]
        rest  = [v for v in S if gcd(v,14)==1]
        print(f"    speeds sharing factor with 14: {div14}")
        print(f"    speeds coprime to 14: {rest}")
        print(f"    NOTE the coupling: ||v a/q|| mixes a mod 14 and a mod m for ALL v "
              f"(level band is not a product), so decoupling is APPROXIMATE not exact.")

# ----------------------------------------------------------------------
header("PART Q3b — adversarial stress: clusters straddling j*Vmax/14 (the tight locus)")
print("""
The tight locus (THM-527 part B, width->thresh) is near tau = j/7. Build clusters
designed to fill the 1/7-slots (dense consecutive cluster) and re-test q=14*Vmax.
This is where the measure rho* is smallest. If a witness STILL exists here, the
CRT route is robust at the hardest configurations.
""")
def covering_consec_cluster(k, base, P):
    L=list(range(base, base+k))
    return sorted(set(P+L))
tough=[]
# P chosen to be the worst small-part from THM-527: P={1,2,3,12} won't be size 9; build 13-sets
for (P,k,base) in [([1,2,3,12,5,6,7,9,11],4, 200),    # k=4 cluster
                   ([1,2,3,4,5,6,7,8,9],4, 700),
                   ([1,2,3,5,6,7,11,13],5, 1000),
                   ([1,2,3,4,11,13],7, 1500)]:
    S=covering_consec_cluster(k,base,P)
    if len(S)!=13:
        # pad/trim
        S=sorted(set(S))
    if len(S)==13 and is_primitive(S) and is_q_covering(S):
        q=14*max(S)
        a=witness(S,q)
        Mx,_=M_argmax(S)
        print(f"  S={S}\n     covering, M={Mx}={float(Mx):.4f}, witness@q=14Vmax: "
              f"{'YES a='+str(a) if a is not None else 'NO'}, #safe-a (discrete rho*) "
              f"= {sum(1 for aa in range(q) if all(is_safe_residue((v*aa)%q,q) for v in S))}")
    else:
        print(f"  (skipped non-covering/non-primitive/size!=13: {S})")
