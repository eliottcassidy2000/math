from fractions import Fraction as F
import random
random.seed(20260707)

def distZ(q):  # distance of Fraction q to nearest integer, exact
    r = q - int(q)
    if r < 0: r += 1
    return min(r, 1-r)

def M_at(v, t):  # min_i ||v_i t||
    return min(distZ(F(vi)*t) for vi in v)

UNITS=[1,3,5,9,11,13]  # (Z/14)*

def witness_ge_1_14(v, d, L, search=400):
    """Does some c in (Z/14)* admit delta with all ||v_i (t_c+delta)|| >= 1/14?
    t_c = c/(14 d L). Search delta in a small symmetric range (step ~ 1/(14 d L^2))."""
    base = F(1,14*d*L)
    step = F(1, 14*d*L*L)   # O(1/L^2) scale
    best=F(0); bestc=None
    for c in UNITS:
        tc = c*base
        for m in range(-search, search+1):
            t = tc + m*step
            mm = M_at(v, t)
            if mm >= F(1,14):
                return True, c, mm
            if mm>best: best,bestc=mm,c
    return False, bestc, best

print("=== INDEPENDENT verify + STRESS of klein-S152 conjugate witness (opus-S131) ===")
print("    escape family v_i = a_i + L*d*i (i=1..13); claim: some c in (Z/14)* gives M>=1/14\n")
fails=[]; tested=0
configs=[(1,500,1),(1,3600,3),(3,3600,3),(1,500,6),(2,4000,6),(1,20000,10)]
# a-patterns: random, adversarial (alternating extreme, all-max, binding-targeted)
def apats(A):
    yield [random.randint(-A,A) for _ in range(13)]
    yield [A*(-1)**i for i in range(13)]          # alternating extreme
    yield [A]*13                                   # all max
    yield [(-A if i%3==0 else A) for i in range(13)]
    yield [random.choice([-A,A]) for _ in range(13)]  # extreme random
for (d,L,A) in configs:
    for _ in range(6):
        for a in apats(A):
            v=[a[i]+L*d*(i+1) for i in range(13)]
            if len(set(v))<13: continue
            ok,c,mm=witness_ge_1_14(v,d,L); tested+=1
            if not ok: fails.append((d,L,A,a,float(mm)))
print(f"tested: {tested} escape families (random + adversarial a-patterns)")
print(f"families with NO c giving M>=1/14 (within delta-search): {len(fails)}")
for f in fails[:8]: print("   FAIL:", f[:3], "a=",f[3], "bestM=",f[4])
print(f"\nVERDICT: {'conjugate witness HOLDS on all tested (incl adversarial) => mechanism robust' if not fails else 'FOUND failures -- mechanism may have gaps at these a; INVESTIGATE'}")
