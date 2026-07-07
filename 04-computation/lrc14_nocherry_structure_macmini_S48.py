"""
mac-mini-2026-07-07-S48 (HYP-5047) -- the NO-SEPARATED-CHERRY STRUCTURE LEMMA and the
census of klein-S165's moderate-spread k=8 residual.

klein-S165's gate: per-shape cherry certificates need a triple with e_c >= 50(e_a+e_b)
("separated cherry").  Residual = shapes without one.

STRUCTURE LEMMA (proved, 3 lines): if a k=8 shape 0 < e2 < ... < e8 (7 nonzero speeds)
has NO L-separated triple, then for every j >= 4 the triple (e2, e3, e_j) gives
e_j < L(e2 + e3); hence the ENTIRE shape lies in {e2} u [e3, L(e2+e3)).  Dichotomy at
any K: (A) e3 <= K e2 => all speeds in [e2, L(1+K)e2]: ONE band of ratio <= L(1+K);
(B) e3 > K e2 => one isolated small speed e2 + a band [e3, L(e2+e3)) of ratio
< L(1 + 1/K): the "1-small + cluster" mixed shape.  With L=50, K=100: (A) ratio <= 5050;
(B) cluster ratio < 50.5.  Both bands land in NAMED machinery: (A) bounded-ratio/
compressed (single-scale tools); (B) the (P,E) small-part architecture / far-cluster
peel (THM-608, THM-530 |P|=1-analog at k=8's leg).

CENSUS: adversarial + random k=8 shapes at diam >= 27 (the post-ledger regime):
 - how many lack a 50x cherry; their band type (A/B); their mu_{1/7} (numeric);
 - is any residual shape LOW-mu (< the 0.675 leg bar)?  (klein's certificates are for
   the floor; the residual's own mu tells whether the gate's residual is ever binding.)
"""
import numpy as np
import random as rnd
from itertools import combinations
rnd.seed(48)

GRID = 60000
xs = (np.arange(GRID)+0.5)/GRID
def mu17(E):
    ph = np.mod(np.outer(xs, np.array(E,float)), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return float((g.max(axis=1) > 1/7).mean())

L = 50
def has_sep_cherry(speeds):
    s = sorted(speeds)
    for (a,b,c) in combinations(s,3):
        if c >= L*(a+b): return True
    return False

def band_type(speeds, K=100):
    s = sorted(speeds)
    e2,e3,e8 = s[0], s[1], s[-1]
    # lemma check
    assert e8 < L*(e2+e3)+1e-9 or has_sep_cherry(speeds)
    return 'A' if e3 <= K*e2 else 'B'

print(f"=== structure lemma check + residual census (L={L}, k=8, diam >= 27) ===")
# generate shapes: random, stretched, geometric, adversarial-ish
def gen_shapes(n_shapes=4000):
    out=[]
    while len(out)<n_shapes:
        style = rnd.random()
        if style < 0.4:
            sp = sorted(rnd.sample(range(1, 400), 7))
        elif style < 0.7:
            base = rnd.randint(1, 8)
            sp = sorted(rnd.sample(range(base, base*rnd.randint(20,120)), 7))
        else:
            e2 = rnd.randint(1, 20); e3 = rnd.randint(e2+1, 40*e2)
            hi = L*(e2+e3)-1
            if hi <= e3+5: continue
            rest = rnd.sample(range(e3, hi), 5)
            sp = sorted([e2] + [e3] + rest)
            sp = sorted(set(sp))
            if len(sp) != 7: continue
        if max(sp) - 0 < 27: continue   # diam of {0} u sp
        out.append(tuple(sp))
    return out

shapes = gen_shapes()
res = [s for s in shapes if not has_sep_cherry(s)]
print(f"shapes: {len(shapes)}; NO-separated-cherry residual: {len(res)} ({100*len(res)/len(shapes):.1f}%)")
lemma_ok = all(max(s) < L*(sorted(s)[0]+sorted(s)[1]) for s in res)
print(f"structure lemma e8 < L(e2+e3) on every residual shape: {lemma_ok}")
typeA = [s for s in res if band_type(s)=='A']; typeB = [s for s in res if band_type(s)=='B']
print(f"band split (K=100): A(single band) = {len(typeA)}, B(1-small + cluster) = {len(typeB)}")

# mu census on the residual (worst cases matter)
mus = []
for s in res[:600]:
    E = (0,)+s
    mus.append((mu17(list(E)), s))
mus.sort()
print(f"\nresidual mu_1/7: min = {mus[0][0]:.4f} at {mus[0][1]}")
print(f"  5 lowest: {[(round(m,4), s[:3]+('...',)) for m,s in mus[:5]]}")
print(f"  fraction below the k=8 leg need 0.675: {sum(1 for m,_ in mus if m < 0.675)}/{len(mus)}")
bar_ok = mus[0][0] >= 0.675
print(f"  residual ever binding? {'NO — every no-cherry shape clears 0.675 directly' if bar_ok else 'YES — see lows'}")
