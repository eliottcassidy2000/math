"""kps-S124: work the final rung -- exact-load vs live-ruler-supply over the covering residual.
THM-681: W0<=0.08 (<=1 global exact relation) => live via restricted-sumset; W0>0.08 => Freiman ladder.
klein-S222: absolute |OffLine| DIVERGES (signed). So check the DICHOTOMY empirically over covering sets:
does every not-near-AP covering family still have a LIVE tall pair-sum ruler (LM>0)? and the E3 relation."""
from itertools import combinations
from math import gcd
def is_cov(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def prim(S):
    g=0
    for s in S:
        g=gcd(g,s)
        if g==1: break
    return g==1
def E3(S):  # Schur triples a+b=c (ordered a,b), the exact relation content
    Sset=set(S); return sum(1 for a in S for b in S if a+b in Sset)
def n_doubling(S):  # v_b = 2 v_a
    Sset=set(S); return sum(1 for a in S if 2*a in Sset)
def LM(S,q):  # # live multipliers p in {1..q-1}: all residues (v*p mod q) in band [q/14,13q/14]
    return sum(1 for p in range(1,q) if all(q<=14*((v*p)%q)<=13*q for v in S))
def tall_rulers(S):
    Vmax=max(S); return sorted({a+b for i,a in enumerate(S) for b in S[i:] if a+b>Vmax})
def min_live_over_tall(S):
    best=0; anylive=0
    for q in tall_rulers(S):
        lm=LM(S,q)
        if lm>0: anylive+=1; best=max(best,lm)
    return anylive, best  # (# live tall rulers, max LM)

cov=[list(S) for S in combinations(range(1,19),13) if is_cov(S) and prim(S)]
print(f"covering primitive [1,18]: {len(cov)} families")
# global exact relation count = "exact load" proxy (Schur triples + doublings)
zero_live=0; minlivecount=99; e3max=0; hardest=[]
for S in cov:
    exact = E3(S) + n_doubling(S)   # global exact unit relations (t=0)
    nlive, maxlm = min_live_over_tall(S)
    if nlive==0: zero_live+=1; hardest.append((S,exact))
    minlivecount=min(minlivecount,nlive)
    e3max=max(e3max,E3(S))
    # record the ones with >=2 exact relations (final-rung candidates)
print(f"families with ZERO live tall ruler (a-priori supply FAILS): {zero_live}")
print(f"min # live tall rulers over all families: {minlivecount}")
print(f"max E3 over covering [1,18]: {e3max}")
# the final rung: >=2 global exact relations. how many, and are they all still live?
rung=[S for S in cov if E3(S)+n_doubling(S)>=2]
rung_live=[S for S in rung if min_live_over_tall(S)[0]>0]
print(f"\nFINAL-RUNG (>=2 global exact relations): {len(rung)} families; all live? {len(rung_live)}=={len(rung)}")
# does high E3 correlate with FEWER live rulers (but still >0)?
import statistics
by_e3={}
for S in cov:
    e=E3(S); by_e3.setdefault(e,[]).append(min_live_over_tall(S)[0])
print("\nE3 -> (avg #live tall rulers, min):")
for e in sorted(by_e3)[:12]:
    vals=by_e3[e]; print(f"  E3={e:3d}: n={len(vals):3d} avg_live={statistics.mean(vals):.1f} min_live={min(vals)}")
