"""What EXACTLY is the 1.04 miss at the final row?  Report k, want, cap, want-cap."""
from math import comb
from fractions import Fraction as F
from itertools import product
import amm12592_gamma35_beam_deathstar as beam
from final_target import viol, supnorm

def reach_final(R=128, D0=3, bw=2500, ctrl=2, span=2, dedup=30, tail=48):
    d = [(3 * (R + i)) // 5 + D0 for i in range(R)]
    states = [([], beam.qpow(R - 1))]
    opts = [None] + list(range(-span, span + 1))
    for i in range(R - 1):
        late = (i >= R - 1 - tail); nxt = []
        for acc, sig in states:
            for tg in product(opts, repeat=ctrl):
                if tg[0] not in (1, -1): continue
                r = beam.step(sig, d[i], tg)
                if r is None: continue
                de, ns = r
                if not ns or abs(ns[0]) != 1: continue
                if viol(ns, R - 1 - i) > 1: continue
                key = (supnorm(ns), len(ns)) if late else (len(ns), sum(abs(c) for c in ns[:6]))
                nxt.append((key, acc + [de], ns))
        if not nxt: return None, d, i
        nxt.sort(key=lambda s: s[0])
        seen, uniq = set(), []
        for v, a, s in nxt:
            k = tuple(s[:dedup])
            if k in seen: continue
            seen.add(k); uniq.append((a, s)); 
            if len(uniq) >= bw: break
        states = uniq
    return states, d, None

states, d, died = reach_final()
dfin = d[127]
print(f"final-row degree budget d = {dfin}  (D0=3)")
best = None
for acc, sig in states:
    if len(sig) - 1 > dfin: continue
    res = list(sig) + [0] * (dfin + 2); rec = []
    for k in range(dfin, -1, -1):
        cap = comb(dfin, k); want = res[dfin - k]
        rec.append((k, want, cap))
        if abs(want) > cap or (want - cap) % 2:
            if best is None or abs(want) / cap < best[0]:
                best = (abs(want) / cap, k, want, cap, 'capacity' if abs(want) > cap else 'parity')
            break
        for t, c in enumerate(beam.basis_poly(dfin, k)): res[t] -= want * c
    else:
        if not beam.trim(res): print("SOLVED!"); break
if best:
    ratio, k, want, cap, mode = best
    print(f"best survivor fails at k={k} ({mode}):")
    print(f"   want = {want}")
    print(f"   cap  = binom({dfin},{k}) = {cap}")
    print(f"   |want|/cap = {ratio:.6f}   want - cap = {want - cap}   |want| - cap = {abs(want)-cap}")
    print(f"   parity ok: {(want - cap) % 2 == 0}")
    print(f"   how many extra degrees would fix it: binom({dfin}+j,{k}) >= |want| needs j =", end=" ")
    j = 0
    while comb(dfin + j, k) < abs(want): j += 1
    print(j)
