"""d_i = floor(3(R+i)/5) + e_i with e_i BOUNDED.  Since T(n) = n+1+floor(3n/5)+e(n)
with e bounded still gives C = lim T(n)/n = 8/5, a non-constant bounded offset is
fully legitimate -- and it is a degree of freedom the fixed-D0 runs never used.
D0=3 reaches the final row and misses at k=22 by ONE degree; D0=4 stalls at row 90.
So try e = 3 everywhere with a +1 bump on the last few rows."""
from math import comb
from fractions import Fraction as F
from itertools import product
import amm12592_gamma35_beam_deathstar as beam
from final_target import viol, supnorm
import time, json

def run_e(R, evec, bw=2500, ctrl=2, span=2, dedup=30, tail=48):
    d = [(3 * (R + i)) // 5 + evec[i] for i in range(R)]
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
        if not nxt: return None, f"died row {i}", d
        nxt.sort(key=lambda s: s[0])
        seen, uniq = set(), []
        for v, a, s in nxt:
            k = tuple(s[:dedup])
            if k in seen: continue
            seen.add(k); uniq.append((a, s))
            if len(uniq) >= bw: break
        states = uniq
    dfin = d[R - 1]; best = None
    for acc, sig in states:
        if len(sig) - 1 > dfin: continue
        res = list(sig) + [0] * (dfin + 2); de = [0] * (dfin + 1); ok = True; w = 0.0
        for k in range(dfin, -1, -1):
            cap = comb(dfin, k); want = res[dfin - k]
            w = max(w, abs(want) / cap)
            if abs(want) > cap or (want - cap) % 2: ok = False; break
            de[k] = want
            for t, c in enumerate(beam.basis_poly(dfin, k)): res[t] -= want * c
        if ok and not beam.trim(res): return acc + [de], "SOLVED", d
        if best is None or w < best: best = w
    return None, f"final row; best |want|/cap = {best:.4g}", d

def verify_e(R, sol, d):
    acc = [0] * (4 * R + 8)
    for i, de in enumerate(sol):
        for k, v in enumerate(de):
            assert abs(v) <= comb(d[i], k) and (v - comb(d[i], k)) % 2 == 0
            if v:
                for t, c in enumerate(beam.basis_poly(d[i], k)): acc[i + t] += v * c
    return beam.trim(acc) == beam.trim(beam.qpow(R - 1))

R = 128
trials = []
for bump, nlast in [(1,1),(1,2),(1,4),(1,8),(2,1),(2,4),(1,16),(2,16)]:
    e = [3]*R
    for j in range(R-nlast, R): e[j] = 3 + bump
    trials.append((f"e=3, last {nlast} rows +{bump}", e))
for name, e in trials:
    t0=time.time()
    sol, msg, d = run_e(R, e)
    ok = verify_e(R, sol, d) if sol else False
    print(f"{name:28s}: {msg}  verify={ok}  ({time.time()-t0:.0f}s)", flush=True)
    if ok:
        print(f"*** R=128 CLOSES at gamma=3/5 with bounded offset ({name}) -> C = 8/5 for n <= 255 ***")
        json.dump({"R":R,"e":e,"sol":sol}, open("sol128.json","w")); break
