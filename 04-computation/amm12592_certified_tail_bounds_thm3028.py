"""Referee for THM-3028: (SUP), (PROP), the float trap, filter soundness, D0 table."""
from math import comb
from fractions import Fraction as F
from itertools import product
import random
import amm12592_gamma35_beam_deathstar as beam

def poly_of(delta, d):
    p = [0] * (d + 1)
    for k, v in enumerate(delta):
        if v:
            for t, c in enumerate(beam.basis_poly(d, k)): p[t] += v * c
    return p

def ev_exact(p, t):
    v = F(0)
    for c in reversed(p): v = v * t + c
    return v

def ev_float(p, t):
    v = 0.0
    for c in reversed(p): v = v * t + c
    return v

TS = [F(i, 64) for i in range(64)]

print("(SUP) every admissible block has |Delta(t)| <= 1 on [0,1]")
random.seed(3); worst = F(0)
for _ in range(2000):
    d = random.randint(1, 12)
    delta = [comb(d, k) - 2 * random.randint(0, comb(d, k)) for k in range(d + 1)]
    p = poly_of(delta, d)
    for t in TS: worst = max(worst, abs(ev_exact(p, t)))
print(f"   2000 random admissible blocks (d<=12): max |Delta(t)| = {worst} "
      f"(<= 1: {worst <= 1}); bound attained: {worst == 1}")

print("\n(PROP) |sigma_{j-1}(t)| <= (1-t^{R-j})/(1-t) -- EXACT, on true solution paths")
def path_ratio(R, D0=0, use_float=False):
    sol, _ = beam.solve(R, beam=250, ctrl=2, span=2, D0=D0)
    assert sol and beam.verify(R, sol, D0=D0)
    d = [(3 * (R + i)) // 5 + D0 for i in range(R)]
    sig = beam.qpow(R - 1); worst = F(0) if not use_float else 0.0
    for i in range(R):
        rows = R - i
        for t in TS:
            b = F(1) if t == 0 else (1 - t ** rows) / (1 - t)
            r = (abs(ev_exact(sig, t)) / b) if not use_float else (abs(ev_float(sig, float(t))) / float(b))
            worst = max(worst, r)
        res = list(sig) + [0] * (d[i] + 4)
        for k, v in enumerate(sol[i]):
            if v:
                for tt, c in enumerate(beam.basis_poly(d[i], k)): res[tt] -= v * c
        assert res[0] == 0
        sig = beam.trim(res[1:])
    return worst

for R in (16, 32, 64):
    w = path_ratio(R)
    print(f"   R={R:3d}: exact max violation ratio along the TRUE path = {w}  "
          f"(satisfies and SATURATES (PROP): {w == 1})")

print("\n(FLOAT) the same quantity in double precision is meaningless")
print(f"   R=64: exact = {path_ratio(64)}   float = {path_ratio(64, use_float=True):.6f}")
print("   => all (SUP)/(PROP) testing must be exact.")

print("\nFILTER SOUNDNESS + D0 death-row table (exact arithmetic throughout)")
TP = [F(1, 2), F(1, 4), F(3, 4), F(1, 8)]
def viol(sig, rows):
    if not sig: return F(0)
    D = len(sig) - 1; w = F(0)
    for t in TP:
        a, b = t.numerator, t.denominator
        num = sum(c * a ** k * b ** (D - k) for k, c in enumerate(sig) if c)
        r = F(abs(num), b ** D) / ((1 - t ** rows) / (1 - t))
        if r > w: w = r
    return w

def run(R, D0=0, bw=400, ctrl=2, span=2, dedup=30):
    d = [(3 * (R + i)) // 5 + D0 for i in range(R)]
    states = [([], beam.qpow(R - 1))]; opts = [None] + list(range(-span, span + 1))
    for i in range(R - 1):
        nxt, best = [], None
        for acc, sig in states:
            for tg in product(opts, repeat=ctrl):
                if tg[0] not in (1, -1): continue
                r = beam.step(sig, d[i], tg)
                if r is None: continue
                de, ns = r
                if not ns or abs(ns[0]) != 1: continue
                v = viol(ns, R - 1 - i)
                if best is None or v < best: best = v
                if v > 1: continue
                nxt.append(((len(ns), sum(abs(c) for c in ns[:6])), acc + [de], ns))
        if not nxt: return None, i, (float(best) if best is not None else None)
        nxt.sort(key=lambda s: s[0])
        seen, uniq = set(), []
        for v, a, s in nxt:
            k = tuple(s[:dedup])
            if k in seen: continue
            seen.add(k); uniq.append((a, s))
            if len(uniq) >= bw: break
        states = uniq
    dfin = d[R - 1]
    for acc, sig in states:
        if len(sig) - 1 > dfin: continue
        res = list(sig) + [0] * (dfin + 2); de = [0] * (dfin + 1); ok = True
        for k in range(dfin, -1, -1):
            cap = comb(dfin, k); want = res[dfin - k]
            if abs(want) > cap or (want - cap) % 2: ok = False; break
            de[k] = want
            for t, c in enumerate(beam.basis_poly(dfin, k)): res[t] -= want * c
        if ok and not beam.trim(res): return acc + [de], None, None
    return None, R, None

for R in (16, 32, 64):
    sol, died, _ = run(R)
    print(f"   soundness R={R:3d}: {'SOLVES through the filter' if sol else f'died row {died}'}"
          f"  verify={beam.verify(R, sol) if sol else False}")
print("   R=128 death row vs constant slack D0 (C = 8/5 for every fixed D0):")
for D0 in (0, 1, 2, 4, 8):
    sol, died, best = run(128, D0=D0)
    print(f"      D0={D0}: {'SOLVED' if sol else f'died row {died} of 128, best achievable ratio {best:.4g}'}")
