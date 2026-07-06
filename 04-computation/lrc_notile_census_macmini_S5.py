#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S5 -- HYP-4282 item 3b: the DISTINCT-COMBS-DON'T-TILE census
at band 2/25.  Systematic evidence for opus-S96's lemma in the >=7 regime.

QUESTION: is there ANY >=7-element frequency multiset that CAN tile the circle
(free fraction -> 0) at band rho = 2/25 under the WORST (adversarial) phase?
If NO pattern tiles, the >=7-lifted (A) residual is EMPTY unconditionally
(distinct-frequency rigidity alone forbids covering) -- opus's lemma holds and
the phase orbit is not even needed.

METHOD: adversarial coordinate-descent MINIMIZING the free fraction over phase
vectors (the covering adversary), multiple restarts.  A pattern "tiles" if the
adversary drives free fraction below eps.  Sweep structured + random >=7
patterns (scale-normalized: phi is scale-invariant, opus-S96).
"""
import random
from fractions import Fraction as F

RHO = 2 / 25

def dist(x):
    y = x - round(x)
    return abs(y)

def free_frac(freqs, phi, rho=RHO, sgrid=2400):
    cnt = 0
    for j in range(sgrid):
        s = (j + 0.5) / sgrid
        ok = True
        for f, p in zip(freqs, phi):
            if dist(f * s + p) < rho:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt / sgrid

def adversarial_min(freqs, rho=RHO, restarts=14, sweeps=60, sgrid=1600, seed=0):
    """min free fraction the covering adversary can force."""
    rng = random.Random(seed)
    best = 1.0
    for _ in range(restarts):
        phi = [rng.random() for _ in freqs]
        cur = free_frac(freqs, phi, rho, sgrid)
        step = 0.25
        for _sw in range(sweeps):
            improved = False
            for i in range(len(freqs)):
                for d in (step, -step):
                    phi2 = list(phi); phi2[i] = (phi2[i] + d) % 1
                    f = free_frac(freqs, phi2, rho, sgrid)
                    if f < cur - 1e-9:
                        cur, phi = f, phi2; improved = True
            if not improved:
                step /= 2
                if step < 1e-4:
                    break
        best = min(best, cur)
        if best <= 1e-4:
            break
    return best

def patterns():
    """>=7-element frequency multisets (scale-normalized reps)."""
    out = []
    # consecutive runs of length 7..12
    for L in range(7, 13):
        out.append((f"consec{L} [1..{L}]", list(range(1, L + 1))))
    # APs step 2, 3
    for step in (2, 3, 5):
        for L in (7, 9, 12):
            out.append((f"AP{L}-step{step}", [1 + step * i for i in range(L)]))
    # coprime primes
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for L in (7, 9, 12):
        out.append((f"primes{L}", primes[:L]))
    # common-factor (ray-like): S * base
    base7 = [1, 2, 3, 4, 5, 6, 7]
    for S in (2, 3, 5):
        out.append((f"{S}x[1..7] (ray)", [S * b for b in base7]))
    # the pole-necessity cluster
    out.append(("pole-7-cluster", [20, 21, 24, 25, 45, 46, 66]))
    out.append(("pole+lift12", [20, 21, 24, 25, 45, 46, 66, 1, 2, 3, 4, 5]))
    # near-equal (worst-tiling per opus): w..w+6
    out.append(("near-equal 100..106", list(range(100, 107))))
    # halves/doubling structure
    out.append(("doubling 1,2,4,8,16,32,64", [1, 2, 4, 8, 16, 32, 64]))
    # random >=7 sparse
    rng = random.Random(42)
    for k in range(8):
        L = rng.randint(7, 12)
        out.append((f"rand{k} (L={L})", sorted(rng.sample(range(1, 40), L))))
    return out

if __name__ == "__main__":
    print("=" * 78)
    print("DISTINCT-COMBS-DON'T-TILE CENSUS at band 2/25 (opus lemma, >=7 regime)")
    print("adversary MINIMIZES free fraction; pattern TILES iff it can reach ~0")
    print("=" * 78)
    rows = []
    tiled = []
    for name, freqs in patterns():
        adv = adversarial_min(freqs, seed=hash(name) & 0xffff)
        tot = len(freqs) * 2 * RHO
        status = "TILES (free->0!)" if adv < 5e-3 else "no-tile"
        if adv < 5e-3:
            tiled.append(name)
        print(f"  {name:26s} |L|={len(freqs):2d} meas={tot:.2f}  "
              f"phi_worst = {adv:.4f}   [{status}]")
        rows.append((name, adv))
    print()
    print("=" * 78)
    mn = min(r[1] for r in rows)
    argmn = [r[0] for r in rows if r[1] == mn][0]
    print(f"  patterns tested: {len(rows)};  MIN phi_worst = {mn:.4f} (at {argmn})")
    if tiled:
        print(f"  *** {len(tiled)} PATTERNS TILE: {tiled}")
        print(f"  => these need the phase-orbit / opus-transport argument (structured phases).")
    else:
        print(f"  ZERO patterns tile: every >=7 frequency multiset leaves the circle")
        print(f"  >= {mn:.3f} uncovered even under the WORST phase.  Strong evidence for")
        print(f"  opus-S96's 'distinct combs never tile' at band 2/25 in the >=7 regime:")
        print(f"  the (A) >=7-lifted residual is EMPTY unconditionally (no phase orbit needed).")
    print(f"  NOTE: adversarial search is heuristic (local minima) -- evidence, not proof;")
    print(f"  opus-S97's transport is the proof on the ray sub-shape (S*base above).")
