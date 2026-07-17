#!/usr/bin/env python3
"""death-star-2026-07-16-S33 (HYP-7161): referee for LRCChainDichotomy.lean — the
split-search dichotomy wiring the S19 chain engine (cite_chain_lonely) into the grand
assembly's ResidualObligation.

THE LEAN THEOREM (lonely_or_denseCore): for sorted positive w_0 <= ... <= w_12, either
  (SPLIT) some k in {0..12} has: all consecutive ratios >= 3 at pair indices >= k, and
          entry fee OK: k = 0 free (B=1), else 21*(k+1)*w_{k-1} <= 2*(13-k)*w_k
          [Fin-12 pair index k' = k-1: 21*(k'+2)*w_{k'} <= 2*(12-k')*w_{k'+1}]
  (DENSE) some pair j has ratio < 3, all pairs above j have ratio >= 3, and the entry
          fee FAILS at every split above j: 2*(12-k)*w_{k+1} < 21*(k+2)*w_k for k >= j.

REFEREE: (1) exhaustiveness+exclusivity of the constructed certificate on random
tuples (the Lean proof picks j* = LAST bad ratio; check the certificate verbatim);
(2) the entry-constant table c_m = 21(m+1)/(2(13-m)) for split position m; (3) the
dense-core cap: above the dense pair the total growth is bounded by prod of the
per-step caps — quantify the middle band; (4) MEASURE: fraction of random
residual-shaped families (compressed, gapped) that fall DENSE vs SPLIT — how much of
the residual the new wire machine-closes."""
import random
from fractions import Fraction as Fr

def split_works(w, k):
    """Lean's split at position k (cited 0..k-1, chain k..12)."""
    # ratios >= 3 at all pair indices j >= k (pairs (j, j+1), j in [0..11])
    for j in range(k, 12):
        if w[j + 1] < 3 * w[j]:
            return False
    if k == 0:
        return True  # B = 1: 26*w0 >= 21 always (w0 >= 1)
    return 21 * (k + 1) * w[k - 1] <= 2 * (13 - k) * w[k]

def dense_certificate(w):
    """The Lean proof's construction: j* = last bad-ratio pair; check the three clauses."""
    bad = [j for j in range(12) if w[j + 1] < 3 * w[j]]
    if not bad:
        return None
    js = max(bad)
    c1 = w[js + 1] < 3 * w[js]
    c2 = all(3 * w[k] <= w[k + 1] for k in range(js + 1, 12))
    c3 = all(2 * (12 - k) * w[k + 1] < 21 * (k + 2) * w[k] for k in range(js, 12))
    return (js, c1, c2, c3)

def referee_dichotomy(trials=200000, seed=33):
    rnd = random.Random(seed)
    ok = True
    n_split = n_dense = 0
    for _ in range(trials):
        style = rnd.random()
        if style < 0.4:   # generic magnitudes
            w = sorted(rnd.randint(1, 10**6) for _ in range(13))
        elif style < 0.7: # geometric-ish chains (near the dichotomy boundary)
            base = rnd.randint(1, 50)
            w = [base]
            for _ in range(12):
                w.append(max(w[-1] + 1, int(w[-1] * rnd.uniform(1.2, 4.5))))
        else:             # clustered (dense blocks + jumps)
            w = [rnd.randint(1, 20)]
            for _ in range(12):
                w.append(w[-1] * rnd.choice([1, 1, 1, 2, 3, 7, 15, 200]) + rnd.randint(0, 3))
            w = sorted(w)
        splits = [k for k in range(13) if split_works(w, k)]
        cert = dense_certificate(w)
        if splits:
            n_split += 1
            # Lean picks a specific k; existence is what the theorem needs
        else:
            n_dense += 1
            if cert is None or not (cert[1] and cert[2] and cert[3]):
                ok = False
                print(f"  FAIL: no split and no dense certificate: {w} cert={cert}")
    print(f"referee 1 (dichotomy exhaustive, {trials} tuples): "
          f"{'PASS' if ok else 'FAIL'}  [split {n_split} | dense {n_dense}]")
    return ok

def entry_constants():
    print("referee 2: entry constants c_m = 21(m+1)/(2(13-m)) (split position m, B = w_{m-1})")
    row = []
    for m in range(1, 13):
        c = Fr(21 * (m + 1), 2 * (13 - m))
        row.append(f"c_{m}={c}={float(c):.3g}")
    print("  " + "  ".join(row))
    # dense-core growth cap above the dense pair at position j: prod_{k=j+1..12} c_k
    print("  dense-core top/dense-pair cap prod c_m (m = j+1..12):")
    for j in range(0, 12):
        cap = 1
        for m in range(j + 1, 13):
            cap *= Fr(21 * (m + 1), 2 * (13 - m))
        print(f"    j={j}: w12/w_j < {float(cap):.4g}")

def residual_shaped_fraction(trials=100000, seed=133):
    """Random families conditioned residual-shaped (positive, distinct, compressed:
    every speed <= 13x some other; gapped: max > 13*min; max >= 23): SPLIT vs DENSE."""
    rnd = random.Random(seed)
    n = n_split = 0
    while n < trials:
        w = sorted(rnd.randint(1, 10**5) for _ in range(13))
        if len(set(w)) < 13 or w[12] <= 13 * w[0] or w[12] < 23:
            continue
        # compressed = no 13x-dominant runner: top <= 13 * second
        if w[12] > 13 * w[11]:
            continue
        n += 1
        if any(split_works(w, k) for k in range(13)):
            n_split += 1
    print(f"referee 4 (residual-shaped random, {trials}): split-closed {n_split} "
          f"({100*n_split/trials:.2f}%), dense-core {trials-n_split} "
          f"({100*(trials-n_split)/trials:.2f}%)")
    print("  NOTE: uniform-random magnitudes rarely admit ratio-3 tails; the wire's value")
    print("  is closing EVERY lacunary/mixed family exactly, not bulk percentage.")

def lacunary_check():
    """Families with a ratio-3 tail anywhere are all split-closed (the wire's guarantee)."""
    rnd = random.Random(233)
    ok = True
    for _ in range(50000):
        k = rnd.randint(0, 12)
        base = sorted(rnd.randint(1, 40) for _ in range(k))
        w = list(base)
        prev = (base[-1] if base else 0)
        # entry: for k>=1 need 21(k+1)w_{k-1} <= 2(13-k)w_k
        need = 1 if k == 0 else -(-21 * (k + 1) * prev // (2 * (13 - k)))
        cur = max(need, prev * 3, 1) + rnd.randint(0, 5)
        for _ in range(13 - k):
            w.append(cur)
            cur = cur * 3 + rnd.randint(0, cur // 2)
        w = w[:13]
        if not any(split_works(w, kk) for kk in range(13)):
            ok = False
            print(f"  FAIL planted-lacunary not closed: k={k} w={w}")
    print(f"referee 3 (planted citable tails all split-closed, 50000): {'PASS' if ok else 'FAIL'}")

if __name__ == "__main__":
    referee_dichotomy()
    entry_constants()
    lacunary_check()
    residual_shaped_fraction()
