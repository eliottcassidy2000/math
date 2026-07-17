#!/usr/bin/env python3
"""death-star-2026-07-16-S34 (HYP-7171): referee for LRCDenseCoreRelationTrap.lean and
LRCPairLiftDichotomy.lean.

TRAP (A1/A2): on ChainDenseCore families (ratios >= 3 at all pair indices > j),
  A1: no vanishing relation with top support t >= j+2 and below-mass sum|c_i| <= 2;
  A2: no unit-coefficient vanishing relation with top support t >= j+4.
  Both follow from the proved inequalities 2*w[t-1] < w[t] (t >= j+2) and
  sum_{i<t} w[i] < w[t] (t >= j+4); referee checks the inequalities AND brute-forces
  actual relations on random dense-core families.

PAIR LIFT: lonely_or_tripleCore -- singles split (THM-937) OR pair-at-dense split
  (entry (13-j)w[j] >= 13(j+1)w[j-1], free at j=0; tail w[j+2] >= 21*w[j+1], vacuous
  at j=11) OR TripleDenseCore. Referee: exhaustiveness + the shrink quantified
  (dense-core families now closed by the dodge)."""
import random
from itertools import combinations, product

def dense_core_certificate(w):
    bad = [j for j in range(12) if w[j + 1] < 3 * w[j]]
    if not bad:
        return None
    js = max(bad)
    if not all(3 * w[k] <= w[k + 1] for k in range(js + 1, 12)):
        return None
    if not all(2 * (12 - k) * w[k + 1] < 21 * (k + 2) * w[k] for k in range(js, 12)):
        return None
    return js

def random_dense_core(rnd):
    """build a random ChainDenseCore family: dense pair at j, controlled ladder above
    (each step ratio in [3, entry-fail cap)), arbitrary below."""
    for _ in range(400):
        j = rnd.randint(0, 11)
        w = sorted(rnd.randint(1, 60) for _ in range(j))  # below the pair
        base = (w[-1] if w else 1) * rnd.randint(1, 4) + rnd.randint(0, 9)
        w.append(base)                                     # w[j]
        w.append(base + max(1, int(base * rnd.uniform(0.2, 1.9))))  # dense pair partner < 3x
        if w[j + 1] >= 3 * w[j]:
            continue
        ok = True
        for k in range(j + 1, 12):
            # ratio in [3, cap_k) where entry fails: 2(12-k)w[k+1] < 21(k+2)w[k]
            capnum, capden = 21 * (k + 2), 2 * (12 - k)
            lo = 3 * w[-1]
            hi = (capnum * w[-1]) // capden  # strict <
            if lo > hi - 1:
                ok = False
                break
            w.append(rnd.randint(lo, max(lo, hi - 1)))
        if not ok or len(w) != 13:
            continue
        if dense_core_certificate(w) == j and len(set(w)) == 13:
            return j, w
    return None

def referee_trap(trials=4000, seed=34):
    rnd = random.Random(seed)
    okA1 = okA2 = okineq = True
    n = 0
    while n < trials:
        got = random_dense_core(rnd)
        if got is None:
            continue
        j, w = got
        n += 1
        # proved inequalities
        for t in range(j + 2, 13):
            if not 2 * w[t - 1] < w[t]:
                okineq = False
        for t in range(j + 4, 13):
            if not sum(w[:t]) < w[t]:
                okineq = False
        # brute A1: top t >= j+2, below-mass <= 2 (patterns: one below +-1/+-2, or two below +-1)
        for t in range(j + 2, 13):
            for s in range(t):
                for c in (1, -1, 2, -2):
                    if w[t] + c * w[s] == 0 or -w[t] + c * w[s] == 0:
                        okA1 = False
            for s1, s2 in combinations(range(t), 2):
                for c1, c2 in product((1, -1), repeat=2):
                    if w[t] + c1 * w[s1] + c2 * w[s2] == 0:
                        okA1 = False
        # brute A2 (sampled): unit coeffs, top t >= j+4
        for t in range(j + 4, 13):
            for _ in range(200):
                cs = [rnd.choice((-1, 0, 1)) for _ in range(t)]
                if w[t] + sum(c * x for c, x in zip(cs, w[:t])) == 0:
                    okA2 = False
    print(f"trap referee ({trials} dense-core families): inequalities "
          f"{'PASS' if okineq else 'FAIL'}; A1 brute {'PASS' if okA1 else 'FAIL'}; "
          f"A2 sampled {'PASS' if okA2 else 'FAIL'}")

def singles_split_works(w, k):
    for jj in range(k, 12):
        if w[jj + 1] < 3 * w[jj]:
            return False
    if k == 0:
        return True
    return 21 * (k + 1) * w[k - 1] <= 2 * (13 - k) * w[k]

def pair_split_works(w, j):
    """pair at dense position j: entry + tail-jump + ladder above j+2."""
    if not w[j + 1] < 3 * w[j]:
        return False
    if j >= 1 and not 13 * (j + 1) * w[j - 1] <= (13 - j) * w[j]:
        return False
    if j + 2 <= 12:
        if not 21 * w[j + 1] <= w[j + 2]:
            return False
        for k in range(j + 2, 12):
            if w[k + 1] < 3 * w[k]:
                return False
    return True

def triple_core(w):
    js = dense_core_certificate(w)
    if js is None:
        return None
    tail_fail = (js + 2 <= 12) and (w[js + 2] < 21 * w[js + 1])
    entry_fail = (js >= 1) and ((13 - js) * w[js] < 13 * (js + 1) * w[js - 1])
    return (js, tail_fail, entry_fail) if (tail_fail or entry_fail) else None

def referee_pairlift(trials=200000, seed=134):
    rnd = random.Random(seed)
    ok = True
    n_s = n_p = n_t = 0
    for _ in range(trials):
        style = rnd.random()
        if style < 0.35:
            w = sorted(rnd.randint(1, 10**6) for _ in range(13))
        elif style < 0.6:
            base = rnd.randint(1, 40)
            w = [base]
            for _ in range(12):
                w.append(max(w[-1] + 1, int(w[-1] * rnd.uniform(1.5, 30))))
        else:
            got = random_dense_core(rnd)
            if got is None:
                continue
            w = got[1]
        if any(singles_split_works(w, k) for k in range(13)):
            n_s += 1
            continue
        js = dense_core_certificate(w)
        if js is None:
            ok = False
            print(f"  FAIL: no split, no dense cert: {w}")
            continue
        if pair_split_works(w, js):
            n_p += 1
            continue
        if triple_core(w) is None:
            ok = False
            print(f"  FAIL: dense, dodge fails, no triple cert: {w}")
        else:
            n_t += 1
    tot = n_s + n_p + n_t
    print(f"pair-lift referee: {'PASS' if ok else 'FAIL'} "
          f"[singles {n_s} | pair-dodge {n_p} | triple-core {n_t}] of {tot}")
    if n_p + n_t:
        print(f"  the dodge closes {n_p}/{n_p+n_t} = {100*n_p/(n_p+n_t):.1f}% of "
              f"the former dense core in this mix")

def planted_dodges(trials=20000, seed=234):
    """families built to pass the pair split must all close."""
    rnd = random.Random(seed)
    ok = True
    n = 0
    while n < trials:
        j = rnd.randint(0, 11)
        w = sorted(rnd.sample(range(1, 50), j)) if j else []
        base = (w[-1] if w else 1) * (13 * (j + 1) // max(1, 13 - j) + 1) + rnd.randint(0, 5)
        w.append(base)
        w.append(base + rnd.randint(1, max(1, 2 * base - 1)))
        if w[j + 1] >= 3 * w[j]:
            continue
        cur = 21 * w[j + 1] + rnd.randint(0, 8)
        for _ in range(11 - j):
            w.append(cur)
            cur = 3 * cur + rnd.randint(0, cur)
        w = w[:13]
        if len(w) != 13 or len(set(w)) != 13:
            continue
        n += 1
        if not (any(singles_split_works(w, k) for k in range(13))
                or pair_split_works(w, dense_core_certificate(w) if dense_core_certificate(w) is not None else j)):
            ok = False
            print(f"  FAIL planted dodge not closed: j={j} {w}")
    print(f"planted pair-dodge families ({trials}): {'PASS' if ok else 'FAIL'}")

if __name__ == "__main__":
    referee_trap()
    referee_pairlift()
    planted_dodges()
