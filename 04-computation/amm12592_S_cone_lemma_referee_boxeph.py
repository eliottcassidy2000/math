"""LANE E1 -- machine certificates for the S-cone one-step lemmas + the
self-propagating entry condition (H4*).

All exact int; no floats in any decision.  Deterministic pseudorandomness
(fixed-seed xorshift, no random module surprises across versions).

Certified claims (each = exact algebra re-checked on many random states):

  L1 (magnitude closed form).  For all-negative junk j = -a and the
     nearest-point clamp in boxes [-2C(d'-1,t), +2C(d'-1,t-1)], one row of
     the T6 flow equals  a'_t = max(0, (K_delta a)_t - 2C(d'-1,t)) for
     t >= 1 and a'_0 = max(0, a_0 - 2);  in particular junk stays <= 0.
     (Checked against an independent clamp implementation.)

  L2 (monotonicity / comparison).  a <= b cellwise  =>  Phi(a) <= Phi(b)
     cellwise, where Phi is the L1 magnitude map (same d, delta).

  L3 (H4* self-propagation).  If  2*(2a_{t-1} + a_{t-2}) <= 2C(d-1,t) for
     all t in [2, m+2]  and  a_0 <= d - 1, then after one step (either
     delta):
       (i)   every cell t in [1, m] is non-increasing; a_0' = a_0 - 2 or 0;
       (ii)  support does not advance (no junk above m);
       (iii) every live cell t in [2, m] loses >= C(d-1, t):
             a'_t <= max(0, a_t - C(d-1,t));
       (iv)  the same inequality family holds at the next row (caps grown).

  L4 (cell-0 clock + cell-1 layer).  a'_0 = max(0, a_0 - 2) always; if
     a_0 > d - 1 (layer), cell-1 gain is <= 2*(a_0 - d) on delta=1 rows
     and cell 1 is non-increasing on delta=0 rows provided a_0 <= 2(d-1).

Output: amm12592_S_cone_lemma_referee_boxeph.{out,json}
"""
import json, os
from math import comb

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "05-knowledge", "results")

_state = [88172645463325252]
def rnd(n):
    """deterministic xorshift64* in [0, n)."""
    x = _state[0]
    x ^= (x << 13) & ((1 << 64) - 1)
    x ^= x >> 7
    x ^= (x << 17) & ((1 << 64) - 1)
    _state[0] = x
    return (x * 2685821657736338717 % (1 << 64)) % n


def step_reference(j, d_new, delta):
    """Independent implementation: transport + nearest-point clamp on the
    signed junk (the engine's semantics, written differently)."""
    w = {}
    K = (1, 1) if delta == 0 else (1, 2, 1)
    for t, v in j.items():
        for s, ks in enumerate(K):
            w[t + s] = w.get(t + s, 0) + ks * v
    jn = {}
    for t, v in sorted(w.items()):
        if v == 0:
            continue
        lo = -2 * comb(d_new - 1, t)
        hi = 2 * comb(d_new - 1, t - 1) if t >= 1 else 0
        u = min(hi, max(lo, v))
        if u != v:
            jn[t] = v - u
    return jn


def step_magnitude(a, d_new, delta):
    """The L1 closed form on magnitudes."""
    an = {}
    K = (1, 1) if delta == 0 else (1, 2, 1)
    load = {}
    for t, v in a.items():
        for s, ks in enumerate(K):
            load[t + s] = load.get(t + s, 0) + ks * v
    for t, v in load.items():
        if t == 0:
            x = max(0, v - 2)
        else:
            x = max(0, v - 2 * comb(d_new - 1, t))
        if x:
            an[t] = x
    return an


def random_state(d, m, big):
    """random all-negative junk with support in [0, m]; magnitudes up to
    big * cap."""
    j = {}
    for t in range(m + 1):
        if rnd(4) == 0 and t != 0:
            continue
        cap = 2 * comb(d - 1, t) if t >= 1 else 2
        v = rnd(max(2, big * cap)) + 1
        j[t] = -2 * ((v + 1) // 2)      # even, negative
    return j


def main():
    rep = {"L1": 0, "L2": 0, "L3": 0, "L4": 0}
    fails = []
    # ---------------- L1: closed form == reference clamp (600 trials)
    for k in range(600):
        d = 30 + rnd(200)
        m = 1 + rnd(min(20, d // 3))
        delta = rnd(2)
        j = random_state(d, m, 3)
        a = {t: -v for t, v in j.items()}
        r1 = step_reference(j, d + delta, delta)
        r2 = step_magnitude(a, d + delta, delta)
        ok = (all(v < 0 for v in r1.values()) and
              {t: -v for t, v in r1.items()} == r2)
        rep["L1"] += ok
        if not ok:
            fails.append(("L1", k, d, delta))
    # ---------------- L2: monotonicity (600 trials)
    for k in range(600):
        d = 30 + rnd(200)
        m = 1 + rnd(min(20, d // 3))
        delta = rnd(2)
        a = {t: -v for t, v in random_state(d, m, 3).items()}
        b = {t: v + 2 * rnd(v + 5) for t, v in a.items()}   # b >= a
        ra = step_magnitude(a, d + delta, delta)
        rb = step_magnitude(b, d + delta, delta)
        ok = all(ra.get(t, 0) <= rb.get(t, 0) for t in set(ra) | set(rb))
        rep["L2"] += ok
        if not ok:
            fails.append(("L2", k, d, delta))
    # ---------------- L3: H4* one-step consequences + self-propagation
    for k in range(600):
        d = 60 + rnd(400)
        m = 2 + rnd(min(25, d // 6))
        delta = rnd(2)
        # build a state satisfying H4* and a_0 <= d-1
        a = {0: 2 * rnd(d // 2)}
        for t in range(1, m + 1):
            capm = comb(d - 1, t)
            # choose a_t freely but then enforce H4* by rejection at t+1;
            # simplest: enforce 4*a_t <= C(d-1, t+1) and 2*a_t <= C(d-1,t+2)
            lim = min(capm * 4, comb(d - 1, t + 1) // 4,
                      comb(d - 1, t + 2) // 2)
            if lim < 2:
                a[t] = 0
            else:
                a[t] = 2 * rnd(lim // 2 + 1)
        # verify H4* holds (construction guarantees; assert exactly)
        H4 = all(2 * (2 * a.get(t - 1, 0) + a.get(t - 2, 0))
                 <= 2 * comb(d - 1, t) for t in range(2, m + 3))
        if not (H4 and a[0] <= d - 1):
            continue
        an = step_magnitude(a, d + delta, delta)
        ok = True
        # (i) non-increasing on [1, m]; cell 0 clock
        ok &= all(an.get(t, 0) <= a.get(t, 0) for t in range(1, m + 1))
        ok &= an.get(0, 0) == max(0, a[0] - 2)
        # (ii) no junk above m
        ok &= all(t <= m for t in an)
        # (iii) live cells in [2, m] lose >= C(d-1,t)
        ok &= all(an.get(t, 0) <= max(0, a.get(t, 0) - comb(d - 1, t))
                  for t in range(2, m + 1))
        # (iv) H4* next row (caps at d+delta)
        ok &= all(2 * (2 * an.get(t - 1, 0) + an.get(t - 2, 0))
                  <= 2 * comb(d + delta - 1, t) for t in range(2, m + 3))
        rep["L3"] += ok
        if not ok:
            fails.append(("L3", k, d, delta, m))
    # ---------------- L4: layer bookkeeping (600 trials)
    for k in range(600):
        d = 100 + rnd(400)
        delta = rnd(2)
        c0 = 2 * rnd(40)
        a = {0: d + c0 - ((d + c0) % 2), 1: 2 * rnd(3 * d)}
        an = step_magnitude(a, d + delta, delta)
        ok = an.get(0, 0) == max(0, a[0] - 2)
        if delta == 1:
            ok &= an.get(1, 0) <= a[1] + 2 * max(0, a[0] - d)
        else:
            if a[0] <= 2 * (d - 1):
                ok &= an.get(1, 0) <= a[1]
        rep["L4"] += ok
        if not ok:
            fails.append(("L4", k, d, delta))
    lines = []
    for L in ("L1", "L2", "L3", "L4"):
        lines.append(f"[{'PASS' if not [f for f in fails if f[0]==L] else 'FAIL'}] "
                     f"{L}: {rep[L]} trials ok")
    lines.append(f"fails: {fails[:10]}")
    out = "\n".join(lines)
    print(out)
    open(os.path.join(RESULTS, "amm12592_S_cone_lemma_referee_boxeph.out"),
         "w").write(out + "\n")
    json.dump({"rep": rep, "fails": fails[:50]},
              open(os.path.join(RESULTS,
                   "amm12592_S_cone_lemma_referee_boxeph.json"), "w"))


if __name__ == "__main__":
    main()
