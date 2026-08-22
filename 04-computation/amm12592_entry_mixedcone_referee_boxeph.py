"""LANE F1 -- machine certificates for the ENTRY-side lemmas:

  EN-D  (signed-part decoupling): for ANY mixed-sign junk j = p - n
        (p_t = j_t^+, n_t = j_t^-), one post-feed row of rule A obeys
          n'_t <= max(0, (K_delta n)_t - 2C(d'-1,t))   (t >= 1)
          n'_0 <= max(0, n_0 - 2)
          p'_t <= max(0, (K_delta p)_t - 2C(d'-1,t-1)) (t >= 1)
          p'_0 <= p_0
        cellwise (junk' = load - clamp, boxes [-2C(d'-1,t), +2C(d'-1,t-1)],
        cell 0 box [-2, 0]).

  EN-M  (monotonicity of the positive majorant flow Phi_P): p <= q cellwise
        implies Phi_P(p) <= Phi_P(q) cellwise, where
        Phi_P(p)_0 = p_0, Phi_P(p)_t = max(0, (K p)_t - 2C(d'-1,t-1)).

  EN-C  (mixed-cone one-step, the S-cone-fc+/- induction step): if at
        degree D0 the parts satisfy
          (G1) supports in [0, m], m + 2 < D0, p_0 = 0,
          (G2) n_0 <= D0 - 1,
          (G3-) 2n_{t-1} + n_{t-2} <= 2C(D0-1, t)   for t in [2, m+2],
          (G3+) 2p_{t-1} + p_{t-2} <= 2C(D0-1, t-1) for t in [2, m+2],
        then after one row at degree D' = D0 + delta (either delta):
          (a) supports stay in [0, m]; no junk beyond m (front frozen);
          (b) n'_t <= n_t and p'_t <= p_t for all t >= 1; p'_0 = 0;
              n'_0 = max(0, n_0 - 2);
          (c) G3- and G3+ still hold (capref frozen at D0);
          (d) clock decays: t in [2,m]:
                n'_t <= max(0, n_t - (2C(D'-1,t)   - 2C(D0-1,t)))
                p'_t <= max(0, p_t - (2C(D'-1,t-1) - 2C(D0-1,t-1)))
              t = 1:
                n'_1 <= max(0, n_1 - max(0, 2(D'-1) - (1+delta) n_0))
                p'_1 <= max(0, p_1 - 2).

600 exact pseudorandom trials per lemma, INDEPENDENT clamp implementation
(no imports from the engines).  Exact ints only; no numpy/sympy/floats.
Output -> 05-knowledge/results/amm12592_entry_mixedcone_referee_boxeph.{out,json}
"""
import json, os, random
from math import comb

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "05-knowledge", "results")
rng = random.Random(0xF1E17)


def clamp_row(j, d, delta):
    """Independent implementation: one post-feed row of rule A.
    Returns junk' dict.  Boxes: t=0 -> [-2,0]; t>=1 ->
    [-2C(d'-1,t), +2C(d'-1,t-1)], d' = d + delta."""
    dp = d + delta
    K = (1, 1) if delta == 0 else (1, 2, 1)
    w = {}
    for t, v in j.items():
        for s, k in enumerate(K):
            w[t + s] = w.get(t + s, 0) + k * v
    jn = {}
    for t, v in w.items():
        if t == 0:
            lo, hi = -2, 0
        else:
            lo, hi = -2 * comb(dp - 1, t), 2 * comb(dp - 1, t - 1)
        u = min(hi, max(lo, v))
        if v != u:
            jn[t] = v - u
    return jn


def parts(j, tmax):
    p = [max(0, j.get(t, 0)) for t in range(tmax + 1)]
    n = [max(0, -j.get(t, 0)) for t in range(tmax + 1)]
    return p, n


def kervec(a, delta, tmax):
    K = (1, 1) if delta == 0 else (1, 2, 1)
    out = [0] * (tmax + 1)
    for t in range(tmax + 1):
        for s, k in enumerate(K):
            if 0 <= t - s < len(a):
                out[t] += k * a[t - s]
    return out


def ev(x):
    return x - (x % 2)


def trial_END():
    d = rng.randrange(20, 300)
    delta = rng.randrange(2)
    dp = d + delta
    m = rng.randrange(3, max(4, d // 2))
    j = {}
    for t in range(m + 1):
        if rng.random() < 0.75:
            mag = rng.randrange(0, 4 * comb(d, min(t + 1, d))) if t else \
                rng.randrange(0, 3 * d)
            sgn = -1 if (t == 0 or rng.random() < 0.5) else 1
            v = ev(sgn * mag)
            if v:
                j[t] = v
    jn = clamp_row(j, d, delta)
    tmax = m + 3
    p, n = parts(j, tmax)
    pn, nn = parts(jn, tmax)
    Kn = kervec(n, delta, tmax)
    Kp = kervec(p, delta, tmax)
    for t in range(tmax + 1):
        if t == 0:
            if not (nn[0] <= max(0, n[0] - 2) and pn[0] <= p[0]):
                return False
        else:
            if not nn[t] <= max(0, Kn[t] - 2 * comb(dp - 1, t)):
                return False
            if not pn[t] <= max(0, Kp[t] - 2 * comb(dp - 1, t - 1)):
                return False
    return True


def trial_ENM():
    d = rng.randrange(20, 200)
    delta = rng.randrange(2)
    dp = d + delta
    m = rng.randrange(3, max(4, d // 3))
    q = [ev(rng.randrange(0, 4 * comb(d, min(t + 1, d)) + 2))
         for t in range(m + 1)]
    p = [ev(rng.randrange(0, qq + 1)) if qq else 0 for qq in q]
    tmax = m + 3
    Kp = kervec(p, delta, tmax)
    Kq = kervec(q, delta, tmax)
    for t in range(1, tmax + 1):
        Pp = max(0, Kp[t] - 2 * comb(dp - 1, t - 1))
        Pq = max(0, Kq[t] - 2 * comb(dp - 1, t - 1))
        if Pp > Pq:
            return False
    return p[0] <= q[0]


def gen_cone_state(D0, m):
    """Random mixed-sign state satisfying G1/G2/G3+/- at degree D0."""
    n = [0] * (m + 3)
    p = [0] * (m + 3)
    n[0] = ev(rng.randrange(0, D0))              # G2, sign at 0 forced -
    for s in range(1, m + 1):
        if rng.random() < 0.5:                   # negative cell
            ub = (2 * comb(D0 - 1, s + 1) - n[s - 1]) // 2
            if s == m:
                ub = min(ub, 2 * comb(D0 - 1, m + 2))
            n[s] = ev(rng.randrange(0, max(1, ub + 1)))
        else:                                    # positive cell
            ub = (2 * comb(D0 - 1, s) - p[s - 1]) // 2
            if s == m:
                ub = min(ub, 2 * comb(D0 - 1, m + 1))
            p[s] = ev(rng.randrange(0, max(1, ub + 1)))
    j = {}
    for t in range(m + 1):
        v = p[t] - n[t]
        if v:
            j[t] = v
    return j, p, n


def check_G3(p, n, D0, m):
    for t in range(2, m + 3):
        nt1 = n[t - 1] if t - 1 < len(n) else 0
        nt2 = n[t - 2] if t - 2 < len(n) else 0
        pt1 = p[t - 1] if t - 1 < len(p) else 0
        pt2 = p[t - 2] if t - 2 < len(p) else 0
        if 2 * nt1 + nt2 > 2 * comb(D0 - 1, t):
            return False
        if 2 * pt1 + pt2 > 2 * comb(D0 - 1, t - 1):
            return False
    return True


def trial_ENC():
    D0 = rng.randrange(40, 300)
    m = rng.randrange(3, min(D0 - 3, 2 * D0 // 3))
    delta = rng.randrange(2)
    Dp = D0 + delta
    j, p, n = gen_cone_state(D0, m)
    assert check_G3(p, n, D0, m) and n[0] <= D0 - 1 and p[0] == 0
    jn = clamp_row(j, D0, delta)
    tmax = m + 3
    pn, nn = parts(jn, tmax)
    # (a) support frozen
    if any(t > m for t in jn):
        return False
    # (b) cellwise non-increase
    if pn[0] != 0 or nn[0] != max(0, n[0] - 2):
        return False
    for t in range(1, m + 1):
        if nn[t] > n[t] or pn[t] > p[t]:
            return False
    # (c) G3 propagation at frozen capref
    if not check_G3(pn, nn, D0, m):
        return False
    # (d) clock decays
    dec1n = max(0, 2 * (Dp - 1) - (1 + delta) * n[0])
    if nn[1] > max(0, n[1] - dec1n):
        return False
    if pn[1] > max(0, p[1] - 2):
        return False
    for t in range(2, m + 1):
        if nn[t] > max(0, n[t] - (2 * comb(Dp - 1, t) - 2 * comb(D0 - 1, t))):
            return False
        if pn[t] > max(0, p[t] -
                       (2 * comb(Dp - 1, t - 1) - 2 * comb(D0 - 1, t - 1))):
            return False
    return True


def trial_ENG():
    """EN-G relaxed growth law: for ANY all-negative state, one row:
    a'_t <= a_t + [spill_t - 2C(d'-1,t)]^+ (t >= 1); a'_0 = max(0,a_0-2)."""
    d = rng.randrange(30, 300)
    delta = rng.randrange(2)
    dp = d + delta
    m = rng.randrange(3, max(4, d // 2))
    j = {}
    for t in range(m + 1):
        if t == 0 or rng.random() < 0.85:
            mag = ev(rng.randrange(0, 6 * comb(d, min(t + 1, d)) + 2))
            if mag:
                j[t] = -mag
    jn = clamp_row(j, d, delta)
    a = {t: -v for t, v in j.items()}
    if any(v > 0 for v in jn.values()):
        return False
    an = {t: -v for t, v in jn.items()}
    if an.get(0, 0) != max(0, a.get(0, 0) - 2):
        return False
    for t in range(1, m + 3):
        spill = 2 * a.get(t - 1, 0) + a.get(t - 2, 0)
        bound = a.get(t, 0) + max(0, spill - 2 * comb(dp - 1, t))
        if an.get(t, 0) > bound:
            return False
    return True


def trial_ENDR():
    """EN-DR drift: C(D0-1+s, t)*DK >= C(D0-1, t)*(DK + t*s) exactly."""
    D0 = rng.randrange(40, 400)
    t = rng.randrange(1, D0 // 3)
    s = rng.randrange(0, D0 // 2)
    DK = D0 + s + rng.randrange(0, D0)
    return comb(D0 - 1 + s, t) * DK >= comb(D0 - 1, t) * (DK + t * s)


def main():
    N = 600
    res = {}
    lines = []
    for name, fn in (("EN_D_decoupling", trial_END),
                     ("EN_M_monotone", trial_ENM),
                     ("EN_C_mixedcone_step", trial_ENC),
                     ("EN_G_growthlaw", trial_ENG),
                     ("EN_DR_capdrift", trial_ENDR)):
        ok = sum(1 for _ in range(N) if fn())
        res[name] = {"trials": N, "pass": ok, "all_pass": ok == N}
        line = f"[{'PASS' if ok == N else 'FAIL'}] {name}: {ok}/{N}"
        lines.append(line)
        print(line, flush=True)
    res["ALL"] = all(v["all_pass"] for v in res.values() if isinstance(v, dict))
    lines.append(f"OVERALL: {'ALL-PASS' if res['ALL'] else 'FAIL'}")
    print(lines[-1], flush=True)
    json.dump(res, open(os.path.join(
        RESULTS, "amm12592_entry_mixedcone_referee_boxeph.json"), "w"),
        indent=1)
    with open(os.path.join(
            RESULTS, "amm12592_entry_mixedcone_referee_boxeph.out"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
