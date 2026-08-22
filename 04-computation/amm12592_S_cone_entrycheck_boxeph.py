"""LANE E1 -- sharp ENTRY* certificate checker on saved feed-end snapshots.

Reads amm12592_S_cone_feedend_R*_D0*_boxeph.json and evaluates, all exact:

  H1: junk all-negative (and support in [0, m], m + 2 < d);
  H2: abar_t <= (Lam-1) * 2C(d-1,t) for all t >= 1   (Lam = 2^11),
      where abar_1 = a_1 + G_L (layer slack), abar_t = a_t otherwise,
      G_L = sum_k 2*max(0, a_0 - 2k - d)  (exact finite sum);
  H3: a_0 <= d + C0 (C0 = 64) and a_0 <= R - 2;
  H4*: for all t in [2, m+2]:  2*(2*abar_{t-1} + abar_{t-2}) <= 2C(d-1,t)
      (the self-propagating half-cap condition; L3-certified);
  BUD: i_pf <= (R-2)/2  and  i_pf + max(ceil(a0/2), L + max(K1R, K2R))
       <= R - 2, with
       K2R = max_t ceil(abar_t / C(d-1,t))            (cells >= 2),
       K1R = ceil((n+1)/(1-g_hi)), n = ceil(abar_1/(d - C0 - 2))  (cell 1),
       L   = layer length = max(0, ceil((a0 - (d-1))/2)).

Prints a margins table; writes amm12592_S_cone_entrycheck_boxeph.{out,json}.
Margins reported as bit-length differences (display only; decisions exact).
"""
import json, os, glob, re
from fractions import Fraction
from math import comb

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "05-knowledge", "results")
LAM = 1 << 11
C0 = 64
G_HI = Fraction(156759, 262144)     # certified upper sandwich for gamma*


def check(path):
    S = json.load(open(path))
    R, D0, i_pf, d = S["R"], S["D0"], S["i_postfeed"], S["d_feedend"]
    j = {int(t): int(v) for t, v in S["junk"].items()}
    a = {t: -v for t, v in j.items()}
    m = max(a) if a else -1
    a0 = a.get(0, 0)
    H1 = all(v > 0 for v in a.values()) and m + 2 < d
    # layer slack
    G_L = 0
    k = 0
    while a0 - 2 * k - d > 0:
        G_L += 2 * (a0 - 2 * k - d)
        k += 1
    L = max(0, -((-(a0 - (d - 1))) // 2)) if a0 > d - 1 else 0
    abar = dict(a)
    abar[1] = a.get(1, 0) + G_L
    H3 = a0 <= d + C0 and a0 <= R - 2
    # H2 + K2R
    H2 = True
    K2R = 0
    worstA = (None, None)
    for t in range(1, m + 1):
        c = comb(d - 1, t)
        v = abar.get(t, 0)
        if v > (LAM - 1) * 2 * c:
            H2 = False
        if v:
            K2R = max(K2R, -((-v) // c))
            rb = v.bit_length() - (2 * c).bit_length()
            if worstA[0] is None or rb > worstA[0]:
                worstA = (rb, t)
    # H4*
    H4 = True
    worstC = (None, None)
    for t in range(2, m + 3):
        lhs = 2 * (2 * abar.get(t - 1, 0) + abar.get(t - 2, 0))
        rhs = 2 * comb(d - 1, t)
        if lhs > rhs:
            H4 = False
        if lhs:
            rb = lhs.bit_length() - rhs.bit_length()
            if worstC[0] is None or rb > worstC[0]:
                worstC = (rb, t)
    # budgets
    n1 = -((-abar[1]) // (d - C0 - 2)) if abar[1] > 0 else 0
    K1R = int(-((-(n1 + 1)) // (1 - G_HI)))
    ext = L + max(K1R, K2R)
    drain = -((-a0) // 2)
    BUD = (2 * i_pf <= R - 2) and (i_pf + max(drain, ext) <= R - 2)
    PASS = H1 and H2 and H3 and H4 and BUD
    return {"R": R, "D0": D0, "i_pf": i_pf, "d_fe": d, "m": m, "a0": a0,
            "a0_minus_d": a0 - d, "G_L": G_L, "L": L,
            "H1": H1, "H2": H2, "H3": H3, "H4star": H4, "BUD": BUD,
            "worstA_bits_t": worstA, "worstH4_bits_t": worstC,
            "K1R": K1R, "K2R": K2R,
            "thm_capture_ub": i_pf + max(drain, ext),
            "PASS": PASS}


def main():
    rows = []
    for p in sorted(glob.glob(os.path.join(
            RESULTS, "amm12592_S_cone_feedend_R*_boxeph.json")),
            key=lambda q: (json.load(open(q))["R"], json.load(open(q))["D0"])):
        rows.append(check(p))
    hdr = (f"{'R':>6} {'D0':>5} {'i_pf':>5} {'d_fe':>6} {'m':>4} "
           f"{'a0-d':>5} {'L':>2} {'H1':>3} {'H2':>3} {'H3':>3} {'H4*':>4} "
           f"{'BUD':>4} {'H4bits@t':>9} {'K2R':>5} {'capUB':>6} {'PASS':>5}")
    lines = [hdr]
    for r in rows:
        lines.append(
            f"{r['R']:>6} {r['D0']:>5} {r['i_pf']:>5} {r['d_fe']:>6} "
            f"{r['m']:>4} {r['a0_minus_d']:>5} {r['L']:>2} "
            f"{str(r['H1'])[0]:>3} {str(r['H2'])[0]:>3} {str(r['H3'])[0]:>3} "
            f"{str(r['H4star'])[0]:>4} {str(r['BUD'])[0]:>4} "
            f"{str(r['worstH4_bits_t'][0])+'@'+str(r['worstH4_bits_t'][1]):>9} "
            f"{r['K2R']:>5} {r['thm_capture_ub']:>6} {str(r['PASS']):>5}")
    out = "\n".join(lines)
    print(out)
    open(os.path.join(RESULTS, "amm12592_S_cone_entrycheck_boxeph.out"),
         "w").write(out + "\n")
    json.dump(rows, open(os.path.join(
        RESULTS, "amm12592_S_cone_entrycheck_boxeph.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
