#!/usr/bin/env python3
"""Part 1: last-digit (mod 10) and mod-3 transitions of consecutive primes up to 10^11.

Reproduces Lemke Oliver--Soundararajan, "Unexpected biases in the distribution of consecutive
primes", PNAS 113(31) (2016) E4446-E4454 (arXiv:1603.03720v4, read), and compares the data with
their Main Conjecture
    pi(x; q, (a,b)) = li(x)/phi(q)^2 * (1 + c1 * loglog x/log x + c2(q;(a,b))/log x + O((log x)^-7/4)),
    c1(q;(a,b)) = 1/2 - (phi(q)/2) [a = b].
c2 is taken from the paper: (1.1) for a = b; the explicit q = 5 formula of Section 5 (last digit:
primes > 5 mod 10 <-> mod 5) and the q = 3 case of (2.23) for a != b.

Needs the C helper procgen_repunit_20260924_primes.c (segmented sieve, ~90 s CPU for 10^11, < 2 MB).
Checks raise on failure.  Session collatz-procgen-20260922, lane procgen_repunit_20260924.
"""
import os, sys, subprocess, tempfile, math, time
from math import gcd, log
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix="procgen_repunit_primes_")
BIN = os.path.join(TMP, "primes")
XMAX = int(float(os.environ.get("PR_PRIMES_XMAX", "1e11")))
mp.mp.dps = 30


def parse(fn):
    snaps = []
    cur = None
    for line in open(fn):
        t = line.split()
        if t[0] == "SNAP":
            cur = {"kind": t[1], "x": float(t[2]), "pi": int(t[3]), "np": int(t[4]),
                   "plast": int(t[5]), "pnext": int(t[6]), "T": np.zeros((120, 120), dtype=np.int64)}
        elif t[0] == "T":
            cur["T"][int(t[1]), int(t[2])] = int(t[3])
        elif t[0] == "END":
            snaps.append(cur)
        elif t[0] == "DONE":
            print("   ", line.strip())
    return snaps


def agg(T, q):
    res = [r for r in range(q) if gcd(r, q) == 1]
    idx = {r: [i for i in range(120) if i % q == r] for r in res}
    return {(a, b): int(T[np.ix_(idx[a], idx[b])].sum()) for a in res for b in res}


def phi(q):
    return sum(1 for r in range(1, q + 1) if gcd(r, q) == 1)


def li2(x):
    return mp.li(x) - mp.li(2)


# ------------------------------------------------------------------ LO-S constants
def A_q_chi(q, chi, P=10**7):
    """A_{q,chi} = prod_{p|q}(1 - chi(p)/p) prod_{p not | q}(1 - (1-chi(p))^2/(p-1)^2); chi: function on ints."""
    s = np.ones(P + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(P ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    primes = np.nonzero(s)[0]
    tot = mp.mpc(1)
    logsum = 0j
    for p in primes:
        p = int(p)
        c = chi(p)
        if q % p == 0:
            tot *= (1 - mp.mpc(c) / p)
        else:
            f = 1 - (1 - complex(c)) ** 2 / (p - 1) ** 2
            if p < 1000:
                tot *= mp.mpc(f)
            else:
                logsum += np.log(f)
    return complex(tot) * np.exp(logsum)


def los_constants():
    out = {}
    # q = 3
    chi3 = lambda n: [0, 1, -1][n % 3]
    L0 = -sum(a * chi3(a) for a in range(1, 3)) / 3            # L(0, chi_-3) = 1/3
    L1 = complex(mp.dirichlet(1, [0, 1, -1]))                  # pi/(3 sqrt 3)
    A3 = A_q_chi(3, chi3)
    C3 = L0 * L1 * A3
    c2 = {}
    for a in (1, 2):
        for b in (1, 2):
            if a == b:
                c2[(a, b)] = (2 / 2) * log(3 / (2 * math.pi)) + log(2 * math.pi) / 2 - (2 / 2) * log(3) / 2
            else:  # q prime: 1/2 log(2pi/q) + q/phi(q) sum_{chi != chi0} C chi(b-a) + (chi(b)-chi(a))/phi(q)
                c2[(a, b)] = 0.5 * log(2 * math.pi / 3) + 1.5 * (C3 * (chi3(b - a) + (chi3(b) - chi3(a)) / 2)).real
    out[3] = c2
    # q = 5 (last digit): complex characters, chi(2) = i
    tab = {1: 1, 2: 1j, 4: -1, 3: -1j, 0: 0}
    chi5 = lambda n: tab[n % 5]
    L0 = -sum(a * chi5(a) for a in range(1, 5)) / 5
    L1 = complex(mp.dirichlet(1, [0, 1, 1j, -1j, -1]))
    A5 = A_q_chi(5, chi5)
    cb = lambda z: z.conjugate() if isinstance(z, complex) else z
    c2 = {}
    for a in range(1, 5):
        for b in range(1, 5):
            if a == b:
                c2[(a, b)] = 1.5 * log(5 / (2 * math.pi))
            else:
                z = L0 * L1 * A5 * (cb(complex(chi5(b - a))) + (cb(complex(chi5(b))) - cb(complex(chi5(a)))) / 4)
                c2[(a, b)] = 0.5 * log(2 * math.pi / 5) + 2.5 * z.real
    out[5] = c2
    out["A5"], out["L1_5"], out["L0_5"], out["A3"], out["C3"] = A5, L1, L0, A3, C3
    # A_{12,chi} for chi mod 3 (the paper prints A_{12,chi} ~ 1.036): p | 12 factors (1 - chi(p)/p)
    out["A12"] = A_q_chi(12, chi3)
    return out


def main():
    t0 = time.time()
    print("=" * 100)
    print("PART 1. Consecutive primes: last digit (mod 10) and mod 3 transitions (Lemke Oliver--Soundararajan)")
    print("=" * 100)
    subprocess.run(["cc", "-O3", "-o", BIN, os.path.join(HERE, "procgen_repunit_20260924_primes.c"), "-lm"], check=True)
    fn = os.path.join(TMP, "primes.txt")
    t1 = time.time()
    subprocess.run([BIN, str(XMAX), fn], check=True)
    print(f"    sieve to {XMAX:.0e}: {time.time() - t1:.0f} s wall")
    S = parse(fn)
    X = {s["x"]: s for s in S if s["kind"] == "X"}
    N = {int(s["x"]): s for s in S if s["kind"] == "N"}

    # ---------------------------------------------------------------- reproductions
    print("\n(1a) Exact reproductions of the published counts")
    pis = {1e6: 78498, 1e7: 664579, 1e8: 5761455, 1e9: 50847534, 1e10: 455052511, 1e11: 4118054813}
    for x, v in pis.items():
        if x in X:
            assert X[x]["pi"] == v, (x, X[x]["pi"])
    print("    pi(10^6..10^11) = " + ", ".join(f"{X[x]['pi']:,}" for x in pis if x in X) + "  (standard values; checked)")
    los3 = {(1, 1): 215873, (1, 2): 283957, (2, 1): 283957, (2, 2): 216213}
    m3 = agg(N[1000002]["T"], 3)
    print(f"    first 10^6 primes > 3 (pairs n = 3..10^6+2), mod 3: {m3}")
    assert m3 == los3
    print("      = LO-S p.2 table exactly.")
    los10 = {(1, 1): 4623042, (1, 3): 7429438, (1, 7): 7504612, (1, 9): 5442345, (3, 1): 6010982, (3, 3): 4442562,
             (3, 7): 7043695, (3, 9): 7502896, (7, 1): 6373981, (7, 3): 6755195, (7, 7): 4439355, (7, 9): 7431870,
             (9, 1): 7991431, (9, 3): 6372941, (9, 7): 6012739, (9, 9): 4622916}
    if 100000004 in N:
        m10 = agg(N[100000004]["T"], 10)
        m10[(7, 1)] -= 1      # remove the pair (p_4, p_5) = (7, 11): the published window is n = 5..10^8+4
        print(f"    first 10^8 primes >= 11 (pairs n = 5..10^8+4), mod 10: equal to LO-S p.2 table: {m10 == los10}")
        assert m10 == los10
        m10b = agg(N[100000003]["T"], 10)
        diff = {k: m10b[k] - los10[k] for k in los10 if m10b[k] != los10[k]}
        print(f"      (with the window n = 4..10^8+3, i.e. first 10^8 primes > 5, two cells differ by one: {diff})")
    # 4-significant-digit tables of Section 5 (Actual rows)
    pub = {  # (q, x): {(a,b): value}
        (3, 1e9): {(1, 1): 1.132e7, (1, 2): 1.411e7}, (3, 1e10): {(1, 1): 1.024e8, (1, 2): 1.251e8},
        (3, 1e11): {(1, 1): 9.347e8, (1, 2): 1.124e9},
        (4, 1e9): {(1, 1): 1.141e7, (1, 3): 1.401e7}, (4, 1e10): {(1, 1): 1.032e8, (1, 3): 1.244e8},
        (4, 1e11): {(1, 1): 9.412e8, (1, 3): 1.118e9},
        (8, 1e9): {(1, 1): 2.356e6, (1, 3): 3.496e6, (1, 5): 3.351e6, (1, 7): 3.508e6},
        (8, 1e10): {(1, 1): 2.170e7, (1, 3): 3.101e7, (1, 5): 2.988e7, (1, 7): 3.117e7},
        (8, 1e11): {(1, 1): 2.010e8, (1, 3): 2.787e8, (1, 5): 2.696e8, (1, 7): 2.802e8},
        (5, 1e9): {(1, 1): 2.328e6, (1, 2): 3.842e6, (1, 3): 3.796e6, (1, 4): 2.745e6, (2, 1): 3.244e6,
                   (2, 2): 2.228e6, (2, 3): 3.444e6, (3, 1): 3.047e6, (3, 2): 3.595e6, (4, 1): 4.092e6},
    }
    nok = 0
    for (q, x), d in pub.items():
        if x not in X:
            continue
        M = agg(X[x]["T"], q)
        for k, v in d.items():
            got = M[k]
            ok = abs(got - v) <= 0.5 * 10 ** (math.floor(math.log10(v)) - 3) + 1
            assert ok, (q, x, k, got, v)
            nok += 1
    print(f"    {nok} four-significant-digit 'Actual' entries of LO-S Section 5 (q = 3, 4, 8 at 10^9..10^11;"
          f" q = 5 at 10^9): all reproduced.")
    if 1e11 in X:
        s = X[1e11]
        M8 = agg(s["T"], 8)
        bnd = (s["plast"] % 8, s["pnext"] % 8)
        exact8 = {(1, 3): 278676326, (3, 5): 278696997, (5, 7): 278692843, (7, 1): 278681776}
        both = {k: M8[k] - (1 if k == bnd else 0) for k in exact8}
        print(f"    mod 8 at 10^11 (Conjecture 1.6 example): ours (p_n <= x) = { {k: M8[k] for k in exact8} }")
        print(f"      boundary pair ({s['plast']}, {s['pnext']}) = {bnd} mod 8; counting pairs with both primes <= x"
              f" gives {both} = LO-S exactly: {both == exact8}")
        assert both == exact8

    # ---------------------------------------------------------------- constants
    print("\n(1b) Main Conjecture constants (formulas as printed in arXiv:1603.03720v4)")
    K = los_constants()
    print(f"    q = 3: c1 = -1/2 (a=b), +1/2 (a!=b); c2 = {{(a,a): {K[3][(1,1)]:+.6f}, (a,b): {K[3][(1,2)]:+.6f}}}"
          f"  [= -/+ (1/2) log(2 pi/3); C_(3,chi) = L(0)L(1)A = {K['C3'].real:.6f}]")
    assert abs(K[3][(1, 2)] - 0.5 * log(2 * math.pi / 3)) < 1e-12
    print(f"    q = 5: L(0,chi) = {K['L0_5']:.6f}, L(1,chi) = {K['L1_5']:.6f}, A_(5,chi) = {K['A5']:.6f}")
    dig = {1: 1, 2: 7, 3: 3, 4: 9}
    print("    c2(5; (a,b)) as last digits (rows a, cols b = 1 3 7 9):")
    for a in (1, 3, 7, 9):
        ra = {1: 1, 3: 3, 7: 2, 9: 4}[a]
        print(f"      {a}: " + " ".join(f"{K[5][(ra, {1: 1, 3: 3, 7: 2, 9: 4}[b])]:+.4f}" for b in (1, 3, 7, 9)))
    for a in range(1, 5):
        rs = sum(K[5][(a, b)] for b in range(1, 5))
        cs = sum(K[5][(b, a)] for b in range(1, 5))
        assert abs(rs) < 1e-6 and abs(cs) < 1e-6, (a, rs, cs)
    print("    consistency: every row and column of c2(5;.) sums to 0 (<1e-6), as pi(x;q,a) ~ li(x)/phi(q) requires;")
    for (a, b) in ((1, 2), (1, 3), (2, 4)):
        lhs = K[5][(a, b)] + K[5][(b, a)]
        d = (b - a) % 5
        rhs = log(2 * math.pi) - 4 * log(5) / 4        # (1.2): q/(q,b-a) = 5, Lambda(5)/phi(5)
        assert abs(lhs - rhs) < 1e-6
    print("    (1.2) c2(a,b)+c2(b,a) = log(2 pi) - log 5 checked; c2(a,b) = c2(-b,-a) checked:",
          all(abs(K[5][(a, b)] - K[5][((-b) % 5, (-a) % 5)]) < 1e-9 for a in range(1, 5) for b in range(1, 5)))
    print(f"    A_(12,chi_-3) = {K['A12'].real:.4f} (the paper prints ~ 1.036)")
    assert abs(K["A12"].real - 1.036) < 0.001
    # (5.1) values printed by the paper for q = 3, 4
    p51 = {(3, 1e9): (1.156e7, 1.387e7), (3, 1e10): (1.042e8, 1.233e8), (3, 1e11): (9.488e8, 1.110e9),
           (3, 1e12): (8.712e9, 1.009e10), (4, 1e9): (1.164e7, 1.378e7), (4, 1e10): (1.049e8, 1.226e8),
           (4, 1e11): (9.547e8, 1.104e9), (4, 1e12): (8.760e9, 1.004e10)}
    for (q, x), (vd, vo) in p51.items():
        L = mp.log(x)
        f = mp.log(2 * mp.pi * L / q) / (2 * L)
        pd, po = li2(x) / 4 * (1 - f), li2(x) / 4 * (1 + f)
        assert abs(float(pd) / vd - 1) < 1e-3 and abs(float(po) / vo - 1) < 1e-3, (q, x, float(pd), float(po))
    print("    the paper's (5.1) predictions for q = 3, 4 at 10^9..10^12 are reproduced to 4 digits (li from 2).")

    # ---------------------------------------------------------------- matrices and decay
    print("\n(1c) Transition matrices P(next = b | current = a) and the decay of the bias")
    xs = sorted(x for x in X if x >= 1e6 - 1)
    for q, lab in ((10, "last digit"), (3, "mod 3")):
        res = [r for r in range(q) if gcd(r, q) == 1]
        for x in [v for v in xs if abs(math.log10(v) - round(math.log10(v))) < 1e-9]:
            M = agg(X[x]["T"], q)
            print(f"    x = 10^{round(math.log10(x))}, {lab}:")
            for a in res:
                n = sum(M[(a, b)] for b in res)
                print(f"      {a}: " + " ".join(f"{M[(a, b)] / n:.4f}" for b in res) + f"   (n = {n:,})")
    # predicted versus actual relative deviations delta = phi^2 N / li(x) - 1
    print("\n    delta(a,b;x) = phi(q)^2 pi(x;q,(a,b))/li(x) - 1: actual vs Main Conjecture c1*LL/L + c2/L")
    for q, qc in ((10, 5), (3, 3)):
        res = [r for r in range(q) if gcd(r, q) == 1]
        ph = phi(q)
        tomod = (lambda r: r % 5) if q == 10 else (lambda r: r)
        print(f"    q = {q}:  x        LL/L     diag: actual  predicted   off-diag: actual  predicted   max|act-pred|")
        fitrows = []
        for x in xs:
            M = agg(X[x]["T"], q)
            L, LL = math.log(x), math.log(math.log(x))
            li = float(li2(x))
            act = {k: ph * ph * v / li - 1 for k, v in M.items()}
            pred = {}
            for (a, b) in M:
                c1 = 0.5 - (ph / 2) * (a == b)
                pred[(a, b)] = c1 * LL / L + K[qc][(tomod(a), tomod(b))] / L
            da = np.mean([act[(a, a)] for a in res])
            dp = np.mean([pred[(a, a)] for a in res])
            oa = np.mean([act[(a, b)] for a in res for b in res if a != b])
            op = np.mean([pred[(a, b)] for a in res for b in res if a != b])
            mx = max(abs(act[k] - pred[k]) for k in M)
            fitrows.append((x, L, LL, act, pred))
            print(f"      {x:9.3g}  {LL / L:.4f}   {da:+.4f}   {dp:+.4f}        {oa:+.4f}   {op:+.4f}       {mx:.4f}")
        # how much of the deviation the two-term Main Conjecture explains, and the size of what is left
        print(f"    q = {q}: ratio actual/predicted (diagonal mean, off-diagonal mean) and residual/(LL/L)^2:")
        line = []
        for x, L, LL, act, pred in fitrows[::2]:
            da = np.mean([act[(a, a)] for a in res]); dp = np.mean([pred[(a, a)] for a in res])
            oa = np.mean([act[(a, b)] for a in res for b in res if a != b])
            op = np.mean([pred[(a, b)] for a in res for b in res if a != b])
            line.append(f"{x:.2g}: {da / dp:.3f}, {oa / op:.3f}, {(da - dp) / (LL / L) ** 2:+.2f}")
        print("      " + "; ".join(line))
        # the per-cell ordering: does the sign/ordering of c2 predict the ordering of the off-diagonal cells?
        x, L, LL, act, pred = fitrows[-1]
        off = [(a, b) for a in res for b in res if a != b]
        ra = sorted(off, key=lambda k: act[k]); rp = sorted(off, key=lambda k: pred[k])
        from scipy.stats import spearmanr
        rho = spearmanr([act[k] for k in off], [pred[k] for k in off]).correlation if len(off) > 2 else float("nan")
        print(f"    q = {q} at x = {x:.2g}: Spearman correlation of actual vs predicted off-diagonal cells = {rho:.3f};"
              f" most favoured actual {ra[-1]}, predicted {rp[-1]}; least favoured actual {ra[0]}, predicted {rp[0]}")
        # decay of the row TV distance from uniform
        print(f"    q = {q}: max row total-variation distance from the uniform row, by x:")
        line = []
        for x, L, LL, act, pred in fitrows:
            M = agg(X[x]["T"], q)
            tv = 0
            for a in res:
                n = sum(M[(a, b)] for b in res)
                tv = max(tv, 0.5 * sum(abs(M[(a, b)] / n - 1 / ph) for b in res))
            line.append(f"{x:.2g}: {tv:.4f}")
        print("      " + ", ".join(line))
        # second eigenvalue of the empirical transition matrix
        line = []
        for x, L, LL, act, pred in fitrows[::4]:
            M = agg(X[x]["T"], q)
            P = np.array([[M[(a, b)] / sum(M[(a, c)] for c in res) for b in res] for a in res])
            ev = sorted(np.abs(np.linalg.eigvals(P)), reverse=True)
            line.append(f"{x:.2g}: {ev[1]:.4f}")
        print(f"    q = {q}: |lambda_2| of the empirical matrix: " + ", ".join(line))
        # the stationary law (the Brouwer/Perron fixed point of the empirical matrix) at the largest x
        x, L, LL, act, pred = fitrows[-1]
        M = agg(X[x]["T"], q)
        P = np.array([[M[(a, b)] / sum(M[(a, c)] for c in res) for b in res] for a in res])
        w, V = np.linalg.eig(P.T)
        v = np.real(V[:, np.argmin(np.abs(w - 1))]); v = v / v.sum()
        print(f"    q = {q} at x = {x:.2g}: stationary law of the empirical matrix = "
              + ", ".join(f"{a}: {vv:.6f}" for a, vv in zip(res, v)) + f"  (uniform = {1 / ph:.6f})")
    print(f"\nPart 1 done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
