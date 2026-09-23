#!/usr/bin/env python3
"""HYP-9129 numerics (procgen 2026-09-23).  Everything in this script is NUMERICAL; nothing here is a proof.

Asymptotic model (see amm12592_procgen_20260923_hyp9128_proof.md, section 6).
A ratio-B super-block [N, BN) with handoff state S_N in Z[w] is given by:
  * S_N(0) = +-1 and deg S_N <= (B-2)N/2;
  * the fold target F_0(u) = ((1+u)/2)^(N-1) S_N((1-u^2)/4).
Its zero-counting measure mu (divided by N) has mass <= M = (B-2)/2. The Mahler constraint int log|zeta| dmu <= 0
follows from S(0) = +-1 and an integral leading coefficient. Then (1/N) log|S_N(w)| ~ U_mu(w) = int log|1 - w/zeta| dmu.

Modes:
 * family: the fold-sufficient threshold of the golden-zero families mu = kappa delta_{-1} + (M-kappa) Unif(|w|=1).
   This is the measure-level version of the explicit HYP-9128 family (kappa = 1/16). "Fold-sufficient" is the
   criterion of the HYP-9128 proof, taken at exponential order with circles only:
     [V1] levels v in [0, B/c-1], sigma = min((c-1)(1+v), B-1-v),
          min over circles enclosing 0,1 of max(log|(1+u)/2| + U_mu(w) - sigma log|u| - v log|(u-1)/2|) <= 0;
     [V2M] level-0 middle regime with the pairing factor g1^min(tau, s0-tau);
     [V2T] level-0 top regime, on circles enclosing -1, 0, 1;
     [Bm] level-0 bottom regime, a Cauchy bound for [z^r] S(w(z)) (1+z)^{-(2-c)N} against log C(R_0, r) at r = rho N;
     [B1] r = 1: |s' - (2-c)| <= 0.95 (c-1), where s' = int Re(-1/zeta) dmu.
 * lp --dict atoms|cyclo|factor: H2, minimising the fold-sufficient threshold over mu by alternating LPs. The
   extremal circles are fixed, an LP is solved in the dictionary weights, and the circles are re-selected. The
   bracket is bisected. The dictionaries are:
     atoms: conjugate-pair atoms on a polar grid, plus Unif(|w|=1). This is the real-measure relaxation and is
       not realizable.
     cyclo: the 2-power cyclotomic factors 1 + w^(2^(k-1)), k <= 6, plus Unif(|w|=1). These are realizable:
       every mixture is the limit of prod_k (1 + w^(2^(k-1)))^(e_k), which lies in Z[w], has S(0) = 1 and
       S == 1 + w^N (mod 2).
     factor: integer factors P(w) of degree <= 4 with P(0) = 1 and P == (1+w)^deg (mod 2), plus Unif(|w|=1),
       with the rest completed by 2-power cyclotomics. These are realizable, and Mahler is automatic.
 * h4: the necessary LP. It uses the single-block evaluation inequality (Long's Thm 11.1 generalised to
   palindromic P_+ = S(w)(1+x)^((B-2)N)):
     U_mu(w) <= min over both preimages x of w of
                max_v [v log|x| + alpha(v) log|1+x| + rho(v) log(1+|x|)] - (B-2) log|1+x|.
   It is imposed on a test grid, with an atom grid for mu. Zeros placed exactly at test points are not excluded,
   so the LP is a heuristic and not a certificate.
Usage:  python3 amm12592_procgen_20260923_hyp9128_h2.py family --B 4
        python3 amm12592_procgen_20260923_hyp9128_h2.py lp --dict atoms --B 4
        python3 amm12592_procgen_20260923_hyp9128_h2.py h4 --B 2 4 8
"""
from __future__ import annotations
import argparse, itertools, math, sys, time
import numpy as np
from scipy.optimize import linprog

EPS = 1e-300
ETA = 3 / 128


def Hent(x):
    return -x * math.log(x) - (1 - x) * math.log(1 - x)


def circles(encl, TH):
    X0 = np.linspace(-0.8, 1.8, 27)
    RHO = np.exp(np.linspace(np.log(0.05), np.log(12), 44))
    C = [(x0, r) for x0 in X0 for r in RHO if abs(x0) < 0.97 * r and (1 not in encl or abs(1 - x0) < 0.9 * r)
         and (-1 not in encl or abs(-1 - x0) < 0.9 * r)]
    U = np.array([x0 + r * np.exp(1j * TH) for x0, r in C])
    return dict(W=(1 - U ** 2) / 4, LA=np.log(np.abs((1 + U) / 2) + EPS), LL=np.log(np.abs(U) + EPS),
                LM=np.log(np.abs((U - 1) / 2) + EPS))


class Dictionary:
    """Basis measures of unit mass. Each is uniform on a finite root set, or 'U1' = Unif(|w| = 1)."""

    def __init__(self, items, labels):
        self.items, self.labels = items, labels
        self.K = len(items)
        self.mahler = np.array([0.0 if isinstance(it, str) else float(np.mean(np.log(np.abs(it)))) for it in items])
        self.sprime = np.array([0.0 if isinstance(it, str) else float(np.mean(np.real(-1 / it))) for it in items])

    def U_one(self, j, W):
        it = self.items[j]
        if isinstance(it, str):
            return np.maximum(0.0, np.log(np.abs(W) + EPS))
        return np.mean([np.log(np.abs(1 - W / z) + EPS) for z in it], axis=0)

    def U_mix(self, W, x):
        out = np.zeros(W.shape)
        for j in np.nonzero(x > 1e-12)[0]:
            out += x[j] * self.U_one(j, W)
        return out

    def U_all(self, Wrow):
        return np.stack([self.U_one(j, Wrow) for j in range(self.K)], axis=-1)


def atoms_dictionary():
    zr = np.array([0.5, 0.7, 0.85, 0.93, 0.97, 1.0, 1.03, 1.08, 1.18, 1.4, 2.0, 3.0])   # dense near |zeta| = 1
    za = np.linspace(0, math.pi, 23)
    items, labels = [], []
    for r in zr:
        for a in za:
            z = r * np.exp(1j * a)
            items.append(np.array([z, np.conj(z)]) if abs(z.imag) > 1e-12 else np.array([complex(z.real, 0)]))
            labels.append(f"{z:.2f}")
    items.append("U1")
    labels.append("Unif|w|=1")
    return Dictionary(items, labels)


def cyclo_dictionary(kmax=6):
    items, labels = [], []
    for k in range(1, kmax + 1):
        n = 2 ** k
        items.append(np.array([np.exp(2j * math.pi * j / n) for j in range(n) if math.gcd(j, n) == 1]))
        labels.append(f"1+w^{2 ** (k - 1)}")
    items.append("U1")
    labels.append("Unif|w|=1")
    return Dictionary(items, labels)


def factor_dictionary():
    items, labels = [], []

    def add(coefs):
        r = np.roots(coefs[::-1])
        if len(r) == 0 or np.any(np.abs(r) < 0.05):
            return
        items.append(r)
        labels.append("+".join(f"{c}w^{k}" for k, c in enumerate(coefs) if c))
    for a in [-2, -1, 0, 1]:
        add([1, 2 * a + 1])
    for a in range(-4, 5):
        for b in [-9, -7, -5, -3, -1, 1, 3, 5, 7, 9]:
            add([1, 2 * a, b])
    for a, b, c_ in itertools.product(range(-2, 3), repeat=3):
        for d in [-5, -3, -1, 1, 3, 5]:
            add([1, 2 * a, 2 * b, 2 * c_, d])
    items.append("U1")
    labels.append("Unif|w|=1")
    return Dictionary(items, labels)


class Model:
    def __init__(self, B, ntheta=161):
        self.B, self.M = B, (B - 2) / 2
        TH = np.linspace(0, math.pi, ntheta)
        self.F01, self.F0, self.FT = circles({0, 1}, TH), circles({0}, TH), circles({-1, 0, 1}, TH)
        SZ = np.exp(np.linspace(np.log(0.002), np.log(0.9), 40))
        ZZ = np.array([s * np.exp(2j * TH) for s in SZ])
        self.WZ = ZZ / (1 + ZZ) ** 2
        self.L1Z = np.log(np.abs(1 + ZZ) + EPS)
        self.LSZ = np.log(SZ)[:, None] * np.ones(ZZ.shape[1])[None, :]

    def select(self, c, Ufun):
        """best circle per condition for the current measure; returns [(W_row, rhs_row)]: need U_mu(W_row) <= rhs_row."""
        B, s0 = self.B, c - 1
        g1 = 1 - 2 * (ETA / s0) * (1 - ETA / s0)
        rows = []
        F = self.F01
        Um = Ufun(F["W"])
        for v in np.linspace(0, B / c - 1, 26):
            sig = min(s0 * (1 + v), B - 1 - v)
            rhs = -(F["LA"] - sig * F["LL"] - v * F["LM"])
            ci = int(np.argmin(np.max(Um - rhs, axis=1)))
            rows.append((F["W"][ci], rhs[ci]))
        F = self.F0
        Um = Ufun(F["W"])
        for tau in np.linspace(0, s0, 24):
            pair = max(min(tau, s0 - tau), 0) * math.log(g1)
            rhs = -(F["LA"] - tau * F["LL"] + pair)
            ci = int(np.argmin(np.max(Um - rhs, axis=1)))
            rows.append((F["W"][ci], rhs[ci]))
        F = self.FT
        Um = Ufun(F["W"])
        rhs = -(F["LA"] - s0 * F["LL"])
        ci = int(np.argmin(np.max(Um - rhs, axis=1)))
        rows.append((F["W"][ci], rhs[ci]))
        UmZ = Ufun(self.WZ)
        for rho in [0.004, 0.008, 0.012, 0.016, 0.020, ETA]:
            rhs = (2 - c) * self.L1Z + rho * self.LSZ + s0 * Hent(rho / s0) - 1e-9
            si = int(np.argmin(np.max(UmZ - rhs, axis=1)))
            rows.append((self.WZ[si], rhs[si]))
        return rows

    def check(self, c, Ufun, sprime, mahler):
        rows = self.select(c, Ufun)
        worst = max(float(np.max(Ufun(W) - h)) for W, h in rows)
        return worst <= 0 and abs(sprime - (2 - c)) <= 0.95 * (c - 1) and mahler <= 1e-9, worst


def bisect(feas, lo, hi, steps, log):
    ok, info = feas(hi)
    log(f"  c = {hi:.4f}: {'feasible' if ok else 'infeasible'} (bracket top)")
    if not ok:
        return None, None
    best = info
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        ok, info = feas(mid)
        log(f"  c = {mid:.4f}: {'feasible' if ok else 'infeasible'}")
        if ok:
            hi, best = mid, info
        else:
            lo = mid
    return hi, best


def run_family(B, kappas, lo, hi, steps, ntheta):
    model = Model(B, ntheta)
    M = model.M
    print(f"=== [family] B = {B}: fold-sufficient threshold of mu = kappa delta_(-1) + (M - kappa) Unif(|w|=1), M = {M} ===")
    for kap in kappas:
        U = lambda W, kap=kap: kap * np.log(np.abs(1 + W) + EPS) + (M - kap) * np.maximum(0, np.log(np.abs(W) + EPS))
        feas = lambda c: (model.check(c, U, kap * 1.0, 0.0)[0], None)
        th, _ = bisect(feas, lo, hi, steps, lambda s: None)
        print(f"  kappa = {kap:.4f}: threshold ~ {th if th is None else round(th, 4)}", flush=True)


def run_lp(B, dname, lo, hi, steps, ntheta, iters=10, keep_every=6, topk=16, full_rows=False):
    model = Model(B, ntheta)
    D = {"atoms": atoms_dictionary, "cyclo": cyclo_dictionary, "factor": factor_dictionary}[dname]()
    M, K = model.M, D.K
    print(f"=== [lp] B = {B}, dictionary '{dname}' ({K} basis measures), theta grid {ntheta}: alternating LP ===", flush=True)
    jU = D.labels.index("Unif|w|=1")
    j1 = [j for j, l in enumerate(D.labels) if l in ("1+w^1", "1w^0+1w^1", "(-1+0j)", "-1.00+0.00j")]
    x0 = np.zeros(K)
    x0[jU] = 0.9375 * M
    if j1:
        x0[j1[0]] = 0.0625 * M
    else:
        x0[jU] = M

    def lp(c, x):
        rows = model.select(c, lambda W: D.U_mix(W, x))
        A, b = [], []
        for W, h in rows:
            cur = D.U_mix(W, x) - h
            keep = np.zeros(len(h), bool)          # cutting-plane rows: a coarse subgrid plus the near-active points
            keep[::keep_every] = True                # (the full grid is re-checked by model.check before accepting)
            keep[np.argsort(cur)[-topk:]] = True
            if full_rows:
                keep[:] = True
            A.append(D.U_all(W[keep]))
            b.append(h[keep])
        A = np.vstack(A)
        b = np.concatenate(b)
        Aub = np.hstack([A, np.ones((A.shape[0], 1))])
        Aub = np.vstack([Aub, np.r_[D.sprime, 0], np.r_[-D.sprime, 0], np.r_[D.mahler, 0]])
        b = np.concatenate([b, [(2 - c) + 0.95 * (c - 1), -((2 - c) - 0.95 * (c - 1)), 0.0]])
        res = linprog(np.r_[np.zeros(K), -1.0], A_ub=Aub, b_ub=b, A_eq=np.r_[np.ones(K), 0][None, :], b_eq=[M],
                      bounds=[(0, None)] * K + [(None, 0.2)], method="highs")
        if res.status != 0:
            return None, x
        return -res.fun, res.x[:-1]

    def feas(c, xs=[x0]):
        x = xs[0].copy()
        for _ in range(iters):
            t, xn = lp(c, x)
            if t is None:
                return False, x
            x = xn
            if t >= 0:
                ok, worst = model.check(c, lambda W: D.U_mix(W, x), float(D.sprime @ x), float(D.mahler @ x))
                if ok:
                    xs[0] = x
                    return True, x
        return False, x

    th, x = bisect(feas, lo, hi, steps, lambda s: print(s, flush=True))
    if th is None:
        print(f"  no feasible point at the top of the bracket {hi}")
        return
    top = np.argsort(-x)[:8]
    print(f"B = {B}, dictionary '{dname}': alternating-LP fold-sufficient threshold ~ {th:.4f}; mass {x.sum():.3f}, "
          f"Mahler {D.mahler @ x:+.4f}, s' = {D.sprime @ x:.3f}; heaviest: "
          + ", ".join(f"{x[j]:.3f}*[{D.labels[j]}]" for j in top if x[j] > 1e-6), flush=True)


def run_h4(Bs):
    def Phi(x, c, B, nv=400):
        vs = np.linspace(0, B - 1, nv)
        rd = np.minimum(c * (1 + vs), B) - (1 + vs)
        un = B - np.minimum(c * (1 + vs), B)
        lx, l1, lp_ = np.log(np.abs(x)), np.log(np.abs(1 + x)), np.log(1 + np.abs(x))
        return np.max(vs[None, :] * lx[:, None] + un[None, :] * l1[:, None] + rd[None, :] * lp_[:, None], axis=1) - (B - 2) * l1

    zr = np.exp(np.linspace(np.log(0.08), np.log(12), 26))
    za = np.linspace(0, math.pi, 25)
    Z = np.array([r * np.exp(1j * a) for r in zr for a in za])
    tr = np.exp(np.linspace(np.log(0.013), np.log(9.7), 40))
    ta = np.linspace(0.011, math.pi - 0.011, 33)
    Wt = np.array([r * np.exp(1j * a) for r in tr for a in ta] + list(-tr) + list(tr * 1.0001))
    A = 0.5 * (np.log(np.abs(1 - Wt[:, None] / Z[None, :])) + np.log(np.abs(1 - Wt[:, None] / np.conj(Z)[None, :])))
    bq = 2 - 1 / Wt
    disc = np.sqrt(bq ** 2 - 4 + 0j)
    x1, x2 = (-bq + disc) / 2, (-bq - disc) / 2
    print("=== [h4] necessary LP: single-block evaluation at both preimages + mass + Mahler (heuristic) ===")
    for B in Bs:
        M, K = (B - 2) / 2, len(Z)

        def slack(c):
            rhs = np.minimum(Phi(x1, c, B), Phi(x2, c, B))
            Aub = np.vstack([np.hstack([A, np.ones((len(Wt), 1))]), np.r_[np.ones(K), 0], np.r_[np.log(np.abs(Z)), 0]])
            res = linprog(np.r_[np.zeros(K), -1.0], A_ub=Aub, b_ub=np.concatenate([rhs, [M, 0.0]]),
                          bounds=[(0, None)] * K + [(None, 1.0)], method="highs")
            return (None, None) if res.status != 0 else (-res.fun, res.x)
        lo, hi = 1.2, 1.7
        for _ in range(16):
            mid = 0.5 * (lo + hi)
            t, x = slack(mid)
            if t is not None and t >= 0:
                hi = mid
            else:
                lo = mid
        t, x = slack(hi)
        top = np.argsort(-x[:-1])[:5]
        print(f"  B = {B}: LP necessary threshold ~ {hi:.4f}; heaviest zeros "
              + ", ".join(f"{x[k]:.3f}@{Z[k]:.2f}" for k in top if x[k] > 1e-6), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["family", "lp", "h4"])
    ap.add_argument("--B", type=int, nargs="+", default=[4])
    ap.add_argument("--dict", default="atoms", choices=["atoms", "cyclo", "factor"])
    ap.add_argument("--lo", type=float, default=None)
    ap.add_argument("--hi", type=float, default=None)
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument("--ntheta", type=int, default=161)
    ap.add_argument("--full-rows", action="store_true", help="LP on the full theta grid (no cutting-plane subsampling)")
    args = ap.parse_args()
    t0 = time.time()
    if args.mode == "family":
        for B in args.B:
            run_family(B, [0.0, 1 / 32, 1 / 16, 0.09, 1 / 8], args.lo or 1.50, args.hi or 1.64, args.steps, args.ntheta)
    elif args.mode == "lp":
        for B in args.B:
            run_lp(B, args.dict, args.lo or 1.40, args.hi or 1.62, args.steps, args.ntheta, full_rows=args.full_rows)
    else:
        run_h4(args.B)
    print(f"[{time.time() - t0:.0f}s]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
