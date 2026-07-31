#!/usr/bin/env python3
"""Lane F2 (G4): policy-independence of the band freeze -- REVISED.

Frame: THM-2966 spine normal form; depth law d_m = floor(g1*m/g2) + D0;
0-side cell (m,k) at monomial (z,o) = (m+d_m-k, k+1), 1-side mirrored;
doubled deficits delta integer, |delta| <= binom(d_m,k), delta == binom
mod 2 (Lucas).  PROVED necessary envelope (lane D):

    |D_m(p)| <= (p^{m+1} + q^{m+1}) / 2   for every row m, every p.  (E_m)

Scaled at bias p = n/dd and level Lmax:  S_i = 2 D b^{Lmax} integer;
Env_m = (n^{m+1} + (dd-n)^{m+1}) dd^{Lmax-m-1}.

New input fact (lane C2/G3): DFS shows the 9-bias envelope constraint set
at (gamma,D0)=(1/2,0) is FEASIBLE to M=10 (>> band birth m*=4), so lane
D's greedy freeze is a POLICY artifact at small M.  This script:

  witness      -- deterministic DFS witness extractor (choice list, not just
                  verdict) + full anatomy: per-row rho profile at 9 biases,
                  binding bias, saturation profile delta/cap, merged-ledger
                  band trajectory (band vs cone value split at p=1/2),
                  divergence-from-greedy, dyadic checkpoint (beta_M == 0 at
                  M+1=2^r, THM-2976 T1) cross-check.
  beam         -- anticipatory policy family: beam search over per-cell
                  transport choices with the PROVED envelope as pruning
                  oracle.  W=1 ~ lookahead-0 rho-greedy; W large ~ deep
                  anticipation; scoring 'max' (min worst-bias rho) or 'sum'
                  (LP-guided weighted-residual Lyapunov).  Survival M and
                  the pinch trajectory best-rho(m) are the outputs.  A beam
                  state that completes row M IS an exact witness for the
                  9-bias envelope set at M (independently re-verified).
  prefixgreedy -- replay a witness prefix of P rows, then continue with
                  lane D's zeroing greedy; measure first envelope-violation
                  row as a function of P (delay-vs-reroute law).
  symdfs       -- mirror-symmetric restricted DFS (delta equal on the two
                  sides of every (m,k)): search space ~ sqrt, exhaustion
                  affordable further out; symmetric infeasibility is NOT
                  general infeasibility but kills the symmetric corridor.

Exact integer arithmetic for every asserted quantity; floats only for
scores/printing.  RULES: no git, outputs under assigned filenames only.
"""

import argparse
import random
import sys
import time
from fractions import Fraction as Fr
from math import comb

sys.setrecursionlimit(200000)

BIASES9 = [(1, 2), (1, 3), (2, 3), (2, 5), (3, 5), (1, 4), (3, 4),
           (1285, 2181), (8847357, 11821757)]


def require(c, msg):
    if not c:
        raise RuntimeError(msg)


# ---------------------------------------------------------------- structures
def build(g1, g2, D0, Mmax, biases):
    d = lambda m: (g1 * m) // g2 + D0
    A = lambda m: m + d(m) + 1
    Lmax = A(Mmax)
    rows = []           # rows[m-1] = list of cells (m, side, k, z, o, cap, par, w)
    for m in range(1, Mmax + 1):
        dm = d(m)
        cells = []
        for k in range(dm + 1):
            cap = comb(dm, k)
            for side, (z, o) in enumerate(((m + dm - k, k + 1),
                                           (k + 1, m + dm - k))):
                w = tuple(n**z * (dd - n)**o * dd**(Lmax - z - o)
                          for (n, dd) in biases)
                cells.append((m, side, k, z, o, cap, cap % 2, w))
        rows.append(cells)
    Env = [tuple((n**(m + 2) + (dd - n)**(m + 2)) * dd**(Lmax - m - 2)
                 for (n, dd) in biases) for m in range(Mmax)]
    return d, A, Lmax, rows, Env


def choices_of(cap, par):
    return [v for v in range(-cap, cap + 1) if (v - par) % 2 == 0]


def verify_assignment(g1, g2, D0, M, biases, choice, tag):
    """Independent exact replay: box, parity, envelope at every row-end."""
    d, A, Lmax, rows, Env = build(g1, g2, D0, M, biases)
    nP = len(biases)
    S = [0] * nP
    idx = 0
    rho_hist = []
    for m0 in range(M):
        for (m, side, k, z, o, cap, par, w) in rows[m0]:
            v = choice[idx]; idx += 1
            require(abs(v) <= cap and (v - par) % 2 == 0,
                    f"{tag}: box/parity fail at ({m},{side},{k})")
            if v:
                for i in range(nP):
                    S[i] += v * w[i]
        rhos = []
        for i in range(nP):
            require(abs(S[i]) <= Env[m0][i],
                    f"{tag}: envelope fail row {m0+1} bias {biases[i]}")
            rhos.append(abs(S[i]) / Env[m0][i] if Env[m0][i] else 0.0)
        rho_hist.append(rhos)
    require(idx == len(choice), f"{tag}: choice length mismatch")
    print(f"[{tag}] RE-VERIFIED EXACTLY: M={M} rows all satisfy the PROVED "
          f"envelope at {nP} biases (box+parity+E_m).", flush=True)
    return rho_hist


# ---------------------------------------------------------------- DFS witness
def dfs_witness(g1, g2, D0, Mtry, biases, node_cap, order="abs", tag="dfs"):
    d, A, Lmax, rows, Env = build(g1, g2, D0, Mtry, biases)
    nP = len(biases)
    R = [tuple(sum(c[5] * c[7][i] for c in rows[m]) for i in range(nP))
         for m in range(Mtry)]
    allowed = []
    for m in range(Mtry):
        best = None
        for mp in range(m, Mtry):
            tot = tuple(Env[mp][i] + sum(R[r][i] for r in range(m + 1, mp + 1))
                        for i in range(nP))
            best = tot if best is None else tuple(min(a, b)
                                                  for a, b in zip(best, tot))
        allowed.append(best)
    suffix = []
    for m in range(Mtry):
        sw = [tuple(0 for _ in range(nP))]
        for c in reversed(rows[m]):
            sw.append(tuple(sw[-1][i] + c[5] * c[7][i] for i in range(nP)))
        sw.reverse()
        suffix.append(sw)

    state = {"nodes": 0}

    def dfs(m, ci, S, hist):
        state["nodes"] += 1
        if state["nodes"] > node_cap:
            raise TimeoutError
        cells = rows[m]
        if ci == len(cells):
            for i in range(nP):
                if abs(S[i]) > Env[m][i] or abs(S[i]) > allowed[m][i]:
                    return None
            if m + 1 == Mtry:
                return list(hist)
            return dfs(m + 1, 0, S, hist)
        rest = suffix[m][ci]
        for i in range(nP):
            if abs(S[i]) > rest[i] + allowed[m][i]:
                return None
        (_, _, _, _, _, cap, par, w) = cells[ci]
        vals = choices_of(cap, par)
        if order == "abs":
            vals.sort(key=abs)
        elif order == "desc":
            vals.sort(key=lambda v: -abs(v))
        for v in vals:
            S2 = tuple(S[i] + v * w[i] for i in range(nP)) if v else S
            hist.append(v)
            r = dfs(m, ci + 1, S2, hist)
            if r is not None:
                return r
            hist.pop()
        return None

    try:
        w = dfs(0, 0, tuple(0 for _ in range(nP)), [])
    except TimeoutError:
        print(f"[{tag}] NODE-CAP {node_cap} at M={Mtry}", flush=True)
        return None, state["nodes"]
    print(f"[{tag}] M={Mtry}: {'WITNESS' if w else 'INFEASIBLE(exhausted)'} "
          f"nodes={state['nodes']}", flush=True)
    return w, state["nodes"]


# ---------------------------------------------------------------- ledger (laneD)
class Ledger:
    """Merged-cell descending ledger (verbatim mechanics from laneD)."""
    def __init__(self, g1, g2, D0):
        self.g1, self.g2, self.D0 = g1, g2, D0
        self.b = {}
        self.A = None
        self.rows = []
    def d(self, m):
        return (self.g1 * m) // self.g2 + self.D0
    def Ad(self, m):
        return m + self.d(m) + 1
    def cells(self, m):
        dm, A = self.d(m), self.Ad(m)
        cel = {}
        for k in range(dm + 1):
            c = comb(dm, k)
            for o in (k + 1, A - 1 - k):
                cap, par = cel.get(o, (0, 0))
                cel[o] = (cap + c, (par + c) % 2)
        return cel
    def descend(self, t):
        if t == 0 or not self.b:
            return
        row = [comb(t, i) for i in range(t + 1)]
        nb = {}
        for o, c in self.b.items():
            for i, r in enumerate(row):
                nb[o + i] = nb.get(o + i, 0) + c * r
        self.b = {o: c for o, c in nb.items() if c}
    def step(self, m, policy):
        Anew = self.Ad(m)
        if self.A is not None:
            self.descend(Anew - self.A)
        self.A = Anew
        rec = []
        for o, (cap, par) in sorted(self.cells(m).items()):
            delta = policy(self, m, o, cap, par)
            require(abs(delta) <= cap and (delta - par) % 2 == 0,
                    "box/parity violation in ledger step")
            nv = self.b.get(o, 0) + delta
            if nv:
                self.b[o] = nv
            else:
                self.b.pop(o, None)
            rec.append((o, cap, par, delta))
        self.rows.append((m, Anew, rec))
    def val(self, p):
        q = 1 - p
        return sum(c * p**(self.A - o) * q**o
                   for o, c in self.b.items()) / 2
    def dval(self, p):
        q = 1 - p
        tot = Fr(0)
        for o, c in self.b.items():
            z = self.A - o
            t = Fr(0)
            if z:
                t += z * p**(z - 1) * q**o
            if o:
                t -= o * p**z * q**(o - 1)
            tot += c * t
        return tot / 2
    def band_split(self, m, p):
        """(band value, cone value) at p after row m has been placed."""
        q = 1 - p
        dm = self.d(m)
        band = Fr(0)
        cone = Fr(0)
        for o, c in self.b.items():
            t = Fr(c) * p**(self.A - o) * q**o / 2
            if dm + 2 <= o <= m - 1:
                band += t
            else:
                cone += t
        return band, cone


def greedy_policy(led, m, o, cap, par):
    t = -led.b.get(o, 0)
    t = max(-cap, min(cap, t))
    if (t - par) % 2:
        b0 = led.b.get(o, 0)
        cand = [x for x in (t - 1, t + 1) if -cap <= x <= cap]
        t = min(cand, key=lambda x: (abs(b0 + x), abs(x)))
    return t


def merged_from_choice(g1, g2, D0, M, choice):
    """(m,o) -> merged doubled deficit, from a per-cell choice list."""
    d = lambda m: (g1 * m) // g2 + D0
    md = {}
    idx = 0
    for m in range(1, M + 1):
        dm = d(m)
        A = m + dm + 1
        for k in range(dm + 1):
            for side, o in enumerate((k + 1, A - 1 - k)):
                md[(m, o)] = md.get((m, o), 0) + choice[idx]
                idx += 1
    require(idx == len(choice), "merged: choice length mismatch")
    return md


def envelope_check(led, biases):
    """Exact E_m check at current last row; returns (ok, rho floats)."""
    m = led.rows[-1][0]
    rhos = []
    ok = True
    for (n, dd) in biases:
        p = Fr(n, dd)
        v = led.val(p)
        env = (p**(m + 1) + (1 - p)**(m + 1)) / 2
        rhos.append(float(abs(v) / env))
        if abs(v) > env:
            ok = False
    return ok, rhos


def dense_check(g1, g2, D0, M, choice, N, tag):
    """Exact E_m check at every row for the dense grid p = j/N (all j),
    on top of the 9 search biases: ledger replay, Fraction arithmetic."""
    md = merged_from_choice(g1, g2, D0, M, choice)
    led = Ledger(g1, g2, D0)
    grid = [Fr(j, N) for j in range(1, N)]
    worst = 0.0
    worst_at = None
    viol = 0
    for m in range(1, M + 1):
        led.step(m, lambda led, mm, o, cap, par: md.get((mm, o), 0))
        rmax, pmax = 0.0, None
        for p in grid:
            v = led.val(p)
            env = (p**(m + 1) + (1 - p)**(m + 1)) / 2
            r = abs(v) / env
            if r > 1:
                viol += 1
            if float(r) > rmax:
                rmax, pmax = float(r), p
        if rmax > worst:
            worst, worst_at = rmax, (m, pmax)
        dv = led.dval(Fr(1, 2))
        dscale = float(dv) * 2.0**(m + 1) / (m + 1)
        print(f"  [dense {tag}] m={m:3d} max_grid_rho={rmax:.4f} at p={pmax} "
              f"D'(1/2)={float(dv):+.3e} [x2^(m+1)/(m+1)={dscale:+.3f}]",
              flush=True)
    print(f"[dense {tag}] grid p=j/{N}: violations={viol}; worst rho="
          f"{worst:.4f} at (m,p)={worst_at}", flush=True)
    return worst, viol


def antisym_defect(g1, g2, D0, M, choice):
    """Fraction of (m,k) pairs with delta_side1 != -delta_side0."""
    d = lambda m: (g1 * m) // g2 + D0
    idx = 0
    bad = tot = 0
    for m in range(1, M + 1):
        for k in range(d(m) + 1):
            v0, v1 = choice[idx], choice[idx + 1]
            idx += 2
            tot += 1
            if v1 != -v0:
                bad += 1
    return bad, tot


# ---------------------------------------------------------------- beta_M clock
def beta_poly(M, dM):
    A = M + dM + 1
    par = {}
    for o in range(A + 1):
        v = 1 if (o & A) == o else 0
        if o <= dM and (o & dM) == o:
            v ^= 1
        j = o - (M + 1)
        if 0 <= j <= dM and (j & dM) == j:
            v ^= 1
        if v:
            par[o] = 1
    return par


def clock_report(g1, g2, D0, Mhi):
    d = lambda m: (g1 * m) // g2 + D0
    print(f"beta_M clock (gamma={g1}/{g2}, D0={D0}): "
          "M with beta_M == 0 (parity-free checkpoints):", flush=True)
    zeros = []
    for M in range(1, Mhi + 1):
        par = beta_poly(M, d(M))
        if not par:
            zeros.append(M)
        # THM-2976 T1 cross-check: beta == 0 iff M+1 is a power of two
        dyadic = (M + 1) & M == 0
        require(bool(par) != dyadic, f"clock drift at M={M}")
    print(f"  zeros = {zeros}  (== M+1 power of two: THM-2976 T1 PASS)",
          flush=True)
    # A-clock (anti-checkpoints): levels with A_m = 2^r - 1.  Claim (PROVED:
    # Lucas all-ones + the correction (1+x^{M+1})(1+x)^{d_M} is supported on
    # o <= d_M and o >= M+1, disjoint from the band [d_M+2, M-1]): EVERY band
    # position is forced odd there.  Machine-checked below.
    ac = []
    for m in range(1, Mhi + 1):
        Am = m + d(m) + 1
        if (Am + 1) & Am == 0:          # A_m = 2^r - 1
            par = beta_poly(m, d(m))
            band = list(range(d(m) + 2, m))
            require(all(par.get(o, 0) == 1 for o in band),
                    f"A-clock all-odd band drift at m={m}")
            ac.append((m, Am, len(band)))
    print(f"  anti-checkpoints (A_m = 2^r - 1, band ALL forced odd -- "
          f"machine-checked): {ac}", flush=True)
    return zeros


# ---------------------------------------------------------------- witness mode
def mode_witness(args):
    g1, g2, D0, M = args.g1, args.g2, args.D0, args.M
    biases = BIASES9
    print(f"=== witness anatomy: gamma={g1}/{g2} D0={D0} M={M} ===",
          flush=True)
    clock_report(g1, g2, D0, max(M + 2, 12))
    choice, nodes = dfs_witness(g1, g2, D0, M, biases, args.node_cap,
                                order=args.order, tag="witness-dfs")
    if choice is None:
        print("no witness within cap; abort anatomy", flush=True)
        return
    rho_hist = verify_assignment(g1, g2, D0, M, biases, choice, "witness")

    d, A, Lmax, rows, Env = build(g1, g2, D0, M, biases)
    # ---- per-row anatomy: rho profile + binding bias + saturation
    print("--- per-row anatomy (witness) ---", flush=True)
    idx = 0
    sat_all = []
    for m0 in range(M):
        deltas = []
        for c in rows[m0]:
            deltas.append((c, choice[idx])); idx += 1
        sat = [abs(v) / c[5] for (c, v) in deltas if c[5] > 0]
        sat_all.append(sum(sat) / len(sat))
        rhos = rho_hist[m0]
        b = max(range(len(biases)), key=lambda i: rhos[i])
        print(f"  m={m0+1:3d} maxrho={rhos[b]:.4f} binding_bias="
              f"{biases[b][0]}/{biases[b][1]} mean|delta|/cap="
              f"{sat_all[-1]:.3f} rho(1/2)={rhos[0]:.4f}", flush=True)

    # ---- ledger trajectory: band vs cone split, vs greedy
    md = merged_from_choice(g1, g2, D0, M, choice)
    led_w = Ledger(g1, g2, D0)
    replay = lambda led, m, o, cap, par: md.get((m, o), 0)
    led_g = Ledger(g1, g2, D0)
    print("--- ledger trajectories: witness (replay) vs laneD greedy ---",
          flush=True)
    p12 = Fr(1, 2)
    first_viol_g = None
    div_row = None
    for m in range(1, M + 1):
        led_w.step(m, replay)
        led_g.step(m, greedy_policy)
        bw, cw = led_w.band_split(m, p12)
        bg, cg = led_g.band_split(m, p12)
        okw, rw = envelope_check(led_w, biases)
        okg, rg = envelope_check(led_g, biases)
        require(okw, f"witness ledger violates envelope at m={m}?!")
        if not okg and first_viol_g is None:
            first_viol_g = m
        if div_row is None:
            recw = {(o, dl) for (o, c, p_, dl) in led_w.rows[-1][2]}
            recg = {(o, dl) for (o, c, p_, dl) in led_g.rows[-1][2]}
            if recw != recg:
                div_row = m
        dmm = led_w.d(m)
        bandpos = f"[{dmm+2},{m-1}]" if m - 1 >= dmm + 2 else "empty"
        print(f"  m={m:3d} band{bandpos} witness: D(1/2)={float(bw+cw):+.3e} "
              f"band={float(bw):+.3e} cone={float(cw):+.3e} maxrho="
              f"{max(rw):.4f} | greedy: D(1/2)={float(bg+cg):+.3e} "
              f"band={float(bg):+.3e} {'VIOL' if not okg else 'ok  '}",
              flush=True)
    print(f"witness-vs-greedy: first divergent row = {div_row}; greedy first "
          f"envelope violation at m = {first_viol_g} (band birth m* = "
          f"{(D0+2)/(1-g1/g2):.1f})", flush=True)

    # band coefficient vectors at final level
    dmm = led_w.d(M)
    band_w = {o: led_w.b.get(o, 0) for o in range(dmm + 2, M)}
    band_g = {o: led_g.b.get(o, 0) for o in range(dmm + 2, M)}
    print(f"final band coeffs (o in [{dmm+2},{M-1}]): witness={band_w}",
          flush=True)
    print(f"                                          greedy ={band_g}",
          flush=True)
    beta = beta_poly(M, dmm)
    for o in range(dmm + 2, M):
        require((band_w.get(o, 0) % 2 == beta.get(o, 0)) and
                (band_g.get(o, 0) % 2 == beta.get(o, 0)),
                "forced band parity drift")
    print("forced band parity (beta_M restriction): PASS for both", flush=True)
    return choice


# ---------------------------------------------------------------- beam mode
def mode_beam(args):
    g1, g2, D0, M, W = args.g1, args.g2, args.D0, args.M, args.width
    biases = list(BIASES9)
    if args.extra:
        for tok in args.extra.split(","):
            n, dd = tok.split("/")
            biases.append((int(n), int(dd)))
        print(f"augmented bias set (+{args.extra}): {len(biases)} biases",
              flush=True)
    scoring = args.scoring
    rng = random.Random(args.seed)
    d, A, Lmax, rows, Env = build(g1, g2, D0, M, biases)
    nP = len(biases)
    print(f"=== beam: gamma={g1}/{g2} D0={D0} Mtarget={M} W={W} "
          f"scoring={scoring} seed={args.seed} rand_keep={args.rand_keep} ===",
          flush=True)
    fEnvRow = None

    # state: (S tuple, hist link, hmax float); hist link = (v, parent)
    states = [(tuple(0 for _ in range(nP)), None, 0.0)]
    survived = 0
    best_final = None
    for m0 in range(M):
        t0 = time.time()
        cells = rows[m0]
        sw = [tuple(0 for _ in range(nP))]
        for c in reversed(cells):
            sw.append(tuple(sw[-1][i] + c[5] * c[7][i] for i in range(nP)))
        sw.reverse()
        fEnvRow = [float(Env[m0][i]) for i in range(nP)]
        for ci, (m, side, k, z, o, cap, par, w) in enumerate(cells):
            rest = sw[ci + 1]
            new = {}
            for (S, hist, hmax) in states:
                for v in choices_of(cap, par):
                    S2 = tuple(S[i] + v * w[i] for i in range(nP)) if v else S
                    dead = False
                    for i in range(nP):
                        if abs(S2[i]) > rest[i] + Env[m0][i]:
                            dead = True
                            break
                    if dead:
                        continue
                    if S2 not in new:
                        new[S2] = ((v, hist), hmax)
            if not new:
                print(f"[beam] DIED in row {m0+1} at cell {ci} "
                      f"(k={k} side={side}); survived through row {survived}",
                      flush=True)
                return survived, best_final
            scored = []
            frest = [float(rest[i]) for i in range(nP)]
            for S2, (hist2, hmax) in new.items():
                if scoring == "max":
                    sc = max(float(abs(S2[i])) / (frest[i] + fEnvRow[i])
                             for i in range(nP))
                else:
                    sc = sum(float(abs(S2[i])) / (frest[i] + fEnvRow[i])
                             for i in range(nP)) / nP
                scored.append((sc, S2, hist2, hmax))
            scored.sort(key=lambda t: t[0])
            keep = scored[:W]
            if args.rand_keep and len(scored) > W:
                extra = rng.sample(scored[W:], min(args.rand_keep,
                                                   len(scored) - W))
                keep = keep + extra
            states = [(S2, hist2, hmax) for (_, S2, hist2, hmax) in keep]
        # row end: all states satisfy |S| <= Env exactly (rest = 0 above)
        end = []
        for (S, hist, hmax) in states:
            rho = max(float(abs(S[i])) / fEnvRow[i] for i in range(nP))
            end.append((rho, S, hist, max(hmax, rho)))
        end.sort(key=lambda t: t[0])
        states = [(S, hist, hm) for (_, S, hist, hm) in end]
        survived = m0 + 1
        bi = max(range(nP), key=lambda i: abs(end[0][1][i]) / fEnvRow[i])
        print(f"[beam] row {m0+1:3d} done: states={len(states)} "
              f"bestrho={end[0][0]:.6f} binding={biases[bi][0]}/"
              f"{biases[bi][1]} min_hmax={min(hm for *_, hm in end):.6f} "
              f"t={time.time()-t0:.1f}s", flush=True)
        best_final = end[0]
    # reconstruct + verify best final state
    rho, S, hist, hmax = best_final
    choice = []
    while hist is not None:
        v, hist = hist
        choice.append(v)
    choice.reverse()
    verify_assignment(g1, g2, D0, M, biases, choice, f"beam-W{W}")
    bad, tot = antisym_defect(g1, g2, D0, M, choice)
    print(f"[beam] SURVIVED all {M} rows; best final rho={rho:.6f} "
          f"hmax={hmax:.6f}; antisym defect {bad}/{tot} cells", flush=True)
    if args.dump:
        print("choice=", choice, flush=True)
    if args.dense:
        dense_check(g1, g2, D0, M, choice, args.dense, f"beam-W{W}")
    return survived, best_final


# ---------------------------------------------------------------- prefix mode
def mode_prefixgreedy(args):
    g1, g2, D0 = args.g1, args.g2, args.D0
    biases = BIASES9
    Mw = args.M
    print(f"=== prefix-greedy: gamma={g1}/{g2} D0={D0} witness M={Mw}, "
          f"greedy continuation to M={args.Mcont} ===", flush=True)
    choice, nodes = dfs_witness(g1, g2, D0, Mw, biases, args.node_cap,
                                order=args.order, tag="prefix-dfs")
    md = merged_from_choice(g1, g2, D0, Mw, choice) if choice else {}
    for P in args.prefixes:
        require(P <= Mw or not choice, "prefix longer than witness")
        led = Ledger(g1, g2, D0)
        first_viol = None
        res = []
        for m in range(1, args.Mcont + 1):
            if m <= P:
                led.step(m, lambda led, mm, o, cap, par: md.get((mm, o), 0))
            else:
                led.step(m, greedy_policy)
            ok, rhos = envelope_check(led, biases)
            if not ok and first_viol is None:
                first_viol = m
            if m in (args.Mcont, 20):
                res.append((m, float(led.val(Fr(1, 2)))))
        tail = " ".join(f"D(1/2)@M={m}:{v:+.3e}" for m, v in res)
        print(f"  P={P:3d}: first envelope violation at m={first_viol} "
              f"(P + {first_viol - P if first_viol else '?'}) {tail}",
              flush=True)


# ---------------------------------------------------------------- symdfs mode
def mode_symdfs(args):
    g1, g2, D0, Mtry = args.g1, args.g2, args.D0, args.M
    biases = BIASES9
    d = lambda m: (g1 * m) // g2 + D0
    A = lambda m: m + d(m) + 1
    Lmax = A(Mtry)
    nP = len(biases)
    print(f"=== symdfs (mirror-symmetric corridor): gamma={g1}/{g2} D0={D0} "
          f"M={Mtry} node_cap={args.node_cap} ===", flush=True)
    rows = []
    for m in range(1, Mtry + 1):
        dm = d(m)
        cells = []
        for k in range(dm + 1):
            cap = comb(dm, k)
            w0 = tuple(n**(m + dm - k) * (dd - n)**(k + 1)
                       * dd**(Lmax - A(m)) for (n, dd) in biases)
            w1 = tuple(n**(k + 1) * (dd - n)**(m + dm - k)
                       * dd**(Lmax - A(m)) for (n, dd) in biases)
            w = tuple(a + b for a, b in zip(w0, w1))
            cells.append((cap, cap % 2, w))
        rows.append(cells)
    Env = [tuple((n**(m + 2) + (dd - n)**(m + 2)) * dd**(Lmax - m - 2)
                 for (n, dd) in biases) for m in range(Mtry)]
    R = [tuple(sum(c[0] * c[2][i] for c in rows[m]) for i in range(nP))
         for m in range(Mtry)]
    allowed = []
    for m in range(Mtry):
        best = None
        for mp in range(m, Mtry):
            tot = tuple(Env[mp][i] + sum(R[r][i] for r in range(m + 1, mp + 1))
                        for i in range(nP))
            best = tot if best is None else tuple(min(a, b)
                                                  for a, b in zip(best, tot))
        allowed.append(best)
    suffix = []
    for m in range(Mtry):
        sw = [tuple(0 for _ in range(nP))]
        for c in reversed(rows[m]):
            sw.append(tuple(sw[-1][i] + c[0] * c[2][i] for i in range(nP)))
        sw.reverse()
        suffix.append(sw)
    state = {"nodes": 0}

    def dfs(m, ci, S, hist):
        state["nodes"] += 1
        if state["nodes"] > args.node_cap:
            raise TimeoutError
        cells = rows[m]
        if ci == len(cells):
            for i in range(nP):
                if abs(S[i]) > Env[m][i] or abs(S[i]) > allowed[m][i]:
                    return None
            if m + 1 == Mtry:
                return list(hist)
            return dfs(m + 1, 0, S, hist)
        rest = suffix[m][ci]
        for i in range(nP):
            if abs(S[i]) > rest[i] + allowed[m][i]:
                return None
        cap, par, w = cells[ci]
        vals = choices_of(cap, par)
        if args.order == "abs":
            vals.sort(key=abs)
        else:
            vals.sort(key=lambda v: -abs(v))
        for v in vals:
            S2 = tuple(S[i] + v * w[i] for i in range(nP)) if v else S
            hist.append(v)
            r = dfs(m, ci + 1, S2, hist)
            if r is not None:
                return r
            hist.pop()
        return None

    t0 = time.time()
    try:
        wtn = dfs(0, 0, tuple(0 for _ in range(nP)), [])
        verdict = "SYM-WITNESS" if wtn else "SYM-INFEASIBLE (exhausted)"
    except TimeoutError:
        wtn, verdict = None, "NODE-CAP"
    print(f"[symdfs] M={Mtry}: {verdict} nodes={state['nodes']} "
          f"t={time.time()-t0:.1f}s", flush=True)
    if wtn:
        # expand symmetric choice to the full per-cell list and re-verify
        full = []
        idx = 0
        for m in range(1, Mtry + 1):
            for k in range(d(m) + 1):
                full.extend([wtn[idx], wtn[idx]])
                idx += 1
        verify_assignment(g1, g2, D0, Mtry, biases, full, "symdfs")
        print("sym witness (per (m,k), both sides equal):", wtn, flush=True)
    return wtn


# ---------------------------------------------------------------- asymdfs mode
def mode_asymdfs(args):
    """Antisymmetric (complementary-mirror) restricted DFS: one variable v
    per (m,k) with delta_side0 = v, delta_side1 = -v; weight w0 - w1.
    D(1/2) == 0 identically in this class.  Exhaustion gives an in-class
    THEOREM (no complementary-mirror scheme passes the envelope at M)."""
    g1, g2, D0, Mtry = args.g1, args.g2, args.D0, args.M
    biases = BIASES9
    d = lambda m: (g1 * m) // g2 + D0
    A = lambda m: m + d(m) + 1
    Lmax = A(Mtry)
    nP = len(biases)
    print(f"=== asymdfs (complementary-mirror corridor): gamma={g1}/{g2} "
          f"D0={D0} M={Mtry} node_cap={args.node_cap} ===", flush=True)
    rows = []
    for m in range(1, Mtry + 1):
        dm = d(m)
        cells = []
        for k in range(dm + 1):
            cap = comb(dm, k)
            w0 = tuple(n**(m + dm - k) * (dd - n)**(k + 1)
                       * dd**(Lmax - A(m)) for (n, dd) in biases)
            w1 = tuple(n**(k + 1) * (dd - n)**(m + dm - k)
                       * dd**(Lmax - A(m)) for (n, dd) in biases)
            w = tuple(a - b for a, b in zip(w0, w1))
            cells.append((cap, cap % 2, w))
        rows.append(cells)
    Env = [tuple((n**(m + 2) + (dd - n)**(m + 2)) * dd**(Lmax - m - 2)
                 for (n, dd) in biases) for m in range(Mtry)]
    R = [tuple(sum(c[0] * abs(c[2][i]) for c in rows[m]) for i in range(nP))
         for m in range(Mtry)]
    allowed = []
    for m in range(Mtry):
        best = None
        for mp in range(m, Mtry):
            tot = tuple(Env[mp][i] + sum(R[r][i] for r in range(m + 1, mp + 1))
                        for i in range(nP))
            best = tot if best is None else tuple(min(a, b)
                                                  for a, b in zip(best, tot))
        allowed.append(best)
    suffix = []
    for m in range(Mtry):
        sw = [tuple(0 for _ in range(nP))]
        for c in reversed(rows[m]):
            sw.append(tuple(sw[-1][i] + c[0] * abs(c[2][i])
                            for i in range(nP)))
        sw.reverse()
        suffix.append(sw)
    state = {"nodes": 0}

    def dfs(m, ci, S, hist):
        state["nodes"] += 1
        if state["nodes"] > args.node_cap:
            raise TimeoutError
        cells = rows[m]
        if ci == len(cells):
            for i in range(nP):
                if abs(S[i]) > Env[m][i] or abs(S[i]) > allowed[m][i]:
                    return None
            if m + 1 == Mtry:
                return list(hist)
            return dfs(m + 1, 0, S, hist)
        rest = suffix[m][ci]
        for i in range(nP):
            if abs(S[i]) > rest[i] + allowed[m][i]:
                return None
        cap, par, w = cells[ci]
        vals = choices_of(cap, par)
        if all(x == 0 for x in w):
            vals = [vals[0]]        # weightless (colliding) cell: canonical
        elif args.order == "abs":
            vals.sort(key=abs)
        else:
            vals.sort(key=lambda v: -abs(v))
        for v in vals:
            S2 = tuple(S[i] + v * w[i] for i in range(nP)) if v else S
            hist.append(v)
            r = dfs(m, ci + 1, S2, hist)
            if r is not None:
                return r
            hist.pop()
        return None

    t0 = time.time()
    try:
        wtn = dfs(0, 0, tuple(0 for _ in range(nP)), [])
        verdict = "ASYM-WITNESS" if wtn else "ASYM-INFEASIBLE (exhausted)"
    except TimeoutError:
        wtn, verdict = None, "NODE-CAP"
    print(f"[asymdfs] M={Mtry}: {verdict} nodes={state['nodes']} "
          f"t={time.time()-t0:.1f}s", flush=True)
    if wtn:
        full = []
        idx = 0
        for m in range(1, Mtry + 1):
            for k in range(d(m) + 1):
                full.extend([wtn[idx], -wtn[idx]])
                idx += 1
        verify_assignment(g1, g2, D0, Mtry, biases, full, "asymdfs")
        if args.dense:
            dense_check(g1, g2, D0, Mtry, full, args.dense, "asymdfs")
    return wtn


# ---------------------------------------------------------------- selftest
def selftest():
    # cross-check scaled-S vs Ledger.val on a small witness (M=5, gamma=1/2)
    g1, g2, D0, M = 1, 2, 0, 5
    biases = BIASES9
    choice, _ = dfs_witness(g1, g2, D0, M, biases, 10**6, tag="selftest-dfs")
    require(choice is not None, "selftest: no witness at M=5?!")
    d, A, Lmax, rows, Env = build(g1, g2, D0, M, biases)
    nP = len(biases)
    S = [0] * nP
    idx = 0
    for m0 in range(M):
        for c in rows[m0]:
            v = choice[idx]; idx += 1
            for i in range(nP):
                S[i] += v * c[7][i]
    md = merged_from_choice(g1, g2, D0, M, choice)
    led = Ledger(g1, g2, D0)
    for m in range(1, M + 1):
        led.step(m, lambda led, mm, o, cap, par: md.get((mm, o), 0))
    for i, (n, dd) in enumerate(biases):
        p = Fr(n, dd)
        require(led.val(p) == Fr(S[i], 2 * dd**Lmax),
                f"selftest: S vs ledger value mismatch at bias {n}/{dd}")
    print("SELFTEST PASS: scaled S == 2 D b^Lmax == descended ledger value "
          "at all 9 biases (exact)", flush=True)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    for name in ("witness", "beam", "prefixgreedy", "symdfs", "asymdfs"):
        p = sub.add_parser(name)
        p.add_argument("--g1", type=int, default=1)
        p.add_argument("--g2", type=int, default=2)
        p.add_argument("--D0", type=int, default=0)
        p.add_argument("--M", type=int, default=10)
        p.add_argument("--node-cap", type=int, default=20_000_000)
        p.add_argument("--order", default="abs")
        if name == "beam":
            p.add_argument("--width", type=int, default=256)
            p.add_argument("--scoring", default="max",
                           choices=("max", "sum"))
            p.add_argument("--seed", type=int, default=0)
            p.add_argument("--rand-keep", type=int, default=32)
            p.add_argument("--dump", action="store_true")
            p.add_argument("--dense", type=int, default=0)
            p.add_argument("--extra", default="")
        if name == "prefixgreedy":
            p.add_argument("--Mcont", type=int, default=30)
            p.add_argument("--prefixes", type=int, nargs="+",
                           default=[0, 2, 4, 6, 8, 10])
        if name == "asymdfs":
            p.add_argument("--dense", type=int, default=0)
    sp = sub.add_parser("selftest")
    cp = sub.add_parser("clock")
    cp.add_argument("--g1", type=int, default=1)
    cp.add_argument("--g2", type=int, default=2)
    cp.add_argument("--D0", type=int, default=0)
    cp.add_argument("--M", type=int, default=100)
    args = ap.parse_args()
    if args.mode == "selftest":
        selftest()
    elif args.mode == "clock":
        clock_report(args.g1, args.g2, args.D0, args.M)
    elif args.mode == "witness":
        mode_witness(args)
    elif args.mode == "beam":
        mode_beam(args)
    elif args.mode == "prefixgreedy":
        mode_prefixgreedy(args)
    elif args.mode == "symdfs":
        mode_symdfs(args)
    elif args.mode == "asymdfs":
        mode_asymdfs(args)


if __name__ == "__main__":
    main()
