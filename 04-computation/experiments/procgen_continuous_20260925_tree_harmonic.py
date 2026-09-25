#!/usr/bin/env python3
"""procgen_continuous_20260925_tree_harmonic.py

Lane: natural boundaries, Mahler/Cobham and harmonic trees
(session collatz-procgen-20260922, 2026-09-25).

Part C: random walk, branching number, harmonic measure and Tutte's barycentric embedding
on the backward (inverse) Collatz trees.

  TC1  PROVED lower bound for the branching number (worst-case ladder flow).  For q = 3 (either
       sign) every non-multiple of 3 has a backward tree with br >= lambda_3 = u^(-1/2),
       2u^3 + u^2 = 1 (lambda_3 = 1.2334); for 5n+1, br >= lambda_5 = u^(-1/4),
       2u^5 + u^4 + u^3 + u^2 = 1 (lambda_5 = 1.1230).  So simple random walk is transient on
       every such tree (Lyons 1990).  The script (a) solves for the roots, (b) checks the worst
       case over all ladder types, and (c) checks the flow bound g_D(root; lambda) >= phi on the
       actual truncated trees for D up to 46 (60 for 5n+1).
  TC2  Level sizes |T_n| and growth ratios (FINITE-EXACT) for the 3n+1 tree of 1, the 3n-1 trees
       of 1, 5, 17, the 5n+1 tree of 1, and a Haar-model multitype Galton-Watson tree.
  TC3  Branching-number profile: g_D(root; lambda) (max-flow strength with capacities
       lambda^-|e|) on a lambda grid for D = 30, 38, 46 (FINITE-EXACT numbers; the estimate
       br ~ gr ~ 4/3 is numerical).
  TC4  Harmonic measure of simple random walk (root reflecting, exit at depth D = 46): level
       entropies H_n, entropy dimension in the 2-adic metric, E-move frequency along
       harmonic rays, compared with log2 of the growth (dimension drop).
  TC5  Tutte's barycentric embedding: leaves at depth D placed at exp(2 pi i 0.m1 m2 ... mD)
       (m = 1 for an odd move), interior vertices at the average of their neighbours; the
       embedding equals the harmonic-measure average; convergence in D is reported.
  TC6  Transport (PROVED; checked): the T_+ tree of x is the negative of the T_- tree of -x.
  TC7  Cheeger constant 0 (PROVED; bare doubling rays).

Every check raises on failure.  Peak memory about 600 MB; runtime about 2 min.
"""
import sys
import time
import math

import numpy as np

T0 = time.time()


def check(cond, msg):
    if not cond:
        raise AssertionError("CHECK FAILED: " + msg)


def mem(tag):
    import resource
    print(f"[mem] {tag}: max RSS so far {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} MB",
          file=sys.stderr)


def hdr(s):
    mem("before " + s[:4])
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


def bisect_root(f, a, b, it=200):
    fa = f(a)
    for _ in range(it):
        m = (a + b) / 2
        fm = f(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return (a + b) / 2


# ---------------------------------------------------------------------------------------------
def build_tree(q, b, root, D):
    """Backward tree of `root` under T_{q,b}: children of x are 2x and (2x-b)/q when integral;
    the root is never re-entered (this cuts the root cycle).  Returns per-level arrays."""
    vals = [np.array([root], dtype=np.int64)]
    par = [np.array([-1], dtype=np.int64)]
    isE = [np.array([0], dtype=np.int8)]
    for n in range(D):
        X = vals[-1]
        check(int(X.max()) < 2**62, "tree value overflow guard")
        dch = 2 * X
        dpar = np.arange(X.size, dtype=np.int64)
        kd = dch != root
        dch, dpar = dch[kd], dpar[kd]
        m = ((2 * X - b) % q) == 0
        ech = (2 * X[m] - b) // q
        epar = np.nonzero(m)[0].astype(np.int64)
        ke = (ech * (1 if root > 0 else -1) >= 1) & (ech != root)
        ech, epar = ech[ke], epar[ke]
        vals.append(np.concatenate([dch, ech]))
        par.append(np.concatenate([dpar, epar]))
        isE.append(np.concatenate([np.zeros(dch.size, np.int8), np.ones(ech.size, np.int8)]))
    return vals, par, isE


def build_gw(D, seed):
    """Haar model: types = residues mod 3.  Type 1 -> one D-child of type 2; type 2 -> D-child of
    type 1 and E-child of uniform type in {0,1,2}; type 0 -> one D-child of type 0 (bare ray).
    Root type 1 (like the root 1)."""
    rng = np.random.default_rng(seed)
    typ = [np.array([1], dtype=np.int8)]
    par = [np.array([-1], dtype=np.int64)]
    isE = [np.array([0], dtype=np.int8)]
    for n in range(D):
        t = typ[-1]
        dtype_ = np.where(t == 0, 0, 3 - t).astype(np.int8)     # D swaps 1<->2, fixes 0
        dpar = np.arange(t.size, dtype=np.int64)
        m = t == 2
        epar = np.nonzero(m)[0].astype(np.int64)
        etype = rng.integers(0, 3, size=epar.size).astype(np.int8)
        typ.append(np.concatenate([dtype_, etype]))
        par.append(np.concatenate([dpar, epar]))
        isE.append(np.concatenate([np.zeros(dpar.size, np.int8), np.ones(epar.size, np.int8)]))
    return typ, par, isE


def flow_strength(par, D, lam):
    """Max-flow strength g_D(root) with edge capacities lam^-|e| (normalized): g(leaf at D) = 1,
    g(v) = lam^-1 sum_children min(1, g(c))."""
    g = np.ones(par[D].size)
    for n in range(D, 0, -1):
        s = np.bincount(par[n], weights=np.minimum(1.0, g), minlength=par[n - 1].size)
        g = s / lam
    return float(g[0])


LIVE_RES = {3: {1, 2}, 5: {1, 2, 3, 4}, 7: {1, 2, 4}}


def flow_min_live(vals, par, D, lam, q, excluded):
    """min over live nodes (not multiples of q, not in `excluded`) at depth <= D-1 of g_D(v).
    `excluded` = the root and the nodes T(root)/2^k, whose ladder contains the root (the rung that
    closes the root cycle is cut from the tree); they keep positive flow but not the uniform bound."""
    g = np.ones(par[D].size)
    mn = np.inf
    exc = np.array(sorted(excluded), dtype=np.int64)
    for n in range(D, 0, -1):
        s = np.bincount(par[n], weights=np.minimum(1.0, g), minlength=par[n - 1].size)
        g = s / lam
        live = np.isin(vals[n - 1] % q, np.array(sorted(LIVE_RES[q]))) & ~np.isin(vals[n - 1], exc)
        if live.any():
            mn = min(mn, float(g[live].min()))
    return mn, float(g[0])


def harmonic(par, isE, D):
    """Unit current / harmonic measure from the root to the level-D boundary (unit conductances).
    Returns per-level theta arrays and per-level E-count arrays."""
    w = [None] * (D + 1)
    Cn = np.full(par[D].size, np.inf)
    for n in range(D, 0, -1):
        wn = np.where(np.isinf(Cn), 1.0, Cn / (1.0 + Cn))      # series: edge 1 + subtree 1/C
        w[n] = wn
        Cn = np.bincount(par[n], weights=wn, minlength=par[n - 1].size)
    theta = [np.array([1.0])]
    ecnt = [np.array([0], dtype=np.int16)]
    # recompute C per level top-down needs the sums; store them
    sums = [None] * (D + 1)
    for n in range(D, 0, -1):
        sums[n - 1] = np.bincount(par[n], weights=w[n], minlength=par[n - 1].size)
    for n in range(1, D + 1):
        p = par[n]
        th = theta[-1][p] * w[n] / sums[n - 1][p]
        theta.append(th)
        ecnt.append((ecnt[-1][p] + isE[n]).astype(np.int16))
    return theta, ecnt


def tutte(par, isE, D):
    """Barycentric (harmonic) embedding: leaves at depth D at exp(2 pi i ang), ang = 0.m1...mD;
    returns (h(root), Dirichlet energy, mean |h| by level)."""
    ang = [np.array([0.0])]
    for n in range(1, D + 1):
        ang.append(ang[-1][par[n]] + isE[n] * 2.0**(-n))
    alpha = np.exp(2j * np.pi * ang[D])
    beta = np.zeros(par[D].size)
    A = [None] * (D + 1)
    B = [None] * (D + 1)
    A[D], B[D] = alpha, beta
    for n in range(D, 0, -1):
        k = np.bincount(par[n], minlength=par[n - 1].size).astype(float)
        sa = np.bincount(par[n], weights=A[n].real, minlength=par[n - 1].size) + \
            1j * np.bincount(par[n], weights=A[n].imag, minlength=par[n - 1].size)
        sb = np.bincount(par[n], weights=B[n], minlength=par[n - 1].size)
        if n - 1 == 0:
            A[0] = sa / (k - sb)
            B[0] = np.zeros(1)
        else:
            den = k + 1.0 - sb
            A[n - 1] = sa / den
            B[n - 1] = 1.0 / den
    h = [A[0]]
    energy = 0.0
    radii = [float(abs(A[0][0]))]
    for n in range(1, D + 1):
        hn = A[n] + B[n] * h[-1][par[n]]
        energy += float(np.sum(np.abs(hn - h[-1][par[n]])**2))
        radii.append(float(np.mean(np.abs(hn))))
        h.append(hn)
    return complex(A[0][0]), energy, radii


def main():
    # -----------------------------------------------------------------------------------------
    hdr("TC1  PROVED lower bound for the branching number (worst-case ladder flow)")
    u3 = bisect_root(lambda u: 2 * u**3 + u**2 - 1, 0.1, 0.99)
    lam3 = u3**-0.5
    u5 = bisect_root(lambda u: 2 * u**5 + u**4 + u**3 + u**2 - 1, 0.1, 0.99)
    lam5 = u5**-0.25
    print(f"  q=3: u3 = {u3:.10f} (2u^3+u^2=1), lambda_3 = u3^(-1/2) = {lam3:.10f}")
    print(f"  q=5: u5 = {u5:.10f} (2u^5+u^4+u^3+u^2=1), lambda_5 = u5^(-1/4) = {lam5:.10f}")

    # generic worst-case ladder sum: chain positions h0 + o*j (o = ord_q 2), rung residues r0 + c*j mod q,
    # a rung is live iff its residue lies in LIVE[q] (residues whose doubling orbit meets the legal class)
    LIVE = {3: {1, 2}, 5: {1, 2, 3, 4}, 7: {1, 2, 4}}
    ORD = {3: 2, 5: 4, 7: 3}
    INC = {3: 1, 5: 3, 7: 1}

    def ladder_sum(lam, q, h0, r0, J=3000):
        tot = 0.0
        for j in range(J):
            if (r0 + INC[q] * j) % q in LIVE[q]:
                tot += lam ** -(h0 + ORD[q] * j + 1)
        return tot

    def worst(lam, q):
        return min((ladder_sum(lam, q, h0, r0), h0, r0) for h0 in range(ORD[q]) for r0 in range(q))
    lam7 = bisect_root(lambda l: worst(l, 7)[0] - 1, 1.0001, 1.5)
    for q, lam in ((3, lam3), (5, lam5), (7, lam7)):
        w = worst(lam, q)
        print(f"  q={q}: at lambda = {lam:.10f}, min over ladder types (h0, first-rung residue r0) of")
        print(f"        sum_(live rungs) lambda^-(distance) = {w[0]:.12f} at (h0, r0) = ({w[1]}, {w[2]})")
        check(abs(w[0] - 1) < 1e-8, f"worst ladder sum is 1 at lambda_{q}")
        check(worst(lam * 1.001, q)[0] < 1, "the bound is sharp for the worst type")
    print(f"  q=7 (7n+1): live residues {{1,2,4}} mod 7, legal every 3 doublings, rungs 8p+1 = p+1 mod 7:")
    print(f"        lambda_7 = {lam7:.10f} (numerical root; the same induction applies)")
    print("  Ladder facts used (PROVED, Prop. 11 of the mod-192 note and its 5n+1 analogue):")
    print("    q=3: E is legal on the doubling chain of x at positions h0, h0+2, h0+4, ...; the rungs")
    print("         satisfy S(p) = 4p+1 = p+1 mod 3 (b=+1) or 4p-1 = p-1 mod 3 (b=-1), so exactly one")
    print("         rung in three consecutive is a multiple of 3 (a bare ray).")
    print("    q=5: legal positions h0, h0+4, ...; rungs S(p) = 16p+3 = p+3 mod 5: one dead rung in five.")
    print("  With phi <= lambda-1, induction from the truncation depth up gives g_D(v) >= phi for every")
    print("  live v and every D, hence br >= lambda_q; Lyons (1990) Thm: br > 1 => SRW transient.")
    print("  Percolation corollary: p_c = 1/br <= 1/lambda_3 = %.4f on every live 3n+-1 tree." % (1 / lam3))

    specs = [("3n+1, root 1", 3, 1, 1, 46), ("3n-1, root 1", 3, -1, 1, 46), ("3n-1, root 5", 3, -1, 5, 44),
             ("3n-1, root 17", 3, -1, 17, 42), ("5n+1, root 1", 5, 1, 1, 60), ("7n+1, root 1", 7, 1, 1, 60)]
    trees = {}
    for name, q, b, root, D in specs:
        vals, par, isE = build_tree(q, b, root, D)
        trees[name] = (q, b, root, D, vals, par, isE)
    print("  flow check on the actual truncated trees at lambda = lambda_q, phi = min(1, lambda_q - 1):")
    for name, (q, b, root, D, vals, par, isE) in trees.items():
        lam = {3: lam3, 5: lam5, 7: lam7}[q]
        phi = min(1.0, lam - 1)
        t = (q * root + b) // 2 if root % 2 else root // 2
        excluded = {root}
        while t % 1 == 0 and t >= 1:
            excluded.add(t)
            if t % 2:
                break
            t //= 2
        mn, groot = flow_min_live(vals, par, D, lam, q, excluded)
        print(f"    {name:<15s} D={D}: min over live nodes (excluding {sorted(excluded)}) of g_D(v; lambda_q) = "
              f"{mn:.4f} (>= phi = {phi:.4f}); g_D(root) = {groot:.4f} > 0")
        check(mn >= phi - 1e-12, "uniform flow lower bound on live nodes")
        check(groot > 0.05, "positive flow at the root")

    # -----------------------------------------------------------------------------------------
    hdr("TC6  Transport (PROVED): the T_+ tree of 1 is the negative of the T_- tree of -1")
    vp, pp, ep = build_tree(3, 1, 1, 30)
    vm, pm, em = build_tree(3, -1, -1, 30)
    for n in range(31):
        check(np.array_equal(np.sort(vp[n]), np.sort(-vm[n])), f"level {n} transport")
    print("  levels 0..30 agree as sets after x -> -x (so every graph invariant of backward trees is")
    print("  the same on the two sheets for mirror roots; the sheets differ only in WHICH roots are positive).")
    del vp, pp, ep, vm, pm, em

    # -----------------------------------------------------------------------------------------
    hdr("TC2  Level sizes and growth (FINITE-EXACT)")
    gw = {}
    for seed in (1, 2):
        typ, par, isE = build_gw(46, seed)
        gw[f"Haar-model GW #{seed}"] = (par, isE, 46)
    rows = []
    for name, (q, b, root, D, vals, par, isE) in trees.items():
        sizes = [v.size for v in vals]
        rows.append((name, D, sizes))
    for name, (par, isE, D) in gw.items():
        rows.append((name, D, [p.size for p in par]))
    for name, D, sizes in rows:
        ratios = [sizes[i + 1] / sizes[i] for i in range(D - 6, D)]
        grn = sizes[D] ** (1 / D)
        print(f"  {name:<18s} |T_D| = {sizes[D]:>9d} (D={D}); last ratios " +
              " ".join(f"{r:.4f}" for r in ratios) + f"; |T_D|^(1/D) = {grn:.4f}")
    print("  => growth 4/3 per level for 3n+-1 (the Haar mean E[N_d] = (4/3)^d), 6/5 for 5n+1, 8/7 for 7n+1.")

    # -----------------------------------------------------------------------------------------
    hdr("TC3  Branching-number profile: g_D(root; lambda) for several D (FINITE-EXACT)")
    lams = [1.20, 1.2334, 1.26, 1.29, 1.31, 1.32, 1.33, 1.3333, 1.34, 1.36, 1.40]
    for name, (q, b, root, D, vals, par, isE) in list(trees.items()) + [(k, (3, 1, 1, v[2], None, v[0], v[1]))
                                                                        for k, v in gw.items()]:
        if q != 3:
            continue
        Ds = [d for d in (30, 36, 42, 46) if d <= D]
        print(f"  {name}:  lambda -> g_D(root) for D in {Ds}")
        for lam in lams:
            vs = [flow_strength(par, d, lam) for d in Ds]
            print(f"     lambda={lam:<7} " + "  ".join(f"{v:.3e}" for v in vs))
    name = "5n+1, root 1"
    q, b, root, D, vals, par, isE = trees[name]
    print(f"  {name}:  lambda -> g_D(root) for D in (40, 50, 60)")
    for lam in (1.10, 1.1230, 1.15, 1.18, 1.19, 1.20, 1.21, 1.25):
        vs = [flow_strength(par, d, lam) for d in (40, 50, 60)]
        print(f"     lambda={lam:<7} " + "  ".join(f"{v:.3e}" for v in vs))
    print("  Reading: g_D(root) stabilizes for lambda below about 4/3 (6/5 for 5n+1) and decays in D")
    print("  above it, consistent with br = gr = 4/3 (resp. 6/5); only br >= 1.2334 (1.1230) is PROVED.")

    # -----------------------------------------------------------------------------------------
    hdr("TC4  Harmonic measure of SRW (root reflecting, boundary at depth D)")
    allrows = [(name, t[5], t[6], t[3], t[4]) for name, t in trees.items()] + \
              [(k, v[0], v[1], v[2], None) for k, v in gw.items()]
    for name, par, isE, D, vals in allrows:
        theta, ecnt = harmonic(par, isE, D)
        tot = float(theta[D].sum())
        check(abs(tot - 1) < 1e-9, "harmonic measure has mass 1")
        H = [float(-(t[t > 0] * np.log(t[t > 0])).sum()) for t in theta]
        n1, n2 = D // 3, (2 * D) // 3
        dimH = (H[n2] - H[n1]) / ((n2 - n1) * math.log(2))
        fE = float((theta[n2] * ecnt[n2]).sum()) / n2
        sizes = [p.size for p in par]
        dimB = math.log2(sizes[D] / sizes[D - 10]) / 10
        extra = ""
        if vals is not None:
            ls = float((theta[n2] * np.log(vals[n2].astype(float))).sum()) / n2
            extra = f"; E_theta[log x_n]/n = {ls:.4f}"
        print(f"  {name:<18s} D={D}: entropy dimension (2-adic) on levels [{n1},{n2}] = {dimH:.4f}; "
              f"log2(growth) = {dimB:.4f}; E-frequency on harmonic rays = {fE:.4f}{extra}")
        del theta, ecnt
    print("  Dimension drop: harmonic measure lives on a set of 2-adic dimension well below the")
    print("  boundary's log2(4/3) = %.4f (cf. Lyons-Pemantle-Peres 1995 for Galton-Watson trees)." % math.log2(4 / 3))

    # -----------------------------------------------------------------------------------------
    hdr("TC5  Tutte's barycentric embedding (boundary on the unit circle at the move-word angle)")
    for name, par, isE, D, vals in allrows:
        out = []
        for d in sorted(set([D - 8, D - 4, D])):
            hroot, en, radii = tutte(par, isE, d)
            out.append((d, hroot, en, radii))
        s = "; ".join(f"D={d}: h(root)={hr.real:+.5f}{hr.imag:+.5f}i, energy={en:.5f}" for d, hr, en, _ in out)
        print(f"  {name:<18s} {s}")
        rad = out[-1][3]
        print(f"      mean |h| by depth (D={out[-1][0]}): " + " ".join(f"{r:.3f}" for r in rad[::6]))
        d1 = abs(out[1][1] - out[0][1])
        d2 = abs(out[2][1] - out[1][1])
        print(f"      |h_D(root) - h_(D-4)(root)| = {d2:.2e} (previous step {d1:.2e})")
        check(d2 <= d1 * 1.05 + 1e-9 and d2 < 0.02, "Tutte root position settles in D (decreasing increments)")
    print("  The embedding of each vertex is the harmonic-measure average of the boundary angles, so its")
    print("  convergence in D is transience made visible; it is sign-blind by TC6.")

    # -----------------------------------------------------------------------------------------
    hdr("TC7  Cheeger constant (PROVED 0)")
    print("  Every multiple of q in a backward tree has only the doubling child, so the tree contains bare")
    print("  rays 3m, 6m, 12m, ...; a segment of L vertices has boundary 2, so inf |dS|/|S| = 0 on every tree")
    print("  (3n+-1, 5n+1).  Anchored expansion and br > 1 carry the transience instead.")

    print()
    print(f"ALL CHECKS PASSED  ({time.time() - T0:.1f}s)")
    mem("end")


if __name__ == "__main__":
    main()
