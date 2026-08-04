#!/usr/bin/env python3
"""AMM 12592 -- ANGLE 5: hybrid endgame algebra at the gamma* floor profile (boxeph).

Problem (*) at epoch [R,2R):  sum_{i=0}^{R-1} x^i Delta_i(x) = (1-x)^{R-1},
Delta_i = sum_k delta_{i,k} B_{d_i,k},  B_{d,k}(x) = x^{d-k}(1-x)^k, with the
Lucas box |delta_{i,k}| <= C(d_i,k) and parity delta_{i,k} == C(d_i,k) (mod 2),
d_i = floor(gamma*(R+i)) via the GS rational proxy (exact for m <= 512).

This file attacks the known failure locus of every sweep/beam (the END) by
characterizing EXACTLY which residuals the last L rows can absorb.

============================  THE ENDGAME ALGEBRA  ============================

S_L(d_a,...,d_z) := { sigma : sigma = sum_{j=0}^{L-1} x^j Delta_j', each block
admissible at its degree }.  Membership tests implemented here:

L = 1 (exact, complete).  {B_{d,k}} is a basis of Z[x]_{<=d}, triangular
against {x^j} by lowest term, so the decode is FORCED cell by cell
(base.final_decode).  sigma in S_1 iff deg <= d and the unique coordinates
pass box+parity.

DUAL FORMULA.  From x^j = x^j (x+(1-x))^{d-j} = sum_m C(d-j,m) B_{d,m):
        b_m(sigma') = sum_j s'_j C(d-j, m)          (deg sigma' <= d).
This is the O(d^2) "dual transform" used everywhere below.

L = 2, BAND REDUCTION.  sigma = A + x*B, A at degree da, B at degree db,
e := db-da >= 0.  The top cell a_da is forced ( = sigma(0) = +-1 ).  Fold it:
sigma~ = sigma - sigma(0)*(1-x)^da.  Then with beta := dual transform of
(sigma~)/x at db, and unknowns a_0..a_{da-1}:
        b_m  =  beta_m - sum_k C(e+1, m-k) a_k ,    m = 0..db.
Proof: column of a_k is the decode of B_{da,k}/x, and by Vandermonde
sum_u (-1)^u C(k,u) C(N-u,m) = C(N-k, m-k) with N = e+1+k the column collapses
to the (e+2)-wide band C(e+1, m-k).  (Verified exactly for a grid of (da,e)
in --selftest; that computational lemma is all the search relies on.)

So 2-row absorbability is a STAIRCASE integer program with bandwidth e+1:
  a_m in A_m := [-C(da,m), C(da,m)] with parity C(da,m) mod 2,
  beta_m - sum_w C(e+1,w) a_{m-w} in B_m := [-C(db,m), C(db,m)] parity likewise,
plus deg sigma <= db+1 and sigma(0) = +-1.

e = 0 (equal degrees): the band is (a_m + a_{m-1}) -- a CHAIN.  Forward
interval+parity DP:  U_m = A_m  cap  (beta_m - B_m - U_{m-1}).
COMPLETE, because (i) all parities are forced, so after normalizing interval
endpoints onto the parity progression, plain interval arithmetic is exact;
(ii) the Minkowski sum of two step-2 progressions is the full step-2
progression of the endpoint sums; (iii) constraint m touches only a_{m-1},a_m,
so forward reachability intervals are exactly the projections of the feasible
set, and backward greedy extraction (choose a_m, then a_{m-1} in
U_{m-1} cap (beta_m - B_m - a_m), nonempty by construction) always lands a
witness.  Decision + construction in O(da) interval ops after the transform.

e = 1: band (a_m + 2a_{m-1} + a_{m-2}).  Implemented as a SANDWICH:
  necessity: the same forward DP with window intervals U_{m-1}, U_{m-2}
     is a RELAXATION (true reachable sets are subsets of the tracked
     intervals), so first-gap or parity-precheck failure PROVES sigma not in
     S_2;
  sufficiency #1 ("downgrade"): solve at (da,da) by the complete chain DP,
     then lift the second block da -> db by B_{da,k} = sum_j C(1,j) B_{da+1,k+j}
     (multiplicative lift, preserves the polynomial and admissibility);
  sufficiency #2: budgeted DFS over the relaxed DP intervals with exact
     look-back constraints + final exact verification;
  sufficiency #3: steered-step grid on row a + exact decode of row b.
Everything that claims success is verified exactly (poly identity + boxes +
parity) before being believed.  A complete e=1 decision would track the exact
reachable sets of (a_{m-1}, a_m); by induction these are the lattice points of
a convex polygon on the fixed parity sublattice (the new coordinate pair
depends only on (a_m, a_{m+1}) while the discarded coordinate enters only
through a 1-D fiber-interval nonemptiness condition), so a piecewise-linear
polygon DP would close the gap -- left as the upgrade path; the sandwich plus
final exact verification is what the hunt needs.

L = 3: semi-decision both ways: necessity by the L=2 relaxation after one
forced step; sufficiency by exhaustive steered grid on the first row followed
by the L=2 solver (this is the C2 completion below).

================================  THE HUNT  ===================================

solve_endgame() = the R=64-winning l1deg beam (rank (deg, full L1), THM-3029
recipe) with the 2-row grid endgame REPLACED by S_2 machinery:
  C1: for every beam state at row R-3, absorb2_solve(sigma, d_{R-2}, d_{R-1});
  C2: for every beam state at row R-4, wide steered grid on row R-3 then
      absorb2_solve  (= L=3 sufficiency).
Optional rank 'a2': after the l1deg pre-sort, re-rank the top `prefilter`
states by (deg, S_2-deficit, L1), where the deficit is the exact total gap
mass of the interval DP against the next two profile degrees -- "distance to
the absorbable set" steering.

Search negatives are NEVER infeasibility evidence (THM-3029).  Exact integer
arithmetic only; every witness is re-verified independently before saving.

Usage:
  --selftest            algebra + DP-vs-bruteforce validation (fast)
  --r64check            S_2 on the real R=64 witness endgame + full R=64 pipeline
  --hunt --R 128 ...    the hunt (see --help)
  --verify-existing     exact re-verification of any R=128 witness JSONs present
"""
import argparse, json, os, random, sys, time
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import amm12592_r64_floor_solver_boxeph as base

RESULTS = os.path.normpath(os.path.join(HERE, "..", "05-knowledge", "results"))

# ---------------------------------------------------------------------------
# exact combinatorial kit
# ---------------------------------------------------------------------------
_PC = [[1]]
def combrow(n):
    while len(_PC) <= n:
        r = _PC[-1]
        _PC.append([1] + [r[i] + r[i + 1] for i in range(len(r) - 1)] + [1])
    return _PC[n]

def dual_transform(sp, db):
    """b_m = sum_j sp[j] * C(db-j, m); requires len(sp) <= db+1."""
    beta = [0] * (db + 1)
    for j, c in enumerate(sp):
        if c:
            row = combrow(db - j)
            for m in range(db - j + 1):
                beta[m] += c * row[m]
    return beta

def block_poly(delta, d):
    acc = [0] * (d + 1)
    for k, v in enumerate(delta):
        if v:
            t = base.qk(k)
            off = d - k
            for s in range(k + 1):
                acc[off + s] += v * t[s]
    return base.trim(acc)

def verify_pair(sigma, A, da, B, db):
    """Exact: A admissible at da, B admissible at db, A + x*B == sigma."""
    if not base.admissible(A, da) or not base.admissible(B, db):
        return False
    pa = block_poly(A, da)
    pb = block_poly(B, db)
    acc = [0] * max(len(pa), len(pb) + 1, 1)
    for j, c in enumerate(pa):
        acc[j] += c
    for j, c in enumerate(pb):
        acc[j + 1] += c
    return base.trim(acc) == base.trim(list(sigma))

def pnorm(lo, hi, r):
    if lo % 2 != r:
        lo += 1
    if hi % 2 != r:
        hi -= 1
    return lo, hi

def cdiv(a, b):
    return -((-a) // b)

# ---------------------------------------------------------------------------
# S_2 machinery
# ---------------------------------------------------------------------------
def absorb2_dp(sigma, da, db):
    """Forward interval+parity DP for sigma = A + x*B at degrees (da, db).

    Returns dict(deficit, pfails, U, beta, s0, e).  deficit == 0 means the
    necessary conditions all pass (and for e == 0 that is COMPLETE:
    feasibility holds and absorb2_extract will construct a witness).
    deficit > 0 is an exact refutation ONLY through its first gap / parity
    fail; the accumulated total is a steering score, not a certificate."""
    e = db - da
    band = combrow(e + 1)
    sigma = base.trim(list(sigma))
    out = dict(deficit=0, pfails=0, U=None, beta=None, s0=0, e=e)
    D = 0
    if not sigma:
        out["deficit"] = 1
        return out
    s0 = sigma[0]
    if abs(s0) != 1:
        D += 4 + 2 * abs(abs(s0) - 1)
        s0 = 1 if s0 >= 0 else -1
    out["s0"] = s0
    if len(sigma) - 1 > db + 1:
        D += sum(abs(c) for c in sigma[db + 2:])
        sigma = sigma[:db + 2]
    # fold forced top cell of A
    st = list(sigma) + [0] * (max(da + 1, len(sigma)) - len(sigma))
    t = base.qk(da)
    for j in range(da + 1):
        st[j] -= s0 * t[j]
    if st[0] != 0:
        D += abs(st[0])
        st[0] = 0
    sp = st[1:]
    beta = dual_transform(sp, db)
    out["beta"] = beta
    rowA, rowB = combrow(da), combrow(db)
    parA = [rowA[k] & 1 for k in range(da)]
    # exact parity precheck (necessary regardless of e)
    pfails = 0
    for m in range(db + 1):
        pa = 0
        for w in range(e + 2):
            k = m - w
            if 0 <= k <= da - 1 and (band[w] & 1) and parA[k]:
                pa ^= 1
        if ((beta[m] & 1) ^ pa) != (rowB[m] & 1):
            pfails += 1
    out["pfails"] = pfails
    D += 2 * pfails
    # forward interval DP
    U = {}
    for w in range(1, e + 2):
        U[-w] = (0, 0)
    for m in range(db + 1):
        lo = beta[m] - rowB[m]
        hi = beta[m] + rowB[m]
        for w in range(1, e + 2):
            ul, uh = U.get(m - w, (0, 0))
            lo -= band[w] * uh
            hi -= band[w] * ul
        if m <= da - 1:
            capA = rowA[m]
            nlo, nhi = max(lo, -capA), min(hi, capA)
            if nlo > nhi:
                D += nlo - nhi          # raw magnitude gap (parity separate)
                mid = min(max(-capA, lo if lo > capA else -capA if hi < -capA else 0), capA)
                # clamp to nearest admissible point of A_m toward the need
                if lo > capA:
                    mid = capA
                elif hi < -capA:
                    mid = -capA
                mlo, mhi = pnorm(mid, mid, parA[m])
                if mlo > mhi:
                    mlo = mhi = mid - 1 if mid - 1 >= -capA else mid + 1
                U[m] = (mlo, mlo)
            else:
                plo, phi = pnorm(nlo, nhi, parA[m])
                if plo > phi:
                    D += 2
                    U[m] = (nlo, nlo)
                else:
                    U[m] = (plo, phi)
        else:
            if lo > 0:
                D += lo
            elif hi < 0:
                D += -hi
    out["deficit"] = D
    out["U"] = U
    return out

def absorb2_extract(sigma, da, db, r, budget=4000):
    """Backward extraction from the DP of absorb2_dp (deficit must be 0).
    e = 0: guaranteed (complete chain).  e >= 1: budgeted DFS, exactness by
    the final verify_pair.  Returns (A, B) or None."""
    e = db - da
    band = combrow(e + 1)
    rowA, rowB = combrow(da), combrow(db)
    beta, U, s0 = r["beta"], r["U"], r["s0"]
    nodes = [0]

    def cands(m, sol):
        lo, hi = U[m]
        mp = m + e + 1
        if mp <= db:
            l2 = beta[mp] - rowB[mp]
            h2 = beta[mp] + rowB[mp]
            for w in range(0, e + 1):
                k = mp - w
                v = sol.get(k, 0) if 0 <= k <= da - 1 else 0
                l2 -= band[w] * v
                h2 -= band[w] * v
            lo, hi = max(lo, l2), min(hi, h2)
        rpar = rowA[m] & 1
        lo, hi = pnorm(lo, hi, rpar)
        if lo > hi:
            return []
        slo, shi = lo, hi
        for mp2 in range(m + 1, min(m + e, db) + 1):
            l3 = beta[mp2] - rowB[mp2]
            h3 = beta[mp2] + rowB[mp2]
            cm = band[mp2 - m]
            for w in range(0, e + 2):
                k = mp2 - w
                if k == m:
                    continue
                if k > m:
                    v = sol.get(k, 0) if 0 <= k <= da - 1 else 0
                    l3 -= band[w] * v
                    h3 -= band[w] * v
                elif 0 <= k <= da - 1:
                    ul, uh = U[k]
                    l3 -= band[w] * uh
                    h3 -= band[w] * ul
            slo = max(slo, cdiv(l3, cm))
            shi = min(shi, h3 // cm)
        slo, shi = pnorm(slo, shi, rpar)
        out = []
        def near0(a, b):
            if a > b:
                return None
            if a <= 0 <= b:
                return 0 if 0 % 2 == rpar else (1 if 1 <= b else -1)
            return a if a > 0 else b
        for v in ([near0(slo, shi), slo, shi] if slo <= shi else []) + [near0(lo, hi), lo, hi]:
            if v is not None and lo <= v <= hi and v not in out:
                out.append(v)
            if len(out) >= 4:
                break
        return out

    def dfs(m, sol):
        nodes[0] += 1
        if nodes[0] > budget:
            return None
        if m < 0:
            return dict(sol)
        for v in cands(m, sol):
            sol[m] = v
            got = dfs(m - 1, sol)
            if got is not None:
                return got
            del sol[m]
        return None

    got = dfs(da - 1, {})
    if got is None:
        return None
    A = [got[k] for k in range(da)] + [s0]
    B = []
    for m in range(db + 1):
        v = beta[m]
        for w in range(e + 2):
            k = m - w
            if 0 <= k <= da - 1:
                v -= band[w] * A[k]
        B.append(v)
    if not verify_pair(sigma, A, da, B, db):
        return None
    return A, B

STATS = {"unknown": 0, "unknown_samples": [], "methods": {}}

def _hit(m):
    STATS["methods"][m] = STATS["methods"].get(m, 0) + 1

def absorb2_solve(sigma, da, db, scan_grid=(3, 6), budget=4000):
    """Constructive S_2 membership: returns (A, B, method) or None.
    Complete for e = db-da = 0; sandwich (sound both ways, possibly
    incomplete in the middle) for e >= 1."""
    sigma = base.trim(list(sigma))
    if not sigma or abs(sigma[0]) != 1 or len(sigma) - 1 > db + 1:
        return None
    e = db - da
    if e == 0:
        r = absorb2_dp(sigma, da, db)
        if r["deficit"] != 0:
            return None                     # complete refutation
        w = absorb2_extract(sigma, da, db, r)
        if w:
            _hit("dp-e0")
            return (w[0], w[1], "dp-e0")
        return None                          # should not happen; stay safe
    # e >= 1 sandwich
    if len(sigma) - 1 <= da + 1:
        r0 = absorb2_dp(sigma, da, da)
        if r0["deficit"] == 0:
            w = absorb2_extract(sigma, da, da, r0)
            if w:
                A, B = w
                B2 = base.lift_block(B, da, db)
                if verify_pair(sigma, A, da, B2, db):
                    _hit("dp-downgrade-lift")
                    return (A, B2, "dp-downgrade-lift")
    r1 = absorb2_dp(sigma, da, db)
    if r1["deficit"] != 0:
        return None                          # sound relaxation refutation
    w = absorb2_extract(sigma, da, db, r1, budget=budget)
    if w:
        _hit("dp-dfs")
        return (w[0], w[1], "dp-dfs")
    ec, es = scan_grid
    eopts = [None] + list(range(-es, es + 1))
    for tg in product(eopts, repeat=ec):
        if tg[0] not in (1, -1):
            continue
        rr = base.step(sigma, da, tg)
        if rr is None:
            continue
        de, ns = rr
        fin = base.final_decode(ns, db)
        if fin is not None and verify_pair(sigma, de, da, fin, db):
            _hit("scan")
            return (de, fin, "scan")
    STATS["unknown"] += 1
    if len(STATS["unknown_samples"]) < 50:
        STATS["unknown_samples"].append({"sigma": sigma, "da": da, "db": db})
    return None

def absorb2_deficit(sigma, da, db):
    return absorb2_dp(sigma, da, db)["deficit"]

# ---------------------------------------------------------------------------
# the hunt: l1deg beam + S_2 endgame (C1) + steered-grid * S_2 (C2)
# ---------------------------------------------------------------------------
def solve_endgame(d, beam=400, ctrl=2, span=2, seed=None, rand_frac=0.0,
                  dedup=999, rank="l1deg", a2_from=64, prefilter=1200,
                  c2_ctrl=3, c2_span=6, skip_c2=False, log=print, ckpt=None,
                  late_from=None, late_ctrl=None, late_span=None,
                  eprune=False, ebranch=False, ecap=None):
    R = len(d)
    rng = random.Random(seed)
    states = [([], base.qpow(R - 1))]
    states_prev = None
    t0 = time.time()
    for i in range(R - 2):
        c_i, s_i = ctrl, span
        if late_from is not None and i >= late_from:
            c_i = late_ctrl if late_ctrl is not None else ctrl
            s_i = late_span if late_span is not None else span
        opts = [None] + list(range(-s_i, s_i + 1))
        nxt = []
        for acc, sig in states:
            for tg in product(opts, repeat=c_i):
                if tg[0] not in (1, -1):
                    continue
                r = base.step(sig, d[i], tg)
                if r is None:
                    continue
                de, ns = r
                if not ns or abs(ns[0]) != 1:
                    continue
                variants = [(de, ns)]
                if ebranch and de and de[0] in (1, -1):
                    # the forced k=0 cell (B_{d,0} = x^d) admits BOTH signs;
                    # flipping delta_{i,0}: v -> -v only shifts the residual's
                    # x^{d} coefficient by +2v, i.e. ns[d-1] += 2v.  This is
                    # the ONLY control surface of the E-walk.
                    v = de[0]
                    ns2 = list(ns) + [0] * (max(0, d[i] - len(ns)))
                    ns2[d[i] - 1] += 2 * v
                    ns2 = base.trim(ns2)
                    if ns2 and abs(ns2[0]) == 1:
                        de2 = list(de)
                        de2[0] = -v
                        variants.append((de2, ns2))
                lim = R - 1 - i
                if ecap is not None:
                    lim = min(lim, ecap)
                for dev, nsv in variants:
                    if eprune and abs(sum(nsv)) > lim:
                        # E-wall (exact necessity for lim = R-1-i): sigma_i(1)
                        # walks by the forced delta_{j,0} = +-1 each remaining
                        # row and must end at 0.  ecap further restricts to
                        # the empirical witness corridor (heuristic).
                        continue
                    nxt.append((acc + [dev], nsv))
        if not nxt:
            return None, f"died at row {i}", states
        nxt.sort(key=lambda st: (len(st[1]), sum(abs(x) for x in st[1])))
        seen, uniq = set(), []
        for a, sg in nxt:
            key = tuple(sg[:dedup])
            if key in seen:
                continue
            seen.add(key)
            uniq.append((a, sg))
        if rank == "a2" and i >= a2_from and i + 2 <= R - 1:
            da2, db2 = d[i + 1], d[i + 2]
            pool = uniq[:max(prefilter, beam)]
            scored = sorted(
                (((len(sg), absorb2_deficit(sg, da2, db2),
                   sum(abs(x) for x in sg)), a, sg) for a, sg in pool),
                key=lambda z: z[0])
            uniq = [(a, sg) for _, a, sg in scored] + uniq[len(pool):]
        n_top = beam if rand_frac <= 0 else max(1, int(beam * (1 - rand_frac)))
        keep = uniq[:n_top]
        if rand_frac > 0 and len(uniq) > n_top:
            pool2 = uniq[n_top:]
            keep = keep + rng.sample(pool2, min(beam - len(keep), len(pool2)))
        states = keep[:beam]
        if i == R - 4:
            states_prev = [(list(a), list(sg)) for a, sg in states]
        best = states[0][1]
        bd = absorb2_deficit(best, d[i + 1], d[i + 2]) if i + 2 <= R - 1 else -1
        log(f"  row {i:3d}: st={len(states):5d} bestdeg={len(best)-1:4d} "
            f"bestL1={sum(abs(x) for x in best)} def2={bd} "
            f"t={time.time()-t0:7.1f}s")
        if ckpt and (i % 8 == 0 or i >= R - 6):
            with open(ckpt, "w") as f:
                json.dump({"row": i, "elapsed": time.time() - t0,
                           "residuals": [sg for _, sg in states[:24]]}, f)
    # C1: complete-ish 2-row absorption on every sigma_{R-3}
    log(f"  C1: absorb2_solve on {len(states)} states at "
        f"(d[{R-2}],d[{R-1}])=({d[R-2]},{d[R-1]})")
    for acc, sig in states:
        w = absorb2_solve(sig, d[R - 2], d[R - 1])
        if w:
            return acc + [w[0], w[1]], f"SOLVED-C1({w[2]})", states
    if skip_c2 or states_prev is None:
        return None, "exhausted (C1)", states
    # C2: wide steered grid on row R-3, then S_2 on the rest
    log(f"  C2: grid ctrl={c2_ctrl} span={c2_span} on {len(states_prev)} "
        f"states at row {R-3} (d={d[R-3]})")
    eopts = [None] + list(range(-c2_span, c2_span + 1))
    tried = 0
    for acc, sig in states_prev:
        for tg in product(eopts, repeat=c2_ctrl):
            if tg[0] not in (1, -1):
                continue
            r = base.step(sig, d[R - 3], tg)
            if r is None:
                continue
            de, ns = r
            if not ns or abs(ns[0]) != 1:
                continue
            tried += 1
            w = absorb2_solve(ns, d[R - 2], d[R - 1])
            if w:
                return acc + [de, w[0], w[1]], f"SOLVED-C2({w[2]})", states
    return None, f"exhausted (C1+C2, {tried} C2 leaves)", states

# ---------------------------------------------------------------------------
# selftest: the computational lemmas + DP vs brute force
# ---------------------------------------------------------------------------
def _rand_block(dd, rng):
    return [rng.choice(range(-base.comb(dd, k) if False else -__import__('math').comb(dd, k),
                             __import__('math').comb(dd, k) + 1, 2))
            if __import__('math').comb(dd, k) % 2 == 0 else
            rng.choice(range(-__import__('math').comb(dd, k),
                             __import__('math').comb(dd, k) + 1, 2))
            for k in range(dd + 1)]

def rand_block(dd, rng):
    from math import comb
    out = []
    for k in range(dd + 1):
        cap = comb(dd, k)
        lo = -cap if (cap - (-cap)) % 2 == 0 else -cap
        vals = list(range(-cap, cap + 1, 2)) if (cap % 2 == 0) else list(range(-cap, cap + 1, 2))
        # parity must equal cap mod 2: range from -cap step 2 hits exactly those
        out.append(rng.choice(vals))
    return out

def brute_absorb2(sigma, da, db):
    """Enumerate all admissible A at da; decode the rest at db.  Exact ground
    truth for small da."""
    from math import comb
    sigma = base.trim(list(sigma))
    if not sigma or abs(sigma[0]) != 1 or len(sigma) - 1 > db + 1:
        return None
    ranges = [range(-comb(da, k), comb(da, k) + 1, 2) for k in range(da + 1)]
    for A in product(*ranges):
        if A[da] != sigma[0]:
            continue
        pa = block_poly(list(A), da)
        res = list(sigma) + [0] * (max(0, len(pa) - len(sigma)))
        for j, c in enumerate(pa):
            res[j] -= c
        if res[0] != 0:
            continue
        fin = base.final_decode(base.trim(res[1:]), db)
        if fin is not None:
            return list(A), fin
    return None

def selftest():
    from math import comb
    print("== SELFTEST: endgame algebra ==")
    rng = random.Random(20260803)
    # 1. dual formula
    ok = True
    for db in (4, 7, 11):
        for _ in range(20):
            b = rand_block(db, rng)
            sp = block_poly(b, db)
            sp = sp + [0] * (db + 1 - len(sp))
            ok &= dual_transform(sp, db) == b
    print(f"dual formula b_m = sum_j s_j C(d-j,m): {'PASS' if ok else 'FAIL'}")
    assert ok
    # 2. band identity: decode column of B_{da,k}/x at db equals C(e+1, m-k)
    ok = True
    for da, e in [(5, 0), (6, 0), (6, 1), (7, 1), (7, 2), (9, 0), (9, 1)]:
        db = da + e
        for k in range(da):
            col = block_poly([0] * k + [1], da)     # B_{da,k}
            colx = col[1:]                          # /x  (k<da => divisible)
            colx = colx + [0] * (db + 1 - len(colx))
            got = dual_transform(colx, db)
            want = [comb(e + 1, m - k) if 0 <= m - k <= e + 1 else 0
                    for m in range(db + 1)]
            if got != want:
                ok = False
                print(f"  band MISMATCH da={da} db={db} k={k}: {got} vs {want}")
    print(f"band identity C(e+1, m-k) for k < da (grid of (da,e)): "
          f"{'PASS' if ok else 'FAIL'}")
    assert ok
    # 3. e=0 completeness vs brute force at (5,5)
    da = db = 5
    agree = pos = neg = 0
    for trial in range(60):
        if trial % 2 == 0:
            A = rand_block(da, rng)
            B = rand_block(db, rng)
            pa, pb = block_poly(A, da), block_poly(B, db)
            sig = [0] * max(len(pa), len(pb) + 1)
            for j, c in enumerate(pa):
                sig[j] += c
            for j, c in enumerate(pb):
                sig[j + 1] += c
            sig = base.trim(sig)
        else:
            sig = [rng.choice((1, -1))] + [rng.randint(-6, 6)
                                           for _ in range(rng.randint(2, db + 1))]
            sig = base.trim(sig)
            if not sig or abs(sig[0]) != 1:
                continue
        bt = brute_absorb2(sig, da, db)
        got = absorb2_solve(sig, da, db)
        if (bt is None) == (got is None):
            agree += 1
            if bt is not None:
                pos += 1
            else:
                neg += 1
        else:
            print(f"  e=0 DISAGREE sigma={sig} brute={'feas' if bt else 'infeas'} "
                  f"dp={'feas' if got else 'infeas'}")
    print(f"e=0 DP vs brute force at (5,5): agree={agree} (pos={pos} neg={neg})"
          f" {'PASS' if agree == pos + neg and agree > 0 else 'FAIL'}")
    assert agree == pos + neg and pos > 0 and neg > 0
    # 4. e=1 sandwich at (5,6)
    da, db = 5, 6
    nec_ok = True
    suff = miss = tot_feas = 0
    meth = {}
    for trial in range(60):
        if trial % 2 == 0:
            A = rand_block(da, rng)
            B = rand_block(db, rng)
            pa, pb = block_poly(A, da), block_poly(B, db)
            sig = [0] * max(len(pa), len(pb) + 1)
            for j, c in enumerate(pa):
                sig[j] += c
            for j, c in enumerate(pb):
                sig[j + 1] += c
            sig = base.trim(sig)
        else:
            sig = [rng.choice((1, -1))] + [rng.randint(-8, 8)
                                           for _ in range(rng.randint(2, db + 1))]
            sig = base.trim(sig)
        if not sig or abs(sig[0]) != 1:
            continue
        bt = brute_absorb2(sig, da, db)
        r = absorb2_dp(sig, da, db) if len(sig) - 1 <= db + 1 else None
        got = absorb2_solve(sig, da, db)
        if bt is not None:
            tot_feas += 1
            if r is None or r["deficit"] != 0:
                nec_ok = False
                print(f"  e=1 NECESSITY VIOLATION sigma={sig}")
            if got is not None:
                suff += 1
                meth[got[2]] = meth.get(got[2], 0) + 1
            else:
                miss += 1
        else:
            if got is not None:
                print(f"  e=1 FALSE POSITIVE (verify_pair let a bad pair "
                      f"through?!) sigma={sig}")
                nec_ok = False
    print(f"e=1 sandwich at (5,6): feasible={tot_feas} constructed={suff} "
          f"missed-by-sandwich={miss} methods={meth} "
          f"necessity={'PASS' if nec_ok else 'FAIL'}")
    assert nec_ok and tot_feas > 0 and suff > 0
    print("== SELFTEST PASSED ==")

def r64check():
    print("== R64 CHECK: S_2 on the real witness endgame + full pipeline ==")
    R = 64
    d = base.prof(R, *base.GS, 0)
    with open(os.path.join(HERE, "amm12592_floor_witness_R64.json")) as f:
        w = json.load(f)
    blocks = w["blocks"]
    a, b = base.verify_witness(R, blocks, d)
    assert a and b
    pa = block_poly(blocks[R - 2], d[R - 2])
    pb = block_poly(blocks[R - 1], d[R - 1])
    sig = [0] * max(len(pa), len(pb) + 1)
    for j, c in enumerate(pa):
        sig[j] += c
    for j, c in enumerate(pb):
        sig[j + 1] += c
    sig = base.trim(sig)
    print(f"sigma_61 (= Delta_62 + x Delta_63): deg={len(sig)-1} "
          f"L1={sum(abs(x) for x in sig)}  degrees=({d[R-2]},{d[R-1]}) "
          f"e={d[R-1]-d[R-2]}")
    got = absorb2_solve(sig, d[R - 2], d[R - 1])
    assert got is not None, "S_2 failed on a residual KNOWN to be absorbable"
    print(f"absorb2_solve: SUCCESS via {got[2]}; pair verifies="
          f"{verify_pair(sig, got[0], d[R-2], got[1], d[R-1])}")
    nb = blocks[:R - 2] + [got[0], got[1]]
    a2, b2 = base.verify_witness(R, nb, d)
    print(f"spliced full witness verifies: admissible={a2} identity={b2}")
    assert a2 and b2
    print("full R=64 pipeline (l1deg beam 400 + S_2 completion):")
    sol, msg, _ = solve_endgame(d, beam=400, ctrl=2, span=2, rank="l1deg",
                                log=lambda s: None)
    print(f"  -> {msg}")
    assert sol is not None
    va, vb = base.verify_witness(R, sol, d)
    print(f"  pipeline witness: admissible={va} identity={vb}")
    assert va and vb
    print("== R64 CHECK PASSED ==")

def verify_existing():
    print("== EXACT RE-VERIFICATION of R=128 witness JSONs present ==")
    R = 128
    d = base.prof(R, *base.GS, 0)
    import glob
    for f in sorted(glob.glob(os.path.join(HERE, "amm12592_floor_witness_R128*.json"))):
        try:
            w = json.load(open(f))
            if w.get("R") != 128:
                continue
            blocks = w["blocks"]
            if w.get("profile") != d:
                print(f"{os.path.basename(f)}: PROFILE MISMATCH")
                continue
            a, b = base.verify_witness(R, blocks, d)
            sat = max((abs(v) / base.comb(d[i], k) if base.comb(d[i], k) else 0.0)
                      for i, de in enumerate(blocks) for k, v in enumerate(de)
                      if v)
            print(f"{os.path.basename(f)}: admissible={a} identity={b} "
                  f"max box saturation={sat:.3e}")
        except Exception as ex:
            print(f"{os.path.basename(f)}: ERROR {ex}")

def hunt(args):
    R = args.R
    d = base.prof(R, *base.GS, 0)
    tag = args.tag or f"{args.rank}b{args.beam}"
    log_path = os.path.join(RESULTS, f"amm12592_r128_endgame_beam_{tag}.log")
    ckpt = os.path.join(RESULTS, f"amm12592_r128_endgame_ckpt_{tag}.json")
    lf = open(log_path, "a", buffering=1)
    def log(s):
        print(s, flush=True)
        lf.write(s + "\n")
    log(f"HUNT R={R} floor profile, rank={args.rank} beam={args.beam} "
        f"ctrl={args.ctrl} span={args.span} a2_from={args.a2_from} "
        f"prefilter={args.prefilter} seed={args.seed} rand_frac={args.rand_frac} "
        f"c2=({args.c2_ctrl},{args.c2_span}) eprune={args.eprune} "
        f"ebranch={args.ebranch} ecap={args.ecap}")
    sol, msg, states = solve_endgame(
        d, beam=args.beam, ctrl=args.ctrl, span=args.span, seed=args.seed,
        rand_frac=args.rand_frac, dedup=args.dedup, rank=args.rank,
        a2_from=args.a2_from, prefilter=args.prefilter,
        c2_ctrl=args.c2_ctrl, c2_span=args.c2_span, skip_c2=args.skip_c2,
        log=log, ckpt=ckpt, late_from=args.late_from,
        late_ctrl=args.late_ctrl, late_span=args.late_span,
        eprune=args.eprune, ebranch=args.ebranch, ecap=args.ecap)
    log(f"  -> {msg}")
    log(f"  absorb2 method hits: {STATS['methods']}  "
        f"unknown(relaxed-pass, unconstructed): {STATS['unknown']}")
    if states:
        dump = os.path.join(RESULTS, f"amm12592_r128_endgame_states_{tag}.json")
        with open(dump, "w") as f:
            json.dump({"tag": tag, "msg": msg,
                       "residuals": [sg for _, sg in states],
                       "top_accs": [acc for acc, _ in states[:5]]}, f)
        log(f"  final residuals dumped: {dump}")
    if STATS["unknown_samples"]:
        up = os.path.join(RESULTS, f"amm12592_r128_endgame_unknowns_{tag}.json")
        with open(up, "w") as f:
            json.dump(STATS["unknown_samples"], f)
        log(f"  unknown residual samples: {up}")
    if sol is None:
        return 1
    va, vb = base.verify_witness(R, sol, d)
    log(f"  EXACT VERIFY: admissible={va} identity={vb}")
    assert va and vb
    out = args.out or os.path.join(HERE, "amm12592_floor_witness_R128_endgame.json")
    with open(out, "w") as f:
        json.dump({"R": R, "profile": d, "blocks": sol, "verified": True,
                   "source_label": f"gamma* floor (ANGLE 5 endgame algebra: "
                                   f"l1deg beam + S_2 interval-decode "
                                   f"completion; {msg})",
                   "search": {"rank": args.rank, "beam": args.beam,
                              "ctrl": args.ctrl, "span": args.span,
                              "seed": args.seed, "rand_frac": args.rand_frac,
                              "a2_from": args.a2_from,
                              "c2": [args.c2_ctrl, args.c2_span]}}, f)
    log(f"  WITNESS WRITTEN: {out}")
    return 0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--r64check", action="store_true")
    ap.add_argument("--verify-existing", action="store_true", dest="verify_existing")
    ap.add_argument("--hunt", action="store_true")
    ap.add_argument("--R", type=int, default=128)
    ap.add_argument("--beam", type=int, default=400)
    ap.add_argument("--ctrl", type=int, default=2)
    ap.add_argument("--span", type=int, default=2)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--rand-frac", type=float, default=0.0, dest="rand_frac")
    ap.add_argument("--dedup", type=int, default=999)
    ap.add_argument("--rank", default="l1deg", choices=["l1deg", "a2"])
    ap.add_argument("--a2-from", type=int, default=64, dest="a2_from")
    ap.add_argument("--prefilter", type=int, default=1200)
    ap.add_argument("--c2-ctrl", type=int, default=3, dest="c2_ctrl")
    ap.add_argument("--c2-span", type=int, default=6, dest="c2_span")
    ap.add_argument("--skip-c2", action="store_true", dest="skip_c2")
    ap.add_argument("--late-from", type=int, default=None, dest="late_from")
    ap.add_argument("--late-ctrl", type=int, default=None, dest="late_ctrl")
    ap.add_argument("--late-span", type=int, default=None, dest="late_span")
    ap.add_argument("--eprune", action="store_true")
    ap.add_argument("--ebranch", action="store_true")
    ap.add_argument("--ecap", type=int, default=None)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    rc = 0
    if args.selftest:
        selftest()
    if args.r64check:
        r64check()
    if args.verify_existing:
        verify_existing()
    if args.hunt:
        rc = hunt(args)
    sys.exit(rc)

if __name__ == "__main__":
    main()
