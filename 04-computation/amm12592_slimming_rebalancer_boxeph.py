"""AMM 12592 / epoch-closure ANGLE 2: split-jump chip-firing rebalancer + reseed.

boxeph 2026-08-03 (post THM-3302).  Problem (*) at epoch [R,2R), gamma* floor
profile d_i = floor(gamma*(R+i)); R = 8,16,32,64 closed; target R = 128.

CELL <-> MONOMIAL BOOKKEEPING (exact).  Row i cell k carries B_{d_i,k} at
shift i, i.e. the monomial  x^z (1-x)^o  with  z = i + d_i - k,  o = k.
Every cell of row i sits on the anti-diagonal z + o = s_i := i + d_i, and
s_i is STRICTLY increasing (d_i nondecreasing), so each monomial (z,o) has
AT MOST ONE host cell: hosting is injective.  The Pascal split-jump
    x^z (1-x)^o = x^z (1-x)^{o+1} + x^{z+1} (1-x)^o
iterates to the exact r-step expansion (multiply by ((1-x)+x)^r = 1):
    x^z (1-x)^k = sum_{t=0}^r C(r,t) x^{z+t} (1-x)^{k+r-t},
and ALL terms on the right live in the single row i2 with s_{i2} = s_i + r.
Hence the complete split-jump move set between rows is

    MOVE(i,k -> i2, amt):   delta[i][k]      -= amt
                            delta[i2][k+r-t] += amt * C(r,t)   (t = 0..r),
    r = s_{i2} - s_i,  legal iff  i < i2 <= z = i + d_i - k  (equivalently
    k + r <= d_{i2}); preserves the epoch identity for ANY integer amt and
    preserves the parity condition delta == binom (mod 2) for EVEN amt.
Consequences derived and used here:
  * top cells k = d_i (z = i) admit NO move in either direction: the forced
    +-1 spine is rigid, exactly the Steer chain of the doubling map;
  * cells k = d_i - 1 can only move one row (i2 = i+1 when s-step is 1);
  * bulk cells (small k) can teleport as far as row z = i + d_i - k.

CONTENTS
  1. moves engine + selftest (random even moves preserve identity/parity);
  2. metrics: interior saturation (cap > 1 cells), max cell, backbone c-mass;
  3. greedy chip-firing slimmer (objective: max interior saturation, then
     backbone deviation), admissibility preserved at every step;
  4. carve (one-level, D_j cached) + generalized carry sweep over arbitrary
     E rows, with OPTIONAL IN-SWEEP TELEPORT REPAIR: when the leftover cell
     |P[k]| exceeds a bound, the even excess is re-expressed forward into
     E[i2] (i2 <= z - buffer) by the exact r-step expansion -- this is the
     rebalancer running DURING the sweep;
  5. drivers: --selftest --metrics --slim --sweep --sweep32x4 --hunt-slim.

Exact arithmetic only (python int / Fraction); floats only in reporting.
"""
import argparse, json, os, sys, time
from fractions import Fraction
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import amm12592_r64_floor_solver_boxeph as S   # import-safe engine

GS = S.GS

def prof(R, D0=0):
    return [(GS[0] * (R + i)) // GS[1] + D0 for i in range(R)]

def conv(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        if u:
            for j, v in enumerate(b):
                r[i + j] += u * v
    return r

lift = S.lift_block
admissible = S.admissible

def ballot(d):
    return [comb(d - 1, k) - (comb(d - 1, k - 1) if k else 0) for k in range(d + 1)]

def block_poly(delta, d):
    poly = [0] * (d + 1)
    for k, v in enumerate(delta):
        if v:
            t = S.qk(k)
            off = d - k
            for s2, c in enumerate(t):
                poly[off + s2] += v * c
    return poly

def rows_poly(rows, degs, R2):
    acc = [0] * (R2 + max(degs) + 2)
    for j, cells in enumerate(rows):
        for t, c in enumerate(block_poly(cells, degs[j])):
            acc[j + t] += c
    return S.trim(acc)

def rows_identity(rows, degs, R2):
    return rows_poly(rows, degs, R2) == S.trim(S.qpow(R2 - 1))

# ---------------------------------------------------------------------------
# 1. moves engine
# ---------------------------------------------------------------------------
def move_ok(d, i, k, i2):
    R = len(d)
    return (0 <= i < i2 < R and 0 <= k <= d[i] and i2 <= i + d[i] - k)

def apply_move(B, d, i, k, i2, amt):
    """transfer amt of monomial x^{i+d_i-k}(1-x)^k from cell (i,k) into its
    r-step expansion in row i2.  Identity-preserving for any integer amt;
    parity-preserving for even amt."""
    assert move_ok(d, i, k, i2), (i, k, i2)
    r = (i2 + d[i2]) - (i + d[i])
    B[i][k] -= amt
    for t in range(r + 1):
        B[i2][k + r - t] += amt * comb(r, t)

def selftest(seed=20260803):
    import random
    rng = random.Random(seed)
    W = {w['R']: w for w in json.load(open(os.path.join(
        HERE, 'amm12592_floor_witnesses_R8_R16_R32.json')))}
    ok_all = True
    for R in (8, 16, 32):
        d = W[R]['profile']
        B = [list(r) for r in W[R]['blocks']]
        B0 = [list(r) for r in B]
        n = 0
        for _ in range(4000):
            i = rng.randrange(R - 1)
            k = rng.randrange(d[i] + 1)
            z = i + d[i] - k
            hi = min(R - 1, z)
            if hi <= i:
                continue
            i2 = rng.randrange(i + 1, hi + 1)
            amt = 2 * rng.randint(-3, 3)
            if amt == 0:
                continue
            apply_move(B, d, i, k, i2, amt)
            n += 1
        ident = rows_identity(B, d, R)
        par = all((v - comb(d[i], k)) % 2 == 0 for i in range(R)
                  for k, v in enumerate(B[i]))
        # exact reversibility on a fresh copy
        Br = [list(r) for r in B0]
        apply_move(Br, d, 1, 0, min(R - 1, 1 + d[1]), 6)
        apply_move(Br, d, 1, 0, min(R - 1, 1 + d[1]), -6)
        rev = (Br == B0)
        print(f'  selftest R={R}: {n} random even moves; identity={ident} '
              f'parity={par} reversible={rev}')
        ok_all &= ident and par and rev
    print('  SELFTEST', 'PASS' if ok_all else 'FAIL')
    return ok_all

# ---------------------------------------------------------------------------
# 2. metrics
# ---------------------------------------------------------------------------
def metrics(R, d, B):
    sat2 = Fraction(0)
    where = None
    maxcell = 0
    nfull = 0
    for i in range(R):
        for k, v in enumerate(B[i]):
            cap = comb(d[i], k)
            if cap > 1:
                s = Fraction(abs(v), cap)
                if s > sat2:
                    sat2, where = s, (i, k)
                if abs(v) == cap:
                    nfull += 1
            maxcell = max(maxcell, abs(v))
    midc = 0
    for i in range(1, R - 1):
        bal = ballot(d[i])
        midc = max(midc, max(abs(a - b) for a, b in zip(B[i], bal)))
    return {'sat2': float(sat2), 'sat2_at': where, 'maxcell': maxcell,
            'nfull_capgt1': nfull, 'midc': midc}

def print_metrics(tag, R, d, B):
    m = metrics(R, d, B)
    print(f'  {tag}: interior sat = {m["sat2"]:.6g} at {m["sat2_at"]}, '
          f'#full(cap>1) = {m["nfull_capgt1"]}, max|cell| = {m["maxcell"]:.4g}, '
          f'mid max|c| = {m["midc"]:.4g}')
    return m

# ---------------------------------------------------------------------------
# 3. greedy chip-firing slimmer
# ---------------------------------------------------------------------------
def slim_greedy(R, d, B, budget_s=600, max_iter=200000, topn=8, verbose=True):
    """Greedy descent on (max interior saturation, backbone deviation) using
    only legal even split-jump moves; admissibility preserved throughout."""
    B = [list(r) for r in B]
    caps = [[comb(d[i], k) for k in range(d[i] + 1)] for i in range(R)]

    def satcells():
        out = []
        for i in range(R):
            for k, v in enumerate(B[i]):
                if caps[i][k] > 1 and v:
                    out.append((Fraction(abs(v), caps[i][k]), i, k))
        out.sort(reverse=True)
        return out

    t0 = time.time()
    applied = 0
    stall = 0
    log = []
    while time.time() - t0 < budget_s and applied < max_iter:
        cells = satcells()
        if not cells:
            break
        gsat = cells[0][0]
        best = None   # (new_local_sat, -amt, descriptor)
        for (s0, i, k) in cells[:topn]:
            v = B[i][k]
            sg = 1 if v > 0 else -1
            z = i + d[i] - k
            # DOWN moves: (i,k) -> row i2, amt = sg*2m
            for i2 in range(i + 1, min(R - 1, z) + 1):
                r = (i2 + d[i2]) - (i + d[i])
                # m bounded by source (do not overshoot 0) and each target
                # staying strictly below s0 and inside its box
                mmax = abs(v) // 2
                okm = mmax
                for t in range(r + 1):
                    c2 = comb(r, t)
                    kk = k + r - t
                    cap2 = caps[i2][kk]
                    cur = B[i2][kk]
                    # target after: cur + sg*2m*c2 ; require |.| <= cap2 and
                    # Fraction(|.|,cap2) < s0 (if cap2>1) -- solve for m
                    lim_hard = cap2
                    lim_soft = (s0.numerator * cap2 - 1) // s0.denominator \
                        if cap2 > 1 else cap2
                    lim = min(lim_hard, lim_soft)
                    # |cur + sg*2m*c2| <= lim  =>  m <= (lim - sg*cur)/ (2 c2)
                    room = lim - sg * cur
                    okm = min(okm, room // (2 * c2) if room >= 0 else -1)
                    if okm <= 0:
                        break
                if okm <= 0:
                    continue
                m = okm
                # local objective after move
                news = Fraction(abs(v) - 2 * m, caps[i][k])
                for t in range(r + 1):
                    kk = k + r - t
                    if caps[i2][kk] > 1:
                        news = max(news, Fraction(
                            abs(B[i2][kk] + sg * 2 * m * comb(r, t)),
                            caps[i2][kk]))
                if news < s0 and (best is None or (news, -m) < best[:2]):
                    best = (news, -m, ('dn', i, k, i2, sg * 2 * m))
            # UP moves: pattern in row i pulled back to (i0,k0)
            for i0 in range(max(0, i - 12), i):
                r0 = (i + d[i]) - (i0 + d[i0])
                for k0 in range(max(0, k - r0), min(k, d[i0] - (i - i0)) + 1):
                    if not move_ok(d, i0, k0, i):
                        continue
                    tstar = k0 + r0 - k
                    if not (0 <= tstar <= r0):
                        continue
                    cstar = comb(r0, tstar)
                    # reverse move: apply_move(i0,k0,i, -sg*2m) adds +sg*2m at
                    # (i0,k0) and subtracts sg*2m*C(r0,t) from row i cells
                    mmax = abs(v) // (2 * cstar)
                    okm = mmax
                    cap0 = caps[i0][k0]
                    lim0 = min(cap0, (gsat.numerator * cap0) //
                               gsat.denominator if cap0 > 1 else cap0)
                    room0 = lim0 - sg * B[i0][k0]
                    okm = min(okm, room0 // 2 if room0 >= 0 else -1)
                    if okm <= 0:
                        continue
                    for t in range(r0 + 1):
                        kk = k0 + r0 - t
                        cap2 = caps[i][kk]
                        lim = min(cap2, (s0.numerator * cap2 - 1) //
                                  s0.denominator if cap2 > 1 else cap2)
                        c2 = comb(r0, t)
                        room = lim + sg * B[i][kk]   # |cur - sg*2m*c2| <= lim
                        okm = min(okm, room // (2 * c2) if room >= 0 else -1)
                        if okm <= 0:
                            break
                    if okm <= 0:
                        continue
                    m = okm
                    news = Fraction(0)
                    if cap0 > 1:
                        news = Fraction(abs(B[i0][k0] + sg * 2 * m), cap0)
                    for t in range(r0 + 1):
                        kk = k0 + r0 - t
                        if caps[i][kk] > 1:
                            news = max(news, Fraction(
                                abs(B[i][kk] - sg * 2 * m * comb(r0, t)),
                                caps[i][kk]))
                    if news < s0 and (best is None or (news, -m) < best[:2]):
                        best = (news, -m, ('up', i0, k0, i, -sg * 2 * m))
            if best is not None:
                break   # repair the worst repairable cell first
        if best is None:
            stall += 1
            break
        _, _, mv = best
        _, a, b, c, amt = mv
        apply_move(B, d, a, b, c, amt)
        applied += 1
        if verbose and applied % 200 == 0:
            cs = satcells()
            print(f'    slim iter {applied}: worst sat = {float(cs[0][0]):.6g} '
                  f'at ({cs[0][1]},{cs[0][2]})  t={time.time()-t0:.0f}s',
                  flush=True)
            log.append((applied, float(cs[0][0])))
    return B, applied, log

# ---------------------------------------------------------------------------
# 4. carve + generalized sweep with teleport repair
# ---------------------------------------------------------------------------
def carve(Rs, dsrc, blocks, dt):
    """one-level differenced carve E_j = D_j - D_{j-1} lifted to dt (D_j cached)."""
    R2 = 2 * Rs
    Dm, Ddeg = [None] * (R2 - 1), [0] * (R2 - 1)
    for j in range(R2 - 1):
        lo, hi = max(0, j - Rs + 1), min(j, Rs - 1)
        ddm = max(dsrc[i] + dsrc[j - i] for i in range(lo, hi + 1))
        cells = [0] * (ddm + 1)
        for i in range(lo, hi + 1):
            pr = conv(blocks[i], blocks[j - i])
            dd = dsrc[i] + dsrc[j - i]
            pl = lift(pr, dd, ddm) if dd < ddm else pr
            for k, v in enumerate(pl):
                cells[k] += v
        Dm[j], Ddeg[j] = cells, ddm
    E = []
    for j in range(R2):
        dj = dt[j]
        cells = [0] * (dj + 1)
        if j <= R2 - 2:
            assert Ddeg[j] <= dj, (j, Ddeg[j], dj)
            src = lift(Dm[j], Ddeg[j], dj) if Ddeg[j] < dj else Dm[j]
            for k, v in enumerate(src):
                cells[k] += v
        if j >= 1:
            assert Ddeg[j - 1] <= dj
            src = lift(Dm[j - 1], Ddeg[j - 1], dj) if Ddeg[j - 1] < dj else Dm[j - 1]
            for k, v in enumerate(src):
                cells[k] -= v
        E.append(cells)
    return E

def rows_overshoot(dt, E):
    ov = Fraction(0)
    at = None
    for j, cells in enumerate(E):
        for k, v in enumerate(cells):
            s = Fraction(abs(v), comb(dt[j], k))
            if s > ov:
                ov, at = s, (j, k)
    return ov, at

def pclamp(v, cap):
    w = max(-cap, min(cap, v))
    if (w - cap) % 2:
        w = w - 1 if w - 1 >= -cap else w + 1
    return w

class Steer:
    """tau-aware top-boundary prescription (verbatim math from the doubling
    module); memo cleared by the caller whenever E rows change (teleport)."""
    def __init__(self, E, dt, R2, L):
        self.E, self.dt, self.R2, self.L = E, dt, R2, L
        self.memo = {}
    def rho(self, j, m, depth):
        if (j, m) in self.memo:
            return self.memo[(j, m)]
        E, dt = self.E, self.dt
        if j + 1 >= self.R2:
            val = 0
        else:
            d = dt[j]; dn = dt[j + 1]; r = dn - (d - 1)
            fixed = 0
            for u in range(max(0, m - r), m):
                fixed += comb(r, m - u) * self.rho(j, u, depth)
            base = E[j + 1][dn - m] + fixed
            if m == 0:
                val = min((1 - base, -1 - base), key=abs)
            elif depth < self.L:
                b = 0 if comb(dn, dn - m) % 2 == 0 else 1
                val = (self.rho(j + 1, m - 1, depth + 1) + b) - base
            else:
                val = pclamp(base, comb(dn, dn - m)) - base
        self.memo[(j, m)] = val
        return val

def sweep_E(R2, dt, E, W=10, L=10, teleport=False, tele_bound_shift=0,
            tele_frac=Fraction(1, 4), tele_buffer=8, verbose=False,
            debug_fallback=None):
    """generalized carry sweep over arbitrary rows E (identity assumed:
    sum_j x^j E_j == q^{R2-1}).  With teleport=True, leftover bulk cells
    beyond ~ box scale are re-expressed forward into E via exact split-jump
    expansions (identity preserved; parity not needed on E/P)."""
    E = [list(r) for r in E]
    st = Steer(E, dt, R2, L)
    K, dK = None, None
    sol, nfall, maxcarry, ntele, telemass_log = [], 0, 0, 0, 0
    for j in range(R2):
        d = dt[j]
        X = list(E[j])
        if K is not None:
            if dK > d:
                return None, f'carry degree {dK} > row degree {d} at row {j}', nfall, maxcarry, ntele
            KL = lift(K, dK, d) if dK < d else list(K)
            for k, v in enumerate(KL):
                X[k] += v
        if abs(X[d]) != 1:
            return None, f'row {j}: forced top cell = {X[d]}', nfall, maxcarry, ntele
        cells = [0] * (d + 1)
        cells[d] = X[d]
        P = [0] * d
        presc = {}
        if j + 1 < R2:
            for m in range(min(W, d)):
                presc[d - 1 - m] = st.rho(j, m, 0)
        for k in range(d - 1, -1, -1):
            cap = comb(d, k)
            if k in presc:
                e = X[k] - presc[k]
                if abs(e) <= cap and (e - cap) % 2 == 0:
                    cells[k] = e
                    P[k] = presc[k]
                    continue
                nfall += 1
                if debug_fallback is not None:
                    debug_fallback.append(
                        (j, k, d - k, 'parity' if abs(e) <= cap else 'mag',
                         e, cap))
            v = pclamp(X[k], cap)
            cells[k] = v
            P[k] = X[k] - v
        sol.append(cells)
        if teleport and j + 1 < R2:
            moved = False
            for k in range(d):
                b0 = comb(d - 1, k)
                bound = (b0 << tele_bound_shift if tele_bound_shift >= 0
                         else max(1, b0 >> (-tele_bound_shift)))
                if abs(P[k]) <= bound:
                    continue
                keep = max(-bound, min(bound, P[k]))
                if (P[k] - keep) % 2:
                    keep += 1 if keep < P[k] else -1
                A = P[k] - keep           # even excess to teleport
                z = j + d - k             # monomial exponent of x
                i2hi = min(R2 - 1, z - tele_buffer)
                for i2 in range(j + 2, i2hi + 1):   # NEAREST first: far
                    if not A:                       # teleports pay 2^r mass
                        break
                    r = (i2 + dt[i2]) - (j + d)
                    if r < 1 or r > 24 or k + r > dt[i2]:
                        continue
                    # chunk limited by tele_frac of target boxes
                    chunk = abs(A) // 2
                    sg = 1 if A > 0 else -1
                    for t in range(r + 1):
                        kk = k + r - t
                        cap2 = comb(dt[i2], kk)
                        lim = (tele_frac.numerator * cap2) // tele_frac.denominator
                        room = lim - sg * E[i2][kk]
                        c2 = comb(r, t)
                        chunk = min(chunk, room // (2 * c2) if room > 0 else 0)
                        if chunk <= 0:
                            break
                    if chunk > 0:
                        amt = sg * 2 * chunk
                        for t in range(r + 1):
                            E[i2][k + r - t] += amt * comb(r, t)
                        A -= amt
                        ntele += 1
                        telemass_log = max(telemass_log, abs(amt))
                        moved = True
                P[k] = keep + A           # whatever could not be teleported
            if moved:
                st.memo.clear()
        if P:
            maxcarry = max(maxcarry, max(abs(x) for x in P))
        if verbose and j % 16 == 0:
            top = max((abs(x) for x in P), default=0)
            print(f'    sweep row {j:3d}: nfall={nfall} maxP={top:.3g} ntele={ntele}',
                  flush=True)
        K, dK = P, d - 1
        if j == R2 - 1 and any(P):
            return None, 'final carry nonzero', nfall, maxcarry, ntele
    return sol, 'CLOSED', nfall, maxcarry, ntele

def band_vacuum(R2, dt, E, M0=40, spread=3, keep=Fraction(1, 1),
                verbose=False):
    """Clear the near-top band of the carve rows E by exact split-jump
    UP-moves: excess in cell (j, dt_j - m) beyond its box is gathered into
    earlier rows' deeper cells (i0, k-r), i0 = j-1..j-spread, which pushes it
    DEEPER below the top (depth grows by j-i0) where capacities are larger;
    the side splash -A*C(r,t) lands in deeper cells of row j and is handled
    by the later depth passes.  Any integer amounts are legal on E rows
    (parity is the sweep's job).  MEASURED NECESSITY: the 32->64 carve that
    closed had top cells +-1..3 and band overshoot <= 1.5, while the raw
    64->128 carve has the whole band at 11-21x capacity -- the Steer chain
    cannot thread that.  Identity is preserved move-by-move (exact)."""
    E = [list(r) for r in E]
    stuck = 0
    stuckmass = 0
    for m in range(M0 + 1):
        for j in range(R2):
            k = dt[j] - m
            if k < 0:
                continue
            cap = comb(dt[j], k)
            lim = max(1, (keep.numerator * cap) // keep.denominator)
            v = E[j][k]
            kp = max(-lim, min(lim, v))
            A = v - kp
            if A == 0:
                continue
            for i0 in range(j - 1, max(-1, j - 1 - spread), -1):
                if A == 0:
                    break
                r = (j + dt[j]) - (i0 + dt[i0])
                k0 = k - r
                if k0 < 0 or not move_ok(dt, i0, k0, j):
                    continue
                cap0 = comb(dt[i0], k0)
                lim0 = max(1, (keep.numerator * cap0) // keep.denominator)
                sg = 1 if A > 0 else -1
                room = lim0 - sg * E[i0][k0]
                a = min(abs(A), room) if room > 0 else 0
                if a <= 0:
                    continue
                # up-move: E[i0][k0] += sg*a ; E[j][k0+r-t] -= sg*a*C(r,t)
                apply_move(E, dt, i0, k0, j, -sg * a)
                A -= sg * a
            if A:
                stuck += 1
                stuckmass = max(stuckmass, abs(A))
        if verbose:
            ov = Fraction(0)
            for j in range(R2):
                k = dt[j] - m
                if k >= 0:
                    ov = max(ov, Fraction(abs(E[j][k]), comb(dt[j], k)))
            print(f'    bandvac depth {m:2d}: residual band overshoot '
                  f'{float(ov):.3g} stuck={stuck}', flush=True)
    return E, stuck, stuckmass

def vacuum_fix(R2, dt, E, M0=45, passes=40, verbose=False):
    """iterated r<=2 band vacuum: excess at (j, dt_j - m) beyond its box is
    gathered ONLY into row j-1 (r = s_j - s_{j-1} in {1,2}; splash at most
    2A into (j,k-1) + A into (j,k-2), re-cleaned by later depth passes /
    iterations).  Runs to fixpoint or `passes`.  Identity exact throughout."""
    E = [list(r) for r in E]
    for p in range(passes):
        moved = 0
        for m in range(M0 + 1):
            for j in range(1, R2):
                k = dt[j] - m
                if k < 0:
                    continue
                cap = comb(dt[j], k)
                v = E[j][k]
                kp = max(-cap, min(cap, v))
                A = v - kp
                if A == 0:
                    continue
                i0 = j - 1
                r = (j + dt[j]) - (i0 + dt[i0])
                k0 = k - r
                if k0 < 0 or not move_ok(dt, i0, k0, j):
                    continue
                cap0 = comb(dt[i0], k0)
                sg = 1 if A > 0 else -1
                room = cap0 - sg * E[i0][k0]
                a = min(abs(A), room) if room > 0 else 0
                if a <= 0:
                    continue
                apply_move(E, dt, i0, k0, j, -sg * a)
                moved += 1
        ov, at = rows_overshoot(dt, E)
        if verbose:
            print(f'    vacuum_fix pass {p:2d}: moves={moved} '
                  f'overshoot={float(ov):.4g} at {at}', flush=True)
        if moved == 0:
            break
    return E

def sweep_E2(R2, dt, E, tele_rmax=24, tele_frac=Fraction(1, 8),
             tele_buffer=8, debug=None):
    """Engine 2: deterministic sweep with an EXACT 2-row top planner.
    The old Steer rho chain assumed its own prescriptions hold; with
    fallbacks the 'fixed' terms are wrong and tops break by +-2/4 (measured).
    Here everything is computed from ACTUAL state at emission time:
      next top      = E[j+1][dn]   + P[d-1]                      (exact)
      next (dn-1)   = E[j+1][dn-1] + P[d-2] + r2 P[d-1]          (exact)
      next (dn-2)   = E[j+1][dn-2] + P[d-3] + r2 P[d-2] + C(r2,2) P[d-1]
    with r2 = dn - (d-1).  The planner enumerates the small grid of carry
    tops (P[d-1] targeting +-1 next top; |P[d-2]|,|P[d-3]| <= 3), requires
    exact emission feasibility at cells d-1..d-3 (box + parity), and scores
    1-row-lookahead plannability of row j+1.  All other cells pclamp; bulk
    leftovers beyond box scale teleport NEAREST-FIRST (r <= tele_rmax --
    far teleports multiply mass by 2^r and are always net losers)."""
    E = [list(r) for r in E]
    K, dK = None, None
    sol, nfall, maxcarry, ntele = [], 0, 0, 0
    for j in range(R2):
        d = dt[j]
        X = list(E[j])
        if K is not None:
            if dK > d:
                return None, f'carry degree {dK} > {d} at row {j}', nfall, maxcarry, ntele
            KL = lift(K, dK, d) if dK < d else list(K)
            for k, v in enumerate(KL):
                X[k] += v
        if abs(X[d]) != 1:
            return None, f'row {j}: forced top cell = {X[d]}', nfall, maxcarry, ntele
        cells = [0] * (d + 1)
        cells[d] = X[d]
        P = [0] * d
        planned = {}
        if j + 1 < R2 and d >= 3:
            dn = dt[j + 1]
            r2 = dn - (d - 1)
            En = E[j + 1]
            Enn_top = E[j + 2][dt[j + 2]] if j + 2 < R2 else None
            best = None
            for s in (X[d], -X[d]):        # prefer keeping the sign pattern
                p1 = s - En[dn]
                e1 = X[d - 1] - p1
                c1 = comb(d, d - 1)
                if abs(e1) > c1 or (e1 - c1) % 2:
                    continue
                for p2 in (0, 1, -1, 2, -2, 3, -3):
                    e2 = X[d - 2] - p2
                    c2 = comb(d, d - 2)
                    if abs(e2) > c2 or (e2 - c2) % 2:
                        continue
                    for p3 in (0, 1, -1, 2, -2, 3, -3):
                        e3 = X[d - 3] - p3
                        c3 = comb(d, d - 3)
                        if abs(e3) > c3 or (e3 - c3) % 2:
                            continue
                        # 1-row lookahead: can row j+1 be planned?
                        ok2 = True
                        if Enn_top is not None:
                            Xn1 = En[dn - 1] + p2 + r2 * p1
                            cn1 = dn
                            feas = False
                            for s2 in (1, -1):
                                q1 = s2 - Enn_top
                                en1 = Xn1 - q1
                                if abs(en1) <= cn1 and (en1 - cn1) % 2 == 0:
                                    feas = True
                                    break
                            ok2 = feas
                        cost = (0 if ok2 else 1, abs(p2) + abs(p3),
                                abs(e1) + abs(e2) + abs(e3))
                        if best is None or cost < best[0]:
                            best = (cost, p1, p2, p3)
                    if best and best[0][0] == 0 and best[0][1] == 0:
                        break
                if best and best[0][0] == 0 and best[0][1] == 0:
                    break
            if best is not None:
                _, p1, p2, p3 = best
                planned = {d - 1: p1, d - 2: p2, d - 3: p3}
                if best[0][0] != 0:
                    nfall += 1
                    if debug is not None:
                        debug.append((j, 'lookahead-infeasible'))
            else:
                nfall += 1
                if debug is not None:
                    debug.append((j, f'no feasible top plan, X top3 = '
                                  f'{X[d-1:d+1]}'))
        for k in range(d - 1, -1, -1):
            cap = comb(d, k)
            if k in planned:
                cells[k] = X[k] - planned[k]
                P[k] = planned[k]
                continue
            v = pclamp(X[k], cap)
            cells[k] = v
            P[k] = X[k] - v
        sol.append(cells)
        if j + 1 < R2:
            moved = False
            for k in range(d - 3 if d >= 3 else d):
                b0 = comb(d - 1, k) if k <= d - 1 else 1
                bound = max(1, b0 >> 2)
                if abs(P[k]) <= bound:
                    continue
                keep = max(-bound, min(bound, P[k]))
                if (P[k] - keep) % 2:
                    keep += 1 if keep < P[k] else -1
                A = P[k] - keep
                z = j + d - k
                for i2 in range(j + 2, min(R2 - 1, z - tele_buffer) + 1):
                    if A == 0:
                        break
                    r = (i2 + dt[i2]) - (j + d)
                    if r < 1 or r > tele_rmax or k + r > dt[i2]:
                        continue
                    chunk = abs(A) // 2
                    sg = 1 if A > 0 else -1
                    for t in range(r + 1):
                        kk = k + r - t
                        cap2 = comb(dt[i2], kk)
                        lim = (tele_frac.numerator * cap2) // tele_frac.denominator
                        room = lim - sg * E[i2][kk]
                        c2 = comb(r, t)
                        chunk = min(chunk, room // (2 * c2) if room > 0 else 0)
                        if chunk <= 0:
                            break
                    if chunk > 0:
                        amt = sg * 2 * chunk
                        for t in range(r + 1):
                            E[i2][k + r - t] += amt * comb(r, t)
                        A -= amt
                        ntele += 1
                        moved = True
                P[k] = keep + A
        if P:
            maxcarry = max(maxcarry, max(abs(x) for x in P))
        K, dK = P, d - 1
        if j == R2 - 1 and any(P):
            return None, 'final carry nonzero', nfall, maxcarry, ntele
    return sol, 'CLOSED', nfall, maxcarry, ntele

def _presc_with_override(E, dt, R2, j, W, L, override):
    """Steer prescriptions for row j, with the m=0 target sign of row j+1
    forced to `override` (+1/-1); deeper rows use the canonical min-abs
    choice.  Fresh memo (exact recomputation each call)."""
    memo = {}
    def rho(jj, m, depth):
        if (jj, m) in memo:
            return memo[(jj, m)]
        if jj + 1 >= R2:
            val = 0
        else:
            d = dt[jj]; dn = dt[jj + 1]; r = dn - (d - 1)
            fixed = 0
            for u in range(max(0, m - r), m):
                fixed += comb(r, m - u) * rho(jj, u, depth)
            base = E[jj + 1][dn - m] + fixed
            if m == 0:
                if jj == j and override is not None:
                    val = override - base
                else:
                    val = min((1 - base, -1 - base), key=abs)
            elif depth < L:
                b = 0 if comb(dn, dn - m) % 2 == 0 else 1
                val = (rho(jj + 1, m - 1, depth + 1) + b) - base
            else:
                val = pclamp(base, comb(dn, dn - m)) - base
        memo[(jj, m)] = val
        return val
    d = dt[j]
    return {d - 1 - m: rho(j, m, 0) for m in range(min(W, d))}

def sweep_beam(R2, dt, E, W=10, L=8, beam=64, verbose=False):
    """Beam sweep: branch on the +-1 top-target of each next row (the only
    self-consistent dial in the Steer chain -- E edits get re-planned
    around, measured).  State = (carry, nfall, pain, emitted-rows); rank by
    (nfall, near-top leftover pain).  Returns (sol, msg)."""
    states = [(None, None, 0, 0, [])]    # K, dK, nfall, pain, rows
    for j in range(R2):
        d = dt[j]
        nxt = []
        for K, dK, nfall, pain, rows in states:
            X = list(E[j])
            if K is not None:
                if dK > d:
                    continue
                KL = lift(K, dK, d) if dK < d else list(K)
                for k, v in enumerate(KL):
                    X[k] += v
            if abs(X[d]) != 1:
                continue
            opts = [None] if j + 1 >= R2 else [1, -1]
            for ov in opts:
                presc = ({} if ov is None else
                         _presc_with_override(E, dt, R2, j, W, L, ov))
                cells = [0] * (d + 1)
                cells[d] = X[d]
                P = [0] * d
                nf2, pain2 = nfall, pain
                for k in range(d - 1, -1, -1):
                    cap = comb(d, k)
                    if k in presc:
                        e = X[k] - presc[k]
                        if abs(e) <= cap and (e - cap) % 2 == 0:
                            cells[k] = e
                            P[k] = presc[k]
                            continue
                        nf2 += 1
                    v = pclamp(X[k], cap)
                    cells[k] = v
                    P[k] = X[k] - v
                for m in range(min(14, d)):
                    c = abs(P[d - 1 - m])
                    cc = comb(d - 1, d - 1 - m)
                    if c > cc:
                        pain2 += (c // max(1, cc)).bit_length()
                nxt.append((P, d - 1, nf2, pain2, rows + [cells]))
        if not nxt:
            return None, f'beam died at row {j}'
        nxt.sort(key=lambda s: (s[2], s[3]))
        seen, keep = set(), []
        for st in nxt:
            key = tuple(st[0][-30:])
            if key in seen:
                continue
            seen.add(key)
            keep.append(st)
            if len(keep) >= beam:
                break
        states = keep
        if verbose and (j % 16 == 0 or j > R2 - 4):
            print(f'    beamsweep row {j:3d}: states={len(states)} '
                  f'best nfall={states[0][2]} pain={states[0][3]}', flush=True)
    for K, dK, nfall, pain, rows in states:
        if not any(K):
            return rows, f'CLOSED (nfall={nfall}, pain={pain})'
    best = states[0]
    return None, (f'final carry nonzero in all {len(states)} states; best '
                  f'nfall={best[2]} pain={best[3]} '
                  f'carryL1~1e{len(str(sum(abs(x) for x in best[0])))-1}')

def repair_sweep(R2, dt, E, W=10, L=8, max_repairs=400, verbose=True):
    """THE REBALANCER RUNNING DURING THE SWEEP (surgical form).  Iterate:
    run the deterministic Steer sweep; on the FIRST prescription fallback
    (j,k) -- overflow e beyond the box by a small even amount -- apply one
    exact split-jump gather that shifts E[j][k] by -sign(e)*(overflow),
    donating to an earlier row's deeper cell (i0, k-r) with hard-cap room
    checks on donor and splash cells; re-run.  Identity preserved by
    construction at every step; the final witness is verified exactly."""
    E = [list(r) for r in E]
    history = []
    last_first = None
    stall = 0
    for it in range(max_repairs):
        dbg = []
        sol, msg, nfall, mc, nt = sweep_E(R2, dt, E, W=W, L=L, teleport=False,
                                          debug_fallback=dbg)
        if sol is not None:
            return sol, E, history, f'CLOSED after {it} repairs'
        if not dbg:
            return None, E, history, f'no fallback to repair; {msg}'
        j, k, depth, kind, e, cap = dbg[0]
        if kind == 'mag':
            over = abs(e) - cap
            shift = -(1 if e > 0 else -1) * (over + (over % 2))
        else:
            shift = 2 if e < 0 else -2   # parity: nudge by odd? e-parity is
            # fixed by X-parity; a 2-shift changes magnitude side.
        done = False
        for i0 in range(j - 1, max(-1, j - 9), -1):
            r = (j + dt[j]) - (i0 + dt[i0])
            k0 = k - r
            if k0 < 0 or not move_ok(dt, i0, k0, j):
                continue
            cap0 = comb(dt[i0], k0)
            if abs(E[i0][k0] - shift) > cap0:
                continue
            ok = True
            for t in range(1, r + 1):
                kk = k0 + r - t
                if abs(E[j][kk] + shift * comb(r, t)) > comb(dt[j], kk) * 2:
                    ok = False
                    break
            if not ok:
                continue
            apply_move(E, dt, i0, k0, j, -shift)
            history.append((it, j, k, kind, e, cap, i0, k0, shift))
            done = True
            break
        if not done:
            return None, E, history, (f'no donor for repair at row {j} cell '
                                      f'{k} ({kind}, e={e}, cap={cap}); {msg}')
        if verbose and (it < 12 or it % 25 == 0):
            print(f'    repair {it}: {kind} at ({j},{k}) depth {depth} '
                  f'e={e} cap={cap} shift={shift} donor=({i0},{k0}); '
                  f'sweep had died: {msg}', flush=True)
        if dbg[0][:2] == last_first:
            stall += 1
            if stall > 6:
                return None, E, history, (f'stalled repairing ({j},{k}); {msg}')
        else:
            stall = 0
        last_first = dbg[0][:2]
    return None, E, history, f'repair budget exhausted; {msg}'

def carve_pow(R0, dsrc, blocks, s):
    """s-fold carve (from the doubling module, verbatim math): source rows are
    2^s-fold products, target rows the (2^s-1)-th backward difference."""
    R = (2 ** s) * R0
    dt = prof(R, 0)
    rows, degs = [list(b) for b in blocks], list(dsrc)
    n = R0
    for lev in range(s):
        nrows, ndegs = [], []
        for a in range(2 * n - 1):
            lo, hi = max(0, a - n + 1), min(a, n - 1)
            dcap = max(degs[i] + degs[a - i] for i in range(lo, hi + 1))
            cells = [0] * (dcap + 1)
            for i in range(lo, hi + 1):
                pl = lift(conv(rows[i], rows[a - i]), degs[i] + degs[a - i], dcap)
                for k, v in enumerate(pl):
                    cells[k] += v
            nrows.append(cells)
            ndegs.append(dcap)
        rows, degs, n = nrows, ndegs, 2 * n - 1
    T = 2 ** s - 1
    E = []
    for j in range(R):
        dj = dt[j]
        cells = [0] * (dj + 1)
        for t in range(T + 1):
            a = j - t
            if a < 0 or a >= len(rows):
                continue
            c = comb(T, t) * (-1) ** t
            pl = lift(rows[a], degs[a], dj) if degs[a] < dj else rows[a]
            for k, v in enumerate(pl):
                cells[k] += c * v
        E.append(cells)
    return dt, E

# ---------------------------------------------------------------------------
# 5. drivers
# ---------------------------------------------------------------------------
def load_witness(path):
    w = json.load(open(path))
    return w['R'], w['profile'], w['blocks'], w

def save_witness(path, R, d, B, note):
    a = all(admissible(B[i], d[i]) for i in range(R))
    b = rows_identity(B, d, R)
    assert a and b, f'REFUSING to save unverified witness: adm={a} id={b}'
    json.dump({'R': R, 'profile': d, 'blocks': B, 'verified': True,
               'source': note}, open(path, 'w'))
    print(f'  SAVED verified witness -> {path}')
    return True

def cmd_metrics(args):
    for path in args.witness:
        R, d, B, w = load_witness(path)
        a = all(admissible(B[i], d[i]) for i in range(R))
        b = rows_identity(B, d, R)
        print(f'{os.path.basename(path)}: R={R} admissible={a} identity={b}')
        print_metrics('metrics', R, d, B)
        if args.carve and R <= 64:
            dt = prof(2 * R, 0)
            E = carve(R, d, B, dt)
            idok = rows_identity(E, dt, 2 * R)
            ov, at = rows_overshoot(dt, E)
            print(f'  carve {R}->{2*R}: identity={idok} overshoot={float(ov):.6g} at {at}')

def cmd_slim(args):
    R, d, B, w = load_witness(args.witness[0])
    print(f'slimming {os.path.basename(args.witness[0])} '
          f'(budget {args.budget}s, topn {args.topn})')
    m0 = print_metrics('before', R, d, B)
    B2, applied, log = slim_greedy(R, d, B, budget_s=args.budget,
                                   topn=args.topn)
    m1 = print_metrics('after ', R, d, B2)
    a = all(admissible(B2[i], d[i]) for i in range(R))
    b = rows_identity(B2, d, R)
    par = all((v - comb(d[i], k)) % 2 == 0 for i in range(R)
              for k, v in enumerate(B2[i]))
    print(f'  moves applied = {applied}; exact recheck: admissible={a} '
          f'identity={b} parity={par}')
    if a and b and args.out:
        save_witness(args.out, R, d, B2,
                     f'split-jump greedy slimming of {os.path.basename(args.witness[0])}'
                     f' ({applied} moves); rebalancer=amm12592_slimming_rebalancer_boxeph.py')
    dt = prof(2 * R, 0)
    E = carve(R, d, B2, dt)
    ov, at = rows_overshoot(dt, E)
    print(f'  carve {R}->{2*R} overshoot after slimming: {float(ov):.6g} at {at}'
          f'  (before: recompute with --metrics --carve)')

def _sweep_and_verify(R2, dt, E, W, L, teleport, tag, outpath=None):
    t0 = time.time()
    sol, msg, nfall, mc, ntele = sweep_E(R2, dt, E, W=W, L=L, teleport=teleport)
    el = time.time() - t0
    if sol is None:
        print(f'  {tag}: {msg}  (nfall={nfall} maxcarry~1e{len(str(mc))-1} '
              f'ntele={ntele})  [{el:.1f}s]', flush=True)
        return False
    a = all(admissible(sol[i], dt[i]) for i in range(R2))
    b = rows_identity(sol, dt, R2)
    print(f'  {tag}: CLOSED nfall={nfall} ntele={ntele}; verify adm={a} id={b} '
          f'[{el:.1f}s]', flush=True)
    if a and b and outpath:
        save_witness(outpath, R2, dt, sol, tag)
        return True
    return a and b

def cmd_sweep(args):
    R, d, B, w = load_witness(args.witness[0])
    R2 = 2 * R
    dt = prof(R2, 0)
    print(f'carve {R}->{R2} from {os.path.basename(args.witness[0])}')
    t0 = time.time()
    E = carve(R, d, B, dt)
    idok = rows_identity(E, dt, R2)
    ov, at = rows_overshoot(dt, E)
    print(f'  carve identity={idok} overshoot={float(ov):.6g} at {at} '
          f'[{time.time()-t0:.1f}s]')
    assert idok
    if args.bandvac >= 0:
        t0 = time.time()
        kf = Fraction(args.keep[0], args.keep[1])
        E, stuck, sm = band_vacuum(R2, dt, E, M0=args.bandvac,
                                   spread=args.spread, keep=kf, verbose=True)
        idok = rows_identity(E, dt, R2)
        ov, at = rows_overshoot(dt, E)
        print(f'  after bandvac(M0={args.bandvac}, spread={args.spread}, '
              f'keep={kf}): identity={idok} overshoot={float(ov):.6g} at {at} '
              f'stuck={stuck}({sm}) [{time.time()-t0:.1f}s]')
        assert idok
    for W in args.W:
        for L in args.L:
            for tp in ([False, True] if args.teleport == 'both'
                       else [args.teleport == 'on']):
                tag = f'W={W} L={L} teleport={tp}'
                if _sweep_and_verify(R2, dt, E, W, L, tp, tag, args.out):
                    return

def cmd_sweep32x4(args):
    W_ = {w['R']: w for w in json.load(open(os.path.join(
        HERE, 'amm12592_floor_witnesses_R8_R16_R32.json')))}
    w32 = W_[32]
    print('double-fold carve 32 x 2^2 -> 128 (T=3 differencing)')
    t0 = time.time()
    dt, E = carve_pow(32, w32['profile'], w32['blocks'], 2)
    idok = rows_identity(E, dt, 128)
    ov, at = rows_overshoot(dt, E)
    print(f'  identity={idok} overshoot={float(ov):.6g} at {at} '
          f'[{time.time()-t0:.1f}s]')
    assert idok
    if args.bandvac >= 0:
        t0 = time.time()
        kf = Fraction(args.keep[0], args.keep[1])
        E, stuck, sm = band_vacuum(128, dt, E, M0=args.bandvac,
                                   spread=args.spread, keep=kf, verbose=True)
        idok = rows_identity(E, dt, 128)
        ov, at = rows_overshoot(dt, E)
        print(f'  after bandvac(M0={args.bandvac}, spread={args.spread}, '
              f'keep={kf}): identity={idok} overshoot={float(ov):.6g} at {at} '
              f'stuck={stuck}({sm}) [{time.time()-t0:.1f}s]')
        assert idok
    for W in args.W:
        for L in args.L:
            for tp in ([False, True] if args.teleport == 'both'
                       else [args.teleport == 'on']):
                tag = f'32x4 W={W} L={L} teleport={tp}'
                if _sweep_and_verify(128, dt, E, W, L, tp, tag, args.out):
                    return

def bern_cells(sig):
    """exact Bernstein decode of sig at its own degree (no box test)."""
    m = len(sig) - 1
    res = list(sig)
    de = [0] * (m + 1)
    for k in range(m, -1, -1):
        v = res[m - k]
        de[k] = v
        if v:
            t = S.qk(k)
            off = m - k
            for s2 in range(k + 1):
                res[off + s2] -= v * t[s2]
    return de

S.RANKS['cellmax'] = lambda sg: (len(sg), max((abs(c) for c in bern_cells(sg)),
                                              default=0))
S.RANKS['cellsum'] = lambda sg: (len(sg), sum(abs(c) for c in bern_cells(sg)))

def cmd_hunt_slim(args):
    """beam hunt at R with rank = Bernstein-cell roughness of the residual:
    slim witnesses are trajectories whose residuals stay Bernstein-smooth."""
    R = args.R
    tgt = prof(R, 0)
    print(f'hunt-slim R={R} rank={args.rank} beam={args.beam}')
    t0 = time.time()
    sol, msg = S.solve_prof(tgt, beam=args.beam, ctrl=2, span=2, seed=None,
                            rand_frac=0.0, dedup=999, end_ctrl=3, end_span=6,
                            rank=args.rank, verbose=True)
    print(f'  -> {msg} [{time.time()-t0:.1f}s]')
    if sol is None:
        return
    a, b = S.verify_witness(R, sol, tgt)
    print(f'  verify: admissible={a} identity={b}')
    if a and b and args.out:
        save_witness(args.out, R, tgt, sol,
                     f'direct beam rank={args.rank} beam={args.beam} (slim hunt)')
        m = print_metrics('hunted', R, tgt, sol)
        if R <= 64:
            dt = prof(2 * R, 0)
            E = carve(R, tgt, sol, dt)
            ov, at = rows_overshoot(dt, E)
            print(f'  carve {R}->{2*R} overshoot: {float(ov):.6g} at {at}')

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('selftest')
    p = sub.add_parser('metrics')
    p.add_argument('witness', nargs='+')
    p.add_argument('--carve', action='store_true')
    p = sub.add_parser('slim')
    p.add_argument('witness', nargs=1)
    p.add_argument('--budget', type=int, default=600)
    p.add_argument('--topn', type=int, default=8)
    p.add_argument('--out', default=None)
    p = sub.add_parser('sweep')
    p.add_argument('witness', nargs=1)
    p.add_argument('--W', type=int, nargs='+', default=[12, 25, 40])
    p.add_argument('--L', type=int, nargs='+', default=[12, 24])
    p.add_argument('--teleport', choices=['on', 'off', 'both'], default='both')
    p.add_argument('--bandvac', type=int, default=-1)
    p.add_argument('--spread', type=int, default=3)
    p.add_argument('--keep', type=int, nargs=2, default=[1, 1])
    p.add_argument('--out', default=None)
    p = sub.add_parser('sweep32x4')
    p.add_argument('--W', type=int, nargs='+', default=[12, 25, 40])
    p.add_argument('--L', type=int, nargs='+', default=[12, 24])
    p.add_argument('--teleport', choices=['on', 'off', 'both'], default='both')
    p.add_argument('--bandvac', type=int, default=-1)
    p.add_argument('--spread', type=int, default=3)
    p.add_argument('--keep', type=int, nargs=2, default=[1, 1])
    p.add_argument('--out', default=None)
    p = sub.add_parser('hunt-slim')
    p.add_argument('--R', type=int, default=64)
    p.add_argument('--beam', type=int, default=150)
    p.add_argument('--rank', default='cellmax', choices=sorted(S.RANKS))
    p.add_argument('--out', default=None)
    args = ap.parse_args()
    if args.cmd == 'selftest':
        sys.exit(0 if selftest() else 1)
    {'metrics': cmd_metrics, 'slim': cmd_slim, 'sweep': cmd_sweep,
     'sweep32x4': cmd_sweep32x4, 'hunt-slim': cmd_hunt_slim}[args.cmd](args)

if __name__ == '__main__':
    main()
