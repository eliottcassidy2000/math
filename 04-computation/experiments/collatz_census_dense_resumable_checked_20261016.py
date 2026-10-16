#!/usr/bin/env python3
"""Exact Collatz parity census on a dense log-uniform grid, resumable across octave boundaries -- CHECKED version.

Identical arithmetic to collatz_census_dense_resumable_20261016.py, plus redundancy that the S38 audit (MISTAKE-623, report
E6) found missing: every accumulator that feeds a published quantity now has a checksum partner computed in the same
loop -- the reverse-weighted sub-block moment sum (n_hi - n) t next to the forward one (their sum must be (2^17 - 1) Sigma t),
and Sigma (hp + 1), Sigma (hp + 1)^2 next to Sigma hp, Sigma hp^2 -- and failing blocks are recomputed.  After the grid pass,
the grid values must satisfy the Lipschitz bounds |Delta S| <= Delta n, |Delta H| <= 2 Delta n and the prefix sum of H^2 must
be non-decreasing; --check-SH now recomputes S, H from the resume state's class counts and aborts on a mismatch.

Same arithmetic as collatz_census_extension_2p40_20261015.py (whose numba kernels tval and grid_pass are reused):
parity bits below 2^LOW from a table (collatz_census_extension_2p37_20261015.py), on-the-fly descent with the exact
8-step jump T^8(256h + l) = 3^(a(l)) h + T^8(l) above, single uint64 steps above 2^64/26 and exact two-limb arithmetic
beyond the uint64 guard.  t(n) = (-1)^tau(n), S(N) = sum_(n<=N) t(n), H(N) = sum_(n<=N) (3[n = 2 mod 3] - 1) t(n).

New here:
  * PER grid points per octave, N_j = round(2^(10 + j/PER)) (deduplicated), restricted to [2^K0, 2^LOGMAX); PER may be
    non-integer (an irrational density avoids the exact aliasing of the Mercator family, THM-4659 statement 5);
  * resume: with --resume STATE.json the heavy pass covers only the blocks of [2^K0, 2^LOGMAX), starting from the
    exact state at 2^K0 written by a previous run (class counts modulo 243 of all n < 2^K0, H(2^K0 - 1), and the exact
    sums of H and H^2 over n < 2^K0);
  * the final state at 2^LOGMAX is written to OUT_state.json (exact integers) for the next resume;
  * outputs: OUT_full.npz (class sums modulo 243 at every grid point; large, kept out of the repository),
    OUT.npz (grid, S, class sums modulo 3, H, exact prefix sums of H and H^2 as strings and floats: committed),
    OUT_moments.npz (signed moments of every 2^17 sub-block of the processed range), OUT.json (summary, exact octave
    sums of H and H^2 for the processed octaves >= 20).
Checks: the jump table against single steps; 20000 random n in [max(2^LOW, 2^K0), 2^LOGMAX) against Python
trajectories; on a fresh run every grid point shared with the 64-per-octave 2^40 census must agree exactly
(--check64 census40.npz); on resume the state at 2^K0 must reproduce S(2^K0) and H(2^K0) if given (--check-SH S H).
Session opus-2026-10-16-S38.  Usage:
  python collatz_census_dense_resumable_checked_20261016.py LOGMAX LOW PER OUTPREFIX [--resume STATE.json]
         [--check64 census40.npz] [--check-SH S H]
"""
import sys, os, json, time, random, argparse
import importlib.util
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('c40', os.path.join(HERE, 'collatz_census_extension_2p40_20261015.py'))
c40 = importlib.util.module_from_spec(spec)
_argv = sys.argv
sys.argv = ['x', '30', '26']
spec.loader.exec_module(c40)
sys.argv = _argv

ap = argparse.ArgumentParser()
ap.add_argument('logmax', type=int)
ap.add_argument('low', type=int)
ap.add_argument('per', type=float)
ap.add_argument('out')
ap.add_argument('--resume', default=None)
ap.add_argument('--check64', default=None)
ap.add_argument('--check-SH', nargs=2, type=int, default=None)
args = ap.parse_args()
M = 243
BLOCK = 1 << 20
SUB = 1 << 17


@nb.njit(parallel=True)
def heavy_range(P, b0, b1, B, Sb, Mod, low, A, Cc, jguard, guard):
    nb_ = b1 - b0
    cnt = np.zeros((nb_, Mod), dtype=np.int64)
    nsub = B // Sb
    mom = np.zeros((nb_ * nsub, 2), dtype=np.int64)
    momr = np.zeros(nb_ * nsub, dtype=np.int64)
    hs1 = np.zeros(nb_, dtype=np.int64)
    hs2 = np.zeros(nb_, dtype=np.int64)
    hs1b = np.zeros(nb_, dtype=np.int64)
    hs2b = np.zeros(nb_, dtype=np.int64)
    hend = np.zeros(nb_, dtype=np.int64)
    mx = np.zeros(nb_, dtype=np.int64)
    bad = np.zeros(nb_, dtype=np.int64)
    tott = np.zeros(nb_, dtype=np.int64)
    for i in nb.prange(nb_):
        b = b0 + i
        lo = b * B
        n0 = lo if lo > 0 else 1
        r = n0 % Mod
        r3 = n0 % 3
        hp = 0
        s1 = 0
        s2 = 0
        s1b = 0
        s2b = 0
        m = 0
        tt = 0
        for n in range(n0, lo + B):
            t = c40.tval(P, n, low, A, Cc, jguard, guard)
            if t == 0:
                bad[i] = n
                t = 1
            tt += t
            cnt[i, r] += t
            k = (n - lo) >> 17
            mom[i * nsub + k, 0] += t
            mom[i * nsub + k, 1] += (n - lo - (k << 17)) * t
            momr[i * nsub + k] += ((Sb - 1) - (n - lo - (k << 17))) * t
            if r3 == 2:
                hp += 2 * t
            else:
                hp -= t
            s1 += hp
            s2 += hp * hp
            s1b += hp + 1
            s2b += (hp + 1) * (hp + 1)
            ah = hp if hp >= 0 else -hp
            if ah > m:
                m = ah
            r += 1
            if r == Mod:
                r = 0
            r3 += 1
            if r3 == 3:
                r3 = 0
        hs1[i] = s1
        hs2[i] = s2
        hs1b[i] = s1b
        hs2b[i] = s2b
        hend[i] = hp
        mx[i] = m
        tott[i] = tt
    return cnt, mom, hs1, hs2, hend, mx, bad, tott, momr, hs1b, hs2b


@nb.njit(parallel=True)
def block_invariants(cnt, mom, hend, tott, nsub, momr, hs1, hs2, hs1b, hs2b, b0, B, Sb):
    # per block: sum of class counts = sum of t (register) = sum of sub-block moments; 3 C(2 mod 3) - S = H increment
    nb_ = cnt.shape[0]
    badb = np.zeros(nb_, dtype=np.int64)
    for i in nb.prange(nb_):
        s = 0
        c2 = 0
        for r in range(cnt.shape[1]):
            s += cnt[i, r]
            if r % 3 == 2:
                c2 += cnt[i, r]
        sm = 0
        for k in range(nsub):
            sm += mom[i * nsub + k, 0]
        if s != tott[i] or sm != tott[i] or 3 * c2 - s != hend[i]:
            badb[i] = 1
        for k in range(nsub):
            if mom[i * nsub + k, 1] + momr[i * nsub + k] != (Sb - 1) * mom[i * nsub + k, 0]:
                badb[i] = 1
        size = B if (b0 + i) > 0 else B - 1
        if hs1b[i] != hs1[i] + size or hs2b[i] != hs2[i] + 2 * hs1[i] + size:
            badb[i] = 1
    return badb


if __name__ == '__main__':
    t0 = time.time()
    LOGMAX, LOW, PER, OUT = args.logmax, args.low, args.per, args.out
    A, Cc = c40.jump_tables()
    rng = random.Random(16101610)
    for _ in range(2000):
        x = rng.randrange(1, 1 << 50)
        y = x
        for _ in range(8):
            y = y // 2 if y % 2 == 0 else (3 * y + 1) // 2
        assert y == int(A[x & 255]) * (x >> 8) + int(Cc[x & 255])
    if args.resume:
        st = json.load(open(args.resume))
        K0 = int(st['logmax'])
        start_counts = np.array([int(x) for x in st['class_counts_below']], dtype=object)
        h_start = int(st['H_last'])
        cum1_start = int(st['sumH_below'])
        cum2_start = int(st['sumH2_below'])
        if args.check_SH:
            cc = [int(x) for x in st['class_counts_below']]
            t0n = -1 if K0 % 2 else 1
            cc[(1 << K0) % M] += t0n
            S0 = sum(cc)
            H0r = 3 * sum(cc[2::3]) - S0
            ok = (S0 == args.check_SH[0] and H0r == args.check_SH[1] and int(st['S_at_2k']) == S0
                  and int(st['H_at_2k']) == H0r and int(st['H_last']) == H0r - (2 * t0n if ((1 << K0) % 3) == 2 else -t0n))
            print(f'resume state recomputed from its class counts: S = {S0}, H = {H0r}; matches the given values: {ok}',
                  flush=True)
            if not ok:
                raise SystemExit('resume state does not reproduce the given S, H')
    else:
        K0 = 0
        start_counts = np.zeros(M, dtype=object)
        h_start = cum1_start = cum2_start = 0
    P, patched = c40.ext.parity_bits(LOW)
    low = 1 << LOW
    print(f'parity table to 2^{LOW}: {time.time() - t0:.1f} s; {len(patched)} words recomputed exactly', flush=True)
    bad = 0
    lo_chk = max(low, 1 << K0)
    for _ in range(20000):
        n = rng.randrange(lo_chk, 1 << LOGMAX)
        if c40.tval(P, n, low, A, Cc, c40.JGUARD, c40.GUARD) != 1 - 2 * c40.ext.tau_parity(n):
            bad += 1
    print(f'spot check of 20000 random n in [2^{max(LOW, K0)}, 2^{LOGMAX}): {bad} mismatches', flush=True)
    if bad:
        raise SystemExit('parity mismatch')
    b0 = (1 << K0) // BLOCK if K0 >= 20 else 0
    b1 = (1 << LOGMAX) // BLOCK
    if K0 and (1 << K0) % BLOCK:
        raise SystemExit('resume point must be a multiple of the block size')
    t1 = time.time()
    cnt, mom, hs1, hs2, hend, mxb, badb, tott, momr, hs1b, hs2b = heavy_range(P, b0, b1, BLOCK, SUB, M, low, A, Cc, c40.JGUARD, c40.GUARD)
    if badb.any():
        raise SystemExit('128-bit overflow')
    # invariant check; recompute any block that fails it (and report), up to 3 attempts
    for attempt in range(3):
        inv = block_invariants(cnt, mom, hend, tott, BLOCK // SUB, momr, hs1, hs2, hs1b, hs2b, b0, BLOCK, SUB)
        badi = np.nonzero(inv)[0]
        print(f'block invariant check (attempt {attempt}): {len(badi)} inconsistent blocks', flush=True)
        if len(badi) == 0:
            break
        for i in badi:
            b = b0 + int(i)
            c_, m_, h1_, h2_, he_, mx_, bd_, tt_, mr_, h1b_, h2b_ = heavy_range(P, b, b + 1, BLOCK, SUB, M, low, A, Cc, c40.JGUARD, c40.GUARD)
            print(f'   block {b}: recomputed (old tott {int(tott[i])}, new {int(tt_[0])}; old hend {int(hend[i])}, new {int(he_[0])})',
                  flush=True)
            cnt[i] = c_[0]
            mom[int(i) * (BLOCK // SUB):(int(i) + 1) * (BLOCK // SUB)] = m_
            hs1[i], hs2[i], hend[i], mxb[i], tott[i] = h1_[0], h2_[0], he_[0], mx_[0], tt_[0]
            momr[int(i) * (BLOCK // SUB):(int(i) + 1) * (BLOCK // SUB)] = mr_
            hs1b[i], hs2b[i] = h1b_[0], h2b_[0]
    else:
        raise SystemExit('block invariants still failing after recomputation')
    mx = int(mxb.max())
    if BLOCK * mx * mx >= 2 ** 63:
        raise SystemExit('possible int64 overflow in block sums of hp^2')
    print(f'heavy pass over blocks [{b0}, {b1}): {time.time() - t1:.1f} s; max within-block |hp| = {mx}', flush=True)
    nb_ = b1 - b0
    tot_py = sum(int(x) for x in tott)
    start = np.zeros((nb_, M), dtype=np.int64)
    base = np.array([int(x) for x in start_counts], dtype=np.int64)
    start[0] = base
    if nb_ > 1:
        start[1:] = base + np.cumsum(cnt[:-1], axis=0)
    H0 = [0] * nb_
    cum1 = [0] * (nb_ + 1)
    cum2 = [0] * (nb_ + 1)
    cum1[0], cum2[0] = cum1_start, cum2_start
    h = h_start
    for i in range(nb_):
        b = b0 + i
        H0[i] = h
        size = BLOCK if b > 0 else BLOCK - 1
        sH = size * h + int(hs1[i])
        sH2 = size * h * h + 2 * h * int(hs1[i]) + int(hs2[i])
        cum1[i + 1] = cum1[i] + sH
        cum2[i + 1] = cum2[i] + sH2
        h += int(hend[i])
    raw = sorted(set(int(round(2 ** (10 + j / PER))) for j in range(0, int((LOGMAX - 10) * PER) + 1)))
    grid = np.array([g for g in raw if (1 << K0) <= g < (1 << LOGMAX)], dtype=np.int64)
    bidx = grid // BLOCK
    gblocks = np.unique(bidx)
    glo = np.searchsorted(bidx, gblocks, side='left').astype(np.int64)
    ghi = np.searchsorted(bidx, gblocks, side='right').astype(np.int64)
    t2 = time.time()
    snapC, snapH, snapP1, snapP2 = c40.grid_pass(P, gblocks, BLOCK, M, start[gblocks - b0].copy(),
                                                 np.array([H0[int(b) - b0] for b in gblocks], dtype=np.int64), grid,
                                                 glo, ghi, low, A, Cc, c40.JGUARD, c40.GUARD)
    print(f'grid pass over {len(gblocks)} blocks ({len(grid)} grid points): {time.time() - t2:.1f} s', flush=True)
    pre1, pre2 = [], []
    for i in range(len(grid)):
        b = int(bidx[i])
        hb = H0[b - b0]
        lo = b * BLOCK
        cntn = int(grid[i]) - max(lo, 1) + 1
        pre1.append(cum1[b - b0] + cntn * hb + int(snapP1[i]))
        pre2.append(cum2[b - b0] + cntn * hb * hb + 2 * hb * int(snapP1[i]) + int(snapP2[i]))
    # Lipschitz and monotonicity checks of the grid values (catch corrupted grid-pass output)
    Sg = snapC.sum(axis=1).astype(np.int64)
    gd = np.diff(grid)
    if len(grid) > 1:
        if np.any(np.abs(np.diff(Sg)) > gd) or np.any(np.abs(np.diff(snapH)) > 2 * gd):
            raise SystemExit('grid values violate the Lipschitz bounds')
        if any(pre2[i + 1] < pre2[i] for i in range(len(pre2) - 1)):
            raise SystemExit('prefix sum of H^2 decreases')
    print(f'grid Lipschitz and monotonicity checks passed on {len(grid)} points', flush=True)
    final_counts = [int(x) for x in (start[-1] + cnt[-1])]
    if sum(final_counts) != sum(int(x) for x in start_counts) + tot_py:
        raise SystemExit('final class counts disagree with the summed block totals')
    Nmax = 1 << LOGMAX
    tN = -1 if LOGMAX % 2 else 1
    S_last = sum(final_counts)
    H_last = 3 * sum(final_counts[2::3]) - S_last
    fc = list(final_counts)
    fc[Nmax % M] += tN
    S_N = sum(fc)
    H_N = 3 * sum(fc[2::3]) - S_N
    if H_last != h:
        print(f'WARNING: H from final class counts {H_last} differs from the block chain {h}', flush=True)
        raise SystemExit('state inconsistency')
    state = dict(logmax=LOGMAX, class_counts_below=[str(x) for x in final_counts], H_last=str(h),
                 sumH_below=str(cum1[-1]), sumH2_below=str(cum2[-1]), S_at_2k=S_N, H_at_2k=H_N)
    json.dump(state, open(OUT + '_state.json', 'w'), indent=1)
    S = snapC.sum(axis=1)
    C3 = np.stack([snapC[:, a::3].sum(axis=1) for a in range(3)], axis=1)
    msg = f'S(2^{LOGMAX}) = {S_N}, H(2^{LOGMAX}) = {H_N}'
    if args.resume:
        msg += f'; resumed from 2^{K0} with S = {st["S_at_2k"]}, H = {st["H_at_2k"]}'
    print(msg, flush=True)
    checks = {}
    if args.check_SH and args.resume:
        checks['resume_state_matches'] = True      # verified before the heavy pass
    if args.check64:
        with np.load(args.check64) as Z:
            g0 = np.array(Z['grid'])
            c0 = np.array(Z['classsums'])
            h0 = np.array(Z['H'])
            q2 = [int(x) for x in Z['preH2_str']]
        sel = (g0 >= (1 << K0)) & (g0 < Nmax)
        g0, c0, h0 = g0[sel], c0[sel], h0[sel]
        q2 = [q for q, s_ in zip(q2, sel) if s_]
        idx = np.searchsorted(grid, g0)
        okidx = (idx < len(grid))
        same = bool(okidx.all() and np.array_equal(grid[idx], g0) and np.array_equal(snapC[idx], c0)
                    and np.array_equal(snapH[idx], h0) and all(pre2[int(i)] == q for i, q in zip(idx, q2)))
        checks['agrees_with_64_per_octave_census'] = same
        checks['shared_points'] = int(len(g0))
        print(f'{len(g0)} shared grid points with the 64-per-octave census agree exactly: {same}', flush=True)
    octs = {}
    for k in range(max(20, K0), LOGMAX):
        i0, i1 = (1 << k) // BLOCK - b0, (1 << (k + 1)) // BLOCK - b0
        octs[str(k)] = (str(cum1[i1] - cum1[i0]), str(cum2[i1] - cum2[i0]))
    np.savez_compressed(OUT + '_full.npz', grid=grid, classsums=snapC, H=snapH)
    np.savez_compressed(OUT + '.npz', grid=grid, S=S, C3=C3, H=snapH,
                        preH1=np.array([float(x) for x in pre1]), preH2=np.array([float(x) for x in pre2]),
                        preH1_str=np.array([str(x) for x in pre1]), preH2_str=np.array([str(x) for x in pre2]),
                        per=PER, logmax=LOGMAX, k0=K0)
    np.savez_compressed(OUT + '_moments.npz', moments=mom, block0=np.int64(b0), N=np.int64(Nmax), B=np.int64(SUB))
    json.dump(dict(logmax=LOGMAX, k0=K0, low=LOW, per=PER, n_grid=int(len(grid)), S=S_N, H=H_N, t_N=tN, checks=checks,
                   octave_sums_H_H2=octs, seconds=time.time() - t0), open(OUT + '.json', 'w'), indent=1)
    print(f'wrote {OUT}.npz, {OUT}_full.npz, {OUT}_moments.npz, {OUT}.json, {OUT}_state.json; total {time.time() - t0:.1f} s',
          flush=True)
