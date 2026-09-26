#!/usr/bin/env python3
"""procgen_gates_20260925_stats.py -- PART B of the cycle-gate equidistribution lane.

For every clock (p, a) with p <= 30 (q = 3) and p <= 24 (q = 5) the full multiset
{c_w mod M : w in W(p,a)}, M = |2^p - q^a|, is generated (meet-in-the-middle, streamed) and:
  N (zeros), collisions Coll = sum_r n_r^2 against the uniform model C + C(C-1)/M, star
  discrepancy D* of {c_w / M mod 1}, near-integer counts lambda(delta), the fraction of words with
  |x_w| < 1, and the perigee statistics of the rational cycles (eligibility: least point >= 1);
  for M <= 2^23 the whole spectrum S(h), h mod M (FFT of the residue histogram): max_{h != 0}|S(h)|,
  the arithmetic form of the top frequencies (h = u 2^-j q^k mod M with |u| minimal), the share of
  the Parseval energy on the structured set Sigma_3 = {+-u 2^-j q^k : 1 <= u <= 3, j <= p, k <= a};
  for M > 2^23, S(h) by the DP on Sigma_3 and on 1000 random h.
Every internal consistency check raises on failure.  stdout = results, stderr = timing/memory.
"""
import math
import os
import random
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from procgen_gates_20260925_core import (S_abs_all_dp, S_dp, cmin_cmax, gate, lean_malloc,  # noqa: E402
                                         lyndon_count, report_mem, residue_chunks)

FFT_MAX = 1 << 20      # numpy rfft needs ~150 B/point for non-smooth lengths: keep below 700 MB
FULLDP_MAX = 1 << 23   # full spectrum by the blocked DP for near-critical clocks (amp >= 1.5)
BUCKET_CAP = 12_000_000
CHUNK = 1 << 18
NCO = 1 << 16          # coarse residue bins used to cut the sort into buckets
DELTAS = (0.1, 0.03, 0.01, 0.003)
U_STRUCT = 3
RNG = random.Random(20260925)


def coarse_bin(r, M):
    """bin index of residues r in [0, M) among NCO equal bins (float; the same function is used in
    both passes, so only consistency matters)."""
    b = np.floor(r.astype(np.float64) * (NCO / M)).astype(np.int64)
    np.clip(b, 0, NCO - 1, out=b)
    return b


def runs_stats(arr, block=1 << 22):
    """sum of squared run lengths and max run length of a sorted 1-d array, block-wise."""
    n = len(arr)
    coll = 0
    mx = 0
    cur_val = None
    cur_len = 0
    for b0 in range(0, n, block):
        blk = arr[b0:b0 + block]
        # boundaries inside the block
        neq = np.flatnonzero(blk[1:] != blk[:-1]) + 1
        starts = np.concatenate(([0], neq))
        ends = np.concatenate((neq, [len(blk)]))
        lens = ends - starts
        # merge first run with the carried run
        if cur_val is not None and blk[0] == cur_val:
            lens[0] += cur_len
        elif cur_val is not None:
            coll += cur_len * cur_len
            mx = max(mx, cur_len)
        # all runs but the last are complete
        if len(lens) > 1:
            done = lens[:-1].astype(np.int64)
            coll += int((done * done).sum())
            mx = max(mx, int(done.max()))
        cur_val = blk[-1]
        cur_len = int(lens[-1])
    if cur_val is not None:
        coll += cur_len * cur_len
        mx = max(mx, cur_len)
    return coll, mx


def disc_block(arr, M, offset, n):
    """contribution to the star discrepancy of the sorted block arr (global ranks offset..)."""
    i = np.arange(offset, offset + len(arr), dtype=np.float64)
    u = arr.astype(np.float64) / M
    return float(max(((i + 1) / n - u).max(), (u - i / n).max()))


def structured_set(p, a, M, q, U=U_STRUCT):
    """Sigma_U folded to representatives in [1, M//2]."""
    inv2 = pow(2, -1, M)
    S = set()
    for j in range(p + 1):
        b = pow(inv2, j, M)
        for k in range(a + 1):
            g = b * pow(q, k, M) % M
            for u in range(1, U + 1):
                h = u * g % M
                h = min(h, M - h)
                if h:
                    S.add(h)
    return np.array(sorted(S), dtype=np.int64)


def decompose(h, p, a, M, q):
    """(u, j, k) with h = u 2^-j q^k mod M and |u| minimal (u centred), over 0<=j<=p, 0<=k<=a."""
    best = None
    qinv = pow(q, -1, M)
    p2 = [pow(2, j, M) for j in range(p + 1)]
    qi = [pow(qinv, k, M) for k in range(a + 1)]
    for j in range(p + 1):
        hj = h * p2[j] % M
        for k in range(a + 1):
            u = hj * qi[k] % M
            if u > M // 2:
                u -= M
            if best is None or abs(u) < abs(best[0]):
                best = (u, j, k)
    return best


def elig_class(p, a, q):
    """classification of a clock by the proven perigee bounds (least point x_min of a rational cycle):
    dyadic G > 0:  x_min <= 1/(2^(p/a) - q)  and  x_min >= q^(a-1)/G;
    q-adic G < 0:  y_min <= 1/(q - 2^(p/a)),  y_min >= q^(a-1)/|G|  and  y_min >= c_min/|G|.
    returns 'none' (no rational cycle has least point >= 1), 'all', or 'transition'."""
    G = gate(p, a, q)
    M = abs(G)
    if G > 0:
        if 2 ** p > (q + 1) ** a:
            return "none"
        if q ** (a - 1) >= G:
            return "all"
        return "transition"
    if 2 ** p < (q - 1) ** a:
        return "none"
    if q ** (a - 1) >= M or cmin_cmax(p, a, q)[0] >= M:
        return "all"
    return "transition"


def clock_stats(p, a, q, do_perigee, do_fft):
    G = gate(p, a, q)
    M = abs(G)
    C = math.comb(p, a)
    exact_ok = cmin_cmax(p, a, q)[1] < (1 << 60)
    single = C <= BUCKET_CAP
    dt = np.uint32 if M < (1 << 32) else np.int64
    store = np.empty(C, dtype=dt) if single else None
    hist = np.zeros(M, dtype=np.int64) if (do_fft and M <= FFT_MAX) else None
    N = 0
    below1 = 0
    near = np.zeros(len(DELTAS), dtype=np.int64)
    elig = elig_prim = prim = 0
    k0 = 0
    coarse = np.zeros(NCO, dtype=np.int64)
    for r, c in residue_chunks(p, a, q, max_chunk=CHUNK, want_float=not exact_ok, want_exact=exact_ok):
        n = len(r)
        x = c.astype(np.float64) / M
        dist = np.minimum(r, M - r).astype(np.float64) / M
        big = x >= 0.5
        for di, d in enumerate(DELTAS):
            near[di] += int(np.count_nonzero(big & (dist <= d)))
        below1 += int(np.count_nonzero(x < 1.0))
        N += int(np.count_nonzero(r == 0))
        coarse += np.bincount(coarse_bin(r, M), minlength=NCO)
        if hist is not None:
            hist += np.bincount(r, minlength=M)
        if single:
            store[k0:k0 + n] = r
        k0 += n
        if do_perigee:
            assert exact_ok
            cc = c.copy()
            cmn = c.copy()
            first = np.zeros(n, dtype=np.int16)
            for s in range(1, p + 1):
                odd = (cc & 1).astype(bool)
                cc = np.where(odd, (q * cc + G) >> 1, cc >> 1)
                np.minimum(cmn, cc, out=cmn)
                newly = (first == 0) & (cc == c)
                first[newly] = s
            assert np.array_equal(cc, c)
            e = cmn >= M
            pr = first == p
            elig += int(np.count_nonzero(e))
            elig_prim += int(np.count_nonzero(e & pr))
            prim += int(np.count_nonzero(pr))
        del x, dist, big
    assert k0 == C
    if do_perigee:
        assert prim == p * lyndon_count(p, a)
    # sorted statistics
    BLK = 1 << 22
    if single:
        store.sort()
        coll, mx = runs_stats(store)
        dstar = max(disc_block(store[b0:b0 + BLK], M, b0, C) for b0 in range(0, C, BLK))
        nzero = int(np.count_nonzero(store[:max(1, N)] == 0))
        assert nzero == N or N == 0
        del store
    else:
        # adaptive buckets from the coarse histogram: each bucket holds <= 8M residues when possible
        edges_idx = [0]
        acc = 0
        for i in range(NCO):
            if acc + coarse[i] > 4_000_000 and acc > 0:
                edges_idx.append(i)
                acc = 0
            acc += int(coarse[i])
        edges_idx.append(NCO)
        coll = mx = 0
        dstar = 0.0
        off = 0
        for b in range(len(edges_idx) - 1):
            i0, i1 = edges_idx[b], edges_idx[b + 1]
            want = int(coarse[i0:i1].sum())
            assert want <= 40_000_000
            arr = np.empty(want, dtype=dt)
            kk = 0
            for r, _ in residue_chunks(p, a, q, max_chunk=CHUNK, want_float=False):
                cb_ = coarse_bin(r, M)
                sel = r[(cb_ >= i0) & (cb_ < i1)]
                arr[kk:kk + len(sel)] = sel
                kk += len(sel)
            assert kk == want
            arr.sort()
            if len(arr):
                cb, mb = runs_stats(arr)
                coll += cb
                mx = max(mx, mb)
                dstar = max(dstar, max(disc_block(arr[b0:b0 + BLK], M, off + b0, C) for b0 in range(0, len(arr), BLK)))
                off += len(arr)
            del arr
        assert off == C
    res = dict(p=p, a=a, q=q, G=G, M=M, C=C, N=N, coll=coll, maxload=mx, dstar=dstar, below1=below1,
               near=near.tolist(), elig=elig if do_perigee else None, elig_prim=elig_prim if do_perigee else None)
    # spectrum: FFT of the histogram (M <= 2^20) or the blocked DP (near-critical, M <= 2^23)
    absS = None
    if hist is not None and C >= 200:
        F = np.fft.rfft(hist.astype(np.float64))
        del hist
        absS = np.abs(F)
        del F
        res["spec"] = "fft"
    elif do_fft and C >= 200 and M <= FULLDP_MAX and max(2 ** p, q ** a) / M >= 1.5:
        absS = S_abs_all_dp(p, a, q)
        res["spec"] = "dp-all"
    if absS is not None:
        # Parseval check: sum_{h} |S|^2 = M * coll
        tot = absS[0] ** 2 + 2.0 * float((absS[1:] ** 2).sum())
        assert abs(tot - M * coll) <= 1e-6 * M * coll + 1e-3
        e_nonzero = M * coll - C * C
        top = np.argsort(absS[1:])[::-1][:6] + 1
        res["maxS"] = float(absS[top[0]])
        res["top"] = [(int(h), float(absS[h]), decompose(int(h), p, a, M, q)) for h in top]
        Sig = structured_set(p, a, M, q)
        e_struct = 2.0 * float((absS[Sig] ** 2).sum())
        res["struct_frac"] = e_struct / e_nonzero if e_nonzero > 0 else float("nan")
        res["struct_size_frac"] = len(Sig) / (M // 2)
        mask = np.ones(len(absS), dtype=bool)
        mask[0] = False
        mask[Sig] = False
        res["max_unstruct"] = float(absS[mask].max()) if mask.any() else float("nan")
        multi = []
        for U in (1, 3, 10, 30):
            SU = structured_set(p, a, M, q, U)
            multi.append((U, 2.0 * float((absS[SU] ** 2).sum()) / e_nonzero if e_nonzero > 0 else float("nan"),
                          len(SU) / (M // 2)))
            if U == 30:
                m30 = np.ones(len(absS), dtype=bool)
                m30[0] = False
                m30[SU] = False
                res["max_unstruct30"] = float(absS[m30].max()) if m30.any() else float("nan")
        res["multi"] = multi
        res["S1"] = float(absS[1])
        # baseline: minimal |u| for random frequencies
        rnd = [decompose(RNG.randrange(1, M), p, a, M, q)[0] for _ in range(40)]
        res["rand_u_median"] = float(np.median(np.abs(rnd)))
    elif C >= 200:
        Sig = structured_set(p, a, M, q)
        vals = np.abs(S_dp(p, a, list(map(int, Sig)), q))
        hr = [RNG.randrange(1, M) for _ in range(1000)]
        vr = np.abs(S_dp(p, a, hr, q))
        res["maxS_struct"] = float(vals.max())
        res["h_struct"] = int(Sig[int(vals.argmax())])
        res["maxS_rand"] = float(vr.max())
        res["rms_rand"] = float(np.sqrt((vr ** 2).mean()))
        res["h_rand_dec"] = decompose(hr[int(vr.argmax())], p, a, M, q)
        res["h_struct_dec"] = decompose(int(Sig[int(vals.argmax())]), p, a, M, q)
        res["S1"] = float(abs(S_dp(p, a, [1], q)[0]))
    return res


def amp(p, a, q):
    return max(2 ** p, q ** a) / abs(gate(p, a, q))


def fmt_row(r):
    p, a, q, M, C, N = r["p"], r["a"], r["q"], r["M"], r["C"], r["N"]
    side = "+" if r["G"] > 0 else "-"
    exp_pairs = C * (C - 1) / M
    ex = r["coll"] - C
    lam = [r["near"][i] / (2 * DELTAS[i] * C) for i in range(len(DELTAS))]
    el = "" if r["elig"] is None else f"{r['elig'] / C:6.3f}"
    if "maxS" in r:
        spec = f"{r['maxS'] / C:7.4f} {r['max_unstruct'] / math.sqrt(C * math.log(M)):5.2f} {r['struct_frac']:6.3f}"
    elif "maxS_struct" in r:
        spec = f"{r['maxS_struct'] / C:7.4f} {r['maxS_rand'] / math.sqrt(C * math.log(M)):5.2f}{'dp':>7}"
    else:
        spec = f"{'':>7} {'':>5} {'':>6}"
    s1 = f"{r['S1'] / C:6.3f}" if "S1" in r else f"{'':>6}"
    return (f"{str((p, a)):>8}{side} {M:>11} {C:>10} {N:>3} {C / M:9.3g} {amp(p, a, q):8.3g} "
            f"{ex:>9} {exp_pairs:10.4g} {r['dstar']:7.4f} {r['below1'] / C:6.3f} "
            f"{lam[0]:5.2f} {lam[2]:5.2f} {lam[3]:5.2f} {el} {s1} {spec}")


HEADER = ("   (p,a)s           M          C   N       C/M      amp  Coll-C   C(C-1)/M      D*  x<1   "
          "lam.1 lam.01 lam.003 elig |S1|/C  maxS/C  mU/rt  Estr")


def run(q, pmax, perigee_all_upto):
    print("=" * 130)
    print(f"B1 (q = {q}): per-clock statistics, p <= {pmax}, clocks with C(p,a) >= 50 printed "
          f"(s = sign of G: + dyadic, - {q}-adic)")
    print("  amp = max(2^p, q^a)/|G| (distance to the critical line); Coll-C = excess coincident pairs;")
    print("  C(C-1)/M = its uniform-model mean; lam.d = #{|x|>=1/2, ||x|| <= d}/(2 d C); elig = share of words")
    print("  whose rational cycle has least point >= 1; maxS/C = max_{h!=0}|S(h)|/C (FFT) or max over Sigma_3")
    print("  (dp); mU/rt = max over h outside Sigma_3 (FFT) or over 1000 random h (dp), divided by sqrt(C ln M);")
    print("  Estr = share of the energy sum_{h!=0}|S(h)|^2 carried by Sigma_3.")
    print("=" * 130)
    print(HEADER)
    rows = []
    t0 = time.time()
    for p in range(1, pmax + 1):
        for a in range(1, p + 1):
            G = gate(p, a, q)
            M = abs(G)
            C = math.comb(p, a)
            if M >= (1 << 62) or M == 1:
                continue
            if p <= perigee_all_upto:
                dop = True
            else:  # only the transition window needs computing (bounds checked by perigee_theorems)
                dop = elig_class(p, a, q) == "transition"
            dop = dop and cmin_cmax(p, a, q)[1] < (1 << 60)
            r = clock_stats(p, a, q, do_perigee=dop, do_fft=True)
            rows.append(r)
            if C >= 50:
                print(fmt_row(r), flush=True)
        print(f"[q={q} p={p} done {time.time() - t0:.1f}s]", file=sys.stderr, flush=True)
    return rows


def perigee_theorems(rows, q):
    """check the proven classification elig_class on every clock where the perigees were computed."""
    cnt = {"none": 0, "all": 0, "transition": 0}
    for r in rows:
        if r["elig"] is None:
            continue
        cl = elig_class(r["p"], r["a"], q)
        if cl == "none":
            assert r["elig"] == 0, r
        elif cl == "all":
            assert r["elig"] == r["C"], r
        cnt[cl] += 1
    print(f"  perigee bounds checked (q = {q}) on every clock with computed perigees: 'none' (least point < 1 "
          f"forced) {cnt['none']} clocks, 'all' (least point >= 1 forced) {cnt['all']}, transition "
          f"{cnt['transition']} (computed).")


def summaries(rows, q):
    print()
    print(f"B2 (q = {q}): spectrum structure where the whole spectrum was computed (FFT for M <= 2^20; blocked DP for")
    print("  amp >= 1.5 and M <= 2^23), C >= 200")
    print("  top frequencies h (folded to [1, M/2]) written h = u 2^-j q^k mod M with |u| minimal; "
          "median |u| of 40 random h for scale")
    shown = 0
    for r in rows:
        if "top" not in r:
            continue
        if not (r["C"] >= 2000 and (amp(r["p"], r["a"], q) > 3 or r["p"] in (16, 20, 23))):
            continue
        shown += 1
        tops = "; ".join(f"h={h} |S|/C={s / r['C']:.3f} u={u:+d},j={j},k={k}" for h, s, (u, j, k) in r["top"][:4])
        print(f"  ({r['p']},{r['a']}) M={r['M']} C={r['C']} amp={amp(r['p'], r['a'], q):.3g} "
              f"rand|u|~{r['rand_u_median']:.0f}: {tops}")
    print()
    print(f"B2b (q = {q}): multi-scale structure.  Energy share E_U of sum_(h!=0)|S(h)|^2 on Sigma_U and the share")
    print("  of frequencies in Sigma_U (clocks with FFT, C >= 2000); mU30 = max outside Sigma_30 / sqrt(C ln M)")
    print(f"  {'(p,a)':>8} {'amp':>7} {'E_1':>6} {'E_3':>6} {'E_10':>6} {'E_30':>6} {'|S_30|/(M/2)':>12} {'mU30':>6}")
    for r in rows:
        if "multi" not in r or r["C"] < 2000:
            continue
        m = r["multi"]
        print(f"  {str((r['p'], r['a'])):>8} {amp(r['p'], r['a'], q):7.3g} {m[0][1]:6.3f} {m[1][1]:6.3f} "
              f"{m[2][1]:6.3f} {m[3][1]:6.3f} {m[3][2]:12.2e} "
              f"{r['max_unstruct30'] / math.sqrt(r['C'] * math.log(r['M'])):6.2f}")
    print()
    print(f"B2c (q = {q}): clocks without FFT (M > 2^23, C >= 200): Sigma_3 max and the random-h maximum, with the")
    print("  (u, j, k) form of both argmaxes; typical |u| for a random h is ~ M/(2(p+1)(a+1))")
    for r in rows:
        if "maxS_struct" not in r:
            continue
        u1, j1, k1 = r["h_struct_dec"]
        u2, j2, k2 = r["h_rand_dec"]
        print(f"  ({r['p']},{r['a']}) M={r['M']} C={r['C']} amp={amp(r['p'], r['a'], q):.3g}: "
              f"Sigma_3 max {r['maxS_struct'] / r['C']:.4f} C at u={u1:+d},j={j1},k={k1}; random max "
              f"{r['maxS_rand'] / math.sqrt(r['C']):.1f} sqrt(C) (rms {r['rms_rand'] / math.sqrt(r['C']):.2f} sqrt(C)) at "
              f"u={u2:+d},j={j2},k={k2} (typical {r['M'] / (2 * (r['p'] + 1) * (r['a'] + 1)):.0f})")
    # aggregate by amplification
    print()
    print(f"B3 (q = {q}): aggregates by distance to the critical line (clocks with C >= 200)")
    bins = [(1.0, 1.05), (1.05, 1.5), (1.5, 3.0), (3.0, 10.0), (10.0, 1e9)]
    print(f"  {'amp bin':>14} {'#':>4} {'med maxS/C':>11} {'med |S1|/C':>11} {'med Estr':>9} {'med mU/rt':>10} "
          f"{'med D*':>8} {'top-1 |u|<=3':>12} {'sum(Coll-C)':>12} {'sum C(C-1)/M':>13}")
    for lo, hi in bins:
        sel = [r for r in rows if r["C"] >= 200 and lo <= amp(r["p"], r["a"], q) < hi]
        if not sel:
            continue
        ms = [r.get("maxS", r.get("maxS_struct", float("nan"))) / r["C"] for r in sel]
        s1 = [r["S1"] / r["C"] for r in sel if "S1" in r]
        es = [r["struct_frac"] for r in sel if "struct_frac" in r]
        mu = [r["max_unstruct"] / math.sqrt(r["C"] * math.log(r["M"])) for r in sel if "max_unstruct" in r]
        ds = [r["dstar"] for r in sel]
        t1 = [abs(r["top"][0][2][0]) <= 3 for r in sel if "top" in r]
        ex = sum(r["coll"] - r["C"] for r in sel)
        er = sum(r["C"] * (r["C"] - 1) / r["M"] for r in sel)
        med = lambda v: float(np.median(v)) if v else float("nan")  # noqa: E731
        print(f"  [{lo:5.2f},{hi:7.3g}) {len(sel):>4} {med(ms):11.4f} {med(s1):11.4f} {med(es):9.3f} {med(mu):10.3f} "
              f"{med(ds):8.4f} {sum(t1):>5}/{len(t1):<6} {ex:>12} {er:13.1f}")


def random_baseline():
    """uniform model: C iid residues mod M; max_{h!=0}|S| / sqrt(C ln M) and D* for comparison."""
    print()
    print("B4: uniform-model baseline (C iid uniform residues mod M, same C and M as the clock)")
    rng = np.random.default_rng(20260925)
    for (p, a) in [(16, 10), (18, 11), (19, 12), (20, 12), (21, 13), (22, 14)]:
        M = abs(gate(p, a))
        C = math.comb(p, a)
        vals = []
        ds = []
        cs = []
        for _ in range(3):
            r = rng.integers(0, M, size=C)
            hist = np.bincount(r, minlength=M).astype(np.float64)
            absS = np.abs(np.fft.rfft(hist))
            vals.append(float(absS[1:].max()) / math.sqrt(C * math.log(M)))
            rs = np.sort(r)
            ds.append(disc_block(rs, M, 0, C))
            cs.append(runs_stats(rs)[0] - C)
        print(f"  ({p},{a}) M={M} C={C}: max|S|/sqrt(C ln M) = {np.mean(vals):.3f} (3 draws), "
              f"D* = {np.mean(ds):.5f}, Coll-C = {np.mean(cs):.0f} (mean C(C-1)/M = {C * (C - 1) / M:.0f})")


def main():
    lean_malloc()
    t = time.time()
    rows3 = run(3, 30, perigee_all_upto=26)
    perigee_theorems(rows3, 3)
    summaries(rows3, 3)
    rows5 = run(5, 24, perigee_all_upto=22)
    perigee_theorems(rows5, 5)
    summaries(rows5, 5)
    random_baseline()
    # machine-readable dump for part C (eligible primitive necklaces per clock)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scratch", "procgen_gates",
                       "partB_elig.tsv")
    try:
        with open(out, "w") as fh:
            for r in rows3 + rows5:
                if r["elig_prim"] is not None:
                    fh.write(f"{r['q']}\t{r['p']}\t{r['a']}\t{r['elig_prim'] // r['p']}\t{r['elig']}\t{r['C']}\n")
    except OSError:
        pass
    report_mem("end")
    print(f"[B total {time.time() - t:.1f}s]", file=sys.stderr)


if __name__ == "__main__":
    main()
