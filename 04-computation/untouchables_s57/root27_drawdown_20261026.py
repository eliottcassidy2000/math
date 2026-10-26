#!/usr/bin/env python3
"""Pre-registered R-P1 and the drawdown law on the root-27 (Antihydra) ray (opus-2026-10-26-S57).

Section of the root-27 packet: T(z) = (3z + c)/2, c = z mod 2, credit q' = q + 3c - 1, source (q, z) = (1, 7); the machine
halts iff q reaches 0.  At odd z the V-return consumes ell = v_2(3z+1) symbols (one 1, ell-1 zeros), q' = q + 3 - ell.
R-P1 (pre-registered, commit 568cb689a4): over the first J = 10^5 returns, the largest deficit excursion
max_{i<j} (K_j - K_i - 3(j-i)), K_J = sum of ell, lies in [15, 35].
Also: the future-drop law.  Under Haar measure on the source the carry bits are fair and independent (root-27 (C3)),
so from any symbol time the walk ever drops by >= d with probability phi^-d (p = 1/2 + p^3/2).  We tabulate the share
of symbol times t <= T/2 whose future drop D_t = q_t - min_{u>t} q_u is >= d, against phi^-d.
Exact integers throughout: z is a Python int; 64 carry bits are read from z mod 2^64 (they depend only on it), then
z <- (3^64 z + C_64)/2^64 exactly (the division is asserted exact).
Usage: python root27_drawdown_20261026.py [symbols]
"""
import sys, os, json, math, time

HERE = os.path.dirname(os.path.abspath(__file__))
PHI = (1 + 5 ** 0.5) / 2
P3_64 = 3 ** 64
MASK = (1 << 64) - 1


def carry_word(z, nsym, q0=1, stop_at_halt=True):
    """Return list of carry bits and credit after each symbol; stops when q would be blocked (q = 1, c = 0)."""
    bits = []
    qs = [q0]
    q = q0
    done = 0
    while done < nsym:
        w = z & MASK
        C = 0
        blk = []
        for i in range(64):
            c = w & 1
            if stop_at_halt and q == 1 and c == 0:
                # the blocked symbol: the macro would halt here; replay exactly to this point and stop
                for j, cj in enumerate(blk):
                    pass
                z_part = z
                for cj in blk:
                    assert (z_part & 1) == cj
                    z_part = (3 * z_part + cj) >> 1
                return bits, qs, z_part, True
            blk.append(c)
            C += c * (1 << i) * 3 ** (63 - i)
            w = (3 * w + c) >> 1
            q += 3 * c - 1
            bits.append(c)
            qs.append(q)
        num = P3_64 * z + C
        assert num & MASK == 0
        z = num >> 64
        done += 64
    return bits, qs, z, False


def returns_from_bits(bits):
    """ell sequence of V-returns: each return = a 1 followed by its run of 0s."""
    ells = []
    i = 0
    n = len(bits)
    while i < n:
        assert bits[i] == 1
        j = i + 1
        while j < n and bits[j] == 0:
            j += 1
        if j == n:
            break
        ells.append(j - i)
        i = j
    return ells


def max_excursion(ells):
    best, run_min, cum = 0, 0, 0
    # deficit walk D_J = K_J - 3J; excursion = max_{i<j} D_j - D_i
    lo = 0
    for e in ells:
        cum += e - 3
        best = max(best, cum - lo)
        lo = min(lo, cum)
    return best


if __name__ == '__main__':
    nsym = int(float(sys.argv[1])) if len(sys.argv) > 1 else 2 * 10 ** 6
    out = {}
    t0 = time.time()
    # positive controls (hostiles): sources that halt
    for src in (15, 5, 9, 13, 27):
        b, qs, zf, halted = carry_word(src, 4096, q0=1)
        out[f'control_source_{src}'] = {'halted': halted, 'macros_before_halt': len(b),
                                        'tree_word': ''.join(map(str, b)) + '0' if halted else None}
    assert out['control_source_15']['tree_word'] == '1111000010000010100000', out['control_source_15']
    print('controls (q0 = 1):', out, flush=True)
    bits, qs, zf, halted = carry_word(7, nsym, q0=1)
    assert not halted
    assert bits[:4] == [1, 1, 1, 0], bits[:8]
    ells = returns_from_bits(bits)
    print(f'source 7: {len(bits)} symbols, {len(ells)} V-returns, min credit {min(qs[1:])}, final credit {qs[-1]}, '
          f'final z has {zf.bit_length()} bits ({time.time() - t0:.0f}s)', flush=True)
    rows = {}
    for J in (10 ** 3, 10 ** 4, 10 ** 5, 3 * 10 ** 5, 10 ** 6):
        if J <= len(ells):
            ex = max_excursion(ells[:J])
            rows[J] = {'max_deficit_excursion': ex, 'haar_scale_lnJ_over_lnphi': math.log(J) / math.log(PHI)}
            print(f'   J = {J:>8}: max deficit excursion {ex}   (ln J / ln phi = {math.log(J) / math.log(PHI):.1f})', flush=True)
    out['excursions'] = rows
    rp1 = rows.get(10 ** 5, {}).get('max_deficit_excursion')
    out['R-P1'] = {'value': rp1, 'window': [15, 35], 'pass': rp1 is not None and 15 <= rp1 <= 35}
    print('R-P1 (J = 1e5, window [15, 35]):', out['R-P1'], flush=True)
    # future-drop law over symbol times t <= T/2
    T = len(qs) - 1
    fut_min = [0] * (T + 1)
    m = qs[T]
    for t in range(T, -1, -1):
        fut_min[t] = m
        m = min(m, qs[t])
    half = T // 2
    hist = {}
    for t in range(half):
        d = qs[t] - fut_min[t]
        if d > 0:
            hist[d] = hist.get(d, 0) + 1
    tab = []
    tail = 0
    dmax = max(hist) if hist else 0
    cnt_ge = {}
    for d in range(dmax, 0, -1):
        tail += hist.get(d, 0)
        cnt_ge[d] = tail
    for d in range(1, min(dmax, 30) + 1):
        share = cnt_ge[d] / half
        tab.append((d, cnt_ge[d], share, PHI ** (-d)))
        if d <= 24 or d == dmax:
            print(f'   drop >= {d:>2}: share {share:.6f}   phi^-d {PHI ** (-d):.6f}   ratio {share / PHI ** (-d):.3f}', flush=True)
    out['future_drop_law'] = {'symbol_times': half, 'rows': tab, 'max_drop_seen': dmax}
    out['symbols'] = len(bits)
    out['returns'] = len(ells)
    out['min_credit'] = min(qs[1:])
    json.dump(out, open(os.path.join(HERE, 'root27_drawdown.json'), 'w'), indent=1)
    print(f'done ({time.time() - t0:.0f}s)')
