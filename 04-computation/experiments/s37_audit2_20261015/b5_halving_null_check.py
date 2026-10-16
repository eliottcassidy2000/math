#!/usr/bin/env python3
"""S37 audit 2, script b5: the halving null of collatz_bias_null_calibration_2p37_20261015.py.

 (A) The claims of HYP-9284 / note section 6 read off the session's JSON (12 replicates + Collatz row) and NPZ:
     ranges, "below every replicate" statements, energy ratios, the pre-registration's "4.6-22.9".
 (B) Properties of the session's null: t(2m) = -t(m) by construction (checked), odd signs: mean, lag correlations.
 (C) An independent null implementation (different 64-bit mixer, sequential exact pass) to 2^30, 4 seeds, plus the
     session's splitmix64 null re-implemented for seed 7919 (row 1 of the session NPZ) to check that the session's per-
     octave level energies are reproduced; statistics on octaves 14..29 compared with Collatz on the same octaves
     (Collatz values from the verified 2^37 census NPZ): window energies, profile sd and one-octave autocorrelation,
     level energies, naive t-counts over the 16 octaves 14..29.
Usage: python b5_halving_null_check.py null.json null.npz census_ext.npz
"""
import sys, json, math, time
import numpy as np
import numba as nb

J = json.load(open(sys.argv[1]))
NZ = np.load(sys.argv[2])
CEN = np.load(sys.argv[3])
c = J['collatz']
nl = J['null']
print(f'(A) session JSON: {len(nl)} null replicates; insample {J["insample"][0]}..{J["insample"][-1]}, out {J["out"]}')


def rng3(vals):
    v = sorted(vals)
    return f'{v[0]:.3g}-{np.median(v):.3g}-{v[-1]:.3g}'
rows = [
    ('naive |t|>5 / >8', f"{c['n_t5']} / {c['n_t8']}", rng3([x['n_t5'] for x in nl]) + ' / ' + rng3([x['n_t8'] for x in nl])),
    ('max naive |t|', f"{c['max_abs_t']:.2f}", rng3([x['max_abs_t'] for x in nl])),
    ('one sign in 19 octaves', c['n_same_sign_19'], rng3([x['n_same_sign_19'] for x in nl])),
    ('top5 signs kept of 15', c['top5']['signs_kept_of_15'], rng3([x['top5']['signs_kept_of_15'] for x in nl])),
    ('top5 passes all', c['top5']['passes_all'], sum(x['top5']['passes_all'] for x in nl)),
    ('all-channel slope', f"{c['all_follow']['slope']:.3f}", rng3([x['all_follow']['slope'] for x in nl])),
    ('lag-1 acf by level', [round(a, 3) for a in c['acf_lag1_by_level']],
     f"{min(min(x['acf_lag1_by_level']) for x in nl):.3f}..{max(max(x['acf_lag1_by_level']) for x in nl):.3f}"),
    ('profile sd', f"{c['profile']['sd']:.3f}", rng3([x['profile']['sd'] for x in nl])),
    ('profile acf 1 oct', f"{c['profile']['acf_1oct']:.3f}", rng3([x['profile']['acf_1oct'] for x in nl])),
    ('window mean in', f"{c['window']['mean_in']:.3f}", rng3([x['window']['mean_in'] for x in nl])),
    ('window mean out', f"{c['window']['mean_out']:.3f}", rng3([x['window']['mean_out'] for x in nl])),
    ('window max in', f"{c['window']['max_in']:.3f}", rng3([x['window']['max_in'] for x in nl])),
    ('window max out', f"{c['window']['max_out']:.3f}", rng3([x['window']['max_out'] for x in nl])),
]
for lab, a, b in rows:
    print(f'  {lab:24s} Collatz {a} | null {b}')
ins = list(range(14, 33)); ks = list(range(14, 37))
ii = [ks.index(k) for k in ins]
LE = NZ['level_energy']           # (13, 23 octaves, 5 levels); row 0 = Collatz
assert bool(NZ['first_is_collatz'])
lev_in = LE[:, ii, :].mean(axis=1)
print('  level energies 1..5 in sample (14..32): Collatz ' + ', '.join(f'{x:.3f}' for x in lev_in[0]) +
      ' | null medians ' + ', '.join(f'{np.median(lev_in[1:, q]):.3f}' for q in range(5)) +
      ' | null ranges ' + ', '.join(f'{lev_in[1:, q].min():.2f}-{lev_in[1:, q].max():.2f}' for q in range(5)))
oi = np.array([x['energy_ratio_out_in_by_level'] for x in nl]); le = np.array([x['energy_ratio_late_early_by_level'] for x in nl])
print('  out/in (34..36 over 14..32) by level: Collatz ' + ', '.join(f'{x:.3f}' for x in c['energy_ratio_out_in_by_level']) +
      ' | null min ' + ', '.join(f'{oi[:, q].min():.3f}' for q in range(5)) +
      ' | Collatz below every replicate: ' + str([bool(c['energy_ratio_out_in_by_level'][q] < oi[:, q].min()) for q in range(5)]))
print('  late/early (30..36 over 14..29): Collatz ' + ', '.join(f'{x:.3f}' for x in c['energy_ratio_late_early_by_level']) +
      ' | null min ' + ', '.join(f'{le[:, q].min():.3f}' for q in range(5)) +
      ' | Collatz below every replicate: ' + str([bool(c['energy_ratio_late_early_by_level'][q] < le[:, q].min()) for q in range(5)]))
wi = np.array([x['window']['mean_in'] for x in nl]); wo = np.array([x['window']['mean_out'] for x in nl])
print(f"  ratio null/Collatz window mean, in sample: {wi.min() / c['window']['mean_in']:.1f}-"
      f"{np.median(wi) / c['window']['mean_in']:.1f}-{wi.max() / c['window']['mean_in']:.1f}; out of sample: "
      f"{wo.min() / c['window']['mean_out']:.1f}-{np.median(wo) / c['window']['mean_out']:.1f}-{wo.max() / c['window']['mean_out']:.1f}")
print(f"  replicates with in-sample window mean > 2: {int((wi > 2).sum())} of {len(wi)}; out-of-sample window means sorted: "
      + ', '.join(f'{x:.2f}' for x in sorted(wo)))
# section 4 table of the note: level energies per octave group
LEc = LE[0]
for grp in ((14, 19), (20, 25), (26, 29), (30, 33), (34, 36)):
    sel = [ks.index(k) for k in range(grp[0], grp[1] + 1)]
    print(f'  Collatz level energies octaves {grp[0]}-{grp[1]}: ' + ', '.join(f'{LEc[sel, q].mean():.2f}' for q in range(5)))

# (B) properties of the session's null
@nb.njit
def splitmix_sign(m, seed):
    z = np.uint64(m) ^ seed
    z = z + np.uint64(0x9E3779B97F4A7C15)
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    z = z ^ (z >> np.uint64(31))
    return 1 - 2 * int(z & np.uint64(1))


@nb.njit
def my_sign(m, seed):
    # independent mixer (Stafford variant 13 style constants differ from splitmix64), middle bit
    z = (np.uint64(m) * np.uint64(0xD1B54A32D192ED03)) ^ seed
    z = (z ^ (z >> np.uint64(32))) * np.uint64(0xAEF17502108EF2D9)
    z = (z ^ (z >> np.uint64(29))) * np.uint64(0x9E3779B97F4A7C15)
    z = z ^ (z >> np.uint64(32))
    return 1 - 2 * int((z >> np.uint64(37)) & np.uint64(1))


@nb.njit
def tnull(n, seed, kind):
    v = 0
    m = n
    while (m & 1) == 0:
        m >>= 1
        v += 1
    s = splitmix_sign(m, seed) if kind == 0 else my_sign(m, seed)
    return s if (v & 1) == 0 else -s


@nb.njit
def props(seed, kind, M):
    tot = 0
    c2 = 0
    c3 = 0
    viol = 0
    for k in range(M):
        m = 2 * k + 1
        a = tnull(m, seed, kind)
        tot += a
        c2 += a * tnull(m + 2, seed, kind)
        c3 += a * tnull((3 * m + 1) // 2 if ((3 * m + 1) // 2) % 2 == 1 else (3 * m + 1) // 2, seed, kind)
        viol += tnull(2 * m, seed, kind) != -a
        viol += tnull(4 * m, seed, kind) != a
    return tot, c2, c3, viol
for kind, lab in ((0, 'session splitmix64'), (1, 'audit mixer')):
    tot, c2, c3, viol = props(np.uint64(7919), kind, 1 << 24)
    M = 1 << 24
    print(f'(B) {lab}, seed 7919, odd m < 2^25: mean s {tot / M:+.5f} (sd 1/sqrt(M) = {1 / math.sqrt(M):.5f}); '
          f'corr(s(m), s(m+2)) {c2 / M:+.5f}; corr(s(m), t((3m+1)/2)) {c3 / M:+.5f}; violations of t(2m) = -t(m), '
          f't(4m) = t(m): {viol}')

# (C) sequential exact pass to 2^LOG for a null (kind, seed) or Collatz values from the census
LOG = 30
grid = CEN['grid'].astype(np.int64)
g30 = grid[grid < (1 << LOG)]


@nb.njit
def null_pass(seed, kind, grid, LOGN):
    G = grid.shape[0]
    C = np.zeros(243, dtype=np.int64)
    snapC = np.zeros((G, 243), dtype=np.int64)
    snapH = np.zeros(G, dtype=np.int64)
    snapP2 = np.zeros(G, dtype=np.float64)
    S = 0
    C2 = 0
    p2 = 0.0
    gi = 0
    r = 1
    for n in range(1, 1 << LOGN):
        t = tnull(n, seed, kind)
        C[r] += t
        S += t
        if r % 3 == 2:
            C2 += t
        h = 3 * C2 - S
        p2 += float(h) * float(h)
        while gi < G and grid[gi] == n:
            snapC[gi, :] = C
            snapH[gi] = h
            snapP2[gi] = p2
            gi += 1
        r += 1
        if r == 243:
            r = 0
    return snapC, snapH, snapP2


def stats(grid, C, H, P2, k0=14, k1=30):
    N = grid.astype(np.float64)
    u = np.log2(N)
    cols, levs = [], []
    for q in range(1, 6):
        m, mp_ = 3 ** q, 3 ** (q - 1)
        Cq = np.stack([C[:, r::m].sum(axis=1) for r in range(m)], axis=1)
        Cp = np.stack([C[:, r::mp_].sum(axis=1) for r in range(mp_)], axis=1) if q > 1 else C.sum(axis=1)[:, None]
        for r in range(m):
            cols.append((3 * Cq[:, r] - Cp[:, r % mp_]) / np.sqrt(6 * N / m))
            levs.append(q)
    Gm = np.stack(cols, axis=1); levs = np.array(levs)
    ks = list(range(k0, k1))
    OM = np.stack([Gm[(u >= k) & (u < k + 1)].mean(axis=0) for k in ks])
    msq = np.stack([(Gm[(u >= k) & (u < k + 1)] ** 2).mean(axis=0) for k in ks])
    lev_ms = np.stack([msq[:, levs == q].mean(axis=1) for q in range(1, 6)], axis=1)
    mi = OM.mean(axis=0); si = OM.std(axis=0, ddof=1)
    ti = mi / (si / math.sqrt(len(ks)))
    g = H / np.sqrt(2 * N)
    sel = (u >= k0) & (u < k1)
    x = g[sel] - g[sel].mean()
    acf = float((x[64:] * x[:-64]).mean() / x.var())
    jj = np.arange(len(grid) - 64)
    E = (P2[jj + 64] - H[jj + 64] ** 2) - (P2[jj] - H[jj] ** 2)
    e = E / N[jj] ** 2
    uY = np.log2(N[jj])
    w = (uY >= k0) & (uY < k1 - 1 + 1e-9)
    return dict(n_t5=int((abs(ti) > 5).sum()), sd=float(g[sel].std()), acf=acf, wmean=float(e[w].mean()),
                wmax=float(e[w].max()), lev=[float(lev_ms[:, q].mean()) for q in range(5)], lev_ms=lev_ms)
# Collatz from the census NPZ (exact), same octaves 14..29
selg = grid < (1 << LOG)
Cc = CEN['classsums'][selg].astype(np.float64)
Hc = CEN['H'][selg].astype(np.float64)
P2c = np.array([float(x) for x in CEN['preH2_str'][selg]])
stc = stats(g30, Cc, Hc, P2c)
print(f"(C) octaves 14..29 (windows with start in [2^14, 2^29]): Collatz: naive |t|>5 (16 octaves) {stc['n_t5']}, profile sd "
      f"{stc['sd']:.3f}, acf(1 oct) {stc['acf']:.3f}, window mean {stc['wmean']:.3f}, max {stc['wmax']:.3f}, level energies "
      + ', '.join(f'{x:.2f}' for x in stc['lev']))
t0 = time.time()
C1, H1, P1 = null_pass(np.uint64(7919), 0, g30, LOG)
st1 = stats(g30, C1.astype(np.float64), H1.astype(np.float64), P1)
sess = NZ['level_energy'][1][:16]           # session replicate seed 7919, octaves 14..29
print(f'  session null seed 7919 re-implemented ({time.time() - t0:.0f} s): per-octave level energies 14..29 max |audit - session| = '
      f'{np.abs(st1["lev_ms"] - sess).max():.2e}; window mean {st1["wmean"]:.3f}, sd {st1["sd"]:.3f}, acf {st1["acf"]:.3f}')
for sd in (11, 22, 33, 44):
    C2_, H2_, P2_ = null_pass(np.uint64(sd * 1000003), 1, g30, LOG)
    st = stats(g30, C2_.astype(np.float64), H2_.astype(np.float64), P2_)
    print(f"  audit null seed {sd}: naive |t|>5 {st['n_t5']}, profile sd {st['sd']:.3f}, acf {st['acf']:.3f}, window mean "
          f"{st['wmean']:.3f} (x{st['wmean'] / stc['wmean']:.1f} Collatz), max {st['wmax']:.3f}, level energies "
          + ', '.join(f'{x:.2f}' for x in st['lev']))
