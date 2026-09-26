#!/usr/bin/env python3
"""
procgen_cube_20260925_run.py -- the strategy cube of Althofer's 3n+-1 game.
Session collatz-procgen-20260922, cube lane, 2026-09-25.

Runs everything and prints the .out, in this order:
  A    engine build and parameters
  B    levels 1-4: exact reference (procgen_cube_20260925_core.py) against the C engine
  C    level 5: all 65536 strategies (engine) + exact Karp densities re-derived in Python
  C'   certificates: (i) descent + finite check, (ii) trap expansion, exact integers
  D    classification tables per level (and the full level-3 table)
  C''  budget sensitivity of (iii) at level 5 (OPEN strategies re-searched with starts <= 2^20)
  E    the negation involution nu
  F    structure of the OPEN set (densities, drift, exceptional dimension)
  F2   where the (iii) witnesses live; near-critical giant cycles
  G    the provability boundary (Theorems A, B, C) checked on every strategy
  H    Collatz's neighbourhood (levels 2-5, exhaustive)
  H2   Hamming-ball searches around Collatz (levels 3-7)
  H3   exact distances delta_k (k <= 9), eps_k (k <= 8) by lazy MaxSAT (procgen_cube_20260925_boundary.py)
  F3   exceptional dimension: Perron formula vs exact counts (from procgen_cube_20260925_macro.py)
  H4   single-flip neighbours of Collatz, levels 2-9 (macro)
  I    microcosm -> macrocosm: sampled levels 6-12, drift concentration (I3), k = infinity model (I4) (macro)

Every check() raises on failure.  Heavy work is in the C engine; each exact claim is re-derived
in Python with integers/Fractions.  One process at a time; peak memory well below 700 MB.

Usage:  python3 procgen_cube_20260925_run.py > 05-knowledge/results/procgen_cube_20260925.out
Env:    CUBE_FAST=1 skips the slow parts of section I (for testing).
"""
import os
import sys
import math
import time
import subprocess
import gzip
from fractions import Fraction
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SCR = os.path.join(ROOT, 'scratch', 'procgen_cube', 'strategy_cube')   # caches; not to be committed
sys.path.insert(0, HERE)
import procgen_cube_20260925_core as core  # noqa: E402

T0 = time.time()
CRIT = math.log(2) / math.log(3)          # the critical odd density log_3 2
LOG2, LOG3 = math.log(2), math.log(3)
N5 = 131072                                # cycle-search budget at levels <= 5
CAP_EXP = 60                               # value cap 2^60 (engine)


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


def hdr(t):
    print()
    print("=" * 100)
    print(t)
    print("=" * 100)


# ============================================================================ engine
ENGINE = os.path.join(SCR, 'engine')


def build_engine():
    os.makedirs(SCR, exist_ok=True)
    src = os.path.join(HERE, 'procgen_cube_20260925_engine.c')
    subprocess.run(['cc', '-O3', '-march=native', '-o', ENGINE, src, '-lm'], check=True)


def parse_line(line):
    parts = line.split()
    d = {}
    for p in parts[1:]:
        key, _, val = p.partition('=')
        d[key] = val
    out = {'label': d.get('mask', d.get('sample'))}
    out['k'] = int(d['k'])
    out['mumax'] = float(d['mumax'])
    out['mumin'] = float(d['mumin'])
    out['ci'] = int(d['ci'])
    out['cii'] = int(d['cii'])
    out['nscc'] = int(d['nscc'])
    out['nbot'] = int(d['nbot'])
    out['npos'] = int(d['npos'])
    bots = []
    for b in d['bots'].split(';'):
        f = b.split(',')
        bots.append({'minnode': int(f[0]), 'size': int(f[1]), 'mn': float(f[2]), 'mx': float(f[3]),
                     'podd': float(f[4]), 'drift': float(f[5]), 'mass': float(f[6])})
    out['bots'] = bots
    out['lz'] = int(d['lz'])
    out['bad'] = [int(x) for x in d['bad'].split(',')]
    out['D'] = float(d['D'])
    if 'N' in d:
        out['N'] = int(d['N'])
        out['esc'] = int(d['esc'])
        out['unres'] = int(d['unres'])
        out['ncyc'] = int(d['ncyc'])
        cyc = []
        if d['cyc']:
            for c in d['cyc'].split(';'):
                mn, ln, od, bs = c.split(':')
                cyc.append((int(mn), int(ln), int(od), int(bs)))
        out['cyc'] = cyc
    return out


def run_engine(args, stdin_text=None, env_extra=None):
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    r = subprocess.run([ENGINE] + [str(a) for a in args], input=stdin_text, capture_output=True, text=True,
                       env=env, check=True)
    recs = []
    words = {}
    for line in r.stdout.splitlines():
        if line.startswith('S '):
            recs.append(parse_line(line))
        elif line.startswith('W '):
            _, lab, w = line.split()
            words[lab.split('=')[1]] = w
    return recs, words


# ============================================================================ exact per-strategy analysis
def exact_analysis(k, mask):
    """Exact (integer/Fraction) G-level data: densities, bottom SCCs, classes (i)/(ii)."""
    succ = core.graph(k, mask)
    comps = core.tarjan_scc(succ)
    bots = core.bottom_sccs(succ, comps)
    K = 1 << k
    rmax = core.karp_density(succ, range(K), True)
    rmin = core.karp_density(succ, range(K), False)
    bdat = []
    for c in bots:
        bmx = core.karp_density(succ, c, True)
        bmn = core.karp_density(succ, c, False)
        bdat.append({'nodes': c, 'rmax': bmx, 'rmin': bmn})
    ci = core.density_vs_critical(rmax) < 0
    cii = any(core.density_vs_critical(b['rmin']) > 0 for b in bdat)
    return {'succ': succ, 'comps': comps, 'bots': bdat, 'rmax': rmax, 'rmin': rmin, 'ci': ci, 'cii': cii}


def verify_cycle(k, mask, mn, ln, od):
    sig = core.sig_table(k, mask)
    x = mn
    odd = 0
    m = mn
    for _ in range(ln):
        odd += x & 1
        x = core.T(x, k, sig)
        m = min(m, x)
    return x == mn and m == mn and odd == od


def fr(ap):
    return f"{ap[0]}/{ap[1]}"


# ============================================================================ section A
def section_A():
    hdr("A. Setup")
    build_engine()
    print("  C engine built: scratch/procgen_cube/strategy_cube/engine  (source 04-computation/experiments/procgen_cube_20260925_engine.c)")
    print("  strategy sigma of level k: odd residues mod 2^k -> {+1,-1};  T(n) = n/2 (even), (3n + sigma(n mod 2^k))/2 (odd)")
    print("  mask bit i <-> residue 2i+1, bit 1 = minus.  Strings list sigma(1), sigma(3), ..., sigma(2^k-1).")
    print("  u/d word: at r, 'd' iff 4 | 3r + sigma(r) (>= 2 halvings follow), 'u' iff exactly one.")
    print(f"  critical odd density c = log_3 2 = {CRIT:.12f}; a cycle (a odd steps, length p) is EXPANDING iff 3^a > 2^p iff a/p > c")
    print(f"  cycle-search budget at levels <= 5: positive starts 1..{N5}, value cap 2^{CAP_EXP}, Brent detection")
    at = core.amin_table(401)
    for j in range(401):
        check(3 ** at[j] > 2 ** j and (at[j] == 0 or 3 ** (at[j] - 1) <= 2 ** j), "amin table")
    print("  exact table amin(j) = least a with 3^a > 2^j checked for j <= 401 (the engine uses the same table)")


# ============================================================================ section B
def section_B(store):
    hdr("B. Levels 1-4: the exact Python reference against the C engine")
    NB = {1: 20000, 2: 20000, 3: 20000, 4: 3000}
    for k in (1, 2, 3, 4):
        recs, _ = run_engine(['all', k, NB[k]], env_extra={'CUBE_DIM': '2'})
        nm = 1 << (1 << (k - 1))
        check(len(recs) == nm, f"engine record count level {k}")
        nmis = Counter()
        for mask, e in enumerate(recs):
            check(int(e['label'], 16) == mask, "mask order")
            ex = exact_analysis(k, mask)
            succ = ex['succ']
            # (1) classes (i)/(ii): exact Karp densities == engine float Karp == simple-cycle enumeration
            cyc = core.simple_cycles(succ)
            has_exp = any(core.expanding(*core.cycle_ap(c)) for c in cyc)
            check(ex['ci'] == (not has_exp), f"Theorem A test by cycle enumeration k={k} mask={mask}")
            check(ex['ci'] == bool(e['ci']), f"engine ci k={k} mask={mask}")
            unif = []
            for b in ex['bots']:
                cc = core.simple_cycles(succ, b['nodes'])
                unif.append(all(core.expanding(*core.cycle_ap(c)) for c in cc))
            check(ex['cii'] == any(unif) == bool(e['cii']), f"engine cii k={k} mask={mask}")
            check(abs(e['mumax'] - (Fraction(*ex['rmax']) * LOG3 - LOG2)) < 1e-9, "mumax value")
            check(abs(e['mumin'] - (Fraction(*ex['rmin']) * LOG3 - LOG2)) < 1e-9, "mumin value")
            # (2) bottom SCCs: exact stationary law vs engine
            check(len(ex['bots']) == e['nbot'], "bottom count")
            eb = sorted(e['bots'], key=lambda b: b['minnode'])
            xb = sorted(ex['bots'], key=lambda b: b['nodes'][0])
            npos = 0
            for b1, b2 in zip(eb, xb):
                check(b1['size'] == len(b2['nodes']) and b1['minnode'] == b2['nodes'][0], "bottom nodes")
                pi = core.stationary(succ, b2['nodes'])
                podd = sum(v for s, v in pi.items() if s & 1)
                check(abs(float(podd) - b1['podd']) < 1e-9, f"stationary odd mass k={k} mask={mask}")
                sgn = core.drift_sign_exact(podd)
                check((b1['drift'] > 0) == (sgn > 0), "drift sign")
                npos += sgn > 0
                check(abs(b1['mn'] - (Fraction(*b2['rmin']) * LOG3 - LOG2)) < 1e-9, "bottom min mean")
                check(abs(b1['mx'] - (Fraction(*b2['rmax']) * LOG3 - LOG2)) < 1e-9, "bottom max mean")
                b2['podd'] = podd
            check(npos == e['npos'], "npos")
            # (3) exceptional counts: DP (Python) == engine; brute force for small L
            dp = core.bad_counts_dp(k, mask, 40)
            eng = [e['bad'][L - 1] for L in range(1, 41)]
            check(dp[1:] == eng, f"Bad_L counts k={k} mask={mask}")
            if k <= 3:
                br = core.bad_counts_brute(k, mask, 9)
                check(br == dp[:10], "Bad_L brute force")
            elif mask % 16 == 0:
                br = core.bad_counts_brute(k, mask, 7)
                check(br == dp[:8], "Bad_L brute force (k=4 sample)")
            lz = next((L for L in range(1, 41) if dp[L] == 0), -1)
            if ex['ci'] and lz < 0:
                cert = core.descent_certificate(k, mask, 400)
                lz = cert[0]
            check(lz == e['lz'], "first empty level")
            # (4) cycle search: engine == Python reference (same N, same cap)
            cyc_py, esc_py, unres_py, _ = core.cycle_search(k, mask, NB[k], cap=1 << CAP_EXP)
            cyc_en = {c[0]: (c[1], c[2]) for c in e['cyc']}
            check(cyc_py == cyc_en, f"cycle lists k={k} mask={mask}: {cyc_py} vs {cyc_en}")
            check(esc_py == e['esc'] and unres_py == e['unres'], "escape counts")
            # (5) dimension: D = -inf iff (i); D = 1 iff a positive-drift bottom SCC (numerically)
            if ex['ci']:
                check(e['D'] == float('-inf'), "D for (i)")
            else:
                check(0 <= e['D'] <= 1 + 1e-6, "D range")
                check((abs(e['D'] - 1) < 1e-4) == (npos > 0), f"D = 1 iff positive drift k={k} mask={mask} D={e['D']}")
            nmis['ok'] += 1
            store[(k, mask)] = {'eng': e, 'ex': ex}
        print(f"  level {k}: {nmis['ok']} strategies; for every one of them:")
        print("     (i)-test [engine float Karp] == exact integer Karp on odd densities == simple-cycle enumeration")
        print("     (ii)-test likewise on every bottom SCC; bottom SCCs and stationary odd mass == exact Fractions")
        print("     |Bad_L| (L<=40): engine DP == Python DP" + (" == brute force (L<=9)" if k <= 3 else " (brute force L<=7 on 16 masks)"))
        print(f"     cycle lists and escape counts for starts <= {NB[k]}: engine == Python reference")
        print("     numeric exceptional dimension: -inf exactly on (i); = 1 exactly when a bottom SCC has positive drift")


# ============================================================================ section C
def section_C(store):
    hdr("C. Level 5: all 65536 strategies")
    path = os.path.join(SCR, f'lvl5_N{N5}.txt.gz')
    if not os.path.exists(path):
        t = time.time()
        r = subprocess.run([ENGINE, 'all', '5', str(N5)], capture_output=True, text=True, check=True)
        with gzip.open(path, 'wt') as f:
            f.write(r.stdout)
        print(f"  engine run: {time.time() - t:.0f} s")
    recs = []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('S '):
                recs.append(parse_line(line))
    check(len(recs) == 65536, "level-5 record count")
    print(f"  engine output parsed: {len(recs)} records (cache scratch/procgen_cube/strategy_cube/lvl5_N{N5}.txt.gz, not committed)")
    t = time.time()
    ncyc_ver = 0
    near = 0
    for mask, e in enumerate(recs):
        check(int(e['label'], 16) == mask, "mask order 5")
        ex = exact_analysis(5, mask)
        check(ex['ci'] == bool(e['ci']) and ex['cii'] == bool(e['cii']), f"level-5 classes by exact Karp, mask {mask}")
        check(len(ex['bots']) == e['nbot'], "bottoms 5")
        eb = sorted(e['bots'], key=lambda b: b['minnode'])
        xb = sorted(ex['bots'], key=lambda b: b['nodes'][0])
        for b1, b2 in zip(eb, xb):
            check(b1['size'] == len(b2['nodes']) and b1['minnode'] == b2['nodes'][0], "bottom nodes 5")
            if abs(b1['drift']) < 1e-9:     # exact fallback near the critical density
                near += 1
                pi = core.stationary(core.graph(5, mask), b2['nodes'])
                podd = sum(v for s, v in pi.items() if s & 1)
                check((core.drift_sign_exact(podd) > 0) == (b1['drift'] > 0), "exact drift sign near 0")
            b2['podd'] = b1['podd']
        for (mn, ln, od, bs) in e['cyc']:
            check(verify_cycle(5, mask, mn, ln, od), f"cycle {mn} of mask {mask}")
            ncyc_ver += 1
        check(e['unres'] == 0, "no unresolved starts")
        check(e['ncyc'] == len(e['cyc']), "cycle table overflow")
        if not ex['ci']:
            check(e['lz'] == -1 and e['bad'][39] > 0, "Theorem A: no certificate at L = 40")
        else:
            check(e['lz'] > 0, "(i) has an empty level")
        ex.pop('succ', None)
        ex.pop('comps', None)
        e['bad'] = tuple(e['bad'])
        store[(5, mask)] = {'eng': e, 'ex': ex}
    print(f"  exact integer Karp (max/min odd density, global and per bottom SCC) re-derived for all 65536: classes (i), (ii) agree")
    print(f"  every reported cycle re-verified exactly: {ncyc_ver} cycles; no unresolved start; drift within 1e-9 of 0: {near}")
    print(f"  ({time.time() - t:.0f} s)")


# ============================================================================ certificates for (i) and (ii)
def certificates(store):
    hdr("C'. Certificates: (i) descent (Terras DP + thresholds + finite check), (ii) trap expansion")
    for k in (1, 2, 3, 4, 5):
        keys = [key for key in store if key[0] == k]
        Lhist = Counter()
        n0max = 0
        cyc_i = Counter()
        for key in keys:
            s = store[key]
            if not s['ex']['ci']:
                continue
            mask = key[1]
            L, n0, nst = core.descent_certificate(k, mask, 400)
            check(L == s['eng']['lz'], "L_min = first empty level")
            found = core.cycles_below(k, mask, n0)
            s['cert_i'] = (L, n0, nst, found)
            Lhist[L] += 1
            n0max = max(n0max, n0)
            # the finite check gives ALL cycles; engine (starts <= N) must agree exactly
            eng_cyc = sorted(c[0] for c in s['eng']['cyc'])
            check(sorted(found) == eng_cyc, f"(i) cycles k={k} mask={mask}: {sorted(found)} vs {eng_cyc}")
            cyc_i[tuple(sorted(found))] += 1
        Mhist = Counter()
        Wmax = Fraction(0)
        for key in keys:
            s = store[key]
            if not s['ex']['cii']:
                continue
            mask = key[1]
            certs = []
            for b in s['ex']['bots']:
                if core.density_vs_critical(b['rmin']) > 0:
                    r = core.divergence_certificate(k, mask, b['nodes'])
                    check(r is not None, "(ii) certificate exists")
                    M, W = r
                    certs.append((b['nodes'], M, W))
                    Mhist[M] += 1
                    Wmax = max(Wmax, W)
            s['cert_ii'] = certs
        ni = sum(Lhist.values())
        nii = sum(1 for key in keys if store[key]['ex']['cii'])
        print(f"  level {k}: (i) {ni} strategies, L_min histogram {dict(sorted(Lhist.items()))}, max threshold n0 = {n0max};")
        print(f"            complete cycle lists of the (i) strategies (by min element): {dict(cyc_i)}")
        print(f"            (ii) {nii} strategies, block length M histogram {dict(sorted(Mhist.items()))}, max W = {float(Wmax):.4f}")


# ============================================================================ classification
def classify(store):
    for key, s in store.items():
        e, ex = s['eng'], s['ex']
        extra = sorted(c[0] for c in e['cyc'] if c[0] > 1)
        s['I'] = ex['ci']
        s['II'] = ex['cii']
        s['III'] = len(extra) > 0
        s['extra'] = extra
        lab = []
        if s['I']:
            lab.append('I')
        if s['II']:
            lab.append('II')
        if s['III']:
            lab.append('III')
        s['label'] = '+'.join(lab) if lab else 'IV'
        k, mask = key
        s['sig'] = core.sigma_string(k, mask)
        s['ud'] = core.ud_word(k, mask)
        s['prim'] = core.is_primitive(k, mask)
        s['nu'] = core.negate_mask(k, mask)
        s['selfdual'] = s['nu'] == mask
        s['pos_bottom'] = e['npos'] > 0


def section_D(store):
    hdr("D. Classification per level: (i) PROVED convergent, (ii) PROVED divergent, (iii) extra cycle found, (iv) OPEN")
    print("  (i)  all cycles of G_sigma contracting  <=>  Bad_L empty for some L (certificate + finite check)")
    print("  (ii) a bottom SCC of G_sigma all of whose cycles expand (trap certificate)")
    print(f"  (iii) a positive-integer cycle not through 1 found among starts <= {N5} (verified exactly)")
    print("  (iv) OPEN: none of the above;  a strategy may carry two labels (I+III, II+III)")
    tot = {}
    for k in (1, 2, 3, 4, 5):
        keys = sorted(key for key in store if key[0] == k)
        c_all = Counter(store[key]['label'] for key in keys)
        c_new = Counter(store[key]['label'] for key in keys if store[key]['prim'])
        nall, nnew = len(keys), sum(1 for key in keys if store[key]['prim'])
        tot[k] = c_all
        print(f"  level {k}: {nall} strategies ({nnew} new at this level)")
        for lab in ('I', 'I+III', 'II', 'II+III', 'III', 'IV'):
            print(f"     {lab:7s} all {c_all.get(lab, 0):6d} ({c_all.get(lab, 0) / nall:7.4f})   new {c_new.get(lab, 0):6d}")
        nI = sum(v for lab, v in c_all.items() if 'I' in lab.split('+'))
        nII = sum(v for lab, v in c_all.items() if 'II' in lab.split('+'))
        nIII = sum(v for lab, v in c_all.items() if 'III' in lab.split('+'))
        print(f"     totals: (i) {nI}  (ii) {nII}  (iii) {nIII}  (iv) {c_all.get('IV', 0)}")
    # level-3 full table
    print()
    print("  Level 3, every strategy (rho_min / pi_odd / rho_max per bottom SCC; D = exceptional dimension):")
    print("  sigma  ud    class   rho_max(G)  bottom SCC: size, rho_min, pi_odd, rho_max, drift/step   D       extra cycles (min)  nu-partner")
    for mask in range(16):
        s = store[(3, mask)]
        b = s['ex']['bots'][0]
        eb = s['eng']['bots'][0]
        print(f"  {s['sig']}   {s['ud']}  {s['label']:7s} {fr(s['ex']['rmax']):>6s}      {len(b['nodes']):2d}, {fr(b['rmin']):>5s}, "
              f"{float(b['podd']):.4f}, {fr(b['rmax']):>5s}, {eb['drift']:+.4f}   {s['eng']['D']:+.3f}   {str(s['extra'][:6]):18s}  "
              f"{core.sigma_string(3, s['nu'])}")
    return tot


# ============================================================================ symmetry
def section_E(store):
    hdr("E. The negation involution nu: sigma'(m) = -sigma(-m mod 2^k)")
    # square
    sq = {core.sigma_string(2, m): core.sigma_string(2, core.negate_mask(2, m)) for m in range(4)}
    print(f"  level 2: nu = {sq}")
    check(sq['++'] == '--' and sq['--'] == '++' and sq['+-'] == '+-' and sq['-+'] == '-+',
          "nu swaps Collatz and 3n-1 and fixes the mixed corners")
    print("  => nu exchanges Collatz (++) and 3n-1 (--), and fixes (+-) and (-+)  [CHECKED]")
    print("  In u/d coordinates nu is the reflection r -> -r of the word (tau'(r) = tau(-r)); self-dual <=> tau symmetric.")
    for k in (1, 2, 3, 4, 5):
        keys = [key for key in store if key[0] == k]
        nself = sum(1 for key in keys if store[key]['selfdual'])
        npairs = (len(keys) - nself) // 2
        expect_self = 0 if k == 1 else 2 ** (2 ** (k - 2))
        check(nself == expect_self, "self-dual count")
        # invariance of every G-level quantity under nu
        for key in keys:
            s = store[key]
            t = store[(k, s['nu'])]
            check(store[(k, t['nu'])] is s, "nu involution")
            for fld in ('ci', 'cii'):
                check(s['ex'][fld] == t['ex'][fld], "nu invariance of classes (i),(ii)")
            check(s['ex']['rmax'] == t['ex']['rmax'] and s['ex']['rmin'] == t['ex']['rmin'], "nu invariance of densities")
            check(s['eng']['bad'] == t['eng']['bad'], "nu invariance of |Bad_L|")
            b1 = sorted((b['size'], round(b['podd'], 9), round(b['mn'], 9), round(b['mx'], 9)) for b in s['eng']['bots'])
            b2 = sorted((b['size'], round(b['podd'], 9), round(b['mn'], 9), round(b['mx'], 9)) for b in t['eng']['bots'])
            check(b1 == b2, "nu invariance of bottom SCC data")
            if not math.isinf(s['eng']['D']) and s['eng']['D'] > -1e8 and t['eng']['D'] > -1e8:
                check(abs(s['eng']['D'] - t['eng']['D']) < 1e-6, "nu invariance of D")
        # OPEN closure
        pair_types = Counter()
        self_types = Counter()
        for key in keys:
            s = store[key]
            if s['selfdual']:
                self_types[s['label']] += 1
            elif key[1] < s['nu']:
                t = store[(k, s['nu'])]
                pair_types[tuple(sorted((s['label'], t['label'])))] += 1
        n_open = sum(1 for key in keys if store[key]['label'] == 'IV')
        n_open_partner_open = sum(1 for key in keys if store[key]['label'] == 'IV' and store[(k, store[key]['nu'])]['label'] == 'IV')
        print(f"  level {k}: {len(keys)} strategies = {nself} self-dual + {npairs} dual pairs"
              f"  [formula: 2^(2^(k-2)) self-dual]")
        print(f"     G-level data (classes (i),(ii), densities, bottom SCC data, |Bad_L|, D) identical on every nu-pair [CHECKED]")
        print(f"     self-dual strategies by class: {dict(self_types)}")
        print(f"     dual pairs by (class, class): {dict(sorted(pair_types.items()))}")
        print(f"     OPEN: {n_open}; with OPEN partner (or self-dual): {n_open_partner_open};"
              f" with partner in (iii) [sheet-obstructed]: {n_open - n_open_partner_open};"
              f" OPEN closed under nu: {n_open == n_open_partner_open}")


# ============================================================================ structure of OPEN
def section_F(store):
    hdr("F. Structure of the OPEN set: densities, drift, exceptional dimension")
    for k in (2, 3, 4, 5):
        keys = [key for key in store if key[0] == k]
        opn = [store[key] for key in keys if store[key]['label'] == 'IV']
        nonI = [store[key] for key in keys if not store[key]['I']]
        print(f"  level {k}: OPEN {len(opn)}")
        if not opn:
            continue
        posdrift = sum(1 for s in opn if s['pos_bottom'])
        allneg = len(opn) - posdrift
        print(f"     all bottom SCCs of negative drift (Collatz-like): {allneg};  some positive-drift bottom SCC (5n+1-like): {posdrift}")
        print(f"     exceptional set nonempty for every OPEN strategy (Theorem A; |Bad_40| > 0): "
              f"{all(s['eng']['bad'][39] > 0 for s in opn)}")
        Ds = sorted(s['eng']['D'] for s in opn if not s['pos_bottom'])
        if Ds:
            q = lambda f: Ds[min(len(Ds) - 1, int(f * len(Ds)))]
            print(f"     D over Collatz-like OPEN: min {Ds[0]:.4f}  q25 {q(.25):.4f}  median {q(.5):.4f}  q75 {q(.75):.4f}  max {Ds[-1]:.4f}")
        rmx = Counter(fr(s['ex']['rmax']) for s in opn)
        print(f"     rho_max(G) over OPEN (most common): {rmx.most_common(8)}")
        s1 = sum(1 for s in opn if s['sig'][0] == '-')
        sm1 = sum(1 for s in opn if s['sig'][-1] == '+')
        both = sum(1 for s in opn if s['sig'][0] == '-' or s['sig'][-1] == '+')
        print(f"     OPEN with sigma(1) = -1 (1 is an expanding fixed point): {s1};  with sigma(-1) = +1 (the Collatz -1 loop): {sm1};"
              f"  either: {both};  neither: {len(opn) - both}")
        # drift of (i) versus OPEN (per step, main bottom SCC = largest basin)

        def main_drift(s):
            return max(s['eng']['bots'], key=lambda b: b['mass'])['drift']
        di = [main_drift(store[key]) for key in keys if store[key]['I']]
        do = [main_drift(s) for s in opn if not s['pos_bottom']]
        if di and do:
            print(f"     drift/step of (i): [{min(di):+.4f}, {max(di):+.4f}];  of Collatz-like OPEN: [{min(do):+.4f}, {max(do):+.4f}]"
                  f"  -> overlap: {min(do) < max(di)}")
        dii = [max(b['drift'] for b in store[key]['eng']['bots']) for key in keys if store[key]['II']]
        dpo = [max(b['drift'] for b in s['eng']['bots']) for s in opn if s['pos_bottom']]
        if dii and dpo:
            print(f"     largest bottom drift of (ii): [{min(dii):+.4f}, {max(dii):+.4f}];  of positive-drift OPEN: [{min(dpo):+.4f}, {max(dpo):+.4f}]")
        # the most nearly provable OPEN strategies: least rho_max, least D
        cand = sorted((Fraction(*s['ex']['rmax']), s['eng']['D'], s['sig'], s['ud']) for s in opn if not s['pos_bottom'])
        print(f"     least rho_max(G) among Collatz-like OPEN: {cand[0][0]} (log_3 2 = {CRIT:.4f}); strategies attaining it:"
              f" {sum(1 for c in cand if c[0] == cand[0][0])}, e.g. {cand[0][2]} ({cand[0][3]}) with D = {cand[0][1]:.4f}")
        candD = sorted((s['eng']['D'], s['sig'], s['ud'], fr(s['ex']['rmax'])) for s in opn if not s['pos_bottom'])
        print(f"     least D among Collatz-like OPEN: {candD[0][0]:.4f} at {candD[0][1]} ({candD[0][2]}), rho_max {candD[0][3]}")
        # escapes as evidence
        esc_pos = [s['eng']['esc'] for s in opn if s['pos_bottom']]
        esc_neg = [s['eng']['esc'] for s in opn if not s['pos_bottom']]
        if esc_neg:
            print(f"     escapes past 2^60 among starts <= {N5}: Collatz-like OPEN max {max(esc_neg)}; "
                  f"positive-drift OPEN min {min(esc_pos) if esc_pos else '-'}")


def section_F2(store):
    hdr("F2. Where the (iii) witnesses live: smallest extra cycle, and the near-critical giants")
    for k in (3, 4, 5):
        keys = [key for key in store if key[0] == k]
        mins = [store[key]['extra'][0] for key in keys if store[key]['extra']]
        bins = Counter()
        for m in mins:
            b = '<=10' if m <= 10 else '<=100' if m <= 100 else '<=1000' if m <= 1000 else '<=10^4' if m <= 10 ** 4 else '>10^4'
            bins[b] += 1
        order = ['<=10', '<=100', '<=1000', '<=10^4', '>10^4']
        print(f"  level {k}: strategies with an extra cycle: {len(mins)}; smallest extra minimum: "
              + ", ".join(f"{b}: {bins.get(b, 0)}" for b in order))
        mc = Counter(m for key in keys for m in store[key]['extra'])
        print(f"     most common extra-cycle minima: {mc.most_common(10)}")
    keys = [key for key in store if key[0] == 5]
    big = Counter()
    bigs = defaultdict(set)
    nbig = 0
    for key in keys:
        for (mn, ln, od, bs) in store[key]['eng']['cyc']:
            if ln >= 300:
                big[(ln, od)] += 1
                bigs[(ln, od)].add(key)
                nbig += 1
                check(abs(od / ln - CRIT) < 3e-4, "near-critical giant")
    nstrat = len(set().union(*bigs.values())) if bigs else 0
    close5 = sum(1 for key in keys for (mn, ln, od, bs) in store[key]['eng']['cyc'] if ln >= 300 and abs(od / ln - CRIT) < 1e-5)
    print(f"  level 5: {nbig} cycles of length >= 300 in {nstrat} strategies; every one has |a/p - log_3 2| < 3e-4,"
          f" and {close5} have < 1e-5.")
    print("  (exact identity for any cycle: a log 3 - p log 2 = -sum over its odd points x of log(1 + sigma(x)/(3x)))")
    print("  most frequent (length p, odd steps a):")
    for (ln, od), c in sorted(big.items(), key=lambda x: -x[1])[:12]:
        g = math.gcd(ln, od)
        print(f"     p = {ln:5d}, a = {od:5d}  (= {g} x ({ln // g}, {od // g}))   3^a {'>' if 3 ** od > 2 ** ln else '<'} 2^p   "
              f"{c} cycles in {len(bigs[(ln, od)])} strategies")
    print("  (the convergents of log_2 3 are 1, 2, 3/2, 8/5, 19/12, 65/41, 84/53, 485/306, 1054/665, ...)")


def section_C3(store):
    hdr("C''. Budget sensitivity of (iii) at level 5: OPEN strategies re-searched with starts <= 2^20")
    N2 = 1 << 20
    path = os.path.join(SCR, f'lvl5_open_N{N2}.txt.gz')
    opn = [key[1] for key in sorted(store) if key[0] == 5 and store[key]['label'] == 'IV']
    if not os.path.exists(path):
        t = time.time()
        text = '\n'.join(format(m, 'x') for m in opn) + '\n'
        r = subprocess.run([ENGINE, 'list', '5', str(N2)], input=text, capture_output=True, text=True, check=True,
                           env=dict(os.environ, CUBE_DIM='0'))
        with gzip.open(path, 'wt') as f:
            f.write(r.stdout)
        print(f"  engine run: {time.time() - t:.0f} s")
    with gzip.open(path, 'rt') as f:
        recs = [parse_line(l) for l in f if l.startswith('S ')]
    check(len(recs) == len(opn), "sensitivity record count")
    flipped = []
    for m, e in zip(opn, recs):
        check(int(e['label'], 16) == m, "sensitivity order")
        extra = [c for c in e['cyc'] if c[0] > 1]
        for (mn, ln, od, bs) in extra:
            check(verify_cycle(5, m, mn, ln, od), "sensitivity cycle")
        if extra:
            flipped.append((m, min(c[0] for c in extra), extra))
    print(f"  {len(opn)} OPEN strategies at starts <= {N5}; with starts <= {N2}: {len(flipped)} acquire a (verified) extra cycle")
    for m, mn, extra in sorted(flipped, key=lambda x: x[1])[:10]:
        c = min(extra)
        print(f"     {core.sigma_string(5, m)}: extra cycle min {c[0]}, length {c[1]}, odd steps {c[2]}")
    esc = sum(1 for e in recs if e['esc'] > 0)
    print(f"  OPEN strategies with some start escaping past 2^60 (starts <= {N2}): {esc}")
    return len(flipped)


# ============================================================================ boundary theorems
def section_G(store):
    hdr("G. The provability boundary, checked on every strategy of levels 1-5")
    n = 0
    for key, s in store.items():
        ex, e = s['ex'], s['eng']
        n += 1
        # Theorem A: (i) <=> rho_max < c <=> some Bad_L empty
        check(s['I'] == (core.density_vs_critical(ex['rmax']) < 0) == (e['lz'] > 0), "Theorem A")
        # Theorem B: (ii) <=> some bottom SCC with rho_min > c ; then that SCC has positive drift
        check(s['II'] == any(core.density_vs_critical(b['rmin']) > 0 for b in ex['bots']), "Theorem B")
        if s['II']:
            check(s['pos_bottom'], "(ii) => positive-drift trap")
        # sandwich: rho_min(C) <= pi_odd(C) <= rho_max(C)
        for b in ex['bots']:
            po = b['podd']
            check(float(Fraction(*b['rmin'])) - 1e-12 <= float(po) <= float(Fraction(*b['rmax'])) + 1e-12, "sandwich")
        # (i) => every bottom drift negative, (i) and (ii) exclusive
        if s['I']:
            check(not s['pos_bottom'] and not s['II'], "(i) => negative drift")
    print(f"  checked on {n} strategies:")
    print("   Theorem A:  (i) <=> every cycle of G_sigma has odd density < log_3 2 <=> Bad_L = empty for some L")
    print("   Theorem B:  (ii) <=> some bottom SCC has every cycle of odd density > log_3 2  (=> that trap has positive drift)")
    print("   sandwich:   rho_min(C) <= pi_odd(C) <= rho_max(C) on every bottom SCC C (drift = pi_odd log 3 - log 2)")
    # the obstruction census for non-(i): which expanding cycles exist
    for k in (2, 3, 4, 5):
        keys = [key for key in store if key[0] == k]
        c = Counter()
        for key in keys:
            s = store[key]
            if s['I']:
                continue
            sig = s['sig']
            f1 = sig[0] == '-'
            fm1 = sig[-1] == '+'
            c[(f1, fm1)] += 1
        print(f"  level {k}: non-(i) by the two fixed-point obstructions (sigma(1)=-1, sigma(-1)=+1): "
              f"{ {('1' if a else '.') + ('-1' if b else '.'): v for (a, b), v in sorted(c.items())} }")


# ============================================================================ Collatz neighbourhood
def section_H(store):
    hdr("H. Collatz's neighbourhood: Hamming / Haar distance from all-plus to each class")
    for k in (2, 3, 4, 5):
        H = 1 << (k - 1)
        keys = [key for key in store if key[0] == k]
        best = {}
        for key in keys:
            s = store[key]
            d = bin(key[1]).count('1')
            for cl in ('I', 'II', 'III', 'IV'):
                has = (cl == 'IV' and s['label'] == 'IV') or (cl != 'IV' and cl in s['label'].split('+'))
                if has and (cl not in best or d < best[cl][0]):
                    best[cl] = (d, [])
                if has and best[cl][0] == d:
                    best[cl][1].append(s['sig'])
        singles = [store[(k, 1 << i)] for i in range(H)]
        sc = Counter(s['label'] for s in singles)
        print(f"  level {k} ({H} odd residues):")
        for cl in ('I', 'II', 'III'):
            if cl in best:
                d, lst = best[cl]
                print(f"     nearest ({cl.lower()}): Hamming {d} = Haar {d / H:.4f}; {len(lst)} strategies at that distance, e.g. {lst[:4]}")
        print(f"     single-flip neighbours of Collatz by class: {dict(sc)}")
        c = store[(k, 0)]
        print(f"     Collatz itself: class {c['label']}, rho_max {fr(c['ex']['rmax'])}, rho_min {fr(c['ex']['rmin'])}, D {c['eng']['D']:.6f}")


def section_H2():
    hdr("H2. Hamming-ball searches around Collatz at levels 6 and 7 (engine 'ball' mode, exhaustive within the radius)")
    res = {}
    for (k, r, which) in ((3, 8, 1), (4, 8, 1), (5, 8, 1), (6, 8, 1), (7, 6, 1), (3, 8, 2), (4, 8, 2), (5, 8, 2), (6, 9, 2)):
        path = os.path.join(SCR, f'ball{k}_{which}_r{r}.txt')
        if not os.path.exists(path):
            t = time.time()
            with open(path, 'w') as f:
                subprocess.run([ENGINE, 'ball', str(k), str(r), str(which)], stdout=f, check=True)
        lines = open(path).read().splitlines()
        hits = [l for l in lines if l.startswith('B ')]
        summ = [l for l in lines if l.startswith('BALL ')][0]
        H = 1 << (k - 1)
        if hits:
            d = int(hits[0].split('d=')[1])
            masks = [int(l.split()[1].split('=')[1], 16) for l in hits]
            # exact verification of every hit
            for m in masks:
                ex = exact_analysis(k, m)
                if which == 1:
                    check(ex['ci'], "ball hit is (i)")
                    cert = core.descent_certificate(k, m, 400)
                    check(cert is not None, "ball hit certificate")
                    found = core.cycles_below(k, m, cert[1])
                    check(sorted(found) == [1], "ball hit transitive")
                else:
                    check(ex['cii'], "ball hit is (ii)")
            res[(k, which)] = (d, len(masks))
            ex_str = core.sigma_string(k, masks[0])
            flips = [2 * i + 1 for i in range(H) if (masks[0] >> i) & 1]
            print(f"  level {k}, class ({'i' if which == 1 else 'ii'}): nearest at Hamming {d} = Haar {d / H:.4f} "
                  f"({len(masks)} strategies); e.g. flips at residues {flips} mod {2 ** k}")
            if which == 1:
                print(f"       each hit re-verified exactly: all cycles contracting, descent certificate, finite check -> only the cycle {{1,2}}")
        else:
            res[(k, which)] = (None, r)
            print(f"  level {k}, class ({'i' if which == 1 else 'ii'}): none within Hamming {r} (Haar {r / H:.4f}); {summ.split('tested=')[1]} strategies tested")
    return res


def section_H3():
    hdr("H3. Exact distances from Collatz to the provable classes (lazy-constraint MaxSAT / CP-SAT, solver-certified)")
    import procgen_cube_20260925_boundary as bd
    print("  delta_k = least number of odd residues mod 2^k at which sigma must differ from Collatz for class (i);")
    print("  eps_k   = the same for class (ii).  Lazy no-good clauses are necessary conditions, so the solver optimum is a")
    print("  lower bound, and the first optimum without offending cycles is optimal.  Every optimum is re-verified exactly.")
    print("  Lifting preserves the map, so delta_k / 2^(k-1) and eps_k / 2^(k-1) are non-increasing in k.")
    out_i, out_ii = {}, {}
    for k in (2, 3, 4, 5, 6, 7, 8, 9):
        r = bd.min_flips_to_I_rc2(k, timeout=3000)
        check(r['status'] == 'optimal', f"RC2 optimum (i) at level {k}")
        if k <= 7 and bd.cp_model is not None:
            r2 = bd.min_flips_to_I(k, timeout=1800)
            check(r2['status'] == 'optimal' and r2['delta'] == r['delta'], "CP-SAT agrees with RC2 on delta_k")
        ok, rmax, cert, cycles = bd.verify_I(k, r['mask'])
        check(ok and cert is not None and cycles == [1], "optimal strategy verified: (i), certificate, only the cycle {1,2}")
        H = 1 << (k - 1)
        out_i[k] = r['delta']
        print(f"  (i)  level {k}: delta = {r['delta']:3d} of {H:3d}  Haar {r['delta'] / H:.4f}   [{r['iterations']} iterations, "
              f"{r['clauses']} clauses, {r['seconds']:.0f} s]  witness verified: rho_max {rmax[0]}/{rmax[1]}, L_min {cert[0]}, n0 {cert[1]}")
        if k <= 7:
            print(f"       flips at residues {r['flips']}")
    for k in (2, 3, 4, 5, 6, 7, 8):
        r = bd.min_flips_to_II_rc2(k, timeout=3000)
        check(r['status'] == 'optimal', f"RC2 optimum (ii) at level {k}")
        ex = exact_analysis(k, r['mask'])
        check(ex['cii'], "optimal (ii) strategy verified by exact Karp on its bottom SCCs")
        H = 1 << (k - 1)
        out_ii[k] = r['delta']
        print(f"  (ii) level {k}: eps   = {r['delta']:3d} of {H:3d}  Haar {r['delta'] / H:.4f}   [{r['iterations']} iterations, "
              f"{r['clauses']} clauses, {r['seconds']:.0f} s]  witness verified by exact Karp")
        if k <= 6:
            print(f"       flips at residues {r['flips']}")
    for k in (3, 4, 5, 6, 7, 8):
        check(out_i[k + 1] / 2 ** k <= out_i[k] / 2 ** (k - 1) + 1e-12, "monotone (lift)")
    return out_i, out_ii


def main():
    store = {}
    section_A()
    section_B(store)
    section_C(store)
    certificates(store)
    classify(store)
    section_D(store)
    section_C3(store)
    section_E(store)
    section_F(store)
    section_F2(store)
    section_G(store)
    section_H(store)
    h2 = section_H2()
    h3i, h3ii = section_H3()
    for k in (3, 4, 5, 6):
        check(h3i[k] == h2[(k, 1)][0], "solver delta_k == exhaustive ball distance")
        check(h3ii[k] == (h2[(k, 2)][0] if h2[(k, 2)][0] is not None else h3ii[k]), "solver eps_k == exhaustive ball")
    check(h2[(6, 2)][0] is None and h3ii[6] > h2[(6, 2)][1], "eps_6 consistent with the exhaustive bound")
    check(h2[(7, 1)][0] is None and h3i[7] > h2[(7, 1)][1], "delta_7 consistent with the exhaustive lower bound")
    import procgen_cube_20260925_macro as macro
    macro.section_I(run_engine, check, hdr, store)
    hdr("DONE")
    print(f"  ALL CHECKS PASSED  ({time.time() - T0:.0f} s)")


if __name__ == '__main__':
    main()
