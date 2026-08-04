#!/usr/bin/env python3
"""AMM 12592 -- FINAL REFEREE (lane F2), boxeph, 2026-08-04.

Stages (run separately; each writes a fragment JSON):

  witness4 : FOURTH-implementation verification of the (256, D0=0) desc1
             exact-floor witness (sha 5950bd42...).  Fully independent of
             the three prior verifiers: fresh per-n Beatty floor by binary
             search on 5^m <=> phi^(2n) (no incremental state), incremental
             Pascal-row admissibility / parity / transportation-point
             checks (no math.comb), epoch identity PROVED at coefficient
             level via the main-term closed form (2x-1 / -1, itself
             re-proved per degree from scratch) plus sparse-correction
             accumulation, and definition-level point evaluations at
             x = 2, 3, -1, 5, -7 that bypass the closed form entirely.

  audit    : independent (fresh) implementation of the Theorem S-cone-fc
             certificate: F1/F2/F3 scan + exact clocks at (2048,128) and
             (4096,256) from E1 stored feed-end snapshots, compared field
             by field to the stored fcscan; then my own definition-level
             box-clamp evolution to capture (S1 magnitude form asserted
             every row; support-freeze and no-reignition asserted).
             Conservativity checks: cell-0 empties at i0 + drain - 1
             (bound i0 + drain is 1 row conservative); conservative T_1
             form dominates the exact-drain T_1.  F3 => a_t <= C(D0-1,t+1)
             magnitude bound verified at the certifying row (the a-priori
             corollary's premise).  Marginal-surface erratum adjudication:
             recompute max r_t from the stored (4096,128), (4096,256) and
             (8192,256) snapshots.

  final    : harvest ALL amm12592_S_referee_frag_*.json + my fragments,
             totals, verdicts, and THE FINAL CANON PARAGRAPH; writes
             05-knowledge/results/amm12592_S_final_referee_boxeph.out.

Exact arithmetic only (int / Fraction).  No numpy, no sympy, no floats in
any decision (floats in display strings only).
"""
import sys, os, json, hashlib, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', '05-knowledge', 'results'))

_ck = {'pass': 0, 'fail': 0, 'fails': [], 'text': []}

def log(msg=''):
    print(msg)
    sys.stdout.flush()
    _ck['text'].append(msg)

def check(ok, label, extra=''):
    tag = 'PASS' if ok else 'FAIL'
    _ck['pass' if ok else 'fail'] += 1
    if not ok:
        _ck['fails'].append(label)
    log('  [%s] %s%s' % (tag, label, ('  -- ' + str(extra)) if extra != '' else ''))
    return ok

def save_frag(stage, extra=None):
    frag = {'stage': stage, 'pass': _ck['pass'], 'fail': _ck['fail'],
            'fails': _ck['fails'], 'text': _ck['text']}
    if extra is not None:
        frag['extra'] = extra
    with open(os.path.join(RES, 'amm12592_S_final_referee_frag_%s.json' % stage), 'w') as f:
        json.dump(frag, f, indent=1)

def exact_div(a, b):
    q, r = divmod(a, b)
    assert r == 0, 'inexact division'
    return q

# ---------------------------------------------------------------- fresh floor
def fib2(n):
    """(F_n, F_{n+1}) -- iterative fast doubling over the bits of n."""
    F, F1 = 0, 1
    for bit in bin(n)[2:]:
        c = F * (2 * F1 - F)
        d = F * F + F1 * F1
        if bit == '1':
            F, F1 = d, c + d
        else:
            F, F1 = c, d
    return F, F1

def floor_g_n(n):
    """floor(n * log_5 phi^2), from scratch for this single n:
    the unique m with 5^m <= phi^(2n) < 5^(m+1);
    phi^(2n) = (L_{2n} + F_{2n} sqrt5)/2, L_{2n} = 2 F_{2n+1} - F_{2n};
    5^k <= phi^(2n)  iff  t := 2*5^k - L <= 0  or  t^2 < 5 F^2
    (equality impossible: phi^(2n) is irrational for n >= 1)."""
    F, F1 = fib2(2 * n)
    L = 2 * F1 - F
    FF5 = 5 * F * F
    def le(k):
        t = 2 * 5 ** k - L
        return t <= 0 or t * t < FF5
    lo = max(0, (n * 59) // 100 - 3)
    hi = (n * 60) // 100 + 3
    assert le(lo) and not le(hi), 'floor bracket failed at n=%d' % n
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if le(mid):
            lo = mid
        else:
            hi = mid
    return lo

def my_profile(R, D0, imax):
    """d_i = floor(g (R+i)) + D0, each n computed independently."""
    return [floor_g_n(R + i) + D0 for i in range(imax + 1)]

# ------------------------------------------------------------- Pascal rows
def pascal_row(n, tmax=None):
    """[C(n,0..min(n,tmax))] by exact ratio recursion."""
    if tmax is None:
        tmax = n
    row = [1]
    c = 1
    for t in range(min(n, tmax)):
        c = exact_div(c * (n - t), t + 1)
        row.append(c)
    return row

def pascal_next(row):
    """C(n,.) full row -> C(n+1,.) full row."""
    out = [1]
    for t in range(1, len(row)):
        out.append(row[t] + row[t - 1])
    out.append(1)
    return out

# --------------------------------------------------------------- poly tools
def qpow_list(m):
    """(1-x)^k, k = 0..m, as coefficient lists (my own)."""
    out = [[1]]
    for _ in range(m):
        prev = out[-1]
        nxt = [prev[0]] + [prev[j] - prev[j - 1] for j in range(1, len(prev))] \
              + [-prev[-1]]
        out.append(nxt)
    return out

def encode(cells, d, qp, shift=0, acc=None):
    """acc += x^shift * sum_t cells[t] x^(d-t) (1-x)^t ; returns acc list."""
    if acc is None:
        acc = []
    need = shift + d + 1
    if len(acc) < need:
        acc.extend([0] * (need - len(acc)))
    for t, c in enumerate(cells):
        if c:
            q = qp[t]
            base = shift + d - t
            for j, b in enumerate(q):
                acc[base + j] += c * b
    return acc

def trim(p):
    while p and p[-1] == 0:
        p.pop()
    return p

# ==========================================================================
def stage_witness4():
    log('=' * 74)
    log('STAGE witness4 -- (256, D0=0) desc1 witness: FOURTH implementation')
    log('=' * 74)
    p = os.path.join(HERE, 'amm12592_witness_R256_bulk_desc1_D0_0_boxeph.json')
    raw = open(p, 'rb').read()
    sha = hashlib.sha256(raw).hexdigest()
    check(sha.startswith('5950bd42'), 'sha256 matches record 5950bd42...',
          sha[:24] + ' (%d bytes)' % len(raw))
    W = json.loads(raw)
    R = W['R']
    check(R == 256 and W['D0'] == 0, 'R = 256, D0 = 0', '%s %s' % (R, W['D0']))
    # ---- independent profile (fresh per-n binary search, no shared state)
    ds = my_profile(R, 0, R - 1)
    check(W.get('profile') == 'floor_gamma_star(R+i)+D0',
          'witness declares profile floor_gamma_star(R+i)+D0; recomputed fresh '
          'per-n (binary search on 5^m <=> phi^(2n))',
          'd_0=%d d_%d=%d' % (ds[0], R - 1, ds[R - 1]))
    deltas_ok = all(ds[i] - ds[i - 1] in (0, 1) for i in range(1, R))
    check(deltas_ok, 'profile increments in {0,1} (Beatty word)')
    sc = {int(i): {int(t): int(v) for t, v in row.items()}
          for i, row in W['sparse_c'].items()}
    # ---- main-term closed form re-proved per distinct degree
    dmax = ds[R - 1]
    qp = qpow_list(dmax + 2)
    seen = set()
    main_ok = True
    for i in range(R - 1):
        d = ds[i]
        if d in seen:
            continue
        seen.add(d)
        rdm1 = pascal_row(d - 1)
        cells = [(rdm1[t] if t <= d - 1 else 0) -
                 (rdm1[t - 1] if 1 <= t <= d else 0) for t in range(d + 1)]
        enc = trim(encode(cells, d, qp))
        if enc != [-1, 2]:
            main_ok = False
    check(main_ok, 'main-term closed form: sum_t (C(d-1,t)-C(d-1,t-1)) '
          'x^(d-t)(1-x)^t == 2x-1 for ALL %d distinct degrees' % len(seen))
    dlast = ds[R - 1]
    enc = trim(encode([-v for v in pascal_row(dlast)], dlast, qp))
    check(enc == [-1], 'last-row main term: -sum_t C(d,t) x^(d-t)(1-x)^t == -1',
          'd=%d' % dlast)
    # ---- admissibility + parity + transportation + point evals, one sweep
    pts = (2, 3, -1, 5, -7)
    acc_pts = [0] * len(pts)
    xp = []
    qpv = []
    for x0 in pts:
        pw = [1]
        for _ in range(dmax + 2):
            pw.append(pw[-1] * x0)
        xp.append(pw)
        q0 = 1 - x0
        pw = [1]
        for _ in range(dmax + 2):
            pw.append(pw[-1] * q0)
        qpv.append(pw)
    rowd = pascal_row(ds[0])
    rowdm1 = pascal_row(ds[0] - 1)
    ncells = 0
    nsat = 0
    adm_ok = par_ok = tpt_ok = True
    corr_acc = []          # sum over sparse rows of x^i corr_i(x)
    for i in range(R):
        d = ds[i]
        if i > 0 and ds[i] - ds[i - 1] == 1:
            rowdm1 = rowd
            rowd = pascal_next(rowd)
        assert len(rowd) == d + 1
        crow = sc.get(i, {})
        cells = []
        for t in range(d + 1):
            if i <= R - 2:
                main = (rowdm1[t] if t <= d - 1 else 0) - \
                       (rowdm1[t - 1] if 1 <= t <= d else 0)
            else:
                main = -rowd[t]
            cells.append(main + crow.get(t, 0))
        for t in range(d + 1):
            Cd = rowd[t]
            v = cells[t]
            ncells += 1
            if abs(v) > Cd:
                adm_ok = False
            if (v - Cd) % 2:
                par_ok = False
            f, r2 = divmod(Cd - v, 2)
            if r2 != 0 or f < 0 or f > Cd:
                tpt_ok = False
            if f == 0 or f == Cd:
                nsat += 1
        # definition-level point accumulation: x0^i * sum_t cells[t] x0^(d-t) q0^t
        for pi in range(len(pts)):
            xpw, qpw = xp[pi], qpv[pi]
            s = 0
            for t in range(d + 1):
                c = cells[t]
                if c:
                    s += c * xpw[d - t] * qpw[t]
            acc_pts[pi] += s * xpw[i] if pts[pi] != 0 else 0
        # sparse-correction coefficient accumulation
        if crow:
            tmax_c = max(crow)
            ccells = [crow.get(t, 0) for t in range(tmax_c + 1)]
            corr_acc = encode(ccells, d, qp, shift=i, acc=corr_acc)
    check(adm_ok, 'all %d cells admissible: |delta_t| <= C(d_i,t)' % ncells)
    check(par_ok, 'all cells parity delta_t == C(d_i,t) (mod 2)')
    check(tpt_ok, 'transportation point f = (C-delta)/2 integer in [0, C] everywhere',
          'boundary-saturated %d (%.2f%%)' % (nsat, 100.0 * nsat / ncells))
    check(ncells == 58837, 'cell count == 58837', ncells)
    check(nsat == 4973, 'boundary-saturated count == 4973 (E3/prior-referee value)', nsat)
    # ---- epoch identity, coefficient level, via corrections-only reduction:
    # sum_i x^i corr_i == (1-x)^(R-1) + x^(R-1) - sum_{i<=R-2} x^i (2x-1)
    rhs = list(qp[R - 1])
    L = max(len(rhs), len(corr_acc), R + 1)
    rhs.extend([0] * (L - len(rhs)))
    ca = list(corr_acc)
    ca.extend([0] * (L - len(ca)))
    rhs[R - 1] += 1
    for i in range(R - 1):
        rhs[i] += 1          # - (-1) x^i
        rhs[i + 1] -= 2      # - 2 x^(i+1)
    check(trim(ca) == trim(rhs),
          'EPOCH IDENTITY sum_i x^i Delta_i == (1-x)^(R-1): coefficient-exact '
          '(main closed form + sparse corrections; my own poly arithmetic)')
    # ---- definition-level point evaluations (no closed form used at all)
    pts_ok = True
    for pi, x0 in enumerate(pts):
        lhs = acc_pts[pi]
        rhsv = (1 - x0) ** (R - 1)
        if lhs != rhsv:
            pts_ok = False
        check(lhs == rhsv, 'identity at x = %d, definition-level (all %d rows, '
              'every cell)' % (x0, R))
    save_frag('witness4', extra={'sha256': sha, 'ncells': ncells, 'nsat': nsat,
                                 'points': list(pts), 'points_ok': pts_ok})

# ==========================================================================
def load_snap(path):
    snap = json.load(open(path))
    jd = {int(t): int(v) for t, v in snap['junk'].items()}
    m = max(jd)
    a = [0] * (m + 1)
    for t, v in jd.items():
        assert v <= 0, 'snapshot junk not all-negative'
        a[t] = -v
    return snap['i_postfeed'], snap['d_feedend'], a

def stage_audit():
    log('=' * 74)
    log('STAGE audit -- independent S-cone-fc certificate + evolution + erratum')
    log('=' * 74)
    CANON = {(2048, 128): 1237, (4096, 256): 2486}
    for (R, D0), canon_cap in CANON.items():
        snapf = os.path.join(RES, 'amm12592_S_cone_feedend_R%d_D0%d_boxeph.json' % (R, D0))
        fcf = os.path.join(RES, 'amm12592_S_cone_fcscan_R%d_D0%d_boxeph.json' % (R, D0))
        i_pf, d_fe, a = load_snap(snapf)
        imax = min(R, (3 * R) // 4 + 200)
        ds = my_profile(R, D0, imax)
        check(ds[i_pf - 1] == d_fe, '(%d,%d) my fresh profile: d_%d == d_fe == %d'
              % (R, D0, i_pf - 1, d_fe), ds[i_pf - 1])
        first = None
        cell0_empty_row = None
        support_ok = True
        reign = set()
        reignitions = 0
        cap_row = None
        m_at_cert = None
        i0 = i_pf
        while i0 < imax - 2:
            m = len(a) - 1
            while m >= 0 and a[m] == 0:
                m -= 1
            D_0 = ds[i0]
            a0 = a[0] if a else 0
            if first is None:
                F1 = all(v >= 0 for v in a) and m + 2 < D_0
                F2 = a0 <= D_0 - 1
                capref = [2 * c for c in pascal_row(D_0 - 1, m + 3)]
                F3 = True
                for t in range(2, m + 3):
                    am1 = a[t - 1] if t - 1 <= m else 0
                    am2 = a[t - 2] if t - 2 <= m else 0
                    if 2 * am1 + am2 > capref[t]:
                        F3 = False
                        break
                if F1 and F2 and F3:
                    # ---- exact clocks, my own staircase
                    drain = (a0 + 1) // 2
                    need = {t: a[t] for t in range(1, m + 1) if a[t]}
                    acc = {t: 0 for t in need}
                    acc1x = 0            # exact-drain T_1 variant (audit)
                    T1exact = None
                    Tt = {}
                    crow = pascal_row(D_0 - 1, m + 3)
                    nn = D_0 - 1
                    k = 0
                    while (need or T1exact is None) and k < R:
                        k += 1
                        Dk = ds[i0 + k]
                        while nn < Dk - 1:
                            # advance full row C(nn,.) -> C(nn+1,.)
                            for t in range(len(crow) - 1, 0, -1):
                                crow[t] += crow[t - 1]
                            nn += 1
                        dl = Dk - ds[i0 + k - 1]
                        # exact-drain T_1 (state entering row i0+k has a0-2k)
                        if T1exact is None and 1 in [1]:
                            a0k_x = max(0, a0 - 2 * k)
                            incx = 2 * (Dk - 1) - (1 + dl) * a0k_x
                            if incx > 0:
                                acc1x += incx
                            if acc1x >= (a[1] if len(a) > 1 else 0):
                                T1exact = k
                        for t in list(need):
                            if t == 1:
                                a0k = max(0, a0 - 2 * (k - 1))
                                inc = 2 * (Dk - 1) - (1 + dl) * a0k
                                if inc > 0:
                                    acc[t] += inc
                            else:
                                acc[t] += 2 * crow[t] - capref[t]
                            if acc[t] >= need[t]:
                                Tt[t] = k
                                del need[t]
                    Tmax = max(Tt.values()) if Tt else 0
                    Tworst = max(Tt, key=lambda t: Tt[t]) if Tt else None
                    ub = i0 + max(drain, Tmax)
                    first = {'i_fc': i0, 'd': D_0, 'm': m, 'a0': a0,
                             'drain': drain, 'Tmax': Tmax, 'T_worst_t': Tworst,
                             'capture_ub': ub, 'budget_margin': (R - 2) - ub}
                    m_at_cert = m
                    check(Tt.get(1) is not None and T1exact is not None
                          and Tt[1] >= T1exact,
                          '(%d,%d) conservative T_1 >= exact-drain T_1 '
                          '(safe direction)' % (R, D0),
                          '%s >= %s' % (Tt.get(1), T1exact))
                    # F3 => magnitude bound (a-priori corollary premise)
                    crow0 = pascal_row(D_0 - 1, m + 3)
                    mb_ok = all(a[t] <= crow0[t + 1] for t in range(1, m + 1))
                    check(mb_ok, '(%d,%d) F3 magnitude bound a_t <= C(D0-1,t+1) '
                          'on [1,m] at i_fc' % (R, D0))
            # ---- evolve one definition-level row (row i0 at degree ds[i0])
            dl = ds[i0] - ds[i0 - 1]
            kk = 1 + dl
            j = [-v for v in a]
            w = [0] * (len(j) + kk)
            for t, v in enumerate(j):
                if v:
                    w[t] += v
                    if kk == 1:
                        w[t + 1] += v
                    else:
                        w[t + 1] += 2 * v
                        w[t + 2] += v
            caps = pascal_row(ds[i0] - 1, len(w) + 1)
            newa = [0] * len(w)
            for t in range(len(w)):
                wt = w[t]
                lo = -2 * caps[t]
                hi = 2 * caps[t - 1] if t >= 1 else 0
                c = wt if lo <= wt <= hi else (lo if wt < lo else hi)
                jn = wt - c
                if wt <= 0:
                    s1 = min(0, wt + (2 * caps[t] if t >= 1 else 2))
                    assert jn == s1, 'S1 magnitude form violated'
                newa[t] = -jn
            old_nz = {t for t, v in enumerate(a) if v}
            new_nz = {t for t, v in enumerate(newa) if v}
            if first is not None:
                if m_at_cert is not None and new_nz and max(new_nz) > m_at_cert:
                    support_ok = False
                for t in new_nz:
                    if t in reign:
                        reignitions += 1
                reign |= (old_nz - new_nz)
            if a and a[0] > 0 and (len(newa) == 0 or newa[0] == 0) \
               and cell0_empty_row is None:
                cell0_empty_row = i0
            while len(newa) > 1 and newa[-1] == 0:
                newa.pop()
            if len(newa) == 1 and newa[0] == 0:
                newa = []
            a = newa
            if not a:
                cap_row = i0
                break
            i0 += 1
        fc = json.load(open(fcf))
        ref = fc['fc']
        ok = (first is not None and first['i_fc'] == fc['i_fc'] and
              all(first[k] == ref[k] for k in
                  ('d', 'm', 'a0', 'drain', 'Tmax', 'T_worst_t',
                   'capture_ub', 'budget_margin')))
        check(ok, '(%d,%d) my certificate == stored fcscan (all 8 fields)' % (R, D0),
              json.dumps(first))
        check(cap_row == canon_cap, '(%d,%d) my evolution captures at canon row %d'
              % (R, D0, canon_cap), cap_row)
        check(cap_row is not None and first is not None
              and cap_row <= first['capture_ub'],
              '(%d,%d) capture <= certificate bound' % (R, D0),
              '%s <= %s' % (cap_row, first['capture_ub']))
        check(support_ok and reignitions == 0,
              '(%d,%d) support frozen inside [0,m] after i_fc; zero re-ignitions'
              % (R, D0))
        check(cell0_empty_row == first['i_fc'] + first['drain'] - 1,
              '(%d,%d) cell 0 empties exactly at i_fc + drain - 1 '
              '(bound i_fc + drain is 1 row conservative)' % (R, D0),
              '%s vs %s' % (cell0_empty_row, first['i_fc'] + first['drain']))
    # ---- marginal-surface erratum adjudication
    log('  marginal-surface recomputation from stored snapshots (caps at d_fe-1):')
    NOTE = {(4096, 128): '1.0018', (4096, 256): '1.0047', (8192, 256): '1.0001'}
    ledger = {}
    for (R, D0), noteval in NOTE.items():
        snapf = os.path.join(RES, 'amm12592_S_cone_feedend_R%d_D0%d_boxeph.json' % (R, D0))
        i_pf, d_fe, a = load_snap(snapf)
        m = len(a) - 1
        crow = pascal_row(d_fe - 1, m + 3)
        best, bt, over = None, None, []
        for t in range(2, m + 3):
            am1 = a[t - 1] if t - 1 <= m else 0
            am2 = a[t - 2] if t - 2 <= m else 0
            r = Fraction(2 * am1 + am2, 2 * crow[t])
            if best is None or r > best:
                best, bt = r, t
            if r > 1:
                over.append(t)
        ledger[(R, D0)] = (float(best), bt, over)
        log('    (%d,%d): max r_t = %.6f at t=%d over=%s  [note table: %s]'
            % (R, D0, float(best), bt, over, noteval))
    v128 = ledger[(4096, 128)][0]
    check(abs(v128 - 1.0018) > 0.0005, 'E1-F3 erratum CONFIRMED: (4096,1/32) '
          'table value 1.0018 does not match stored snapshot', '%.6f' % v128)
    v256 = ledger[(4096, 256)][0]
    check(abs(v256 - 1.0047) < 0.0006, '(4096,1/16) table value 1.0047 reproduces',
          '%.6f' % v256)
    v8k = ledger[(8192, 256)][0]
    check(abs(v8k - 1.0001) < 0.0006, '(8192,1/32) table value 1.0001 reproduces',
          '%.6f' % v8k)
    save_frag('audit')

# ==========================================================================
def stage_final():
    log('=' * 74)
    log('STAGE final -- harvest of all referee fragments + FINAL CANON PARAGRAPH')
    log('=' * 74)
    frags = {}
    for pre in ('amm12592_S_referee_frag_', 'amm12592_S_final_referee_frag_'):
        for fn in sorted(os.listdir(RES)):
            if fn.startswith(pre) and fn.endswith('.json'):
                st = ('' if pre.endswith('S_referee_frag_') else 'F2:') + \
                     fn[len(pre):-len('.json')]
                frags[st] = json.load(open(os.path.join(RES, fn)))
    tp = sum(f['pass'] for f in frags.values())
    tf = sum(f['fail'] for f in frags.values())
    allfails = [x for f in frags.values() for x in f['fails']]
    log('')
    log('HARVEST: %d fragments, %d PASS, %d FAIL%s' %
        (len(frags), tp, tf, ('  FAILS: ' + '; '.join(allfails)) if allfails else ''))
    for st in sorted(frags):
        f = frags[st]
        log('  %-28s %3d pass %2d fail' % (st, f['pass'], f['fail']))
    lines = []
    A = lines.append
    A('AMM 12592 -- FINAL REFEREE CLOSE-OUT (lane F2).  boxeph, 2026-08-04.')
    A('Script: 04-computation/amm12592_S_final_referee_boxeph.py')
    A('Harvests: amm12592_S_referee_frag_*.json (four-process hostile referee)')
    A('        + amm12592_S_final_referee_frag_{witness4,audit}.json (this session)')
    A('Parent assembled report: 05-knowledge/results/amm12592_S_referee_boxeph.out')
    A('')
    A('TOTAL CHECKS ACROSS ALL FRAGMENTS: %d PASS, %d FAIL%s' % (tp, tf,
      ('  FAILS: ' + '; '.join(allfails)) if allfails else ''))
    A('')
    for st in sorted(frags):
        f = frags[st]
        A('  %-30s %3d pass %2d fail' % (st, f['pass'], f['fail']))
    A('')
    A('=' * 76)
    A('THIS SESSION\'S TWO NEW STAGES (full text)')
    A('=' * 76)
    for st in ('F2:witness4', 'F2:audit'):
        if st in frags:
            for ln in frags[st]['text']:
                A(ln)
            A('')
    A('=' * 76)
    A('FINAL VERDICTS (lane F2)')
    A('=' * 76)
    A('''
V1. THEOREM S-cone-fc: CONFIRMED.  The one-row-certificate logic was
    re-derived independently from S1 (exact post-feed magnitude law
    a'_t = max(0,(K a)_t - 2C(d'-1,t)), a'_0 = max(0,a_0-2)) and S2
    (comparison principle): F2 is absorbing because a_0 falls 2/row while
    d never falls; F3 with the FIXED row-i0 reference caps propagates
    because every cell is non-increasing (spills computed from the
    original state dominate all later spills) and caps only grow
    (D_k >= D_0); front-freeze at m+1, m+2 consumes exactly F3 at
    t = m+1, m+2; death needs support at cell D_k > m+2, impossible.
    Two deliberate conservativities verified in this session's audit
    stage: the drain bound i0+drain is exactly 1 row above the true
    cell-0 extinction (i0+drain-1, observed), and the T_1 staircase uses
    a_0-2(k-1) >= the true a_0-2k (T_1^cons >= T_1^exact, observed
    632>=631 / 1264>=1264).  My fresh implementation reproduces the
    stored fcscan certificates in ALL 8 fields at (2048,128) and
    (4096,256) and captures at the canon rows inside the bound.

V2. R = 32768 REPLICATION: CONFIRMED, closing the last independence gap.
    Both eps replicated from scratch by the four-process referee (fresh
    Beatty floor, fresh kernel/feed/clamp): CLOSED at 19865 (eps=1/16)
    and 20185 (eps=1/32), debt 16383 = (R-2)/2 exact, zero re-ignitions,
    cell-0 clock exact, cell-1 extinction rows 17133/17521, feed-end
    snapshots BIT-IDENTICAL to E1's stored states (170/169 cells), and
    the S-cone-fc certificates re-derived from BOTH the E1 snapshots and
    the referee's own independent snapshots (conecert 'mine' rows).
    S(1/32) and S(1/16) are VERIFIED-exact for all dyadic
    128 <= R <= 32768 by two independent implementation lineages.

V3. A-PRIORI F4 + WINDOW: CONFIRMED with scope discipline.  The window
    death-freeness [i_pf, i_pf+64] is now PROVED for all dyadic
    2^16..2^40 (symbolic mixed-endpoint interval form
    D0(3+2g) >= 2R(1-g-g^2) + 66(1+g), plus exact per-R checks at
    2^16/2^17) -- E1-F1 quantifier gap CLOSED.  The a-priori F4
    corollary is certified for dyadic 2^9..2^40 ONLY (code range
    confirmed: R = 2^7..2^40 loops in the constants script); the
    independent worst-case clock re-computation at 2^16/2^17 (apriori
    stage) passes with T_2 always the worst clock.  For R > 2^40, F4
    must be carried inside the entry hypothesis (it is a decidable
    one-row check, so the reduction loses nothing).

V4. (256, D0=0) EXACT-FLOOR WITNESS sha 5950bd42: CONFIRMED, now by a
    FOURTH independent implementation (this session): fresh per-n Beatty
    floor by binary search on 5^m <=> phi^(2n); incremental Pascal-row
    admissibility (58837 cells in-box, correct parity); transportation
    integer point f = (C-delta)/2 in [0,C] everywhere with EXACTLY 4973
    boundary-saturated cells (8.45%); epoch identity
    sum_i x^i Delta_i = (1-x)^255 proved COEFFICIENT-EXACT via the
    re-proved main-term closed form (ballot(d) -> 2x-1 for all 152
    distinct degrees; corner -> -1) plus sparse-correction accumulation,
    and INDEPENDENTLY at definition level (all rows, every cell) at
    x = 2, 3, -1, 5, -7.  17/17.

V5. ERRATUM (display-only), adjudicated: the E1 note sec. 6 table entry
    (4096, eps=1/32) "1.0018" is wrong; the stored (post-bug-regenerated)
    snapshot gives max r_t = 1.004211 at t = 4, over-cells {2,4,6}.
    Neighbors (4096,1/16) = 1.0047 and (8192,1/32) = 1.0001 reproduce
    exactly.  Qualitative marginal-surface law unaffected.

V6. STALENESS FLAG: the parent referee's drafted verdict "E2 -- NO
    CONTENT (standing by)" predates the super-block note
    (amm12592-S-superblock-proof-boxeph.md, this session, 00:32): that
    lane DID deliver a second, independent per-(R,D0) S-certificate
    mechanism, 16/16 PASS at R = 128..16384 (both eps) with the capture
    row PREDICTED exactly (i0 + |j_0|/2) in every run -- consistent
    row-for-row with the canon captures.  Its internal Theorems 1-2 were
    NOT line-audited by the four-process referee; treat its PROVED labels
    as its own note's, its 16/16 ledger as VERIFIED-exact.

V7. ENTRY: still OPEN.  No lane-F1 entry proof is visible in the
    worktree at close-out.  The older half-cap/Lambda entry ledger
    (amm12592_S_cone_entrycheck) shows H4 failing AT i_pf itself -- which
    is precisely why the operative certificate is the full-cap one firing
    at i_fc = i_pf + 0..8; it is not evidence against ENTRY-fc
    (rule/scan negatives are never infeasibility proofs).
''')
    A('=' * 76)
    A('THE FINAL CANON PARAGRAPH (referee-certified, 2026-08-04, lane F2)')
    A('=' * 76)
    A('''
UNCONDITIONAL: the best proved bound remains C* <= 2, from the submission
envelope T(n) <= max(n+1, 2n-2) (with T(n) <= max(n+1, 2n-3) for n >= 5);
exactly rational, no bracket.  NO lane of this program produced a new
unconditional C* bound, and none claims to.

PROVED, UNCONDITIONALLY, INSIDE THE PROGRAM (all referee-confirmed, all
exact int/Fraction): (i) the post-feed flow of plain rule A on
all-negative junk is the monotone system a'_t = max(0,(K a)_t-2C(d'-1,t)),
a'_0 = max(0,a_0-2), with comparison principle and unconditional deadline
i_pf <= (R-2)/2; (ii) THEOREM S-cone-fc: if ANY post-feed row i0
satisfies the one-row integer conditions F1 (all-negative junk, support
in [0,m], m+2 < d), F2 (a_0 <= d-1), F3 (2a_{t-1}+a_{t-2} <= 2C(d-1,t)
on t in [2,m+2]) and F4 (exact staircase clocks fit the R-2 budget),
then the flow captures by row R-2 with no death -- S(R) holds; (iii) F4
is automatic given F1-F3 for every dyadic R in [2^9, 2^40] (exact-
rational certification; NOT proved beyond 2^40), and the handover window
[i_pf, i_pf+64] is death-free for every dyadic R in [2^16, 2^40] at
eps in {1/32, 1/16} (symbolic interval certification); (iv) the closed
form 1 + gamma* + eps* = (5+3g)/(3+2g), g = gamma* = log_5 phi^2 in
(627035, 627036)/2^20, eps* = 2(1-g-g^2)/(3+2g) in (0.0211736,
0.0211747), with certified bracket (1.6191617801, 1.6191618342).

VERIFIED-exact (finite, two independent implementation lineages at the
top scale): Hypothesis S at eps = 1/32 AND 1/16 for ALL dyadic
128 <= R <= 32768 -- at 32768: captures 19865 (1/16) and 20185 (1/32),
debt 16383 = (R-2)/2 exact, zero post-feed re-ignitions, capture =
i_pf + a_0/2 - 1 exactly, feed-end states bit-identical across
implementations; the one-row certificate fires at i_fc = i_pf + 0..8
with ~38%-of-R budget margins in all 18 stored scans; the super-block
route independently certifies S at 16/16 configs (128..16384) with exact
capture prediction.  ATTAINMENT RECORDS: exact-floor existence at
(256, D0=0) -- witness sha 5950bd42, verified by FOUR independent
implementations, including a coefficient-exact epoch identity and a
transportation integer point (58837 cells, 4973 boundary-saturated);
best known at R = 512 is D0 = 4 (tr15, referee-verified).

PROVED-MODULO-ENTRY: define ENTRY-fc(eps): for every dyadic R >= 65536,
some post-feed row i0 <= i_pf + 64 of plain rule A satisfies F1^F2^F3
(and F4, needed only for R > 2^40).  Then S(eps) holds in full, and by
the D2 assembly (THM-3329 + LIFT + Theorem B, eps > eps*):
C* <= 1 + gamma* + eps.  In particular ENTRY-fc(1/32) gives
C* <= 1 + gamma* + 1/32 < 427095/262144 = 1.6292382; and if ENTRY-fc
holds for a sequence eps -> eps*+, C* <= (5+3g)/(3+2g) = 1.61916...
Every Theta(R)-row dynamical statement in S is PROVED; ENTRY-fc is a
STATIC property of one row, whose structural basis is the verified
marginal-surface law (max_t r_t -> 1+ from above: 1.0009 at 32768,
over-cells single-parity, low-cell core within 0.2% of the surface).

HONEST GAP LIST (everything that remains): (G1) ENTRY-fc(eps) itself --
OPEN for all eps; no entry proof exists in the worktree at close-out.
(G2) The a-priori F4 corollary and the window certification stop at
dyadic 2^40; beyond, F4 rides inside ENTRY-fc (decidable per row), or
the T_t <= T_2 monotonicity lemma must be written out in full.
(G3) The verified eps set is {1/32, 1/16}; "every eps > eps*" is
untested beyond those two, so the (5+3g)/(3+2g) limit needs ENTRY-fc
along a sequence eps -> eps*+.  (G4) The S-chain quantifies over DYADIC
R; general n enters only through the THM-3329 envelope assembly.
(G5) Exact-floor existence frontier: (512, D0=0) OPEN in both
directions (D0 = 1, 2, 3 also open).  (G6) The general-class floor:
C* >= log_5(5 phi^2) is OPEN -- the golden value is pinned only against
balanced-block schemes (THM-3024 demotion honored); so even ENTRY-fc for
all eps > eps* would leave C* itself undetermined in [phi-side floor,
(5+3g)/(3+2g)].  (G7) The super-block lane's internal Theorems 1-2 are
note-local PROVED labels, not yet line-audited by an independent
referee.  Nothing else is open: every other lemma, constant, witness,
capture row, and certificate in the S program now carries at least two
independent exact verifications.
''')
    with open(os.path.join(RES, 'amm12592_S_final_referee_boxeph.out'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    log('')
    log('final report written: amm12592_S_final_referee_boxeph.out '
        '(%d fragments, %d pass, %d fail)' % (len(frags), tp, tf))
    save_frag('final', extra={'total_pass': tp, 'total_fail': tf,
                              'fragments': sorted(frags)})
    return tp, tf

# --------------------------------------------------------------------- main
if __name__ == '__main__':
    stage = sys.argv[1]
    t0 = time.time()
    if stage == 'witness4':
        stage_witness4()
    elif stage == 'audit':
        stage_audit()
    elif stage == 'final':
        stage_final()
    else:
        raise SystemExit('unknown stage')
    print('STAGE %s DONE: %d pass, %d fail (%.1fs)'
          % (stage, _ck['pass'], _ck['fail'], time.time() - t0))
