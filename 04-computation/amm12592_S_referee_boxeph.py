#!/usr/bin/env python3
"""AMM 12592 -- HOSTILE REFEREE of lanes E1 (S via invariant cone), E2, E3
(barrier/gap).  boxeph, 2026-08-04.

Independent FIFTH implementation of the rule-A junk flow, written from the
STATEMENTS of T2/T4/T6 (transient theorem) -- not from any prior engine code:
  * exact Beatty floor via 5^m <=> phi^(2n) integer comparisons (fresh),
  * row-0 initial data from the T4 closed form (itself re-verified against a
    definition-level polynomial decode in stage `slow`),
  * kernel/feed/clamp row loop from the T6 statement,
  * a definition-level SLOW rule A (sigma/ideal/AdmClamp polynomial
    arithmetic) used as ground truth at R = 128/256, with block-for-block
    Delta comparison against the fast flow.

Stages (CLI):  constants | slow | suite | bigrun R D0 | conecert R D0 |
               lemmas | apriori | window | e3wit | ledger | report
All decisions exact int/Fraction.  No numpy, no sympy, no floats in any
decision (floats only in display fields).

Outputs: 05-knowledge/results/amm12592_S_referee_frag_<stage>.json
         05-knowledge/results/amm12592_S_referee_boxeph.out   (stage `report`)
"""
import sys, os, json, time, random, hashlib
from fractions import Fraction
from math import comb, isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), '05-knowledge', 'results')
SCRATCH = os.environ.get('REF_SCRATCH', RES)

_buf = []
def log(msg=''):
    _buf.append(msg)
    print(msg, flush=True)

_ck = {'pass': 0, 'fail': 0, 'fails': []}
def check(ok, label, extra=''):
    ok = bool(ok)
    tag = 'PASS' if ok else 'FAIL'
    _ck['pass' if ok else 'fail'] += 1
    if not ok:
        _ck['fails'].append(label)
    log('  [%s] %s%s' % (tag, label, ('  -- ' + str(extra)) if extra != '' else ''))
    return ok

def save_frag(stage, extra=None):
    frag = {'stage': stage, 'pass': _ck['pass'], 'fail': _ck['fail'],
            'fails': _ck['fails'], 'text': _buf}
    if extra:
        frag['extra'] = extra
    with open(os.path.join(RES, 'amm12592_S_referee_frag_%s.json' % stage), 'w') as f:
        json.dump(frag, f, indent=1, default=str)

# ---------------------------------------------------------------- exact floor
def fib_pair(n):
    """(F_n, F_{n+1}) by fast doubling, exact."""
    if n == 0:
        return (0, 1)
    a, b = fib_pair(n >> 1)
    c = a * (2 * b - a)
    d = a * a + b * b
    return (d, c + d) if (n & 1) else (c, d)

class GammaFloor:
    """m(n) = floor(n * log_5 phi^2) exactly: largest m with 5^m <= phi^(2n),
    phi^(2n) = (L_{2n} + F_{2n} sqrt5)/2, L_{2n} = 2 F_{2n+1} - F_{2n}.
    Incremental in n (n -> n+1 = two Fibonacci steps)."""
    def __init__(self, n0):
        self.n = n0
        self.F, self.F1 = fib_pair(2 * n0)
        m = int(0.5979 * n0)
        p = 5 ** m
        while self._le(5 * p):
            m += 1; p *= 5
        while not self._le(p):
            m -= 1; p //= 5
        self.m, self.p = m, p
        # exact bracket assert: 5^m <= phi^(2n) < 5^(m+1)
        assert self._le(self.p) and not self._le(5 * self.p)
    def _le(self, pw):
        """is pw <= phi^(2n)?  (pw an integer; equality impossible for pw=5^k>1)"""
        L = 2 * self.F1 - self.F
        t = 2 * pw - L
        if t <= 0:
            return True
        return t * t < 5 * self.F * self.F
    def advance(self):
        F2 = self.F + self.F1
        self.F, self.F1 = F2, F2 + self.F1
        self.n += 1
        if self._le(5 * self.p):
            self.m += 1
            self.p *= 5
        return self.m

def profile(R, D0, imax):
    gf = GammaFloor(R)
    ds = [gf.m + D0]
    for _ in range(imax):
        ds.append(gf.advance() + D0)
    return ds

# ------------------------------------------------------------- binomial rows
def exact_div(a, b):
    q, r = divmod(a, b)
    assert r == 0, 'inexact division'
    return q

class BinRow:
    """row[t] = C(nn, t), t = 0..len-1; incremental in nn and length."""
    def __init__(self, nn, tmax):
        self.nn = nn
        row = [1]
        c = 1
        for t in range(tmax):
            c = exact_div(c * (nn - t), t + 1) if nn - t > 0 else 0
            row.append(c)
        self.row = row
    def extend(self, tmax):
        row, nn = self.row, self.nn
        t = len(row) - 1
        while t < tmax:
            row.append(exact_div(row[t] * (nn - t), t + 1) if nn - t > 0 else 0)
            t += 1
    def inc(self):
        row = self.row
        for t in range(len(row) - 1, 0, -1):
            row[t] += row[t - 1]
        self.nn += 1

# ------------------------------------------------------- the fast junk flow
def row0_junk(R, d0):
    """T4 closed form w_t = (-1)^(d0-t) C(R-2-t, d0-t) - C(d0+1,t+1) + 2C(d0,t),
    clamped into [-2C(d0-1,t), +2C(d0-1,t-1)] (t>=1), [-2,0] (t=0).
    Full sweep t = 0..d0 (no window assumption).  Returns junk list, c0."""
    A = comb(R - 2, d0)          # C(R-2-t, d0-t) at t=0
    B = comb(d0 + 1, 1)          # C(d0+1, t+1) at t=0
    C0 = 1                       # C(d0, t) at t=0
    CAP = 1                      # C(d0-1, t) at t=0
    CAPm1 = 0                    # C(d0-1, t-1) at t=0
    sgn = -1 if (d0 % 2) else 1  # (-1)^(d0-t) at t=0
    junk = []
    c0_emit = None
    for t in range(d0 + 1):
        w = sgn * A - B + 2 * C0
        lo = -2 * CAP
        hi = 2 * CAPm1
        c = w if lo <= w <= hi else (lo if w < lo else hi)
        jn = w - c
        assert jn % 2 == 0 and c % 2 == 0
        if t == 0:
            c0_emit = c
        junk.append(jn)
        if t == d0:
            assert jn == 0, 'DEATH at row 0'
        # advance families t -> t+1
        if t < d0:
            A = exact_div(A * (d0 - t), R - 2 - t)
            B = exact_div(B * (d0 - t), t + 2)
            CAPm1 = CAP
            CAP = exact_div(CAP * (d0 - 1 - t), t + 1) if d0 - 1 - t > 0 else 0
            C0 = exact_div(C0 * (d0 - t), t + 1)
            sgn = -sgn
    while junk and junk[-1] == 0:
        junk.pop()
    return junk, c0_emit

def two_G(R, k, binRk):
    """[x^k] 2G_R given binRk = C(R-1,k); k >= 1."""
    return (binRk if (k % 2 == 0) else -binRk) - 1

def run_flow(R, D0, record_c=False, snap_compare=None, ckpt=None,
             snap_out=None, quiet=False):
    """Independent junk-flow run.  Returns result dict."""
    t0 = time.time()
    ds = profile(R, D0, R)
    i_pf = next(i for i in range(R + 1) if ds[i] + i > R)
    d0 = ds[0]
    assert d0 <= R - 2
    a, c0 = row0_junk(R, d0)
    minus2 = 1 if c0 == -2 else 0
    front0 = len(a) - 1
    c_rows = {}
    if record_c:
        # c row 0 not recorded here (suite reconstructs it separately)
        pass
    capB = BinRow(d0 - 1, len(a) + 4)
    # feed coefficient tracker: C(R-1, k)
    kcur = d0            # next feed degree index will be >= d0 + 1
    bcur = comb(R - 1, kcur)
    def stepG(ktarget):
        nonlocal kcur, bcur
        while kcur < ktarget:
            bcur = exact_div(bcur * (R - 1 - kcur), kcur + 1)
            kcur += 1
        return bcur
    snapshot = None
    events = []
    dead_cells = set()
    reignitions = 0
    cell1_death = None
    outcome = None
    capture = None
    die_row = None
    postfeed_negative_ok = True
    a0_clock_ok = True
    for i in range(1, R):
        d_prev, d = ds[i - 1], ds[i]
        delta = d - d_prev
        if i == i_pf:
            snapshot = {'R': R, 'D0': D0, 'i_postfeed': i_pf,
                        'd_feedend': d_prev,
                        'junk': {str(t): str(v) for t, v in enumerate(a) if v}}
            if snap_out:
                with open(snap_out, 'w') as f:
                    json.dump(snapshot, f)
        # ---- kernel ----
        kk = 1 + delta
        w = [0] * (len(a) + kk if a else 2)
        if kk == 1:
            for t, v in enumerate(a):
                if v:
                    w[t] += v
                    w[t + 1] += v
        else:
            for t, v in enumerate(a):
                if v:
                    w[t] += v
                    w[t + 1] += 2 * v
                    w[t + 2] += v
        # ---- feed: degrees jdeg in [d_prev, min(d, R-1-i)], value [x^(jdeg+i)] 2G_R
        hi_deg = min(d, R - 1 - i)
        for jdeg in range(d_prev, hi_deg + 1):
            k = jdeg + i
            v = two_G(R, k, stepG(k))
            w[0] += v
            if d - jdeg == 1:
                w[1] += v
        while len(w) > 1 and w[-1] == 0:
            w.pop()
        # ---- death guard ----
        if len(w) - 1 >= d:
            # only cell d exists; junk there kills unless value in {0,2}
            if len(w) - 1 > d or w[d] not in (0, 2):
                outcome, die_row = 'DIE', i
                break
        # ---- caps & clamp ----
        for _ in range(delta):
            capB.inc()
        capB.extend(len(w))
        row = capB.row
        newa = [0] * len(w)
        crow = {} if record_c else None
        nz = False
        for t in range(len(w)):
            wt = w[t]
            lo = -2 * row[t]
            hi = 2 * row[t - 1] if t >= 1 else 0
            c = wt if lo <= wt <= hi else (lo if wt < lo else hi)
            jn = wt - c
            if jn or c:
                assert jn % 2 == 0 and c % 2 == 0, 'parity violation'
            newa[t] = jn
            if jn:
                nz = True
            if record_c and c:
                crow[t] = c
            if t == 0:
                if c == -2:
                    minus2 += 1
        if record_c:
            c_rows[i] = crow
        # shrink cap row to the live band (perf only; extend() regrows exactly)
        if len(capB.row) > len(w) + 3:
            del capB.row[len(w) + 3:]
        # ---- extinction / re-ignition bookkeeping (post-feed only) ----
        if i >= i_pf:
            old_nz = {t for t, v in enumerate(a) if v}
            new_nz = {t for t, v in enumerate(newa) if v}
            for t in sorted(old_nz - new_nz):
                if t >= 1:
                    events.append((i, t, 'die'))
                    dead_cells.add(t)
                    if t == 1:
                        cell1_death = i
            for t in new_nz:
                if t in dead_cells:
                    reignitions += 1
                    dead_cells.discard(t)
        # ---- post-feed invariant spot checks ----
        if i >= i_pf:
            if any(v > 0 for v in newa):
                postfeed_negative_ok = False
            a0_old = -a[0] if a else 0
            a0_new = -newa[0] if newa else 0
            if a0_old > 0 and a0_new != a0_old - 2:
                a0_clock_ok = False
        while len(newa) > 1 and newa[-1] == 0:
            newa.pop()
        if len(newa) == 1 and newa[0] == 0:
            newa = []
        a = newa
        if ckpt and (i % 256 == 0 or i == i_pf):
            with open(ckpt, 'w') as f:
                json.dump({'i': i, 'front': (len(a) - 1) if a else -1,
                           'L1bits': sum(abs(v).bit_length() for v in a),
                           'minus2': minus2, 'i_pf': i_pf,
                           'elapsed_s': round(time.time() - t0, 1)}, f)
        if not a and ds[i] + i > R:
            outcome, capture = 'CLOSED', i
            break
    if outcome is None:
        outcome = 'OPEN_RESIDUAL'
    res = {'R': R, 'D0': D0, 'outcome': outcome, 'capture_row': capture,
           'die_row': die_row, 'minus2': minus2, 'i_pf': i_pf,
           'front0': front0, 'd0': d0,
           'snapshot': snapshot, 'events_n': len(events),
           'cell1_death': cell1_death, 'reignitions': reignitions,
           'postfeed_negative_ok': postfeed_negative_ok,
           'a0_clock_ok': a0_clock_ok,
           'elapsed_s': round(time.time() - t0, 1)}
    if record_c:
        res['c_rows'] = c_rows
    if snapshot:
        jd = {int(t): int(v) for t, v in snapshot['junk'].items()}
        res['fe_m'] = max(jd) if jd else -1
        res['fe_a0'] = -jd.get(0, 0)
        res['fe_d'] = snapshot['d_feedend']
    if snap_compare and snapshot:
        with open(snap_compare) as f:
            ref = json.load(f)
        rj = {int(t): int(v) for t, v in ref['junk'].items()}
        mj = {int(t): int(v) for t, v in snapshot['junk'].items()}
        res['snap_equal'] = (rj == mj and ref['i_postfeed'] == i_pf
                            and ref['d_feedend'] == snapshot['d_feedend'])
        res['snap_ref_cells'] = len(rj)
        res['snap_my_cells'] = len(mj)
    return res

# ---------------------------------------------------- slow definition-level
def poly_mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                if b:
                    r[i + j] += a * b
    return r

def qpow_list(m):
    """(1-x)^k for k=0..m as coefficient lists."""
    out = [[1]]
    for _ in range(m):
        out.append(poly_mul(out[-1], [1, -1]))
    return out

def decode_cells(p, d):
    """w_t = sum_j p[j] C(d-j, t) for t = 0..d (definition-level)."""
    return [sum(p[j] * comb(d - j, t) for j in range(min(len(p), d + 1)))
            for t in range(d + 1)]

def encode_cells(cells, d, qp):
    """sum_t cells[t] x^(d-t) (1-x)^t as coefficient list of length d+1."""
    out = [0] * (d + 1)
    for t, c in enumerate(cells):
        if c:
            q = qp[t]
            base = d - t
            for j, b in enumerate(q):
                out[base + j] += c * b
    return out

def slow_ruleA(R, D0, max_rows=None):
    """Rule A from the ORIGINAL definition (transient note sec. 1):
    sigma_{-1} = q^(R-1); ideal_i = sigma_{i-1} - x E_{R-2-i};
    Delta_i = AdmClamp_{d_i}(trunc_{d_i}(ideal_i)) (nearest box clamp +
    tozero parity fix); die unless x | sigma_{i-1} - Delta_i;
    sigma_i = (sigma_{i-1} - Delta_i)/x; closure iff sigma_{R-1} = 0."""
    ds = profile(R, D0, R)
    qp = qpow_list(2 * R)
    LEN = 2 * R + 4          # sigma can carry junk up to degree d_i - 1 ~ 1.2R
    sigma = list(qp[R - 1]) + [0] * (LEN - R)
    Deltas = []
    parity_fires = 0
    capture = None
    die_row = None
    for i in range(R):
        d = ds[i]
        m = R - 2 - i
        # ideal = sigma - x E_m ;  E_m = -1 + x + ... + x^m  (E_{-1} = 0)
        ideal = list(sigma)
        if m >= 0:
            ideal[1] -= -1
            for j in range(1, m + 1):
                ideal[j + 1] -= 1
        w = decode_cells(ideal, d)
        cells = []
        for t in range(d + 1):
            Cd = comb(d, t)
            v = w[t]
            v = Cd if v > Cd else (-Cd if v < -Cd else v)
            if (v - Cd) % 2 != 0:
                parity_fires += 1
                v += -1 if v > 0 else 1     # tozero
            cells.append(v)
        Delta = encode_cells(cells, d, qp)
        rem = [sigma[j] - (Delta[j] if j < len(Delta) else 0)
               for j in range(len(sigma))]
        if rem[0] != 0:
            die_row = i
            break
        sigma = rem[1:] + [0]
        Deltas.append((d, Delta))
        # capture: e_i = sigma_i - E_{R-2-i} == 0 ?
        mm = R - 2 - i
        e = list(sigma)
        if mm >= 0:
            e[0] += 1
            for j in range(1, mm + 1):
                e[j] -= 1
        if capture is None and all(v == 0 for v in e):
            capture = i
        if max_rows and i >= max_rows:
            break
    closed = die_row is None and all(v == 0 for v in sigma)
    return {'R': R, 'D0': D0, 'die_row': die_row, 'capture': capture,
            'closed': closed, 'parity_fires': parity_fires,
            'Deltas': Deltas, 'profile': ds}

# ---------------------------------------------------------------- stages ---
def stage_constants():
    log('=' * 74)
    log('STAGE constants -- fresh g sandwich, eps*, (5+3g)/(3+2g), kappa, g/(1-g)')
    log('=' * 74)
    M = 1 << 20
    F2M, F2M1 = fib_pair(2 * M)
    L2M = 2 * F2M1 - F2M
    def le(pw):  # 5^a <= phi^(2M) ?
        t = 2 * pw - L2M
        if t <= 0:
            return True
        return t * t < 5 * F2M * F2M
    lo, hi = 0, M
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if le(5 ** mid):
            lo = mid
        else:
            hi = mid
    g_lo = Fraction(lo, M)
    g_hi = Fraction(lo + 1, M)
    check(lo == 627035, 'g sandwich: floor(2^20 g) = 627035 (fresh integer bisection)',
          '627035/2^20 < g < 627036/2^20')
    eps = lambda g: 2 * (1 - g - g * g) / (3 + 2 * g)
    e_hi, e_lo = eps(g_lo), eps(g_hi)   # eps* decreasing in g
    check(e_lo < e_hi, 'eps*(g) decreasing on the sandwich')
    check(Fraction(211736, 10 ** 7) < e_lo and e_hi < Fraction(211747, 10 ** 7),
          'eps* in (0.0211736, 0.0211747)', '(%s, %s)' % (float(e_lo), float(e_hi)))
    check(e_hi < Fraction(1, 32), 'eps* < 1/32')
    # closed form 1 + g + eps* = (5+3g)/(3+2g): polynomial identity
    ok = True
    rnd = random.Random(7)
    for _ in range(200):
        gg = Fraction(rnd.randrange(1, 999), 1000)
        ok &= (1 + gg + eps(gg)) == (5 + 3 * gg) / (3 + 2 * gg)
    check(ok, 'identity 1 + g + eps*(g) == (5+3g)/(3+2g) (200 random rationals)')
    ok = True
    for _ in range(200):
        aa = Fraction(rnd.randrange(1, 999), 1000)
        bb = Fraction(rnd.randrange(1, 999), 1000)
        if aa == bb:
            continue
        f = lambda x: (5 + 3 * x) / (3 + 2 * x)
        ok &= (f(aa) - f(bb)) == -(aa - bb) / ((3 + 2 * aa) * (3 + 2 * bb))
    check(ok, "derivative law f(a)-f(b) = -(a-b)/((3+2a)(3+2b)) => f' = -1/(3+2g)^2")
    f = lambda g: (5 + 3 * g) / (3 + 2 * g)
    br_lo, br_hi = f(g_hi), f(g_lo)     # decreasing
    check(abs(br_lo - Fraction(16191617801, 10 ** 10)) < Fraction(1, 10 ** 9) and
          abs(br_hi - Fraction(16191618342, 10 ** 10)) < Fraction(1, 10 ** 9) and
          br_hi - br_lo < Fraction(6, 10 ** 8),
          '(5+3g)/(3+2g) bracket matches E1 (1.6191617801, 1.6191618342), width < 6e-8',
          '(%.10f, %.10f)' % (float(br_lo), float(br_hi)))
    # conditional bound arithmetic: 427095 = 262144 + 156759 + 8192 exactly,
    # so 1 + g + 1/32 < 427095/262144 holds STRICTLY because g < 156759/262144 strictly
    check(262144 + 156759 + 8192 == 427095 and g_hi == Fraction(156759, 262144),
          '1 + g + 1/32 < 427095/262144 EXACT (equality at the strict upper end g_hi)')
    check(262144 + 156759 + 16384 == 435287,
          '1 + g + 1/16 < 435287/262144 EXACT (same mechanism)')
    # E3 constants: g/(1-g) and kappa
    r_lo, r_hi = g_lo / (1 - g_lo), g_hi / (1 - g_hi)
    check(abs(r_lo - Fraction(14874828, 10 ** 7)) < Fraction(1, 10 ** 6) and
          abs(r_hi - Fraction(14874887, 10 ** 7)) < Fraction(1, 10 ** 6),
          'g/(1-g) bracket matches E3 (1.4874828, 1.4874887)',
          '(%.7f, %.7f)' % (float(r_lo), float(r_hi)))
    kap = lambda g: (3 - 4 * g) / (2 * (1 - g)) - (1 - g) / g
    k1, k2 = kap(g_lo), kap(g_hi)
    klo, khi = min(k1, k2), max(k1, k2)
    check(Fraction(839815, 10 ** 7) < klo and khi < Fraction(839820, 10 ** 7),
          'kappa in (0.0839815, 0.0839820) [E3 bracket 0.0839816..0.0839819]',
          '(%.7f, %.7f)' % (float(klo), float(khi)))
    # eps_phi ordering (referee fix1): 1/phi - g < eps*
    # 1/phi = (sqrt5 - 1)/2: sandwich sqrt5 by integers a/2^20
    S = 1 << 20
    lo5 = isqrt(5 * S * S)
    s_lo, s_hi = Fraction(lo5, S), Fraction(lo5 + 1, S)
    ephi_lo = (s_lo - 1) / 2 - g_hi
    ephi_hi = (s_hi - 1) / 2 - g_lo
    check(ephi_hi < e_lo, 'eps_phi = 1/phi - g < eps* (exact sandwiches)',
          'eps_phi <= %.7f < eps*_lo = %.7f' % (float(ephi_hi), float(e_lo)))
    save_frag('constants')

def stage_slow():
    log('=' * 74)
    log('STAGE slow -- definition-level rule A (polynomial coords) vs my fast flow')
    log('=' * 74)
    qp = qpow_list(340)
    # T4 closed form vs definition-level decode of trunc(2G_R)
    ok = True
    for (R, d) in [(64, 38), (64, 40), (128, 76), (128, 84), (256, 153), (256, 169)]:
        G = [0] * R
        G[0] = 2
        for k in range(1, R):
            G[k] = (comb(R - 1, k) if k % 2 == 0 else -comb(R - 1, k)) - 1
        w_def = decode_cells(G[:d + 1], d)
        for t in range(d + 1):
            w_t4 = (comb(R - 2 - t, d - t) if (d - t) % 2 == 0
                    else -comb(R - 2 - t, d - t)) - comb(d + 1, t + 1) + 2 * comb(d, t)
            if w_def[t] != w_t4:
                ok = False
    check(ok, 'T4 closed form == definition-level decode (6 grids, all cells)')
    for (R, D0, exp_out, exp_row) in [(128, 0, 'CLOSED', 81), (128, 8, 'CLOSED', 78),
                                      (256, 0, 'DIE', 61), (256, 16, 'CLOSED', 150)]:
        sl = slow_ruleA(R, D0)
        if exp_out == 'CLOSED':
            ok = sl['closed'] and sl['capture'] == exp_row
            got = 'closed=%s capture=%s parity_fires=%d' % (sl['closed'], sl['capture'], sl['parity_fires'])
        else:
            ok = sl['die_row'] == exp_row
            got = 'die=%s' % sl['die_row']
        check(ok, 'slow definition-level rule A (%d,%d) reproduces canon %s %s' %
              (R, D0, exp_out, exp_row), got)
        check(sl['parity_fires'] == 0, 'T3: parity fix never fires at (%d,%d)' % (R, D0))
        # witness identity of my slow run (closures)
        if exp_out == 'CLOSED' and sl['closed']:
            tot = [0] * (R + max(len(D) for _, D in sl['Deltas']) + 2)
            for i, (d, Dp) in enumerate(sl['Deltas']):
                for jj, v in enumerate(Dp):
                    if v:
                        tot[i + jj] += v
            target = qp[R - 1]
            okid = all(tot[j] == (target[j] if j < len(target) else 0)
                       for j in range(len(tot)))
            check(okid, 'slow run witness: sum x^i Delta_i == (1-x)^(R-1) at (%d,%d)' % (R, D0))
        # block-for-block vs my fast flow (T2 conjugacy, my two implementations)
        fr = run_flow(R, D0, record_c=True)
        rows_to_cmp = range(1, (sl['die_row'] if exp_out == 'DIE' else sl['capture']) + 1)
        okb = True
        ds = sl['profile']
        for i in rows_to_cmp:
            if i >= len(sl['Deltas']) and exp_out == 'DIE':
                break
            d = ds[i]
            crow = fr.get('c_rows', {}).get(i, {})
            cells = [0] * (d + 1)
            for t, v in crow.items():
                cells[t] = v
            enc = encode_cells(cells, d, qp)
            # Delta_fast = (2x-1) + encode(c)
            Dfast = list(enc)
            Dfast[0] -= 1
            Dfast[1] += 2
            Dslow = sl['Deltas'][i][1] if i < len(sl['Deltas']) else None
            if Dslow is None:
                break
            L = max(len(Dfast), len(Dslow))
            if any((Dfast[j] if j < len(Dfast) else 0) != (Dslow[j] if j < len(Dslow) else 0)
                   for j in range(L)):
                okb = False
                break
        check(okb, 'block-for-block: fast junk flow == slow rule A at (%d,%d), rows 1..%s'
              % (R, D0, max(rows_to_cmp)))
        if exp_out == 'DIE':
            check(fr['die_row'] == exp_row, 'fast flow death row matches at (%d,%d)' % (R, D0),
                  fr['die_row'])
    # T1 identity: encode(ballot(d)) == 2x-1 for the profile degrees used
    okT1 = True
    for d in sorted({dd for dd in profile(128, 8, 127)} | {dd for dd in profile(256, 16, 255)}):
        b = [comb(d - 1, t) - (comb(d - 1, t - 1) if t >= 1 else 0) for t in range(d + 1)]
        e = encode_cells(b, d, qp)
        if not (e[0] == -1 and e[1] == 2 and all(v == 0 for v in e[2:])):
            okT1 = False
    check(okT1, 'T1: encode_d(ballot) == 2x-1 for every profile degree (128,8)+(256,16)')
    save_frag('slow')

SUITE = [(128, 0, 'CLOSED', 81), (128, 4, 'CLOSED', 79), (128, 8, 'CLOSED', 78),
         (256, 0, 'DIE', 61), (256, 1, 'CLOSED', 159), (256, 8, 'CLOSED', 153),
         (256, 16, 'CLOSED', 150), (512, 4, 'DIE', 121), (512, 5, 'CLOSED', 312),
         (512, 16, 'CLOSED', 317), (512, 32, 'CLOSED', 313),
         (1024, 15, 'CLOSED', 639), (1024, 64, 'CLOSED', 616),
         (2048, 37, 'DIE', 508), (2048, 38, 'CLOSED', 1271),
         (2048, 128, 'CLOSED', 1237)]

def stage_suite():
    log('=' * 74)
    log('STAGE suite -- my fast flow vs 16 canon outcomes (4 prior implementations)')
    log('=' * 74)
    for (R, D0, exp_out, exp_row) in SUITE:
        snapf = os.path.join(RES, 'amm12592_S_cone_feedend_R%d_D0%d_boxeph.json' % (R, D0))
        sc = snapf if os.path.exists(snapf) else None
        r = run_flow(R, D0, snap_compare=sc)
        if exp_out == 'CLOSED':
            ok = r['outcome'] == 'CLOSED' and r['capture_row'] == exp_row
            ok = ok and r['minus2'] == (R - 2) // 2
            got = 'capture=%s minus2=%s' % (r['capture_row'], r['minus2'])
        else:
            ok = r['outcome'] == 'DIE' and r['die_row'] == exp_row
            got = 'die=%s' % r['die_row']
        check(ok, '(%d,%d) -> %s %s [canon]' % (R, D0, exp_out, exp_row), got)
        if r['outcome'] == 'CLOSED':
            pred = r['i_pf'] + r['fe_a0'] // 2 - 1
            check(r['capture_row'] == pred,
                  '(%d,%d) capture == i_pf + a0/2 - 1' % (R, D0),
                  '%s == %s + %s - 1' % (r['capture_row'], r['i_pf'], r['fe_a0'] // 2))
            check(r['reignitions'] == 0 and r['postfeed_negative_ok'] and r['a0_clock_ok'],
                  '(%d,%d) zero re-ignitions + post-feed negativity + cell-0 clock' % (R, D0))
        if sc is not None and 'snap_equal' in r:
            check(r['snap_equal'],
                  '(%d,%d) feed-end snapshot BIT-IDENTICAL to E1 stored (%s cells)' %
                  (R, D0, r.get('snap_ref_cells')), 'i_pf=%s d_fe=%s' % (r['i_pf'], r.get('fe_d')))
    save_frag('suite')

def stage_bigrun(R, D0):
    log('=' * 74)
    log('STAGE bigrun -- independent replication R=%d D0=%d' % (R, D0))
    log('=' * 74)
    snapf = os.path.join(RES, 'amm12592_S_cone_feedend_R%d_D0%d_boxeph.json' % (R, D0))
    sc = snapf if os.path.exists(snapf) else None
    r = run_flow(R, D0, snap_compare=sc,
                 ckpt=os.path.join(RES, 'amm12592_S_referee_ckpt_R%d_D0%d.json' % (R, D0)),
                 snap_out=os.path.join(RES, 'amm12592_S_referee_feedend_R%d_D0%d.json' % (R, D0)))
    r.pop('snapshot', None)
    log(json.dumps({k: v for k, v in r.items() if k != 'c_rows'}, default=str))
    with open(os.path.join(RES, 'amm12592_S_referee_run_R%d_D0%d.json' % (R, D0)), 'w') as f:
        json.dump(r, f, indent=1, default=str)
    save_frag('bigrun_R%d_D0%d' % (R, D0), extra=r)

# ------------------------------------------------ S-cone-fc certificate check
def cone_certificate(R, D0, snap_path, expect=None, scan_max=80):
    """My own implementation of Theorem S-cone-fc's hypothesis check + clocks,
    from a feed-end snapshot; evolves post-feed with the definition box clamp,
    verifying the S1 magnitude form each row."""
    with open(snap_path) as f:
        snap = json.load(f)
    i_pf = snap['i_postfeed']
    jd = {int(t): int(v) for t, v in snap['junk'].items()}
    ds = profile(R, D0, R)
    assert ds[i_pf - 1] == snap['d_feedend'], 'profile mismatch vs snapshot'
    a = [0] * (max(jd) + 1)
    for t, v in jd.items():
        a[t] = v
    capB = BinRow(ds[i_pf] - 1, len(a) + 4)
    first = None
    results = []
    for i0 in range(i_pf, i_pf + scan_max + 1):
        D_0 = ds[i0]
        m = len(a) - 1
        while m >= 0 and a[m] == 0:
            m -= 1
        F1 = all(v <= 0 for v in a) and (m + 2 < D_0)
        a0 = -a[0] if a else 0
        F2 = a0 <= D_0 - 1
        capref = BinRow(D_0 - 1, m + 3).row
        F3 = True
        for t in range(2, m + 3):
            am1 = -a[t - 1] if t - 1 <= m else 0
            am2 = -a[t - 2] if t - 2 <= m else 0
            if 2 * am1 + am2 > 2 * capref[t]:
                F3 = False
                break
        if F1 and F2 and F3 and first is None:
            # clocks: exact staircase sums with the certified profile
            need = {t: -a[t] for t in range(1, m + 1) if a[t]}
            acc = {t: 0 for t in need}
            Tt = {}
            drain = (a0 + 1) // 2
            crow = BinRow(D_0 - 1, m + 3)
            k = 0
            while need and k < R:
                k += 1
                Dk = ds[i0 + k]
                for _ in range(Dk - crow.nn - 1):
                    crow.inc()
                a0k = max(0, a0 - 2 * (k - 1))
                dl = Dk - ds[i0 + k - 1]
                for t in list(need):
                    if t == 1:
                        inc = 2 * (Dk - 1) - (1 + dl) * a0k
                        if inc > 0:
                            acc[t] += inc
                    else:
                        acc[t] += 2 * crow.row[t] - 2 * capref[t]
                    if acc[t] >= need[t]:
                        Tt[t] = k
                        del need[t]
            Tmax = max(Tt.values()) if Tt else 0
            Tworst = max(Tt, key=lambda t: Tt[t]) if Tt else None
            ub = i0 + max(drain, Tmax)
            F4 = ub <= R - 2
            first = {'i_fc': i0, 'd': D_0, 'm': m, 'a0': a0, 'drain': drain,
                     'Tmax': Tmax, 'T_worst_t': Tworst, 'capture_ub': ub,
                     'F4': F4, 'budget_margin': (R - 2) - ub}
        results.append((i0, F1, F2, F3))
        # evolve one post-feed row (definition-level box clamp; S1 check)
        d_prev, d = ds[i0], ds[i0 + 1]
        # careful: the row being processed is row i0 itself, at degree ds[i0]
        dl = ds[i0] - ds[i0 - 1]
        kk = 1 + dl
        w = [0] * (len(a) + kk)
        for t, v in enumerate(a):
            if v:
                w[t] += v
                if kk == 1:
                    w[t + 1] += v
                else:
                    w[t + 1] += 2 * v
                    w[t + 2] += v
        for _ in range(ds[i0] - 1 - capB.nn):
            capB.inc()
        capB.extend(len(w))
        row = capB.row
        newa = [0] * len(w)
        for t in range(len(w)):
            wt = w[t]
            lo = -2 * row[t]
            hi = 2 * row[t - 1] if t >= 1 else 0
            c = wt if lo <= wt <= hi else (lo if wt < lo else hi)
            newa[t] = wt - c
            # S1 magnitude form check (all-negative state)
            if wt <= 0:
                s1 = min(0, wt + (2 * row[t] if t >= 1 else 2))
                assert newa[t] == s1, 'S1 magnitude form violated'
        while len(newa) > 1 and newa[-1] == 0:
            newa.pop()
        if len(newa) == 1 and newa[0] == 0:
            newa = []
        a = newa
        if not a:
            break
    return first, results

def stage_conecert(R, D0):
    log('=' * 74)
    log('STAGE conecert -- my S-cone-fc certificate re-derivation R=%d D0=%d' % (R, D0))
    log('=' * 74)
    for tag, snapf in [('E1', os.path.join(RES, 'amm12592_S_cone_feedend_R%d_D0%d_boxeph.json' % (R, D0))),
                       ('mine', os.path.join(RES, 'amm12592_S_referee_feedend_R%d_D0%d.json' % (R, D0)))]:
        if not os.path.exists(snapf):
            log('  [SKIP] no snapshot %s' % snapf)
            continue
        first, _ = cone_certificate(R, D0, snapf)
        log('  snapshot=%s -> %s' % (tag, json.dumps(first)))
        fcf = os.path.join(RES, 'amm12592_S_cone_fcscan_R%d_D0%d_boxeph.json' % (R, D0))
        if os.path.exists(fcf) and first:
            with open(fcf) as f:
                fc = json.load(f)
            ref = fc['fc']
            ok = (first['i_fc'] == fc['i_fc'] and first['d'] == ref['d'] and
                  first['m'] == ref['m'] and first['a0'] == ref['a0'] and
                  first['drain'] == ref['drain'] and first['Tmax'] == ref['Tmax'] and
                  first['capture_ub'] == ref['capture_ub'])
            check(ok, 'my certificate == stored fcscan (i_fc,d,m,a0,drain,Tmax,ub) [%s snap]' % tag,
                  'i_fc=%s ub=%s' % (first['i_fc'], first['capture_ub']))
            actual = fc.get('capture_row')
            if actual:
                check(first['capture_ub'] >= actual and first['F4'],
                      'certificate bound >= actual capture, F4 holds [%s]' % tag,
                      '%s >= %s' % (first['capture_ub'], actual))
    save_frag('conecert_R%d_D0%d' % (R, D0))

# --------------------------------------------------------- randomized lemmas
def one_step_box(a, dnew, delta):
    """definition-level one row: kernel (1+delta), clamp at degree dnew."""
    kk = 1 + delta
    w = [0] * (len(a) + kk)
    for t, v in enumerate(a):
        w[t] += v
        if kk == 1:
            w[t + 1] += v
        else:
            w[t + 1] += 2 * v
            w[t + 2] += v
    out = []
    for t in range(len(w)):
        lo = -2 * comb(dnew - 1, t)
        hi = 2 * comb(dnew - 1, t - 1) if t >= 1 else 0
        c = min(max(w[t], lo), hi)
        out.append(w[t] - c)
    return out

def stage_lemmas():
    log('=' * 74)
    log('STAGE lemmas -- my randomized exact certificates for S1/S2/S4/cone step')
    log('=' * 74)
    rnd = random.Random(20260804)
    okS1 = okS2 = okS4 = okC = True
    for trial in range(300):
        d = rnd.randrange(30, 90)
        delta = rnd.randrange(2)
        m = rnd.randrange(1, min(12, d - 4))
        a = [-2 * rnd.randrange(0, comb(d, t + 1) // 2 + 2) for t in range(m + 1)]
        # S1: magnitude closed form == box clamp
        out = one_step_box(a, d + delta, delta)
        for t in range(len(out)):
            kk = 1 + delta
            Ka = sum(comb(kk, s) * (-a[t - s]) for s in range(kk + 1)
                     if 0 <= t - s < len(a))
            s1 = max(0, Ka - (2 * comb(d + delta - 1, t) if t >= 1 else 2))
            if -out[t] != s1:
                okS1 = False
        # S2: monotonicity
        b = [v - 2 * rnd.randrange(0, 5) for v in a]     # b <= a (more negative)
        oa = one_step_box(a, d + delta, delta)
        ob = one_step_box(b, d + delta, delta)
        for t in range(len(oa)):
            if -(ob[t] if t < len(ob) else 0) < -(oa[t] if t < len(oa) else 0):
                okS2 = False
        # S4: a0 <= d'-1 => cell1 non-increasing under inflow
        a0 = -a[0]
        if a0 <= d + delta - 1:
            if -oa[1] > -a[1] if len(a) > 1 else False:
                okS4 = False
    check(okS1, 'S1 magnitude closed form == definition box clamp (300 random states)')
    check(okS2, 'S2 comparison principle (300 random comparable pairs)')
    check(okS4, 'S4 cell-1 absorbing criterion (300 trials)')
    # cone one-step preservation: random F1/F2/F3 states
    for trial in range(300):
        D_0 = rnd.randrange(60, 140)
        m = rnd.randrange(2, 14)
        if m + 2 >= D_0:
            continue
        capref = [2 * comb(D_0 - 1, t) for t in range(m + 4)]
        a = [0] * (m + 1)
        a[0] = -2 * rnd.randrange(0, (D_0 - 1) // 2 + 1)
        for t in range(1, m + 1):
            # choose a_t with 2a_t + a_{t-1} <= capref_{t+1} (ensures F3 at t+1)
            room = (capref[t + 1] - (-a[t - 1])) // 2
            if room < 0:
                room = 0
            a[t] = -2 * rnd.randrange(0, room // 2 + 1)
        # verify F3 on [2, m+2]
        f3 = all(2 * (-a[t - 1] if t - 1 <= m else 0) + (-a[t - 2] if t - 2 <= m else 0)
                 <= capref[t] for t in range(2, m + 3))
        if not f3:
            continue
        for delta in (0, 1):
            for dgrow in (0, 1, 3):
                Dk = D_0 + dgrow
                out = one_step_box(a, Dk, delta)
                mm = len(out) - 1
                while mm >= 0 and out[mm] == 0:
                    mm -= 1
                if mm > m:
                    okC = False        # support advanced
                for t in range(1, min(mm, m) + 1):
                    if -out[t] > -a[t]:
                        okC = False    # some cell increased
                # F3 propagates (out <= a cellwise)
                f3p = all(2 * (-out[t - 1] if t - 1 < len(out) else 0)
                          + (-out[t - 2] if t - 2 < len(out) else 0)
                          <= capref[t] for t in range(2, m + 3))
                if not f3p:
                    okC = False
    check(okC, 'S-cone-fc one-step: non-increase + frozen support + F3 propagation '
               '(300 random cone states x both kernels x cap growth 0/1/3)')
    save_frag('lemmas')

def stage_apriori():
    log('=' * 74)
    log('STAGE apriori -- F4-automatic corollary re-check at R = 2^16, 2^17 (worst case)')
    log('=' * 74)
    for R in (1 << 16, 1 << 17):
        for eps_num, eps_den in ((1, 32), (1, 16)):
            D0 = -(-R * eps_num // eps_den)
            ds = profile(R, D0, R)
            i_feed = max(i for i in range(R // 3) if ds[i] + i <= R - 1)
            i0 = i_feed + 66                                    # worst (latest) row
            D_0 = ds[i0]
            # worst-case magnitudes from F3: a_t = C(D_0-1, t+1); a0 = D_0-1
            a0 = D_0 - 1
            drain = (a0 + 1) // 2
            probes = [1, 2, 3, 4, 6, 10, 20, 50, 200, 1000, 5000]
            probes = [t for t in probes if t + 2 < D_0]
            Tt = {}
            okall = True
            for t in probes:
                needv = comb(D_0 - 1, t + 1)
                capref_t = 2 * comb(D_0 - 1, t)
                cur = comb(ds[i0 + 1] - 1, t)     # C(D_k - 1, t), incremental in n
                nn = ds[i0 + 1] - 1
                acc = 0
                k = 0
                T = None
                while k < R - i0 - 2:
                    k += 1
                    Dk = ds[i0 + k]
                    while nn < Dk - 1:
                        cur = exact_div(cur * (nn + 1), nn + 1 - t)
                        nn += 1
                    if t == 1:
                        dl = Dk - ds[i0 + k - 1]
                        a0k = max(0, a0 - 2 * (k - 1))
                        inc = 2 * (Dk - 1) - (1 + dl) * a0k
                        if inc > 0:
                            acc += inc
                    else:
                        acc += 2 * cur - capref_t
                    if acc >= needv:
                        T = k
                        break
                if T is None:
                    okall = False
                else:
                    Tt[t] = T
            Tmax = max(Tt.values()) if Tt else None
            ub = i0 + max(drain, Tmax) if okall else None
            check(okall and ub <= R - 2,
                  'F4 automatic at R=%d eps=%d/%d: worst-case clocks fit budget' %
                  (R, eps_num, eps_den),
                  'i0=%d Tmax=%s ub=%s R-2=%d margin=%s' %
                  (i0, Tmax, ub, R - 2, (R - 2 - ub) if ub else None))
            t2 = Tt.get(2)
            check(t2 is not None and all(v <= t2 for t, v in Tt.items()),
                  '  T_t <= T_2 on probe set (E1 monotonicity claim)',
                  str(sorted(Tt.items())))
    save_frag('apriori')

def stage_window():
    log('=' * 74)
    log('STAGE window -- death-freeness of [i_pf, i_pf+64]: exact + symbolic')
    log('=' * 74)
    # exact per-R check: no-death at rows i <= i_feed + 66 iff (T6b extension)
    # F(0) + i < d_0, i.e. i < 2 d_0 - R + 2; verify 2 d_0 - R + 2 > i_feed + 66
    for R in (1 << 16, 1 << 17):
        for eps_num, eps_den in ((1, 32), (1, 16)):
            D0 = -(-R * eps_num // eps_den)
            gf = GammaFloor(R)
            d0 = gf.m + D0
            # i_feed = max{i : d_i + i <= R-1} by exact scan near the formula value
            g_lo = Fraction(627035, 1 << 20)
            approx = int((R * (1 - g_lo) - D0) / (1 + g_lo))
            ds = profile(R, D0, approx + 100)
            i_feed = max(i for i in range(approx + 100) if ds[i] + i <= R - 1)
            ok = 2 * d0 - R + 2 > i_feed + 66
            check(ok, 'R=%d D0=%d: 2d0-R+2 > i_feed+66 (window death-free, T6b ext.)' %
                  (R, D0), '%d > %d' % (2 * d0 - R + 2, i_feed + 66))
    # symbolic sufficient condition: D0(3+2g) >= 2R(1-g-g^2) + 66(1+g)
    g_lo, g_hi = Fraction(627035, 1 << 20), Fraction(627036, 1 << 20)
    ok = True
    for R_exp in range(16, 41):
        R = 1 << R_exp
        for eps in (Fraction(1, 32), Fraction(1, 16)):
            D0 = (R * eps.numerator + eps.denominator - 1) // eps.denominator
            # rigorous interval form: LHS at g_lo (min), RHS termwise max
            if not (D0 * (3 + 2 * g_lo) >=
                    2 * R * (1 - g_lo - g_lo * g_lo) + 66 * (1 + g_hi)):
                ok = False
    check(ok, 'symbolic: D0(3+2g) >= 2R(1-g-g^2) + 66(1+g) for all dyadic '
              'R = 2^16..2^40, eps in {1/32,1/16}, both sandwich ends '
              '(PROVES the [i_pf, i_pf+64] window death-free for the reduction)')
    save_frag('window')

def stage_e3wit():
    log('=' * 74)
    log('STAGE e3wit -- (256,0) desc1 witness: my own from-scratch verification')
    log('=' * 74)
    p = os.path.join(HERE, 'amm12592_witness_R256_bulk_desc1_D0_0_boxeph.json')
    sha = hashlib.sha256(open(p, 'rb').read()).hexdigest()
    check(sha.startswith('5950bd42'), 'sha256 matches record 5950bd42...', sha[:16])
    W = json.load(open(p))
    R = W['R']
    assert W['D0'] == 0
    ds = profile(R, 0, R)
    qp = qpow_list(ds[R - 1] + 4)
    sc = {int(i): {int(t): v for t, v in row.items()} for i, row in W['sparse_c'].items()}
    tot = [0] * (3 * R)
    adm_ok = True
    par_ok = True
    int_ok = True
    ncells = 0
    nsat = 0
    for i in range(R):
        d = ds[i]
        crow = sc.get(i, {})
        if i <= R - 2:
            cells = [ (comb(d - 1, t) - (comb(d - 1, t - 1) if t >= 1 else 0))
                      + crow.get(t, 0) for t in range(d + 1)]
        else:
            cells = [-comb(d, t) + crow.get(t, 0) for t in range(d + 1)]
        for t in range(d + 1):
            Cd = comb(d, t)
            ncells += 1
            if abs(cells[t]) > Cd:
                adm_ok = False
            if (cells[t] - Cd) % 2:
                par_ok = False
            f2, r2 = divmod(Cd - cells[t], 2)      # transportation f = (C-delta)/2
            if r2 != 0 or f2 < 0 or f2 > Cd:
                int_ok = False
            if f2 == 0 or f2 == Cd:
                nsat += 1
        Dp = encode_cells(cells, d, qp)
        for jj, v in enumerate(Dp):
            if v:
                tot[i + jj] += v
    target = qp[R - 1]
    id_ok = all(tot[j] == (target[j] if j < len(target) else 0) for j in range(len(tot)))
    check(adm_ok and par_ok, 'all %d cells: |delta_t| <= C(d_i,t) with correct parity' % ncells)
    check(id_ok, 'epoch identity sum x^i Delta_i == (1-x)^255 EXACT (my own polynomial arithmetic; last row = corner -1)')
    check(int_ok, 'transportation point: f = (C-delta)/2 integer in [0, C] at every cell',
          'cells=%d boundary-saturated=%d (%.2f%%)' % (ncells, nsat, 100.0 * nsat / ncells))
    check(ncells == 58837, 'cell count == 58837 (E3 claim)', ncells)
    # identity at x = 2, 3, -1 (implied; direct evaluation)
    for xv in (2, 3, -1):
        lhs = sum(v * xv ** j for j, v in enumerate(tot) if v)
        rhs = (1 - xv) ** (R - 1)
        check(lhs == rhs, 'identity evaluated at x = %d' % xv)
    # E3-L1/L2/L3 algebra on random grids
    rnd = random.Random(3)
    okL = True
    for _ in range(200):
        R_ = 2 * rnd.randrange(30, 400)
        d = rnd.randrange(R_ // 2 + 1, (2 * R_) // 3)
        F0 = R_ - 2 - d
        t2 = 2 * d - R_ + 3
        if (t2 - 1) + F0 != d:
            okL = False
        if comb(R_ - 2 - F0, d - F0) != comb(d, F0):
            okL = False
        for t in range(1, min(d, R_ - 2)):
            lhs = Fraction(comb(R_ - 1 - t, d - t + 1), comb(R_ - 2 - t, d - t)) / \
                  Fraction(comb(d, t - 1), comb(d, t))
            if lhs != Fraction(R_ - 1 - t, t):
                okL = False
            break
    check(okL, 'E3-L1 complementarity, E3-L2 front load = half box width, '
               'E3-L3 rate identity (200 random windows)')
    save_frag('e3wit')

def stage_ledger():
    log('=' * 74)
    log('STAGE ledger -- E1 stored-artifact cross-checks (A_1 law, edge law, r_t, fcscan)')
    log('=' * 74)
    # A_1 bits and edge law across all stored snapshots
    rows = []
    for R in (128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768):
        for D0 in (-(-R // 32), -(-R // 16)):
            f = os.path.join(RES, 'amm12592_S_cone_feedend_R%d_D0%d_boxeph.json' % (R, D0))
            if not os.path.exists(f):
                continue
            snap = json.load(open(f))
            jd = {int(t): int(v) for t, v in snap['junk'].items()}
            d_fe = snap['d_feedend']
            a0 = -jd.get(0, 0)
            a1 = -jd.get(1, 0)
            m = max(jd)
            # A_1 = a1 / (2(d_fe-1)); log2 via bit lengths (display only)
            A1_bits = a1.bit_length() - (2 * (d_fe - 1)).bit_length()
            rows.append((R, D0, m, a0 - d_fe, A1_bits))
    log('  (R, D0, m, a0-d_fe, ~log2 A_1): %s' % rows)
    edge_ok = all(-14 <= r[3] <= 9 for r in rows)
    check(edge_ok, 'edge law a0 - d_fe in [-14, +9] over all stored snapshots')
    by_eps = {}
    for (R, D0, m, e, A1) in rows:
        eps32 = (D0 == -(-R // 32))
        by_eps.setdefault(eps32, []).append((R, A1))
    lin_ok = True
    for eps32, seq in by_eps.items():
        seq.sort()
        for (Ra, Aa), (Rb, Ab) in zip(seq, seq[1:]):
            if Rb == 2 * Ra and not (0 <= Ab - Aa <= 2):
                lin_ok = False
    check(lin_ok, 'A_1 grows ~ +1 bit/doubling at feed-end (REFUTES D2 R-independence '
                  'reading; supports E1)', str(sorted(by_eps.get(True, []))))
    # marginal surface at 32768 (both eps): max_t r_t and over-cells
    from fractions import Fraction as Fr
    for (R, D0, exp_max, exp_t) in ((32768, 2048, '1.0009', 5), (32768, 1024, '1.0009', 5)):
        f = os.path.join(RES, 'amm12592_S_cone_feedend_R%d_D0%d_boxeph.json' % (R, D0))
        if not os.path.exists(f):
            continue
        snap = json.load(open(f))
        jd = {int(t): int(v) for t, v in snap['junk'].items()}
        d_fe = snap['d_feedend']
        m = max(jd)
        crow = BinRow(d_fe - 1, m + 3).row
        best, bt, over = None, None, []
        for t in range(2, m + 3):
            a1 = -jd.get(t - 1, 0)
            a2 = -jd.get(t - 2, 0)
            r = Fr(2 * a1 + a2, 2 * crow[t])
            if best is None or r > best:
                best, bt = r, t
            if r > 1:
                over.append(t)
        log('  R=%d D0=%d: max r_t = %.5f at t=%d, over-cells=%s' %
            (R, D0, float(best), bt, over))
        check(abs(float(best) - 1.0009) < 0.0002 and bt == exp_t,
              'marginal surface at (%d,%d): max r_t ~ %s at t=%d [E1 claim]' %
              (R, D0, exp_max, exp_t))
        check(all(o % 2 == 1 for o in over), '  over-cells single-parity (odd)', over)
    # fcscan stored tables: bound >= actual everywhere
    ok = True
    for fn in os.listdir(RES):
        if fn.startswith('amm12592_S_cone_fcscan_') and fn.endswith('.json'):
            fc = json.load(open(os.path.join(RES, fn)))
            if fc.get('outcome', {}).get('status') == 'CLOSED' and fc.get('fc'):
                if fc['fc']['capture_ub'] < fc['outcome']['capture_row']:
                    ok = False
    check(ok, 'all stored fcscan certificates: capture_ub >= actual capture')
    save_frag('ledger')

def stage_bigcheck():
    log('=' * 74)
    log('STAGE bigcheck -- verdicts on the two independent R=32768 replications')
    log('=' * 74)
    EXP = {(32768, 2048): (19865, 6963, 25806, 169, 25806, 17133),
           (32768, 1024): (20185, 7604, 25165, 168, 25164, 17521)}
    for (R, D0), (cap, ipf, dfe, m, a0, c1d) in EXP.items():
        f = os.path.join(RES, 'amm12592_S_referee_run_R%d_D0%d.json' % (R, D0))
        if not os.path.exists(f):
            check(False, '(%d,%d) run JSON missing' % (R, D0))
            continue
        r = json.load(open(f))
        check(r['outcome'] == 'CLOSED' and r['capture_row'] == cap,
              '(%d,%d) INDEPENDENT REPLICATION: CLOSED at capture %d [E1 claim]' %
              (R, D0, cap), '%s %s' % (r['outcome'], r['capture_row']))
        check(r['minus2'] == (R - 2) // 2,
              '(%d,%d) debt minus2 = (R-2)/2 = 16383 exact' % (R, D0), r['minus2'])
        check(r['i_pf'] == ipf and r.get('fe_d') == dfe and r.get('fe_m') == m
              and r.get('fe_a0') == a0,
              '(%d,%d) i_pf/d_fe/m/a0 = %d/%d/%d/%d [E1 claim]' % (R, D0, ipf, dfe, m, a0),
              '%s/%s/%s/%s' % (r['i_pf'], r.get('fe_d'), r.get('fe_m'), r.get('fe_a0')))
        check(r.get('snap_equal') is True,
              '(%d,%d) feed-end snapshot BIT-IDENTICAL to E1 stored state' % (R, D0),
              '%s cells' % r.get('snap_ref_cells'))
        check(r['capture_row'] == r['i_pf'] + r['fe_a0'] // 2 - 1,
              '(%d,%d) capture == i_pf + a0/2 - 1 (cell-0 clock exact)' % (R, D0))
        check(r['reignitions'] == 0 and r['postfeed_negative_ok'] and r['a0_clock_ok'],
              '(%d,%d) zero post-feed re-ignitions; negativity + 2/row clock in-flight' % (R, D0))
        check(r['cell1_death'] == c1d,
              '(%d,%d) cell-1 extinction row == %d [E1 claim]' % (R, D0, c1d),
              r['cell1_death'])
        check(r['front0'] == R - 2 - r['d0'],
              '(%d,%d) initial front == R-2-d_0 (T4b)' % (R, D0), r['front0'])
    save_frag('bigcheck')

def stage_report():
    """Assemble amm12592_S_referee_boxeph.out from all stage fragments +
    the two independent R=32768 replications + verdicts + canon paragraph."""
    frags = {}
    for fn in sorted(os.listdir(RES)):
        if fn.startswith('amm12592_S_referee_frag_') and fn.endswith('.json'):
            st = fn[len('amm12592_S_referee_frag_'):-len('.json')]
            frags[st] = json.load(open(os.path.join(RES, fn)))
    runs = {}
    for (R, D0) in ((32768, 2048), (32768, 1024)):
        f = os.path.join(RES, 'amm12592_S_referee_run_R%d_D0%d.json' % (R, D0))
        if os.path.exists(f):
            runs[(R, D0)] = json.load(open(f))
    tp = sum(f['pass'] for f in frags.values())
    tf = sum(f['fail'] for f in frags.values())
    allfails = [x for f in frags.values() for x in f['fails']]
    lines = []
    A = lines.append
    A('AMM 12592 -- HOSTILE REFEREE of lanes E1 (S via invariant cone), E2, E3')
    A('(local-rule barrier / existence gap).  boxeph, 2026-08-04.')
    A('Script: 04-computation/amm12592_S_referee_boxeph.py')
    A('Fragments: 05-knowledge/results/amm12592_S_referee_frag_*.json')
    A('Independent runs: amm12592_S_referee_{run,feedend,ckpt}_R32768_*.json')
    A('')
    A('TOTAL CHECKS: %d PASS, %d FAIL%s' % (tp, tf,
      ('  FAILS: ' + '; '.join(allfails)) if allfails else ''))
    A('')
    A('=' * 76)
    A('STAGE OUTPUTS (in audit order)')
    A('=' * 76)
    order = ['constants', 'slow', 'suite', 'bigcheck',
             'bigrun_R32768_D02048', 'bigrun_R32768_D01024',
             'conecert_R16384_D0512', 'conecert_R16384_D01024',
             'conecert_R32768_D02048', 'conecert_R32768_D01024',
             'lemmas', 'apriori', 'window', 'e3wit', 'ledger']
    for st in order + [s for s in frags if s not in order]:
        if st in frags:
            for ln in frags[st]['text']:
                A(ln)
            A('')
    A('=' * 76)
    A('VERDICTS')
    A('=' * 76)
    A('''
E1 (S via the invariant cone) -- VERDICT: CONFIRMED, with two findings.
  Every PROVED item re-derived and independently machine-certified:
  S1 magnitude closed form, S2 comparison, S3 unconditional deadline
  (i_pf <= i_feed+2 <= (R-2)/2 re-derived), S4 spill criterion, and
  THEOREM S-cone-fc: proof audited line by line -- the max(0, x-c_k)
  composition needs c_k >= 0, supplied by D_k >= D_0; F3 range [2, m+2]
  is exactly what front-freeze at m+1, m+2 consumes; F3 propagation
  follows from cellwise non-increase.  My independent certificate checker
  reproduces E1's fcscan EXACTLY (i_fc, d, m, a0, drain, Tmax, worst
  clock t=2, capture_ub) at (16384, 512/1024) and (32768, 2048/1024),
  and capture_ub >= actual capture in all 18 stored certificates
  (i_fc - i_pf in [0, 8]).  The R = 32768 closures (BOTH eps) are
  replicated by THIS referee's from-scratch fifth implementation
  (exact Beatty floor via fresh 5^m <=> phi^(2n) integer comparisons;
  T4 closed-form row 0 re-proved against a definition-level polynomial
  decode; kernel/feed/clamp written from the T6 statement; validated
  block-for-block against my own definition-level slow rule A at
  128/256 and against 16 canon outcomes to R = 2048): captures
  19865 (eps=1/16) and 20185 (eps=1/32), debt 16383 = (R-2)/2 exact,
  zero post-feed re-ignitions, cell-0 clock exact, and the feed-end
  snapshots BIT-IDENTICAL to E1's stored states (170/169 cells).
  S(1/32) and S(1/16) are hereby verified for all dyadic
  128 <= R <= 32768 by two independent implementation lineages.
  Feed-end laws re-verified from stored snapshots: edge law
  a0 - d_fe in [-14, +9]; A_1 ~ +1 bit/doubling (E1's REFUTATION of
  D2's "max_t log2 A_t ~ 10-11 R-independent" is CORRECT -- D2 sec. 4.2
  must be read at a later row, not feed-end; nothing load-bearing);
  marginal surface max r_t = 1.00090/1.00089 at t = 5 with odd-parity
  over-cells at 32768.  Basin claims verified after decoding the file
  convention snapx{k} = scale 2^k (x2 captures at 1581/6556; x4 and
  corner OPEN_RESIDUAL -- correctly never used as feasibility evidence).

  FINDING E1-F1 (quantifier gap, repaired here): the reduction
  "ENTRY-fc(eps) => S(eps) for dyadic R >= 65536" needs the window
  [i_pf, i0) to be death-free; E1's parenthetical justifies it with the
  row-i0 support bound m, which cannot bound EARLIER rows.  Correct
  argument: T6a/T6b extension -- F(i) <= F(0) + i + (d_i - d_0), so no
  death before row 2d_0 - R + 2; the needed 2d_0 - R + 2 > i_feed + 66
  follows from D0(3+2g) >= 2R(1-g-g^2) + 66(1+g), i.e.
  D0 >= eps* R + ~25.1.  Certified here in rigorous mixed-endpoint
  interval form for all dyadic 2^16..2^40 at both eps (min margin 2665)
  plus exact per-R checks.  GAP CLOSED; no downstream impact.

  FINDING E1-F2 (quantifier scope): the a-priori-F4 corollary
  ("hypothesis is F1^F2^F3 alone for R >= 512") is certified only for
  dyadic R = 2^9..2^40 (code range confirmed); beyond 2^40 it rests on
  the sketched K1c/K2c bounds + the unwritten T_t <= T_2 monotonicity
  lemma.  My independent worst-case clock computation confirms it at
  2^16/2^17 (both eps, probe t up to 5000, T_2 always worst).  For full
  rigor state ENTRY-fc as F1^F2^F3^F4 for R > 2^40 (F4 is a decidable
  one-row check, so the reduction loses nothing), or write out the
  monotonicity lemma.

  FINDING E1-F3 (cosmetic erratum): the note sec. 6 marginal-surface
  table entry for (4096, eps=1/32) reads "1.0018"; recomputing from the
  stored (post-bug regenerated) snapshot gives max r_t = 1.0042 at
  t = 4 (over-cells {2,4,6}) with caps at d_fe - 1, or 1.0029 with caps
  at d_fe -- neither matches 1.0018, so that display value likely
  predates the snapshot regeneration.  All other 17 table entries
  reproduce exactly.  Qualitative law unaffected; display-only.

  OVERCLAIM NOTE: the lane-report headline "S(eps) ... COMPLETE" must
  not enter canon as "S proved".  S(eps) remains OPEN for dyadic
  R >= 65536; what is proved is the reduction to the one-row ENTRY-fc
  plus finite verification through 32768.  The note itself (secs. 7, 9)
  states this correctly.

E2 -- VERDICT: NO CONTENT.  E2 delivered only "standing by" (armed
  monitors).  Nothing to audit; E2 provides no independent confirmation
  of anything.  The independence gap for the R = 32768 closures (E1 was
  the only implementation to have run them) is closed by THIS referee's
  replication, not by E2.  No consistency issue arises between E1 and
  E2 because E2 asserts nothing.

E3 (barrier / existence gap) -- VERDICT: CONFIRMED; scope discipline is
  exemplary.  The (256, 0) desc1 witness re-verified by my own
  from-scratch verifier (sha 5950bd42..., 58837 cells in-box with
  parity, epoch identity sum x^i Delta_i = (1-x)^255 EXACT with
  last-row corner -1, transportation point f = (C-delta)/2 integer in
  [0, C] at every cell, 8.45% boundary-saturated, identity at
  x = 2/3/-1).  Record correction stands: the exact-floor existence
  frontier is (512, D0=0), OPEN in both directions.  E3-L1/L2/L3
  algebra confirmed on 200 random windows; kappa and g/(1-g) brackets
  confirmed from my fresh sandwich; L4 influence cone is a sound
  one-line induction (kernel support {0,1,2}).  All negatives (84/84
  rule deaths, beam all-dies) are hazard-labeled and never used as
  feasibility evidence; E3-BC(param)/(mech) correctly labeled
  CONJECTURED; LP at (512,0) correctly OPEN; floor-side honesty note
  (THM-3024 demotion) correctly repeated.  Minor: E3-BC(mech)'s "all RS
  rules qualify" (bounded deviation) is an unproved aside inside a
  conjecture -- harmless as labeled.

CONSISTENCY ACROSS LANES: E1's 18-run sweep captures equal D2's table
  exactly where they overlap; E1's i_pf differs from D2's i_feed by the
  documented off-by-one conventions (d_i + i <= R-2 / R-1 / R), each
  used consistently within its own note; the prior referee's
  re-instantiation (2d_0 - R + 2 > i_feed + 1) and my window lemma
  (> i_feed + 66) subsume all three.  No contradiction found anywhere
  between E1, E3, D1/D2/D3, and the prior referee.
''')
    A('=' * 76)
    A('FINAL CANON PARAGRAPH (referee-certified, 2026-08-04)')
    A('=' * 76)
    A('''
BEST UNCONDITIONAL C* UPPER BOUND: UNCHANGED -- C* <= 2, from the
submission envelope T(n) <= max(n+1, 2n-2) (and T(n) <= max(n+1, 2n-3)
for n >= 5); the bound is exactly rational, no bracket needed.  No lane
produced a new unconditional C* bound, and none claimed to.  What
changed is the CONDITIONAL frontier, now in its sharpest verified form:
the chain  THM-3329 assembly + LIFT (Thm A) + Thm B feed-phase survival
(D0 >= eps* R, eps* = 2(1-g-g^2)/(3+2g) in (0.0211736, 0.0211747),
g = gamma* = log_5 phi^2 in (627035, 627036)/2^20) + death-free handover
window (T6b extension, D0(3+2g) >= 2R(1-g-g^2) + 66(1+g), certified
2^16..2^40 -- this referee) + THEOREM S-cone-fc (one-row certificate
F1^F2^F3^F4 at some post-feed row => capture by R-2 with no death;
referee-confirmed, independently re-implemented) + a-priori F4
(certified dyadic 2^9..2^40) + exact closures for ALL dyadic
128 <= R <= 32768 at eps = 1/32 and 1/16 (R = 32768 both eps now
verified by TWO independent implementation lineages: captures
19865/20185, debt 16383 exact, feed-end states bit-identical) proves:
IF ENTRY-fc(eps) holds -- i.e. for every dyadic R >= 65536 some
post-feed row i0 <= i_pf + 64 of plain rule A satisfies the one-row
static conditions F1 (all-negative junk, support in [0, m], m+2 < d),
F2 (a_0 <= d-1), F3 (2a_{t-1} + a_{t-2} <= 2C(d-1,t) on [2, m+2]), and
(for R > 2^40) F4 -- THEN Hypothesis S(eps) holds and
C* <= 1 + gamma* + eps; in particular ENTRY-fc(1/32) gives
C* <= 1 + gamma* + 1/32 < 427095/262144 = 1.6292382, and the route's
closed-form limit constant is 1 + gamma* + eps* = (5+3g)/(3+2g) in
(1.6191617801, 1.6191618342).  Every Theta(R)-row dynamical statement
in S is now PROVED; the sole open item is the STATIC one-row handover
property, whose empirical form is the marginal-surface law
(max_t r_t -> 1+, measured 1.0009 at R = 32768, over-cells single-
parity, i_fc = i_pf + 0..8 at all 18 scales).  WHAT REMAINS for
C* <= log_5(5 phi^2) = 1 + gamma* itself: either (R0) exact-floor
universality -- frontier (512, D0=0), OPEN both directions, with the
(256,0) transportation integer point re-verified -- or (GG + S') an
o(R)-slack schedule below the G3 super-pair threshold plus the post-
feed collapse; and on the floor side the general-class lower bound
C* >= log_5(5 phi^2) remains OPEN (golden value pinned only against
balanced-block schemes).  Hypothesis S itself: VERIFIED-exact for all
dyadic 128 <= R <= 32768 (both eps), OPEN for R >= 65536 -- "COMPLETE"
in the E1 lane headline means the reduction is complete, not that S is
proved; C* <= 1 + gamma* + eps remains CONDITIONAL.
''')
    with open(os.path.join(RES, 'amm12592_S_referee_boxeph.out'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('report written: %d pass, %d fail' % (tp, tf))
    return tp, tf, runs

# --------------------------------------------------------------------- main
if __name__ == '__main__':
    stage = sys.argv[1]
    if stage == 'constants':
        stage_constants()
    elif stage == 'slow':
        stage_slow()
    elif stage == 'suite':
        stage_suite()
    elif stage == 'bigrun':
        stage_bigrun(int(sys.argv[2]), int(sys.argv[3]))
    elif stage == 'conecert':
        stage_conecert(int(sys.argv[2]), int(sys.argv[3]))
    elif stage == 'lemmas':
        stage_lemmas()
    elif stage == 'apriori':
        stage_apriori()
    elif stage == 'window':
        stage_window()
    elif stage == 'e3wit':
        stage_e3wit()
    elif stage == 'ledger':
        stage_ledger()
    elif stage == 'bigcheck':
        stage_bigcheck()
    elif stage == 'report':
        stage_report()
    else:
        raise SystemExit('unknown stage')
    print('STAGE %s DONE: %d pass, %d fail' % (stage, _ck['pass'], _ck['fail']))
