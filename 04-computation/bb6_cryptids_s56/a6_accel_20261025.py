#!/usr/bin/env python3
"""A6 exact two-level accelerated simulator with certified translated-cycle jumps (S56).

Machine A6 = 1RB0LF_1RC1RA_0RD0RB_1LE---_0LA1LE_1RD1LA.  Phase P(v0, w0) = J(0, v0, w0) and its end cases are
in a6_phase_map_20261025.py (verified against the paper's rule system for every W < 20000).

Jump principle (three-point lemma).  Fix a finite label sequence L = (l_1..l_p) of phases (each label = depth j and
end case).  The set D_L of states from which the next p phases carry exactly the labels L is (convex polyhedron) cap
(coset of a lattice) -- depth and end-case conditions are linear inequalities, parities are congruences -- and on D_L
the p-step map G_L is affine.  If b, b+D, b+K*D all lie in D_L with G_L(b)=b+D, G_L(b+D)=b+2D, then for every
0 <= k <= K:  b+k*D lies in D_L and G_L(b+k*D) = b+(k+1)*D.  (Convexity gives the polyhedron along the segment; b and
b+D in one coset puts D in the lattice; the linear part fixes D.)  The largest such K is found by exponential +
binary search; if K = 2^E holds with E far above the bit sizes involved, no inequality decreases along the ray, so
the chain is infinite: a certified non-returning sweep (a translated cycler at the phase level).
Level 1 applies this to phases inside a sweep (states (v0, w0)); level 2 to whole sweeps (states W, labels = the
tokenised phase-label sequence of each sweep, which fixes the composite affine map).
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a6_phase_map_20261025 import phase


class Capture(Exception):
    """Certified infinite chain inside a sweep: args = (base state, translation, label tuple, exponent E)."""


class Halt(Exception):
    pass


def big_exponent(*xs, p=1):
    bits = max(int(abs(x)).bit_length() for x in xs) if xs else 1
    return max(1000, 4 * (bits + 2 * p) + 200)


# ---------------------------------------------------------------- level 1: phases inside a sweep
def _run_labels(s, labels):
    """Apply phases from state s; return the final state if the labels match exactly and all results are 'P'."""
    v0, w0 = s
    for L in labels:
        if v0 < 0 or w0 < 0:
            return None
        res, j, lab = phase(v0, w0)
        if (j, lab) != L or res[0] != 'P':
            return None
        v0, w0 = res[1], res[2]
    return (v0, w0)


def _valid1(b, D, labels, K):
    s = (b[0] + K * D[0], b[1] + K * D[1])
    return _run_labels(s, labels) == (b[0] + (K + 1) * D[0], b[1] + (K + 1) * D[1])


def _max_valid(valid, E):
    """valid(0), valid(1) known true; the valid set is an interval [0, Kmax]. Return Kmax, or None if valid(2^E)."""
    if valid(1 << E):
        return None
    lo, hi = 1, 2
    while valid(hi):
        lo, hi = hi, 2 * hi
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if valid(mid):
            lo = mid
        else:
            hi = mid
    return lo


STATS = {'l1_jumps': 0, 'l1_jumped_phases': 0, 'l2_jumps': 0, 'l2_jumped_sweeps': 0, 'phases': 0}


def sweep(W, cap=2 * 10 ** 6, pmax=4000, gram=3):
    """Sweep from J(0,0,W). Returns (W', tokens) or raises Halt(tokens) / Capture(cert)."""
    xs = [(0, W)]          # states x_0, x_1, ...
    labs = []              # labs[i] = label of the phase x_i -> x_{i+1}
    tokens = []
    occ = {}
    n_ph = 0
    while True:
        v0, w0 = xs[-1]
        res, j, lab = phase(v0, w0)
        n_ph += 1
        STATS['phases'] += 1
        L = (j, lab)
        if res[0] == 'HALT':
            tokens.append(L)
            raise Halt(tuple(tokens))
        if res[0] == 'ret':
            tokens.append(L)
            return res[1], tuple(tokens)
        labs.append(L)
        tokens.append(L)
        xs.append((res[1], res[2]))
        if n_ph > cap:
            raise RuntimeError(('too many phases', xs[-1]))
        n = len(labs)   # x_n is the newest state
        if n < 2 * gram:
            continue
        key = tuple(labs[-gram:])
        lst = occ.setdefault(key, [])
        found = False
        for i in reversed(lst):
            p = n - i
            if p > pmax or 2 * p > n:
                break
            xn, xp, x2p = xs[n], xs[n - p], xs[n - 2 * p]
            D = (xn[0] - xp[0], xn[1] - xp[1])
            if D != (xp[0] - x2p[0], xp[1] - x2p[1]) or D == (0, 0):
                continue
            if labs[n - p:n] != labs[n - 2 * p:n - p]:
                continue
            if any(isinstance(l[0], str) for l in labs[n - 2 * p:n]):
                continue          # windows containing earlier jump tokens are not used at level 1
            cyc = tuple(labs[n - p:n])
            b = x2p
            # constraint values along the cycle are integer-affine in k with |value at k=0| < 2^(2*maxbits+O(1))
            maxbits = max(max(abs(a).bit_length(), abs(c).bit_length()) for a, c in xs[n - 2 * p:n + 1])
            E = max(1000, 3 * maxbits + 300)
            Kmax = _max_valid(lambda K: _valid1(b, D, cyc, K), E)
            if Kmax is None:
                raise Capture((b, D, cyc, E))
            # jump: from b, (Kmax+1) periods lead to b + (Kmax+1) D
            new = (b[0] + (Kmax + 1) * D[0], b[1] + (Kmax + 1) * D[1])
            STATS['l1_jumps'] += 1
            STATS['l1_jumped_phases'] += (Kmax + 1) * p
            # rewrite history: drop the 2p phases after b, record a token
            del xs[n - 2 * p + 1:]
            del labs[n - 2 * p:]
            del tokens[len(tokens) - 2 * p:]
            tokens.append(('cyc', cyc, Kmax + 1))
            xs.append(new)
            labs.append(('cyc', cyc, Kmax + 1))
            occ = {}
            found = True
            break
        if not found:
            lst.append(n)
            if len(lst) > 8:
                del lst[0]


# ---------------------------------------------------------------- level 2: whole sweeps
def _run_sweeps(W, toks):
    for T in toks:
        if W < 0:
            return None
        try:
            W2, t = sweep(W)
        except (Halt, Capture, RuntimeError):
            return None
        if t != T:
            return None
        W = W2
    return W


def orbit(W0=0, max_sweeps=10 ** 9, pmax=3000, gram=2, report=None, stop=None):
    """Iterate sweeps from W0 with level-2 jumps. Yields nothing; returns (status, W, n_sweeps, info)."""
    Ws = [W0]
    toks = []
    occ = {}
    nsw = 0
    while nsw < max_sweeps:
        W = Ws[-1]
        if stop is not None and stop(W, nsw):
            return 'stopped', W, nsw, None
        try:
            W2, t = sweep(W)
        except Halt as h:
            return 'HALT', W, nsw, h.args[0]
        except Capture as c:
            return 'CAPTURE', W, nsw, c.args[0]
        nsw += 1
        toks.append(t)
        Ws.append(W2)
        if report is not None:
            report(W, W2, t, nsw)
        n = len(toks)
        if len(Ws) > 20000:          # keep memory bounded
            cut = 10000
            del Ws[:cut]
            del toks[:cut]
            occ = {}
            continue
        if n < 2 * gram:
            continue
        key = tuple(toks[-gram:])
        lst = occ.setdefault(key, [])
        found = False
        for i in reversed(lst):
            p = n - i
            if p > pmax or 2 * p > n:
                break
            D = Ws[n] - Ws[n - p]
            if D != Ws[n - p] - Ws[n - 2 * p] or D == 0:
                continue
            if toks[n - p:n] != toks[n - 2 * p:n - p]:
                continue
            cyc = tuple(toks[n - p:n])
            b = Ws[n - 2 * p]
            maxbits = max(abs(x).bit_length() for x in Ws[n - 2 * p:n + 1])
            E = max(1000, 3 * maxbits + 300)

            def valid(K):
                return _run_sweeps(b + K * D, cyc) == b + (K + 1) * D
            Kmax = _max_valid(valid, E)
            if Kmax is None:
                return 'INFINITE-SWEEP-CYCLE', b, nsw, (b, D, cyc, E)
            m = min(Kmax + 1, (max_sweeps - nsw + 2 * p) // p)   # never overshoot max_sweeps
            new = b + m * D
            STATS['l2_jumps'] += 1
            STATS['l2_jumped_sweeps'] += m * p - 2 * p
            nsw += m * p - 2 * p
            del Ws[n - 2 * p + 1:]
            del toks[n - 2 * p:]
            toks.append(('cyc2', len(cyc), m))
            Ws.append(new)
            occ = {}
            found = True
            break
        if not found:
            lst.append(n)
            if len(lst) > 8:
                del lst[0]
    return 'max', Ws[-1], nsw, None


if __name__ == '__main__':
    import time
    from a6_fast_20261025 import sweep as sweep_fast, TooManyPhases
    # self-test 1: level-1 sweep equals the fast sweep on W < 20000 (returning sweeps) and flags the known chains
    bad = []
    caps = []
    for W in range(20000):
        try:
            a = sweep(W)[0]
        except Halt:
            a = 'HALT'
        except Capture:
            caps.append(W)
            continue
        try:
            b = sweep_fast(W, cap=200000)[0]
        except TooManyPhases:
            b = 'INF'
        if a != b:
            bad.append((W, a, b))
    print('level-1 sweep == fast sweep for W < 20000:', not bad, bad[:3])
    print('certified captures W < 20000:', len(caps), caps[:25])
    # self-test 2: orbit with level-2 jumps vs plain fast sweeps for the first 3e6 sweeps
    t0 = time.time()
    W = 0
    plain = {}
    for i in range(3 * 10 ** 6):
        W = sweep_fast(W)[0]
    print('plain fast orbit after 3e6 sweeps: W =', W, f'({time.time() - t0:.0f}s)')
    t0 = time.time()
    st = orbit(0, max_sweeps=3 * 10 ** 6)
    print('accelerated orbit after', st[2], 'sweeps: status', st[0], 'W =', st[1], f'({time.time() - t0:.0f}s)', STATS)
