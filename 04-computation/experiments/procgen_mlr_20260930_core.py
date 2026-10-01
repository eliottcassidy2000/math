#!/usr/bin/env python3
"""procgen_mlr_20260930_core.py -- shared library for the "multiplicative lonely runner" lane (mlr).

Session collatz-procgen-20260922, lane "mlr", 2026-09-30.

Objects
  * modulus D coprime to 6; runner box B(J,K) = {2^j 3^k : 0 <= j < J, 0 <= k < K};
  * a time h in Z/D; the runner (j,k) sits at 2^j 3^k h / D mod 1; its centred residue
    r(j,k) in (-D/2, D/2] is the representative of 2^j 3^k h mod D (negative exponents = inverses);
  * near set at integer threshold R (delta = R/D): {(j,k) : |r(j,k)| < R};
  * L(D,J,K) = max_h min_{(j,k) in B} ||2^j 3^k h/D||   (maximal loneliness, all h != 0);
  * I(t) = inf_{j,k >= 0} ||2^j 3^k t||  (the x2x3 "lonely spectrum" of a real t).

Everything that is labelled PROVED/FINITE-EXACT in the report is computed here in exact integer or
Fraction arithmetic; float is used only for summaries and for the gate exponential sums (which are
certified by the gates lane, procgen_gates_20260925_core.py, imported read-only).
"""
import bisect
import math
import os
import resource
import sys
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

LOG2, LOG3 = math.log(2.0), math.log(3.0)


class CheckFailed(AssertionError):
    pass


def check(cond, msg):
    if not cond:
        raise CheckFailed(msg)


def mem_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 * 1024) if sys.platform == "darwin" else r / 1024


def lean_malloc():
    """macOS: re-exec once with MallocLargeCache=0 (returns freed large blocks to the OS)."""
    if sys.platform == "darwin" and os.environ.get("MallocLargeCache") != "0":
        env = dict(os.environ)
        env["MallocLargeCache"] = "0"
        os.execve(sys.executable, [sys.executable, "-u"] + sys.argv, env)


# ----------------------------------------------------------------------------------------------
# elementary arithmetic
# ----------------------------------------------------------------------------------------------
def is_prime(n):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0:
            return n == sp
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for b in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(b, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def next_prime(n):
    while not is_prime(n):
        n += 1
    return n


def centred(x, D):
    x %= D
    return x - D if x > D // 2 else x


def mult(j, k, D):
    """2^j 3^k mod D for integers j, k of any sign (D coprime to 6)."""
    a = pow(2, j, D) if j >= 0 else pow(pow(2, -1, D), -j, D)
    b = pow(3, k, D) if k >= 0 else pow(pow(3, -1, D), -k, D)
    return a * b % D


def smooth_upto(X):
    """all 3-smooth integers 2^a 3^b <= X, sorted."""
    out = []
    a = 1
    while a <= X:
        b = a
        while b <= X:
            out.append(b)
            b *= 3
        a *= 2
    return sorted(out)


def sigma3(X):
    """#{(a,b) >= 0 : 2^a 3^b < X}  (exact, X a Fraction/int/float > 0)."""
    if X <= 1:
        return 0
    c = 0
    b = 1
    while b < X:
        a = b
        while a < X:
            c += 1
            a *= 2
        b *= 3
    return c


# ----------------------------------------------------------------------------------------------
# runner boxes: best numerator, near counts, loneliness
# ----------------------------------------------------------------------------------------------
def box_multipliers(D, J, K, k0=0):
    """2^j 3^k mod D for 0 <= j < J, k0 <= k < K (k may be negative: inverse powers)."""
    return np.array([mult(j, k, D) for j in range(J) for k in range(k0, K)], dtype=np.int64)


def scan_box(D, J, K, R_list=(), k0=0, hs=None):
    """For every h (or the given hs): best[h] = min over the box of |r(j,k)| and, for each integer
    threshold R in R_list, cnt[R][h] = #{(j,k) in box : |r(j,k)| < R}.  O(|box| * #h) int64 ops."""
    G = box_multipliers(D, J, K, k0)
    h = np.arange(D, dtype=np.int64) if hs is None else np.asarray(hs, dtype=np.int64)
    best = np.full(len(h), D, dtype=np.int64)
    cnt = {R: np.zeros(len(h), dtype=np.int32) for R in R_list}
    check(D < (1 << 31), "scan_box: D too large for int64 products")
    for g in G:
        r = (g * h) % D
        c = np.minimum(r, D - r)
        np.minimum(best, c, out=best)
        for R in R_list:
            cnt[R] += (c < R)
    return best, cnt


def loneliness(D, J, K):
    """L(D,J,K) as an exact Fraction (max over h != 0 of min over box of |r|/D) and an argmax."""
    best, _ = scan_box(D, J, K)
    best[0] = -1
    h = int(best.argmax())
    return F(int(best[h]), D), h


# ----------------------------------------------------------------------------------------------
# Triangle lemma machinery (Z^2, exact residues)
# ----------------------------------------------------------------------------------------------
def near_components(D, h, R, j_range, k_range):
    """Near set of h at threshold R on the rectangle j_range x k_range of Z^2 (exponents of any
    sign), its 4-connected components, and per component: the meet (min j, min k), the residue
    there, and whether the component equals the full up-right triangle of its meet (within the
    rectangle).  Returns (res dict, comps list)."""
    res = {}
    for j in j_range:
        for k in k_range:
            res[(j, k)] = centred(h * mult(j, k, D), D)
    near = {P for P, r in res.items() if abs(r) < R}
    comps = []
    seen = set()
    for P in sorted(near):
        if P in seen:
            continue
        stack = [P]
        comp = []
        seen.add(P)
        while stack:
            Q = stack.pop()
            comp.append(Q)
            j, k = Q
            for N in ((j + 1, k), (j - 1, k), (j, k + 1), (j, k - 1)):
                if N in near and N not in seen:
                    seen.add(N)
                    stack.append(N)
        comps.append(comp)
    return res, near, comps


def verify_triangle_lemma(D, h, R, jlo, jhi, klo, khi):
    """(Triangle Lemma holds for 1 <= R <= D/4.)  Check, on a window W = [jlo,jhi) x [klo,khi) of Z^2 large enough to contain every component
    whose meet lies in the inner window (margin handled by caller), the claims of the Triangle Lemma:
      (b) each component meeting the inner window has a unique minimum P0 and equals
          T(P0) = {P0 + (x,y) : 2^x 3^y |u| < R}, with exact residues r = 2^x 3^y u;
      (c) the apex numerator u = r(P0) is coprime to 6;
      (d) the shadows Sh(P0) = {P0 + (x,y) : R <= 2^x 3^y |u| < D/2} carry exact residues, are pairwise
          disjoint and disjoint from the near set.
    Returns number of components checked."""
    jr, kr = range(jlo, jhi), range(klo, khi)
    res, near, comps = near_components(D, h, R, jr, kr)
    shadow_owner = {}
    ncheck = 0
    for comp in comps:
        j0 = min(P[0] for P in comp)
        k0 = min(P[1] for P in comp)
        # skip components touching the window boundary (may be truncated)
        if j0 == jlo or k0 == klo or max(P[0] for P in comp) >= jhi - 1 or max(P[1] for P in comp) >= khi - 1:
            continue
        P0 = (j0, k0)
        check(P0 in near, f"meet not near D={D} h={h} R={R} comp={comp[:4]}")
        for (P, Q) in ((P, (P[0] + 1, P[1])) for P in comp):
            if Q in near:
                check(res[Q] == 2 * res[P], "adjacent near points without exact x2 relation")
        for (P, Q) in ((P, (P[0], P[1] + 1)) for P in comp):
            if Q in near:
                check(res[Q] == 3 * res[P], "adjacent near points without exact x3 relation")
        u = res[P0]
        check(u % 2 != 0 and u % 3 != 0, f"apex numerator not 6-free D={D} h={h} u={u}")
        T = set()
        x = 0
        while (1 << x) * abs(u) < R:
            y = 0
            while (1 << x) * 3 ** y * abs(u) < R:
                Q = (j0 + x, k0 + y)
                T.add(Q)
                if Q in res:
                    check(res[Q] == (1 << x) * 3 ** y * u, "inexact residue in triangle")
                y += 1
            x += 1
        check(set(comp) == T, f"component != triangle D={D} h={h} R={R} P0={P0}")
        # shadow
        x = 0
        while (1 << x) * abs(u) < D / 2:
            y = 0
            while (1 << x) * 3 ** y * abs(u) < D / 2:
                v = (1 << x) * 3 ** y * abs(u)
                if v >= R:
                    Q = (j0 + x, k0 + y)
                    if Q in res:
                        check(res[Q] == (1 << x) * 3 ** y * u, "inexact residue in shadow")
                        check(Q not in near, "shadow point is near")
                        check(Q not in shadow_owner, f"shadows overlap D={D} h={h} at {Q}")
                        shadow_owner[Q] = P0
                y += 1
            x += 1
        ncheck += 1
    return ncheck


def subgroup23(D):
    """sorted array of the elements of <2,3> in (Z/D)^*."""
    seen = np.zeros(D, dtype=bool)
    seen[1 % D] = True
    frontier = [1 % D]
    while frontier:
        nf = []
        for x in frontier:
            for g in (2, 3):
                y = (x * g) % D
                if not seen[y]:
                    seen[y] = True
                    nf.append(y)
        frontier = nf
    return np.nonzero(seen)[0].astype(np.int64)


def cosets23(D):
    """list of (representative h, coset array) for the cosets of H = <2,3> in (Z/D)^*."""
    H = subgroup23(D)
    done = np.zeros(D, dtype=bool)
    out = []
    for h in range(1, D):
        if done[h] or math.gcd(h, D) != 1:
            continue
        cos = (h * H) % D
        done[cos] = True
        out.append((h, cos))
    return H, out


def coset_crowding_bound(D, R, m):
    """PROVED upper bound for the near fraction of a coset whose least |element| is m:
    max over 6-free u in [m, R) of sigma3(R/u) / sigma3(D/(2u))."""
    best = F(0)
    for u in range(m, R):
        if u % 2 == 0 or u % 3 == 0:
            continue
        num = sigma3(F(R, u))
        den = sigma3(F(D, 2 * u))
        if den > 0:
            best = max(best, F(num, den))
    return best


# ----------------------------------------------------------------------------------------------
# discrete counting criterion (6-free reduction of the union bound)
# ----------------------------------------------------------------------------------------------
def bad_set_bound(J, K, R):
    """F_{J,K}(R) = sum_{1<=u<R, (u,6)=1} 2 [JK + K A_u + J B_u + tau^+_u], where
    A_u = max{a : 2^a u < R}, B_u = max{b : 3^b u < R}, tau^+_u = #{a,b >= 1 : 2^a 3^b u < R}.
    PROVED: #{h != 0 : some runner of the box has 0 < |r| < R} <= F_{J,K}(R)."""
    tot = 0
    for u in range(1, R):
        if u % 2 == 0 or u % 3 == 0:
            continue
        A = 0
        while (u << (A + 1)) < R:
            A += 1
        B = 0
        while u * 3 ** (B + 1) < R:
            B += 1
        tp = 0
        a = 1
        while (u << a) * 3 < R:
            b = 1
            while (u << a) * 3 ** b < R:
                tp += 1
                b += 1
            a += 1
        tot += 2 * (J * K + K * A + J * B + tp)
    return tot


def certified_R(D, J, K):
    """largest integer R >= 1 with F_{J,K}(R) < D - 1 (then L(D,J,K) >= R/D); 0 if none."""
    lo = 0
    R = 1
    while bad_set_bound(J, K, R) < D - 1:
        lo = R
        R = R * 2
        if R > D:
            break
    hi = R
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if bad_set_bound(J, K, mid) < D - 1:
            lo = mid
        else:
            hi = mid
    return lo


def bad_set_exact(D, J, K, R):
    """exact #{h != 0 : min over box of |r| < R} (for checking bad_set_bound)."""
    best, _ = scan_box(D, J, K)
    best[0] = D
    return int((best < R).sum())


# ----------------------------------------------------------------------------------------------
# continuous loneliness: exact interval arithmetic
# ----------------------------------------------------------------------------------------------
def allowed_intersect(intervals, v, lo, hi):
    """intersect a sorted list of disjoint closed intervals [a,b] (Fractions) with
    {t : v t mod 1 in [lo, hi]} = union_m [(m+lo)/v, (m+hi)/v]."""
    out = []
    for a, b in intervals:
        m = max(int(math.floor(a * v - hi)) - 1, -1)
        while True:
            L = (m + lo) / v
            Rr = (m + hi) / v
            if L > b:
                break
            x, y = max(a, L), min(b, Rr)
            if x <= y:
                out.append((x, y))
            m += 1
    return out


def X_set(speeds, delta):
    """{t in [0,1] : ||v t|| >= delta for all v in speeds} as a list of closed Fraction intervals."""
    iv = [(F(0), F(1))]
    for v in sorted(set(speeds)):
        iv = allowed_intersect(iv, v, delta, 1 - delta)
    return iv


def kappa_exact(V):
    """kappa(V) = sup_t min_{v in V} ||t v|| (exact Fraction) via the classical reduction: the sup is
    attained at t = l/(v+v') (Perarnau-Serra survey eq. (4)); we scan every N = v+v' and all l."""
    V = sorted(set(V))
    Ns = sorted({v + w for v in V for w in V})
    best = F(0)
    arg = None
    Varr = np.array(V, dtype=np.int64)
    for N in Ns:
        l = np.arange(1, N, dtype=np.int64)
        r = (l[:, None] * Varr[None, :]) % N
        m = np.minimum(r, N - r).min(axis=1)
        i = int(m.argmax())
        val = F(int(m[i]), N)
        if val > best:
            best, arg = val, F(int(l[i]), N)
    return best, arg


def I_rat(t):
    """I(t) = inf_{j,k>=0} ||2^j 3^k t|| for a rational t (the x2x3 orbit of t is finite)."""
    r, N = t.numerator % t.denominator, t.denominator
    seen = {r}
    frontier = [r]
    best = N
    while frontier:
        nf = []
        for x in frontier:
            best = min(best, x, N - x)
            for g in (2, 3):
                y = (g * x) % N
                if y not in seen:
                    seen.add(y)
                    nf.append(y)
        frontier = nf
    return F(best, N)


def _stabilisers(N, limit=10 ** 7):
    N0, Np = 1, N
    while Np % 2 == 0:
        Np //= 2
        N0 *= 2
    while Np % 3 == 0:
        Np //= 3
        N0 *= 3
    return N0, Np, [s for s in smooth_upto(limit) if s > 1 and Np > 1 and s % Np == 1]


def local_cover(t0, sigma, w, delta, n_stab=6):
    """PROVE: for all x in (0, w], I(t0 + sigma x) < delta.
    Method: pick s* in <2,3> with s* = 1 mod N' (N' = 6-free part of the denominator N of t0) and use
    multipliers u in N0*<2,3> (N0 = {2,3}-part of N), for which (u s*^i) t0 = u t0 mod 1.  It suffices
    to cover the closed interval [w/s*, w] by the open sets J_u = {x : ||u t0 + sigma u x|| < delta}
    (then x in (0,w] is moved into [w/s*, w] by s*^i).  Exact Fractions throughout.
    Returns (ok, s*, max multiplier u used)."""
    N0, Np, cands = _stabilisers(t0.denominator)
    if Np == 1:
        return False, None, None
    for s in cands[:n_stab]:
        lo_x, hi_x = w / s, w
        U = [N0 * u for u in smooth_upto(int(20 * s / w) + 10)]
        pieces = []
        for u in U:
            y = (u * t0) % 1
            zlo, zhi = u * lo_x, u * hi_x
            if zlo > 10:
                break
            mlo = math.floor(min(y + sigma * zlo, y + sigma * zhi)) - 1
            mhi = math.ceil(max(y + sigma * zlo, y + sigma * zhi)) + 1
            for m in range(mlo, mhi + 1):
                A_, B_ = m - delta - y, m + delta - y
                z1, z2 = (A_, B_) if sigma > 0 else (-B_, -A_)
                x1, x2 = z1 / u, z2 / u
                if max(x1, lo_x) < min(x2, hi_x):
                    pieces.append((x1, x2, u))
        pieces.sort()
        pos, umax, ok = lo_x, 0, True
        while True:
            best = None
            for (p1, p2, u) in pieces:
                if p1 < pos < p2 and (best is None or p2 > best[0]):
                    best = (p2, u)
            if best is None:
                ok = False
                break
            umax = max(umax, best[1])
            if best[0] > hi_x:
                break
            pos = best[0]
        if ok:
            return True, s, umax
    return False, None, None


def certify_spectrum(delta, J, K, Dcen=300, far=F(1, 1000)):
    """Computer-assisted proof that every t in R/Z either has I(t) < delta or is one of the returned
    exceptional rationals (all with I >= delta).  Steps: (1) box certificate X = X_set(box, delta):
    t outside X has some ||v t|| < delta; (2) each component of X is assigned to a centre r/N
    (N <= Dcen) inside it or within `far`; (3) local_cover around every centre on both sides.
    Returns dict with keys ok, exceptional (list of (t0, I(t0))), ncomp, ncentre, unassigned,
    worst_multiplier (max over centres/sides of umax * s* * w * N, used for the discrete corollary)."""
    B = [2 ** j * 3 ** k for j in range(J) for k in range(K)]
    X = X_set(B, delta)
    cen = sorted({F(r, N) for N in range(2, Dcen + 1) for r in range(1, N) if math.gcd(r, N) == 1})
    assign = {}
    unassigned = []
    for (a, b) in X:
        i = bisect.bisect_left(cen, a)
        inside = [c for c in cen[max(0, i - 1):i + 6] if a <= c <= b]
        if inside:
            c = min(inside, key=lambda c: c.denominator)
        else:
            cands = [cen[k] for k in (i - 1, i) if 0 <= k < len(cen)]
            c = min(cands, key=lambda c: min(abs(c - a), abs(c - b)))
            if min(abs(c - a), abs(c - b)) > far:
                unassigned.append((a, b))
                continue
        wl, wr = assign.get(c, (F(0), F(0)))
        assign[c] = (max(wl, c - a), max(wr, b - c))
    ok = not unassigned
    exc = []
    nonexc = []
    worst = 0
    fails = []
    for c in sorted(assign):
        Ic = I_rat(c)
        if Ic >= delta:
            exc.append((c, Ic))
        else:
            nonexc.append((c, Ic))
        wl, wr = assign[c]
        for sigma, w in ((-1, wl), (1, wr)):
            if w <= 0:
                continue
            good, s, umax = local_cover(c, sigma, w, delta)
            if not good:
                ok = False
                fails.append((c, sigma, w))
            else:
                worst = max(worst, umax * s * w * c.denominator)
    return dict(ok=ok, exceptional=exc, nonexceptional=nonexc, ncomp=len(X), ncentre=len(assign),
                unassigned=unassigned, worst_multiplier=worst, fails=fails)


# ----------------------------------------------------------------------------------------------
# gate helpers (read-only use of the gates lane's certified DP)
# ----------------------------------------------------------------------------------------------
def gates_core():
    import procgen_gates_20260925_core as gc
    return gc


def gate_spectrum_fft(p, a):
    """full spectrum S(h), h mod M, via the residue histogram of c_w mod M and one inverse FFT
    (S(h) = sum_r n_r e(hr/M) = M * ifft(n)[h]).  Only for M <= 2^20 (memory)."""
    gc = gates_core()
    M = abs(gc.gate(p, a))
    check(M <= (1 << 20), "gate_spectrum_fft: M too large")
    n = np.zeros(M, dtype=np.float64)
    for blk, _ in gc.residue_chunks(p, a, want_float=False):
        n += np.bincount(blk, minlength=M)
    S = np.fft.ifft(n) * M
    return M, S, n


def gate_spectrum_dp(p, a, block=4096):
    """full spectrum S(h), h mod M, by the gates lane's certified float DP (blocks of h, h <= M/2, then
    conjugate symmetry S(M-h) = conj S(h)), plus the exact residue histogram n_r (for N = n_0 and
    Coll = sum n_r^2).  Memory ~ 40 bytes per residue (no FFT)."""
    gc = gates_core()
    M = abs(gc.gate(p, a))
    n = np.zeros(M, dtype=np.float64)
    for blk, _ in gc.residue_chunks(p, a, want_float=False):
        n += np.bincount(blk, minlength=M)
    S = np.empty(M, dtype=np.complex128)
    H = M // 2 + 1
    for h0 in range(0, H, block):
        h1 = min(H, h0 + block)
        S[h0:h1] = gc.S_dp(p, a, range(h0, h1))
    S[H:] = np.conj(S[1:M - H + 1][::-1])
    return M, S, n


def path_weights(p, a):
    """W[s, m] = P(the one with 3-exponent m = a-1-i sits at position s) for a uniform word."""
    C = math.comb(p, a)
    W = np.zeros((p, a))
    for i in range(a):
        for s in range(i, p - (a - 1 - i)):
            W[s, a - 1 - i] = math.comb(s, i) * math.comb(p - 1 - s, a - 1 - i) / C
    return W
