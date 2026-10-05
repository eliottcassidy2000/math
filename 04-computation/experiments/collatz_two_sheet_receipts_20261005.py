#!/usr/bin/env python3
"""
The two-sheet receipt calculus and the trunk multiplier move (opus, 2026-10-05).

SHEETS.  U_+(u) = oddpart(3u+1) and U_-(u) = oddpart(3u-1) on positive odd u.  An edge (u, s) with sheet s in {+,-}
goes from the label [u] to [U_s(u)] and has value u/U_s(u).  A two-sheet receipt is an integer-valued finite function c
on edges with  c >= 0 on the + sheet (actual Collatz stock)  and  c <= 0 on the - sheet (reversed 3x-1 edges =
multipliers).  Value V(c) = prod value^c, boundary d c = sum c ([u] - [U_s u]).  A receipt for the ratio x/x' has
transport defect  Delta = d c - ([x] - [x']).

MULTIPLIER MOVE.  If m is 3x-1-ROOTED (its - orbit reaches 1; the - map also has the cycles {5,7} and
{17,25,37,55,41,61,91}), the reversed - path P_-(m) has value 1/m and boundary -[m] + [1].  For any + path from m x
to x' the packet  c = (+ path from m x) - P_-(m)  has value x/x' and transport defect exactly the ONE fusion relation
R(m, x) = [m x] + [1] - [m] - [x]  (Proposition A below).  The TRUNK MOVE is the case m = m_j = (2^j + 1)/3, j odd
(3 m_j - 1 = 2^j, so P_-(m_j) is the single edge m_j -> 1), applied to x = -1 mod 2^j, where U_+(m_j x) = x + (x+1)/2^j
exactly (Applegate--Lagarias Lemma 2.2 in accelerated form).

TESTS.
A. exact checks of the move on many (x, j) and of Proposition A on random (m, x, depth);
B. a prefix code of residue classes mod 2^K (K <= Kmax) paying every class by a plain descent or by ONE multiplier
   move with m from a given set, searched exactly on affine classes (descent for all members, AL's worst-case
   criterion); compared for M_minus (3x-1-rooted m <= 100) against M_AL = {5,7,11,13,23,29,43,25,35};
C. universal two-sheet receipts for every odd x <= N built stage by stage from the code: number of fusion relations
   (multiplier events) per x, i.e. the defect rank of universal coverage, for both codes;
D. the trunk move as a dynamical shortcut: for x in the deep cells -1 mod 2^k, does the + orbit of x meet the + orbit
   of m_j x before x descends below itself?
E. the mixed trunk: x = (4^e - 1)/(2^j + 1) (j | e, j odd) has m_j x on the + trunk, so the trunk move lands at 1.
Usage: python collatz_two_sheet_receipts_20261005.py [Kmax=16] [N=65536]
"""
import sys, math, time, random
from fractions import Fraction
from collections import Counter, defaultdict

def oddpart(n):
    while n % 2 == 0:
        n //= 2
    return n

def Up(u): return oddpart(3 * u + 1)
def Um(u): return oddpart(3 * u - 1) if u > 1 else 1
def v2(n): return (n & -n).bit_length() - 1

def minus_orbit(m, cap=10000):
    """- orbit of m until 1 or a repeat; returns (path, rooted?)"""
    path = [m]; seen = {m}
    while m != 1 and len(path) < cap:
        m = Um(m)
        if m in seen:
            return path, False
        seen.add(m); path.append(m)
    return path, (m == 1)

def minus_rooted(m):
    return minus_orbit(m)[1]

def plus_path(x, steps):
    p = [x]
    for _ in range(steps):
        x = Up(x); p.append(x)
    return p

# ---------------- receipts ----------------
class Receipt:
    def __init__(self):
        self.c = Counter()
    def add_plus_path(self, path, mult=1):
        for u in path[:-1]:
            self.c[(u, '+')] += mult
    def sub_minus_path(self, path):
        for u in path[:-1]:
            self.c[(u, '-')] -= 1
    def value(self):
        v = Fraction(1)
        for (u, s), k in self.c.items():
            img = Up(u) if s == '+' else Um(u)
            v *= Fraction(u, img) ** k
        return v
    def boundary(self):
        b = Counter()
        for (u, s), k in self.c.items():
            img = Up(u) if s == '+' else Um(u)
            b[u] += k; b[img] -= k
        return Counter({k: v for k, v in b.items() if v})
    def valid_sheets(self):
        return all((k >= 0) if s == '+' else (k <= 0) for (u, s), k in self.c.items())

def fusion(a, b):
    r = Counter({a * b: 1, 1: 1}); r[a] -= 1; r[b] -= 1
    return Counter({k: v for k, v in r.items() if v})

def transport_defect(rec, x, xp):
    d = rec.boundary(); d[x] -= 1; d[xp] += 1
    return Counter({k: v for k, v in d.items() if v})

def multiplier_move(x, m, steps):
    """+ path of `steps` odd steps from m x, minus the - path of m to 1 (m must be - rooted)."""
    mp, rooted = minus_orbit(m)
    assert rooted, "multiplier not 3x-1 rooted"
    pp = plus_path(m * x, steps)
    r = Receipt(); r.add_plus_path(pp); r.sub_minus_path(mp)
    return r, pp[-1]

# ---------------- prefix code search ----------------
def class_certificate(s, K, m, pos, max_steps):
    """Class x = s + 2^K t (t >= 0).  Multiply by m at odd index `pos` (pos = 0: at x).  Returns (steps, ratio) if the
    affine image after `steps` odd steps is y = P x + Q with P < 1, the valuation word determined by the K known bits,
    and y < x for the smallest member (hence for all members); else None."""
    P = Fraction(1); Q = Fraction(0); A = 0
    rep = s if s > 1 else s + (1 << K)      # worst case = smallest member > 1
    y = rep
    for i in range(max_steps + 1):
        if i == pos:
            P *= m; Q *= m; y *= m
        if i > 0 and P < 1 and y < rep:
            return (i, P, Q)
        if i == max_steps:
            return None
        a = v2(3 * y + 1)
        if a >= K - A:                       # the word is not determined by the known bits
            # coarse rule: every member has valuation >= K - A here, so U(member) <= (3 member' + 1)/2^(K-A);
            # accept if that bound descends for the smallest member (affine bound with coefficient 3P/2^(K-A))
            Pb = P * 3 / (1 << (K - A)); yb = (3 * y + 1) // (1 << (K - A))
            if i + 1 <= max_steps and Pb < 1 and yb < rep:
                return (i + 1, Pb, Fraction(0))
            return None
        y = (3 * y + 1) >> a
        P = P * 3 / (1 << a); Q = (Q * 3 + 1) / (1 << a); A += a

def build_code(mults, Kmax, max_steps=20, max_pos=3):
    """DFS over odd residue classes; returns dict class(s, K) -> certificate (m, pos, steps) and the uncovered classes."""
    code = {}; uncovered = []
    stack = [(1, 2), (3, 2)]
    while stack:
        s, K = stack.pop()
        cert = None
        c = class_certificate(s, K, 1, 0, max_steps)
        if c:
            cert = (1, 0, c[0], c[1], c[2])
        else:
            for m in mults:
                for pos in range(0, max_pos + 1):
                    c = class_certificate(s, K, m, pos, max_steps)
                    if c:
                        cert = (m, pos, c[0], c[1], c[2]); break
                if cert: break
        if cert:
            code[(s, K)] = cert
        elif K < Kmax:
            stack.append((s, K + 1)); stack.append((s + (1 << K), K + 1))
        else:
            uncovered.append((s, K))
    return code, uncovered

def build_code_plain_first(mults, Kmax, max_steps=20, max_pos=3):
    """refine for plain descent up to Kmax first; multipliers only on the leaves that plain descent cannot pay."""
    code, unc = build_code([], Kmax, max_steps, max_pos)
    still = []
    for (s, K) in unc:
        cert = None
        for m in mults:
            for pos in range(0, max_pos + 1):
                c = class_certificate(s, K, m, pos, max_steps)
                if c:
                    cert = (m, pos, c[0], c[1], c[2]); break
            if cert: break
        if cert: code[(s, K)] = cert
        else: still.append((s, K))
    return code, still

def code_stats(code, uncovered):
    plain = Fraction(0); mult = Fraction(0); unc = Fraction(0); used = Counter()
    for (s, K), (m, pos, steps, P, Q) in code.items():
        d = Fraction(2, 1 << K)
        if m == 1: plain += d
        else: mult += d; used[m] += d
    for (s, K) in uncovered:
        unc += Fraction(2, 1 << K)
    return plain, mult, unc, used

def lookup(code, x):
    """deepest class containing x"""
    K = 2
    while True:
        s = x % (1 << K)
        if (s, K) in code:
            return (s, K), code[(s, K)]
        K += 1
        if K > 40: return None, None

def universal_receipt(x, code, record=False):
    """stage-by-stage descent of x to 1 with the code; returns (number of multiplier events, list of (m, x) fusions)."""
    events = []
    guard = 0
    while x > 1:
        cls, cert = lookup(code, x)
        if cert is None:
            k = v2(x + 1)
            if k >= 11:                        # deep cell -1 mod 2^k: trunk move with j = 1 mod 6, k-10 <= j <= k-5 (AL)
                j = k - 5
                while j % 6 != 1: j -= 1
                mj = ((1 << j) + 1) // 3
                y = Up(mj * x)                 # = x + (x+1)/2^j, exact depth k - j in [5, 10]
                events.append((mj, x, 'trunk'))
                cls, cert = lookup(code, y)
                if cert is None:
                    return None, events
                x = y
                # the next stage descends y below y, and y < (1 + 2^-5) x; continue
            else:
                return None, events
        m, pos, steps, P, Q = cert
        y = x
        for i in range(steps):
            if i == pos: y *= m
            y = Up(y)
        if pos == steps:                      # multiplication at the final index never happens (pos < steps)
            pass
        if m != 1:
            events.append((m, x, 'code'))
        assert y < x, (x, cert, y)
        x = y
        guard += 1
        if guard > 100000: raise RuntimeError
    return len(events), events

def main():
    Kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 65536
    out = []
    def P(*a):
        s = " ".join(str(v) for v in a); print(s, flush=True); out.append(s)
    t0 = time.time()

    # ---- A: exact checks
    P("# A. the trunk move and Proposition A (exact)")
    ok = 0; tot = 0; exact = 0; exact_tot = 0
    for k in range(3, 21):
        for j in range(1, k + 1, 2):
            mj = ((1 << j) + 1) // 3
            for t in range(1, 8):
                x = (t << k) - 1
                xp = oddpart(x + (x + 1) // (1 << j))
                r = Receipt(); r.add_plus_path([mj * x, Up(mj * x)]); r.sub_minus_path([mj, 1])
                tot += 1
                if Up(mj * x) == xp and r.value() == Fraction(x, xp) and transport_defect(r, x, xp) == fusion(mj, x) and r.valid_sheets():
                    ok += 1
                if j < v2(x + 1):
                    exact_tot += 1
                    if Up(mj * x) == x + (x + 1) // (1 << j): exact += 1
    P(f"  trunk move: U_+(m_j x) = oddpart(x + (x+1)/2^j), value x/x', transport defect = R(m_j, x), sheets valid: {ok}/{tot}; "
      f"with j < v_2(x+1) the landing point x + (x+1)/2^j is itself odd: {exact}/{exact_tot}")
    rng = random.Random(5)
    roots = [m for m in range(3, 400, 2) if minus_rooted(m)]
    ok = 0; tot = 0
    for _ in range(400):
        m = rng.choice(roots); x = rng.randrange(3, 10 ** 6, 2); d = rng.randrange(1, 12)
        r, xp = multiplier_move(x, m, d)
        tot += 1
        if r.value() == Fraction(x, xp) and transport_defect(r, x, xp) == fusion(m, x) and r.valid_sheets():
            ok += 1
    P(f"  Proposition A (general - rooted multiplier, any + depth): value x/x' and defect R(m,x): {ok}/{tot}")
    P(f"  3x-1-rooted odd m <= 100: {[m for m in range(3, 101, 2) if minus_rooted(m)]}")
    P(f"  AL multipliers and their - fate: " + ", ".join(f"{m}:{'rooted' if minus_rooted(m) else 'cycle ' + str(minus_orbit(m)[0][-1])}" for m in (5, 7, 11, 13, 23, 29, 43, 25, 35)))
    ex = Receipt(); ex.add_plus_path([21, 1]); ex.sub_minus_path([3, 1])
    P(f"  example x = 7, j = 3, m_3 = 3: 21 -> 1 on +, 3 -> 1 on -: value {ex.value()}, boundary {dict(ex.boundary())}, defect vs [7]-[1]: {dict(transport_defect(ex, 7, 1))} = R(3,7)")
    P(f"  [{time.time()-t0:.0f}s]")

    # ---- B: prefix codes
    P(f"# B. prefix codes to Kmax = {Kmax}, at most one multiplier, inserted at odd index <= 3, <= 20 odd steps")
    M_minus = [m for m in range(3, 101, 2) if minus_rooted(m)]
    M_AL = [5, 7, 11, 13, 23, 29, 43, 25, 35]
    codes = {}
    for name, M, builder in (("plain only", [], build_code), ("M_minus greedy (multiplier tried before refining)", M_minus, build_code),
                             ("M_minus plain-first", M_minus, build_code_plain_first), ("M_AL plain-first", M_AL, build_code_plain_first),
                             ("trunk only {3, 11, 43, 171, 683} plain-first", [3, 11, 43, 171, 683], build_code_plain_first)):
        code, unc = builder(M, Kmax)
        plain, mult, uncd, used = code_stats(code, unc)
        codes[name] = (code, unc)
        P(f"  {name}: classes {len(code)}; odd density plain {float(plain):.4f}, multiplier {float(mult):.4f}, uncovered {float(uncd):.5f} ({len(unc)} classes" + (f", e.g. {unc[:8]}" if unc else "") + ")")
        if used:
            P("     multipliers used (density): " + ", ".join(f"{m}:{float(d):.4f}" for m, d in sorted(used.items(), key=lambda kv: -kv[1])[:12]))
    # the AL-style 12-bit code with M_minus for direct comparison with the paper's 20.26%
    code12, unc12 = build_code_plain_first(M_minus, 12)
    plain, mult, uncd, used = code_stats(code12, unc12)
    P(f"  Kmax = 12 with M_minus: plain {float(plain):.4f}, multiplier {float(mult):.4f}, uncovered {float(uncd):.5f} ({len(unc12)} classes {unc12[:4]}); AL paper: plain 0.7969, multiplier 0.2026, uncovered -1 mod 4096 only")
    P(f"  [{time.time()-t0:.0f}s]")

    # ---- C: universal receipts and their defect rank
    P(f"# C. universal two-sheet receipts for odd x <= {N}: fusion relations per x")
    for name in ("M_minus greedy (multiplier tried before refining)", "M_minus plain-first", "M_AL plain-first", "trunk only {3, 11, 43, 171, 683} plain-first"):
        code, unc = codes[name]
        hist = Counter(); fail = 0; maxev = (0, 0); tot_ev = 0; cnt = 0; trunk_ev = 0
        for x in range(3, N + 1, 2):
            n, ev = universal_receipt(x, code)
            if n is None:
                fail += 1; continue
            hist[n] += 1; tot_ev += n; cnt += 1; trunk_ev += sum(1 for e in ev if e[2] == 'trunk')
            if n > maxev[0]: maxev = (n, x)
        P(f"  {name}: no receipt {fail}; fusion relations per x: mean {tot_ev/cnt:.3f} (of which deep-cell trunk moves {trunk_ev/cnt:.4f}), max {maxev}; distribution " + ", ".join(f"{k}:{hist[k]}" for k in sorted(hist)[:12]))
    P(f"  [{time.time()-t0:.0f}s]")

    # ---- D: the trunk move as a shortcut of x's own orbit
    P("# D. deep cells x = -1 mod 2^k (exact depth k): does the + orbit of x meet the + orbit of m_j x before x descends below itself?")
    for k in (8, 10, 12, 14, 16, 20):
      for j in (3, 5, 7, 9, 11, 13):
        if j > k - 2: continue
        mj = ((1 << j) + 1) // 3
        meet = 0; tot = 0; meet_any = 0
        for t in range(1, 1001, 2):            # x = t 2^k - 1 with t odd: exact depth k
            x = (t << k) - 1
            # orbit of x until below x
            ox = []; y = x
            while y >= x and len(ox) < 5000:
                ox.append(y); y = Up(y)
            ox_set = set(ox[1:])
            # orbit of m_j x until below x
            oy = []; z = mj * x
            while z >= x and len(oy) < 5000:
                z = Up(z); oy.append(z)
            tot += 1
            if any(w in ox_set for w in oy): meet += 1
            # meet anywhere above 1 (including after descent)
            full_x = set(plus_path(x, 3000)) - {1}
            if any(w in full_x for w in plus_path(mj * x, 3000)[1:] if w != 1): meet_any += 1
        P(f"  k={k}, j={j}, m_j={mj}: pairs {tot}; orbits meet before x descends: {meet}; meet anywhere above 1: {meet_any}")
    P(f"  [{time.time()-t0:.0f}s]")


    # ---- F: the sheet commutator and the sibling chain
    P("# F. the sheet commutator U_+(3x) = U_-(U_+(x)) for v_2(3x+1) = 1, and the sibling chain X = 2^c Y + 1")
    ok = 0; tot = 0; bad1 = 0; tot1 = 0
    for x in range(3, 1000001, 2):
        if v2(3 * x + 1) == 1:
            tot += 1
            if Up(3 * x) == Um(Up(x)): ok += 1
        else:
            tot1 += 1
            if Up(3 * x) == Um(Up(x)): bad1 += 1
    P(f"  x = 3 mod 4 (valuation one): identity holds {ok}/{tot}; x = 1 mod 4: holds only {bad1}/{tot1}")
    def fresh_pair(n):
        """sheet pair of n: big B = 2^c S + eps with {B, S} = {U_+(n), U_-(n)}."""
        X, Y = Up(n), Um(n)
        if v2(3 * n + 1) == 1:
            c = v2(3 * n - 1) - 1; assert X == (1 << c) * Y + 1; return X, Y, +1, c     # big = + image
        else:
            c = v2(3 * n + 1) - 1; assert Y == (1 << c) * X - 1; return Y, X, -1, c     # big = - image
    def chain(x, cap=10000):
        """x = 3 mod 4: start from the sheet pair of x1 = U_+(x) (the x side x2 = U_+(x1), the shadow y = U_-(x1)).
        Both sides take + steps.  With B = 2^c S + eps:
          (+1,1) -> (+1, b) with b = v2(3S+1);  (+1,2) -> MERGE (B = 4S+1);  (+1,3) -> (-1, b+1);  (+1,>=4) -> BREAK;
          (-1,1) -> fresh pair of S (U_+(B) = U_-(S));  (-1,>=2) -> BREAK.   Returns (result, steps, entry)."""
        B, S, eps, c = fresh_pair(Up(x))
        entry = f"{'P' if eps > 0 else 'M'}{min(c, 4)}"
        steps = 0
        while steps < cap:
            assert B == (1 << c) * S + eps
            if eps == 1:
                if c == 1:
                    b = v2(3 * S + 1); B, S, c = 3 * S + 2, (3 * S + 1) >> b, b
                elif c == 2:
                    return 'merge', steps, entry
                elif c == 3:
                    b = v2(3 * S + 1); B, S, eps, c = 6 * S + 1, (3 * S + 1) >> b, -1, b + 1
                else:
                    return 'break', steps, entry
            else:
                if c == 1:
                    B, S, eps, c = fresh_pair(S)
                else:
                    return 'break', steps, entry
            steps += 1
        return 'cap', steps, entry
    def truth_meet(x, cap=400):
        """ground truth: do the + orbits of x2 = U_+^2(x) and of y = U_-(U_+(x)) share an element > 1 within cap steps?"""
        x1 = Up(x); y = Um(x1); X = Up(x1)
        ox = set(); z = X
        for _ in range(cap):
            ox.add(z)
            if z == 1: break
            z = Up(z)
        z = y
        for _ in range(cap):
            if z in ox and z != 1: return True
            if z == 1: return False
            z = Up(z)
        return False
    cnt = Counter(); n3 = 0; truth = 0; agree = 0; merge_truth = 0
    for x in range(3, 1000001, 4):
        r, k, how = chain(x); cnt[(r, how)] += 1; n3 += 1
        if x <= 200001:
            tm = truth_meet(x); truth += tm
            if r == 'merge':
                merge_truth += tm
    P(f"  x = 3 mod 4, x <= 10^6 ({n3} values): by entry relation and outcome: " + ", ".join(f"{how}/{r}:{cnt[(r,how)]/n3:.4f}" for (r, how) in sorted(cnt, key=lambda t: (t[1], t[0]))))
    tot_merge = sum(v for (r, how), v in cnt.items() if r == 'merge') / n3
    P(f"    structural merge total {tot_merge:.4f} (automaton prediction: from (+1,1) 1/2, from a fresh pair 1/3); "
      f"ground truth (orbits of x2 and y share an element > 1 within 400 steps, x <= 2 10^5): {truth/50001:.4f}")
    mm = 0; mt = 0
    for x in range(3, 200001, 4):
        r, k, how = chain(x)
        if r == 'merge':
            mt += 1; mm += truth_meet(x)
    P(f"    automaton merges below 2 10^5: {mt}, confirmed by the orbits: {mm}")
    # deep cells: entry is P1 (a = 2), merge decided at step k - 3
    for k in (8, 12, 16, 20, 24):
        m = Counter(); tot = 0; lens = []
        for t in range(1, 4001, 2):
            x = (t << k) - 1
            r, kk, how = chain(x); m[r] += 1; tot += 1
            if r == 'merge': lens.append(kk)
        P(f"  deep cell -1 mod 2^{k} (exact): merge {m['merge']/tot:.4f}, break {m['break']/tot:.4f} ({tot} values; merge step always {min(lens)}..{max(lens)}, i.e. k - 3 = {k-3})")
    # sibling repair of the x3 events of the universal receipts (M_minus plain-first code)
    code, unc = codes["M_minus plain-first"]
    ev3 = 0; rep = 0; evo = 0
    for x in range(3, N + 1, 2):
        n, ev = universal_receipt(x, code)
        if n is None: continue
        for (m, xi, kind) in ev:
            if m == 3:
                ev3 += 1
                if xi % 4 == 3 and chain(xi)[0] == 'merge': rep += 1
            else:
                evo += 1
    P(f"  universal receipts x <= {N}: x3 events {ev3}, of which sibling-repairable (the x-orbit merges with the orbit of U_-(U_+(x))) {rep} ({rep/max(1,ev3):.4f}); other-multiplier events {evo}")
    P(f"  [{time.time()-t0:.0f}s]")

    # ---- E: mixed trunk
    P("# E. mixed trunk x = (4^e - 1)/(2^j + 1), j | e, j odd: m_j x = (4^e - 1)/3 is on the + trunk")
    rows = []
    for j in (3, 5, 7, 9, 11):
        cnt = 0
        for e in range(j, 60, j):
            x = ((1 << (2 * e)) - 1) // ((1 << j) + 1)
            if ((1 << (2 * e)) - 1) % ((1 << j) + 1): continue
            mj = ((1 << j) + 1) // 3
            if x % 2 == 1 and Up(mj * x) == 1 and cnt < 3:
                rows.append((j, e, x, mj, None)); cnt += 1
    P("  (j, e, x, m_j): " + ", ".join(str(r[:4]) for r in rows))
    for (j, e, x, mj, _) in rows[:6]:
        p = [x]; y = x
        while y != 1: y = Up(y); p.append(y)
        P(f"    x={x}: own + orbit length {len(p)-1} odd steps; the trunk move gives a two-sheet receipt with defect R({mj}, {x}) in one step")
    P(f"  [{time.time()-t0:.0f}s total]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
