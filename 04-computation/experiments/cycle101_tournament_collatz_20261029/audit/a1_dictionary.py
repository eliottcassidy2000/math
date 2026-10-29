#!/usr/bin/env python3
"""Independent audit of THM-4849 (A), (B), (B'), (C), (F)-traces, (G) and HYP-9446 (4).

Written from scratch for the audit; it imports nothing from the session's scripts.

Conventions checked (THM-4849 section 0):
  T(x) = x/2 (x even), (3x+1)/2 (x odd), on rationals with odd denominator (2-adic parity);
  U(x) = (3x+1)/2^v; a U-word w = (v_1..v_k) with H = sum v_i has marked point C(w)/(2^H - 3^k),
  C(w) = sum_j 3^(k-1-j) 2^(s_j), s_j = v_1+...+v_j, s_0 = 0.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import comb, gcd
import sys

def T(x):
    # parity of a rational with odd denominator = parity of numerator (denominator odd)
    assert x.denominator % 2 == 1
    if x.numerator % 2 == 0:
        return x / 2
    return (3 * x + 1) / 2

def parity(x):
    return x.numerator % 2

def point_from_parity_word(u):
    """The unique 2-adic (here rational) point whose T-parity vector is the periodic word u (independent
    derivation: x_{j+1} = (3^{u_j} x_j + u_j)/2, so 2^L x_L = 3^k x_0 + sum_j u_j 3^{#ones after j} 2^j)."""
    L = len(u)
    k = sum(u)
    c = 0
    for j in range(L):
        if u[j]:
            ones_after = sum(u[j + 1:])
            c += 3 ** ones_after * 2 ** j
    return Fr(c, 2 ** L - 3 ** k)

def C_of_word(w, a=3):
    k = len(w)
    s = 0
    tot = 0
    for j in range(k):
        tot += a ** (k - 1 - j) * 2 ** s
        s += w[j]
    return tot

def D_of_word(w, a=3):
    return 2 ** sum(w) - a ** len(w)

def U_orbit(x, steps):
    out = [x]
    for _ in range(steps):
        y = 3 * x + 1
        v = 0
        while y.numerator % 2 == 0:
            y /= 2
            v += 1
        out.append(y)
        x = y
    return out

def gaps(S, q):
    s = sorted(S)
    return tuple(((s[(i + 1) % len(s)] - s[i]) % q) or q for i in range(len(s)))

def canon_necklace(S, q):
    """translation class of a subset of Z_q: lexicographically least translate (as sorted tuple)."""
    best = None
    for t in range(q):
        c = tuple(sorted((x + t) % q for x in S))
        if best is None or c < best:
            best = c
    return best

def is_connection_set(S, q):
    S = set(S)
    if 0 in S:
        return False
    for x in range(1, q):
        if (x in S) == ((-x) % q in S):
            return False
    return True

def connection_sets(q):
    k = (q - 1) // 2
    pairs = [(x, q - x) for x in range(1, k + 1)]
    for mask in range(2 ** k):
        yield frozenset(p[0] if (mask >> i) & 1 == 0 else p[1] for i, p in enumerate(pairs))

def units(q):
    return [m for m in range(1, q) if gcd(m, q) == 1]

def rot_word(w, i):
    return tuple(w[i:] + w[:i])

def necklace_of_word(w):
    return min(rot_word(tuple(w), i) for i in range(len(w)))

def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True

LOG = []
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.append(s)

def check(cond, msg):
    if not cond:
        say("FAIL:", msg)
        raise SystemExit(1)

# ----------------------------------------------------------------------------------------------
# (0) conventions: cycle formula vs exact iteration, T-word positions s_0..s_(k-1)
def check_conventions():
    import random
    rnd = random.Random(1)
    n = 0
    for _ in range(400):
        k = rnd.randint(1, 6)
        w = tuple(rnd.randint(1, 6) for _ in range(k))
        H = sum(w)
        D = D_of_word(w)
        if D == 0:
            continue
        x = Fr(C_of_word(w), D)
        # U-orbit with the prescribed valuations must close
        y = x
        for v in w:
            z = 3 * y + 1
            check(z != 0, "zero")
            # exact valuation
            num = z.numerator
            vv = 0
            while num % 2 == 0:
                num //= 2
                vv += 1
            check(vv == v, f"valuation mismatch w={w}")
            y = z / 2 ** v
        check(y == x, f"cycle does not close w={w}")
        # T-cycle of length H with parity ones exactly at s_0..s_(k-1)
        pos = [0]
        for v in w[:-1]:
            pos.append(pos[-1] + v)
        y = x
        par = []
        for _ in range(H):
            par.append(parity(y))
            y = T(y)
        check(y == x, "T-cycle length")
        check([i for i, b in enumerate(par) if b] == pos, "T-word odd positions")
        n += 1
    say(f"(0) conventions: cycle formula C(w)/(2^H-3^k) and T-word odd positions s_0..s_(k-1) verified on {n} random words")

# ----------------------------------------------------------------------------------------------
# (A) dictionary at primes q
def check_dictionary(q):
    k = (q - 1) // 2
    Dq = 2 ** q - 3 ** k
    nsets = 0
    for S in connection_sets(q):
        nsets += 1
        u = [1 if j in S else 0 for j in range(q)]
        z = point_from_parity_word(u)
        # exact T-period q with k odd points, parity vector = u
        y = z
        pts = []
        for j in range(q):
            check(parity(y) == u[j], "parity vector of z_S")
            pts.append(y)
            y = T(y)
        check(y == z, "period")
        check(len(set(pts)) == q, "exact period q")
        odd_pts = [p for p in pts if parity(p) == 1]
        check(len(odd_pts) == k, "k odd points")
        # points are C/D_q
        for p in pts:
            check((p * Dq).denominator == 1, "denominator divides D_q")
        # U-word is the cyclic gap sequence of S
        s = sorted(S)
        first = s[0]
        x0 = pts[first]
        orb = U_orbit(x0, k)
        check(orb[k] == x0, "U-cycle closes")
        g = []
        for i in range(k):
            z1 = 3 * orb[i] + 1
            num = z1.numerator
            vv = 0
            while num % 2 == 0:
                num //= 2
                vv += 1
            g.append(vv)
        check(tuple(g) == gaps(S, q), "U-word = gap sequence")
        # vertices: z_v has parity vector 1_(v+S); T(z_v) = z_(v-1); z_v odd iff -v in S
        zv = {}
        for v in range(q):
            uv = [1 if ((j - v) % q) in S else 0 for j in range(q)]
            zv[v] = point_from_parity_word(uv)
        for v in range(q):
            check(T(zv[v]) == zv[(v - 1) % q], "T(z_v) = z_(v-1)")
            check((parity(zv[v]) == 1) == (((-v) % q) in S), "z_v odd iff v -> 0")
        # multiplier m: w_(mS)_j = w_S(m^-1 j)
        for m in range(1, q):
            mS = frozenset((m * x) % q for x in S)
            check(is_connection_set(mS, q), "mS connection set")
            mi = pow(m, -1, q)
            for j in range(q):
                check((j in mS) == (((mi * j) % q) in S), "decimation")
        # converse: -S is the reversal j -> -j
        nS = frozenset((-x) % q for x in S)
        for j in range(q):
            check((j in nS) == (((-j) % q) in S), "reversal")
    say(f"(A) q={q}: D_q = {Dq}; all {nsets} connection sets: z_S has exact T-period q, k = {k} odd points, "
        f"points in (1/D_q)Z, U-word = gap sequence, T(z_v) = z_(v-1), z_v odd iff v->0, decimation, reversal: OK")

# ----------------------------------------------------------------------------------------------
# (B) base points, counts, Catalan comparison, class/orbit tables
def base_points_of_necklace(rep, q):
    """number of translates of the subset rep that are connection sets"""
    return sum(1 for t in range(q) if is_connection_set({(x + t) % q for x in rep}, q))

def check_two_base_points(q):
    k = (q - 1) // 2
    I = frozenset(range(1, k + 1))
    rot_sets = {frozenset((m * x) % q for x in I) for m in units(q)}
    cycles = {}
    for S in connection_sets(q):
        cycles.setdefault(canon_necklace(S, q), []).append(S)
    nb_hist = {}
    for neck, Ss in cycles.items():
        nb = len(Ss)
        check(nb == base_points_of_necklace(neck, q), "base point count consistency")
        nb_hist[nb] = nb_hist.get(nb, 0) + 1
        is_rot = any(S in rot_sets for S in Ss)
        # reversal-symmetric: necklace of -S equals necklace of S
        S0 = Ss[0]
        rev_sym = canon_necklace({(-x) % q for x in S0}, q) == neck
        check((nb == 2) == is_rot, f"two base points iff rotational, q={q}")
        check((nb == 2) == rev_sym, f"two base points iff achiral, q={q}")
        check(nb <= 2, "at most two")
    ncyc = len(cycles)
    check(ncyc == 2 ** k - k, f"count 2^k-k at q={q}")
    return ncyc, nb_hist

def all_necklaces(q, k):
    seen = set()
    for c in combinations(range(q), k):
        seen.add(canon_necklace(c, q))
    return len(seen)

def catalan(k):
    return comb(2 * k, k) // (k + 1)

def class_table(q):
    """isomorphism classes = multiplier orbits on connection sets; orbit sizes on cycles"""
    k = (q - 1) // 2
    sets = list(connection_sets(q))
    seen = set()
    classes = []
    for S in sets:
        if S in seen:
            continue
        orb = {frozenset((m * x) % q for x in S) for m in range(1, q)}
        seen |= orb
        stab = sum(1 for m in range(1, q) if frozenset((m * x) % q for x in S) == S)
        cyc = {canon_necklace(T_, q) for T_ in orb}
        classes.append((len(orb), stab, len(cyc), sorted(S)))
    return classes

def ham_path_count(adj, n):
    # number of directed Hamiltonian paths (ordered vertex sequences)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for m in range(1, 1 << n):
        for v in range(n):
            c = dp[m][v]
            if c == 0:
                continue
            for w in range(n):
                if not (m >> w) & 1 and adj[v][w]:
                    dp[m | (1 << w)][w] += c
    return sum(dp[(1 << n) - 1])

def circ_adj(S, q):
    return [[1 if ((j - i) % q) in S else 0 for j in range(q)] for i in range(q)]

def tournament_invariants(S, q):
    adj = circ_adj(S, q)
    # 3-cycles, 4-vertex profile via counts of cyclic triangles through... use H and number of 3-cycles + 5-cycles
    H = ham_path_count(adj, q)
    # number of directed 4-cycles and 5-cycles through vertex 0 (cheap invariants)
    def cycles_through0(L):
        cnt = 0
        def rec(path, used):
            nonlocal cnt
            if len(path) == L:
                if adj[path[-1]][0]:
                    cnt += 1
                return
            for w in range(1, q):
                if not used[w] and adj[path[-1]][w]:
                    used[w] = True
                    path.append(w)
                    rec(path, used)
                    path.pop()
                    used[w] = False
        used = [False] * q
        used[0] = True
        rec([0], used)
        return cnt
    return (H, cycles_through0(4), cycles_through0(5))

def brute_iso_classes_q7():
    from itertools import permutations
    sets = list(connection_sets(7))
    def adjm(S):
        return [[1 if ((j - i) % 7) in S else 0 for j in range(7)] for i in range(7)]
    mats = [adjm(S) for S in sets]
    reps = []
    for idx, A in enumerate(mats):
        found = False
        for r in reps:
            B = mats[r]
            for p in permutations(range(7)):
                if all(A[i][j] == B[p[i]][p[j]] for i in range(7) for j in range(7)):
                    found = True
                    break
            if found:
                break
        if not found:
            reps.append(idx)
    return [sorted(sets[r]) for r in reps]

# ----------------------------------------------------------------------------------------------
# (B') composite q
def lexi_structure(S, t, q):
    g = gcd(t, q)
    n = q // g
    sub = {(g * j) % q for j in range(n)}  # <t> = gZ_q
    # every coset other than <t> inside S or missing S
    for a in range(1, g):
        coset = {(a + x) % q for x in sub}
        inter = coset & S
        if inter and inter != coset:
            return False
    target = {(-2 * t * j) % q for j in range(1, (n - 1) // 2 + 1)}
    return (S & sub) == target

def check_composite(q):
    k = (q - 1) // 2
    sets = list(connection_sets(q))
    setset = set(sets)
    pairs = 0
    agree = 0
    for S in sets:
        for t in range(1, q):
            St = frozenset((x + t) % q for x in S)
            a = St in setset
            b = lexi_structure(set(S), t, q)
            check(a == b, f"(B') iff fails at q={q}, S={sorted(S)}, t={t}")
            if a:
                pairs += 1
    cycles = {}
    for S in sets:
        cycles.setdefault(canon_necklace(S, q), []).append(S)
    hist = {}
    for neck, Ss in cycles.items():
        hist[len(Ss)] = hist.get(len(Ss), 0) + 1
    I = frozenset(range(1, k + 1))
    rot = {frozenset((m * x) % q for x in I) for m in units(q)}
    doubly_not_rot = sum(1 for neck, Ss in cycles.items() if len(Ss) >= 2 and not any(S in rot for S in Ss))
    # closed formula: ordered pairs = sum_{n|q, n>1} phi(n) 2^((q/n-1)/2)
    def phi(n):
        return sum(1 for m in range(1, n + 1) if gcd(m, n) == 1)
    formula = sum(phi(n) * 2 ** ((q // n - 1) // 2) for n in range(3, q + 1) if q % n == 0)
    check(pairs == formula, f"pair-count formula q={q}")
    return pairs, len(cycles), hist, doubly_not_rot, formula

# ----------------------------------------------------------------------------------------------
def main():
    check_conventions()
    for q in (3, 5, 7, 11, 13):
        check_dictionary(q)
    say("")
    say("(B) two-base-point lemma (at most two; two iff m*{1..k}; two iff reversal-symmetric) and count 2^k - k:")
    for q in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        ncyc, hist = check_two_base_points(q)
        k = (q - 1) // 2
        msg = f"   q={q}: {ncyc} tournament cycles = 2^{k}-{k}; base points per cycle {dict(sorted(hist.items()))}"
        if q <= 23:
            allc = all_necklaces(q, k)
            check(allc == catalan(k), "Cat_k necklaces")
            msg += f"; all shape-({k},{q}) necklaces = {allc} = Cat_{k}"
        say(msg)
    for k in range(1, 30):
        if catalan(k) == 2 ** k - k:
            say(f"   Cat_k = 2^k - k at k = {k}")
    check(all(catalan(k) > 2 ** k - k for k in range(4, 200)), "Catalan above for k>=4")
    say("   Cat_k > 2^k - k for 4 <= k < 200: OK")
    say("")
    say("(B) class tables (multiplier orbits):")
    for q in (5, 7, 11, 13):
        cl = class_table(q)
        say(f"   q={q}: {len(cl)} classes; (orbit size on connection sets, |stabiliser|, cycles contributed, rep) = "
            + "; ".join(f"({a},{b},{c},{d})" for a, b, c, d in cl))
    # isomorphism classes: brute force at q=7, invariants separate the orbits at q = 11, 13
    reps7 = brute_iso_classes_q7()
    say(f"   q=7 brute-force isomorphism classes over all 7! relabellings: {len(reps7)} -> {reps7}")
    for q in (11, 13):
        cl = class_table(q)
        invs = [tournament_invariants(set(r[3]), q) for r in cl]
        say(f"   q={q}: invariants (H, #4-cycles through 0, #5-cycles through 0) of the {len(cl)} multiplier orbits: {invs}; "
            f"pairwise distinct: {len(set(invs)) == len(invs)}")
    say("")
    say("(B') composite q: S + t connection set <=> lexicographic structure; base points; pair-count formula:")
    for q in (9, 15, 21, 25, 27, 33, 35, 39):
        pairs, ncyc, hist, dnr, formula = check_composite(q)
        say(f"   q={q}: ordered pairs (S,S+t) = {pairs} = formula {formula}; tournament cycles {ncyc}; "
            f"base points per cycle {dict(sorted(hist.items()))}; doubly-based not m*I: {dnr}; max base points {max(hist)}")
    say("")
    # (C) the (3,7) table
    say("(C) the (3,7) family:")
    q, k = 7, 3
    I = {1, 2, 3}
    for S in connection_sets(7):
        pass
    necks = {}
    for S in connection_sets(7):
        necks.setdefault(canon_necklace(S, 7), []).append(sorted(S))
    QR = {1, 2, 4}
    for neck, Ss in sorted(necks.items()):
        S0 = set(Ss[0])
        w = gaps(S0, 7)
        # markings: the three rotations of the word
        ms = []
        for i in range(3):
            wi = rot_word(w, i)
            C = C_of_word(wi)
            D = D_of_word(wi)
            x = Fr(C, D)
            # append the letter 1 and reclose
            wv = wi + (1,)
            Cv, Dv = C_of_word(wv), D_of_word(wv)
            check(Cv == 3 * C + 2 ** 7, "C(wv) = 3C + 2^H")
            red = Fr(Cv, Dv).denominator
            ms.append((wi, C, red))
        nums = [m[1] for m in ms]
        trace = sum(Fr(n, 101) for n in nums)
        Ttrace = None
        # T-trace: all points of the T-cycle
        x = Fr(nums[0], 101)
        y = x
        tt = Fr(0)
        for _ in range(7):
            tt += y
            y = T(y)
        check(y == x, "T-cycle")
        check(tt == 2 * trace + 3, "T-trace = 2 U-trace + k")
        rev = canon_necklace({(-a) % 7 for a in S0}, 7)
        dec2 = canon_necklace({(2 * a) % 7 for a in S0}, 7)
        dec3 = canon_necklace({(3 * a) % 7 for a in S0}, 7)
        letters = w
        locked5 = all((v - 1) % 4 == 0 for v in letters)
        say(f"   necklace {neck}: connection sets {Ss}; word {w}; markings (word, C, reduced denom after appending 1): {ms}; "
            f"U-trace {trace}; T-trace {tt}; reversal -> {rev}; x2 -> {dec2}; x3 -> {dec3}; "
            f"letters = 1 mod 4: {locked5}")
    # the v=1 fibre
    fib = {}
    for w in [(a, b, 7 - a - b) for a in range(1, 6) for b in range(1, 7 - a) if 7 - a - b >= 1]:
        d = Fr(C_of_word(w + (1,)), D_of_word(w + (1,))).denominator
        fib.setdefault(d, []).append(C_of_word(w))
    say(f"   v=1 fibre over all 15 marked words (reduced denominator: numerators C): {dict(sorted(fib.items()))}")
    m5 = sorted(C_of_word(w) for w in [(a, b, 7 - a - b) for a in range(1, 6) for b in range(1, 7 - a) if 7 - a - b >= 1]
                if (Fr(C_of_word(w), 101) + 1).numerator % 5 == 0)
    m25 = [c for c in m5 if (Fr(c, 101) + 1).numerator % 25 == 0]
    m7 = sorted(C_of_word(w) for w in [(a, b, 7 - a - b) for a in range(1, 6) for b in range(1, 7 - a) if 7 - a - b >= 1]
                if (Fr(C_of_word(w), 101) + 1).numerator % 7 == 0)
    say(f"   markings = -1 mod 5: {m5}; also mod 25: {m25}; = -1 mod 7: {m7}")
    say("")
    # (G) 5x+1
    say("(G) 5x+1 cycles and the all-odd criterion on the (3,7) shape:")
    for start in (1, 13, 17):
        x = start
        orb = [x]
        wd = []
        while True:
            y = 5 * x + 1
            v = 0
            while y % 2 == 0:
                y //= 2
                v += 1
            wd.append(v)
            x = y
            if x == start:
                break
            orb.append(x)
        say(f"   5x+1 cycle from {start}: {orb}, word {tuple(wd)}")
    for neck, Ss in sorted(necks.items()):
        w = gaps(set(Ss[0]), 7)
        res = []
        for i in range(3):
            wi = rot_word(w, i)
            c5 = C_of_word(wi, 5)
            d5 = D_of_word(wi, 5)
            res.append((wi, c5, c5 % 3, Fr(c5, d5)))
        allodd = all(v % 2 == 1 for v in w)
        integral = all(r[3].denominator == 1 for r in res)
        check(allodd == integral, "all-odd criterion")
        say(f"   necklace {neck} word {w}: (marking, C_5, C_5 mod 3, point) {res}; all letters odd {allodd}; integral {integral}")
    # interval marking C(1^(k-1), k+2) = 3^k - 2^k
    for kk in range(1, 15):
        w = (1,) * (kk - 1) + (kk + 2,)
        check(C_of_word(w) == 3 ** kk - 2 ** kk, "interval marking")
    say("   C(1^(k-1), k+2) = 3^k - 2^k for k = 1..14: OK; C(1,5,1) =", C_of_word((1, 5, 1)), ", C(5,1,1) =", C_of_word((5, 1, 1)))

if __name__ == "__main__":
    main()
