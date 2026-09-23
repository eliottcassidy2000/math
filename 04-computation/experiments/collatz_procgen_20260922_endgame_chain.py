#!/usr/bin/env python3
"""collatz_procgen_20260922_endgame_chain.py -- exact chain calculus for the Q2 endgame (backward E-game).

Reverse move x -> (2^k x - 1)/3 (legal iff the result is a 3-adic unit; for positive integers: a positive
integer prime to 3).  A path of s moves and K halvings has multiplier R = 2^K/3^s and endpoint R(x - beta).

Canonical escape Psi (note, section 1).  Every 3-adic unit x lies in the class of exactly one of the two
primary hostile points, h(x) = 1 if x = 1 (mod 3) and h(x) = 1/2 if x = 2 (mod 3).  With k = v_3(x - h),
    Psi(x) = rho_h(k) (x - h),  rho_1(k) = 2^{K*(k-1)}/3^k,  rho_{1/2}(k) = 2^{K*(k-1)+1}/3^k,
K*(0)=0, K*(1)=2, K*(s)=K0(s)=floor((s+1) log2 3) for s>=2 (HYP-9122; FINITE-EXACT for s<=6000).
Psi is realised by an explicit legal path of exactly k moves (a minimal loop of length k-1 through 1, first
move +1 for h=1/2, then the move 0).

Parts (all exact; every check is an assert):
  A  multiplier table, the budget constant c* = ln(128/81)/4 (max over all branches of ln(rho)/k)
  B  explicit execution of the escapes with loop words (s <= SMAX) on random integers; digit consumption
  C  optimality inside k moves (brute force, small k)
  D  the frontier recursion w' = (2^(K_t+3) w - 1)/3^j' versus the optimal w' = (2^(K0(k-1)+1) w - 1)/3^k'
  E  the universal landing recursion u_{t+1} = (2^{a_t} u_t - h_{t+1})/3^{k_{t+1}} on random Psi-orbits
  F  thread table: Psi-itineraries of the 7 known threads and the 98-point backward census; the Psi-lift
     lemma checked on random integers near every thread; per-digit link prices
Usage: python3 collatz_procgen_20260922_endgame_chain.py [SMAX] [path to bwd_bad_r41.txt (optional census check)]
"""
import sys, math, random
from fractions import Fraction as F

LN2, LN3 = math.log(2.0), math.log(3.0)
HALF = F(1, 2)

def K0(s): return (3 ** (s + 1)).bit_length() - 1
def Kstar(s): return 0 if s == 0 else (2 if s == 1 else K0(s))
def c_of(s): return F(2 ** K0(s), 3 ** s)
def rho(half, k): return F(2 ** (Kstar(k - 1) + (1 if half else 0)), 3 ** k)
def v3(x):
    x = F(x); n, d = x.numerator, x.denominator; v = 0
    if n == 0: return 10 ** 9
    while n % 3 == 0: n //= 3; v += 1
    while d % 3 == 0: d //= 3; v -= 1
    return v
def res3(x, mod):
    x = F(x); return (x.numerator * pow(x.denominator, -1, mod)) % mod
def psi(x):
    x = F(x); h = F(1) if res3(x, 3) == 1 else HALF
    if x == h: return None
    k = v3(x - h); r = rho(h == HALF, k)
    return (r * (x - h), h, k, r)
def psi_int(m):
    t, half = (m - 1, False) if m % 3 == 1 else (2 * m - 1, True)
    k = 0
    while t % 3 == 0: t //= 3; k += 1
    return t << Kstar(k - 1), half, k
def rev_int(x, k):
    t = (x << k) - 1
    assert t % 3 == 0, (x, k)
    y = t // 3
    assert y >= 1 and y % 3 != 0, ("illegal", x, k, y)
    return y
def fmt(x):
    x = F(x); return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"
def itinerary(x, maxsteps=200):
    x = F(x); P = F(1); D = 0; steps = []
    for _ in range(maxsteps):
        r = psi(x)
        if r is None: return steps, ('1' if x == 1 else '1/2')
        y, h, k, rr = r; P *= rr; D += k
        steps.append((h, k, rr, P, D, y)); x = y
    return steps, None

def find_loops(smax, vmax):
    """minimal loops through 1 (layered DP over values <= vmax)."""
    INF = 10 ** 9; dist = {1: 0}; par = []; out = {0: (0, [])}
    for s in range(1, smax + 1):
        nd = {}; pp = {}
        for v, Kv in dist.items():
            k = 0 if v % 3 == 1 else 1; p = (1 << k) * v
            while p <= 3 * vmax + 1:
                if p % 9 in (4, 7):
                    y = (p - 1) // 3
                    if 1 <= y <= vmax and Kv + k < nd.get(y, INF): nd[y] = Kv + k; pp[y] = (v, k)
                k += 2; p *= 4
        par.append(pp); dist = nd
        if 1 in dist:
            ks = []; v = 1
            for t in range(s - 1, -1, -1):
                u, k = par[t][v]; ks.append(k); v = u
            ks.reverse(); out[s] = (dist[1], ks)
    return out
def escape_word(half, k, loops):
    K, ks = loops[k - 1]; ks = list(ks)
    assert K == Kstar(k - 1), (k, K)
    if half:
        if k == 1: return [1]
        ks[0] += 1
    return ks + [0]

# The 98-point backward census (dyadic points p/2^e, e <= 40, value in [0.3,2], lying in a class of the
# dimension lane's Bad_41; reconstructed from its dump bwd_bad_r41.txt).  Stored as data; section F
# recomputes it from the dump when the dump path is given as the 2nd argument.
CENSUS = """1/2 1580017/2097152 832185350369/1099511627776 6283/8192 3309745115/4294967296 13171577/16777216
52451/65536 27640222387/34359738368 209/256 110155585/134217728 439291/524288
231582132299/274877906944 1753/2048 924291401/1073741824 3691475/4194304 14753/16384
7781751697/8589934592 59/64 31126123/33554432 124585/131072 65740797977/68719476736 499/512
263357891/268435456 1 1055729/1048576 557307443425/549755813888 4235/4096 2236003291/2147483648
8977273/8388608 36067/32768 19050287795/17179869184 145/128 4754357/4194304 76601153/67108864
9958534033/8589934592 19127/16384 10074987409/8589934592 308219/262144 10106274631/8589934592
39629179/33554432 162862655563/137438953472 10230631825/8589934592 40084075/33554432
10361641873/8589934592 83155056665/68719476736 1241/1024 40692061/33554432 159577/131072
655855945/536870912 84086683673/68719476736 84336981449/68719476736 41203819/33554432
331382339/268435456 85331839001/68719476736 41742955/33554432 335021507/268435456 163951/131072
86379919385/68719476736 2642899/2097152 86661504383/68719476736 42463531/33554432
339885395/268435456 166057/131072 1397044423331/1099511627776 87484069913/68719476736
87780719129/68719476736 343979459/268435456 661/512 88959809561/68719476736 348292547/268435456
10657/8192 349451333/268435456 171241/131072 5634268049/4294967296 90201978905/68719476736
90535709273/68719476736 354057155/268435456 91862186009/68719476736 358909379/268435456 43/32
177073/131072 22737515/16777216 93259626521/68719476736 365394563/268435456 93635073185/68719476736
370853315/268435456 715/512 91817/65536 376604099/268435456 378149147/268435456
48560928793/34359738368 384290243/268435456 371/256 196249027/134217728 793585/524288
1675672562339/1099511627776 419868489953/274877906944 1690578594467/1099511627776"""

def census_points(dump=None):
    pts = sorted(F(t) for t in CENSUS.split())
    if dump:
        Mb = 41; MODB = 3 ** (Mb + 1)
        bb = [int(l.split()[0]) for l in open(dump)]
        rec = set()
        for e in range(0, 41):
            q = 2 ** e
            for c in bb:
                p = (c * q) % MODB
                if q * 3 // 10 <= p <= 2 * q and (p % 2 == 1 or e == 0): rec.add(F(p, q))
        assert sorted(rec) == pts, "census mismatch with the dump"
        print(f"  census recomputed from {dump}: identical ({len(pts)} points)")
    return pts

def main():
    SMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    dump = sys.argv[2] if len(sys.argv) > 2 else None
    random.seed(20260922)
    print("=" * 100)
    print("A. Psi multipliers rho_h(k) and the budget constant c* = max ln(rho)/k")
    print("   k   rho_1(k)           rho_1/2(k)          ln(rho_1/2)/k")
    best = (-1, None)
    for k in range(1, 13):
        r1, r2 = rho(False, k), rho(True, k)
        val = math.log(r2) / k
        if val > best[0]: best = (val, k)
        print(f"  {k:2d}   {fmt(r1):>14s}={float(r1):.4f}  {fmt(r2):>14s}={float(r2):.4f}   {val:+.6f}")
    for k in range(13, 6002):
        assert Kstar(k - 1) + 1 - k * math.log2(3) < 1          # rho_1/2 < 2 for k >= 3
        assert math.log(2) / k < best[0]
    for k in range(3, 6002):
        assert rho(False, k) < 1 < rho(True, k) < 2            # 1-escape descends, 1/2-transfer costs in (1,2)
    assert best[1] == 4 and abs(best[0] - math.log(128 / 81) / 4) < 1e-15
    print(f"  c* = max over all branches = ln(128/81)/4 = {best[0]:.9f}, attained only at (h=1/2, k=4);"
          f" rho_1(k) < 1 < rho_1/2(k) < 2 for all 3 <= k <= 6001")

    print("=" * 100)
    print(f"B. explicit escapes with minimal loops (layered DP, s <= {SMAX})")
    loops = find_loops(SMAX, 5000 if SMAX <= 80 else 60000)
    bad = [s for s in range(2, SMAX + 1) if s not in loops or loops[s][0] != K0(s)]
    assert not bad, bad
    print(f"  loops with K = K0(s) found for all 2 <= s <= {SMAX} (HYP-9122 instances), K*(1)=2 (trivial loop)")
    n_exec = 0
    for k in range(1, SMAX + 2):
        for half in (False, True):
            word = escape_word(half, k, loops)
            assert len(word) == k and sum(word) == Kstar(k - 1) + (1 if half else 0)
            for _ in range(6):
                u = random.randrange(1, 10 ** 25)
                if u % 3 == 0: u += 1
                if half:
                    if u % 2 == 0: u += 3                      # need (1 + 3^k w)/2 integral: w odd
                    m = (1 + 3 ** k * u) // 2
                else:
                    m = 1 + 3 ** k * u
                y = m
                for kk in word: y = rev_int(y, kk)
                yy, hh, kv = psi_int(m)
                assert kv == k and hh == half and y == yy
                assert y == (u << Kstar(k - 1)), "landing formula"
                n_exec += 1
    print(f"  {n_exec} escapes executed move by move on random integers (legality asserted): endpoint = Psi(m)"
          f" = 2^K* u exactly, for every k <= {SMAX+1}")
    # digit consumption: on the shell v_3(x-h)=k, x mod 3^(n+k) <-> Psi(x) mod 3^n is a bijection onto units
    for k in range(1, 6):
        for half in (False, True):
            h = HALF if half else F(1)
            for n in range(1, 5):
                img = []
                for r in range(3 ** (n + k)):
                    if r % 3 == 0: continue
                    t = (2 * r - 1) if half else (r - 1)          # r - h up to the unit factor 1/2
                    t %= 3 ** (n + k)
                    if t % 3 ** k != 0 or t % 3 ** (k + 1) == 0: continue   # exact valuation k (known since n >= 1)
                    img.append(res3(rho(half, k) * (F(r) - h), 3 ** n))
                units = [a for a in range(3 ** n) if a % 3]
                assert sorted(img) == units, (k, half, n)
    print("  digit consumption: on each shell {v_3(x-h)=k}, x mod 3^(n+k) determines Psi(x) mod 3^n, uniformly onto"
          " the units (checked k<=5, n<=4): every Psi step consumes exactly k ternary digits")

    print("=" * 100)
    print("C. optimality inside k moves (inherited Prop 2.1 / Lemma 7.2), brute force over all legal paths")
    def paths_min(m, L, kcap):
        best = None; stack = [(m, 0)]
        while stack:
            x, s = stack.pop()
            if s >= 1 and (best is None or x < best[0]): best = (x, s)
            if s == L: continue
            k = 0 if x % 3 == 1 else 1; p = (1 << k) * x
            while k <= kcap:
                if p % 9 in (4, 7): stack.append(((p - 1) // 3, s + 1))
                k += 2; p *= 4
        return best
    nchk = 0
    for k in range(3, 7):
        for half in (False, True):
            for _ in range(20):
                u = random.randrange(1, 10 ** 6)
                if u % 3 == 0: u += 1
                if half and u % 2 == 0: u += 3
                m = (1 + 3 ** k * u) // 2 if half else 1 + 3 ** k * u
                # every endpoint of a path of length <= k is increasing in each move exponent, so moves <= k*2+4 suffice
                bst = paths_min(m, k, 2 * k + 6)
                y, _, _ = psi_int(m)
                if half:
                    assert bst[0] == y and y > m          # least endpoint within k moves = Psi(m) > m
                else:
                    assert bst[0] == y and y < m          # the 1-escape is the least endpoint within k moves
                nchk += 1
    print(f"  {nchk} random shells (k=3..6): least endpoint over all legal paths of length <= k equals Psi(m);"
          f" > m on the 1/2 class, < m on the 1 class")

    print("=" * 100)
    print("D. the frontier-item-1 recursion versus the optimal transfer")
    diff = []; agree = []
    for k in range(4, SMAX + 2):
        Kfront = Kstar(k - 2) + 2          # route 1/2 -(3)-> 1 -> loop of length k-2 -> 0: lands at 2^(K*(k-2)+2) w
        Kopt = Kstar(k - 1)
        assert Kfront >= Kopt
        (agree if Kfront == Kopt else diff).append(k)
        # execute the frontier route on a random integer and check its landing
        w = random.randrange(1, 10 ** 20) * 6 + 1
        m = (1 + 3 ** k * w) // 2
        y = rev_int(m, 3)
        for kk in loops[k - 2][1]: y = rev_int(y, kk)
        y = rev_int(y, 0)
        assert y == w << Kfront
        # the recursion w' = (2^(K_t+3) w - 1)/3^j'  with K_t = K*(k-2)
        t = (y << 1) - 1; j2 = 0
        while t % 3 == 0: t //= 3; j2 += 1
        assert t == ((w << (Kstar(k - 2) + 3)) - 1) // 3 ** j2
    frac = len(diff) / (len(diff) + len(agree))
    print(f"  the frontier route (3, loop(k-2), 0) is legal and lands at 2^(K*(k-2)+2) w (checked, k=4..{SMAX+1});")
    print(f"  it equals the optimal landing 2^(K0(k-1)) w iff K0(k-1) = K0(k-2)+2; it is twice as large for"
          f" {len(diff)} of {len(diff)+len(agree)} k (fraction {frac:.3f}; limit 2 - log2 3 = {2-math.log2(3):.3f}),"
          f" e.g. k = {diff[:8]}")
    print("  optimal recursion on the 1/2 thread:  w' = (2^(K0(k-1)+1) w - 1)/3^k'   (frontier: 2^(K_t+3), K_t=K0(k-2))")

    print("=" * 100)
    print("E. universal landing recursion on random Psi-orbits (big integers)")
    nE = 0
    for _ in range(300):
        m = random.randrange(10 ** 30, 10 ** 60)
        if m % 3 == 0: m += 1
        x = m
        hs = []; ks = []; us = []
        for t in range(12):
            y, half, k = psi_int(x)
            h = HALF if half else F(1)
            u = (F(x) - h) / 3 ** k
            hs.append(h); ks.append(k); us.append(u)
            x = y
        for t in range(11):
            a = Kstar(ks[t] - 1) + (1 if hs[t] == HALF else 0)
            assert us[t + 1] == (2 ** a * us[t] - hs[t + 1]) / 3 ** ks[t + 1]
            nE += 1
    print(f"  {nE} consecutive steps: u_(t+1) = (2^(a_t) u_t - h_(t+1)) / 3^(k_(t+1)) with a_t = K*(k_t - 1) + [h_t = 1/2]")

    print("=" * 100)
    print("F. thread table: Psi-itineraries of the known threads and of the 98-point census")
    pts = census_points(dump)
    known = [F(1, 2), F(43, 32), F(59, 64), F(145, 128), F(209, 256), F(371, 256), F(499, 512)]
    rows = []
    maxperdigit = 0
    for x in pts:
        st, term = itinerary(x)
        assert term is not None
        P = st[-1][3] if st else F(1); D = st[-1][4] if st else 0
        assert all(s[3] > 1 for s in st)                     # every prefix multiplier > 1
        word = " ".join(("1/2" if s[0] == HALF else "1") + f"@{s[1]}" for s in st)
        # link price at precision j > D: P * rho_term(j - D); per-digit sup over j
        sup = 0
        for j in range(D + 1, D + 400):
            price = P * rho(term == '1/2', j - D)
            sup = max(sup, math.log(price) / j)
        maxperdigit = max(maxperdigit, sup)
        rows.append((x, word, D, P, term, sup))
    print(f"  all {len(pts)} census points reach 1 or 1/2 exactly, every prefix multiplier > 1;"
          f" terminals: 1 -> {sum(1 for r in rows if r[4]=='1')}, 1/2 -> {sum(1 for r in rows if r[4]=='1/2')};"
          f" max itinerary length {max(len(r[1].split()) for r in rows)} steps, max digits {max(r[2] for r in rows)}")
    print(f"  max over threads and precisions of ln(link price)/digits = {maxperdigit:.6f} <= c* = {math.log(128/81)/4:.6f}")
    print("  known threads (value, itinerary [class@precision], digits D_h, multiplier P_h, terminal):")
    for x in known:
        r = [r for r in rows if r[0] == x][0]
        print(f"    {fmt(x):>8s} = {float(x):.4f}: {r[1] or '(base point)':>14s}  D={r[2]} P={fmt(r[3])} -> {r[4]}")
    # Psi-lift lemma on random integers near every thread
    nlift = 0
    for x in pts:
        if x == 1: continue
        st, term = itinerary(x); D = st[-1][4] if st else 0
        p, q = x.numerator, x.denominator; e = q.bit_length() - 1
        for j in range(D + 1, D + 25):
            for _ in range(3):
                # m integer with v_3(m - x) = j exactly: 2^e m = p + 3^j W, W = -p 3^-j mod 2^e, 3 !| W
                W0 = (-p * pow(3 ** j, -1, 2 ** e)) % (2 ** e) if e else 0
                W = W0 + (2 ** e) * random.randrange(1, 10 ** 12)
                while W % 3 == 0: W += 2 ** e
                m = (p + 3 ** j * W) >> e
                assert (p + 3 ** j * W) % (2 ** e) == 0 and v3(F(m) - x) == j
                y = m
                for s in st:                                  # the itinerary steps are followed exactly
                    yy, half, k = psi_int(y)
                    assert (HALF if half else F(1)) == s[0] and k == s[1]
                    y = yy
                tau = F(1) if term == '1' else HALF
                assert v3(F(y) - tau) == j - D                 # landing on the terminal thread at precision j - D
                # landing parameter: y - tau = 3^(j-D) * 2^A * W / 2^e * ...  with A = log2 of P_h * 3^D
                A = (st[-1][3] * 3 ** D).numerator.bit_length() - 1 if st else 0
                assert F(y) - tau == F(3 ** (j - D) * (2 ** A) * W, 2 ** e)
                y2, half2, k2 = psi_int(y)                     # the terminal escape at precision j - D
                assert k2 == j - D and half2 == (term == '1/2')
                a_t = Kstar(j - D - 1) + (1 if term == '1/2' else 0)
                assert F(y2) == F(2 ** (A + a_t) * W, 2 ** e)  # link landing  y' = 2^(A_h + a_tau(j-D_h)) W / 2^e
                nlift += 1
    print(f"  Psi-lift lemma: {nlift} random integers m with v_3(m - h) = j > D_h follow h's itinerary exactly and land at"
          f" tau + 3^(j-D_h) 2^(A_h) W / 2^e  (A_h = log2(P_h 3^D_h))")
    print("  link landing y' = 2^(A_h + a_tau(j - D_h)) W / 2^e exactly (a_tau(k) = K*(k-1) + [tau = 1/2]); link multiplier"
          " 2^(A_h + a_tau(j-D_h)) / 3^j")
    print("  => a repeated landing at a thread h'' = p''/2^e'' with precision j'' has W'' = (2^(A_h + a_tau(j-D_h) + e'' - e) W - p'') / 3^j''")

if __name__ == '__main__':
    main()
