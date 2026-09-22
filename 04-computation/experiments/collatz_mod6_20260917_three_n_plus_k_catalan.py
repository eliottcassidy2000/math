#!/usr/bin/env python3
"""
collatz_mod6_20260917_three_n_plus_k_catalan.py
(session collatz-mod6-20260917, LANE = three_n_plus_k_catalan)

Object: the odd-only maps  T_k(n) = (3n+k)/2^{v_2(3n+k)}  on the odd integers of
both signs, k odd, 3 does not divide k.

Sections
  S1  PROVED identities, checked exactly: negation conjugacy T_{-k}(-n)=-T_k(n);
      scaling T_{dk}(dn)=d T_k(n) (d odd); gcd(n,k) is an orbit invariant;
      translation identity T_k = T_{+-1} o (n -> n + (k -+ 1)/3);
      inverse-fibre braid R_k(n)=4n+k with v_3(R_k^t(n)-n)=v_3(t).
  S2  PROVED necklace mechanism for k = 2^K-3^L: rotation identity
      3B(w)+k = 2^{k_1} B(w'), B(w') odd; injectivity of B on words; necklace
      count = cyclic classes of compositions; periodic words = scaled cycles.
      Verified on k = 13, 5, -1, 7, -11, 1.
  S3  Catalan typing: 2^K-3^L = +-1 iff (K,L) in {(1,1),(2,1),(3,2)} (elementary
      proof in the note; finite check here); 3n+1 has one necklace cycle, 3n-1
      has two, and its 7-cycle is sporadic with 139 | B(1,1,1,2,1,1,4).
  S4  FINITE-EXACT census: all odd k, 3 not | k, |k| <= 199, all odd starts
      |n| <= 10^6, escape bound 10^18, step cap 10^5; every cycle typed by
      (word, K, L, D=2^K-3^L, g=gcd(n0,k), m=Dg/k) as necklace / anti-necklace /
      scaled / sporadic.  Cross-checked against F3 and the scaling law.
  S5  3n-5 in detail; rows mod 6 / mod 9 of 3n-5 versus 3n+1; R_k braid law.

All load-bearing arithmetic is exact Python int.  All checks use explicit
`raise` (active under python -O).  RAM << 1 GB, runtime ~ 2-3 minutes.

Inheritance (read, not re-derived):
  05-knowledge/results/arithmetic_braids_20260917_collatz.md   (rows 6j+1->9j+2; R(n)=4n+1;
      cycle gate n = B/(2^K-3^L) for k=1, eq. (B7); three 3n-1 cycles)
  05-knowledge/results/arithmetic_braids_20260917_summand.md   (signed conjugation C_+(-n)=-C_-(n))
  05-knowledge/results/arithmetic_braids_20260917_divisors.md  (F=S+U iff p, p^3, p^2qr)
  05-knowledge/results/collatz_mod6_20260917_row_braid_typing.out (F_p tower, R period 3 mod 6)
  Session-lead facts F1-F4 (F3 = cycle census |n|<=2*10^5, k=+-1,+-5,+-7,+-11,+-13).
"""
import sys, time, array, itertools
from math import gcd, comb
from sympy import totient, mobius, divisors, factorint

T0 = time.time()

def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)

def hr(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)

def v2(y):
    if y == 0:
        raise ValueError("v2(0)")
    return (y & -y).bit_length() - 1

def T(k, n):
    y = 3 * n + k
    return y // (y & -y)

def word_B(w):
    """B(w) = sum_{i=0}^{L-1} 3^{L-1-i} 2^{K_i}, K_0=0, K_i=k_1+..+k_i."""
    L = len(w)
    B = 0
    Ki = 0
    for i in range(L):
        B += 3 ** (L - 1 - i) * (1 << Ki)
        Ki += w[i]
    return B

def compositions(K, L):
    """All compositions of K into L positive parts."""
    if L == 1:
        yield (K,)
        return
    for first in range(1, K - L + 2):
        for rest in compositions(K - first, L - 1):
            yield (first,) + rest

def canon_rot(w):
    L = len(w)
    return min(tuple(w[i:] + w[:i]) for i in range(L))

def necklace_count(K, L):
    """# cyclic classes of compositions of K into L parts (Burnside)."""
    g = gcd(K, L)
    return sum(int(totient(d)) * comb(K // d - 1, L // d - 1) for d in divisors(g)) // L

def lyndon_count(K, L):
    g = gcd(K, L)
    return sum(int(mobius(d)) * comb(K // d - 1, L // d - 1) for d in divisors(g)) // L

def reps(k, Kmax=80, Lmax=50):
    """All (K,L), 1<=K<=Kmax, 1<=L<=Lmax with 2^K-3^L = k."""
    out = []
    for L in range(1, Lmax + 1):
        t = k + 3 ** L
        if t >= 2 and (t & (t - 1)) == 0:
            K = t.bit_length() - 1
            if 1 <= K <= Kmax:
                out.append((K, L))
    return sorted(out)

# ----------------------------------------------------------------------------
hr("S0  Definitions, universe, inheritance")
print("T_k(n) = (3n+k)/2^{v_2(3n+k)},  k odd, 3 not | k, n odd integer of either sign.")
print("Word of a cycle (n_0,..,n_{L-1}): w=(k_1..k_L), k_i = v_2(3 n_{i-1}+k);  K = sum w;  D = 2^K-3^L.")
print("B(w) = sum_{i=0}^{L-1} 3^{L-1-i} 2^{K_i} (K_0=0).  Cycle gate (PROVED below): n_0 * D = k * B(w).")
print("Census universe: k in {odd, 3 not | k, |k|<=199} (134 values), odd starts |n|<=10^6, escape 10^18, cap 10^5.")
print("Inherited: cycle gate (B7) for k=1 and the three 3n-1 cycles [arithmetic_braids_20260917_collatz.md];")
print("  signed conjugation [arithmetic_braids_20260917_summand.md]; F3 census |n|<=2*10^5 (session lead).")

# ----------------------------------------------------------------------------
hr("S1  PROVED identities (checked exactly on a finite box)")
KS = [k for k in range(-199, 200) if k % 2 and k % 3]
check(len(KS) == 134, "134 admissible k")
box = [n for n in range(-3001, 3002, 2)]
# (a) negation conjugacy
for k in KS:
    for n in box:
        if T(-k, -n) != -T(k, n):
            raise RuntimeError("negation conjugacy fails")
print("(S1.1) PROVED  T_{-k}(-n) = -T_k(n)   [v_2(-y)=v_2(y)]; checked all 134 k, odd |n|<=3001.")
# (b) scaling
for k in KS:
    for d in (3, 5, 7, 9, 15, -1, -5):
        for n in box[::7]:
            if T(d * k, d * n) != d * T(k, n):
                raise RuntimeError("scaling fails")
print("(S1.2) PROVED  T_{dk}(dn) = d T_k(n) for odd d (v_2(dy)=v_2(y)); checked d in {3,5,7,9,15,-1,-5}.")
# (c) gcd invariance
for k in KS:
    for n in box[::3]:
        if gcd(abs(T(k, n)), abs(k)) != gcd(abs(n), abs(k)):
            raise RuntimeError("gcd invariance fails")
print("(S1.3) PROVED  gcd(T_k(n),k) = gcd(n,k):  gcd((3n+k)/2^v, k) = gcd(3n+k,k) = gcd(3n,k) = gcd(n,k) (k odd, 3 not | k).")
print("       Hence cycles(T_k) = disjoint union over positive g | k of  g * (primitive cycles of T_{k/g}),")
print("       primitive := gcd(n_0,k)=1.  Every cycle element is odd and not divisible by 3 (T_k(n) = k*2^{-v} mod 3 != 0).")
# (d) translation identity
for k in KS:
    if k % 6 == 1:
        s = (k - 1) // 3
        for n in box[::5]:
            if T(k, n) != T(1, n + s):
                raise RuntimeError("translation +1 fails")
    else:
        s = (k + 1) // 3
        for n in box[::5]:
            if T(k, n) != T(-1, n + s):
                raise RuntimeError("translation -1 fails")
print("(S1.4) PROVED  k = 1 mod 6: T_k(n) = T_1(n + (k-1)/3);   k = 5 mod 6: T_k(n) = T_{-1}(n + (k+1)/3)")
print("       (3n+k = 3(n+(k-1)/3)+1, and (k-1)/3 is even exactly when k = 1 mod 6).  Every 3n+k is 3n+-1 with a")
print("       shifted input; T_k is conjugate (by the shift) to n -> T_{+-1}(n) + shift.  Checked all 134 k.")
# (e) braid law
for k in KS[::9]:
    for n in box[::11]:
        y = 3 * n + k
        for t in range(1, 28):
            m = n
            for _ in range(t):
                m = 4 * m + k
            diff = m - n
            if diff != (4 ** t - 1) * y // 3:
                raise RuntimeError("braid closed form")
            v3 = 0
            while diff % 3 == 0:
                diff //= 3; v3 += 1
            tt = t; vt = 0
            while tt % 3 == 0:
                tt //= 3; vt += 1
            if v3 != vt:
                raise RuntimeError("braid valuation")
            if T(k, m) != T(k, n) or v2(3 * m + k) != v2(y) + 2 * t:
                raise RuntimeError("braid fibre")
print("(S1.5) PROVED  R_k(n)=4n+k: 3R_k(n)+k = 4(3n+k), so T_k(R_k n)=T_k(n) with valuation +2;")
print("       R_k^t(n)-n = (4^t-1)(3n+k)/3 and v_3(R_k^t(n)-n) = v_3(t) (3 not | 3n+k; v_3(4^t-1)=1+v_3(t)).")
print("       Same 3-adic odometer law as R_1 (inherited (B3)); checked t<=27 on a sample of k,n.")

# ----------------------------------------------------------------------------
hr("S2  PROVED necklace mechanism for k = 2^K - 3^L")
# rotation identity, exact, all compositions K<=16
cnt = 0
for K in range(1, 17):
    for L in range(1, K + 1):
        k = (1 << K) - 3 ** L
        for w in compositions(K, L):
            B = word_B(w)
            w2 = w[1:] + w[:1]
            B2 = word_B(w2)
            if B % 2 == 0 or 3 * B + k != (1 << w[0]) * B2 or T(k, B) != B2:
                raise RuntimeError("rotation identity fails at %s" % (w,))
            cnt += 1
print("(S2.1) PROVED  3B(w) + (2^K-3^L) = 2^{k_1} B(w'), w' = left rotation, B(w') odd; hence T_k(B(w)) = B(w')")
print("       with valuation exactly k_1, and B(w) lies on a cycle of T_k following the word w.")
print("       Proof: 3B(w) = 3^L + sum_{i>=1} 3^{L-i} 2^{K_i}; add k; 2^{k_1}B(w') = sum_{i=1}^{L} 3^{L-i}2^{K_i}.")
print("       Checked exactly on all %d compositions with K<=16." % cnt)
# injectivity of B on words with fixed (K,L)
for K in range(1, 15):
    for L in range(1, K + 1):
        seen = {}
        for w in compositions(K, L):
            B = word_B(w)
            if B in seen:
                raise RuntimeError("B not injective")
            seen[B] = w
print("(S2.2) PROVED  B is injective on words of fixed (K,L): the valuation sequence of the T_k-orbit of B(w) is")
print("       w itself (by S2.1), so B(w) determines w.  Hence the cycle through B(w) is exactly {B(rot^j w)},")
print("       its length is the number of distinct rotations of w (= L for aperiodic w, = L/e for w = u^e),")
print("       B-values of distinct necklaces lie on DISTINCT cycles, and distinct (K,L) with the same k give")
print("       distinct cycles (the cycle determines L: if L<L' then 2^K-3^L would divide 2^{K'}-3^{L'}=k with")
print("       cofactor >1 in absolute value).  Checked injectivity K<=14.")
# necklace counts vs enumeration
for K in range(1, 13):
    for L in range(1, K + 1):
        classes = {canon_rot(w) for w in compositions(K, L)}
        lyn = {c for c in classes if len({tuple(c[i:] + c[:i]) for i in range(L)}) == L}
        if len(classes) != necklace_count(K, L) or len(lyn) != lyndon_count(K, L):
            raise RuntimeError("necklace count")
print("(S2.3) FINITE-EXACT  #necklaces(K,L) = (1/L) sum_{d|gcd(K,L)} phi(d) C(K/d-1,L/d-1) and")
print("       #aperiodic = (1/L) sum mu(d) C(K/d-1,L/d-1) agree with enumeration for K<=12.")
print("(S2.4) PROVED  periodic word w = u^e (u of length p=L/e, sum K/e): B_L(u^e) = B_p(u) * (2^K-3^L)/(2^{K/e}-3^p),")
print("       so this necklace cycle is the (k/d)-scaled copy, d = 2^{K/e}-3^p | k, of the necklace cycle of T_d")
print("       (geometric-sum identity; the cofactor (k/d) is odd since both k and d are odd).")
for K in range(2, 17):
    for L in range(2, K + 1):
        for e in divisors(gcd(K, L)):
            if e == 1:
                continue
            p = L // e
            for u in compositions(K // e, p):
                d = (1 << (K // e)) - 3 ** p
                k = (1 << K) - 3 ** L
                if word_B(u * e) * d != word_B(u) * k:
                    raise RuntimeError("periodic scaling")
print("       Checked exactly for all periodic words with K<=16.")

def necklace_cycles(k):
    """All necklace cycles of T_k (n_0 = B(w), reps of k) and anti-necklace cycles (n_0 = -B(w), reps of -k).
    Returns list of (canonical cycle tuple, K, L, word, sign, periodic?)."""
    out = []
    for sign in (1, -1):
        for (K, L) in reps(sign * k):
            for c in sorted({canon_rot(w) for w in compositions(K, L)}):
                cyc = []
                n = sign * word_B(c)
                for _ in range(L):
                    cyc.append(n)
                    n = T(k, n)
                if n != cyc[0]:
                    raise RuntimeError("necklace cycle does not close")
                per = len(set(cyc))
                cyc = cyc[:per]
                mn = min(cyc); r = cyc.index(mn)
                cyc = tuple(cyc[r:] + cyc[:r])
                out.append((cyc, K, L, c, sign, per < L))
    return out

print()
print("Verification on the requested k (necklace cycles n_0=B(w) from k=2^K-3^L; anti-necklace n_0=-B(w) from -k=2^K-3^L):")
for k in (13, 5, -1, 7, -11, 1):
    R = reps(k); Rm = reps(-k)
    print("  k=%4d: reps(k)=%s  reps(-k)=%s  necklace counts %s / %s" % (
        k, R, Rm, [necklace_count(K, L) for (K, L) in R], [necklace_count(K, L) for (K, L) in Rm]))
    for (cyc, K, L, w, sign, per) in necklace_cycles(k):
        print("     (K,L)=(%d,%d) sign=%+d word=%s %s -> cycle %s" % (K, L, sign, w, "PERIODIC(scaled)" if per else "aperiodic", cyc))
print("k=13: 2^8-3^5=13 gives C(7,4)=35 compositions / 5 = 7 necklaces (gcd(8,5)=1): the seven 5-cycles of F3;")
print("      2^4-3=13 gives {1}.  k=-1: (1,1)->{1}, (3,2)->{5,7}.  k=5: (3,1)->{1}, (5,3)->2 three-cycles of F3.")

# ----------------------------------------------------------------------------
hr("S3  Catalan typing of 3n+1 and 3n-1")
sols = [(K, L) for K in range(1, 201) for L in range(1, 131) if abs((1 << K) - 3 ** L) == 1]
check(sols == [(1, 1), (2, 1), (3, 2)], "Catalan +-1 solutions")
print("FINITE-EXACT (K<=200,L<=130) and PROVED (note): 2^K-3^L = +-1 iff (K,L) in {(1,1),(2,1),(3,2)}:")
print("  2^K-3^L=+1: K<=2 forces (2,1); K>=3 needs 3^L = -1 mod 8, but 3^L in {1,3} mod 8.")
print("  2^K-3^L=-1: K=1 gives (1,1); K>=2 forces 3^L = 1 mod 4 so L even, then 2^K=(3^{L/2}-1)(3^{L/2}+1)")
print("     is a product of two powers of two differing by 2, i.e. 2*4, so L=2, K=3.")
print("  CITED (special case of Mihailescu 2004 / classical): the only consecutive perfect powers are 8,9.")
print("Consequences (PROVED via S2):")
print("  3n+1 (k=+1): reps(1)={(2,1)} -> one necklace cycle {1}; reps(-1)={(1,1),(3,2)} -> anti-necklaces {-1},{-5,-7}.")
print("  3n-1 (k=-1): reps(-1)={(1,1),(3,2)} -> necklace cycles {1} and {5,7}; reps(1)={(2,1)} -> anti-necklace {-1}.")
w7 = (1, 1, 1, 2, 1, 1, 4)
B7 = word_B(w7)
D7 = (1 << 11) - 3 ** 7
check(D7 == -139 and B7 == 2363 and B7 % 139 == 0 and B7 // 139 == 17, "139 | B(w7)")
cyc = [17]
for _ in range(7):
    cyc.append(T(-1, cyc[-1]))
check(cyc == [17, 25, 37, 55, 41, 61, 91, 17], "7-cycle of 3n-1")
check([v2(3 * c - 1) for c in cyc[:7]] == list(w7), "word of the 7-cycle")
print("  The 7-cycle 17->25->37->55->41->61->91 of 3n-1 has word %s, K=11, L=7, D=2^11-3^7=%d, B(w)=%d = 139*17." % (w7, D7, B7))
print("  Gate n_0*D = k*B:  17*(-139) = (-1)*2363 = -2363.  D = k*m with m=139 = -D/k, m | B(w): SPORADIC (m != +-1).")
print("PROVED typing:  the 3 positive cycles of 3n-1 = 2 Catalan-necklace cycles ({1} from 2-3=-1, {5,7} from 8-9=-1)")
print("                + 1 sporadic cycle (m=139).  Equivalently the 3 negative cycles of 3n+1 are 2 Catalan anti-necklaces + 1 sporadic.")
print("  Necessity direction (PROVED): a cycle of T_{-1} with n_0 = B(w) (m=1) needs 2^K-3^L = -1, so only (1,1),(3,2):")
print("  NO further necklace cycles exist for k=+-1; any further cycle of 3n+-1 must be sporadic (|m|>1).  Completeness OPEN.")

# ----------------------------------------------------------------------------
hr("S4  FINITE-EXACT census: all odd k, 3 not | k, |k|<=199; odd starts |n|<=10^6")
NMAX = 10 ** 6
ESC = 10 ** 18
CAP = 10 ** 5

def find_cycles(k, N=NMAX):
    UN = -2; ESCAPE = -1
    size = N + 1
    cid = array.array('i', [UN]) * (2 * size)
    cycles = []
    stats = {'escape': 0, 'cap': 0, 'steps': 0, 'starts': 0}
    for a in range(1, N + 1, 2):
        for n in (a, -a):
            stats['starts'] += 1
            path = []; pset = set(); x = n; c = None
            while True:
                if -N <= x <= N:
                    i = x if x > 0 else size - x
                    cc = cid[i]
                    if cc != UN:
                        c = cc; break
                if x in pset:
                    j = path.index(x); cyc = path[j:]
                    mn = min(cyc); r = cyc.index(mn); cyc = tuple(cyc[r:] + cyc[:r])
                    cycles.append(cyc); c = len(cycles) - 1; break
                if abs(x) > ESC:
                    c = ESCAPE; stats['escape'] += 1; break
                if len(path) > CAP:
                    c = ESCAPE; stats['cap'] += 1; break
                pset.add(x); path.append(x)
                y = 3 * x + k; x = y // (y & -y)
            stats['steps'] += len(path)
            for v in path:
                if -N <= v <= N:
                    i = v if v > 0 else size - v
                    cid[i] = c
    return cycles, stats

def type_cycle(k, cyc):
    """Return dict with word, K, L, D, B, g, m, type string."""
    L = len(cyc)
    w = tuple(v2(3 * c + k) for c in cyc)
    K = sum(w)
    D = (1 << K) - 3 ** L
    B = word_B(w)
    n0 = cyc[0]
    if n0 * D != k * B:
        raise RuntimeError("cycle gate fails for k=%d cycle %s" % (k, cyc))
    g = gcd(abs(n0), abs(k))
    if (D * g) % k:
        raise RuntimeError("m not integral")
    m = D * g // k
    if B % m or n0 != g * (B // m):
        raise RuntimeError("m | B / n0 = gB/m fails")
    if m == 1:
        base = "necklace(K=%d,L=%d)" % (K, L)
    elif m == -1:
        base = "anti-necklace(K=%d,L=%d)" % (K, L)
    else:
        base = "SPORADIC(K=%d,L=%d,D=%d,m=%d)" % (K, L, D, m)
    typ = base if g == 1 else ("%d x [T_%d: %s]" % (g, k // g, base))
    return dict(word=w, K=K, L=L, D=D, B=B, g=g, m=m, type=typ, prim=(g == 1), neck=(abs(m) == 1))

census = {}
tstart = time.time()
for k in KS:
    cyc, st = find_cycles(k)
    if st['escape'] or st['cap']:
        raise RuntimeError("escape/cap hit for k=%d: %s" % (k, st))
    census[k] = sorted(cyc, key=lambda c: (len(c), c))
print("census done: %d values of k, %d odd starts each, no escapes, no step-cap hits, %.1f s" % (
    len(KS), 2 * (NMAX // 2 + (NMAX % 2)), time.time() - tstart))

# every necklace/anti-necklace cycle must have been found
for k in KS:
    found = set(census[k])
    for (cyc, K, L, w, sign, per) in necklace_cycles(k):
        if cyc not in found:
            raise RuntimeError("necklace cycle missed by census: k=%d %s" % (k, cyc))
print("(S4.1) every necklace / anti-necklace cycle predicted by S2 was met by the census (all 134 k).")

# negation conjugacy of the census
for k in KS:
    neg = sorted((tuple(sorted([-c for c in cyc], key=lambda x: x)) for cyc in census[k]))
    pos = sorted((tuple(sorted(cyc)) for cyc in census[-k]))
    if neg != pos:
        raise RuntimeError("census not negation-symmetric at k=%d" % k)
print("(S4.2) census(-k) = -census(k) as sets, all 67 pairs (PROVED by S1.1; here checked on the finite universe).")

rows = {}
for k in KS:
    typed = [type_cycle(k, c) for c in census[k]]
    c_all = len(typed)
    c_prim = sum(t['prim'] for t in typed)
    n_scaled = c_all - c_prim
    R = reps(k); Rm = reps(-k)
    N_all = sum(necklace_count(K, L) for (K, L) in R + Rm)
    N_lyn = sum(lyndon_count(K, L) for (K, L) in R + Rm)
    N_prim = sum(1 for t in typed if t['prim'] and t['neck'])
    s_prim = c_prim - N_prim
    n_Dk = sum(1 for t in typed if abs(t['D']) == abs(k))
    if n_Dk != N_all:
        raise RuntimeError("#cycles with D=+-k != necklace count at k=%d (%d vs %d)" % (k, n_Dk, N_all))
    rows[k] = dict(c=c_all, prim=c_prim, scaled=n_scaled, N_all=N_all, N_lyn=N_lyn, N_prim=N_prim,
                   s_prim=s_prim, reps=R + Rm, typed=typed)
# scaling-law consistency: c(k) = sum_{g | k, g>0} c_prim(k/g)
for k in KS:
    s = sum(rows[k // g]['prim'] for g in divisors(abs(k)))
    if s != rows[k]['c']:
        raise RuntimeError("scaling consistency at k=%d" % k)
print("(S4.3) c(k) = sum_{g|k} c_prim(k/g) holds for all 134 k (scaled cycles are exactly the g-multiples).")
# non-primitive necklace cycles with D=k but g>1 (necklace of T_k that is a scaled sporadic of T_{k/g})
odd_neck = [(k, t) for k in KS for t in rows[k]['typed'] if abs(t['D']) == abs(k) and t['g'] > 1]
print("(S4.4) necklace-mechanism cycles (D=+-k) with gcd(B(w),k)=g>1 (they are g x sporadic/necklace cycles of T_{k/g}): %d" % len(odd_neck))
for (k, t) in odd_neck:
    print("       k=%d  word=%s  B=%d  g=%d  -> %s" % (k, t['word'], t['B'], t['g'], t['type']))

# F3 cross-check
F3 = {1: 4, -1: 4, 5: 9, -5: 9, 7: 5, -7: 5, 11: 7, -11: 7, 13: 13, -13: 13}
for k, c in F3.items():
    if rows[k]['c'] != c:
        raise RuntimeError("F3 mismatch at k=%d: %d vs %d" % (k, rows[k]['c'], c))
print("(S4.5) F3 cross-check: c(+-1)=4, c(+-5)=9, c(+-7)=5, c(+-11)=7, c(+-13)=13 reproduced.")
minel13 = sorted(min(c) for c in census[13] if len(c) == 5 and min(c) > 0)
check(minel13 == [211, 227, 251, 259, 283, 287, 319], "F3 5-cycle minima")
check(sorted(c for c in census[5] if len(c) == 3) == [(19, 31, 49), (23, 37, 29)], "F3 3-cycles of k=5")
check((5, 7) in census[-1], "F3 (5,7)")
print("       seven 5-cycle minima of k=13 = %s; k=5 3-cycles (19,31,49),(23,37,29); k=-1 has (5,7)." % minel13)

print()
print("TABLE S4-A  (k>0 only; k<0 is the negative by S1.1)   c=#cycles, prim=#primitive (gcd(n0,k)=1), scaled=c-prim,")
print("  N_all = #necklaces over reps of k and -k (= #cycles with D=+-k), N_lyn = aperiodic ones, N_prim = aperiodic with gcd(B,k)=1,")
print("  s_prim = prim - N_prim = #primitive SPORADIC cycles (|m|>1).  reps = (K,L) with 2^K-3^L = +-k.")
print("   k |  c prim scaled | N_all N_lyn N_prim | s_prim | reps (k) ; reps (-k)")
for k in KS:
    if k < 0:
        continue
    r = rows[k]
    flag = "  <-- SPORADIC" if r['s_prim'] > 0 else ""
    print(" %3d | %2d  %2d   %2d   |  %2d    %2d    %2d  |   %2d   | %s ; %s%s" % (
        k, r['c'], r['prim'], r['scaled'], r['N_all'], r['N_lyn'], r['N_prim'], r['s_prim'], reps(k), reps(-k), flag))

print()
print("TABLE S4-B  all primitive SPORADIC cycles (|m|>1), k>0 (negate for -k): k, L, K, D=2^K-3^L, m=D/k, n_0, min element, cycle sign")
tot_spor = 0
for k in KS:
    if k < 0:
        continue
    for t in rows[k]['typed']:
        if t['prim'] and not t['neck']:
            tot_spor += 1
            cyc = census[k][rows[k]['typed'].index(t)]
            print("   k=%3d  L=%2d  K=%2d  D=%7d  m=%6d  B=%10d  n0=%8d  %s  word=%s" % (
                k, t['L'], t['K'], t['D'], t['m'], t['B'], cyc[0], "neg" if cyc[0] < 0 else "pos",
                t['word'] if t['L'] <= 20 else str(t['word'][:20]) + "..."))
print("   total primitive sporadic cycles for k>0: %d;  k values with at least one: %s" % (
    tot_spor, [k for k in KS if k > 0 and rows[k]['s_prim'] > 0]))
print("   k values (k>0) with NO primitive cycle at all beyond necklaces and with N_prim=0: %s" % (
    [k for k in KS if k > 0 and rows[k]['prim'] == 0]))

print()
top = sorted((k for k in KS if k > 0), key=lambda k: (-rows[k]['c'], k))[:12]
print("TABLE S4-C  the k>0 with the most cycles and where they come from")
print("   k |  c | prim | scaled | N_all | s_prim | #divisors g of k | factorization | why")
for k in top:
    r = rows[k]
    why = []
    if r['N_all'] >= 2:
        why.append("%d necklace cycles from reps %s" % (r['N_all'], r['reps']))
    if r['scaled'] >= 2:
        why.append("%d scaled from divisors %s" % (r['scaled'], divisors(k)[1:]))
    if r['s_prim'] >= 1:
        why.append("%d sporadic" % r['s_prim'])
    print(" %3d | %2d |  %2d  |   %2d   |  %2d   |   %2d   | %d | %s | %s" % (
        k, r['c'], r['prim'], r['scaled'], r['N_all'], r['s_prim'], len(divisors(k)), dict(factorint(k)), "; ".join(why)))

# stability probe: bigger universe for the top k and the requested k
print()
print("Stability probe: re-run the census with |n| <= 4*10^6 for k in {1,5,7,11,13} u top-4:")
for k in sorted(set([1, 5, 7, 11, 13] + top[:4])):
    cyc4, st = find_cycles(k, 4 * 10 ** 6)
    if st['escape'] or st['cap']:
        raise RuntimeError("escape at 4e6")
    same = sorted(cyc4, key=lambda c: (len(c), c)) == census[k]
    print("   k=%3d: c=%d at 10^6, c=%d at 4*10^6, identical cycle sets: %s" % (k, rows[k]['c'], len(cyc4), same))
    if not same:
        print("      NEW cycles at 4*10^6: %s" % [c for c in cyc4 if c not in set(census[k])])

# ----------------------------------------------------------------------------
hr("S5  3n-5 in detail; rows mod 6 / mod 9; the braid R_k")
k = -5
print("Divisors of |k|=5: g in {1,5}.  cycles(T_{-5}) = prim(T_{-5}) u 5*prim(T_{-1});  prim(T_{-1}) = all 4 cycles of 3n-1.")
print("reps(-5) = %s (necklace, n_0=+B);  reps(+5) = %s (anti-necklace, n_0=-B)." % (reps(-5), reps(5)))
print("  (2,2): compositions of 2 into 2 parts: (1,1) only, PERIODIC = (1)^2: B=3+2=5, cycle {5} = 5 x {1}, {1} the necklace")
print("         of T_{-1} from 2-3=-1.  So the (2,2) necklace IS the scaled cycle 5x{1}: counted once, as scaled (S2.4).")
print("  (3,1): word (3): B=1 -> anti-necklace {-1}:  3(-1)-5 = -8 -> -1.")
print("  (5,3): two aperiodic necklaces (1,1,3),(1,2,2) -> anti-necklaces {-19,-31,-49}, {-23,-29,-37}.")
print()
print("  All 9 cycles of 3n-5 (canonical rotation from the minimum), typed:")
for cyc in census[-5]:
    t = type_cycle(-5, cyc)
    print("   %-95s L=%2d K=%2d D=%6d g=%d m=%5d  %s" % (str(cyc) if len(cyc) <= 7 else str(cyc[:5])[:-1] + ", ...)",
                                                            t['L'], t['K'], t['D'], t['g'], t['m'], t['type']))
r5 = rows[-5]
print("  Count: 9 = %d scaled (5x{1},5x{5,7},5x{17..91},5x{-1}) + %d primitive necklace/anti-necklace ({-1},{-19,..},{-23,..})"
      " + %d primitive sporadic (two negative 17-cycles, D=2^27-3^17 = %d = -5*%d)." % (
          r5['scaled'], r5['N_prim'], r5['s_prim'], (1 << 27) - 3 ** 17, ((1 << 27) - 3 ** 17) // -5))
for cyc in census[-5]:
    if len(cyc) == 17:
        t = type_cycle(-5, cyc)
        check(t['D'] == (1 << 27) - 3 ** 17 and t['m'] == -t['D'] // 5 and t['B'] % t['m'] == 0, "17-cycle typing")
print("  Both 17-cycles have K=27 (2^27=134217728 vs 3^17=129140163, ratio 1.0393): word sum 27, so each is a 'near-Catalan'")
print("  pair (27,17) with D=%d = -5 x %d, and %d | B(w) in both cases." % ((1 << 27) - 3 ** 17, -((1 << 27) - 3 ** 17) // 5, -((1 << 27) - 3 ** 17) // 5))

print()
print("Rows mod 6 -> mod 9 of the single-halving maps F_k(n)=(3n+k)/2, n odd:")
for kk in (1, -5, 7, 13, -11):
    imgs = []
    for r in (1, 3, 5):
        vals = {((3 * (6 * j + r) + kk) // 2) % 9 for j in range(0, 60)}
        if len(vals) != 1:
            raise RuntimeError("row image not single class")
        imgs.append(vals.pop())
    print("   k=%3d: 6j+1 -> %d mod 9,  6j+3 -> %d mod 9,  6j+5 -> %d mod 9   (shift vs k=1: %d)" % (
        kk, imgs[0], imgs[1], imgs[2], (imgs[0] - 2) % 9))
print("PROVED: F_k(6j+r) = 9j + (3r+k)/2, so the image row of source row r is (3r+k)/2 mod 9; the three images are")
print("  always {2,5,8} mod 9 (a permutation of the k=1 rows), shifted by (k-1)/2 mod 9: k=-5 shifts by -3 = 6,")
print("  i.e. rows (1,3,5) -> (8,2,5) instead of (2,5,8): F_{-5}(n) = F_1(n) - 3 = F_1(n-2) (translation S1.4).")
print("  Powers of two in the image rows: 2^e = 2 mod 3 needs e odd, and the row is 2^e mod 9 in {2,8,5} for")
print("  e = 1,3,5 mod 6 -- identical for every k = 1 mod 6 up to the row relabelling.")
print("Inverse-fibre braid R_k(n) = 4n+k: PROVED in S1.5 to have the SAME 3-adic period law v_3(R_k^t n - n)=v_3(t) as R_1;")
print("  in particular R_k has period exactly 3^s on the odd classes mod 2*3^s, and period 3 on rows mod 6, for every k.")
print("  Mod 6: R_k adds 4n+k-n = 3n+k = 3+k mod 6 to an odd n: k=1 -> +4 (1->5->3->1), k=-5=1 mod 6 -> +4 too.")

print()
print("Elapsed %.1f s.  No check failed." % (time.time() - T0))
