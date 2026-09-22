#!/usr/bin/env python3
"""collatz_mod6_20260917_extended_collatz_scc_audit_proof-audit.py

INDEPENDENT proof-audit recomputation for lane extended_collatz_scc
(session collatz-mod6-20260917).  Imports nothing from the explorer's script.

E  : n -> n/2 (n even), n -> 3n+1 (all n).   E_- : 3n-1.   E_p : pn+1.

Audit items (numbered A1..A9):
  A1  leaf identity / transient forest (C1, C2), incl. 'at most v_2 halvings'
  A2  drift table, hostile 3-adic family, carry formula (4.1) of the note
      (the note's exponent of 3 is reversed; corrected form verified)
  A3  Q2 greedy to 10^6, Q1 to 10^6 (C5) -- own implementation
  A4  SCC data by an independent iterative Kosaraju (C6)
  A5  word determination, Markov chain, E[2^K_J], f_J, g_J, tail bound (C7,C8)
  A6  cycle census: own DFS (nodes<=2000,len<=40) AND a global word-enumeration
      census for every length <= 24 (no node bound) -- tests 'unique 16-cycle'
      and 'h = ceil(a log2 3)' beyond the explorer's universe (C9)
  A7  signed side (C10): leaf identity, {4,11}, k-word of the rescue path,
      Q2_-/Q1_- to 10^6, note's explicit walks (repetition!), Kosaraju, census
  A8  pn+1 hierarchy (C11): closed forms, general supercriticality p>=5,
      Chernoff certificate, E_5/E_7 controls
  A9  misc note statements: 31 -> 47 (not 121), (1,4,2) parent pair
All checks use explicit raise (survive python -O).
"""
import sys
import time
import hashlib
from fractions import Fraction
from collections import deque

T0 = time.time()


def chk(c, msg):
    if not c:
        raise AssertionError(msg)


def v2(n):
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    return c


def banner(s):
    print()
    print("#" * 76)
    print(s)
    print("#" * 76)


# ---------------------------------------------------------------------------
banner("A1  leaf identity (C1) and transient forest (C2)")
# ---------------------------------------------------------------------------
V = 10 ** 5
cnt_new = 0
for v in range(1, V + 1):
    # predecessors in E: 2v always; n with 3n+1=v iff v = 1 mod 3 and v >= 4
    p3 = (v - 1) // 3 if (v % 3 == 1 and v >= 4) else None
    new = p3 is not None and p3 % 2 == 0     # even -> 3n+1 is NOT a Collatz arrow
    if new:
        cnt_new += 1
        j = (v - 1) // 6
        chk(v % 6 == 1 and v >= 7 and p3 == 2 * j, "C1 target/pred at v=%d" % v)
        n0 = (4 * v - 1) // 3
        chk((4 * v - 1) % 3 == 0 and n0 % 2 == 1 and n0 % 8 == 1, "C1 n0 at v=%d" % v)
        chk(3 * n0 + 1 == 4 * v and (n0 - 1) // 4 == p3 and (n0 - 1) % 4 == 0, "C1 R^{-1} at v=%d" % v)
        # least odd T-predecessor: no odd n with oddpart(3n+1)=v and 3n+1 < 4v
        for n in range(1, n0, 2):
            x = 3 * n + 1
            while x % 2 == 0:
                x //= 2
            chk(x != v, "odd T-predecessor smaller than n0 at v=%d" % v)
        # (B1) at index -2 not integral
        chk((v - 4) % 12 != 0, "index -2 integral at v=%d" % v)
    else:
        chk(not (v % 6 == 1 and v >= 7), "missing new arrow at v=%d" % v)
        if v % 6 in (3, 5):
            chk(p3 is None, "3/5 mod 6 target with 3n+1 predecessor")
        if v % 6 == 4:
            chk(p3 is not None and p3 % 2 == 1, "4 mod 6 target: predecessor odd (ordinary Collatz)")
chk(cnt_new == len(range(7, V + 1, 6)), "count of new arrows")
print("C1 CONFIRMED on v <= %d: %d new arrows = #{v = 1 mod 6, 7<=v<=%d}; n0=(4v-1)/3 is the least odd"
      % (V, cnt_new, V))
print("   T-predecessor (brute-forced), n0 = 1 mod 8, 3n0+1 = 4v, (n0-1)/4 = (v-1)/3; index -2 never integral.")
# transient forest: for every n in 3Z, every maximal E-path leaves 3Z after AT MOST v_2(n) halvings
for n in range(3, 3 * 10 ** 4, 3):
    chk((3 * n + 1) % 3 == 1 and (3 * n - 1) % 3 == 2, "3n+-1 in 3Z?!")
    x = n
    halvings = 0
    while x % 2 == 0:
        x //= 2
        halvings += 1
        chk(x % 3 == 0, "halving left 3Z?!")
    chk(halvings == v2(n), "v2")
    # from x (odd multiple of 3) the only arrow is 3x+1, outside 3Z
    chk((3 * x + 1) % 3 != 0, "?")
# an E-path can leave 3Z EARLIER than v_2(n) halvings: e.g. 12 -> 37 immediately
chk((3 * 12 + 1) % 3 == 1 and v2(12) == 2, "?")
print("C2 CONFIRMED (n <= 3*10^4): no arrow enters 3Z from outside; inside 3Z only halvings; singleton SCCs.")
print("   WORDING: a path from 3a leaves 3Z after AT MOST v_2(3a) halvings (12 -> 37 leaves after 0),")
print("   the note's 'after v_2(3a) halvings' should read 'after at most v_2(3a) halvings'.")

# ---------------------------------------------------------------------------
banner("A2  drift table (C3), 3-adic hostile family (C4), carry formula (4.1)")
# ---------------------------------------------------------------------------
U9 = (1, 2, 4, 5, 7, 8)
KMIN = {}
KCL = {}
NXT = {}
for r in U9:
    ks = sorted(k for k in range(6) if (pow(2, k, 9) * r) % 9 in (4, 7))
    KCL[r] = ks
    KMIN[r] = ks[0]
    NXT[r] = ((pow(2, ks[0], 9) * r - 1) // 3) % 3
chk(KMIN == {1: 2, 2: 1, 4: 0, 5: 3, 7: 0, 8: 1}, "k_min table")
chk(KCL == {1: [2, 4], 2: [1, 3], 4: [0, 2], 5: [3, 5], 7: [0, 4], 8: [1, 5]}, "admissible classes")
chk(NXT == {1: 1, 2: 1, 4: 1, 5: 1, 7: 2, 8: 2}, "result mod 3")
for r in U9:
    chk(all(k % 2 == (1 if r % 3 == 2 else 0) for k in KCL[r]), "k parity")
# brute: for m<=2000, 3 not| m, k<=12: (2^k m-1)/3 integral & not 0 mod 3  iff 2^k m in {4,7} mod 9
for m in range(1, 2001):
    if m % 3 == 0:
        continue
    for k in range(13):
        y = (m << k) - 1
        ok = (y % 3 == 0) and ((y // 3) % 3 != 0)
        chk(ok == (((m << k) % 9) in (4, 7)), "admissibility criterion m=%d k=%d" % (m, k))
        chk(ok == ((k % 6) in KCL[m % 9]), "class criterion")
print("C3 CONFIRMED: drift table, parity rule, admissible classes, result classes (brute m<=2000, k<=12).")


def G(x):
    k = KMIN[x % 9]
    return ((x << k) - 1) // 3, k


# hostile family: v_3(m-1)=j exactly (m = 1 mod 3^j, m != 1 mod 3^(j+1)), several m per j
for j in range(1, 9):
    for m in (3 ** j + 1, 2 * 3 ** j + 1, 3 ** j * 4 + 1, 3 ** j * 7 + 1):
        if (m - 1) % 3 ** (j + 1) == 0:
            continue
        x = m
        word = []
        for _ in range(j):
            x, k = G(x)
            word.append(k)
        chk(word == [2] * (j - 1) + [0], "word at j=%d m=%d: %s" % (j, m, word))
        chk(x * 3 ** j == 4 ** (j - 1) * (m - 1), "closed form j=%d m=%d" % (j, m))
        chk((x < m) == (j <= 4) or (j >= 5 and m * (4 ** (j - 1) - 3 ** j) > 4 ** (j - 1)), "descent iff")
        if j <= 4:
            chk(x < m, "j<=4 should descend")
        else:
            chk(x > m, "j>=5 should grow (m=%d)" % m)
chk(all(4 ** (j - 1) < 3 ** j for j in range(1, 5)) and all(4 ** (j - 1) > 3 ** j for j in range(5, 40)), "4^(j-1) vs 3^j")
path = [244]
for _ in range(5):
    path.append(G(path[-1])[0])
chk(path == [244, 325, 433, 577, 769, 256], "244 family")
chk(G(256)[0] == 85 and 577 * 4 == 2308, "244 continues to 85 with E-node peak 2308")
print("C4 CONFIRMED: word 2^(j-1) 0, image 4^(j-1)(m-1)/3^j, threshold j>=5; 244->...->256->85, E-node peak 2308.")

# carry formula: note (4.1) claims B_i = sum_{l=0}^{i-1} 3^{i-1-l} 2^{K_i-K_{l+1}}.
# Derivation: the '-1' at step l (1-based) is divided by 3 at step l and multiplied by 2^{k} at
# each later step, then everything is scaled by 3^i:  contributes 3^{l-1} 2^{K_i-K_l}.
def carry_note(ks):
    i = len(ks)
    K = [0]
    for k in ks:
        K.append(K[-1] + k)
    return sum(3 ** (i - 1 - l) * 2 ** (K[i] - K[l + 1]) for l in range(i))


def carry_true(ks):
    i = len(ks)
    K = [0]
    for k in ks:
        K.append(K[-1] + k)
    return sum(3 ** (l - 1) * 2 ** (K[i] - K[l]) for l in range(1, i + 1))


bad_note = 0
for m in range(1, 3000):
    if m % 3 == 0:
        continue
    x = m
    ks = []
    for i in range(1, 7):
        x, k = G(x)
        ks.append(k)
        K = sum(ks)
        chk(3 ** i * x == 2 ** K * m - carry_true(ks), "corrected carry formula m=%d i=%d" % (m, i))
        if 3 ** i * x != 2 ** K * m - carry_note(ks):
            bad_note += 1
chk(bad_note > 0, "the note's (4.1) would have to fail somewhere")
print("NOTE (4.1) REFUTED as written: B_i = sum_{l=0}^{i-1} 3^{i-1-l} 2^{K_i-K_{l+1}} fails in %d of the" % bad_note)
print("   (m<=3000, i<=6) instances; e.g. word (2,0): note gives 3*2^0+1=4, true B_2 = 2^0+3 = 4 (coincidence),")
print("   word (0,2): note gives 3*4+1=13, true B_2 = 4+3 = 7.  Correct form: B_i = sum_{l=1}^{i} 3^{l-1} 2^{K_i-K_l}.")
chk(carry_note([0, 2]) == 13 and carry_true([0, 2]) == 7, "example")
print("   (Positivity B_i > 0, the only property used in Theorems 4.3/4.4, is unaffected.)")

# ---------------------------------------------------------------------------
banner("A3  Q2 greedy to 10^6 and Q1 to 10^6 (C5)")
# ---------------------------------------------------------------------------
NQ = 10 ** 6
t = time.time()
steps_to_1 = [0] * (NQ + 1)
arrows_to_1 = [0] * (NQ + 1)
peak_to_1 = [0] * (NQ + 1)
peak_to_1[1] = 1
stat = {"below": 0, "one": 0}
worst = {"peak": (1, 1), "steps": (0, 1), "arrows": (0, 1)}
for m in range(2, NQ + 1):
    if m % 3 == 0:
        continue
    x = m
    s = a = 0
    pk = m
    while x >= m:
        k = KMIN[x % 9]
        y = x << k
        if y > pk:
            pk = y
        x = (y - 1) // 3
        s += 1
        a += k + 1
        chk(s < 10 ** 5, "cap at m=%d" % m)
        if x == 1:
            break
        chk(x != m, "greedy cycle at m=%d" % m)
    chk(x < m and x % 3 != 0 and x >= 1, "descent target at m=%d" % m)
    stat["one" if (x == 1 and s == 1) else "below"] += 1
    steps_to_1[m] = s + steps_to_1[x]
    arrows_to_1[m] = a + arrows_to_1[x]
    peak_to_1[m] = max(pk, peak_to_1[x])
    if peak_to_1[m] > worst["peak"][0]:
        worst["peak"] = (peak_to_1[m], m)
    if steps_to_1[m] > worst["steps"][0]:
        worst["steps"] = (steps_to_1[m], m)
    if arrows_to_1[m] > worst["arrows"][0]:
        worst["arrows"] = (arrows_to_1[m], m)
chk(stat["below"] + stat["one"] == 666666, "count of non-multiples of 3 in [2,10^6]")
chk(worst["peak"] == (150994948, 797162), "worst peak %s" % (worst["peak"],))
chk(150994948 == 9 * 2 ** 24 + 4, "9*2^24+4")
chk(worst["steps"] == (74, 984104), "worst steps %s" % (worst["steps"],))
chk(worst["arrows"] == (172, 984104), "worst arrows %s" % (worst["arrows"],))
print("C5 CONFIRMED (Q2): every m<=10^6, 3 not|m, reduces to 1 by the greedy map alone; worst overall")
print("   E-node peak 150994948 = 9*2^24+4 at m=797162; 74 compound moves / 172 arrows at m=984104. (%.1fs)"
      % (time.time() - t))
# note: 'peak' counts intermediate doubled values 2^k x, i.e. E-nodes of the forward path 1 -> m.
# m=2 and m=4 hit 1 in one move; every other m descends strictly below itself first
chk(G(2)[0] == 1 and G(4)[0] == 1, "m=2,4 hit 1")
t = time.time()
for n in range(2, NQ + 1):
    x = n
    while x >= n:
        x = x // 2 if x % 2 == 0 else 3 * x + 1
print("C5 CONFIRMED (Q1): deterministic Collatz descends below n for 2<=n<=10^6 (%.1fs)." % (time.time() - t))

# ---------------------------------------------------------------------------
banner("A4  SCC data by iterative Kosaraju (C6)")
# ---------------------------------------------------------------------------


def kosaraju(N, up):
    succ = [[] for _ in range(N + 1)]
    pred = [[] for _ in range(N + 1)]
    for n in range(1, N + 1):
        if n % 2 == 0:
            succ[n].append(n // 2)
            pred[n // 2].append(n)
        u = up(n)
        if 1 <= u <= N:
            succ[n].append(u)
            pred[u].append(n)
    order = []
    seen = [False] * (N + 1)
    for r in range(1, N + 1):
        if seen[r]:
            continue
        seen[r] = True
        st = [(r, 0)]
        while st:
            n, i = st[-1]
            if i < len(succ[n]):
                st[-1] = (n, i + 1)
                w = succ[n][i]
                if not seen[w]:
                    seen[w] = True
                    st.append((w, 0))
            else:
                st.pop()
                order.append(n)
    comp = [0] * (N + 1)
    ncomp = 0
    sizes = []
    for r in reversed(order):
        if comp[r]:
            continue
        ncomp += 1
        comp[r] = ncomp
        sz = 1
        st = [r]
        while st:
            n = st.pop()
            for w in pred[n]:
                if not comp[w]:
                    comp[w] = ncomp
                    sz += 1
                    st.append(w)
        sizes.append(sz)
    return comp, sizes


EXPECT = {10 ** 3: (854, 147, 520, 202, [31, 47, 55, 71, 73]),
          10 ** 4: (8249, 1752, 4915, 1737, [383, 511, 575, 608, 667]),
          10 ** 5: (82924, 17077, 49590, 17764, [1535, 2047, 2207, 2287, 2303])}
for N in (10 ** 3, 10 ** 4, 10 ** 5):
    t = time.time()
    comp, sizes = kosaraju(N, lambda n: 3 * n + 1)
    c1 = comp[1]
    giant = sum(1 for n in range(1, N + 1) if comp[n] == c1)
    nontriv = sum(1 for s in sizes if s > 1)
    chk(nontriv == 1 and giant == max(sizes), "one nontrivial SCC containing 1 at N=%d" % N)
    chk(all(sizes[comp[n] - 1] == 1 for n in range(3, N + 1, 3)), "3Z singleton at N=%d" % N)
    out = [n for n in range(1, N + 1) if n % 3 and comp[n] != c1]
    half = sum(1 for n in out if n <= N // 2)
    e = EXPECT[N]
    chk((len(sizes), giant, len(out), half, out[:5]) == e, "SCC data at N=%d: %s" % (N, (len(sizes), giant, len(out), half, out[:5])))
    print("  N=%d: #SCC=%d giant=%d outsiders=%d (<=N/2: %d) smallest=%s  (%.1fs)" % (N, len(sizes), giant, len(out), half, out[:5], time.time() - t))
    if N == 10 ** 5:
        for n, pk_exp in zip(out[:5], (118096, 1276936, 190996, 250504, 118096)):
            x = n
            pk = n
            while x != 1:
                x = x // 2 if x % 2 == 0 else 3 * x + 1
                pk = max(pk, x)
            chk(pk == pk_exp, "Collatz peak of %d = %d" % (n, pk))
print("C6 CONFIRMED: Kosaraju reproduces all SCC counts, giant sizes, outsider counts, smallest outsiders, peaks.")

# ---------------------------------------------------------------------------
banner("A5  word determination, chain law, E[2^K_J], tail bound (C7, C8)")
# ---------------------------------------------------------------------------


def word(m, J):
    ks = []
    x = m
    for _ in range(J):
        x, k = G(x)
        ks.append(k)
    return ks, x


# exhaustive: for J<=6, all units mod 3^(J+1) and lifts by 3^(J+1) t, t<=8: same word; and
# exact sharpness: for every J<=8 the pair (1, 1+3^J) has different J-th letters.
for J in range(1, 7):
    mod = 3 ** (J + 1)
    for m in range(1, mod):
        if m % 3 == 0:
            continue
        w0, x0 = word(m, J)
        for tt in (1, 2, 5, 8):
            w1, x1 = word(m + mod * tt, J)
            chk(w1 == w0 and x1 % 3 == x0 % 3, "word not determined mod 3^(J+1)")
for J in range(1, 9):
    chk(word(1, J)[0] == [2] * J and word(1 + 3 ** J, J)[0] == [2] * (J - 1) + [0], "sharpness J=%d" % J)
print("C7 CONFIRMED: first J letters and m_J mod 3 depend only on m mod 3^(J+1) (exhaustive J<=6); sharp.")

# chain enumeration, exact rationals
FJ_EXP = [Fraction(2, 3), Fraction(8, 9), Fraction(25, 27), Fraction(26, 27), Fraction(236, 243),
          Fraction(239, 243), Fraction(241, 243), Fraction(2173, 2187), Fraction(19609, 19683), Fraction(58868, 59049)]
GJ_EXP = [Fraction(2, 3), Fraction(5, 6), Fraction(43, 54), Fraction(73, 81), Fraction(214, 243),
          Fraction(686, 729), Fraction(2126, 2187), Fraction(12647, 13122), Fraction(19339, 19683), Fraction(6413, 6561)]
t = time.time()
for J in range(1, 11):
    mod = 3 ** (J + 1)
    tot = 0
    e2 = 0
    fcount = 0
    gcount = 0
    lawJ = {}
    for m in range(1, mod):
        if m % 3 == 0:
            continue
        tot += 1
        x = m
        K = 0
        hit = False
        for i in range(1, J + 1):
            k = KMIN[x % 9]
            K += k
            x = ((x << k) - 1) // 3
            if 2 ** K < 3 ** i:
                hit = True
            if i == J:
                lawJ[x % 9] = lawJ.get(x % 9, 0) + 1
        fcount += hit
        gcount += 2 ** K < 3 ** J
        e2 += 2 ** K
    chk(tot == 2 * 3 ** J, "units")
    chk(Fraction(e2, tot) == 3 * Fraction(7, 3) ** (J - 1), "E[2^K_J] at J=%d" % J)
    chk(Fraction(fcount, tot) == FJ_EXP[J - 1], "f_J at J=%d: %s" % (J, Fraction(fcount, tot)))
    chk(Fraction(gcount, tot) == GJ_EXP[J - 1], "g_J at J=%d: %s" % (J, Fraction(gcount, tot)))
    chk(1 - Fraction(fcount, tot) <= Fraction(7, 9) ** (J - 1), "tail bound J=%d" % J)
    # law of m_J mod 9 over units mod 3^(J+1): 2/9 on {1,4,7}, 1/9 on {2,5,8}
    chk({r: Fraction(c, tot) for r, c in lawJ.items()} ==
        {1: Fraction(2, 9), 4: Fraction(2, 9), 7: Fraction(2, 9), 2: Fraction(1, 9), 5: Fraction(1, 9), 8: Fraction(1, 9)},
        "stationary law at step J=%d" % J)
print("C8 CONFIRMED: E[2^K_J] = 3(7/3)^(J-1), f_J, g_J tables, 1-f_J <= (7/9)^(J-1), law 2/9,1/9 for J<=10 (%.1fs)."
      % (time.time() - t))
# the Markov chain transition: exact conditional law of r_{i+1} given the whole past = uniform on coset
mod = 3 ** 5
for m in range(1, mod):
    if m % 3 == 0:
        continue
    x = m
    for i in range(3):
        c = NXT[x % 9]
        lifts = set()
        for tt in range(3):
            y = m + tt * 3 ** (i + 2)
            z = y
            for _ in range(i + 1):
                z = ((z << KMIN[z % 9]) - 1) // 3
            lifts.add(z % 9)
        chk(lifts == {c, c + 3, c + 6}, "three lifts do not biject onto the coset")
        x = ((x << KMIN[x % 9]) - 1) // 3
M = [[Fraction(5, 3), Fraction(1, 3)], [Fraction(10, 3), Fraction(2, 3)]]
chk(M[0][0] * M[1][1] - M[0][1] * M[1][0] == 0 and M[0][0] + M[1][1] == Fraction(7, 3), "tilt spectrum")
uA = Fraction(4 + 2 + 1 + 8, 6)
uB = Fraction(1 + 2, 6)
chk(uA == Fraction(5, 2) and uB == Fraction(1, 2), "initial vector")
chk(uA * (M[0][0] + M[0][1]) + uB * (M[1][0] + M[1][1]) == 7, "uM1=7")
chk(Fraction(2, 9) * 2 + Fraction(1, 9) * 1 + Fraction(1, 9) * 3 + Fraction(1, 9) * 1 == 1, "E_pi[k]=1")
print("   Markov structure (three lifts biject onto the coset), tilted matrix, u M 1 = 7, E_pi[k]=1: CONFIRMED.")
# mismatch check J=6, m<=10^5 (prediction vs actual descent)
mism = 0
for m in range(2, 10 ** 5 + 1):
    if m % 3 == 0:
        continue
    x = m
    K = 0
    pd = ad = False
    for i in range(1, 7):
        k = KMIN[x % 9]
        K += k
        x = ((x << k) - 1) // 3
        pd |= 2 ** K < 3 ** i
        ad |= x < m
    chk(not (pd and not ad), "prediction => descent violated at m=%d" % m)
    mism += pd != ad
chk(mism == 0, "mismatches J=6: %d" % mism)
print("   J=6, m<=10^5: zero mismatches between residue prediction and actual descent: CONFIRMED.")

# ---------------------------------------------------------------------------
banner("A6  cycle census (C9): own DFS (nodes<=2000, len<=40) + GLOBAL word census len<=24")
# ---------------------------------------------------------------------------


def census_dfs(LIM, MAXLEN, up):
    found = []
    for s in range(1, LIM + 1):
        # iterative DFS over simple paths from s using nodes in [s, LIM]
        st = [(s, 0)]
        path = [s]
        onp = {s}
        while st:
            n, i = st[-1]
            nb = []
            if n % 2 == 0:
                nb.append(n // 2)
            u = up(n)
            if u <= LIM:
                nb.append(u)
            if i >= len(nb):
                st.pop()
                path.pop()
                onp.discard(n)
                continue
            st[-1] = (n, i + 1)
            w = nb[i]
            if w == s:
                found.append(tuple(path))
                continue
            if w < s or w in onp or len(path) >= MAXLEN:
                continue
            st.append((w, 0))
            path.append(w)
            onp.add(w)
    return found


def stats(c, sgn):
    L = len(c)
    a = sum(1 for i in range(L) if c[(i + 1) % L] == 3 * c[i] + sgn)
    h = L - a
    e = sum(1 for i in range(L) if c[(i + 1) % L] == 3 * c[i] + sgn and c[i] % 2 == 0)
    return a, h, e


t = time.time()
CY = census_dfs(2000, 40, lambda n: 3 * n + 1)
hist = {}
ah = {}
for c in CY:
    hist[len(c)] = hist.get(len(c), 0) + 1
    a, h, e = stats(c, 1)
    ah[(a, h)] = ah.get((a, h), 0) + 1
    chk(2 ** h > 3 ** a and all(x % 3 for x in c), "cycle equation / 3Z")
    chk((e >= 1) == (c != (1, 4, 2)), "even->3n+1 usage")
    chk(h == min(tt for tt in range(100) if 2 ** tt > 3 ** a), "h = ceil(a log2 3)")
chk(len(CY) == 74, "74 cycles, got %d" % len(CY))
chk(sorted(hist.items()) == [(3, 1), (8, 1), (13, 6), (16, 1), (21, 2), (26, 11), (34, 9), (39, 43)], "histogram")
chk(sorted(ah) == [(1, 2), (3, 5), (5, 8), (6, 10), (8, 13), (10, 16), (13, 21), (15, 24)], "(a,h) set")
chk((2, 7, 22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4) in CY, "user's 16-cycle")
print("C9 CONFIRMED (bounded universe): 74 cycles, histogram, (a,h) values, h=ceil(a log2 3) (%.2fs)." % (time.time() - t))

# GLOBAL census by words: a cycle of length L with a '3n+1' steps and h halvings has
# n0 (2^h - 3^a) = B, B determined by the word; enumerate every word, test validity.
t = time.time()
GLOBAL = {}
LMAX = 24
for L in range(1, LMAX + 1):
    seen_sets = set()
    for mask in range(1, 1 << L):        # bit set = '3n+1' step at that position
        a = bin(mask).count("1")
        h = L - a
        if 2 ** h <= 3 ** a:
            continue
        # forward accumulation: state (3^a n0 + B)/2^hh
        B = 0
        hh = 0
        for i in range(L):
            if (mask >> i) & 1:
                B = 3 * B + 2 ** hh
            else:
                hh += 1
        D = 2 ** h - 3 ** a
        if B % D:
            continue
        n0 = B // D
        if n0 <= 0:
            continue
        x = n0
        nodes = [x]
        ok = True
        for i in range(L):
            if (mask >> i) & 1:
                x = 3 * x + 1
            else:
                if x % 2:
                    ok = False
                    break
                x //= 2
            nodes.append(x)
        if not ok or x != n0:
            continue
        cyc = nodes[:-1]
        if len(set(cyc)) != L:
            continue           # not simple (a shorter cycle traversed several times)
        key = frozenset(cyc)
        if key in seen_sets:
            continue
        seen_sets.add(key)
        GLOBAL.setdefault(L, []).append((a, h, min(cyc), max(cyc), tuple(cyc)))
for L in sorted(GLOBAL):
    lst = GLOBAL[L]
    ahset = sorted({(a, h) for a, h, *_ in lst})
    big = sum(1 for a, h, mn, mx, c in lst if mx > 2000)
    hforced = all(h == min(tt for tt in range(100) if 2 ** tt > 3 ** a) for a, h, *_ in lst)
    print("  L=%2d: %3d cycles globally, (a,h)=%s, with max node > 2000: %d, h=ceil(a log2 3) for all: %s, max node overall %d"
          % (L, len(lst), ahset, big, hforced, max(mx for *_, mx, c in lst)))
# compare against bounded census for L <= 24 (there L in {3,8,13,16,21})
for L in (3, 8, 13, 16, 21):
    bounded = sum(1 for c in CY if len(c) == L)
    glob = len(GLOBAL.get(L, []))
    print("  L=%d: bounded(<=2000) census %d vs global %d" % (L, bounded, glob))
chk(len(GLOBAL.get(16, [])) >= 1, "16-cycle")
print("  GLOBAL length-16 cycles: %s" % [c for *_, c in GLOBAL[16]])
print("  (%.1fs)" % (time.time() - t))

# ---------------------------------------------------------------------------
banner("A7  signed side E_- (C10)")
# ---------------------------------------------------------------------------
for v in range(1, V + 1):
    p3 = (v + 1) // 3 if (v + 1) % 3 == 0 else None
    new = p3 is not None and p3 % 2 == 0
    if new:
        chk(v % 6 == 5, "new minus target")
        n0 = (4 * v + 1) // 3
        chk((4 * v + 1) % 3 == 0 and n0 % 8 == 7 and 3 * n0 - 1 == 4 * v and 4 * p3 - 1 == n0, "minus leaf identity v=%d" % v)
        for n in range(1, n0, 2):
            x = 3 * n - 1
            while x % 2 == 0:
                x //= 2
            chk(x != v, "smaller odd T_- predecessor at v=%d" % v)
    else:
        chk(v % 6 != 5, "missing minus arrow")
KM = {}
for r in U9:
    ks = [k for k in range(6) if (pow(2, k, 9) * r) % 9 in (2, 5)]
    KM[r] = min(ks)
    chk(KM[r] == KMIN[9 - r], "negation conjugacy")
chk(KM == {1: 1, 2: 0, 4: 3, 5: 0, 7: 1, 8: 2}, "minus table")
print("C10 CONFIRMED: minus leaf identity (n0 = 7 mod 8 least odd T_- predecessor, brute-forced), k_-(r)=k_+(9-r).")


def Gm(x):
    k = KM[x % 9]
    return ((x << k) + 1) // 3, k


chk(Gm(4) == (11, 3) and Gm(11) == (4, 0), "{4,11} 2-cycle")
p = [242]
for _ in range(5):
    p.append(Gm(p[-1])[0])
chk(p == [242, 323, 431, 575, 767, 256], "minus hostile")
# rescue path of 4 and its k-word
resc = [4, 11, 59, 20, 7, 5, 2]
kw = []
for x, y in zip(resc, resc[1:]):
    ks = [k for k in range(12) if ((x << k) + 1) % 3 == 0 and ((x << k) + 1) // 3 == y]
    chk(len(ks) == 1, "step %d->%d" % (x, y))
    kw.append(ks[0])
chk(kw == [3, 4, 0, 0, 1, 0], "k-word of the rescue path is %s" % kw)
print("   NOTE k-word of 4->11->59->20->7->5->2 is (3,4,0,0,1,0), NOT (3,4,0,1,1,0) as the draft states (20->7 has k=0).")


def bfs_below(m, sign, kcap=40, budget=200000):
    par = {m: None}
    dq = deque([m])
    n = 0
    while dq:
        x = dq.popleft()
        n += 1
        if n > budget:
            return None
        for k in range(kcap + 1):
            y = (x << k) - sign
            if y % 3:
                continue
            z = y // 3
            if z <= 0 or z % 3 == 0 or z in par:
                continue
            par[z] = x
            if z < m:
                pth = [z]
                while par[pth[-1]] is not None:
                    pth.append(par[pth[-1]])
                return pth[::-1]
            dq.append(z)
    return None


t = time.time()
tp = [0] * (NQ + 1)
ts = [0] * (NQ + 1)
tp[1] = 1
wp = (1, 1)
ws = (0, 1)
nongreedy = []
for m in range(2, NQ + 1):
    if m % 3 == 0:
        continue
    x = m
    s = 0
    pk = m
    cyc = False
    while x >= m:
        k = KM[x % 9]
        y = x << k
        pk = max(pk, y)
        x = (y + 1) // 3
        s += 1
        if x == m or s > 10 ** 5:
            cyc = True
            break
    if cyc:
        pth = bfs_below(m, -1)
        chk(pth is not None, "Q2_- fails at %d" % m)
        nongreedy.append((m, pth))
        x = pth[-1]
        pk = max(pth)
        s = len(pth) - 1
    tp[m] = max(pk, tp[x])
    ts[m] = s + ts[x]
    if tp[m] > wp[0]:
        wp = (tp[m], m)
    if ts[m] > ws[0]:
        ws = (ts[m], m)
chk(nongreedy == [(4, [4, 11, 59, 20, 7, 5, 2])], "non-greedy starts: %s" % nongreedy[:3])
chk(wp == (150994940, 797161) and 150994940 == 9 * 2 ** 24 - 4, "Q2_- peak %s" % (wp,))
chk(ws == (82, 919795), "Q2_- max moves %s" % (ws,))
print("   Q2_- to 10^6 CONFIRMED: only m=4 needs a non-greedy move; peak 9*2^24-4 at 797161; 82 moves at 919795 (%.1fs)."
      % (time.time() - t))
# Q1_-: deterministic descent except 5, 17
t = time.time()
nondesc = []
for n in range(2, NQ + 1):
    x = n
    s = 0
    while x >= n:
        x = x // 2 if x % 2 == 0 else 3 * x - 1
        s += 1
        if x == n or s > 10 ** 5:
            nondesc.append(n)
            break
chk(nondesc == [5, 17], "non-descending 3n-1 starts: %s" % nondesc[:5])
print("   Q1_- to 10^6 CONFIRMED: deterministic 3n-1 descends except n=5,17 (%.1fs)." % (time.time() - t))


def walk_ok(w, sgn):
    return all(y == 3 * x + sgn or (x % 2 == 0 and y == x // 2) for x, y in zip(w, w[1:]))


W5 = [5, 14, 7, 20, 10, 29, 86, 43, 128, 64, 32, 16, 8, 4, 2, 1]
W17in = [1, 2, 5, 14, 7, 20, 59, 176, 88, 44, 22, 11, 32, 16, 8, 4, 11, 32, 16, 8, 23, 68, 34, 17]
W17out = [17, 50, 25, 74, 37, 110, 55, 164, 82, 41, 122, 61, 182, 91, 272, 136, 68, 34, 101, 302, 151, 452, 226,
          113, 338, 169, 506, 253, 758, 379, 1136, 568, 284, 142, 71, 212, 106, 53, 158, 79, 236, 118, 59, 176,
          88, 44, 22, 11, 32, 16, 8, 4, 2, 1]
chk(walk_ok(W5, -1) and walk_ok(W17in, -1) and walk_ok(W17out, -1), "note's walks are legal E_- walks")
chk(len(W5) - 1 == 15 and len(W17in) - 1 == 23 and len(W17out) - 1 == 53, "arrow counts")
chk(len(set(W17in)) < len(W17in), "the 1->17 walk repeats nodes")
W17simple = [1, 2, 5, 14, 7, 20, 59, 176, 88, 44, 22, 11, 32, 16, 8, 23, 68, 34, 17]
chk(walk_ok(W17simple, -1) and len(set(W17simple)) == len(W17simple) and len(W17simple) - 1 == 18, "simple 1->17 path")
chk(len(set(W17out)) == len(W17out) and len(set(W5)) == len(W5), "other two are simple")
print("   NOTE: the displayed 1->17 'path' (23 arrows) is a WALK repeating 11,32,16,8 (detour 8->4->11);")
print("   the simple path 1->2->5->14->7->20->59->176->88->44->22->11->32->16->8->23->68->34->17 has 18 arrows.")
t = time.time()
comp, sizes = kosaraju(10 ** 5, lambda n: 3 * n - 1)
c1 = comp[1]
giant = sum(1 for n in range(1, 10 ** 5 + 1) if comp[n] == c1)
chk(len(sizes) == 82979 and giant == 17022 and sum(1 for s in sizes if s > 1) == 1, "E_- SCC data (%d,%d)" % (len(sizes), giant))
chk(comp[5] == c1 and comp[17] == c1, "cycles in giant")
print("   E_-|[1,10^5] Kosaraju: 82979 SCCs, giant 17022, one nontrivial: CONFIRMED (%.1fs)." % (time.time() - t))
CYM = census_dfs(2000, 40, lambda n: 3 * n - 1)
hm = {}
dets = []
for c in CYM:
    hm[len(c)] = hm.get(len(c), 0) + 1
    a, h, e = stats(c, -1)
    chk(2 ** h < 3 ** a and all(x % 3 for x in c), "E_- cycle equation")
    if e == 0:
        dets.append(c)
chk(len(CYM) == 70 and sorted(hm.items()) == [(2, 1), (5, 2), (15, 2), (18, 12), (23, 4), (28, 1), (31, 2), (36, 46)], "E_- census")
chk(sorted(dets) == sorted([(1, 2), (5, 14, 7, 20, 10),
                            (17, 50, 25, 74, 37, 110, 55, 164, 82, 41, 122, 61, 182, 91, 272, 136, 68, 34)]), "deterministic E_- cycles")
print("   E_- census: 70 cycles, histogram, 2^h<3^a, deterministic = three known cycles: CONFIRMED.")

# ---------------------------------------------------------------------------
banner("A8  pn+1 hierarchy (C11)")
# ---------------------------------------------------------------------------


def ordp(p):
    o = 1
    while pow(2, o, p) != 1:
        o += 1
    return o


def pinfo(p):
    o = ordp(p)
    e = ((pow(2, o, p * p) - 1) // p) % p
    chk(e != 0, "Wieferich")
    reach = sorted({pow(2, t, p) for t in range(o)})          # residues s with 2^k s = 1 solvable <=> s in <2>
    reach = [s for s in range(1, p) if any(pow(2, t, p) * s % p == 1 for t in range(o))]
    k0 = {s: min(t for t in range(o) if pow(2, t, p) * s % p == 1) for s in reach}
    # brute minimal admissible k for every unit r mod p^2 with s in reach, and next class
    k = {}
    c = {}
    for r in range(1, p * p):
        if r % p == 0 or (r % p) not in k0:
            continue
        for t in range(0, o * p + 2):
            y = pow(2, t, p * p) * r % (p * p)
            if y % p == 1 and ((y - 1) // p) % p != 0:
                k[r] = t
                c[r] = ((y - 1) // p) % p
                break
        s = r % p
        a = ((pow(2, k0[s], p * p) * r - 1) // p) % p
        chk(k[r] == k0[s] + o * (a == 0) and c[r] == (a if a else e), "k/c formula p=%d r=%d" % (p, r))
    for s in reach:
        law = {}
        for t in range(p):
            law[c[s + p * t]] = law.get(c[s + p * t], 0) + 1
        chk(law == {cc: 1 + (cc == e) for cc in range(1, p)}, "class law p=%d s=%d" % (p, s))
    return o, e, reach, k0, k, c


ROWS = {}
for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 37):
    o, e, reach, k0, k, c = pinfo(p)
    prim = (o == p - 1)
    if prim:
        Ek = sum(Fraction(1 + (s == e), p) * k0[s] for s in reach) + Fraction(o, p)
        rho = Fraction(sum(2 ** k0[s] for s in reach) + 2 ** (o + k0[e]), p)
        chk(Ek == Fraction(p - 1, 2) + Fraction(k0[e], p), "E_pi[k] closed form p=%d" % p)
        # rank-one tilted matrix eigenvalue = rho  (v^T u)
        ROWS[p] = (o, e, Ek, 2 ** Ek.numerator < p ** Ek.denominator, rho, rho < p)
        # general inequalities
        chk(rho >= Fraction(2 ** p - 1, p), "rho >= (2^p-1)/p")
        if p >= 5:
            chk(rho > p and 2 ** p - 1 > p * p, "mean-supercritical for p>=5")
        if p >= 7:
            chk(2 ** Ek.numerator > p ** Ek.denominator and 2 ** (p - 1) >= p ** 2, "log-supercritical for p>=7")
    else:
        ROWS[p] = (o, e, None, None, None, None, reach, Fraction(p - 1 - len(reach), p) if False else Fraction(sum(1 for t in range(p) if c[reach[0] + p * t] not in k0), p))
for p, row in ROWS.items():
    print("  p=%2d: %s" % (p, row))
chk(ROWS[3][2] == 1 and ROWS[3][4] == Fraction(7, 3), "p=3")
chk(ROWS[5][2] == Fraction(11, 5) and ROWS[5][4] == Fraction(47, 5), "p=5")
chk(ROWS[11][2] == Fraction(61, 11) and ROWS[11][4] == Fraction(66559, 11), "p=11")
chk(ROWS[13][2] == Fraction(86, 13) and ROWS[13][4] == Fraction(1052671, 13), "p=13")
chk(ROWS[19][2] == Fraction(176, 19) and ROWS[19][4] == Fraction(8650751, 19), "p=19")
chk([p for p in ROWS if ROWS[p][2] is not None and ROWS[p][5]] == [3], "mean-subcritical only p=3")
chk([p for p in ROWS if ROWS[p][2] is not None and ROWS[p][3]] == [3, 5], "log-subcritical only p=3,5")
print("C11 CONFIRMED: k/c formulas, class law, E_pi[k], rho_p for p in {3,5,11,13,19} (+23,29,37 as extra).")
print("   STRENGTHENING (elementary): rho_p >= (2^p-1)/p > p for every p>=5 with 2 a primitive root (non-Wieferich),")
print("   and E_pi[k] >= (p-1)/2 >= log2 p for p>=7; so mean-subcriticality <=> p=3 and log-subcriticality <=> p in {3,5}")
print("   hold for ALL such primes, not only the tested list.")
lam = Fraction(21, 20)
o, e, reach, k0, k, c = pinfo(5)
rl = (sum(lam ** k0[s] for s in reach) + lam ** (o + k0[e])) / 5
chk(rl == Fraction(17876501, 16000000) and rl ** 10 < lam ** 23 and 2 ** 23 < 5 ** 10, "p=5 certificate")
print("   p=5 Chernoff certificate rho(21/20)^10 < (21/20)^23, 2^23 < 5^10: CONFIRMED (rate %.5f)." % float((rl ** 10 / lam ** 23) ** Fraction(1, 10)))
# E[lambda^K_J] enumeration vs rank-one formula, p=5, J<=3
for lamb in (Fraction(2), lam):
    S = reach
    Mx = [[lamb ** k0[s] / 5 * (1 + (lamb ** o if s2 == e else 0)) for s2 in S] for s in S]
    for J in range(1, 4):
        vec = [Fraction(1)] * 4
        for _ in range(J):
            vec = [sum(Mx[i][j] * vec[j] for j in range(4)) for i in range(4)]
        pred = sum(vec) / 4
        acc = Fraction(0)
        tot = 0
        for m in range(1, 5 ** (J + 1)):
            if m % 5 == 0:
                continue
            tot += 1
            x = m
            K = 0
            for _ in range(J):
                kk = k[x % 25]
                K += kk
                x = ((x << kk) - 1) // 5
            acc += lamb ** K
        chk(acc / tot == pred, "E[lambda^K_J] p=5")
print("   E[lambda^K_J] = (1/4) 1^T M^J 1 for p=5, J<=3, lambda in {2,21/20}: CONFIRMED.")
# E_5 greedy cycles
cyc5 = []
for m in range(2, 10 ** 5 + 1):
    if m % 5 == 0:
        continue
    x = m
    s = 0
    while x >= m:
        kk = k[x % 25]
        x = ((x << kk) - 1) // 5
        s += 1
        if x == m:
            cyc5.append(m)
            break
        chk(s < 10 ** 4, "cap")
chk(cyc5 == [13, 17], "E_5 greedy cycle minima %s" % cyc5)
c13 = [13]
while True:
    c13.append(((c13[-1] << k[c13[-1] % 25]) - 1) // 5)
    if c13[-1] == 13:
        break
chk(c13[:-1] == [13, 83, 33], "13-cycle")
# deterministic 5n+1 cycles through 13 and 17
for start in (13, 17):
    x = start
    for _ in range(50):
        x = x // 2 if x % 2 == 0 else 5 * x + 1
        if x == start:
            break
    chk(x == start, "5n+1 cycle at %d" % start)
# E_7: reachable residues from 1 are in <2> = {1,2,4}; residue 3 unreachable (m=3 refutes Q2_7)
seen = {1}
st = [1]
while st:
    x = st.pop()
    for y in ([x // 2] if x % 2 == 0 else []) + ([7 * x + 1] if 7 * x + 1 <= 10 ** 4 else []):
        if y not in seen:
            seen.add(y)
            st.append(y)
chk({x % 7 for x in seen} == {1, 2, 4}, "E_7 residues")
print("   E_5 greedy cycles {13,83,33},{17,27,43} (reversed 5n+1 cycles), E_7 residues {1,2,4}: CONFIRMED.")

# ---------------------------------------------------------------------------
banner("A9  miscellaneous note statements")
# ---------------------------------------------------------------------------
# 2-adic block 31 = 2^5-1 (L=4): odd T-orbit is 31 -> 47 -> 71 -> 107 -> 161 = 2*3^4-1, NOT 31 -> 121


def Todd(n):
    x = 3 * n + 1
    while x % 2 == 0:
        x //= 2
    return x


orb = [31]
for _ in range(4):
    orb.append(Todd(orb[-1]))
chk(orb == [31, 47, 71, 107, 161] and 161 == 2 * 3 ** 4 - 1, "2-adic block")
print("   NOTE: '31 -> 121 -> ...' is wrong; T-orbit of 31 = 2^5-1 is 31 -> 47 -> 71 -> 107 -> 161 = 2*3^4-1.")
# (1,4,2): the arrow 1 -> 4 is 1 + (2*1+1) = 1 + 3, parents (1,3) distinct; the diagonal (1,1)->2 is the
# shortcut arrow 1 -> 2, which is NOT an E-arrow (2 is reached from 4 by halving).
chk(3 * 1 + 1 == 4 and 1 + 3 == 4 and 2 * 1 + 1 != 1, "?")
print("   NOTE: 'E's cycle (1,4,2) uses the diagonal parent pair (1,1)' is wrong in E: 1->4 = 1+3 has distinct")
print("   parents (1,3); the diagonal (1,1)->2 belongs to the shortcut map, and 1->2 is not an E-arrow.")
# 'index -2 is never integral': (v/4-1)/3 for odd v; indices <= -2 all non-integral: (4^{j+1} v - 1)/3 with j<=-2
for v in range(7, 1000, 6):
    for j in (-2, -3, -4):
        chk((Fraction(4) ** (j + 1) * v - 1) / 3 != int((Fraction(4) ** (j + 1) * v - 1) / 3), "index j integral?")
print("   (B1) at indices <= -2 non-integral for v = 1 mod 6 < 1000: CONFIRMED.")

print()
with open(__file__, "rb") as fh:
    print("audit source sha256: %s" % hashlib.sha256(fh.read()).hexdigest())
print("python %s; total time %.1fs" % (sys.version.split()[0], time.time() - T0))
print("ALL AUDIT CHECKS PASSED")
