#!/usr/bin/env python3
"""
collatz_mod6_20260917_three_adic_g_map.py  (LANE = three_adic_g_map)

The greedy 3-adic inverse map of firsthand fact F1:

    G(m) = (2^k m - 1)/3,  k >= 0 minimal with 2^k m = 4 or 7 (mod 9),
    domain  {m in Z : 3 does not divide m}  (we study m >= 1 and m <= -1).

The G-orbit of m reversed is a path from 1 to m in the user's extended graph E
(odd n -> 3n+1, even n -> 3n+1 and n -> n/2), because m' = G(m) satisfies
3m'+1 = 2^k m.

Sections
  S1  exact transfer law at every 3-adic level, Markov structure (PROVED + enumeration J<=7);
      S1.3 conjugacy of (Z_3^x, G) to a 6-state SFT, k-word as a 2-block recoding, entropy log 3
  S2  invariant measure, stationary law, ergodicity, exact drift (PROVED + enumeration)
  S3  Terras-type density theorem, d_J exact for J<=12, Collatz comparison;
      S3.5 exact exponential moment for every tilt, sharp Chernoff exponent, proved sigma=sigma_res threshold
  S4  G_k for odd k, 3 not | k: E[j]=1 for every k, cycle gate, census |k|<=49, duality, scaling lemma
  S5  hostiles: 3^j+1, the exact peak bound, worst m<=10^7, (8,5,1^n) family
  S6  reframe theorem (Q1, Q2, one SCC), typed analogies incl. the blueprint affine word model
All checks are explicit `check(...)` calls (active under python -O).
"""
import sys, time, hashlib
from fractions import Fraction as Fr
from math import log, gcd
import numpy as np

T0 = time.time()
FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        raise RuntimeError("CHECK FAILED: " + msg)


def banner(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


# ----------------------------------------------------------------------------
# basic maps
# ----------------------------------------------------------------------------
def kmin(m):
    """greedy exponent for k=1: minimal k>=0 with 2^k m = 4 or 7 mod 9 (depends on m mod 9)."""
    r = m % 9
    k = 0
    while (r % 3 != 1) or (r == 1):
        r = (2 * r) % 9
        k += 1
        if k > 6:
            raise RuntimeError("no k for residue %d" % (m % 9))
    return k


KT = {r: kmin(r) for r in (1, 2, 4, 5, 7, 8)}
KT_ARR = np.zeros(9, dtype=np.int64)
for r, k in KT.items():
    KT_ARR[r] = k


def G(m):
    k = KT[m % 9]
    return ((m << k) - 1) // 3, k


def jtab_k(k):
    """greedy exponent table for G_k: minimal j with 2^j m = k+3 or k+6 mod 9."""
    tab = {}
    t1, t2 = (k + 3) % 9, (k + 6) % 9
    for r in (1, 2, 4, 5, 7, 8):
        x, j = r, 0
        while x != t1 and x != t2:
            x = (2 * x) % 9
            j += 1
            if j > 6:
                raise RuntimeError("no j")
        tab[r] = j
    return tab


def Gk_factory(k):
    tab = jtab_k(k)

    def Gk(m):
        j = tab[m % 9]
        return (m * (1 << j) - k) // 3, j

    return Gk


# ============================================================================
banner("S0  Definitions, universe, inheritance")
# ============================================================================
print("G(m) = (2^k m - 1)/3, k minimal with 2^k m in {4,7} mod 9; domain: 3 not | m.")
print("k-table (depends only on m mod 9):", KT)
check(KT == {1: 2, 2: 1, 4: 0, 5: 3, 7: 0, 8: 1}, "k-table of F1")
print("Inheritance: 05-knowledge/results/collatz_mod6_20260917_extended_collatz_scc.out (S2, S4, S6, S7),")
print("             05-knowledge/results/arithmetic_braids_20260917_{collatz,summand,divisors}.md")
print("Firsthand facts F1-F4 of the session lead are verified where touched, never re-derived silently.")
for m in (1, 2, 4, 5, 7, 8, 10, 28, 244):
    g, k = G(m)
    check(3 * g + 1 == (m << k), "E-path identity 3G(m)+1 = 2^k m at m=%d" % m)
    check(g % 3 != 0, "G(m) is a 3-adic unit at m=%d" % m)
print("PROVED (S0.1): 3G(m)+1 = 2^k m, so G-orbit reversed is an E-path 1 -> ... -> m (checked samples).")
print("PROVED (S0.2): G(1)=1 (k=2) is a fixed point; 3 never divides G(m) by the exclusion 2^k m != 1 mod 9.")

# ============================================================================
banner("S1  Exact transfer law at every 3-adic level; Markov structure of the residue chain")
# ============================================================================
print("""THEOREM S1.1 (PROVED).  (a) k(m) depends only on m mod 9 (the defining congruence is mod 9).
(b) Let J>=2 and let a be a unit class mod 3^J.  Its three lifts a+t*3^J (t=0,1,2) mod 3^(J+1)
    all have the same k (J+1>=2 so the lift determines m mod 9), and
        G(a + t 3^J) = (2^k a - 1)/3 + t 2^k 3^(J-1)  ==  G(a) + t 2^k 3^(J-1)   (mod 3^J).
    Since 2^k is a unit mod 3, t -> t 2^k is a bijection of Z/3, so the three lifts of a map
    BIJECTIVELY onto the three lifts (mod 3^J) of the class G(a) mod 3^(J-1).  In particular
    G(m) mod 3^J is a function of m mod 3^(J+1), and G maps every cylinder a+3^(J+1)Z_3 (J>=1)
    affinely and bijectively onto the cylinder G(a)+3^J Z_3, scaling the 3-adic metric by 3.
(c) Consequently, for m Haar-distributed on Z_3^x (or on any class mod 3^J, J>=1), the residue
    sequence X_n = G^n(m) mod 3^J is a MARKOV CHAIN with kernel
        P_J(a,b) = 1/3  if  b = G(a) mod 3^(J-1)  (three such unit classes b mod 3^J),  else 0.
    Proof of the Markov property: the history event {X_0=a_0,...,X_n=a_n} is a disjoint union of
    classes D mod 3^(J+n); by (b) iterated n times (all intermediate moduli are >= 3^(J+1) >= 9),
    G^n maps each D affinely onto the full class a_n + 3^J Z_3, pushing Haar|_D to a multiple of
    Haar on that class.  Hence conditionally on the history G^n(m) is Haar on a_n+3^J Z_3, its
    extra digit mod 3^(J+1) is uniform, and by (b) X_(n+1) is uniform on the three lifts of
    G(a_n) mod 3^(J-1).  (For J=1 the target 'mod 3^0' is trivial: P_1(a,1)=2/3, P_1(a,2)=1/3
    because the three lifts of any class mod 3 to mod 9 have G-residues mod 3 equal to 1,1,2.)
(d) Level J=2 is exactly F1: from {1,2,4,5} the next residue is uniform on {1,4,7}
    (G(a) = 1 mod 3), from {7,8} uniform on {2,5,8} (G(a) = 2 mod 3).""")

# exact enumeration of (b) for J<=7.  NOTE (truth, not a weakening): the bijective-lift statement
# needs J>=2 (the three lifts of a class mod 3 to mod 9 have DIFFERENT k, and their images mod 3
# are (1,1,2), not distinct); at J=1 the correct statement is the kernel P_1(a,1)=2/3, P_1(a,2)=1/3.
for J in range(1, 8):
    M1 = 3 ** (J + 1)
    MJ = 3 ** J
    units_J = [a for a in range(MJ) if a % 3]
    for a in units_J:
        imgs = []
        ks = []
        for t in range(3):
            g, k = G(a + t * MJ)
            imgs.append(g % MJ)
            ks.append(k)
        if J == 1:
            check(sorted(imgs) == [1, 1, 2], "J=1: images mod 3 of the three lifts are (1,1,2), a=%d" % a)
            check(ks == ([2, 0, 0] if a == 1 else [1, 3, 1]), "J=1: k along the lifts of a=%d is (2,0,0)/(1,3,1), not constant" % a)
        else:
            check(len(set(ks)) == 1, "lifts share k, J=%d a=%d" % (J, a))
            # distinct, and all congruent to G(a) mod 3^(J-1)
            check(len(set(imgs)) == 3, "bijective lifts J=%d a=%d" % (J, a))
            base = G(a)[0] % (3 ** (J - 1))
            check(all(x % (3 ** (J - 1)) == base for x in imgs), "lifts land over G(a) mod 3^(J-1), J=%d a=%d" % (J, a))
        # also: G(m) mod 3^J depends only on m mod 3^(J+1): test with a second representative
        for t in range(3):
            m1 = a + t * MJ
            m2 = m1 + 5 * M1
            check(G(m1)[0] % MJ == G(m2)[0] % MJ, "G mod 3^J determined by m mod 3^(J+1)")
    if J == 1:
        print("  J=1: both unit classes mod 3: the three lifts mod 9 have NON-constant k ((2,0,0) resp. (1,3,1)) and images (1,1,2) mod 3 -> P_1 = [[2/3,1/3],[2/3,1/3]]  [checked]")
    else:
        print("  J=%d: all %d unit classes mod 3^%d: three lifts (same k) map bijectively onto the three lifts of G(a) mod 3^%d  [checked]"
              % (J, len(units_J), J, J - 1))
# sharpness: mod 3^J does NOT determine G mod 3^J
for J in range(1, 5):
    wit = None
    MJ = 3 ** J
    for a in range(1, MJ):
        if a % 3 == 0:
            continue
        # try both non-trivial lifts a+3^J and a+2*3^J (at J=1, a=1: G(1)=G(4)=1 mod 3 but G(7)=2 mod 3)
        for t in (1, 2):
            if G(a)[0] % MJ != G(a + t * MJ)[0] % MJ:
                wit = (J, a, a + t * MJ, G(a)[0] % MJ, G(a + t * MJ)[0] % MJ)
                break
        if wit is not None:
            break
    print("  sharpness J=%d: witness (a, a+3^J, G mod 3^J) = %s" % (J, wit))
    check(wit is not None, "sharpness witness")
print("FINITE-EXACT (S1.2): transfer law (b) verified for all unit classes at levels J=1..7; modulus 3^(J+1) is sharp.")

# the level-2 kernel explicitly
P2 = {}
for a in (1, 2, 4, 5, 7, 8):
    row = {}
    for t in range(3):
        b = G(a + 9 * t)[0] % 9
        row[b] = row.get(b, 0) + Fr(1, 3)
    P2[a] = row
print("  P_2 (mod 9) rows:")
for a in (1, 2, 4, 5, 7, 8):
    print("    %d -> %s" % (a, dict(sorted(P2[a].items()))))
check(all(set(P2[a]) == {1, 4, 7} for a in (1, 2, 4, 5)), "F1 rows from {1,2,4,5}")
check(all(set(P2[a]) == {2, 5, 8} for a in (7, 8)), "F1 rows from {7,8}")

print("""THEOREM S1.3 (PROVED).  Symbolic conjugacy and entropy.
Let A be the 6x6 0/1 matrix on the unit residues mod 9 with A(a,b)=1 iff b = G(a) (mod 3), i.e.
rows {1,2,4,5} -> columns {1,4,7}, rows {7,8} -> columns {2,5,8}; let Sigma_A be the one-sided SFT.
(i)  The residue coding Phi(m) = (G^n(m) mod 9)_(n>=0) is a homeomorphism Z_3^x -> Sigma_A with
     Phi o G = shift o Phi.  Proof: continuity and shift-equivariance are immediate; injectivity:
     by induction the block (X_0..X_(J-1)) determines m mod 3^(J+1): true for J=1; if m mod 3^(J+1)
     and G(m) mod 3^(J+1) are known, the level-(J+1) transfer law S1.1(b) (three lifts of m mod 3^(J+1)
     have three DISTINCT images mod 3^(J+1)) pins m mod 3^(J+2).  Surjectivity: the admissible blocks
     of length J number 6*3^(J-1) = 2*3^J = the number of unit classes mod 3^(J+1), so the injective
     block map is a bijection at every length, and compactness gives every infinite path.
(ii) The k-word (k_1 k_2 ...) is a 2-block recoding of Phi: k_n determines X_(n-1) except for the
     pairs {4,7} (k=0) and {2,8} (k=1), and k_(n+1) in {0,2} <=> X_n in {1,2,4,5} resolves the pair.
     Hence the k-word map Z_3^x -> {0,1,2,3}^N is injective (a homeomorphism onto a proper SFT), and
     the number of admissible k-words of length J is 4*3^(J-1) (< 2*3^J: a J-letter word determines
     m only mod 3^(J+1) up to the last pair, a (J+1)-letter word determines m mod 3^(J+1) exactly).
(iii) Entropy.  A has constant row sum 3, so its Perron eigenvalue is 3 and h_top(Sigma_A) = log 3;
     the Parry (maximal-entropy) measure has kernel A/3 = P_2 and stationary law = the left Perron
     vector (2,1,2,1,2,1)/9 = pi_2 = mu at level 2 (S2).  At every level J the same holds (P_J = A_J/3,
     left vector pi_J), so mu is the measure of maximal entropy of G and h_mu(G) = h_top(G) = log 3.
     (The k-word shift is conjugate to Sigma_A, hence also of entropy log 3: 4*3^(J-1) words.)
[CITED for the Parry-measure facts: Lind-Marcus, Symbolic Dynamics and Coding, Thm 4.4.4 / Sec. 13.3.]""")
for J in range(1, 8):
    M1 = 3 ** (J + 1)
    paths = {}
    for a in range(M1):
        if a % 3 == 0:
            continue
        x = a
        path = []
        for i in range(J):
            path.append(x % 9)
            x = G(x)[0]
        path = tuple(path)
        check(path not in paths, "residue block of length %d is injective on classes mod 3^%d" % (J, J + 1))
        paths[path] = a
        for i in range(J - 1):
            check(path[i + 1] % 3 == G(path[i])[0] % 3, "admissibility of the coded block")
    check(len(paths) == 2 * 3 ** J == 6 * 3 ** (J - 1), "block count = admissible path count at J=%d" % J)
    # k-word: (J+1)-letter word determines m mod 3^(J+1); J-letter words number 4*3^(J-1)
    words = {}
    M2 = 3 ** (J + 2)
    for a in range(M2):
        if a % 3 == 0:
            continue
        x = a
        w = []
        for i in range(J + 1):
            x, k = G(x)
            w.append(k)
        words.setdefault(tuple(w), set()).add(a % M1)
    check(all(len(v) == 1 for v in words.values()), "(J+1)-letter k-word determines m mod 3^(J+1), J=%d" % J)
    nJ = len(set(w[:J] for w in words))
    check(nJ == 4 * 3 ** (J - 1), "number of J-letter k-words = 4*3^(J-1) at J=%d" % J)
    print("  J=%d: residue blocks <-> unit classes mod 3^%d bijective (%d); k-words of length %d: %d = 4*3^%d; (J+1)-word pins m mod 3^%d  [checked]"
          % (J, J + 1, len(paths), J, nJ, J - 1, J + 1))
# adjacency matrix, Perron data, entropy (exact rationals)
STATES = (1, 2, 4, 5, 7, 8)
A6 = [[1 if b % 3 == G(a)[0] % 3 else 0 for b in STATES] for a in STATES]
check(all(sum(row) == 3 for row in A6), "constant row sum 3")
left = {1: Fr(2, 9), 2: Fr(1, 9), 4: Fr(2, 9), 5: Fr(1, 9), 7: Fr(2, 9), 8: Fr(1, 9)}
for j, b in enumerate(STATES):
    check(sum(left[a] * A6[i][j] for i, a in enumerate(STATES)) == 3 * left[b], "left Perron vector (2,1,2,1,2,1)/9")
h_pi = sum(pi * 3 * Fr(1, 3) for pi in left.values())  # coefficient of log 3 in -sum pi P log P = sum_a pi(a) * log 3
check(h_pi == 1, "h_pi(P_2) = 1 * log 3")
print("  A =", A6)
print("  Perron eigenvalue 3 (row sums), left vector (2,1,2,1,2,1)/9 = pi_2; h_top = h_mu = log 3 = %.6f  [checked]" % log(3))
print("PROVED (S1.4): (Z_3^x, G) is conjugate to the 6-state SFT Sigma_A; the k-word is an injective 2-block recoding;")
print("  entropy log 3 with mu the Parry measure.  FINITE-EXACT: coding bijections verified for J<=7.")

# ============================================================================
banner("S2  Invariant measure, stationary law, ergodicity, exact drift")
# ============================================================================
print("""THEOREM S2.1 (PROVED).  Let w(1)=4/3, w(2)=2/3 and mu = w(m mod 3) * Haar on Z_3^x (Haar normalized
to mass 1 on the units, so each unit class mod 3^(J+1) has Haar mass 1/(2*3^J)).  Then mu is G-invariant.
Proof.  It suffices to check cylinders B = b + 3^J Z_3, J>=1, b a unit.  Preimage classes a mod 3^(J+1)
with G(a) = b mod 3^J satisfy 2^k a = 3b+1 mod 3^(J+1) for k = k(a), i.e. a = 2^(-k)(3b+1) mod 3^(J+1):
at most one class per k, valid iff k(a)=k.  With 3b+1 = 4 mod 9 (b = 1 mod 3): the candidates
a = 4*2^(-k) mod 9 are 4,2,1,5 for k=0,1,2,3 and each has k(a)=k (k>=4 never occurs: k(a)<=3), so
mu(G^-1 B) = [w(1)+w(2)+w(1)+w(2)]/(2*3^J) = 4/(2*3^J) = (4/3)*3/(2*3^J) = mu(B).
With 3b+1 = 7 mod 9 (b = 2 mod 3): candidates 7,8,4,2 for k=0,1,2,3, valid only for k=0 (a=7) and
k=1 (a=8), so mu(G^-1 B) = [w(1)+w(2)]/(2*3^J) = 2/(2*3^J) = (2/3)*3/(2*3^J) = mu(B).   QED
Equivalently w is the Perron left eigenvector of the branch matrix [[2,1],[2,1]] (class 1 mod 3 has
three lifts mod 9 going to classes 1,1,2 mod 3; so does class 2): (w1,w2)[[2,1],[2,1]] = 3 (w1,w2)
gives w1 = 2 w2.  G is a 6-branch piecewise-affine 3-adic expanding map, 4-to-1 over 1+3Z_3 and
2-to-1 over 2+3Z_3.""")

for J in range(1, 8):
    MJ, M1 = 3 ** J, 3 ** (J + 1)
    acc = {}
    for a in range(M1):
        if a % 3 == 0:
            continue
        b = G(a)[0] % MJ
        acc[b] = acc.get(b, 0) + Fr(4, 3) * (1 if a % 3 == 1 else Fr(1, 2))
    for b in range(MJ):
        if b % 3 == 0:
            continue
        w_b = Fr(4, 3) if b % 3 == 1 else Fr(2, 3)
        check(acc[b] == 3 * w_b, "invariance at J=%d b=%d" % (J, b))
    print("  J=%d: sum_{a mod 3^%d : G(a)=b mod 3^%d} w(a) = 3 w(b) for all %d unit classes b  [checked]" % (J, J + 1, J, 2 * 3 ** (J - 1)))
print("FINITE-EXACT (S2.2): mu-invariance verified on all cylinders of level J<=7.")

print("""THEOREM S2.3 (PROVED).  Stationary law and ergodicity at every level.
(i)  pi_J(b) := mu(b+3^J Z_3) = w(b mod 3)/(2*3^(J-1)) is P_J-invariant (this is S2.1 read at level J).
     At J=2: pi_2 = 2/9 on each of 1,4,7 and 1/9 on each of 2,5,8  (= F1).
(ii) Support of P_J^n(a,.) for n < J is EXACTLY the set of 3^n unit classes c = G^n(a) mod 3^(J-n)
     (induction on n with (b): one step from the class of x mod 3^(J-n) fills the class of G(x) mod 3^(J-n-1)).
     For n >= J every entry of P_J^n is positive: after J-1 steps a full class mod 3 is reached, and
     the three sub-classes mod 9 of a class mod 3 have G-images covering both classes mod 3 at every
     depth (G maps each class mod 9 onto a full class mod 3).  So P_J is primitive, hence irreducible
     and aperiodic, with unique stationary law pi_J, for every J>=1.
(iii) Drift: E_pi[k] = 2*(2/9) + 1*(1/9) + 0 + 3*(1/9) + 0 + 1*(1/9) = 1, so
     E_mu[log(2^k/3)] = log 2 - log 3 = log(2/3) = -0.405465.
     The uniform-residue average would give E[k] = 7/6 and the WRONG drift (7/6)log2 - log3 = -0.2899.
(iv) Equivalent phrasing on E-arrows: one G-step reverses one 3n+1 arrow and k halvings; at
     stationarity halvings and triplings are used in ratio 1:1 (forward Collatz uses 2:1).""")

# exact stationary law at J=2 and drift
pi2 = {a: (Fr(2, 9) if a % 3 == 1 else Fr(1, 9)) for a in (1, 2, 4, 5, 7, 8)}
for b in (1, 2, 4, 5, 7, 8):
    s = sum(pi2[a] * P2[a].get(b, 0) for a in pi2)
    check(s == pi2[b], "pi_2 invariant at b=%d" % b)
Ek = sum(pi2[a] * KT[a] for a in pi2)
Ek_unif = Fr(sum(KT.values()), 6)
check(Ek == 1, "E_pi[k] = 1")
print("  E_pi[k] = %s (exact);  uniform-residue E[k] = %s;  drift log(2/3) = %.6f;  naive %.6f"
      % (Ek, Ek_unif, log(2 / 3), float(Ek_unif) * log(2) - log(3)))

# ergodicity/reachability structure enumerated for J<=7
for J in range(1, 8):
    MJ = 3 ** J
    units_J = [a for a in range(MJ) if a % 3]
    succ = {}
    for a in units_J:
        succ[a] = frozenset(G(a + t * MJ)[0] % MJ for t in range(3))
    allset = frozenset(units_J)
    ok_struct = True
    for a in units_J:
        reach = frozenset([a])
        x = a
        for n in range(1, J + 1):
            reach = frozenset().union(*(succ[y] for y in reach))
            x = G(x)[0]
            if n < J:
                mod = 3 ** (J - n)
                pred = frozenset(c for c in units_J if c % mod == x % mod)
                if reach != pred:
                    ok_struct = False
            else:
                if reach != allset:
                    ok_struct = False
    check(ok_struct, "reachability structure at J=%d" % J)
    print("  J=%d: supp P_J^n(a,.) = {c = G^n(a) mod 3^(J-n)} for n<J and = all %d classes for n=J  [checked all a]" % (J, len(units_J)))
print("FINITE-EXACT (S2.4): primitivity with exact support law verified at levels J<=7.")

# ============================================================================
banner("S3  Terras-type theorem for the greedy 3-adic stopping time; exact d_J for J<=12")
# ============================================================================
print("""Notation.  Word k_1 k_2 ... of m; K_i = k_1+...+k_i;  G^i(m) = (2^(K_i) m - B_i)/3^i with
    B_i = sum_{s=1}^{i} 3^(s-1) 2^(K_i - K_s)  > 0   (PROVED by induction: G^(i+1) = (2^(k_(i+1)) G^i - 1)/3).
sigma(m)     := min{i>=1 : G^i(m) < m}           (real stopping time; sigma(1) = infinity),
sigma_res(m) := min{i>=1 : 2^(K_i) < 3^i}        (residue stopping time).

THEOREM S3.1 (PROVED).
(a) sigma(m) <= sigma_res(m) for every m>=2: if 2^(K_i) < 3^i then G^i(m) < 2^(K_i) m / 3^i < m.
(b) If sigma(m) = i < sigma_res(m) then m (2^(K_i) - 3^i) < B_i with 2^(K_i) - 3^i >= 1, hence m < B_i.
    Since the word of length J is a function of m mod 3^(J+1) (S1.1b) and there are finitely many words
    of length J, sup{B_i : i<=J} =: C_J < infinity, and {m : sigma(m)<=J} and {m : sigma_res(m)<=J}
    differ only inside [1, C_J].  Thus sigma(m)<=J is decided by m mod 3^(J+1) for all m > C_J,
    and the two sets have the same natural density.
(c) {m : 3 not| m, sigma_res(m) <= J} is a union of unit classes mod 3^(J+1); the natural density of a
    unit class mod 3^(J+1) relative to the non-multiples of 3 is 1/(2*3^J) = its normalized Haar mass.
    Hence d_J := dens{sigma<=J} = dens{sigma_res<=J} = Haar{sigma_res<=J} = (#good classes)/(2*3^J).
(d) d_J -> 1.  Route 1 (elementary, with rate; inherited from the wave-one lane S4.5): under Haar the
    class-mod-3 chain with weights E[2^k | class] gives the tilted matrix M = [[5/3,1/3],[10/3,2/3]]
    (rank 1, eigenvalue 7/3), so E_Haar[2^(K_J)] = 3 (7/3)^(J-1) and by Markov's inequality
    1 - d_J <= P(2^(K_J) >= 3^J) <= (7/9)^(J-1).
    Route 2 (ergodic): P_2 is primitive (S2.3), so the strong law for ergodic finite Markov chains
    [CITED: Norris, Markov Chains, Thm 1.10.2] gives K_n/n -> E_pi[k] = 1 for mu-a.e. m; since
    2/3 <= dmu/dHaar <= 4/3, the same holds Haar-a.e., so 2^(K_n) < 3^n eventually, i.e. sigma_res < infinity
    Haar-a.e.; d_J = Haar{sigma_res<=J} increases to Haar{sigma_res<infinity} = 1.  Route 2 uses
    the cited LLN; Route 1 is complete and elementary.  No step is OPEN.
(e) OPEN (as for Collatz): stopping time -> convergence to 1 needs 'no positive G-cycle except {1}
    and no divergent positive G-orbit'; FINITE-EXACT to 10^6 (wave one) and to 10^7 here (S5).""")

# exact enumeration mod 3^(J+1) for J<=12 (vectorised, exact integers)
JMAX = 12
KMAX_TAB = [0] + [((3 ** i).bit_length() - 1) for i in range(1, 400)]  # 2^K < 3^i  iff  K <= KMAX_TAB[i]
KMAX_ARR = np.array(KMAX_TAB, dtype=np.int64)
t1 = time.time()
mod = 3 ** (JMAX + 1)
r = np.arange(mod, dtype=np.int64)
r = r[r % 3 != 0]
n_units = len(r)
K = np.zeros(n_units, dtype=np.int64)
done = np.zeros(n_units, dtype=bool)
d_list = []
g_list = []
inherited_f = {1: Fr(2, 3), 2: Fr(8, 9), 3: Fr(25, 27), 4: Fr(26, 27), 5: Fr(236, 243), 6: Fr(239, 243),
               7: Fr(241, 243), 8: Fr(2173, 2187), 9: Fr(19609, 19683), 10: Fr(58868, 59049)}
print("\n  J | #unit classes mod 3^(J+1) | d_J exact | float | 1-d_J | (7/9)^(J-1) | bound holds | g_J=P(2^K_J<3^J) | wave-one f_J")
for i in range(1, JMAX + 1):
    k = KT_ARR[r % 9]
    r = ((r << k) - 1) // 3
    K += k
    cur_mod = 3 ** (JMAX + 1 - i)
    r %= cur_mod
    desc_now = K <= KMAX_TAB[i]
    done |= desc_now
    cnt = int(done.sum())
    dJ = Fr(cnt, 2 * 3 ** i)
    gJ = Fr(int(desc_now.sum()), 2 * 3 ** i)
    # every class mod 3^(i+1) is 3^(JMAX-i) times over-represented among the n_units classes mod 3^(JMAX+1)
    check((cnt * (2 * 3 ** i)) % n_units == 0, "class over-representation is exact at J=%d" % i)
    dJ = Fr(cnt, n_units)
    gJ = Fr(int(desc_now.sum()), n_units)
    d_list.append(dJ)
    g_list.append(gJ)
    bound = Fr(7, 9) ** (i - 1)
    check(1 - dJ <= bound, "Markov bound at J=%d" % i)
    if i in inherited_f:
        check(dJ == inherited_f[i], "d_J agrees with wave-one f_J at J=%d" % i)
    check(dJ >= (d_list[-2] if len(d_list) > 1 else 0), "monotone d_J")
    print("  %2d | %8d | %s | %.6f | %s | %.6f | %s | %s | %s"
          % (i, 2 * 3 ** i, dJ, float(dJ), 1 - dJ, float(bound), 1 - dJ <= bound, gJ,
             inherited_f.get(i, "(new)")))
print("FINITE-EXACT (S3.2): d_J for J<=12 (J=11,12 new; J<=10 equal to the wave-one f_J), monotone, 1-d_J <= (7/9)^(J-1).  time %.1fs" % (time.time() - t1))

# exceptional set check: sigma vs sigma_res on m <= 10^6 (vectorised, exact)
t1 = time.time()
N_SIG = 10 ** 6
m0 = np.arange(1, N_SIG + 1, dtype=np.int64)
m0 = m0[m0 % 3 != 0]
v = m0.copy()
Kc = np.zeros(len(m0), dtype=np.int64)
sig = np.zeros(len(m0), dtype=np.int64)
sigres = np.zeros(len(m0), dtype=np.int64)
active = np.ones(len(m0), dtype=bool)
active[m0 == 1] = False
step = 0
while active.any() and step < 400:
    step += 1
    idx = np.nonzero(active)[0]
    vv = v[idx]
    k = KT_ARR[vv % 9]
    vv = ((vv << k) - 1) // 3
    v[idx] = vv
    Kc[idx] += k
    newres = (Kc[idx] <= KMAX_TAB[step]) & (sigres[idx] == 0)
    sigres[idx[newres]] = step
    newdesc = (vv < m0[idx]) & (sig[idx] == 0)
    sig[idx[newdesc]] = step
    fin = (sig[idx] > 0) & (sigres[idx] > 0)
    active[idx[fin]] = False
check(not active.any(), "all m<=10^6 (m>1) have finite sigma and sigma_res within 400 steps")
mism = np.nonzero((sig != sigres) & (m0 > 1))[0]
check(np.all(sig[m0 > 1] <= sigres[m0 > 1]), "sigma <= sigma_res")
print("\n  m<=10^6, 3 not|m, m>1: max sigma = %d (at m=%d), max sigma_res = %d (at m=%d)"
      % (sig[m0 > 1].max(), m0[m0 > 1][sig[m0 > 1].argmax()], sigres[m0 > 1].max(), m0[m0 > 1][sigres[m0 > 1].argmax()]))
print("  #{m : sigma(m) != sigma_res(m)} = %d;  witnesses (m, sigma, sigma_res): %s"
      % (len(mism), [(int(m0[i]), int(sig[i]), int(sigres[i])) for i in mism[:20]]))
for i in mism:
    check(sig[i] < sigres[i], "mismatch direction")
    # m < B_{sigma} must hold: verify by recomputing B_i exactly
    m = int(m0[i]); i_s = int(sig[i])
    # exact B_i via G^i(m) = (2^K m - B)/3^i
    x = m; Kw = 0
    for s in range(i_s):
        x, kk = G(x); Kw += kk
    B_i = (1 << Kw) * m - (3 ** i_s) * x
    check(B_i > 0 and m < B_i, "exceptional m=%d is below B_sigma=%d" % (m, B_i))
print("  every mismatch has sigma < sigma_res and m < B_sigma (S3.1b)  [checked];  time %.1fs" % (time.time() - t1))
print("FINITE-EXACT (S3.3): the exceptional set {sigma != sigma_res} inside [1,10^6] is as listed; empty means d_J is exact on this range.")

# Collatz comparison: Terras densities mod 2^J
t1 = time.time()
JC = 20
modc = 1 << JC
n = np.arange(modc, dtype=np.int64)
a = np.zeros(modc, dtype=np.int64)
donec = np.zeros(modc, dtype=bool)
# 3^a < 2^i  iff  a <= AMAX[i]  where AMAX[i] = max a with 3^a < 2^i
AMAX = [max([aa for aa in range(0, 2 * JC + 2) if 3 ** aa < 2 ** i] or [-1]) for i in range(0, JC + 1)]
print("\n  Collatz T(n) = n/2 (even), (3n+1)/2 (odd); Terras residue stopping time: first i with 3^(a_i) < 2^i")
print("  J | Terras density F(J) (exact, residues mod 2^J) | float | greedy 3-adic d_J")
terras = []
for i in range(1, JC + 1):
    odd = (n & 1) == 1
    n = np.where(odd, (3 * n + 1) >> 1, n >> 1)
    a += odd
    n &= (1 << (JC - i)) - 1
    donec |= (a <= AMAX[i])
    F = Fr(int(donec.sum()), modc)
    terras.append(F)
    dj = d_list[i - 1] if i <= JMAX else None
    print("  %2d | %s | %.6f | %s" % (i, F, float(F), (str(dj) + " = %.6f" % float(dj)) if dj else "-"))
check(terras[0] == Fr(1, 2) and terras[1] == Fr(3, 4), "Terras F(1)=1/2, F(2)=3/4")
check(all(terras[i] >= terras[i - 1] for i in range(1, len(terras))), "Terras monotone")
print("FINITE-EXACT (S3.4): Terras densities to J=20 (Terras 1976 theorem F(J)->1 is CITED; values computed here).")
print("  Comparison: the 3-adic greedy stopping set reaches 0.99 by J=7 and 0.9990 by J=12, the 2-adic Collatz")
print("  one is only %.4f at J=12 and %.4f at J=20: drift log(2/3)=-0.405 per G-step vs (1/2)log(3/4)=-0.144 per T-step;" % (float(terras[11]), float(terras[19])))
print("  and the 3-adic carry B_i>0 HELPS descent (sigma<=sigma_res), whereas the 2-adic carry hinders it.  time %.1fs" % (time.time() - t1))

print("""THEOREM S3.5 (PROVED).  Exact exponential moments for EVERY tilt; the sharp Chernoff exponent.
Let u = e^theta > 0 and lump the level-2 chain to the class c_n = G^n(m) mod 3 in {1,2}.  Given c_n, the
residue X_n is uniform on the three residues of its class (S1.1c), and the pairs (k(X_n), c_(n+1)) are
   class 1: (2,1),(0,1),(0,2)        class 2: (1,1),(3,1),(1,2)   [= class 1 shifted by one halving].
So the tilted matrix M_u(c,c') = E[u^k ; c_(n+1)=c' | c_n=c] is
   M_u = (1/3) [[u^2+1, 1], [u+u^3, u]] = v w^T,  v=(1,u), w=((u^2+1)/3, 1/3),
of rank ONE for every u, with eigenvalue rho(u) = w.v = (u^2+u+1)/3.  With X_0 Haar (c_0 = 1,2 each 1/2):
   E_Haar[u^(K_J)] = (1/2,1/2) M_u^J (1,1)^T = (u+1)(u^2+2)/6 * ((u^2+u+1)/3)^(J-1)      (exact, all J>=1).
u=2 gives 3*(7/3)^(J-1), the wave-one identity (its 'initial vector' (5/2,1/2) = (1/2,1/2)M_2).
Consequences.  (a) Lambda(theta) := lim (1/J) log E[e^(theta K_J)] = log((u^2+u+1)/3) exactly, so by the
Gartner-Ellis theorem / the LDP for additive functionals of a finite irreducible chain
[CITED: Dembo-Zeitouni, Large Deviations Techniques and Applications, Sec. 3.1 and Thm 2.3.6]
P(K_J >= cJ) = exp(-J I(c) + o(J)) with I(c) = sup_theta (theta c - Lambda(theta)), c > Lambda'(0) = 1.
(b) The residue-stopping tail is exactly this event at c = log_2 3: 1 - g_J = P(2^(K_J) >= 3^J).
The wave-one bound (7/9)^(J-1) is Chernoff at u = 2 (per-step factor rho(2)/3 = 7/9); it is valid but not
sharp: the optimizer solves Lambda'(theta) = c, i.e. u(2u+1)/(u^2+u+1) = log_2 3, while at u=2 the slope is
10/7 = 1.4286 < 1.5850.  The sharp per-step factor is exp(-I(log_2 3)) = rho(u*)/u*^(log_2 3) (below).
(c) Since {sigma_res <= J} contains {2^(K_i) < 3^i for i=J} , 1 - d_J <= 1 - g_J and the prefix tail decays
at least at the sharp rate; the observed ratios (1-d_(J+1))/(1-d_J) for J<=12 are printed for comparison.""")
# exact moment identity, enumerated over units mod 3^(J+1), for u in {2,3,5,7} and J<=8
for u in (2, 3, 5, 7):
    for J in (1, 2, 3, 5, 8):
        M1 = 3 ** (J + 1)
        tot = 0
        cnt = 0
        for a in range(M1):
            if a % 3 == 0:
                continue
            x = a
            K = 0
            for i in range(J):
                x, k = G(x)
                K += k
            tot += u ** K
            cnt += 1
        lhs = Fr(tot, cnt)
        rhs = Fr((u + 1) * (u * u + 2), 6) * Fr(u * u + u + 1, 3) ** (J - 1)
        check(lhs == rhs, "E[u^K_J] exact at u=%d J=%d" % (u, J))
    Mu = [[Fr(u * u + 1, 3), Fr(1, 3)], [Fr(u + u ** 3, 3), Fr(u, 3)]]
    check(Mu[0][0] * Mu[1][1] - Mu[0][1] * Mu[1][0] == 0, "det M_u = 0 at u=%d" % u)
    check(Mu[0][0] + Mu[1][1] == Fr(u * u + u + 1, 3), "trace = rho(u)")
    print("  u=%d: E_Haar[u^K_J] = %s * (%s)^(J-1) verified by enumeration at J in {1,2,3,5,8}; det M_u = 0, rho = %s"
          % (u, Fr((u + 1) * (u * u + 2), 6), Fr(u * u + u + 1, 3), Fr(u * u + u + 1, 3)))
c_thr = log(3) / log(2)
slope2 = Fr(2 * 5, 7)
check(slope2 == Fr(10, 7) and float(slope2) < c_thr, "Lambda'(log 2) = 10/7 < log_2 3: u=2 is not the optimizer")
# optimal u*: (2-c) u^2 + (1-c) u - c = 0
import math
disc = (1 - c_thr) ** 2 + 4 * c_thr * (2 - c_thr)
u_star = ((c_thr - 1) + math.sqrt(disc)) / (2 * (2 - c_thr))
rho_star = (u_star * u_star + u_star + 1) / 3
fac_star = rho_star / u_star ** c_thr
fac_2 = (7 / 3) / 3
check(abs(u_star * (2 * u_star + 1) / (u_star ** 2 + u_star + 1) - c_thr) < 1e-12, "u* solves Lambda' = c")
check(fac_star < fac_2, "sharp factor below 7/9")
# convexity check: the factor rho(u)/u^c is minimized at u*
for uu in (1.5, 2.0, 2.5, 3.0, 3.5):
    check((uu * uu + uu + 1) / 3 / uu ** c_thr >= fac_star - 1e-12, "u* is the minimizer (sampled)")
print("  c = log_2 3 = %.6f;  Lambda'(log 2) = 10/7 = %.4f;  u* = %.6f;  rho(u*) = %.6f;" % (c_thr, 10 / 7, u_star, rho_star))
print("  sharp per-step factor exp(-I(c)) = rho(u*)/u*^c = %.6f  vs wave-one Chernoff 7/9 = %.6f  [checked u* minimizes]" % (fac_star, fac_2))
ratios_d = [float((1 - d_list[i + 1]) / (1 - d_list[i])) for i in range(len(d_list) - 1)]
ratios_g = [float((1 - g_list[i + 1]) / (1 - g_list[i])) for i in range(len(g_list) - 1)]
print("  observed (1-d_(J+1))/(1-d_J), J=1..11: %s" % ["%.3f" % r for r in ratios_d])
print("  observed (1-g_(J+1))/(1-g_J), J=1..11: %s" % ["%.3f" % r for r in ratios_g])
gm_d = (float(1 - d_list[11]) / float(1 - d_list[7])) ** 0.25
print("  geometric-mean ratio of 1-d_J over J=8..12: %.4f (finite-J, lattice effects; asymptotic rate is %.4f)" % (gm_d, fac_star))
print("PROVED (S3.6): exact moment identity for every tilt (rank-one M_u), exact LD rate function I(c) = sup(theta c - log((u^2+u+1)/3));")
print("  the wave-one (7/9)^(J-1) is the u=2 Chernoff bound, valid but not the sharp exponent (%.4f)." % fac_star)

print("""THEOREM S3.7 (PROVED per level).  Exact threshold for sigma = sigma_res.  Fix i>=1 and a word w of length i
realized by a unit class a mod 3^(i+1), with 2^(K_i) > 3^i.  On that class G^i is the affine map
m -> (2^(K_i) m - B_i)/3^i, so G^i(m) < m  iff  m < m*(w) := B_i/(2^(K_i) - 3^i).  Hence a mismatch
sigma(m) = i < sigma_res(m) forces (1) m < m*(w), m in the class a, and (2) no earlier descent: G^l(m) >= m for
all l < i (then G^i(m) < m automatically, so sigma(m) = i, and 2^(K_i) > 3^i gives sigma_res(m) > i).  Both
conditions are decidable by finite enumeration: the candidates m < m*(w) in each class are finitely many
(m* <= 145 for i <= 12).  If no candidate other than m = 1 (word 2^i, m* = 1 exactly) survives (2), then
sigma(m) = sigma_res(m) for every m >= 2 with sigma(m) <= 12, and d_J is the EXACT natural density of
{sigma <= J} with no exceptional set for J <= 12.  Condition (1) alone is NOT sufficient (candidates that
already descended earlier are listed below); the table gives both counts.""")
CANDS = []
t1 = time.time()
JT = 12
mod = 3 ** (JT + 1)
r = np.arange(mod, dtype=np.int64)
r = r[r % 3 != 0]
a0 = r.copy()
K = np.zeros(len(r), dtype=np.int64)
Bv = np.zeros(len(r), dtype=object)  # exact ints
Bv[:] = 0
worst = []
n_bad_total = 0
for i in range(1, JT + 1):
    k = KT_ARR[r % 9]
    r = ((r << k) - 1) // 3
    K += k
    # B_(i) = 2^k B_(i-1) + 3^(i-1)   (from G^(i) = (2^k G^(i-1) - 1)/3)
    Bv = Bv * (1 << 0)
    Bv = np.array([int(b) * (1 << int(kk)) + 3 ** (i - 1) for b, kk in zip(Bv, k)], dtype=object)
    # classes mod 3^(i+1) are represented 3^(JT-i) times; pick the representatives a0 < 3^(i+1)
    sel = a0 < 3 ** (i + 1)
    grow = sel & (K > KMAX_TAB[i])  # 2^K > 3^i
    idx = np.nonzero(grow)[0]
    max_ratio = Fr(0)
    arg = None
    n_cand = 0
    n_bad = 0
    for t in idx:
        a = int(a0[t]); Kt = int(K[t]); Bt = int(Bv[t])
        ratio = Fr(Bt, (1 << Kt) - 3 ** i)
        if ratio > max_ratio:
            max_ratio, arg = ratio, (a, Kt, Bt)
        if a == 1:
            check(ratio == 1 and Kt == 2 * i, "m=1 class: word 2^i, m* = 1")
        # candidates: integers m >= 2 in the class a mod 3^(i+1) with m < m*(w)  (condition (1))
        mm = a
        while mm < ratio:
            if mm >= 2:
                n_cand += 1
                CANDS.append((i, mm, float(ratio)))
                # condition (2): no earlier descent
                x = mm
                early = False
                for l in range(1, i):
                    x = G(x)[0]
                    if x < mm:
                        early = True
                        break
                if not early:
                    x = G(x)[0]
                    check(x < mm, "affine criterion: G^i(m) < m iff m < m*")
                    n_bad += 1
            mm += 3 ** (i + 1)
    n_bad_total += n_bad
    worst.append((i, max_ratio, arg, n_cand, n_bad))
    print("  i=%2d: classes with 2^K>3^i: %6d; max m*(w) = %s = %.4f at (a=%d,K=%d,B=%d); candidates m>=2 below m*: %3d; true mismatches: %d"
          % (i, len(idx), max_ratio, float(max_ratio), arg[0], arg[1], arg[2], n_cand, n_bad))
check(n_bad_total == 0, "no mismatch at any level i<=12")
print("  candidates (level i, m, m*(w)) -- each descended at an earlier step: %s" % CANDS)
for (i_c, m_c, _) in CANDS:
    x = m_c
    first = None
    for l in range(1, i_c):
        x = G(x)[0]
        if x < m_c:
            first = l
            break
    check(first is not None and first < i_c, "candidate m=%d at level %d descended earlier" % (m_c, i_c))
print("PROVED (S3.8): sigma(m) = sigma_res(m) for all m >= 2 with sigma(m) <= 12 (exact threshold enumeration, %.1fs);" % (time.time() - t1))
print("  FINITE-EXACT (S3.3) extends this to all m <= 10^6 (where sigma_res <= 31).  OPEN: all m, all levels.")
del r, a0, K, Bv

# ============================================================================
banner("S4  G_k for odd k, 3 not | k: E[j]=1 for every k; cycle gate; census |k|<=49; duality with T_k")
# ============================================================================
print("""G_k(m) = (2^j m - k)/3, j>=0 minimal with 2^j m in {k+3, k+6} mod 9; domain: all m with 3 not| m.
THEOREM S4.1 (PROVED).  For every k with 3 not| k: the transfer law S1.1 holds verbatim (the map on
a class mod 9 is affine with unit slope 2^j/3 of 3-adic size 3), the invariant density is
w = 4/3 on the class m = k mod 3 and 2/3 on the class m = -k mod 3, and E_pi[j] = 1.
Proof.  Write T_1 = k+3 (result 1 mod 3), T_2 = k+6 (result 2 mod 3) and u = the discrete log with
a = T_1 2^(-u) mod 9 (u in Z/6; 2 generates (Z/9)^x).  T_2/T_1 = (k+6)/(k+3) is 4 = 2^2 if k = 1 mod 3
and 7 = 2^4 if k = 2 mod 3.  Case k = 1 mod 3: j(a) = min(u, u+2 mod 6): u=0..3 -> j=u, target T_1;
u=4,5 -> j=0,1, target T_2; and a mod 3 = (k+3)(-1)^u = (-1)^u.  This is EXACTLY the k=1 table
(u=0,1,2,3,4,5 <-> residues 4,2,1,5,7,8), so the chain, w and E[j]=1 are those of S2.  Case k = 2
mod 3: G_(-k)(-m) = -G_k(m) conjugates to the case -k = 1 mod 3 and swaps the classes mod 3.  QED""")

for kk in (1, 2, 4, 5, 7, 8):
    tab = jtab_k(kk)
    Gk = Gk_factory(kk)
    Pk = {}
    for a in (1, 2, 4, 5, 7, 8):
        row = {}
        for t in range(3):
            b = Gk(a + 9 * t)[0] % 9
            row[b] = row.get(b, 0) + Fr(1, 3)
        Pk[a] = row
    wk = {1: Fr(4, 3) if kk % 3 == 1 else Fr(2, 3), 2: Fr(2, 3) if kk % 3 == 1 else Fr(4, 3)}
    pik = {a: wk[a % 3] / 6 for a in (1, 2, 4, 5, 7, 8)}
    for b in pik:
        check(sum(pik[a] * Pk[a].get(b, 0) for a in pik) == pik[b], "stationary law for k=%d mod 9" % kk)
    Ej = sum(pik[a] * tab[a] for a in pik)
    check(Ej == 1, "E[j]=1 for k=%d mod 9" % kk)
    check(sorted(tab.values()) == [0, 0, 1, 1, 2, 3], "j multiset for k=%d" % kk)
    print("  k = %d mod 9: j-table %s  w(1 mod 3)=%s w(2 mod 3)=%s  E_pi[j]=%s" % (kk, tab, wk[1], wk[2], Ej))
print("FINITE-EXACT+PROVED (S4.2): all six k mod 9 have stationary E[j]=1, drift log(2/3), j-multiset {0,0,1,1,2,3}.")

print("""THEOREM S4.3 (PROVED, cycle gate).  If m_0 -> m_1 -> ... -> m_L = m_0 is a G_k-cycle with word
(j_1..j_L), J = sum j_i, J_i = j_1+..+j_i, then
    m_0 (2^J - 3^L) = k B',   B' = sum_{i=1}^{L} 3^(i-1) 2^(J - J_i)  =  B(j_L, ..., j_1),
where B(w) = sum_{i=0}^{L-1} 3^(L-1-i) 2^(K_i) is the inherited T_k gate polynomial
(T_k cycle: n_0 (2^K - 3^L) = k B(k_1..k_L), arithmetic_braids_20260917_collatz.md / F3).
So the two gates coincide up to REVERSING the word.  Consequences (PROVED):
 (i) sign(m_0) = sign(k) * sign(2^J - 3^L): G_k-cycles with 2^J < 3^L (the typical, drift-side words)
     have sign -sign(k); cycles with sign(k) need the atypical 2^J > 3^L.
 (ii) A G_k-cycle is (reversed) a T_k-cycle iff all its elements are odd (then j_i = v_2(3 m_i + k)).
     A T_k-cycle is (reversed) a G_k-cycle iff every exponent is the greedy minimum for its node.
 (iii) G_(-k)(-m) = -G_k(m) and T_(-k)(-n) = -T_k(n): cycle sets of k and -k are negatives.
There is NO word-level bijection between G_k-cycles and T_k- or T_(-k)-cycles: the census below has
different counts, and the two maps read different digit expansions (3-adic vs 2-adic).""")


def census_G(k, M=10 ** 5, ESC=10 ** 18):
    Gk = Gk_factory(k)
    state = {}
    cycles = []
    starts = [m for m in range(-M, M + 1) if m != 0 and m % 3 != 0]
    starts.sort(key=abs)
    n_escape = 0
    for m in starts:
        if m in state:
            continue
        path = []
        pos = {}
        x = m
        cid = None
        while True:
            if x in state:
                cid = state[x]
                break
            if x in pos:
                cyc = path[pos[x]:]
                cycles.append(cyc)
                cid = len(cycles) - 1
                break
            if abs(x) > ESC:
                cid = -1
                n_escape += 1
                break
            pos[x] = len(path)
            path.append(x)
            x = Gk(x)[0]
        for y in path:
            state[y] = cid
    return cycles, n_escape


def census_T(k, M=2 * 10 ** 5, ESC=10 ** 18):
    state = {}
    cycles = []
    starts = [n for n in range(-M, M + 1) if n % 2 != 0]
    starts.sort(key=abs)
    n_escape = 0
    for n0 in starts:
        if n0 in state:
            continue
        path = []
        pos = {}
        x = n0
        cid = None
        while True:
            if x in state:
                cid = state[x]
                break
            if x in pos:
                cycles.append(path[pos[x]:])
                cid = len(cycles) - 1
                break
            if abs(x) > ESC:
                cid = -1
                n_escape += 1
                break
            pos[x] = len(path)
            path.append(x)
            y = 3 * x + k
            if y == 0:
                cid = -2
                break
            while y % 2 == 0:
                y //= 2
            x = y
        for y in path:
            state[y] = cid
    return cycles, n_escape


def gate_check_G(cyc, k):
    Gk = Gk_factory(k)
    L = len(cyc)
    js = []
    for i in range(L):
        nxt, j = Gk(cyc[i])
        check(nxt == cyc[(i + 1) % L], "cycle consistency")
        js.append(j)
    J = sum(js)
    Ji = 0
    Bp = 0
    for i in range(1, L + 1):
        Ji += js[i - 1]
        Bp += 3 ** (i - 1) * 2 ** (J - Ji)
    check(cyc[0] * (2 ** J - 3 ** L) == k * Bp, "G_k gate for cycle %s k=%d" % (cyc, k))
    # reversed-word identity with the T gate polynomial
    rev = js[::-1]
    Ki = 0
    Bt = 0
    for i in range(L):
        Bt += 3 ** (L - 1 - i) * 2 ** Ki
        Ki += rev[i]
    check(Bt == Bp, "B'(w) = B(reversed w)")
    return js, J, Bp


def gate_check_T(cyc, k):
    L = len(cyc)
    ks = []
    for i in range(L):
        y = 3 * cyc[i] + k
        e = 0
        while y % 2 == 0:
            y //= 2
            e += 1
        check(y == cyc[(i + 1) % L], "T cycle consistency")
        ks.append(e)
    Kt = sum(ks)
    Ki = 0
    B = 0
    for i in range(L):
        B += 3 ** (L - 1 - i) * 2 ** Ki
        Ki += ks[i]
    check(cyc[0] * (2 ** Kt - 3 ** L) == k * B, "T_k gate")
    return ks, Kt, B


t1 = time.time()
KS = [k for k in range(-49, 50) if k % 2 and k % 3]
resG = {}
resT = {}
for k in KS:
    resG[k] = census_G(k)
    resT[k] = census_T(k)
print("  census time %.1fs" % (time.time() - t1))
F3 = {1: 4, 5: 9, 7: 5, 11: 7, 13: 13}
print("\n  k | #G_k cycles (|m|<=1e5, both signs) | #T_k cycles (odd |n|<=2e5, both signs) | #common (as sets) | G_k cycle minima (by |.|) | T_k cycle minima")
common_tab = {}
for k in KS:
    cG, eG = resG[k]
    cT, eT = resT[k]
    check(eG == 0 and eT == 0, "no escapes for k=%d" % k)
    setsG = [frozenset(c) for c in cG]
    setsT = [frozenset(c) for c in cT]
    common = [c for c in setsG if c in setsT]
    common_tab[k] = common
    for c in cG:
        gate_check_G(c, k)
    for c in cT:
        gate_check_T(c, k)
    # (ii): G-cycle is a T-cycle iff all odd
    for c in cG:
        allodd = all(x % 2 for x in c)
        check(allodd == (frozenset(c) in setsT), "criterion (ii) G->T for k=%d cycle %s" % (k, c))
    # (ii): T-cycle is a G-cycle iff greedy exponents (reversed)
    tab = jtab_k(k)
    for c in cT:
        L = len(c)
        greedy = True
        for i in range(L):
            y = 3 * c[i] + k
            e = 0
            while y % 2 == 0:
                y //= 2
                e += 1
            # reversed arrow: c[i] = G_k(c[i+1]) requires exponent e to be greedy at node c[i+1]
            if tab[c[(i + 1) % L] % 9] != e:
                greedy = False
        check(greedy == (frozenset(c) in setsG), "criterion (ii) T->G for k=%d cycle %s" % (k, c))
    if abs(k) in F3:
        check(len(cT) == F3[abs(k)], "F3 count for k=%d" % k)
    mG = sorted([min(c, key=abs) for c in cG], key=abs)
    mT = sorted([min(c, key=abs) for c in cT], key=abs)
    print("  %3d | %2d | %2d | %d | %s | %s" % (k, len(cG), len(cT), len(common), mG, mT))
# negation symmetry
for k in KS:
    if k > 0:
        sG = set(frozenset(-x for x in c) for c in resG[k][0])
        check(sG == set(frozenset(c) for c in resG[-k][0]), "G negation symmetry k=%d" % k)
        sT = set(frozenset(-x for x in c) for c in resT[k][0])
        check(sT == set(frozenset(c) for c in resT[-k][0]), "T negation symmetry k=%d" % k)
print("  negation symmetry (iii) verified for all |k|<=49: cycles(G_-k) = -cycles(G_k), cycles(T_-k) = -cycles(T_k).")
print("  F3 counts for |k| in {1,5,7,11,13} reproduced: %s" % F3)
# specific checks
cG1 = [sorted(c) for c in resG[1][0]]
check(any(set(c) == {-11, -4} for c in resG[1][0]), "G_1 has the negative cycle {-11,-4}")
check(any(set(c) == {1} for c in resG[1][0]), "G_1 fixed point 1")
print("  G_1 cycles: %s  (the task's {-11,-4} check: PASS; word of (-4,-11) is (3,0), gate -4*(2^3-3^2) = 4 = B'(3,0))" % cG1)
posG1 = [c for c in resG[1][0] if min(c) > 0]
check(len(posG1) == 1, "only positive G_1 cycle is {1} for |m|<=1e5")
for k in KS:
    cG, _ = resG[k]
    for c in cG:
        js, J, Bp = gate_check_G(c, k)
        L = len(c)
        sgn = (1 if k > 0 else -1) * (1 if 2 ** J > 3 ** L else -1)
        check(sgn == (1 if c[0] > 0 else -1), "sign law (i) k=%d" % k)
print("  sign law (i) verified on every G_k cycle found.")
print("""THEOREM S4.5 (PROVED, scaling lemma).  For 3 not| k and 3 not| m:  G_k(k m) = k G_1(m), with the same j.
Proof: 2^j k m in {k+3, k+6} mod 9  <=>  2^j m in {1 + 3k^-1, 1 + 6k^-1} mod 9 = {4, 7} mod 9 (k^-1 = 1 or 2
mod 3 permutes the two targets), which is the G_1 condition; then (2^j k m - k)/3 = k (2^j m - 1)/3.  QED
Hence every G_k has the three UNIVERSAL cycles k*{1}, k*{-1}, k*{-4,-11} (scaled from G_1: 3 cycles, S4.4),
exactly as T_k(k n) = k T_1(n) gives T_k the four universal cycles k*{1}, k*{-1}, k*{-5,-7,-10..}, k*{-17,...}
(braids2 signed_cycles: content d = |b|/q; universal <=> q = 1).  PRIMITIVE cycles (not inside kZ) are the
ones with q = |k|; their counts are tabulated below (G: #cycles - 3; T: #cycles - 4).""")
for k in KS:
    Gk = Gk_factory(k)
    for m in range(-200, 201):
        if m % 3 == 0:
            continue
        gk, jk = Gk(k * m)
        g1, j1 = Gk_factory(1)(m)
        check(gk == k * g1 and jk == j1, "scaling lemma at k=%d m=%d" % (k, m))
    setsG = set(frozenset(c) for c in resG[k][0])
    for base in ({1}, {-1}, {-4, -11}):
        check(frozenset(k * x for x in base) in setsG, "universal cycle k*%s present for k=%d" % (sorted(base), k))
print("  scaling lemma verified for all |k|<=49, |m|<=200; universal G_k cycles k*{1}, k*{-1}, k*{-4,-11} present for every k.")
print("\n  k | #G_k | primitive G_k (= #-3) | #T_k | primitive T_k (= #-4) | primitive G_k cycle minima")
for k in KS:
    if k < 0:
        continue
    cG = resG[k][0]
    cT = resT[k][0]
    prim = [c for c in cG if any(x % k != 0 for x in c)] if k != 1 else []
    check(len(prim) == len(cG) - 3, "primitive count = total - 3 at k=%d" % k)
    check(len(cT) >= 4, "T_k has >= 4 cycles")
    print("  %3d | %2d | %2d | %2d | %2d | %s" % (k, len(cG), len(prim), len(cT), len(cT) - 4,
                                                sorted([min(c, key=abs) for c in prim], key=abs)))
print("FINITE-EXACT (S4.6): primitive-cycle table for 0<k<=49 (k<0 by negation); k=1: no primitive cycle on either side.")
# detail table for the G_k cycles with their words and gates, small k
print("\n  G_k cycles with words (k in {1,-1,5,-5,7,-7,11,13}):")
for k in (1, -1, 5, -5, 7, -7, 11, 13):
    for c in resG[k][0]:
        js, J, Bp = gate_check_G(c, k)
        print("    k=%3d  cycle %s  word %s  J=%d L=%d  2^J-3^L=%d  B'=%d  m0=k*B'/(2^J-3^L)=%s" %
              (k, c, js, J, len(c), 2 ** J - 3 ** len(c), Bp, Fr(k * Bp, 2 ** J - 3 ** len(c))))
print("FINITE-EXACT (S4.4): census |k|<=49 above.  Duality verdict: (a) PROVED exact gate identity with word")
print("  reversal; (b) PROVED negation conjugation k<->-k inside each family; (c) REFUTED: any bijection")
print("  G_k-cycles <-> T_k-cycles or <-> T_(-k)-cycles (counts differ already at k=1: #G_1=%d vs #T_1=%d);" % (len(resG[1][0]), len(resT[1][0])))
print("  (d) the exact intersection is given by criterion (ii) (common cycles counted above).")

# ============================================================================
banner("S5  Hostiles: 3^j+1, the exact peak bound, the worst m<=10^7, (8,5,1^n) family")
# ============================================================================
print("""THEOREM S5.1 (PROVED; wave-one S2.2 gives the word, here the exact orbit values).
For j>=1 and m = 3^j + 1:  G^i(m) = 4^i 3^(j-i) + 1 for 0<=i<=j-1 (each such value is 1 mod 9 when
j-i>=2, so k=2 and G = (4m-1)/3 = 4^(i+1) 3^(j-i-1) + 1), and G^(j-1)(m) = 3*4^(j-1)+1 = 4 mod 9,
so k=0 and G^j(m) = 4^(j-1).  Peak/m >= (3*4^(j-1)+1)/(3^j+1) -> infinity like (4/3)^(j-1):
a family with peak/m -> infinity, the exact 3-adic mirror of Collatz's 2^L-1 (T^L(2^L-1) = 3^L-1).""")
for j in range(1, 16):
    m = 3 ** j + 1
    x = m
    for i in range(j - 1):
        check(x == 4 ** i * 3 ** (j - i) + 1, "3^j+1 orbit value")
        x, k = G(x)
        check(k == 2 or (j - i - 1 == 0), "k=2 along 3^j+1")
    if j >= 2:
        check(x == 3 * 4 ** (j - 1) + 1, "penultimate value")
        x, k = G(x)
        check(k == 0 and x == 4 ** (j - 1), "G^j(3^j+1) = 4^(j-1)")
    # full orbit peak
    x = m; peak = m; steps = 0
    while x != 1 and steps < 10000:
        x, k = G(x); steps += 1; peak = max(peak, x)
    check(x == 1, "3^j+1 reaches 1")
    print("  j=%2d m=%9d  G^(j-1)(m)=%12d  ratio=%8.3f  (4/3)^(j-1)=%8.3f  full-orbit peak/m=%8.3f steps-to-1=%d"
          % (j, m, 3 * 4 ** (j - 1) + 1 if j >= 2 else m, (3 * 4 ** (j - 1) + 1) / m if j >= 2 else 1.0, (4 / 3) ** (j - 1), peak / m, steps))

print("""THEOREM S5.2 (PROVED, exact peak bound).  For every m (3 not| m) and every n>=0, along the word of m
    K_n <= 2n + [k_1 = 3],   hence   G^n(m) < 2^(K_n) m/3^n <= 2*(4/3)^n m,  and  <= (4/3)^n m unless m = 5 mod 9.
Proof.  Letters with k=3 are exactly the residues 5 mod 9.  From residue 5 the next residue is in {1,4,7}
(k in {2,0,0}); residue 5 is entered only from {7,8} (k in {0,1}): indeed G(5+9t) = 13+24t = 4+6t mod 9
is never 5, and G(a) = 5 mod 9 forces G(a) = 2 mod 3, i.e. a in {7,8} mod 9.  Pair every non-initial 3
with its predecessor letter (<= 1): the pairs are disjoint (a predecessor of a 3 is never itself a 3) and
each pair sums to <= 4 = 2*2; every unpaired letter is <= 2.  So K_n <= 2n, plus 1 if the first letter
is 3.  QED.   The family 3^j+1 attains K_(j-1) = 2(j-1), so the rate 4/3 per step is optimal, and no
family can grow faster than 2*(4/3)^n in n steps (Collatz's optimum is (3/2)^n from 2^L-1).""")

# vectorised full-orbit statistics for all m <= 10^7 (F1 verification) + the K_n <= 2n+[k1=3] bound
# processed in chunks to keep RAM small (each chunk ~1.7*10^6 starts)
t1 = time.time()
N_ALL = 10 ** 7
NCHUNK = 4
NU = 0
best_steps = (-1, 0)
best_ratio = (0.0, 0, 0)
bound_ok = True
for ci in range(NCHUNK):
    lo = ci * (N_ALL // NCHUNK) + 1
    hi = (ci + 1) * (N_ALL // NCHUNK)
    m0 = np.arange(lo, hi + 1, dtype=np.int64)
    m0 = m0[m0 % 3 != 0]
    NU += len(m0)
    v = m0.copy()
    peak = m0.copy()
    steps = np.zeros(len(m0), dtype=np.int32)
    Kc = np.zeros(len(m0), dtype=np.int32)
    first3 = (m0 % 9 == 5).astype(np.int32)
    active_idx = np.nonzero(m0 != 1)[0]
    step = 0
    while len(active_idx) and step < 1000:
        step += 1
        vv = v[active_idx]
        k = KT_ARR[vv % 9]
        vv = ((vv << k) - 1) // 3
        v[active_idx] = vv
        Kc[active_idx] += k.astype(np.int32)
        if np.any(Kc[active_idx] > 2 * step + first3[active_idx]):
            bound_ok = False
        peak[active_idx] = np.maximum(peak[active_idx], vv)
        steps[active_idx] = step
        active_idx = active_idx[vv != 1]
    check(step < 1000 and len(active_idx) == 0, "every m in chunk %d reaches 1" % ci)
    imax = int(steps.argmax())
    if int(steps[imax]) > best_steps[0]:
        best_steps = (int(steps[imax]), int(m0[imax]))
    ratio = peak / m0
    irat = int(ratio.argmax())
    if float(ratio[irat]) > best_ratio[0]:
        best_ratio = (float(ratio[irat]), int(m0[irat]), int(peak[irat]))
    del m0, v, peak, steps, Kc, first3, ratio, active_idx
check(bound_ok, "K_n <= 2n + [m=5 mod 9] along all orbits m<=10^7")
print("  all %d non-multiples of 3 in [1,10^7] reach 1 under G  [checked]  time %.1fs" % (NU, time.time() - t1))
print("  max steps-to-1 = %d at m=%d   (F1: 93 at 8751065)" % best_steps)
print("  max peak/m = %.3f at m=%d, peak=%d   (F1: 133.03 at 4847486)" % best_ratio)
check(best_steps == (93, 8751065), "F1 max steps")
check(best_ratio[1] == 4847486, "F1 worst peak/m location")
print("  K_n <= 2n + [k_1=3] verified along every orbit (S5.2)  [checked]")
# full orbit of the worst m
m = 4847486
orb = [m]; word = []
x = m
while x != 1:
    x, k = G(x); orb.append(x); word.append(k)
pk = max(orb)
print("  worst m=%d: word (k_i) = %s" % (m, word))
print("    residues mod 9 along the orbit: %s" % [y % 9 for y in orb[:-1]])
print("    orbit: %s" % orb)
print("    peak %d at step %d, peak/m = %.4f, K at peak = %d, 2*(4/3)^%d = %.1f (bound S5.2 holds)" %
      (pk, orb.index(pk), pk / m, sum(word[:orb.index(pk)]), orb.index(pk), 2 * (4 / 3) ** orb.index(pk)))
def v3(x):
    e = 0
    while x % 3 == 0:
        x //= 3
        e += 1
    return e


def run_of_twos(w, start):
    n = 0
    while start + n < len(w) and w[start + n] == 2:
        n += 1
    return n


# 3-adic reason (truth): m = 5 mod 9, so the word starts with a 3 (factor 8/3); then G(m) = 1 mod 3^6
# exactly, forcing 2^5 0 (S5.1); G^7(m) = 5 mod 9 again, then G^8(m) = 1 mod 3^10 exactly, forcing 2^9 0.
# The peak (step 17) is the end of that second run: pattern (7,5,1^9) of S5.3 stacked on (5,1^5).
check(m % 9 == 5 and word[0] == 3, "worst m is 5 mod 9")
check(v3(orb[1] - 1) == 6 and word[1:7] == [2] * 5 + [0], "G(m) = 1 mod 3^6 exactly -> 2^5 0")
check(orb[7] % 9 == 5 and word[7] == 3, "G^7(m) = 5 mod 9")
check(v3(orb[8] - 1) == 10 and word[8:18] == [2] * 9 + [0], "G^8(m) = 1 mod 3^10 exactly -> 2^9 0")
print("    3-adic reason: m = 5 mod 9 (letter 3), v_3(G(m)-1) = %d (run of %d twos then 0), G^7(m) = 5 mod 9 (letter 3),"
      % (v3(orb[1] - 1), run_of_twos(word, 1)))
print("      v_3(G^8(m)-1) = %d (run of %d twos then 0 at step 17 = the peak): word prefix 3 2^5 0 3 2^9 0, K_17 = %d,"
      % (v3(orb[8] - 1), run_of_twos(word, 8), sum(word[:17])))
print("      exact growth to the peak: 2^%d/3^17 = %.2f  vs observed peak/m %.4f (the carry B_17 costs the difference)"
      % (sum(word[:17]), 2 ** sum(word[:17]) / 3 ** 17, pk / m))

print("""S5.3  Residue-5 chains (PROVED structure).  Residue 5 mod 9 (factor 8/3) can never repeat consecutively
and is always preceded by 7 or 8 (factors 1/3, 2/3): the best 5-containing pattern is 8,5,1^n
with product (2/3)(8/3)(4/3)^n = (16/9)(4/3)^n over n+2 steps, i.e. again rate 4/3.  The word
8,5,1,...,1 (n ones) then 0 is realized by EXACTLY ONE class mod 3^(n+3): m = 8 mod 27 (forces
G(m) = 5 mod 9), G(m) = 5+18s, G^2(m) = 13+48s, and 13+48s = 1 mod 3^n (not mod 3^(n+1)) iff
4s = -1 mod 3^(n-1) with the lift condition.  Smallest members and their exact growth:""")
for nn in range(1, 12):
    # find smallest m = 8 mod 27 with G^2(m) = 1 mod 3^nn but not mod 3^(nn+1)
    found = None
    s = 0
    while found is None:
        m = 8 + 27 * s
        g2 = 13 + 48 * s
        if (g2 - 1) % 3 ** nn == 0 and (g2 - 1) % 3 ** (nn + 1) != 0:
            found = m
        s += 1
    m = found
    x = m; word = []; vals = [m]
    for i in range(nn + 2):
        x, k = G(x); word.append(k); vals.append(x)
    exp_word = [1, 3] + [2] * (nn - 1) + [0]
    check(word == exp_word, "8,5,1^n word for n=%d" % nn)
    print("  n=%2d  m=%10d  word %s  value after n+1 steps = %d  ratio %.4f  (16/9)(4/3)^(n-1)=%.4f" %
          (nn, m, word, vals[nn + 1], vals[nn + 1] / m, (16 / 9) * (4 / 3) ** (nn - 1)))
print("FINITE-EXACT+PROVED (S5.4): the (8,5,1^n) family has word 1,3,2^(n-1),0 exactly and grows at rate 4/3;")
print("  no family beats 2*(4/3)^n (S5.2).  Hence sup_m peak(m)/m = infinity (S5.1), but peak(m) <= 2 (4/3)^{steps} m always.")

# ============================================================================
banner("S6  Typed analogies and the three-piece decomposition")
# ============================================================================
print("""A1  Source: Collatz T on Z_2 (Lagarias 1985, CITED: T is Haar-preserving on Z_2 and the parity-vector map
    Q: Z_2 -> Z_2 is a measure-preserving bijection; Terras 1976, CITED: stopping-time densities -> 1).
    Target: G on Z_3^x.  Map: 2-adic digits <-> 3-adic digits, parity word <-> k-word, n mod 2^J <-> m mod 3^(J+1),
    gate n0(2^K-3^L)=kB(w) <-> m0(2^J-3^L)=kB(rev w).  Preserved: 'word of length J is a function of the
    residue at depth J(+1)', cycle gate shape, Terras-type density theorem.  Lost: Haar-invariance and
    bijectivity of the word map onto the FULL shift (G is 4:1 / 2:1 with a.c. invariant density 4/3, 2/3 != Haar;
    the k-word map IS injective (S1.3ii) but its image is a proper SFT of entropy log 3 < log 4, whereas Lagarias'
    parity map is onto the full 2-shift); the direction of the carry (helps here, hinders there).  Sidecar: the rank-1
    tilted matrix [[5/3,1/3],[10/3,2/3]] with eigenvalue 7/3.  Decisive test: E[k]=1 (3-adic) vs E[parity]=1/2 (2-adic).
A2  Source: squarefree density 6/pi^2 = prod_p (1-p^-2) (user's prompt).  Target: the stationary law (4/3, 2/3).
    Map: 'a global density is a product/limit of local Haar masses' (CRT).  Preserved: density of a residue
    set = Haar measure (S3.1c).  Lost: multiplicativity over primes -- G lives at the single prime 3 (2 acts as
    a unit), so there is no Euler product and NO MAP FOUND from 6/pi^2 to any G-quantity; the analogy stops at
    'local densities govern global counts'.
A3  The user's three pieces {0,1,2} mod 3 in E.  PROVED: 3Z is a source forest (nothing enters it: 3n+1 != 0 mod 3
    and n/2 = 0 mod 3 forces n = 0 mod 3), so '1 generates every number' can only mean every non-multiple of 3
    (wave-one Q2, OPEN, FINITE-EXACT to 10^6), the multiples of 3 being attached by a single last move
    3t -> 9t+1 (E-forward) read backwards.  Inside the units, G's level-1 chain is the 2x2 matrix
    [[2/3,1/3],[2/3,1/3]]: from EITHER class the next class is 1 mod 3 with probability 2/3.  The three
    pieces fit as: 0 mod 3 = leaves; 1 mod 3 = the heavy class (density 4/3); 2 mod 3 = the light class (2/3),
    with the transfer law S1.1 gluing the levels.  This is the exact 3-adic content of F1; it is a statement about
    the GREEDY inverse walk, not about Collatz forward dynamics (whose image residues mod 3 are 1,2 with
    probability 1/2 each: (3n+1)/2^k = (-1)^k mod 3).
A4  3n-5 / 3n+k versions (user): S4 shows the 3-adic greedy inverse of 3n+k has the SAME chain and drift for
    every k with 3 not| k, and its typical cycles are small with sign -sign(k).  No map found from the G_k
    census to Bott periodicity or octonions (no 8-fold or 2-fold structure appears in the data; the period 6 of
    2 mod 9 is the only periodicity present).
A5  Source: the blueprint audit's exact affine word model (collatz_blueprint_20260921_affine.md, Sec. 1 and 3):
    chronological D/E words with W(n) = (2^r n - B)/3^m, legality iff 2^r n = B mod 3^m (one legal source class
    mod 3^m), fixed point n = B/(2^r-3^m), and the ternary/dyadic progression bijection s+3^m j <-> W(s)+2^r j.
    Target: the greedy word of G.  Map: letter k>=1 -> the guarded word E o D^(k-1) (3G+1 = 2^k m); letter k=0
    -> the E-graph's extra arrow (even -> 3n+1 reversed), which is OUTSIDE the D/E semigroup (wave-one leaf
    identity, fibre index j=-1).  Preserved: the carry formula (r = K_J, m = J, B = B_J of S3), the legality
    congruence (the word of length J is legal on exactly one class mod 3^J), the fixed-point/cycle gate, and the
    free-semigroup property (distinct k-words are distinct affine maps).  Lost: the dyadic endpoint address
    mod 2^r and the uniform word mass 2^(-r) -- G selects the word from the 3-adic side (m mod 3^(J+1)), so word
    masses are the Markov products 2*3^J P_chain(w) (wave-one S4), not 2^(-r); and the blueprint's 'source class
    mod 3^m' becomes 'source class mod 3^(J+1)': the greedy rule spends one extra ternary digit to choose k.
    Sidecar: the blueprint's spectral hostile (DEDE vs DDEE, same trace/determinant, different carry) is the
    statement that the cycle gate is NOT a function of (K,L) alone -- exactly S4.3's B'(w) dependence.
    Decisive test: 'same words, opposite adic completion' -- the blueprint counts words by dyadic endpoints
    (mass 2^-r, E-count binomial), G counts them by ternary sources (mass 2*3^J P_chain, E[k]=1).

THEOREM S6.1 (PROVED, the reframe).  Let U = {n>=1 : 3 not| n}.  Q1: every n>=1 reaches 1 in E.  Q2: 1 reaches
every m in U in E.  Then  (Q1 and Q2)  <=>  (U is one strongly connected component of E containing 1, and every
multiple of 3 is a transient singleton feeding it).
Proof.  (=>) For m, m' in U, Q1 gives m -> 1 and Q2 gives 1 -> m', so U lies in one SCC; multiples of 3 are
never entered (3n+1 = 1 mod 3; n/2 = 0 mod 3 forces n = 0 mod 3) and leave 3Z after v_2(n) halvings via
3t -> 9t+1 (wave-one S1.3), so they are singleton SCCs.  (<=) Q2 is the reachability 1 -> m inside the SCC;
Q1 for n in U is m -> 1 inside the SCC, and for n = 3t: 3t -> 9t+1 in U -> 1.  QED
Relation to Collatz (SCOPE, stated exactly): Collatz => Q1 (the deterministic path is an E-path); Q1 => Collatz
is NOT established here (E has extra even -> 3n+1 arrows, so reaching 1 in E is a priori weaker); Q2 is the new
question, with the G-certificate: a G-orbit m -> ... -> 1 reversed is an explicit E-path 1 -> m, so
Q2 <= 'every m in U reaches 1 under G' (FINITE-EXACT to 10^7, S5; stopping-time density 1, S3).""")

# ----------------------------------------------------------------------------
banner("Provenance")
src = open(__file__, "rb").read()
print("source sha256:", hashlib.sha256(src).hexdigest())
print("python:", sys.version.split()[0], " numpy:", np.__version__)
print("total time: %.1fs" % (time.time() - T0))
check(not FAILS, "no failures")
print("ALL CHECKS PASSED")
