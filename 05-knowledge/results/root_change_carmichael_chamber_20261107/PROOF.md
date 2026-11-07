# THM-4904: the Carmichael chamber of root-change words, and two root changes to 10^9

kind-pasteur-2026-11-07-S14. Status: **PROVED** (sections 1-4, given Session199's adjoint and Session200's physical
inequalities, both re-verified here); **FINITE-EXACT** (sections 5-6); **CONJECTURE** (HYP-9494, section 7).
Pre-registration: [PREREGISTRATION.md](PREREGISTRATION.md), pushed in 6d0e88f0e8 before the 10^9 output was read.
Literature and repository map: [CARMICHAEL_MAP.md](CARMICHAEL_MAP.md).

## 0. Objects

Words over {O, Z, J} act on the rooted edge t = (D,E,U,d) = (1,0,2,1) by the physical moves (Session200, section 1)

    O(D,E,U,d) = (D,0,2U,d),  Z = (D-E+U,d,2U-d,d),  J = (D-E+U,D-E,D+U,D).

Every move is a literal vertex birth. The value of a word is the first coordinate D of its final state, which is a
Hamiltonian-path count. Write F(s) = D - E + U. Both Z and J give the new first coordinate F(s), so

    value(wZ) = value(wJ) = F(t.w).

Write #J(w) for the number of letters J in w. For K >= 0 and odd N, the **word multiplicity** is

    M_K(N) = [N = 1] + #{wJ : #J(w) <= K-1, F(t.w) = N} + #{wZ : #J(w) <= K, F(t.w) = N}.

Here the term [N = 1] is the empty word, and words ending in O are omitted (O does not change the first coordinate).
This is exactly the count made by the S13 enumerator `root_changes.c` and by `rc_census.c`. A trailing O inside the
last block is the same as no letter; a trailing O before a J matters and is counted. F_<=K is the set of N with
M_K(N) > 0: the H-values reachable from the rooted edge with at most K root changes (HYP-9480 (D)).

The dagger of a word is w^dagger: reverse w and exchange O with Z, fixing J. Session199's physical adjoint (frozen in
Session200 as `dependencies/PHYSICAL_ADJOINT.md`) proves

    F(t.w) = F(t.w^dagger)   for every word w.                              (adj)

`verify_chamber.c` re-checks (adj) on all 7,174,452 words of length 1..14, with any number of letters J. It also
re-checks Session200's physical inequalities

    D - E >= d >= 1,   D <= U,   U >= 2d                                    (phys)

on all 76,333,010 states with F <= 2^20 and at most two letters J. No exceptions were found in either check.

## 1. Twin decomposition

**Theorem 1.** For K >= 0 and odd N, with A_j(N) = #{w : #J(w) <= j, F(t.w) = N} (A_(-1) = 0) and
B_K(N) = #{w : #J(w) = K, F(t.w) = N},

    M_K(N) = [N = 1] + 2 A_(K-1)(N) + B_K(N).

*Proof.* Group the counted words by the prefix w. A prefix with at most K - 1 letters J contributes both wJ and wZ, of
equal value F(t.w). A prefix with exactly K letters J contributes only wZ. QED.

This is the word-level analogue of phi(x) = phi(2x) for odd x: a trivial partner exists unless the prefix has
exhausted its root-change budget.

## 2. Parity law

**Theorem 2.** Let Fix_K(N) = #{w : w = w^dagger, #J(w) = K, F(t.w) = N}. Then

    M_K(N) = [N = 1] + Fix_K(N)   (mod 2).

*Proof.* The dagger is an involution on words and preserves #J. By (adj) it preserves F(t.w). Hence the set counted
by B_K(N) is a union of dagger orbits of size 2, together with the dagger-fixed words. Theorem 1 finishes. QED.

Dagger-fixed words have the form w = A c A^dagger with c in {empty, J}: a middle letter is fixed only if it is J. Thus
#J(w) = 2 #J(A) + [c = J]. The **K-chamber** is

    K = 0:  A A^dagger,          A in {O,Z}*;
    K = 1:  A J A^dagger,        A in {O,Z}*;
    K = 2:  A A^dagger with one J in A, i.e. B J C C^dagger J B^dagger,   B, C in {O,Z}*.

## 3. Singletons live in the chamber

**Theorem 3.** Let N > 1 with M_K(N) = 1. Then:

1. A_(K-1)(N) = 0, that is, N is not in F_<=K-1 (for K >= 1);
2. the unique word is wZ with #J(w) = K and w = w^dagger;
3. for K >= 1, w neither begins nor ends with J. For K = 2 this means w = B J C C^dagger J B^dagger with B nonempty.

*Proof.* (1) and the first half of (2) follow from Theorem 1. The dagger orbit of the unique prefix has size 1, so the
prefix is dagger-fixed. For (3), Session200's Theorem 1 gives F(O(s)) = F(J(s)) for every source s. If w = uJ,
then F(t.uO) = F(t.uJ) = N. The prefix uO has K - 1 letters J, so A_(K-1)(N) >= 1, contradicting (1). If w begins
with J, then w^dagger = w ends with J, and the same argument applies. QED.

This is the analogue of the Carmichael-Klee chamber (a unique preimage is divisible by 4, then by 2^2 3^2 7^2 43^2,
and so on). The analogue of Ford's remark that a unique preimage is divisible by 6 is (1): a singleton value has all of
its words in the exactly-K layer.

## 4. Singletons have density zero, with a sharp square-root bound

**Lemma 4 (growth).** For every word w,

    U(t.w) >= 2 c0^(#J(w)+1) 2^(n(w)),   with c0 = prod_(i>=0) (1 - (1/4)(2/3)^i) = 0.44257...,

where n(w) is the number of letters O and Z. Consequently F(t.w) > U(t.w), and every word of value at most X has
n(w) <= log2 X + (#J(w)+1) log2(1/c0).

*Proof.* Put rho = d/U. By (phys), rho <= 1/2 for every state.
- O: U' = 2U and rho' = rho/2.
- Z: U' = (2 - rho)U and rho' = rho/(2 - rho) <= (2/3)rho.
- J: U' = D + U > U and rho' = D/(D + U) <= 1/2, because D <= U.

Inside a maximal O/Z block the i-th letter therefore has rho <= (1/2)(2/3)^i, and multiplies U by at least
2(1 - (1/4)(2/3)^i). A block of length L multiplies U by at least c0 2^L. There are #J + 1 blocks, U(t) = 2, and J
never decreases U. Finally F = (D - E) + U >= d + U > U. QED.

`verify_chamber.c` checks Lemma 4 on all 76,333,010 enumerated states with no exceptions.

**Theorem 5.** For each K, the number of N <= X with M_K(N) = 1 is at most the number of K-chamber words of value at
most X, which is O_K(sqrt(X) (log X)^floor(K/2)). In particular the singleton values of every bounded-root-change
grammar have density zero.

*Proof.* By Theorem 3 distinct singleton values have distinct chamber words. A chamber word A c A^dagger is determined
by A, which has n(A) = n(w)/2 letters O or Z and floor(K/2) letters J. By Lemma 4, n(A) <= (log2 X)/2 + C_K. The number
of such A is at most sum over a <= (log2 X)/2 + C_K of C(a + floor(K/2), floor(K/2)) 2^a, which is
O_K(sqrt(X) (log X)^floor(K/2)). QED.

This is where the analogy with Ford's dichotomy breaks. For phi, Ford proved that counterexamples either do not exist
or form a positive proportion of totients. For root-change words, a chamber of square-root size confines singletons
unconditionally. The bound is sharp for K = 0. Below 2^20 the 0-chamber has exactly 1,023 = 2^10 - 1 words, of which
500 are singletons. Below 2^26 there are 3,983 singletons, about 0.49 sqrt(X). The per-window counts grow by a factor
close to sqrt(2) per doubling (section 6).

## 5. Two root changes to 10^9 (pre-registered P1-P4)

The run was `rc_census 1000000000 2 512`: 8,019 tasks, 1,845.7 s of enumeration on 12 threads. Its raw output is
[rc_1e9.txt](rc_1e9.txt), and the outcomes, computed by `evaluate_prereg.py`, are in
[PREREG_OUTCOMES.txt](PREREG_OUTCOMES.txt). The pre-registration was pushed (6d0e88f0e8) while the run was in progress
and before any output line was read.

| Prediction | Outcome |
| --- | --- |
| P1: F_<=2 misses exactly the 140 values below 10^8 (largest 3,780,967) and nothing in (10^8, 10^9] | **CONFIRMED** |
| P2: min M_2 >= 10 in the windows from 2^26 up to 10^9 | **FAILED**: minima 10, 8, 8, 14 |
| P3: no odd 1 < N <= 10^9 has M_2(N) = 1 | **CONFIRMED** |
| P4: the minimum over N = 1 (mod 3) is strictly smallest in each window | **FAILED**: (34,10,10) is a tie; (40,12,8) has residue 2 lowest |

P1 is the decisive test of HYP-9480 (D), extended tenfold. P2 and P4 predicted window extremes of a heavy-tailed
statistic, and both failed. The bulk statistics behind P4 hold (section 5.3), but an extreme of a few hundred
million values is the wrong quantity to register. That lesson goes to MISTAKES. The value of the minimum is not
monotone: above 2^26 it stays between 8 and 14 while the window mean exceeds 200 (the 8-bit counter saturates at 255,
so the printed means are lower bounds).

### 5.1 Root-change depth and the uniquely Hamiltonian spectrum below 10^9 (FINITE-EXACT, exploratory)

Let depth(N) be the least number of letters J in a word of value N. `rc_depth.c` (the census with up to seven
levels and full miss lists) gives [depth_K7_1e6.txt](depth_K7_1e6.txt) and
[depth_K3_4e6.txt](depth_K3_4e6.txt). Combined with P1:

| Number of root changes allowed | odd misses below the cap | largest miss |
| --- | --- | --- |
| <= 2 | 140 (cap 10^9) | 3,780,967 |
| <= 3 | 44 (cap 4,000,000, hence also 10^9) | 27,781 |
| <= 4 | 33 (cap 10^6) | 3,523 |
| <= 5 | 31 | 2,647 |
| <= 6 | 31 | 2,647 |
| <= 7 | 29 | 2,647 |

Depth is therefore:
- 3 for 96 values;
- 4 for the 11 values 373, 1811, 2149, 2999, 4411, 4547, 7129, 17383, 19951, 22091, 27781;
- 5 for 667 and 3523;
- 6 for no value;
- 7 for 1877 and 1897;
- undefined for 29 values.

The 29 are exactly Session197's list of values missing from the C = 1 spectrum: 7, 13, 21, 41, 43, 55, 67, 77, 83,
163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569, 583, 637, 821, 1067, 1339, 1861, 2647.

`literal_witnesses.py` checks this table independently, using exhaustive word searches in which F rises by at least
2 per letter.
- It confirms the minimal depth of every listed depth-4, -5 and -7 value, and of three depth-3 values.
- It builds each witness word, closed by J, as a literal tournament with Session200's builder, and counts its
  Hamiltonian paths and cycles by subset DP. Every one has H = N and exactly one Hamiltonian cycle.
- 1877 = ZJJZJJJJJOO and 1897 = ZZJJJJJJJOO each give a 14-vertex uniquely Hamiltonian tournament.
- The 29 values have no word at all, with any number of root changes.

Session197 (`session197/lowcycle/C1_COMPLETE.md`, section 3) classifies the strong C = 1 tournaments by the rooted
grammar. Session199 (`session199/physical-adjoint/PROOF.md`, Theorem 1) proves that the Hamiltonian-path spectrum of
strong tournaments with exactly one Hamiltonian cycle is exactly {F(s) : s in S}, every value with a literal word
closed by J. So:

**Corollary (FINITE-EXACT).** Every odd N with 2647 < N <= 10^9 is the Hamiltonian-path count of a strong tournament
with exactly one Hamiltonian cycle, built with at most seven root changes. Below 10^9 the strong C = 1 spectrum
misses exactly Session197's 29 values. This extends Session197 (to 5,001) and THM-4884 (2b) (to 50,001).

Session199's balanced readout with at most three letters J (`session199/hidden-structure/GAPS_THREE.json`) has 73
gaps up to 10^6, the largest also 27,781. It is a sub-readout, and its gaps contain the 44 listed here.

### 5.2 Near-misses (exploratory)

Above 4,194,304, every odd N below 10^9 has M_2(N) >= 4. `explore/rc_low.c`, run to 2^28, lists the 130 values in
(2^22, 2^28) with M_2 <= 12 ([explore/rc_low_2p28.txt](explore/rc_low_2p28.txt)).
- By residue mod 3, 104 are 1 (mod 3), 26 are 2 (mod 3), and **none** is divisible by 3.
- The minima are 4, attained at 5,661,961, 8,203,243, 11,340,277, 20,195,929 and 21,375,313.
- No binary pattern (popcount, 2-adic valuation of N - 1 or N + 1) stands out.

### 5.3 A 3-adic local factor (exploratory)

With 16-bit counters (`explore/rc_mod3.c`, `explore/rc_v3.c`, X = 2^24), the mean of M_2 in [2^23, 2^24) depends
strongly on the 3-adic valuation of N.

| v3(N) | 0 | 1 | 2 | 3 | 4 | 5 | >= 6 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| mean M_2 | 205.8 | 332.0 | 457.3 | 510.7 | 556.8 | 614.2 | 618.7 |

The ratios to the 3 ∤ N mean rise from 1 to about 3. Among 3 ∤ N, the class N = 1 (mod 3) is lowest: 199.96 against
211.73 for N = 2 (mod 3). This accounts for the mod-3 profile of the misses and near-misses. Every miss of F_<=2 above
10^5 except 136,607 is 1 (mod 3), and no near-miss is 0 (mod 3). A plausible mechanism, not proved here: ordinal sums
multiply Hamiltonian-path counts, and H(C_3) = 3 = F(t), so multiples of 3 gain extra families of representations. A
derivation of the local factor is OPEN.

## 6. The singleton hierarchy (FINITE-EXACT, exploratory)

Singleton values N > 1 (M_K(N) = 1) per dyadic window [2^k, 2^(k+1)), from `rc_census` at X = 2^26 for K = 0, 1 and
X = 10^9 for K = 2:

- K = 0: 1, 1, 0, 1, 2, 0, 1, 2, 3, 5, 8, 9, 10, 19, 30, 41, 61, 62, 98, 147, 234, 297, 418, 597, 854, 1183
  (k = 0..25). The first values are 3, 11, 25, 31, 91, 179, 213, 353, 371, 459, with words Z, OZZ, ZZOOZ, ZOZOZ, ...
  Every one is dagger-fixed.
- K = 1: 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 5, 2, 0, 2, 4, 3, 2, 7, 1, 2, 3, 4 (k = 0..25; the k = 0 entry is
  N = 1). The singletons above 1 begin 49 (ZZJOOZ), 161 (ZZZJOOOZ), 2457, 12553, 16777, 19785, 24217, 26089, 31377.
- K = 2: none above 1 through 10^9 (P3 CONFIRMED).

Odd multiplicities occur exactly where Fix_K(N) is odd (Theorem 2). For K = 2 there are 2,024 such N up to 2^20 and
10,055 up to 2^24 (`explore/rc_par.c`). That is about 2 to 2.5 sqrt(X), within the O(sqrt(X) log X) chamber bound and
of density zero. The small odd values M_2 = 3, 5, 7 occur only below 2^15, because multiplicities grow.

The K = 0 singleton values include Session200's profile singletons 3 (word Z) and 11 (word OZZ), and the diagonal
two-run value 3971 = 63^2 + 2 (word O^5 Z^6). Session200's diagonal family is exactly the K = 0 chamber of
two-run words: the dagger exchanges O^a Z^b with O^b Z^a.

## 7. Conjectures (HYP-9494)

- (A) K = 0 singletons number (c + o(1)) sqrt(X) with c about 0.49. That is about half of the 0-chamber, matching the
  root's Poisson miss rate e^(-lambda/2) = 0.4985 (HYP-9480 (A)).
- (B) K = 1 singletons are infinite but number O(log X), about 3 per dyadic window.
- (C) K >= 2 grammars have no singleton above N = 1. This is the Carmichael conjecture for words with two root
  changes: FINITE-EXACT through 10^9 (P3).
- (D) Two root changes cover every odd N > 3,780,967 (HYP-9480 (D), now FINITE-EXACT through 10^9). The maximal
  root-change depth over the whole strong C = 1 spectrum is 7, attained only at 1877 and 1897.
- (E) The mean of M_2 has a 3-adic local factor that rises with v3(N) towards about 3 times the 3 ∤ N value
  (exploratory, section 5.3).

The heuristic: a chamber word of value N is alone with probability exp(-lambda_K(N)). Here lambda_K is the mean number
of other words. lambda_0 is bounded, lambda_1 grows like (1/2) ln N (one free block), and lambda_2 like c ln^2 N. The
expected singleton counts per dyadic window are then about sqrt(N), O(1), and exp(-c ln^2 N) sqrt(N), which is summable.
