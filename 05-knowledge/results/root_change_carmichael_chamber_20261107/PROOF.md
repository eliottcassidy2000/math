# THM-4904: the Carmichael chamber of root-change words, and two root changes to 10^9

kind-pasteur-2026-11-07-S14. Status:
- **PROVED**, sections 1-4. New for K >= 1; the K = 0 case and two steps are KNOWN and credited (section 4 end).
- **FINITE-EXACT**, sections 5-6.
- **CONJECTURE**, HYP-9494, section 7.
- **Independent audit:** PASS WITH CORRECTIONS, recorded in [audit/AUDIT.md](audit/AUDIT.md) and applied (MISTAKE-733).

Pre-registration: [PREREGISTRATION.md](PREREGISTRATION.md), pushed in 6d0e88f0e8 before the 10^9 output was written.
Literature and repository map: [CARMICHAEL_MAP.md](CARMICHAEL_MAP.md).

## 0. Objects

Words over {O, Z, J} act on the rooted edge t = (D,E,U,d) = (1,0,2,1) by the physical moves (Session200, section 1)

    O(D,E,U,d) = (D,0,2U,d),  Z = (D-E+U,d,2U-d,d),  J = (D-E+U,D-E,D+U,D).

Every move is a literal vertex birth. The value of a word is the first coordinate D of its final state, which is a
Hamiltonian-path count. Write F(s) = D - E + U. Both Z and J give the new first coordinate F(s), so

    value(wZ) = value(wJ) = F(t.w).

Write #J(w) for the number of letters J in w. For K >= 0 and odd N, the **word multiplicity** is

    M_K(N) = [N = 1] + #{wJ : #J(w) <= K-1, F(t.w) = N} + #{wZ : #J(w) <= K, F(t.w) = N}.

Here the term [N = 1] is the empty word, and words ending in O are omitted, since O does not change the first
coordinate. This is exactly the count made by the S13 enumerator `root_changes.c` and by `rc_census.c`. The auditor's
letter-by-letter enumerator, which uses no closed form, reproduces every window statistic of the census at 2^20, 10^8
and 10^9.

F_<=K is the set of N with M_K(N) > 0: the H-values reachable from the rooted edge with at most K root changes
(HYP-9480 (D)). The root-change depth of N is the least K with N in F_<=K.

The dagger of a word is w^dagger: reverse w and exchange O with Z, fixing J. Session199's physical adjoint
(`session199/physical-adjoint/PROOF.md`, section 2, eq. (4)-(5)) proves

    F(t.w) = F(t.w^dagger)   for every word w.                              (adj)

The inequalities

    D - E >= d >= 1,   D <= U,   U >= 2d                                    (phys)

hold for every state (Session200, section 1). The auditor's Appendix A gives a five-line induction: the set
{E >= 0, E <= d, D - E >= d >= 1, D - E + d <= U, U >= 2d} contains t and is closed under O, Z and J.

Finite re-checks, with no exceptions:
- `verify_chamber.c` checks (adj) on all 7,174,452 words of length 1 to 14, and (phys) on all 76,333,010 states with
  F <= 2^20 and at most two letters J;
- the auditor checks (adj) on 64,570,080 words of length up to 16.

## 1. Twin decomposition

**Theorem 1.** For K >= 0 and odd N, put A_j(N) = #{w : #J(w) <= j, F(t.w) = N}, with A_(-1) = 0, and
B_K(N) = #{w : #J(w) = K, F(t.w) = N}. Then

    M_K(N) = [N = 1] + 2 A_(K-1)(N) + B_K(N).

*Proof.* Group the counted words by the prefix w. A prefix with at most K - 1 letters J contributes both wJ and wZ, of
equal value F(t.w). A prefix with exactly K letters J contributes only wZ. QED.

This is the word-level analogue of phi(x) = phi(2x) for odd x: a trivial partner exists unless the prefix has
exhausted its root-change budget.

## 2. Parity law

**Theorem 2.** Let Fix_K(N) = #{w : w = w^dagger, #J(w) = K, F(t.w) = N}. Then

    M_K(N) = [N = 1] + Fix_K(N)   (mod 2).

*Proof.* The dagger is an involution on prefixes and preserves #J. By (adj) it preserves F(t.w); the counted word is
always wZ, and the dagger acts on the prefix w. Hence the set counted by B_K(N) is a union of dagger orbits of size
2, together with the dagger-fixed words. Theorem 1 finishes. QED.

Dagger-fixed words have the form w = A c A^dagger with c in {empty, J}, since a middle letter is fixed only if it is J.
Thus #J(w) = 2 #J(A) + [c = J]. The **K-chamber** is:

    K = 0:  A A^dagger,          A in {O,Z}*;
    K = 1:  A J A^dagger,        A in {O,Z}*;
    K = 2:  A A^dagger with one J in A, i.e. B J C C^dagger J B^dagger,   B, C in {O,Z}*.

**K = 0 is KNOWN.** On J-free prefixes, the dagger is Session206's converse involution
T(2^(n-1) + a) = 2^n - 1 - rev_(n-1)(a) (`session206/root-oz-noncofiniteness/ROOT_OBSTRUCTION.md`, sections 1-2;
THM-4888 (B)). It maps the tail word b_0 ... b_(n-2) Z to dagger(b_0 ... b_(n-2)) Z. Its fixed points are exactly the
0-chamber A A^dagger. So the K = 0 case of Theorem 2 is the orbit decomposition of T. The extension to K >= 1 through
Session199's adjoint is new.

## 3. Singletons live in the chamber

**Theorem 3.** Let N > 1 with M_K(N) = 1. Then:

1. A_(K-1)(N) = 0. That is, N is not in F_<=K-1, and N has root-change depth exactly K.
2. The unique word is wZ with #J(w) = K and w = w^dagger.
3. For K >= 1, w neither begins nor ends with J. For K = 2 this means w = B J C C^dagger J B^dagger with B nonempty.
4. w contains neither of the factors JZ and OJ. Consequently:
   - for K = 1, A ends in Z;
   - for K = 2, B ends in Z, and C is empty or begins with O.

*Proof.*
- (1) and the first half of (2) follow from Theorem 1.
- The dagger orbit of the unique prefix has size 1, so the prefix is dagger-fixed.
- For (3), Session200's Theorem 1 gives F(O(s)) = F(J(s)) for every source s. If w = uJ, then F(t.uO) = F(t.uJ) = N.
  The prefix uO has K - 1 letters J, so A_(K-1)(N) >= 1, contradicting (1). If w begins with J, then w^dagger = w ends
  with J, and the same argument applies.
- For (4), Session205's operator identity JZ = OJ (THM-4888 (D)1) says that uJZv and uOJv reach the same state with the
  same number of letters J. So a prefix containing either factor has a second prefix in the same B_K-layer with the
  same value.

QED.

**Credit.** The J-exclusion in (3) is the argument of Session200 carmichael-fibers, Theorem 3 (a terminal J or O gives
a twin value, and the dagger moves the start to the end). Here it is run inside a root-change budget. The JZ/OJ
remark (4) was suggested by the audit (m2). It is checked on all 520 singletons below 2^20 with K <= 2, and on all
4,122 below 2^26 with K <= 1.

This is the analogue of the Carmichael-Klee chamber: a unique preimage is divisible by 4, then by 2^2 3^2 7^2 43^2,
and so on. Ford remarks that a unique preimage is divisible by 6; the analogue of that remark is (1), which says a
singleton value has all of its words in the exactly-K layer.

## 4. Singletons have density zero, with a square-root chamber bound

**Lemma 4 (KNOWN, Session200 `root-weighted-height/PROOF.md`, section 1).** Put W = U - d. The moves act by

    O: W' = 2W + d,   Z: W' = 2W,   J: W' = W + d,

and W(t) = 1. Hence W(t.w) >= 2^(n(w)) for every word, where n(w) counts the letters O and Z. So

    F(t.w) = (D - E) + U >= d + U = W + 2d >= 2^(n(w)) + 2.

The first version of this packet proved a weaker bound, U >= 2 c0^(#J+1) 2^(n(w)) with c0 = 0.44257, by the ratio
rho = d/U. `verify_chamber.c` still checks that version on 76.3 M states. The auditor found the KNOWN W-bound, which
loses nothing per root change and attains equality at Z^n.

**Theorem 5.** For each K,

    #{1 < N <= X : M_K(N) = 1} <= #{chamber words w with F(t.w) <= X} = O_K(sqrt(X) (log X)^floor(K/2)).

In particular the singleton values of every bounded-root-change grammar have density zero. For K = 0 the chamber
itself has Theta(sqrt(X)) words.

*Proof.* By Theorem 3, distinct singleton values have distinct chamber words. A chamber word A c A^dagger is
determined by A, which has n(A) = n(w)/2 letters O or Z and floor(K/2) letters J. By Lemma 4, n(A) < (1/2) log2 X.
The number of such A is at most 2 sqrt(X) C((1/2) log2 X + floor(K/2), floor(K/2)).

For the lower bound at K = 0: U at most doubles per letter, and F <= 2U. So every A in {O,Z}^a with 2a + 2 <= log2 X
gives a chamber word of value at most X, which is at least sqrt(X)/2 - 1 words. QED.

The K = 0 case of Theorem 5 is again Session206's orbit count for T. New here are K >= 1 and the comparison with
Ford.

**The contrast with phi.** Ford proved that counterexamples to Carmichael's conjecture are either absent or a
positive proportion of totients. Here a square-root chamber confines the singletons unconditionally. Counting
chamber words cannot improve the bound for K = 0. Whether the singletons are themselves of order sqrt(X) is
HYP-9494 (A):
- below 2^20 the 0-chamber has exactly 1,023 = 2^10 - 1 words, and 500 of them are singletons;
- below 2^26 it has 8,191 words, and 4,083 are singletons, which is 0.498 sqrt(X).

## 5. Two root changes to 10^9 (pre-registered P1-P4)

The run was `rc_census 1000000000 2 512`: 8,019 tasks, 1,845.7 s of enumeration on 12 threads. Its raw output is
[rc_1e9.txt](rc_1e9.txt), and the outcomes, computed by `evaluate_prereg.py`, are in
[PREREG_OUTCOMES.txt](PREREG_OUTCOMES.txt). The timing evidence is in the README:
- the run started at about 10:55:34;
- the pre-registration was pushed at 11:00:57;
- the output was written at 11:26:28;
- stdout is printed only after the enumeration ends.

The auditor reran the whole census with an independent segmented enumerator in 875 s and reproduced every outcome.

| Prediction | Outcome |
| --- | --- |
| P1: F_<=2 misses exactly the 140 values below 10^8 (largest 3,780,967) and nothing in (10^8, 10^9] | **CONFIRMED** |
| P2: min M_2 >= 10 in the windows from 2^26 up to 10^9 | **FAILED**: minima 10, 8, 8, 14 |
| P3: no odd 1 < N <= 10^9 has M_2(N) = 1 | **CONFIRMED** |
| P4: the minimum over N = 1 (mod 3) is strictly smallest in each window | **FAILED**: (34,10,10) is a tie; (40,12,8) has residue 2 lowest |

P1 is the decisive test of HYP-9480 (D), extended tenfold. P2 and P4 registered window extremes of a heavy-tailed
statistic, and both failed. The bulk statistics behind P4 hold (section 5.3). Registering an extreme over hundreds of
millions of values was the mistake, recorded in MISTAKE-733.

The minimum is not monotone. Above 2^26 it stays between 8 and 14, while the window mean exceeds 200. The 8-bit
counter saturates at 255, so the printed means are lower bounds. The auditor finds only two values with M_2 <= 12 in
[2^28, 10^9]: 277,438,727 (8 words) and 301,189,927 (12 words).

### 5.1 Root-change depth and the uniquely Hamiltonian spectrum below 10^9 (FINITE-EXACT, exploratory)

`rc_depth.c` is the census with up to seven levels and full miss lists. It gives [depth_K7_1e6.txt](depth_K7_1e6.txt)
and [depth_K3_4e6.txt](depth_K3_4e6.txt). Combined with P1:

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

Root changes are counted in source words. The closing J that turns a source into its C = 1 graph is not counted, which
is Session199's convention, so a depth-7 witness has eight letters J in all.

`literal_witnesses.py` checks this table independently, using exhaustive word searches in which F rises by at least
2 per letter.
- It confirms the minimal depth of every listed depth-4, -5 and -7 value, and of three depth-3 values.
- It builds each witness word, closed by J, as a literal tournament with Session200's builder, and counts its
  Hamiltonian paths and cycles by subset DP. Every one has H = N and exactly one Hamiltonian cycle.
- 1877 = ZJJZJJJJJOO and 1897 = ZZJJJJJJJOO each give a 14-vertex uniquely Hamiltonian tournament.
- The 29 values have no word at all, with any number of root changes.
- The auditor reproduces all of this with its own code, and its builder agrees with the state machine on all 29,524
  words of length at most 9.

The readout "strong C = 1 spectrum = F(S)" is Session198 c1-monotone, section 1 (`session198/c1-monotone/PROOF.md`).
Session199 (`physical-adjoint/PROOF.md`, Theorem 1) proves that the C = 1 spectrum is exactly {F(s) : s in S}, every
value with a literal word closed by J. Session197 (`lowcycle/C1_COMPLETE.md`, section 3) classifies the strong C = 1
tournaments by the rooted grammar. So:

**Corollary (FINITE-EXACT).** Every odd N with 2647 < N <= 10^9 is the Hamiltonian-path count of a strong tournament
with exactly one Hamiltonian cycle, from a source word with at most seven root changes. Below 10^9 the strong C = 1
spectrum misses exactly Session197's 29 values. This extends Session197 (to 5,001) and Session198 c1-monotone,
section 3 (to 50,001, replayed in THM-4884 (2b)).

Prior depth data: Session199 `hidden-structure/PROOF.md`, section 6, tabulated the minimal J cost (root-change depth)
of the balanced sub-readout with at most three J's, to 1,000,001. Its cost 0/1/2/3 counts are 213,873 / 266,697 /
19,137 / 220. It has 73 gaps, and 293 misses with at most two J's, the largest 879,001. The table here is for the full
C = 1 readout, to 10^9.

### 5.2 Near-misses (exploratory)

Above 4,194,304, every odd N below 10^9 has M_2(N) >= 4. `explore/rc_low.c`, run to 2^28, lists the 130 values in
(2^22, 2^28) with M_2 <= 12 ([explore/rc_low_2p28.txt](explore/rc_low_2p28.txt)).
- By residue mod 3, 104 are 1 (mod 3), 26 are 2 (mod 3), and **none** is divisible by 3.
- The minima are 4, attained at 5,661,961, 8,203,243, 11,340,277, 20,195,929 and 21,375,313.
- No binary pattern (popcount, 2-adic valuation of N - 1 or N + 1) stands out.

### 5.3 A 3-adic local factor (exploratory)

With 16-bit counters (`explore/rc_mod3.c`, `explore/rc_v3.c`, X = 2^24; exact, since the maximum M_2 below 2^24 is
3,792), the mean of M_2 in [2^23, 2^24) depends strongly on the 3-adic valuation of N.

| v3(N) | 0 | 1 | 2 | 3 | 4 | 5 | >= 6 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| mean M_2 | 205.8 | 332.0 | 457.3 | 510.7 | 556.8 | 614.2 | 618.7 |

The ratios to the 3 ∤ N mean rise from 1 to about 3. Among 3 ∤ N, the class N = 1 (mod 3) is lowest: 199.96 against
211.73 for N = 2 (mod 3).

This is consistent with the mod-3 profile of the misses and near-misses. Every miss of F_<=2 above 10^5 except
136,607 is 1 (mod 3), and no near-miss is 0 (mod 3). A 6% gap in the class means does not by itself explain a 104:26
split.

One mechanism, an illustration rather than a derivation: Session206 proves that a seed with 3 | gcd(R, d) confines
each parity of its tail to a single class mod 3. The two-root-change seed JJ(t) = (5,2,6,3) is the dilated root
H = 3 H_t(m) + 2 - 2(m mod 2) (THM-4888 (D)3), whose values are 0 or 2 (mod 3) and never 1. A derivation of the local
factor is OPEN.

## 6. The singleton hierarchy (FINITE-EXACT, exploratory)

Singleton values N > 1 (M_K(N) = 1), per dyadic window [2^k, 2^(k+1)):

- **K = 0**, from `rc_census` at X = 2^26 (`explore/rc_K0_2p26_windows.txt`), k = 0..25:

      1, 1, 0, 1, 2, 0, 1, 2, 3, 5, 8, 9, 10, 19, 30, 41, 61, 62, 98, 147, 234, 297, 418, 597, 854, 1183

  The k = 0 entry is N = 1, so there are 4,083 above 1. The first values are 3, 11, 25, 31, 91, 179, 213, 353, 371,
  459. Their counted words are Z, OZZ, ZZOOZ, ZOZOZ, ...; every prefix is dagger-fixed.
- **K = 1**, k = 0..25: 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 5, 2, 0, 2, 4, 3, 2, 7, 1, 2, 3, 4. Again the k = 0
  entry is N = 1. The singletons above 1 begin 49 (counted word ZZJOOZ), 161 (ZZZJOOOZ), 2457, 12553, 16777, 19785,
  24217, 26089, 31377.
- **K = 2**: none above 1 through 10^9 (P3 CONFIRMED).
- **K >= 3**: a singleton must sit at a value of root-change depth exactly K, with a single, dagger-fixed exactly-K
  prefix (Theorem 3 (1)). Below 10^9 this happens once: **M_7(1897) = 1**. The word is ZZJJJJJJJOOZ, the prefix
  A J A^dagger with A = ZZJJJ. The audit found it in this packet's own `depth_K7_1e6.txt` (window [1024, 2048),
  c1 = 1). By contrast, M_7(1877) = 4. No other word-singleton above 1 with 2 <= K exists below 10^9.

Odd multiplicities occur exactly where Fix_K(N) is odd (Theorem 2). For K = 2 there are 2,024 such N up to 2^20 and
10,055 up to 2^24 (`explore/rc_par.c`). That is about 2 to 2.5 sqrt(X), within the O(sqrt(X) log X) chamber bound and
of density zero. The small odd values M_2 = 3, 5, 7 occur only below 2^15, because multiplicities grow.

**The two-run K = 0 chamber and repunits (exploratory, `explore/diag_root.c`).** The two-run dagger-fixed J-free words
are O^a Z^a and Z^a O^a.
- O^a Z^a is Session200's Mersenne diagonal, with value 2 + (2^(a+1) - 1)^2: two plus the square of a binary
  repunit. It begins with O.
- Z^a O^a is the other half, with value (2^a + 1)^2.
- For every a <= 13 both families have odd M_0, as Theorem 2 allows.
- The Mersenne diagonal is a K = 0 singleton for a = 1, 5, 9, 11 (11, 3971, 1046531, 16769027).
- The Fermat-square half is a K = 0 singleton for a = 2, 5, 9, 13 (25, 1089, 263169, 67125249).
- With root changes allowed, Session201 shows that the even-a Mersenne diagonal values have at least four profiles.

## 7. Conjectures (HYP-9494)

- **(A)** K = 0 singletons number (c + o(1)) sqrt(X) with c about 0.50: 4,083 of the 8,191 chamber words below 2^26.
  That is half of the 0-chamber, matching the root's Poisson miss rate e^(-lambda/2) = 0.4985 (HYP-9480 (A)).
- **(B)** K = 1 singletons are infinite but number O(log X), about 3 per dyadic window.
- **(C)** For **K = 2**, no value N > 1 has exactly one word. This is the Carmichael conjecture for words with two root
  changes, FINITE-EXACT through 10^9 (P3).
  - (C) is stated only for K = 2. For larger K, Theorem 3 (1) puts singletons at values of exact depth K, where
    M_K = B_K has no doubled term, and M_7(1897) = 1 shows they occur.
  - The first draft stated (C) for all K >= 2; the audit refuted that (M4).
- **(D)** Two root changes cover every odd N > 3,780,967 (HYP-9480 (D), now FINITE-EXACT through 10^9). The maximal
  root-change depth over the whole strong C = 1 spectrum is 7, attained only at 1877 and 1897.
- **(E)** The mean of M_2 has a 3-adic local factor that rises with v3(N) towards about 3 times the 3 ∤ N value
  (exploratory, section 5.3).

The heuristic: a chamber word of value N is alone with probability exp(-lambda_K(N)), where lambda_K is the mean
number of other words.
- lambda_0 is bounded;
- lambda_1 grows like (1/2) ln N, for one free block;
- lambda_2 grows like c ln^2 N.

The expected singleton counts per dyadic window are then about sqrt(N), O(1), and exp(-c ln^2 N) sqrt(N), which is
summable.
