# Independent audit of THM-5126 / HYP-9658 (packet `kps_s21_two_adic_lines_20261110`)

**Audited state:** commit `7035dd3c6e` (checkpoint 2: cut-matrix form (1'), all 65 lines `c <= 129`, order-ten run,
rewritten HYP-9658). The audit started on checkpoint 1 (`b576542d6c`) and was redone on the current files.
**Auditor:** independent agent, fresh context, own code (`aud_engine.c`, `aud_checks.py`, `aud_compare.py`; nothing
imported from the packet's scripts; only the packet's data files are read, for comparison). Total compute about
4 minutes. No file outside `audit/` was edited; nothing was committed.

## Verdict: PASS WITH CORRECTIONS

- **The mathematics is sound.** Every proof in THM-5126 (A)-(E) was checked line by line and no gap was found that
  affects a statement. The headline - *every odd `N > 1` whose `N - 1` has odd part at most 129 is a
  Hamiltonian-path count, except 7 and 21* - is **confirmed in full by independent recomputation** (all 65 seeds, all
  121 points below the rays), conditional only on the tail law (proved, and tested from the definitions), the
  substitution law (THM-3121, tested literally) and THM-1370.
- **Every census number was reproduced** with an independent class generator and independent counters, except three
  small slips (minor findings 1-3).
- **The corrections are about credit, scope and wording:** two identities typed "new" are already in the repository
  (M1, M2), one citation points to the wrong theorem (M3), the negative statements about "no seed of order <= 10"
  exceed the searched universe (M4; the audit closes most of the gap), and three README sentences are false or
  overstated as written (M5).

The Chan-Pak conjecture remains OPEN; nothing in the packet or in this audit bounds the spectrum in density.

---

## MAJOR findings (false, unproved-but-claimed, or miscredited)

### M1. The cut-matrix form (1') is KNOWN (Session216), not new

Quoted (THM-5126 status): *"new here are the identities c0 = d - delta' = H(B) - W(O) and M = H - E + delta' =
U_B(O), which express the offset through the cut matrix W of THM-4114 / THM-5122."* Also the heading *"Offset formula
(new form)"* and README section 1: *"The law is known; the identities are the new form."*

What is already in the repository, with the dictionary `h = H(B)`, `A = U_B(O)`, `T = W(O) = delta'`, `E0 = d - E`:

- `session216_overnight_constructions_20261107/three-core/PROOF.md` (2026-11-07), section 2: the cut coordinates
  `(h, q, A, G, T, S0, E)` - *"A count ordered spanning path covers (P,Q) with Q nonempty starting in S and P
  optional; ... T count such covers with both paths nonempty, P ending in S and Q starting in S"* - and eq. (3)
  `H = 2 A_L + h_L - T_L` for the closure `L -> v -> x` (`x` beats exactly the cut), eq. (4)
  `F(p, r, b) = (qp + A)(2^(b+1) - rev_b(r)) + h(2r + 1 - f(r)) - pG - h f(p) - T - ...`. At `p = r = 0` this is
  `h - T + A 2^(b+1)`, i.e. exactly (1) + (1') with `t = b + 1`.
- `session216.../star/PROOF.md` section 6, eq. (A): appending an unselected sink gives `alpha' = 2 alpha`,
  `tau' = tau` - the statement that the tail doubles the slope and fixes the offset deficit.
- Restated in `session229_combined_constructions_20261109/dynamic/PROOF.md` (1), THM-4960 (*"The path count is
  2A+h-T"*), `session244.../audit/REPORT.md` (*"For a root attached to a cut, h+A-T-E0 ... The x,y closing gadget
  therefore has count h+2A-T"*), `session245.../interfaces/PROOF.md` section 1 (*"(h+2A-T, h, h+2A, h) ... the extra
  sink doubles A and leaves T unchanged"*: inserting this rooted state in (1) gives `c0 = h - T`, `M = 2A`
  immediately), and THM-5122 / Session246 `switches/PROOF.md` (1): `F(x) = h + 2 a^T x - (1/2) x^T W x`, the case
  `t = 1` in the very notation of (1').

So `delta'` ("crossed covers") is the coordinate `T` of the seven-coordinate interface, and (1') is the star-closure
law of Session216 written with `W` and `U_B`. THM-5126 cites THM-5122 for the matrix `W` but not for the identity,
and does not cite Sessions 216/229/244/245 or THM-4960 at all.

What is new in (A), as far as the audit can tell: the translation to rooted O/Z coordinates
(`H - E + 2d - U = h - T`, `U - d = A`), the reading of `h - T` as the 2-adic offset of the fixed-root ray, the
remarks "`c0` odd" and "`c0 <= H(B)` with equality iff `|O| <= 1`", and the use of `T = h - 1` as a search equation
(unit sets).

Proposed status text: *"(A) KNOWN: the tail law (...) and its cut form (1') - the star-closure law of Session216
(three-core/PROOF.md (3)-(4), star/PROOF.md section 6: `h - T + 2^t A`, with `A = U_B(O)`, `T = W(O)`), restated in
Session229 dynamic (1), THM-4960, Sessions 244-245 and THM-5122 (1). Re-derived here in rooted coordinates
(`delta' = T`). New here: the reading of `h - T` as the 2-adic offset of the ray, its parity and its equality case,
and the unit-offset equation `W(O) = H(B) - 1` as a seed criterion."* Heading: *"Offset formula (KNOWN in cut
coordinates; rooted form)"*. README section 1: *"The law and both identities are known (Session216 star closure,
`h - T + 2^t A`); what is used here is the reading of `h - T` as the offset."*

### M2. The extension-mean law at `k = 1` is KNOWN (THM-4111, THM-4115)

Quoted (status): *"(C) PROVED, elementary (new as an identity in this form; the one-vertex cut form is THM-4114 (32)
and its strict positivity is THM-5122)."*

THM-4111 (2026-08-25) Theorem 1.1, eq. (4): `sum over S of H(T + x_S) = 2^(n-2) ((n+3) H(T) + F_1(T))`, and eq. (6):
`sum_(a != b) Q_T(a,b) = (n-1) H(T) + F_1(T)`. Since `P_2(T) = sum_(a != b) Q_T(a,b)`, this is literally
`2^(n-2) (4H + P_2)`. THM-4115 eq. (5) writes the mean as `H + W/2` with `W = P_2/2`, and eq. (6) gives the exact
variance of the same family. Neither theorem is cited.

Proposed: *"(C) PROVED, elementary. KNOWN at k = 1: THM-4111 (4), (6) and THM-4115 (5) (`P_2 = (n-1)H + F_1`; the
variance is THM-4115 (6)). New here, as far as a canon search shows: the identity for every k and the bound
`P_j >= j! S(n, j)` with its equality case."* The general-`k` identity, the bound and the uniqueness of the minimiser
are correct (checked below).

### M3. The module law is cited to the wrong theorem

Quoted: *"(B) PROVED from (A) and the KNOWN module law (THM-2195)"*, the `related:` entry *"THM-2195-... (module
law)"*, and (B1) *"This is the KNOWN module law"*.

THM-2195 (*"Transitive quotients exactly control universal tournament substitution products"*) is about
arc-reversal distance between substitutions; it contains no statement about Hamiltonian paths. The law
`H(K[M at x]) = sum_j c_j(K, x) P_j(M)` is THM-3121, eq. (4)
(`THM-3121-path-cover-walk-content-substitution-kernel.md`: `H(Q[S_1..S_q]) = sum_c [x^c] W_Q(x) prod_i c_i!
pc_(S_i)(c_i)`), with THM-1975 for the path-cover polynomial. Proposed: replace "THM-2195" by "THM-3121 (4)" in the
status, in `related:` and in (B1).

### M4. Negative statements about seeds of order ten exceed the searched universe (mostly closed by the audit)

Quoted (THM status): *"2,016 of the 2,049 odd c <= 4097 have a seed of order <= 10; the least c without one is
245"*. README section 2: *"2,016 of the 2,049 odd c <= 4097 have a seed of order <= 10; the 33 without are 245,
853, ..."*. HYP-9658 title and (B): *"the 33 odd c <= 4097 without a seed of order <= 10 (least 245)"*.

A seed is a unit-offset root or a doubling vertex. The packet searched unit sets of size 2-4 over bases of order
nine (PREREGISTRATION and the HYP status say so: *"out-degree >= 5 is NOT searched"*) and doubling vertices only to
order eight. The three sentences above carry neither restriction, so "without a seed of order <= 10" is not
established by the packet. The scope note is present in PREREGISTRATION.md, in the HYP status and outcome lines and in
the README table of section 4; it is missing from the THM status, README section 2 and the HYP title / (B).

What the audit adds (own census, `aud_engine census`, `dbl9`):

- **Every** vertex set `O` (all sizes `0..9`) over **every** class of order `<= 9` was tested: there is no unit set of
  size `>= 4` over any base of order `<= 9`. So the list of unit-offset roots of order ten is complete (19,002 unit
  sets of size 2 and 3, plus the fans over the transitive base) and the caveat "out-degree >= 5 not searched" can be
  replaced by "none exists".
- Doubling vertices of order nine (all 1,723,824 class-vertex pairs): 972 with `c_2 > 0`, 160 values, 97 odd `c`.
  No `c <= 4097` gains a seed, so the 33 exceptions stand. **But the "least s0" drops on nine lines**:
  `c = 27: 6 -> 3`, `45: 3 -> 1`, `63: 4 -> 1`, `85: 3 -> 1`, `87: 2 -> 1`, `105: 3 -> 1`, `165: 6 -> 1`,
  `219: 3 -> 1`, `261: 5 -> 1` (each recounted literally, `aud_dbl9_extension.out`). The table in (B3) is correctly
  labelled "least s0 found" and (B3) is unaffected, but HYP-9658's *"The least s0(c) found is at most 8"* and the
  table should say that module seeds were searched to order eight only.
- Doubling vertices of order ten were **not** searched by anyone.

Proposed: *"2,016 of the 2,049 odd c <= 4097 have a seed among the unit-offset roots of order <= 10 (complete: the
audit found no unit set of size >= 4 over a base of order <= 9) and the doubling vertices of order <= 9 (audit);
doubling vertices of order ten were not searched. The least c without such a seed is 245."*

### M5. Three README sentences are false or overstated as written

(a) Section 7: *"The mean prices a family. ... compute the exact mean of the count; if it is at least `2^m / C` the
family cannot cover an interval by itself. (THM-5126 (C).)"* This inference is false. A family of `2^m` members
whose values are exactly the odd numbers `1, 3, ..., 2^(m+1) - 1` has mean `2^m` (`C = 1`) and covers an interval.
A lower bound on the mean bounds the number of members per unit of mean; it does not forbid covering. Proposed:
*"... if it is at least `2^m / C`, the family has at most `C` members per unit of its mean, so no coverage statement
can come from counting members; whether it covers is a collision (second-moment) question."*

(b) Section 3: *"This is the exact form of the observation that one-apex counts cover a positive proportion but not
everything."* The mean law implies neither half. That one-apex counts have positive density is a numerical
observation, not a theorem (the README itself says in section 5 that positive density of the spectrum is unproved),
and "not everything" is the fact that 13 is missed, which the mean does not give. Proposed: *"This is consistent
with the numerical observation that one-apex counts have density about one half and miss 13, 29, 35, ...; the mean
alone implies neither statement."* The theorem file was already softened ("has the least mean"); the README still
says "densest" and "exact form".

(c) Section 6.1: *"Unit offset means the ray points at the 2-adic integer 1, which is a path count only as the empty
tournament."* False: 1 is the path count of every transitive tournament, of every order. The typing "EXACT,
inherited" covers only the first sentence of that item (THM-5096); the gloss about a formal inverse
`Z^(-infinity)` is an analogy. Proposed: *"... the 2-adic integer 1, the path count of the transitive tournaments;
no tournament with a tail of negative length is meant."*

---

## MINOR findings

1. **HYP-9658:** *"8 of the 336 are trivial (M = 0)"*. Own count: **7** of the 336 (the transitive tournament rooted
   at its sink, one per order 2..8); 8 of the 1,669 through order nine. The packet's own root file gives 7.
2. **README section 5:** *"(eight of them for c <= 257 are above the order-nine range)"*. For `129 < c <= 257` there
   are 108 points below the rays and **six** are not values of order `<= 9`: 5217, 5281, 5345, 10689, 4113 (above
   3,357) and **3153 = 1 + 197 * 16**, which is inside the range but is not attained at order nine. "Eight" is right
   only if the two points 5313, 3425 of `c <= 129` are added. Four of the six have order-join certificates
   (`5217 = 3 * 1739`, `5345 = 5 * 1069`, `3153 = 3 * 1051`, `4113 = 3 * 1371`); 5281 lies on the order-nine module
   ray of `c = 165` (M4); 10689 has neither. `c = 245` has no seed.
3. **README section 4:** *"49 of the 320 values have a doubling realiser"* counts the trivial value 1 (every vertex
   of a transitive tournament). 48 values have a doubling realiser with `c_2 > 0`.
4. **Citation "THM-4888 (B)"** (status, README section 0) for the tail law: in THM-4888 the tail formula is in the
   Setting (Session 205); (B) is the involution rigidity. The ray (1) itself is THM-5085 (3), "bridge rays", at
   `x = 2`.
5. **(A), last paragraph:** *"the only unit-offset one-apex roots are the fans (A(u) a power of two)"*. Also
   `A(u) = 0` (the transitive tournament rooted at its sink, `M = 0`).
6. **(A), proof of (1'):** *"W_ef is even and positive (THM-5122)"*. THM-5122 states positivity
   (`W_ij >= 2 H(C - i - j)`); evenness is not stated there. One line suffices for the parity:
   `c0 = H(L o 0) - 2M` is odd by Redei, so `W(O) = H(B) - c0` is even for every `O`, in particular for `|O| = 2`.
   An insertion proof of `W_ef >= 2 H(B - e - f)` is in the appendix below (the audit also checked the inequality on
   every pair of every class of order `<= 9`).
7. **(B1), definition of a doubling vertex versus the test in `dbl.c`** (`a = 2, 3, 4`). The new sentence *"A walk
   visits x at most min(d+(x), d-(x)) + 1 times"* makes the test complete at order `<= 8`. In fact `c_3 = 0` already
   forces `c_j = 0` for all `j >= 4` at every order (appendix), so `a = 2, 3` suffice. The audit computed
   `c_2, ..., c_5` by walk counting and found no disagreement with the `h_2, h_3, h_4` test on any vertex of any class
   of order 3-8.
8. **(E):** *"exact when no walk visits x or y three times"* should read "three or more times". (At order `<= 8`
   exactness on the grid `a, b <= 3`, the test of `modfam.c`, coincides with true exactness: 20,702 triples either
   way, `aud_exact4.out`. So the tuples used by `cover3.py` up to `a, b = 26` are legitimate.)
9. **README section 2:** *"The hostile fact '11 cannot be doubled' is the same statement."* It is a consequence of
   `21 not in S`, not an equivalent statement. Proposed: "... is a consequence."
10. **README section 5:** *"the first statement of this thread that is complete on an infinite set defined by a
    simple arithmetic condition"*. The line `c = 1` (`N - 1` a power of two, `H(R_n) = 2^(n-2) + 1`, THM-1370) is
    such a statement and is old. Proposed: "... extends this from the classical line `c = 1` to the 65 lines
    `c <= 129`." For calibration: all odd `N <= 2^32` other than 7, 21 are already path counts by THM-4875 (D)
    (FINITE-EXACT, audited; its addendum extends this to `2^35` in one run), so the new content of the headline is
    the part of each line above that range, and every point below a ray (largest 5313) was already covered there by
    a different construction.
11. **README section 6.1:** see M5 (c).
12. **README section 6.2:** the identity is true for all `n` and has a three-line proof (appendix); the proposed
    route through Grinberg-Stanley is not needed. It can be typed PROVED. It is linearity of expectation
    (`E[P_j] = n! C(n-1, j-1) 2^(j-n)`); the "odd order" index set is the expansion of `((2+t)/(2-t))^x`.
13. **README section 4:** *"15 has no doubling realiser of order <= 8"* can be stated for every order (appendix).
14. **Counting conventions.** "19,002 unit sets" counts pairs (base class, `O`); "1,333" counts class-vertex pairs.
    They are not the same count (order nine: 1,328 unit sets of size `>= 2` against 1,322 class-vertex pairs), so
    19,002 must not be appended to the sequence `2, 5, 7, 16, 31, 68, 207, 1333` of HYP-9658 (D).
15. **`below_rays.txt` is not checkable from the shipped files.** `verify_lines.py` checks only the arithmetic of the
    certificates; "census" needs the class lists of `classgen.c`, which are not shipped, and the factors 1771, 685
    are not exhibited. The audit's `out/values_witness.txt` gives a witness (order, code) for every value of order
    `<= 9`; adding the witness code to each line of `below_rays.txt` would make the certificate self-contained.
16. **Stale docstrings:** `verify_s21.py` says "THM-5120" and "about 3-5 minutes"; `finalize.py` calls itself
    `finalize2.py`.
17. **Pre-registration.** Both texts are byte-identical in HEAD to the committed versions (`2c27353a14`, 11:59:54;
    `050845a98d`, 12:28:11), both commits are on `origin/main`, and outcomes were appended later (`b576542d6c`,
    `7035dd3c6e`). P1 and P4 are reported as FAILED, faithfully. Only commit order can be checked, not wall-clock
    order of the runs. Q1-Q4 were confirmed with wide margins (11 of 11 against "at least 6"; 510 against "at least
    450").

---

## What was verified independently

### Proofs (read line by line)

- **Tail law (A).** The case split (`Z_1` empty / nonempty; `L_2` empty / nonempty) is disjoint and exhaustive. Every
  junction into `Z` is automatic because `B => Z`; the only conditions are `z_j -> v` (`y_j = 0`) and, when `L_2` is
  empty, `v ->` first vertex of `Z_2`, which depends on `y` only. The conventions `H(empty) = 1` and `L = {v}`
  (`H = E = d = U = 1`) are consistent. Correct.
- **Offset formula.** Covers with `P2 != (v)` split into: `P1` empty (paths from `v`), last vertex of `P1` beating
  `v` (paths with `v` interior), and crossed covers. Correct. Cut form: a crossed cover is an ordered two-path cover
  of `B` ending in `O` and starting in `O`; a cover with `P2 != (v)` is one of `B` whose second path starts in `O`.
  Correct. `c0` odd and `c0 <= H(B)` with equality iff `|O| <= 1`: correct (minor 6).
- **(B1)-(B4).** Correct given THM-3121 (4) and `P_2(TT_a) = 2^a - 2`. The seed of `L_3` (`B` a 3-cycle,
  `W_01 = 0 + 2 = 2`) and of `L_5` were recomputed by hand and by machine.
- **(C).** Junctions lie on distinct pairs, each containing a new vertex; range `1 <= j <= k + 1`; equality "for all
  `j`" is needed (`P_n = n!` always) and is what is stated; the minimiser is unique because the coefficient of
  `P_1` is positive. Correct.
- **(D).** Correct (`(x - 1) + ((x - 1) mod 2) = 2 floor(x/2)`).
- **(E).** The Camion walk is valid: the inserted copy of `y` is separated from the other copy by `x`, and
  `y -> x` excludes `y = q_1`. Correct; it also follows from THM-5124 (2), because a twin of a vertex of `B` cannot
  be a one-triangle vertex of the rigid configuration with hubs in `B`. So the credit to S90 is right.

### Recomputation (own code)

| item | universe | result |
|---|---|---|
| classes of orders 1-9 | own canonical form; counts and `sum n!/|Aut| = 2^C(n,2)` | 1, 1, 2, 4, 12, 56, 456, 6880, 191536; mass formula exact at every order |
| spectrum | every class of order `<= 9`; and all `2^28` labelled tournaments of order 8 without canonical forms | 1, 1, 2, 3, 7, 19, 77, 320, 1520 values; maxima 661, 3357; 7 and 21 absent to order 9; every odd `N <= 609` other than 7, 21 at order `<= 8`; least non-value of order `<= 9`: 2883 |
| identities `c0 = d - delta'`, `M = H - E + delta'`, `c0` odd | 1,782,476 roots (every root, every class, orders 2-9), `delta'` from its definition | no failure |
| cut form (1') | every `(B, O)` with `|B| <= 7` and every unit set over `|B| = 8, 9`, on the literal `L` (82,771 cases); 809 rooted tournaments from the definitions by permutations | no failure |
| `W_ef` even, `>= 2 H(B - e - f)` | every pair of every class of order `<= 9` | no failure |
| tail law | 8,080 triples `(L, v, y)`, state from the definitions | no failure |
| unit-offset roots, orders 2-9 | every class, every root | 2, 5, 7, 16, 31, 68, 207, 1333; equal to the packet's file as multisets of (order, class, H, d, E, U, M), 1,669 rows; `d = E`: 116 of 336 and 197 of 1,333; out-degree `<= 2` to order 7, 3 from order 8 |
| unit sets over bases of order 9 | every `O` of every size | 2,771 pairs, 16,231 triples, nothing larger (19,002); 140 and 1,007 pairs (out-degree, `M`) at orders 8, 9 |
| odd `c` with a seed | tail `<= 9`: 636 (module seeds add none); with order ten: 5,529; `c <= 1025`: 377 -> 510 of 513 (missing 245, 853, 909); `c <= 4097`: 2,016 of 2,049 | all as claimed; the list of 33 agrees |
| `line_seeds.txt` | all 2,016 lines: `(H, d, E, U)`, `c0 = 1`, `M = c 2^(s0-1)`, cut form, `H(L o 0^t)` for `t = 1, 2, 3` literally (orders to 13); module seeds by walk counts and `a = 2..5` | 6,052 literal counts, 0 bad; `(c, s0, order)` equal to the own census for every `c <= 4097` |
| required sample | `c = 3, 5, 7, 13`, ten random lines, the 4 module seeds, **all eleven** order-ten seeds at `t = 1, 2` | pure Python, subset DP and DFS enumeration of the paths: all agree (`aud_checks.out`) |
| `s0` table of (B3) | typed in from the theorem | equals `line_seeds.txt` and the own census (65 entries) |
| `below_rays.txt` | 121 points = `sum (s0 - 1)`; 117 census, 2 products, 2 gaps | every census value has an own witness of order `<= 9`, recounted in Python; `1771`, `685` have witnesses of order 9; the joins of orders 12 and 13 have 5313 and 3425 paths (DP and DFS) |
| doubling vertices, orders 3-8 | every vertex of every class, walk counts `c_1..c_5` and literal `K[TT_a at x]` | 48 nontrivial values (+ the value 1); the list `N = 3 mod 4, N < 200` as claimed; 15 and 11 absent |
| two-module families | 406,412 triples, walk counts, literal spot checks | 3,477 distinct exact tuples; none additive (no triple at all has `de = 0` with `be, ga != 0`); `de > 0` iff both on cycles; weight-three rows `d <= 4` and `d = 7` complete |
| (C) | `k = 1`: all labelled `n <= 6`; `k = 2`: `n <= 5`; `k = 3`: `n <= 4` and 6 of order 5; `k = 4`: `n <= 3`; `k = 5`: `n <= 2`; `k = 6`: `n = 1`; `P_j` also from the definition | no failure; bound and equality case hold |
| (D) | literal counts of all 8,190 words of length `<= 12` | no failure |
| (E) | construction carried out on 17,452 ordered pairs | every walk valid |
| packet scripts | `verify_lines.py`, `verify_s21.py`, `unitpairs.py`, `mean_pi_identity.py` rerun | exit 0; recorded outputs reproduced |

### Scope and wording of the headline (items 4 and 5 of the brief)

- **Headline.** *"every odd N > 1 whose N - 1 has odd part at most 129 is a Hamiltonian-path count, except 7 and
  21"* is precise and verified. It is safe.
- **"Density zero".** Correct and correctly placed (THM status, Scope, README header and section 5): the set has
  about `65 log2 X` members up to `X`. No sentence of the theorem or the hypothesis reads as progress on density or
  completeness; the two README sentences that do are M5.
- **"Equivalent".** *"Chan-Pak is equivalent to: every line L_c is in the spectrum, with one exception on L_3 and
  one on L_5"* is true, but it is a relabelling (`N <-> (c, s)` is a bijection) and carries no content; the README
  says as much ("It reorganises the target").
- **"Isolated".** (B2) is the exact statement `L_3 n S = L_3 - {7}`, `L_5 n S = L_5 - {21}`; verified. The word
  describes the two lines only.
- **Numbers.** 65 of 65 lines, 61 + 4 seeds, 121 = 117 + 2 + 2 points, 636, 1,333, 116 of 336, 197, 406,412, 3,477,
  19,002 = 2,771 + 16,231, 5,529, 2,016 of 2,049, 510 of 513 and the `s0` table are all confirmed. The three that
  are not: "8 of the 336", "eight of them for c <= 257", and "49 of the 320" (minors 1-3).
- **Scope note for the order-ten run.** Present in PREREGISTRATION.md, the HYP status and outcomes, and the README
  table; absent from the THM status, README section 2 and the HYP title and (B) (M4).

### Novelty search

Targeted searches of `01-canon/theorems` and of the cited result folders for: offset / crossed covers / cut
coordinates (found: M1), mean over extensions (found: M2), substitution law (found: M3), "unit offset", "2-adic
line", "doubling vertex", `1 + 3 * 2^s`, weight three, the `(x+1)`-trajectory. No earlier statement was found of:
the 65 complete lines, the isolation of 7 and 21 on `L_3` and `L_5` for all `s`, unit-offset roots as a notion, the
general-`k` mean identity, or the trajectory form of (D). Known neighbours, correctly not claimed: THM-5109 gives
`1 + 2^b + 2^K` for fixed small `b` and all large `K` (columns of the weight-three table; (B4) gives rows);
THM-5085 (3) gives rays with other offsets (`2^K - 159`).

---

## What could not be verified

- Doubling vertices of order ten (9,733,056 classes): not searched by the packet or by the audit (M4).
- THM-1370 is taken as known (its consequence was checked to order nine only).
- The repository is large; the novelty search covers the canon and the folders named in the brief plus those the
  trail led to (Sessions 216, 229, 244, 245, 246; THM-4111/4115, 4960, 3121). Result notes elsewhere may contain
  further precedents.
- `classgen.c`, `seeds.c`, `dbl.c`, `unitsets.c`, `modfam.c` were not rerun; their outputs were compared with the
  own census instead. `finalize.py` cannot be rerun without the unshipped class files.
- "2,256 doubling vertices of orders 3-5" (a labelled count in `verify_s21.py`) was reproduced only by rerunning
  that script.
- Wall-clock order of pre-registration and runs (only commit order).
- README sections 6.4, 6.5 and 8 are analogies and a log; they were read, not tested.

---

## Appendix: short proofs supplied by the audit

**A1. `W_ef >= 2 H(B - e - f)`.** Let `r_1 ... r_k` be a Hamiltonian path of `R = B - e - f`. Insert the block
`e f` or `f e` at a position `i` in `0..k`; the block need not respect the arc between `e` and `f`. The result is a
two-path cover counted by `W_ef` iff the two outer junctions are arcs. Put `a_i = [r_i -> e]`, `c_i = [f -> r_i]`,
`(a_0, c_0) = (1, 0)`, `(a_(k+1), c_(k+1)) = (0, 1)`. Then `e f` fits at `i` iff `a_i = 1` and `c_(i+1) = 1`, and
`f e` fits iff `c_i = 0` and `a_(i+1) = 0`. With the potential `phi(10) = 0`, `phi(11) = phi(00) = 1`,
`phi(01) = 2` on the states `(a_i, c_i)`, each step contributes at least `phi(next) - phi(current)` fits, and the
walk goes from `10` to `01`: at least two fits. Different paths or positions give different covers.

**A2. `c_3 = 0` forces `c_j = 0` for `j >= 4`.** A walk `S_0 x S_1 x S_2 x ...` with three or more visits has two
consecutive loops `S_1, S_2`; the tournament on `{x} u S_1 u S_2` is strong, so it has a Hamiltonian cycle (Camion),
and replacing `x S_1 x S_2 x` by that cycle read from `x` to `x` gives a walk with one visit fewer.

**A3. 15 has no doubling realiser of any order.** `H` is multiplicative over strong components, and for `x` in a
component `Q`, `c_j(K, x) = c_j(Q, x) H(K)/H(Q)`. If `H(K) = 15`, either the non-trivial components are `C_3` and
the strong 4-tournament, where `c_2(K, x)` is `5` or a multiple of `3`, never `7`; or one component has `H = 15`,
and then its order is 5 or 6 because the strong floor is at least 25 from order 7 on (THM-1370). That component is
inside the order-`<= 8` census, which has no doubling vertex with `H = 15`.

**A4. The odd-order mean identity (README 6.2) for all `n`.** An ordered `j`-tuple of disjoint nonempty paths
covering `[n]` is a permutation cut into `j` blocks whose `n - j` inner junctions are arcs, so
`E[P_j] = n! C(n-1, j-1) 2^(j-n)` and
`sum_n E[pi_T(x)] t^n / n! = sum_j C(x, j) (t / (1 - t/2))^j = ((2 + t)/(2 - t))^x`. The right-hand side has
exponential generating function `exp(x sum over odd l of 2^(1-l) t^l / l) = exp(2x artanh(t/2))`, the same function.
Checked in closed form for `n <= 12`, `x <= 5`, by brute force over permutations for `n <= 7` and over all labelled
tournaments for `n <= 4`.

---

## Files in this folder

- `aud_engine.c` - census engine (`gcc -O2 -o aud_engine aud_engine.c`). Modes: `census OUT` (39 s), `seeds FILE`,
  `roots FILE`, `labelled8` (24 s), `mean` (41 s), `oneapex 12`, `dbl9` (75 s), `exact4` (10 s).
- `aud_checks.py` (11 s) - brute force from the definitions: T1-T7. Output `aud_checks.out`.
- `aud_compare.py OUT` - packet claims against the census. Output `aud_compare.out`. Its two FAIL lines are minor
  findings 1 and 3; every other line passes.
- Outputs: `aud_seeds.out`, `aud_labelled8.out`, `aud_mean.out`, `aud_oneapex.out`, `aud_exact4.out`,
  `aud_dbl9_extension.out`.
- `out/`: `census_summary.txt`, `spectrum_order_8.txt`, `spectrum_order_9.txt`, `values_witness.txt` (value, least
  order, code), `roots_unit_offset.txt`, `seed_table_c_le_4097.txt`, `doubling_vertices_orders_3_to_8.txt`,
  `doubling_vertices_order_9.txt`, `exact_tuples_order_le8.txt`. The full list of unit sets (20,677 lines, sha256
  `5071a803f6ff0baa43e5ca4c8f9e316e8e627428c4b9b4717e2db4084c692a35` with LF endings) is regenerated by
  `aud_engine census OUT` and is not stored.
