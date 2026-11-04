# Audit 2 of THM-4868 (G) and (H): admissible Shelor parents, run transitions, and what is built on them

**Scope.** An independent, adversarial audit of the new parts (G) and (H) of THM-4868 (kind-pasteur-2026-11-04-S10), and of
the derivations built on them: the bijection, Garey's count, the coupled recursion, the elimination to the run-length
recurrence, the generating function, and HYP-9458's closed forms. A separate auditor covers (A)-(F). They are touched
here only where (G) and (H) consume them: (C), (E), Lemma 1 and Lemma 2 as inputs, and (D) for the passage to T(n, k).

**Texts audited.**
- the note `05-knowledge/results/uniquely_hamiltonian_noncritical_run_20261104.md`;
- `01-canon/theorems/THM-4868-*.md`;
- HYP-9460, and for context HYP-9458, the S74 note's section 4, and THM-4862 (b).

**Moving text.** The author edited the note, the canon file and both HYP files throughout the audit, apparently while
applying the other auditor's report. Several of my draft findings were fixed while I worked; they are listed as RESOLVED
in section 4. Locations below quote the text. Line numbers refer to the note as it stood at about 13:20 (488 lines).

**Independence.** All code in this folder was written for this audit. No program or output of the author was used as
evidence. The author's `confirm.out` was read once, as a bookkeeping check of the pre-registration record (section 3).

## Verdict: SOUND WITH CORRECTIONS (minor)

The proofs of (G) and (H) are correct, and so is every derivation built on them. Each step was checked line by line, and
I found no logical gap.
- Every finite claim was re-verified with my own code:
  - labelled exhaustive enumeration, with no canonical forms, through order 9;
  - two independent class generators through order 16;
  - (G), (H), the bijection and the key path-count identity on every class of order 4-15, with all children of order
    5-16 measured directly.
- The open corrections are minor:
  - the closed-form induction lacks ranges, and two displayed one-step formulas are false at their first order;
  - the section-8 typing rows are inconsistent with the corrected status block;
  - wording.

## 1. Proof audit, line by line

Notation as in the note: T is uniquely Hamiltonian (UH) of order m, y is its run start, U = T - y, p_0 = y-, and
P_p = C_U - (p -> p+).

**Lemma 3** (hc(U^p) = e_p(U)). Correct for every tournament U with |U| >= 2.
- The cycle enters y only from p and leaves to some z != p. The rest is a Hamiltonian path z ... p, and conversely y
  beats every z != p.
- Verified for every class U of order 4-14 and every vertex p: the exact cycle count of the (m+1)-vertex tournament U^p
  equals e_p(U). At order 15 it was checked in capped form, min(hc, 2) = min(e_p, 2). It was also verified exactly on
  every labelled UH tournament of order <= 9.

**(G), sole entries are admissible.** Correct. The current text separates the two equalities correctly: the in-neighbour
observation, then Lemma 3.

**Base, order 4.** Checked by hand.
- The Hamiltonian paths are 0123, 1230, 1302, 2301 and 3012, so e = (1, 1, 2, 1).
- NC = {1, 2}, y = 1, y- = 0 is a sole entry, and adm = {0, 1, 3}.
- The uniqueness of S_4 is now stated. My enumeration agrees: 24 labelled UH 4-tournaments, one asymmetric class.

**Step, set-up.** Correct.
- y is non-critical, so U is UH by Lemma 1, with C_U = C_T - y + (p_0 -> y+).
- y is a sole entry with in-neighbour p_0, so T = U^(p_0).
- Hence e_(p_0)(U) = hc(T) = 1.

**Classification** e_p(T) = [e_p(U) - hp_U(p_0, p)] + W(p_0, p) for p not in {y, p_0}. Every path is counted exactly once.
- Either y is first: the rest is a Hamiltonian path of U ending at p whose start z satisfies y -> z, that is z != p_0.
- Or y is not first: its predecessor is its only in-neighbour p_0. y is not last, since p != y. The path is A, y, B with
  A ending at p_0 and B non-empty and ending at p, and y -> first(B) holds automatically because p_0 is in A.
- The two cases are disjoint, and each correspondence is bijective.
- For p = p_0, y has no possible predecessor, since p_0 is last, so y is first. For p = y, the path is (path ending at
  p_0), y. Both are right.
- Verified numerically on every class of order 5-15 and every p not in {y, p_0}, with W computed independently as
  sum over S of E_S(p_0) E_(V-S)(p) by a subset DP (`decomp = 0`).

**W >= 1, and the case p+ != p_0.** Correct.
- P_p contains p_0 != p. Cutting it after p_0 gives the two non-empty paths p+ ... p_0 and p_0+ ... p.
- If p+ != p_0, then P_p starts at p+ != p_0 and is counted in the bracket. So e_p(T) >= 2.
- Edge cases: p = y+ (= p_0+ on C_U) lies in this case and is harmless for |U| >= 3. p = p_0- = y-- is exactly the
  excluded case.

**p = p_0-, p_0 a sole entry.** Correct.
- p_0's only in-neighbour is p, which is last, so p_0 is first. y's only in-neighbour is p_0, so y is second.
- The rest is a Hamiltonian path of U - p_0 ending at p_0-, and conversely.
- So e_p(T) = hc(U) = 1, by Lemma 3.
- "A sole entry of U iff of T" is right, because y does not beat p_0.

**p = p_0-, p_0 not a sole entry** (the only use of the induction hypothesis, (G) for U at order m - 1 >= 4). Correct.
- U has one or two sole entries: y_U, and possibly y_U-. Two sole entries are adjacent, and the second cannot be y_U+,
  since that would make y_U critical (Lemma 2). So the admissible non-sole-entries are {y_U-} if y_U- is not a sole entry,
  and {y_U--} if it is.
- Hence p = y_U-- (not admissible, because y_U- is not a sole entry) or p = y_U---.
- y_U--- lies outside {y_U, y_U-, y_U--} exactly when |U| >= 4; at |U| = 3 it would be y_U. So |U| >= 4 is precisely what
  is needed, and it holds because m >= 5.
- "Not admissible, so e_p(U) >= 2" silently uses e_p(U) >= 1, which P_p gives (F8).
- The current text gets W >= h from P_p being one of the h >= 1 paths from p_0 to p, so e_p(T) >= e_p(U) >= 2.
- Measured at orders 5-15: min (W - h) = 0, so W >= h is tight. In this case min e_p(U) = 3 (2 at order 5) and min h = 1.

**Bijection.** Correct.
- Phi: [T] -> [(T - y_T, y_T-)] is well defined, because y_T has an isomorphism-invariant description: the unique
  non-critical sole entry, by (C).
- Psi: [(U, p)] -> [U^p] is well defined.
- Psi(Phi) = id, since T = (T - y_T)^(y_T-).
- Phi(Psi) = id: the new vertex of U^p is a non-critical sole entry (U is strong), hence y of U^p, and its predecessor is
  its only in-neighbour p.
- The count needs distinct vertices of one class to give distinct children. If U^p is isomorphic to U^p', the isomorphism
  sends y to y', so it restricts to an automorphism of U sending p to p'. By (E), p = p'. The note's "By (E), pairs
  (U, p) of a class and a vertex are distinct" says this in compressed form.
- m >= 4 is needed: at m = 3, the three pairs (C_3, p) all give S_4, because C_3 is symmetric.
- Verified: every class of order 5-16 is hit by exactly one admissible pair from the order below, with no misses and no
  repeats. Also verified without canonical forms, on labelled tournaments (section 2).

**Sole-entry rule.** Correct. In U^p, y has in-degree 1, p's in-degree is unchanged, and every other vertex gains 1 and
so has in-degree >= 2.

**Count.** Correct.
- u(m + 1) = sum of |adm(U)| = 2u + b, because |adm| = 2 + [two sole entries].
- b(m + 1) = #{(U, p) : p a sole entry} = u + b, because sole entries are admissible.
- Hence u(m + 1) = 3u(m) - u(m - 1) for m >= 5, with u(4) = 1 and u(5) = 3. So u(m) = F(2m - 6) and b(m) = F(2m - 7).
- Measured: |adm| = 2 for F(2n - 8) classes and 3 for F(2n - 7) at every order 4-15.

**(H): C_T is C_U with y inserted after p.** Correct, and the current text now gives the reason (e_p(U) = 1).

**Fact 1.** Correct.
- T - v = (U - v)^p.
- U - v is strong, and y has an in-neighbour p (since v != p) and out-neighbours (|U - v| >= 3).
- Adding such a vertex preserves strong connectivity.

**Fact 2.** Correct: y is a sole entry with y- = p, so p is critical by Lemma 2.

**Fact 3.** Correct.
- x_U is not admissible. It differs from y_U (r >= 2), and x_U is non-critical while y_U- and, when relevant, y_U-- are
  critical.
- So y -> x_U. Its out-neighbourhood stays {x_U+}, its C_T-successor is still x_U+, and x_U+ is critical in T by Lemma 2.

**Run identification.** Correct. By (C) for T, NC(T) is one contiguous run starting at y. The C_T-segment from y to x_U
is:
- y, y_U+, ..., x_U for p = y_U;
- y, y_U, ..., x_U for p = y_U-;
- y, y_U-, y_U, ..., x_U for p = y_U--.

All of it is non-critical, and x_U+ is critical. So the run is exactly that segment, and no other critical vertex of U
can turn non-critical. This answers "why the run cannot extend elsewhere": by contiguity from y.

**y_U- becomes non-critical when p = y_U--.** Correct. The current text uses the direct Hamiltonian cycle
y -> y_U -> (C_U-path) -> y_U-- -> y of T - y_U-.

**Coupled recursion.** Correct for m >= 4 and every integer r.
- The children of U have runs r(U), r(U) + 1 and r(U) + 2, the last only when U has two sole entries.
- The children with two sole entries are exactly those whose p is a sole entry: p = y_U always, with run r(U); and
  p = y_U- when U has two sole entries, with run r(U) + 1.

**Elimination.** Correct, n >= 6.
- The second line at m = n - 2 (r -> r - 2) and the first at m = n - 2 (r -> r - 1) give
  N_2(n - 1, r - 2) = N(n - 1, r - 1) - N(n - 2, r - 1).
- Substituting into the first line at m = n - 1 gives the recurrence for every integer r.
- The recurrence puts no mass outside 2 <= r <= n - 1 (checked).
- The initial rows N(4, .) = (1, 0) and N(5, .) = (1, 1, 1) follow from the first line at m = 4, with N_2(4, 2) = 1.
- The recurrence rows equal the coupled-system rows at every order 4-18. So do the totals: F(2n - 6), and F(2n - 7) for
  N_2.

**Generating function.** Correct.
- Multiply G by 1 - (1 + 2t)x + t x^2. The x^4 coefficient is t^2, the x^5 coefficient is
  (t^2 + t^3 + t^4) - (1 + 2t) t^2 = t^4 - t^3, and every higher coefficient vanishes.
- sympy solves the coupled system symbolically and returns exactly G; the series agrees through x^18; and G(x, 1) =
  x^4 / (1 - 3x + x^2).
- The sum runs over n >= 4, which is unstated; C_3 is excluded (F4).
- Audit-derived companion, verified: the sum of N_2(n, r) t^r x^n is t^2 x^4 (1 - t x) / (1 - (1 + 2t) x + t x^2).

**Closed forms** (end of the proof of (H), note l.336-340). The formulas are right, but the induction is stated without
ranges, and two displayed steps are false at their first order (F2).
- The T-recurrence follows from (D) and the recurrence for N.
- T(n, 2) = 2 T(n - 1, 2) holds for n >= 6, since T(., 1) = 0. This gives 2^(n-5) for n >= 5.
- "T(n, 3) = 2^(n-7) + 2 T(n - 1, 3)" uses T(n - 2, 2) = 2^(n-7), which needs n >= 7.
  - At n = 6 it gives 5/2, while T(6, 3) = 3.
  - With the base T(6, 3) = 3 it holds for n >= 7.
- "T(n, 4) = n 2^(n-9) + 2 T(n - 1, 4)" is false at n = 7 (23/4 against 6). It holds for n >= 8, from the base
  T(7, 4) = 6.
- The diagonal recurrences hold for n >= 6, from bases at n = 5.
- The final forms hold exactly in HYP-9458's ranges (n >= 5, 6, 7). Verified for every n <= 60 from the recurrence, and
  symbolically: each closed form satisfies the T-recurrence. Below those ranges they fail (n = 4; n = 4, 5; n = 4-6).
- Cross-check of the two independent proofs: S74's audit closed form, from Shelor's grammar,
  T(n, k) = sum over d of C(n - 4, k - 3 + 2d) C(k - 2 + d, d),
  equals the (H) recurrence for every 5 <= n <= 40. So do its next diagonal n(n - 4)(n - 5)/6 and the recurrence's
  N(n, 5).

**Inputs from (A)-(F).** (G) uses:
- (C): y is the unique non-critical sole entry; there are at most two sole entries, and if two, they are y- and y;
  Lemma 2;
- (E), for the count;
- Lemma 1, which uses Camion;
- through (C), Moon 1966.

(H) uses (C) and Lemma 2. The HYP-9458 passage uses (D). I re-checked only the parts of (C) and (E) that (G) and (H)
consume: the run structure, the uniqueness of y, the adjacency of two sole entries, and the rotation argument for (E).
All are correct, and all are verified on every class below (`C-input = 0`, `asym = 0`, `cut = 0`).

## 2. Independent computation (own code, this folder)

**`labelled_audit2.c`: no canonical forms anywhere.**
- Mode A enumerates all 2^C(n,2) labelled tournaments, n <= 7.
- Mode B enumerates all tournaments containing the fixed cycle 0 -> 1 -> ... -> n-1 -> 0, n <= 9. This covers every
  class, by relabelling along its cycle.
- Uniqueness of the cycle is decided by a DFS count. NC comes from direct strong-connectivity tests, and e_p from a
  subset DP.
- On every UH tournament it checks:
  - the (C)-inputs;
  - (G);
  - Lemma 3 for every p (exact DP count on T^p);
  - (H): the child is UH, its run length is r + 0, 1, 2, its run starts at the new vertex with predecessor p, and the
    proof's facts hold. Those facts are: p is critical; x != p; x is still a sole exit; x+ is critical; NC(U) - p stays
    non-critical; and y_U- becomes non-critical for p = y_U--.
  - the labelled bijection identities L(n+1) = (n+1) * sum of |adm| (mode A) and n L_0(n+1) = (n+1) * sum of |adm|
    (mode B).

| run | tournaments | UH | expected | bad counters | sum of \|adm\| | bijection identity |
|---|---:|---:|---|---|---:|---|
| A 4 | 64 | 24 | 1 * 4! | 0 | 72 | 5 * 72 = 360 = L(5) |
| A 5 | 1,024 | 360 | 3 * 5! | 0 | 960 | 6 * 960 = 5,760 = L(6) |
| A 6 | 32,768 | 5,760 | 8 * 6! | 0 | 15,120 | 7 * 15,120 = 105,840 = L(7) |
| A 7 | 2,097,152 | 105,840 | 21 * 7! | 0 | 277,200 | — |
| B 4-9 | 4; 32; 512; 16,384; 1,048,576; 134,217,728 | 4, 15, 48, 147, 440, 1,296 | n F(2n - 6) | 0 | 12, 40, 126, 385, 1,152, 3,393 | n L_0(n+1) = (n+1) * sum of \|adm\| at n = 4-8 |

The labelled run-length rows divided by n! (mode A) or by n (mode B) are the note's section-4 rows, through order 9.

**`classes_audit2.c`: classes by two generators.**

Generation:
- Order 4 is exhaustive.
- **Method A** forms U^p for every class U and every vertex p, and keeps it when an exact DP count of the cycles of the
  (m+1)-vertex tournament gives 1. Its completeness alone would need Shelor or (C).
- **Method B** forms U + v with an arbitrary out-set (2^m - 2 choices) and keeps it when the sum of hp(w, u) over
  w in O and u not in O equals 1.
  - hp comes from start-fixed subset DPs.
  - Its completeness uses only Moon 1966 (with Camion) and the trivial extension lemma, not (C) or (G). T - v strong
    and hc(T) = 1 force hc(T - v) = 1.
  - At orders <= 12 each hit is re-verified by an exact DP count.
- Classes are identified by the minimum rotation word along the unique cycle, a complete invariant for UH tournaments.
  As a sanity check, all classes of order 4-9 also get brute-force canonical forms over all n! relabellings.

Per class of order 4-15, the program computes:
- hc = 1 (DP), and asymmetry (all n rotation words distinct);
- the (C)-inputs;
- (G) via exact e_p;
- Lemma 3 for every p (exact at orders 4-14, capped at 15);
- every admissible child: its run length, the proof's facts, its canonical identity and hit counts, and its own run, N_2
  and cut-arc counts;
- the path-count identity with an independent W;
- the cut arcs, by deleting each C-arc.

| order | method A = method B | count = F(2n - 6) | other checks |
|---:|---|---|---|
| 5-7 | yes | 3, 8, 21 | equal to the exhaustive labelled class sets |
| 8-12 | yes | 55, 144, 377, 987, 2,584 | B hits re-verified by DP: 0 bad |
| 13-15 | yes | 6,765, 17,711, 46,368 | — |
| 16 | yes | 121,393 | A 259 s, B 210 s, on a shared machine; from `classes_audit2_n16_partial_stopped.out` |
| 4-9 | — | — | brute-force canonical forms: all classes distinct |

| parent order n | bad counters (hc, asym, C-input, G, L3, H, facts, canon, decomp, pathsplit, cut) | \|adm\| = 2 / 3 | child deltas 0/1/2 | bijection onto order n + 1 |
|---:|---|---|---|---|
| 4 | all 0 | 0 / 1 | 1 / 1 / 1 | 3 of 3 hit exactly once |
| 5-13 | all 0 | F(2n - 8) / F(2n - 7) | F(2n-6) / F(2n-6) / F(2n-7) | every class hit exactly once |
| 14 | all 0 | 6,765 / 10,946 | 17,711 / 17,711 / 10,946 | 46,368 hit once, 0 missed, 0 repeated |
| 15 (top-order mode) | all 0 (L3 in capped form) | 17,711 / 28,657 | 46,368 / 46,368 / 28,657 | 121,393 admissible pairs give 121,393 distinct children, which equals \|UH(16)\| from the A = B generation; 0 repeated |

**Rows** (`recurrence_audit2.out`, part 5).
- The measured N, N_2 and T rows at orders 4-15, and the children's rows at order 16, equal the recurrence and the
  coupled system.
- They equal the note's section-4 table, and HYP-9458's published T(n, k) rows at orders 5-13 and the full order-16 row
  2048, 8192, 16896, 23552, 24600, 20256, 13504, 7380, 3312, 1210, 352, 78, 12, 1.
- P3-P6 are measured at orders 14, 15 and 16.
- Every measured cut-arc row is the run row reversed, so #cut + r = n + 1 on every class (D), and on every child.
- The recurrence reproduces HYP-9460 P7's rows at orders 14-18. Rows 17 and 18 are arithmetic only; my census reaches
  order 16.

**`recurrence_audit2.py` (sympy).** It checks:
- the elimination;
- the generating function: series, symbolic solve of the coupled system, and t = 1;
- the companion generating function for N_2;
- all quoted tables;
- the closed forms for n <= 60, and symbolically;
- the failure orders of the note's one-step forms.

Output: `ALL OK`.

## 3. Typing

- **"Garey's count re-derived".** Correct as typed: KNOWN, re-derived. The route needs neither Garey nor THM-4862 (b),
  which is re-proved as Lemma 3. It is self-contained modulo two classical citations, Moon 1966 (NC non-empty, used in
  (B) and (C)) and Camion 1959 (Lemma 1). "Self-contained" should name them (F6).
- **"HYP-9458 proved".**
  - The body and canon now type (G) and (H) correctly as a second, Shelor-free proof.
  - The first proof is opus-S74's audit, given Shelor's grammar. Git confirms its order: `a0518f53cd` at 12:51:03 -0700,
    after this session's pre-registration `a8db6218f9` at 12:28:48 -0700.
  - The note's title still reads "(HYP-9458 proved)" (F7).
  - The derivation of the formulas is complete given (D), which the other auditor covers; my census confirms (D)'s count
    through order 16.
- **"Canonical".** Justified. The deletion vertex y is defined isomorphism-invariantly, as the unique non-critical sole
  entry, so T = (T - y)^(y-) commutes with isomorphisms.
- **Priority and dictionary.**
  - The priority paragraph is consistent with git.
  - The "dictionary" (Shelor's j = 0, 1, 2 corresponds to {y_U, y_U-, y_U--}) interprets a CITED grammar whose source was
    not read here.
  - It is exact only with the proviso that j = 2 requires y_U- to be a sole entry; otherwise y_U-- is not admissible
    (F7).
- **Silent reliance.** None on unproved results. The one small unstated step left is e_p(U) >= 1 (F8). S_4's uniqueness
  is now stated.
- **FINITE-EXACT scope.** The status block and the canon now list exactly what each program checks. The section-8 typing
  row "all of the above at orders 4-18 | FINITE-EXACT" still claims more (F1).
- **Las Vergnas.** The canon now types it as corroborating only. The note's section-8 typing row still says "the top
  diagonal uses Las Vergnas" (F5).
- **Bookkeeping, resolved during the audit.** Section 5 of the note, HYP-9460's results and status, and its `note:` field
  were placeholders or stale when I started. All are now filled. As a bookkeeping check only, the author's
  `confirm.out` rows at orders 14-18 equal the pre-registered P7 rows.

## 4. Findings (against the text at about 13:20)

| # | severity | location | problem | fix |
|---|---|---|---|---|
| F1 | minor | note section 8, typing row "all of the above at orders 4-18 \| FINITE-EXACT" (l.472) | The row contradicts the corrected status block (l.40-45) and the canon. (G) and (H) were machine-checked by the author only on parents of order 4-14. The (F) bijection, the rooted count and (B)'s "equivalently" form were not checked by the author's programs. | Use the status block's wording: "(A), (C)-(E), the degree-1 and (F) counts, and one-arc cuts at 4-18; (G) and (H) on parents 4-14". Optionally add "(G) and (H) on parents 4-15, children to 16, by both audits". |
| F2 | minor | note, end of the proof of (H) (l.336-340): "T(n, 3) = 2^(n-7) + 2 T(n - 1, 3); T(n, 4) = n 2^(n-9) + 2 T(n - 1, 4)" | There are no ranges or bases. The T(n, 3) step is false at n = 6 (5/2 against 3), and the T(n, 4) step is false at n = 7 (23/4 against 6). The ranges of the final forms are not stated in section 1 (H) or here. | State: T(n, 2) = 2T(n - 1, 2) for n >= 6; the T(n, 3) step for n >= 7 from T(6, 3) = 3; the T(n, 4) step for n >= 8 from T(7, 4) = 6; the diagonal steps for n >= 6 from n = 5. Hence 2^(n-5) for n >= 5, n 2^(n-7) for n >= 6, and (n^2 + n - 8) 2^(n-10) for n >= 7, which are HYP-9458's ranges. |
| F3 | (resolved) | base of (G) | The uniqueness of S_4 (needed for the induction base and u(4) = b(4) = 1) was unstated. | Fixed during the audit (l.259-261). |
| F4 | minor | canon (H) (l.97): "The children T^y, T^(y-) and T^(y--) have run lengths r, r + 1 and r + 2"; canon (G) (l.95): "bijection from UH(n + 1) onto the admissible pairs of UH(n)"; note (G), l.122-123; the generating function (note l.395, HYP-9460 P7) | T^(y--) is uniquely Hamiltonian, so its run is defined, only when y- is a sole entry. The bijection is between isomorphism classes, with the parent order >= 4; at order 3 the pairs (C_3, p) all give S_4. The generating-function sum is over n >= 4. | Add "(the last when y- is a sole entry)", "isomorphism classes", "m >= 4" and "n >= 4". |
| F5 | remark | note section 8, typing row (D): "the top diagonal uses Las Vergnas 1975 (CITED ...)" (l.467) | (H) proves N(n, 2) = 1 outright: N(n, 2) = N(n - 1, 2) for n >= 6, and N(4, 2) = N(5, 2) = 1. The note's (D) statement and the canon already say so; the typing row does not. | "PROVED by (H); Las Vergnas (CITED) corroborates." |
| F6 | remark | note (G), "by a self-contained route" (l.124); canon (G) (l.96), "re-derived without Douglas's characterization" | The route still uses Moon 1966 (through (B) and (C)) and Camion 1959 (Lemma 1). | "Self-contained modulo Moon 1966 and Camion 1959". |
| F7 | remark | note title "(HYP-9458 proved)"; "Dictionary" (l.39); canon status (l.15-16: "parent letters j = 0, 1, 2 are (G)'s y, y-, y--") | The title reads as a first proof, although the body credits S74's audit with the first. The dictionary lacks the proviso on j = 2. | Retitle as "a second, Shelor-free proof of HYP-9458". Add "j = 2 only when y_U- is a sole entry". |
| F8 | remark | note (G), case p = p_0-, p_0 not a sole entry: "This is not admissible in U (\|U\| >= 4), so e_p(U) >= 2" (l.287) | Uses e_p(U) >= 1 without saying so. (The other two unstated steps I first noted, at "sole entries are admissible" and "C_T is C_U with y inserted after p", were fixed during the audit.) | Add "(e_p(U) >= 1 because P_p ends at p)". |
| F9 | (resolved) | (H), y_U- non-critical | The strong-component detour was unnecessary. | Fixed during the audit: the note now uses the direct Hamiltonian cycle. |
| F10 | (resolved) | HYP-9460 results, status and `note:` field; note section 5 | These were placeholders while orders 14-18 were already claimed FINITE-EXACT. | Fixed during the audit. |

**No finding** on any of the specific points asked:
- the path classification counts each Hamiltonian path of T ending at p exactly once, including p adjacent to y or to p_0;
- W >= 1, W >= h, and the conclusion e_p(T) >= 2;
- the induction uses (G) for U correctly, and |U| >= 4 is exactly what keeps y_U--- outside {y_U, y_U-, y_U--};
- the base case S_4;
- the bijection's canonicity, and its use of asymmetry;
- the sole-entry rule;
- every step of the run-transition proof;
- the elimination, with index ranges, n >= 6;
- the generating function.

## 5. Reproduction (Git Bash, from the worktree root, 12 threads)

```
D=04-computation/experiments/uh_noncritical_run_20261104_audit2_S10
export PATH="/c/Users/Eliott/scoop/apps/gcc/current/bin:$PATH"
gcc -O2 -march=native -fopenmp -o $D/labelled_audit2.exe $D/labelled_audit2.c
gcc -O2 -march=native -fopenmp -o $D/classes_audit2.exe  $D/classes_audit2.c
(for n in 4 5 6 7; do $D/labelled_audit2.exe A $n; done; for n in 4 5 6 7 8 9; do $D/labelled_audit2.exe B $n; done) > $D/labelled_audit2.out   # ~4 s
$D/classes_audit2.exe 15 14 7 12 > $D/classes_audit2_n15.out     # ~1.5 min: generation to 15, parents 4-14, children rows to 15
$D/classes_audit2.exe 15 15 7 12 > $D/classes_audit2_n15top.out  # ~3 min: plus parents of order 15 in top-order mode, children rows at 16
$D/classes_audit2.exe 16 15 7 12                                  # optional, ~8 min generation: order 16 by both methods (A = B = 121,393)
python3 $D/recurrence_audit2.py $D/classes_audit2_n15.out $D/classes_audit2_n15top.out > $D/recurrence_audit2.out
```

**Top-order mode** (the last order NCHECK = NMAX). Order NMAX + 1 is not generated.
- Lemma 3 is checked in capped form, min(hc(T^p), 2) = min(e_p, 2), by a DFS.
- The bijection is certified by two facts: all children's canonical words are distinct, and their number is F(2n - 4),
  which equals the order-16 class count found independently by methods A and B in the partial run.

**Files** (all in this folder):
- `labelled_audit2.c` and `labelled_audit2.out`;
- `classes_audit2.c`; `classes_audit2_n15.out` and `classes_audit2_n15top.out` (final build);
- `classes_audit2_n13.out` and `classes_audit2_n14.out`, from an earlier build without the children rows; their other
  lines have the same content;
- `classes_audit2_n16_partial_stopped.out`: generation to order 16 by both methods, with checks on parents 4-14. It was
  stopped during the order-15 checks to stay within the compute budget; those checks were redone in top-order mode;
- `recurrence_audit2.py` and `recurrence_audit2.out`;
- this `AUDIT.md`;
- `labelled_audit2.exe` and `classes_audit2.exe`: build artifacts. Rebuild them from source; they need not be committed.

Total compute for the audit was about 19 minutes on a shared 12-thread machine.
