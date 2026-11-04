# Independent audit: THM-4866, HYP-9457 and the single-cycle-floor note

**Audited:**
- note `05-knowledge/results/single_cycle_floor_crevasses_lorentzian_joins_20261104.md`;
- `01-canon/theorems/THM-4866-single-cycle-floor-bottom-band-uniquely-hamiltonian-one-arc-lorentzian-joins-and-crevasses.md`;
- `05-knowledge/hypotheses/HYP-9457-single-cycle-floor-band-and-lorentzian-join-crevasses.md`;
- `04-computation/experiments/scf_crevasses_20261104.c` and its three `.out` files.

**State:** worktree `C:/tmp/wt-s73`, branch `opus-2026-11-04-S73`, head `3222ec5189`.

**Auditor:** independent adversarial audit agent, 2026-11-04. All code in this folder was written from scratch and shares
nothing with `scf_crevasses_20261104.c`: its own individualisation-refinement canonical form, its own path and cycle
DPs, |Pf| as sqrt(det) by exact Bareiss elimination rather than Pfaffian expansion, and rooted triples computed directly
from their definitions rather than through the converse. No file outside this folder was edited.

## Verdict: PASS WITH CORRECTIONS

The mathematics checks. Every proof is valid, and every number I recomputed agrees with the note, with one exception:
the interval-chain total is 1,596, not 1,601. Section "Independent reproduction" lists the reproduced numbers.

Three of the headline items, however, are typed as new or open although the repository already holds them:
- the one-arc join law and its Lorentzian form (M1);
- the one-arc cut of every uniquely Hamiltonian tournament, HYP-9457 (b) (M2);
- the uniqueness of the Hamiltonian cycle of interval chains, which is classical (m1).

A further headline statement, "the bottom band equals the uniquely Hamiltonian stratum", is false as worded (M3).

These are retyping, credit and wording corrections. None changes a computed number, so the verdict is not FAIL. All
of M1-M3 must be fixed before THM-4866 is treated as canon.

## MAJOR issues

### M1. The one-arc join law (B) and its Lorentzian reading are earlier repository results, and are not credited

The note (section 4, "The creative grouping: one-arc joins are Lorentzian inner products (PROVED; checked)"),
THM-4866 (B) and the typing table (row "one-arc join law ... PROVED") present the law as this session's result. It is
not.

**(i) Session29, 2026-10-22, independently audited twice.**

`05-knowledge/results/session29_framed_profiles_20261022/session29-profiles/ADDENDUM.md`, section 2, "The seamed
target and a two-block port formula", commit `25a4f6ba21`. Its setting: write the target as Y join X, but reverse b->a
into a->b, where a is the exit of X and b the entry of Y. It proves, for this one-arc join:

    H_target = H_X H_Y - E_Y(b) S_X(a) + (E_X(a) + Q_X(a)) (S_Y(b) + Q_Y^start(b)),

where Q_X(a) = sum over nonempty proper U containing a of E_{X[U]}(a) H(X - U).

- **The dictionary.** Put A = Y, a = b (the entry of Y), B = X and b = a (the exit of X). Then:
  - S_Y(b) + Q_Y^start(b) = R_a(A);
  - E_X(a) + Q_X(a) = R'_b(B);
  - E_Y(b) S_X(a) = e_a(A) s_b(B).

  The formula is THM-4866 (B) term for term.
- **The proof** is the same case split: the path patterns XY, YXY, XYX and YXYX.
- **Checks.** The formula was checked on 324 ten-vertex targets.
- **Audit.** `session29-wiring-audit/AUDIT.md` gave PASS.

**(ii) The Lorentzian form is also Session29's.**
- **ADDENDUM section 5**, "Mirror becomes a genuine Lorentzian adjoint". With J = diag(1, -1, 1) and port vector
  v = (H, S, R), the per-letter transfer matrices satisfy K_tau(t) = J K_t^T J, and the target count of a mirror-fixed
  word is v^T J v = H^2 - S^2 + R^2.
- **`LORENTZ_MOMENT_COROLLARIES.md`** proves signature (2, 1) for every shifted Hankel matrix. Its section 3 is
  "What the Lorentzian analogy preserves". `session29-moment-audit/AUDIT.md` ACCEPTED it, stating that "The Lorentzian
  form describes an exact marked-interface transfer and mirror-adjoint relation".
- **Iterated joins.** Session29's seamed targets (blocks joined by single seams) are iterated one-arc joins. Their 3 x 3
  transfer propagates exactly the triple (h, s, R'). So the note's "interval chains are iterated one-arc joins" is a
  special case of that construction.

**(iii) polarized_path_calculus, 2026-11-03.**

`05-knowledge/results/polarized_path_calculus_20261103/CALCULUS.md`, section C1, commit `9d16cea3bb`. It proves:

    H(X diamond Y) = F_X F_Y + h_X F_Y + F_X h_Y + (h_X - s_X)(h_Y - s_Y)
                   = u_X u_Y + v_X v_Y - w_X w_Y,

in the coordinates u = h + F, v = h - s, w = h. It calls this "a literal polarized quadratic response".

Z = X diamond Y is itself a one-arc join. The only arc leaving X toward Y* + c is a -> c, so Z = J(Y* + c, c; X, a).

**(iv) klein, the same day.**
- **THM-4862 (b)**, commit `74ebbdde57`, 09:21: C(U + a_exit) = s_r(U) and C(U + a_entry) = e_r(U). These are the
  |B| = 1 and |A| = 1 cases of hc(J) = s_a(A) e_b(B), with the same proof.
- **THM-4865 (f)** introduced the realisable-triple set {(h, e_t, B_t)}, which is the note's X, and wrote all three
  sweeps as functionals on it:
  - h + 2B is Q(x, (1, 2, 0));
  - the near-join is Q(x, (1, 2^k, 0));
  - the tower is Q(x, (A_k, D_k, C_k)).

**Required correction.**
- Retype (B) in the THM title and status, in note section 4 and in the typing table as KNOWN in the repository:
  - Session29 ADDENDUM sections 2 and 5;
  - LORENTZ_MOMENT_COROLLARIES;
  - the C1 special case of polarized_path_calculus.
- Record what is genuinely added: the statement for arbitrary rooted pairs in klein's root-response notation, and
  hc(J) = s_a e_b, an elementary generalisation of THM-4862 (b).
- Withdraw the "creative grouping" novelty framing.
- Add the Session29 and polarized_path_calculus packets to the `related` lists of THM-4866 and HYP-9457, and to both
  INDEX entries.

### M2. HYP-9457 (b) is known, and P3's join half was not a prediction; the cut is always a single vertex

**The claim.** HYP-9457 (b) (Open statements; title (b); hypotheses INDEX line 2: "OPEN ... (b) one-arc cut for every
uniquely Hamiltonian tournament") and note section 4 (lines 189-192) say that it is open whether every uniquely
Hamiltonian tournament has a cut crossed backward by exactly one arc. They add that "Douglas's characterization is the
natural place to prove it".

**The repository already records the classical complete constructor.**
- **Session103**, `05-knowledge/results/session103_words_necklaces_20261029/README.md` section 1, 2026-10-29, audited
  twice. It describes Shelor's grammar:
  - "Each new vertex x beats every old vertex except its parent p, which beats x";
  - "Type A has one maximum-outdegree vertex a. Type B has two, a0->a1. These maxima have indegree one."

  It cites the classical bijectivity and completeness of the constructor: M. R. Shelor, *On Tournaments Having A Unique
  Hamiltonian Circuit* (1977), sections 2-5.
- **`operation_closed_synthesis_20261029/THEORY.md`** section II cites the same constructor, "complete for n>=4".
- **`two_entry_funnels_20261103/README.md`** says the Garey/Shelor grammar "was already recovered".

**Consequence.** A vertex v of in-degree 1 is a cut {v} | V - v crossed backward by exactly one arc. So every uniquely
Hamiltonian tournament of order at least 3 has a one-arc cut (C_3 directly). (b) is KNOWN for all n.

**My census agrees and sharpens it.** At every order 3-9, every uniquely Hamiltonian class has a vertex of in-degree 1
and a vertex of out-degree 1. The (in-degree-1 count, out-degree-1 count) histograms are:
- n = 5: (1,1): 1, (2,2): 2;
- n = 6: 2, 1, 1, 4;
- n = 7: 4, 4, 4, 9;
- n = 8: 9, 12, 12, 22;
- n = 9: 22, 33, 33, 56.

These are Shelor's types A/B. At orders 3-9 every class also has an in-degree-1 vertex whose deletion leaves a uniquely
Hamiltonian tournament. In the full order-10 census, all 377 classes have both a vertex of in-degree 1 and a vertex of
out-degree 1.

So the "one-arc join" decomposition of the stratum is always available with |A| = 1 (or |B| = 1). There
x(A, a) = (1, 1, 1), and the form degenerates to H = h(B) + R'_b(B) - s_b(B). That is klein's one-exit law (THM-4862
(a), THM-4865 Lemma 1) on the converse side.

**klein had made the link the same morning (THM-4862 (b), 09:21, before the S73 reservation at 10:57).** The bijection
Phi sends digit-1 ready tournaments onto strong sole-exit tournaments with exactly one Hamiltonian cycle. The class
counts are Garey's F(2n - 6), verified through order 10, and THM-4862 (f) uses that count as a completeness certificate
through order 13.

**Required correction.**
- Close HYP-9457 (b) as KNOWN (Shelor 1977, as recorded in Session103 and operation_closed_synthesis II; Douglas 1970
  for the characterisation).
- State that the cut can always be taken to be a single vertex.
- Delete "This is a two-parameter numeration of the crevasse region, not a one-parameter ladder" (note line 190) and
  "so the bottom band is a Lorentzian numeration" (HYP (b)).
- Downgrade P3's join half, and P1 (see m3), from confirmed predictions to sanity checks.
- Credit THM-4862 (b) and (f) for the sole-exit / uniquely-Hamiltonian correspondence.

### M3. "The bottom band equals the uniquely Hamiltonian stratum" is false as stated

By definition the band [f(n), f_2(n)) contains only uniquely Hamiltonian classes. The stratum U(n) extends far above
the band:

| n | band end | max H in U(n) |
|---:|---:|---:|
| 7 | 41 | 51 |
| 8 | 69 | 107 |
| 9 | 115 | 227 |
| 10 | 205 | 467 |

At order 7, for example, U(7) contains 45, 47 and 51, all above the band [25, 41).

**Where the false or misleading wording occurs:**
- note title ("the bottom band of the strong spectrum is the uniquely Hamiltonian stratum");
- note status, line 25 ("The bottom band equals the uniquely Hamiltonian stratum, and all 1 + ... + 377 of its classes
  are one-arc joins"; the 610 classes are U's, not the band's);
- note line 189 ("the bottom band is the image of the Lorentzian form on pairs ... with s_a = e_b = 1"; that image is the
  value set of all of U, strictly larger than the band from order 7 on);
- HYP-9457 title (b);
- the results INDEX line 145 and the hypotheses INDEX line 2 ("bottom band = uniquely Hamiltonian stratum");
- the note's typing-table row "bottom band = uniquely Hamiltonian".

**Two further points.**
- The inclusion band ⊆ U is definitional. The FINITE-EXACT content is only the value f_2(n) = min over k >= 2 of
  f_k(n), which I confirm.
- The note's own section 3 wording, "entirely uniquely Hamiltonian", is correct.

**Required correction.** Say that the band is the bottom part of the uniquely Hamiltonian stratum, by the definition of
f_2, and that the one-arc-join (Shelor) numeration covers all of U(n), not just the band.

## MINOR issues

**m1. (A): interval chains are classically uniquely Hamiltonian; the typing and attribution need fixing.**
- **The chains are upset tournaments.** THM-4744 (2) shows that every IC(p) has score list (1, 1, 2, ..., n-2, n-2).
  Tournaments with this list are exactly Brualdi-Li's upset tournaments (R. Brualdi and Q. Li, Upsets in round robin
  tournaments, JCT B 35 (1983) 62-77).
- **Poet-Shader** (Short score certificates for upset tournaments, Electron. J. Combin. 5 (1998) #R24, Introduction)
  says: "they each have a unique hamiltonian cycle".
- **Savchenko** (arXiv:2311.13377v1, Introduction) is more nuanced than the note's citation suggests.
  - It credits [3] (Burzio-Demaria) only for the characterisation: equality in Moon's bound at l = 3 holds only for
    strong tournaments with the score list (1, 1, 2, ..., n-2, n-2).
  - The next sentence, "They all also contain exactly one hamiltonian circuit", carries no separate citation.
  - The wider uniquely Hamiltonian class is credited to Douglas [6].
  - The Crossref abstract of Burzio-Demaria says only that they characterise the family and count it.
- **Moon 1972 and Douglas 1970** were not accessible to me. The repository's reading of Moon 1972 (session27 literature
  report; `repunit_rarity_hyp9357_20261021.md`) records only the counts.
- **In-repo precursor.** THM-4655 (1e) (2026-10-14) already states that every strong transitive-plus-one tournament has
  a unique Hamiltonian cycle. That covers NT_n and every IC(p) with at most two blocks.

**Required:**
- Retype the typing-table row "interval chains: one Hamiltonian cycle, all n | PROVED" as KNOWN (Poet-Shader 1998 /
  Brualdi-Li 1983; Savchenko 2023), with the cut proof as an independent proof.
- Rephrase the Burzio-Demaria attribution as "stated in Savchenko (2023), Introduction, after citing Burzio-Demaria for
  the score-list characterisation".
- Cite THM-4655 (1e).

The cut proof works verbatim for all upset chains (parts >= 1). I checked all 2,046 labelled compositions through
order 12: each has exactly one Hamiltonian cycle.

**m2. Count error: "1,601 interval chains of order <= 17".** This appears at note line 94, THM-4866 (A) line 52 and
HYP-9457 line 34. The chains output gives 1, 1, 2, ..., 610 at orders 3-17, which is F(1), ..., F(15). Their sum is
F(17) - 1 = **1,596**. My own chain run gives 1,596, all uniquely Hamiltonian, with H equal to the product.

**m3. P1 and the hc check.**
- **P1 was not a falsifiable test.** It is implied by THM-4744 (4), which gives E(10) = 1 with the unique floor class
  IC(3, 3, 3), together with the cut proof that interval chains are uniquely Hamiltonian, proved before the test.
- **The program never checks hc.** It never tests hc(J) = s_a(A) e_b(B). The "9,000 random joins" of THM (B) "Checks"
  and the typing table check only H.
- **The 9,000 joins are not independent.** They are three runs of 3,000 with the same seed (20261104), differing only
  in the order and class moduli.

I checked both H and hc exhaustively on all 138,905 rooted class pairs with |A| + |B| <= 9, and on 4,000 independent
random labelled joins with sides of order <= 6. There were no mismatches. The |A| = |B| = 1 exception is genuine: there
H = 1 and hc = 0, while s_a e_b = 1.

**m4. (E): credit and scope.**
- **Component pairing.** That an anti-automorphism maps C_i onto C_(m+1-i) is already in THM-4743 (1b), "the
  construction".
- **H = 3 at even order.** "H = 3 never at even order" is immediate from THM-791, which gives H ≡ 1 (mod 4) for every
  tournament with a fixed-point-free involutive anti-automorphism. At even order every self-converse tournament has
  one: by THM-4818 (A), anti-automorphism cycles have length 2 (mod 4), so a suitable odd power of the anti-automorphism
  is such an involution. This is the reading of THM-4744 (3). I checked it on all 2, 12 and 176 self-converse classes of
  orders 4, 6 and 8. With no middle component, H = s^2 ≡ 1 (mod 8).
- **Scope.** THM-4866 (E), lines 107-109, must say "for self-converse T". As printed, "H = 5 is impossible at odd order"
  is false for tournaments in general: NT_4 plus a dominated vertex has order 5 and H = 5.
- **Typing.** The L-theory / Bannwart sentence (note line 246) should be typed ANALOGY.

**m5. P6's proposed mechanism is not supported, and "not on N mod 5 or mod 7" overstates flatness.**
- **The proposed mechanism fails.** Among the 2,732 triples of order <= 7, only 891 (32.6 %, about 1/3) have
  h ≡ 0 (mod 3). Small tournaments do not have H divisible by 3 disproportionately. This contradicts note line 226.
- **A generic local factor exists at every prime.** For triples equidistributed mod p, a nondegenerate bilinear form has
  P(Q = 0) / P(Q = c ≠ 0) = 1 + p/(p^3 - 1). That is 1.115, 1.040 and 1.020 at p = 3, 5, 7. The observed residue-0
  excesses on [25000, 50000) are about 1.05 mod 5 and 1.01 mod 7, the generic size.
- **The extra 3-adic bias has another source.** It comes from an over-representation of triples ≡ (0, 0, 0) (mod 3):
  238 of 2,732, or 8.7 %, against 1/27 = 3.7 %.

Replace the mechanism sentence, or type it HEURISTIC.

**m6. (D): the comet language and the coverage context.**
- **The Hardy-Littlewood sentence is only an analogy.** With a finite piece set, r(N) must vanish beyond the range, so
  "the representation counts grow, as Hardy-Littlewood predicts" (note line 211) is an analogy and should be typed so.
- **The coverage is not new spectrum information.** The repository's banks realise every odd number up to 99,999,999
  except 7 and 21 (HYP-9453 status). klein's sweeps from bases of order <= 10 reach [69, 8,079,323] (THM-4865 (e);
  commit `d0a9ad0c69`, 11:08, before the S73 checkpoint at 11:16). The note quotes only the order <= 9 reach
  [69, 702189] (lines 56 and 231).

Add both facts so that (D) reads as a statement about a family, not about the spectrum.

**m7. The localization is stated more generally than shown.** Note line 143, "the region where crevasses can arise is
exactly the bottom band of each order", is a finite observation. I confirm it through order 10: at n = 4, ..., 10 every
hole of the strong spectrum in [f(n), f(n+1)) lies in [f(n), f_2(n)). At order 10 the 13 holes in [125, 225) are all
below 205. Add "through order 10".

**m8. Terminology clash.** In `h_spectrum_multiplicities_20261021.md` (S50, line 147), "the second floor" means the
second-lowest strong count: 11, 17, 27, 49, 81, 135, 243 at orders 5-11. Here f_2(n) = min{H : hc = 2}. Rename it
(for example, "two-cycle floor") or cross-reference the other usage.

**m9. The Cayley analogy (note section 7).**
- **The direction is reversed.** Ding-Kim's abstract describes Cayley's identities as expressing "principal and
  almost-principal minors of a skew-symmetric matrix" in terms of its Pfaffians. The note has the identities recover
  Pfaffians from minors.
- **A qualifier is missing.** "The deck determines H up to the one spanning term 2hc" holds at odd order only. At even
  order the merged deck determines H (THM-4863 (H)).

The INTERPRETATION typing is honest.

**m10. Wording in section 5.** "X is closed under converse, so both sides draw from the same set" gives the wrong
reason. The y-set {x(B^op, b)} equals the x-set because B -> B^op is a bijection on tournaments. I verified that the
two sets are equal: 42, 254 and 2,732 triples at orders <= 5, 6, 7.

**m11. Further credit.** Two packets deserve a line.
- `polarized_path_calculus_20261103` classifies the strong rooted tournaments with s_root = 1, with class series
  x^3/(1 - 3x + x^2). These are the "s_a = 1" pieces of the uniquely Hamiltonian joins.
- `exclusion_7_21_20261103` re-proves the 7/21 exclusion and gives a marked core-counting refinement through H <= 21.

The deck-form remark in section 1 ("the strong floor is therefore G_min + 2") is a restatement of the ceiling. It is
correctly called equivalent, but "therefore" should go.

**m12. Pre-registration record (checked; fine).**
- **What was committed first.** `4d5099bc95` (11:06:26) commits P1-P6 together with the program. That program is
  byte-identical to the final one (`git diff 4d5099bc95 3222ec5189 -- *.c` is empty).
- **When the outputs appear.** All three `.out` files first appear in `3222ec5189` (11:16:27).
- **Disclosure and timing.** The exploratory runs (census 9, chains 17, joins 5 5, joins 6 6) are disclosed. The
  ten-minute gap is consistent with the 4.5-minute order-10 census.
- **Limits.** Git cannot exclude a run before 11:06. P2, P4, P5 and P6 were genuine tests; P1 and P3 were not (m3, M2).

## Confirmed as correct

- **(A) the cut proof.** Each cut j | j+1 is spanned by exactly one back arc. A Hamiltonian cycle therefore crosses
  each cut exactly once in each direction and uses every back arc. The forward arcs tile [0, N) with forced tails, so
  there is at most one cycle, and Camion supplies one. The Moon-Yang combination gives the single-cycle floor.
- **(B) the proof.** It is correct in every junction case:
  - A_1 cannot end at a, and B_2 cannot start at b, so every A-to-B junction exists;
  - the cycle count uses b -> a exactly once;
  - it fails only for |A| = |B| = 1.

  All the special cases check: B a point gives h - e + R; B transitive rooted at its sink gives h + 2^k R; TT_k with
  TT_m rooted at source and sink gives NT, with 1 + 2^(k+m-2).
- **(C).** All numbers are correct through order 10; see the next section.
- **(D).** All coverage, window-minimum and residue numbers are correct.
- **(E).** The proof is correct: palindromic condensation, multiplicativity, and the parity of |M|. The parity
  exclusions hold for self-converse tournaments. The decomposition was checked on every self-converse class of orders
  3-10.
- **Pfaffian.** |Pf| is an isomorphism invariant, and both |Pf| and H are odd, so the chance level is 1/2. The counts at
  orders 4, 6, 8 and 10 are reproduced independently.
- **The program.** Every reproduced output agrees with my code.

## Independent reproduction (commands and numbers)

Build and run from this folder. Each program is its own code; `audit_joins.c` and `audit_extra.c` include
`audit_census.c` with `NO_MAIN`.

```
gcc -O2 -fopenmp -o audit_census audit_census.c
./audit_census 9 ears   > audit_census_9_ears.out   # full census 3-9 + order 10 by ears (Moon + insertion), ~5 s
./audit_census 10       > audit_census_10.out       # full census 3-10 (own canonical form), 4 min 2 s
gcc -O2 -fopenmp -o audit_joins audit_joins.c && ./audit_joins > audit_joins.out    # ~9 s
gcc -O2 -fopenmp -o audit_extra audit_extra.c && ./audit_extra > audit_extra.out    # ~15 s
```

**The full order-10 census (`audit_census_10.out`, 4 min 2 s on 12 threads, own canonical form)** gives:
- 9,733,056 classes, of which 9,355,949 are strong;
- f(10) = 125, with one floor class, which has hc = 1;
- cycle ceiling: minimum slack 0, one equality class (the floor), no violations;
- min{H : hc >= 2} = 205, with f_2 = 205 attained by 6 classes and f_3 = 225 by 6;
- U(10): 377 classes, 91 distinct H, range [125, 467]. All 377 have a one-arc cut, and all 377 have both an
  in-degree-1 and an out-degree-1 vertex;
- band holes: 127 129 131 133 137 139 141 143 149 151 157 163 173 (13). These are also exactly the strong holes in
  [125, 225);
- 8,784 self-converse classes, H = H(M) s^2 in all of them, and 67 self-converse odd gaps below 200;
- |Pf| ≡ H (mod 4) in 4,868,460 / 9,733,056 classes, and (|Pf| ≡ 1) iff hc even in 4,929,248 / 9,733,056.

Every order-10 number in the note's census output is reproduced.

**Census, orders 3-9 (full), with order 10 by ears and by full census.**

| claim | note / THM | audit |
|---|---|---|
| classes / strong, n = 3..9 | 2/1, 4/1, 12/6, 56/35, 456/353, 6880/6008, 191536/178133 | same (A000568) |
| f(n), n = 3..10 | 3, 5, 9, 15, 25, 45, 75, 125 | same |
| floor classes, all with hc = 1 | 1, 1, 2, 2, 1, 5, 3, 1 | same, all hc = 1 |
| ceiling H >= f + 2(hc - 1) | no violation; equality only at the floor | same, n = 3..10 |
| f_2(n) = min over k >= 2 of f_k(n), n = 5..10 | 13, 23, 41, 69, 115, 205 | same; 205 attained by 6 classes |
| U(n): classes / distinct H / range | 3/2/[9,11], 8/4/[15,23], 21/11/[25,51], 55/18/[45,107], 144/48/[75,227], 377/91/[125,467] | same |
| U(4), U(5), U(6), U(7) | {5}, {9,11}, {15,17,19,23}, {25,...,39,45,47,51} | same |
| one-arc cut in every U-class | 610/610 | same; always a single vertex (M2) |
| band holes | 21 (n = 6); 47, 55, 67 (8); 77, 79, 83 (9); 13 values 127-173 (10) | same; at 10: 127 129 131 133 137 139 141 143 149 151 157 163 173 |
| \|Pf\| ≡ H (mod 4); (\|Pf\| ≡ 1) iff hc even | 4/4, 30/56, 3608/6880; 1/4, 34/56, 3068/6880 | same (Bareiss det) |
| self-converse odd gaps below 200, n = 9, 10 | 38, 67 | 38, 67 |
| \|Pf\| at n = 10 | 4,868,460 / 9,733,056; 4,929,248 / 9,733,056 | same |
| self-converse classes; H = H(M) s^2 | (not stated) | 2, 2, 8, 12, 88, 176, 2752, 8784 (n = 3..10); decomposition holds in all |
| interval chains of order <= 17 | "1,601" | **1,596** (F(1) + ... + F(15)), all hc = 1, H = product |

For the ceiling at order 10, the ears route used a reduction. A Hamiltonian path whose last vertex beats its first
closes to a Hamiltonian cycle, so H >= n hc. A violation or an equality at order 10 therefore needs hc <= 15 and
H <= 153. All strong 10-classes with H <= 205 arise as ears of strong 9-classes with H <= 205, by Moon's strong-deletion
theorem and the insertion lemma, and every uniquely Hamiltonian 10-class arises from a uniquely Hamiltonian 9-class.
The full order-10 census confirms the same numbers.

**Joins** (`audit_joins.out`).
- **Law.** H and hc agree with the law on all 138,905 rooted class pairs with |A| + |B| <= 9, and on 4,000 random
  labelled joins: 0 mismatches.
- **Triple sets.** There are 42, 254 and 2,732 distinct triples at orders <= 5, 6, 7, and the x-set equals the y-set at
  each order.
- **First gap above 21.** It is 47 at orders <= 5, 2,395 at <= 6, and 84,707 at <= 7.
- **Window minima (orders <= 7).** On [128, 256), ..., [32768, 65536) the minima are 34, 40, 34, 48, 94, 160, 198, 152
  and 24.
- **Mod-3 ratios.** They are 1.972, 2.082 and 1.370 for the windows [375, 750), [5000, 10000) and [25000, 50000).
- **Mod 5 and mod 7.** The relative means are 1.039 / ≈ 0.99 for mod 5 and 1.010 / ≈ 1.00 for mod 7.

**Extra** (`audit_extra.out`).
- The degree-1 histograms of M2.
- Shelor deletion property: every U-class at orders 5-9 has an in-degree-1 vertex whose deletion leaves a uniquely
  Hamiltonian tournament.
- All self-converse classes at orders 4, 6 and 8 (2, 12, 176) have H ≡ 1 (mod 4), and H ≡ 1 (mod 8) when there is no
  middle component.
- The order-10 holes in [125, 225) all lie below 205.
- The interval chains: 1,596 through order 17. The upset chains: 2,046 through order 12, all uniquely Hamiltonian.

**Literature read in this audit.**
- **Savchenko**, arXiv:2311.13377v1: Introduction and references, read in full from the PDF.
- **Poet-Shader**, EJC 5 (1998) #R24: Introduction and references, read from the PDF.
- **Burzio-Demaria**, JGT 14 (1990): Crossref abstract only.
- **Ding-Kim**, arXiv:2601.19830: abstract.
- **Not accessed:** Moon 1972, Douglas 1970 and Shelor 1977. Shelor is relied on through the repository's own audited
  record.
