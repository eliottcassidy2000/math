# Independent adversarial audit of kind-pasteur-2026-11-09-S20 (THM-5109, HYP-9656, MISTAKE-750)

**Auditor:** independent subagent, 2026-11-10. All numbers below were recomputed with the audit's own programs in this
folder (`a_*.c`, `a*.py`); the packet's scripts were read but not used to produce any audit number.

## Verdict: PASS WITH CORRECTIONS

The mathematics that is actually proved is correct: the convention map to THM-4879 (A)-(B), the chain law (D), the
Pillai inequality (E1) with an accurate Ellison citation, the classification (E2) and non-chaining (E3) *for new ports
inside the attached gadget*, the tiling (F), the census, and every closure number of P1-P6 reproduce exactly.

The corrections concern priority and scope, not arithmetic:
- the dyadic cone of (C), including the example 69, 133, 261, ..., is already in Session236 and in this machine's own
  S19 (THM-5085 (3));
- several summaries drop the "new port inside the gadget" qualifier, and the unqualified reading is false;
- several summaries misstate (D)/(E) as "only through h + 3^a = 2^t ... near-misses 1, 5, 13";
- HYP-9656 (A) is literally false as written (C(r0) is a set of triples and is always infinite), is quoted as fact in
  five places, and its interpretation (C) wrongly places finite fixed-port banks inside its scope.

Five MAJOR findings, nine MINOR findings.

---

## MAJOR findings

### M1. Priority: the dyadic cone of THM-5109 (C) is KNOWN (Session236, and this machine's own S19)

- **Location.** THM-5109 status ("PROVED (new): (C)-(F)"), title and section (C); PROOF.md section 3; README "Status",
  "The idea and the answer", "Inheritance" and the line "No one had the dyadic-port lane" (README line 45).
- **Problem.** The mechanism "a root response B = 2^j turns the near-join into the pure ray h + 2^K", with the very
  same example family, was published and audited earlier the same day:
  - **Session236**, `05-knowledge/results/session236_chanpak_offsets_orbits_20261109/offsets/PROOF.md` section 2
    (commit 0352fb55fa, 2026-11-09 16:19, audited): "Clock shifts A=2^r allow more rays", with the seeds
    rows [12,29,24,16,1], port 1: H = 5, R = 32, giving H(J_m) = 2^(m+5) + 5, and rows [14,28,24,16,1], port 2:
    H = 9, R = 16, giving 2^(m+4) + 9. The first is exactly THM-5109's "69, 133, 261, ...". It also has the
    unit-coefficient lemma R = 1 iff transitive and rooted at the sink, used in (E2).
  - **S19 (this machine, THM-5085 (3))**, `kps_s19_fixed_port_tails_20261109/PROOF.md` section 2 (b) (a single
    bridge after s sinks gives h0 + 2^s R0, "Session236's ray") and section 3 ("pure dyadic rays": for x = 2 a pure
    dyadic sink ray forces a single bridge and the offset h0). s sinks followed by one bridge is literally the one-arc
    join with TT_(s+1) entered at its sink, i.e. the cone of (C) with p = m.
  - **Ingredients of the source-port families are also inherited:**
    - (H, Q) = (5, 32) for the strong 4-tournament: Session227 `multiroot/PROOF.md` section 5 (2026-11-08);
    - (H, Q) = (19, 256): Session228 `unsaturated_module/PROOF.md` and Session229 `overlap/PROOF.md`;
    - Q multiplicative over ordinal sums: Session227 multiroot section 5 (h = 3^k 5^l, q = 12^k 32^l) and Session229
      `dynamic/PROOF.md` eq. (6) ("q' = q1 q2", the "inherited law");
    - the response of a root that beats every later block is the product of the blocks' Q: Session229
      `overlap/PROOF.md` section 1, eq. (1) (U = d + sum_i alpha_i prod_(j<i) H_j prod_(j>i) Q_j), with
      alpha_i = Q_i - H_i, telescopes to U = prod Q_i.
- **Evidence.** `a14_session236_seeds.out`: the audit's counter gives (h, e, R) = (5, 0, 32) and (9, 2, 16) for the two
  Session236 seeds, and rays 69, 133, 261, 517 and 41, 73, 137, 265. `a15_scope_witnesses.out` (W2): the cone output
  of (1, 0, 4) through TT_2 is literally Session236's second seed, out-masks [14, 28, 24, 16, 1].
- **Correction.**
  - Retype (C) as "KNOWN mechanism (Session236 offsets section 2; THM-5085 (3) and S19 PROOF sections 2 (b) and 3;
    near-join HYP-9453 / THM-4865)".
  - Keep as new only:
    - the exhaustive census of the 11 strong classes of order 3-8 with Q a power of two (H = 217, 223, 599 are new
      data; 5 and 19 are known pairs);
    - the explicit product family 5^a 19^b 217^c 223^d 599^f + 2^K, an immediate corollary of known ingredients.
  - Add Session227/228/229 and Session236 to `related:` and to README "Inheritance". Delete "No one had the
    dyadic-port lane".
  - Record a MISTAKE entry, or extend MISTAKE-750: this is the same failure as MISTAKE-737/749/750, the machine's own
    canon from earlier the same day. The grep that would have found it is `dyadic` or `power of 2` in
    `05-knowledge/results/kps_s19*` and `session236*`.

### M2. Scope: (E2)/(E3) hold only for new ports inside the attached gadget; the unqualified summaries are false

- **Location.**
  - THM-5109 title: "dyadic cone steps are 1 -> 9 and one explicit e-relation family ..., and no two dyadic cone steps
    chain";
  - THM-5109 (E2) (line 110): "every dyadic cone step is one of";
  - README line 27: "no two cone steps ever chain";
  - README (E2): "The dyadic cone steps are exactly";
  - README "What this means": "Pure cone inductions are impossible. This is a theorem for transitive gadgets";
  - THM-5109 "Consequence": "Pure cone inductions ... cannot be run in tournaments";
  - HYP-9640 S20 note: "Pure cones ... cannot be iterated (THM-5109 (E))";
  - both INDEX rows.

  PROOF.md (E2)/(E3) and THM-5109 (E3) are correctly qualified.
- **Problem.** A cone step can re-land on a dyadic port that lies outside the gadget: at an old vertex of T, or at a
  new source placed above N. These steps are not in the (E2) list.
- **Evidence.**
  - Explicit witnesses (`a15_scope_witnesses.out`):
    - **W1.** T = TT_2 rooted at its sink, triple (1, 1, 1). TT_3 entered at p = 3 gives H = 5 = 1 + 2^2. The old
      source of T has (e, R) = (0, 32), a dyadic port. So "1 -> 5" is a dyadic cone step.
    - **W2.** (1, 0, 4) (TT_3 at its source), TT_2 at p = 2, gives H = 9 with (9, 2, 16) at an old vertex.
    - **W3.** (33, 0, 128), TT_2 at p = 2, gives H = 289 with (289, 64, 512) at an old vertex.
  - `a11_cone_reland_anywhere.out` takes all 181 dyadic bases of order <= 8 and 1,086 cone steps with m <= 4:
    - 5 re-landings in the gadget (the classified ones);
    - 6 re-landings at old vertices;
    - 5 outputs with Q(N) a power of two.
    The values involved are 5, 9 and 289.
  - `a12_reroot_chains.out`: even with re-rooting allowed (any vertex, or a new source above N), the search finds no
    second consecutive re-landing. It used bases of order <= 8, TT_m with m <= 4, and outputs of order <= 15. So the
    broad "no two cone steps chain" is supported by evidence but is not proved.
- **Correction.**
  - Add "with the new dyadic port at the exit inside the attached transitive gadget" to the title, to (E2) and to
    every summary.
  - State the re-rooted version as OPEN with this evidence. One-step re-landings outside the gadget exist (W1-W3);
    two-step re-landings have none known at order <= 15.
  - Replace "Pure cone inductions are impossible ... theorem" by the proved statement plus the window bound of M4.

### M3. (D)/(E) misquoted in the canon summaries ("only through h + 3^a = 2^t ... near-misses 1, 5, 13")

- **Location.**
  - HYP-9453 S20 note, lines 40-41: "A cone step lands on a dyadic port again only if h + 3^a = 2^t (THM-5109 (D)). By
    Ellison's bound this happens only at the Pillai near-misses 1, 5, 13";
  - `05-knowledge/hypotheses/INDEX.md` S20 row: "a cone step lands on a dyadic port iff h + 3^(q+j-1) = 2^t (or an
    e-variant); by Ellison only the Pillai near-misses 1, 1, 5, 13 occur";
  - `05-knowledge/results/INDEX.md` S20 row: "cone steps chain only through h + 3^a = 2^t: by Ellison the only Pillai
    near-misses are 1, 1, 5, 13";
  - the softer README line 26-27;
  - the commit message of be5ce9a33c, which cannot be changed.
- **Problem.** Each of these contradicts THM-5109 itself:
  1. Case (D)(ii) with e > 0 is not of the form h + 3^a = 2^t. The realised step 33 -> 289 has h = 33, and 33 + 3^a is
     never a power of two (it is 34, or divisible by 3).
  2. The values 5 and 13 do not occur: (E2) excludes them through THM-4810.
  3. The steps that do occur are 1 -> 9 (two base triples) and the (beta) family, realised by 33 -> 289.
- **Evidence.**
  - `a07_realised.out`: the only realised steps among the 181 dyadic bases of order <= 8, with m <= 20, are
    (1, 0, 2) -> 9, (1, 1, 1) -> 9 and (33, 4, 64) -> 289. The last is brute-forced on class code 23649, order 6,
    root 2: output (289, 64, 512).
  - `a03_census.out`: H = 13 has no dyadic port, and H = 5 has dyadic ports only at j >= 5.
- **Correction (suggested text).** "A cone step re-lands on a dyadic port at the exit inside the gadget iff
  h + 3^(q+j-1) = 2^t with p = 2q + j, or, with exit at the gadget source, 2^(p-1)h - e + 3^(p-2)2^(j+1) is a power of
  two. By Ellison's bound the realisable steps are 1 -> 9 and a family with e = 2^(p-1)h - d_p 2^(j+1),
  d_p = 1, 5, 13 (example 33 -> 289); the near-miss values 5 and 13 carry no suitable port, and no two such steps
  chain."

### M4. HYP-9656 (C) and README: finite banks and fixed-port tails are not bounded-response families; "by (A)" is a non sequitur

- **Location.**
  - HYP-9656 title (C): "finite banks and pure dyadic cones (r0 = 1) cannot reach every large N";
  - HYP-9656 (C), lines 51-53: "Finite banks have bounded root responses ... By (A), neither reaches every large N";
  - README line 103: "This covers finite banks and fixed-port tails";
  - README lines 29-30: "Bounded banks and pure cones cannot do it".
- **Problem.**
  - **Fixed-port tails are outside the closure.** They keep the root in the core, whereas the closure's moves put the
    new root in the gadget.
  - **Their port responses are unbounded.** The odd parts grow without bound (THM-5096 normal form
    R = 2^n R0 - 2dm).
  - **Finite banks are infinite.** The order-7 bank attains every odd value below 10^7 except 7 and 21 (THM-5085 (4)).
    Whether a finite bank misses a positive proportion is HYP-9480 (C), which is OPEN and independent of HYP-9656.
  - **Cones are not covered by (A) either.**
    - Cone outputs h + 2^K have non-dyadic responses, so they are not members of C(1).
    - The bases 5^a 19^b ... are reached only by a root change (a source above X), so they are not in C(1) either.
- **Evidence.**
  - `a13_misc.out` (b): from a singleton core, the word SBSBSSBSBSSSB gives port responses with odd parts 1, 3, 3, 11,
    11, 11, 87, 87, 347, ..., 5551.
  - From a 3-cycle core the odd parts reach 21935.
- **Correction.**
  - Delete the finite-bank and fixed-port-tail clauses, or replace them with: "fixed-port tails are not in the scope
    of (A) (root in the core; unbounded response odd parts); their coverage is HYP-9480 (C)".
  - For cones, replace the heuristic with a theorem that needs no conjecture. Since B >= h, a dyadic port has
    h <= 2^j, so a cone value N = h + 2^K with K >= j + 1 satisfies 2^K < N <= 2^K + 2^(K-1).
  - Hence no cone from any dyadic port of any tournament ever reaches an odd N in (3 * 2^(K-1), 2^(K+1)), which is
    half of every dyadic window: 7, 13, 15, 25-31, 49-63, ... (compare THM-5055's "values just below 2^(n-1) are never
    cones").

### M5. HYP-9656 (A) is false as written, and the conjecture is stated as fact in five places

- **Location.**
  - HYP-9656 (A): "C(r0) is finite for every r0 and every finite library";
  - the HYP-9656 title: "... closures ... are finite";
  - README line 27-28: "every self-feeding one-arc closure whose root responses have odd part at most r0 is finite";
  - README line 103: "Bounded-response families are finite";
  - results INDEX: "Bounded-response self-feeding closures are finite";
  - hypotheses INDEX: "... are finite";
  - HYP-9640 S20 note: "gives finite slices".
- **Problem.**
  - C(r0) is defined as a set of triples, and it is always infinite. Sink padding keeps the odd part, so (1, 0, 2^k)
    is in C(r0) for every k >= 1. The packet's code hides this by capping R at 2 * HMAX + 2. What the data support is
    the finiteness of the value set h(C(r0)).
  - Separately, (A) is a conjecture. HYP-9656's own status (OPEN) types it correctly, but the README (twice), both
    INDEX rows and the HYP-9640 note state it unconditionally.
- **Evidence.**
  - Definition, plus `a10_closure.c`. In exact mode the program stores triples of any R and expands only those with
    R <= HMAX.
  - In `a17_tail.out` the number of P1 triples grows with the cap (4,775 at 10^10, 6,745 at 10^13), while the value
    count stays at 197.
  - The out-of-sample stress runs are consistent with finiteness of the values (`a16_stress.out`, `a17_tail.out`; see
    the item-by-item checks below).
- **Correction.**
  - Restate (A) as: "the value set h(C(r0)) is finite (equivalently, C(r0) is finite modulo sink padding and transitive
    re-padding)".
  - In every summary write "conjecturally finite (HYP-9656 (A), OPEN; saturation evidence P1-P6)".

---

## MINOR findings

1. **Check counts and file references** (PROOF.md section 2 "Checks"; THM-5109 status "4,440 cases").
   - `check_tt.py` runs 24 bases with m <= 7, all (m, p, q): 3,360 cases. Of these, 672 have p = 1, where R' is not
     compared (`Rp = None`). So "3,360 cases with p >= 2 and all q" is wrong: 2,688 cases have p >= 2.
   - R' was compared on 2,688 + 1,080 = 3,768 cases. h' and e' were compared on 4,440.
   - PROOF cites the p = 1 rows as `out/check_tt.out`; they are in `out/check_tt_p1.out`. No script for
     `check_tt_p1.out` is in the packet.
   - **Correction:** fix the counts and the file name, and add the p = 1 script.
   - **Audit:** `a01` checks every row, p = 1 included, on 7,956 cases (39 bases of orders 1-7, m <= 8, N <= 15) with
     0 mismatches. THM-4879 (B), read verbatim through the map s = m - i, S = e, also gives 0 mismatches.

2. **Closure exactness** (HYP-9656 "Exactness"; `rclosure3/4.c` header).
   - The programs drop triples with R > RMAX = 2 * HMAX + 2. "h never decreases" does not justify this.
   - The needed fact: every library matrix has E_u >= 1 and 0 <= s_u <= h_G (a09: min E_u = 1 over all 22,243
     matrices). So on realisable triples h' = h_G h - s_u e + E_u R >= R.
   - Hence triples with R > HMAX yield only values > HMAX (gadget moves) or the same value (sink move).
   - Also, h' >= h is elementary: h' - h >= (h_G - 1)h - h_G h + h = 0. THM-4875 (C) states only the equality case.
   - In practice the cap removed nothing. `a10_closures.out` gives identical triple and value counts in exact mode and
     in cap mode.
   - Compare THM-5070 (4)'s saturation argument (S18).
   - **Correction:** add the E_u >= 1 justification and the one-line monotonicity proof.

3. **(F) typing and precision.**
   - **Typing.** (F) is an elementary reformulation of known exact decoders: Session197
     `session197/path-polynomials/ONE_CORE_BINARY_INVERSE.md` section 3 (the one-core binary normal form), THM-5096
     and THM-5085 (1)/(2).
   - **Criticality.** The "exactly critical / mass one" reading has a quantitative precursor: Session 206's
     zero-root-change mass (ln 2) X and "one free block never suffices", recorded in HYP-9480 (D)/(C).
   - **Precision.** The lattice claim "every point has exactly one predecessor" fails at R' = -c, which has none
     (a08). Predecessors of physical points need not be physical.
   - **Correction:** type (F) "PROVED (elementary; reformulation of KNOWN decoders)". Credit Session197 and Session
     206. Restrict the lattice to R' > 0, or to R' != -c.

4. **"(E3) holds even for abstract triples"** (README).
   - The written proof excludes first steps from h in {5, 13} only through realisability (THM-4810).
   - Abstract first steps exist: (5, 0, 8) -> 133, (13, e, 16) -> 2061 and (13, 0, 32) -> 2061. Their outputs need the
     one-line check that they do not chain.
   - a06: a full abstract sweep with j <= 14 (44.7 million states), all p <= 60 and every exit gives 0 classification
     violations and 0 chains. The sweep with j <= 12 and no parity condition gives the same.
   - **Correction:** add the paragraph, or write "cross-checked for j <= 17".

5. **"Equivalent to cofiniteness"** (README line 110; HYP-9640 note "Cofiniteness of the orbit is equivalent ...").
   - N(r0) -> infinity is equivalent to completeness: the full closure attains every odd n except 7 and 21. A
     cofinite closure with one further miss keeps N(r0) bounded.
   - Credit THM-5070 (4) (S18): the same 264-matrix order-<= 5 library, with no odd-part bound, already attains every
     odd v <= 2 * 10^5 except 7 and 21. This is the r0 = infinity baseline of HYP-9656 (B).

6. **Pre-registration "Already seen: No persistence through the standard near-join. 15 cases tested; none has a
   dyadic port afterwards"** is false in general.
   - P4's own 33 -> 289 is a near-join (TT_3 with p = m = 3) landing on (289, 64, 512).
   - a11 and a15 give re-landings at old vertices.
   - Leave the pre-registration text unedited, but annotate it in README.

7. **Credits for the Q facts** in PROOF section 0 and README "Inheritance".
   - "Q is multiplicative over ordinal sums, because H is" should credit Session229 `dynamic/PROOF.md` eq. (6) and
     Session227 multiroot section 5.
   - "the source port over X has response Q(X)" is Session229's U-formula (`overlap/PROOF.md` section 1) specialised
     to a root that beats every block. See M1.

8. **Wording and numbers.**
   - "exactly 11 strong classes of order <= 8": the trivial order-1 tournament has Q = 2 (a03 lists 12 including it).
     Write "of order 3-8".
   - "roughly r0^1.8": the least-squares exponent over r0 = 63-2047 is 1.76, and the endpoint exponent is 1.75 (a13).
     The out-of-sample N(4095) = 145,433 gives a local exponent of about 1.4. Write "about r0^1.75 so far, possibly
     slowing".
   - P2 outcome "all below 86,371,467": that is the maximum, so write "at most".
   - MISTAKE-750 "printed in HYP-9471's title": the counts 1, 4, 12, 44, 220, 1,660, 20,500 are in HYP-9471's
     exploratory table, not its title. a09 reproduces them per exact order.
   - The Ellison quote is in Waldschmidt section 2.2.2.
   - The `rclosure4.c` header still says rclosure3.

9. **THM-5109 "Consequence (heuristic)".** Replace it with the proved content: (E3) within its scope, and the
   cone-window bound of M4. Mark the remaining heuristic as conditional on HYP-9656 (A).

---

## Item-by-item checks (what was reproduced, and how)

| item | claim | audit result | program / output |
|---|---|---|---|
| 1 | (A)-(B) = THM-4879 (A)-(B) via s = m - i, S = e | CONFIRMED. The table and THM-4879 (B) read verbatim both match brute-force joins on all 7,956 rows (p = 1 included). | a01 |
| 1 | class lists complete | cls_1..8 are pairwise non-isomorphic (canonical forms over all n! labellings) and number 1, 1, 2, 4, 12, 56, 456, 6,880 = A000568 | a02 |
| 2 | Q multiplicative; source response = Q | CONFIRMED (60 + 60 random cases) | a04 (a), (b) |
| 2 | 11 strong classes of order 3-8 with Q = 2^i | CONFIRMED, same codes and H, Q (plus the trivial order-1 class) | a03 |
| 2 | 5 + 2^K, 25 + 2^K, 95 + 2^K, 217/223/599 + 2^K | CONFIRMED by explicit tournaments (orders up to 15) | a04 (c) |
| 2 | census 148 / 69 / 181; odd parts at H = 5, 13 | CONFIRMED, directly from class lists (not from S19 port files); identical least-j list | a03 |
| 3 | (D) (i), (ii), (ii-a) | 2-adic arguments correct. No dyadic exit with q >= p; no step outside the classification | a06 |
| 4 | (E1) solutions; Ellison citation | Exhaustive x <= 400: exactly (1,0), (2,1), (5,3), (8,5). Waldschmidt arXiv:0908.4031 section 2.2.2 (audit's own download) states the bound for x >= 12, x not in {13, 14, 16, 19, 27}, all y. Numerically the inequality fails for x >= 12 exactly at those five x in [12, 20,000]. Constants 0.85573 > 0.855 and 0.63093 < 0.631 correct. | a05 |
| 5 | (E2): B >= h use, h = 5 and h = 13 exclusion, (beta) relation, 33 -> 289 | Proof correct (within scope, M2). THM-4810 (A) does state F_s(5), F_s(13) singletons. C4: Q = 32, B_w = 5, 6, 7, 9. C13 (code 305): Q = 104, B_w = 14, 15, 18, 18, 26, so no dyadic port at H = 13. 33 -> 289 brute-forced. | a03, a07, a13 (a) |
| 6 | (E3), 9-case table, after-(alpha) values | Table entries and 4, -16, 160 recomputed. No chains in the abstract sweep (j <= 14) or the realised search (m <= 20). | a06, a07 |
| 7 | (F) letters, bijectivity, parity | Letters S, B and the post-bridge maps T_j confirmed on 142 + 72 physical tournaments. Unique preimage on 6 x 30,702 lattice points; none at R' = -c. | a08 |
| 8 | P1 | CONFIRMED: 197 values (list identical), max 1,995,807, 4,775 triples | a10 |
| 8 | P2 | CONFIRMED: 12,352 values with HMAX = 10^8 and with 10^9, max 86,371,467 | a10 |
| 8 | prefix table, P3, P5, P6 | CONFIRMED, all ten (values, N(r0)) pairs; 7,587,283 / 26,567,358 triples | a10 |
| 8 | gadget tables | Audit matrices by brute force + linear solve (5 bases, 2 checks each) equal the author's sets: 264 / 1,884 / 22,243; 0 linearity failures. Author tables' SHA-256 match out/gadget_tables.sha256. | a09 |
| 8 | pre-registration timing | origin/main reflog: P1-P4 pushed 21:04:10, outputs 21:06-21:07; P5-P6 pushed 21:13:12, output 21:13:51 | git reflog, file times |
| 9 | priority | see M1, minor 3, 5, 7 | a14 |

**Out-of-sample stress of HYP-9656 (A)-(B).** All runs are in exact mode, with the audit's own closure code.

- **Larger r0 and libraries** (`a16_stress.out`, HMAX = 5 * 10^10):
  - order <= 5, r0 = 4095: 580,505 values, max 1,038,447,593, none in (1.04 * 10^9, 5 * 10^10]. The complete prefix
    is N(4095) = 145,433, so N keeps growing (55,343 -> 145,433).
  - order <= 6, r0 = 255: 104,217 values, max 12,039,412,527, prefix 18,959. The tail is thin: 51 values in
    (10^8, 10^9] and 2 above 10^9. A run capped at 10^9 would have looked saturated while missing 2 values. "Saturation"
    therefore needs a wide margin above the largest value.
  - order <= 7, r0 = 3: 823 values, max 5,987,421, prefix 13.
- **Tail probes to HMAX = 10^13** (`a17_tail.out`, three-word-key variant `a10b_closure.c`).
  - The value counts are unchanged: P1 197, P2 12,352, P3 93,159, and order-<= 6 r0 = 255 104,217.
  - The triple counts grow with the cap (P1: 4,775 triples at 10^10, 6,745 at 10^13) while the values stay fixed.
    This illustrates M5: the set of triples is infinite, the value set apparently is not.

These results are consistent with (A), read as finiteness of the values, and with (B). They are evidence, not proof.

## Files

| file | purpose |
|---|---|
| `a_paths.c`, `a_lib.py` | own Hamiltonian-path counter (H, Q, e_v, R_v, s_v), validated against definition-level brute force; constructions |
| `a_classes.c` | certifies the class lists by canonical forms |
| `a01`-`a08` | Theorem B and the THM-4879 map; census; cones; Pillai/Ellison; abstract (D)/(E) sweep; realised steps; tiling |
| `a09_gadgets.c`, `a09_compare_gadgets.py` | independent gadget matrices and comparison with the author's tables |
| `a10_closure.c`, `a10b_closure.c`, `a10_run_closures.sh`, `a16_stress.sh`, `a17_tail.sh` | independent closures (exact mode and author-cap mode); stress and tail probes |
| `a11`, `a12`, `a15` | re-landing outside the gadget (scope, M2) |
| `a13`, `a14` | misc (H = 5/13 classes, tail responses, exponent fit); Session236 seeds (M1) |
| `run_all.sh` | rebuild and rerun everything |

Binaries and the regenerable data are not kept: `a03_triples_le8.txt` (SHA-256 139d73af...deb6, rewritten by a03 in
about 1 s), the gadget tables and the value lists.
