# Independent audit: THM-4867, HYP-9458 and the uniquely-Hamiltonian one-arc-cut note

**Audited:**
- note `05-knowledge/results/uniquely_hamiltonian_one_arc_cut_shelor_20261104.md`;
- `01-canon/theorems/THM-4867-every-uniquely-hamiltonian-tournament-is-a-one-arc-join-shelor-independent-check-cut-arcs.md`;
- `05-knowledge/hypotheses/HYP-9458-cut-arcs-of-uniquely-hamiltonian-tournaments.md`;
- `04-computation/experiments/uh_one_arc_cut_20261104.c` and `.out`;
- the two INDEX lines added in `addeb7d11b`.

**State:** worktree `C:/tmp/wt-s74`, branch `opus-2026-11-04-S74`, head `addeb7d11b`; `origin/main` at audit time was
`4b6ad434f2` (klein's cross-citation of THM-4867, 12:13, after the S74 checkpoint).

**Auditor:** independent adversarial audit agent, 2026-11-04. All code in this folder was written from scratch and
shares nothing with `uh_one_arc_cut_20261104.c`:
- its own isomorphism test (rotation-minimum chord code along the unique cycle, not refinement + search);
- its own Hamiltonian cycle and path counts (bitmask subset DPs saturated at 2);
- its own cut-arc test (reachability in T - e, not interval scanning);
- three independent generators of U(n): brute force over chord orientations, extension by *arbitrary* out-sets, and
  Shelor's grammar as recorded in session103 section 1;
- for routes (a) and (b), its own labelled enumeration and its own refinement/individualisation canonical form for
  all tournaments through order 9.

No file outside this folder was edited.

## Verdict: PASS WITH CORRECTIONS

The mathematics checks:
- the theorem is correctly typed KNOWN, and its proof (given Shelor's construction) is correct;
- the equivalences and the extension lemma are correct;
- every number in the note, the THM, the HYP and the `.out` is reproduced by independent code.

The corrections are about credit, citation and typing. In addition, this audit **proves HYP-9458** (section A): the
number of cut arcs of a Shelor word is 3 + #{b, e} - #{c}. This gives the whole triangle T(n, k), its bivariate
generating function and every fitted formula, so HYP-9458 must be retyped from "proof OPEN" to PROVED.

No computed number changes, so the verdict is not FAIL. M1 and M2 must be fixed before THM-4867 and HYP-9458 are
treated as canon.

**Suggested ledger entry.** One MISTAKE entry, of the same kind as MISTAKE-644, -677 and -699, would cover:
- same-day repository precursors presented as this session's work (M1);
- the inaccurate "already cites Shelor" list (m4);
- the missed session127 packet (m6).

## A. A result supplied by this audit: HYP-9458 is a theorem

**Setting.** This is the grammar of session103 section 1 and operation_closed_synthesis II:
- The seed is S4 on 0, 1, 2, 3, with cycle 0 -> 1 -> 2 -> 3 -> 0, chords 0 -> 2 and 1 -> 3, and type B with
  (a0, a1) = (0, 1).
- A letter inserts x right after its parent p on the cycle. The new vertex x beats every old vertex except p.
- The **newest vertex** v is a in type A and a1 in type B; at the seed, v = 1. It has in-degree 1, and its parent is
  pred(v). In type B, a0 = pred(v) also has in-degree 1, and it has out-degree n - 2 >= 2.
- The parent is p = pred^j(v):
  - j = 0 for letters b and e ("U");
  - j = 1 for letters a and d ("F");
  - j = 2 for letter c ("D").

For a cut arc e of the cycle C, every set entered only by e is a cyclic interval of C that starts at head(e) and avoids
tail(e). These intervals are nested. Write A(e) for the largest one.

**Lemma 1 (one insertion).** Let T* = T + x with parent p. The cut arcs of T* are (p, x) together with those cut arcs
e of T for which p is in A(e). Moreover A*(p, x) = {x}, and A*(e) = A(e) ∪ {x}.

*Proof.* Let A be a one-arc cut of T*.
- **x is not in A.** Then x enters every vertex of A - {p}. Suppose A - {p} is non-empty. The single entering arc
  lies on C*, so it would have to be x -> succ*(x) = s. Since a one-arc cut is an interval of C*, this forces
  A = {s}. But p -> s also enters {s}, a contradiction. So A = {p}. In that case {p} is a one-arc cut of T entered by
  the same arc, and this arc is already covered by the third case.
- **A starts at x.** If |A| >= 2, then s is in A and p -> s enters A. So A = {x}, which is a one-arc cut because the
  only in-neighbour of x is p.
- **x is in A, and A does not start at x.** Then p is in A, and A' = A - {x} is a cyclic interval of C containing p.
  The arcs entering A in T* are exactly the arcs entering A' in T: the only in-neighbour of x is p, which lies in A',
  and the arcs out of x leave A. Conversely, every one-arc cut A' of T that contains p gives the interval A' ∪ {x}
  of C*.

QED.

**Lemma 2.** In every grammar tournament, the newest vertex v lies in A(e) for every cut arc e.

*Proof.* At the seed the cut arcs are (0, 1), (3, 0) and (2, 3), with A = {1}, {0, 1, 2} and {3, 0, 1}. Each contains
v = 1. For the induction step use Lemma 1: A*(e) = A(e) ∪ {x} and A*(p, x) = {x}. QED.

**Lemma 3 (the count).** Write k for the number of cut arcs. By Lemma 1, k(T*) = 1 + #{e : p ∈ A(e)}. Each A(e) is
an interval that starts at head(e), contains v and avoids tail(e).
- **j = 0.** Every A(e) contains p = v, so k* = k + 1.
- **j = 1.** Exactly one cut arc is lost: (pred(v), v), whose tail is p. It is a cut arc because v has in-degree 1.
  Every other A(e) starts at a head other than v, runs forward to v, and so contains pred(v). So k* = k.
- **j = 2 (letter c, type B only).** Exactly two cut arcs are lost.
  - The first is (pred(a0), a0), whose tail is p. It is a cut arc because a0 has in-degree 1.
  - The second is (a0, v). Its interval A(e) could contain pred^2(v) only by wrapping around to
    A(e) = V - {a0}. That would need out(a0) = 1, but out(a0) = n - 2.

  Every other A(e) starts at a head outside {a0, v} and runs forward to v, so it contains pred^2(v). So k* = k - 1.

The seed has k = 3. **Hence k = 3 + #{b, e} - #{c} for every word.** QED.

**Corollary (all n, given the classical bijectivity of the grammar; unconditionally through order 16, see section D).**
- **Words.** The legal words are the {U, F, D} words in which deleting the F letters leaves no two adjacent D's. In
  each type, U and F each have exactly one letter.
- **Closed form.** With m = n - 4,

      T(n, k) = sum_{d >= 0} C(m, k - 3 + 2d) C(k - 2 + d, d),

      sum_{n >= 4, k} T(n, k) x^(n-4) z^k = z^2 (x + (1 - x) z) / (1 - 2x - x(1 - x) z).

  At z = 1 this is 1/(1 - 3x + x^2), Garey's count.
- **Columns.**
  - T(n, 2) = sum C(m, 2d - 1) = 2^(n-5).
  - T(n, 3) = sum (d + 1) C(m, 2d) = (m + 4) 2^(m-3) = n 2^(n-7) for n >= 6.
  - T(n, 4) = (1/8) sum_{j odd} (j + 1)(j + 3) C(m, j) = (m^2 + 9m + 12) 2^(m-6) = (n^2 + n - 8) 2^(n-10) for n >= 7.
- **Diagonals.** T(n, n-1) = 1, T(n, n-2) = n - 4 and T(n, n-3) = C(n - 3, 2). The next diagonal is
  T(n, n-4) = n(n - 4)(n - 5)/6. So the pattern C(n - 5 + j, j) in HYP-9458's title stops at j = 2.
- **Lower bound.** k >= 2 always, because D's are separated by U's, so #c <= #(b, e) + 1. This proves the note's
  "at least two cut arcs" for every n.

**OEIS (searched 2026-11-04).**
- T(n, 3) = A045891(n - 4).
- T(n, 4) = A158920(n - 5), the binomial transform of the doubled triangular numbers.
- T(n, n-4) = A005581(n - 4).
- The triangle itself is not in OEIS: neither row orientation matched, with or without the order-4 row.

I found no repository precedent for the cut-arc statistic.
- `git grep` for the rows finds only the S74 files.
- `git grep` for "cut arc" finds only the S74 files on this topic. The other hits (THM-602, capacity_flow, session130,
  two_port_seven_core and LRC reflections) use the phrase in unrelated senses.

**Machine check.** `aud_uh` builds every Shelor word of orders 4-20, 9,227,464 tournaments in all. For each it counts
the cut arcs directly, by reachability in T - e, and compares the count with 3 + #U - #D. There are 0 mismatches. The
direct T(n, k) table equals:
- the closed form at every order 4-20;
- the generating function (`aud_gf.py`);
- S74's table at every order 5-16.

## MAJOR issues

### M1. The extension lemma and the "Shelor-free generator" have same-day repository precursors that are not credited

**What the S74 files claim.**
- THM-4867 title, "PROVED here: ... the extension lemma (Moon deletion; cyclic-interval out-neighbourhoods;
  hc(T' + v) = sum P_T'(O -> I))".
- Note status lines 23-25 and section 2.
- The INDEX line: "a Shelor-free generator ... reproduces Garey's F(2n-6) through order 16".

These present both as this session's.

**Precursors in the repository.**
- **klein, `two_branch_rearming_sigma_sweep_20261104/PROOFS.md` section 6, Lemma 5** (commit `a92e0b49be`, 10:05).
  This is before the S74 reservation `87166334f1` (11:46). It reads: "if Y + w has exactly one Hamiltonian cycle,
  b is exactly one nonempty, proper cyclic block". That is steps 2-3 of the extension lemma, with the same proof by
  insertion switches.
- **klein's generators `scripts/gen1.c` / THM-4862 (f).** They are Shelor-free, built on Lemma 5, and certified by
  Garey's count:
  - 144, 377, 987, 2584 and 6765 at orders 9-13;
  - 17,711 at order 14 (THM-4862 (h)).
- **THM-4862 (b): C(U + a_exit) = s_r(U) and C(U + a_entry) = e_r(U).** These are the |O| = 1 and |I| = 1 cases of
  step 4. The second is also exactly the corollary "hc(T) = e_p(T - v)" of the theorem.
- **The S73 audit** (`scf_crevasses_20261104_audit_S73/AUDIT.md`, lines 367-371, `4884563b32`, 11:56). It generated
  U(10) from U(9) "by Moon's strong-deletion theorem and the insertion lemma". It also checked that every U-class at
  orders 3-10 has in- and out-degree-1 vertices, and the Shelor deletion property at orders 5-9 (lines 122-132,
  386-387). S74 read this audit and cites it for MISTAKE-705, but not for these points.
- **klein has already acknowledged the overlap** from the other side (`4b6ad434f2`, 12:13): "opus's extension lemma
  coincides with the block lemma".

**Required correction.**
- Retype the extension lemma as PROVED (elementary), with credit:
  - the interval part is klein's Lemma 5 (PROOFS.md section 6);
  - the |O| = 1 and |I| = 1 cycle counts are THM-4862 (b);
  - Moon-deletion generation of U(n) was used by the S73 audit at order 10.
- State the genuine additions:
  - completeness of the Moon-deletion + interval generator for every n (klein's completeness relied on exhaustive
    heredity through order 11 plus Garey's count);
  - the extension to orders 15 and 16.
- Change "independent, Shelor-free generator through order 16" to "a third Shelor-free generation, after klein's
  (through order 14) and the S73 audit's (order 10)".
- Make this correction in the THM title and status, the note status and section 2, and the results and hypotheses
  INDEX lines.

### M2. HYP-9458 is proved (section A); "proof OPEN" and "FINITE-EXACT fits" must be retyped

**Where the old typing appears.**
- HYP-9458 status, line 4 ("OPEN (proof)").
- Its "Open" paragraph.
- Note section 4, lines 146-147 ("still FINITE-EXACT fits ... a proof would come from a transfer matrix").
- The note's typing-table row "T(n, k) table; fitted forms".
- THM-4867 section "Cut arcs (FINITE-EXACT; HYP-9458)".
- The hypotheses INDEX line ("proof OPEN").

**Required correction.**
- Retype these as PROVED, given the classical bijectivity of Shelor's grammar. Through order 16 the result is
  unconditional, because the grammar's class sets are certified by two Shelor-free generators (section D).
- Record the lemma k = 3 + #{b, e} - #{c}, the closed form and the generating function.
- Record the OEIS identifications of the columns, and that the triangle is absent from OEIS.
- Credit this audit.

The P3-P6 predictions were genuine and correct. They are now consequences of the theorem.

## MINOR issues

**m1. Routes (a) and (b) have no committed program or output.**

Problem:
- THM-4867 lists only `uh_one_arc_cut_20261104.c/.out` under verification. That program computes only route (c)'s
  numbers (`ext2`).
- The FINITE-EXACT numbers of note section 3 (a) and (b) therefore have no artifact: 72/24, 4080/720, 292320/18480;
  2,385 classes; 99,618 classes; 4 pairs.

I reproduced every one of them exactly (section D). Required:
- commit the session's code and outputs, or cite this audit's `aud_routes.c` and its outputs;
- add the classical context to (b): Fraisse-Thomassen (Graphs Combin. 3 (1987) 239-250, doi 10.1007/BF01788546)
  proved that every k-strong tournament minus any k - 1 arcs is Hamiltonian. In particular every 2-strong tournament
  minus any one arc is.
  - So the order-9 witnesses must have vertex-connectivity 1. They do, with cut vertices 5, 4, 3 and 3
    (`aud_badclasses.out`).
  - Each witness has |Aut| = 9 and exactly one bad arc. The 4 (class, arc) pairs are therefore also 4 arc orbits.
- Note, for (a), that the 24 order-4 failures are exactly S4 minus its single non-cut cycle arc. The printed example
  is that digraph, and both completions of it are S4.

**m2. A gap in the reasoning at note line 111.**

The note argues "Since in- and out-degree-1 vertices exist, k >= 2 for n >= 5". That does not follow:
- the arc into an in-degree-1 vertex v and the arc out of an out-degree-1 vertex u coincide when u = p is v's parent;
- a further argument is then needed.

The statement is true. Use section A (#c <= #(b, e) + 1), or add the missing case.

**m3. Citations in the "Background" paragraph (note lines 52-57) and in the KNOWN attributions.**
- **"A 1994 local criterion ..."** has no reference, and I could not locate one.
  - The nearest items found are Demaria-Kiihl's 1989 report on the simple quotients of tournaments with exactly one
    Hamiltonian cycle, and a Demaria-Gianella paper on local characterizations of Hamiltonian tournaments (title
    only).
  - My exhaustive check agrees with the statement. Among all strong classes of orders 7, 8 and 9 (353, 6,008 and
    178,133), hc = 1 holds if and only if every strong 6-vertex subtournament has exactly one Hamiltonian cycle,
    with no exceptions.
  - Give the reference, or type the statement as "reported; consistent with an exhaustive check through order 9".
- **"Moon gave a direct proof"** needs its citation: J. W. Moon, *The number of tournaments with a unique spanning
  cycle*, J. Graph Theory 6 (1982) 303-308. The abstract describes "another proof ... based on a recursive
  construction". klein's README already cites it.
- **The attribution "KNOWN (Shelor 1977; Douglas 1970; Garey 1972)"** appears in the THM title, the THM status and
  note line 182. It credits Garey (a count) and Douglas (a characterization) with the in-degree-1 statement. Attribute
  that statement to Shelor 1977, as recorded in session103 section 1 ("These maxima have indegree one"). Keep Douglas
  and Garey as background.
- **Moon (1966)** should get its full reference: *On subtournaments of a tournament*, Canad. Math. Bull. 9(3) (1966)
  297-301 (vertex-pancyclicity implies a vertex v with T - v strong).

**m4. "The repository already cites Shelor: ... klein's THM-4862" (note line 19; THM status lines 5-6) is inaccurate.**
- At the time of writing, THM-4862 did not mention Shelor. It cites Garey via THM-4846. Shelor entered THM-4862 only
  through klein's cross-citation `4b6ad434f2` (12:13), after the S74 checkpoint.
- Rephrase: THM-4862 (b) and (f) supply the sole-exit bijection, finite-range evidence and the corollary
  hc = e_p(T - v).

**m5. Further credit.**
- **session103's `singularity/FORBIDDEN_WORDS_ARE_STRONG_BRIDGE_CUTS.md`** proves that a Shelor tournament has a
  one-arc cut with both sides strong and of size >= 3 iff its word contains `be` or `ee`. Its
  `languages/ACTUAL_GRAMMAR_REFINEMENTS.md` gives a bivariate occurrence generating function Z(t, z). These are the
  closest prior results on one-arc cuts in Shelor's grammar. Cite them in section 4.
- **`polarized_path_calculus_20261103/FUNNELS.md` F2** classifies the rooted strong tournaments with exactly one
  Hamiltonian path from the root, with class series x^3/(1 - 3x + x^2). The corollary "e_p(T - v) = 1" is the
  converse-side form of that condition, so cite it there. The S73 audit's m11 already pointed to it.
- **klein's HYP-9459** (`a1c6c549f7` 11:54, `1df3d50e93` 11:59) observed that every one-cycle (and two-cycle)
  tournament of order 6-10 has an out-degree-1 vertex. This was concurrent with S74 (reservation 11:46, checkpoint
  12:03). klein has cross-cited THM-4867, so reciprocate.
- **THM-4867 `related`** should list the two_branch packet, session103, polarized_path_calculus and the S73 audit.

**m6. Note section 5 contains a factual error about the repository.**
- Line 158 says a search "found no session" for pi^(n/2)/(n/2)!. But
  `05-knowledge/results/session127_gsma_norm_collapse_20261031/volume/THEORY.md` (Session127, 2026-10-31,
  "Unit-ball volume, algebraic norm, and the rank-two Thue slice") is exactly such a session. It opens: "The formula
  `pi^(n/2)/((n/2)!)` is the volume of an n-dimensional Euclidean unit ball".
- Cite it, and say whether anything in it connects. Nothing in THM-4867 depends on it.

**m7. Section 5 wording and an optional strengthening.**
- **The n = 5 peak.** Among the computed orders 3-10, n = 5 is the unique interior local maximum of the self-converse
  proportion, and the maximum over n >= 4. "True only in the trivial sense" undersells this.
- **A structural reason the owner's polynomials fail.** Every tournament has tr A = tr A^2 = 0, so the coefficients
  of λ^4 and λ^3 vanish, and the coefficient of λ^2 is minus the number of cyclic triples.
  - λ^2(λ^3 - λ - 1) has λ^3-coefficient -1, so it is the characteristic polynomial of no tournament.
  - λ^5 - 2λ^2 - 2λ - 1 would need 2 cyclic triples, while a regular 5-tournament has 5.
  - The claimed spectrum {2, ±i, golden values} has sum of squares 5, not 0.
- **The "irregular" identification.** Three strong 5-classes have scores (3, 2, 2, 2, 1):
  - H = 11 with hc = 1 (uniquely Hamiltonian);
  - H = 13 with hc = 2, polynomial λ^5 - 4λ^2 - 3λ - 2;
  - H = 15 with hc = 3, polynomial λ^5 - 4λ^2 - 3λ - 3.

  The identification with Q_5 therefore rests on H = 15. I did not have the owner's pasted text, so I could not
  check the matrices against it.

**m8. Route (c) wording (note lines 104-107).**
- The numbers 3, 14, 52, ..., 46,124 (and 133,570, 383,243 and 1,091,184 at orders 14-16) count pairs (parent class,
  proper cyclic interval O) with hc = 1 and |O|, |I| >= 2. They are not isomorphism classes of tournaments or of
  (T, v). Say so.
- Give an order-5 witness: S4 + v with N+(v) = {1, 2}, which has hc = 1, T - v = S4 strong, and d+(v) = d-(v) = 2.
  Three of the four length-2 intervals work.

**m9. THM-4867 "Lemma" (line 55) omits n >= 4.** C3 has no non-critical vertex. The note's section 2 states n >= 4.

**m10. Shelor's primary text was not accessible to this audit.**
- The paperzz mirror returns a Cloudflare challenge, which I did not attempt to bypass. A DTIC PDF that a title
  search returned as a candidate gave 403, zbMATH has no record under the author name, and I found no other copy.
- As for the S73 audit, the attribution therefore rests on two records: operation_closed_synthesis II's inspection of
  the mirror, and session103's twice-audited use of the grammar.
- Computationally, the recorded grammar is complete and injective through order 16. Its class sets equal those of two
  Shelor-free generators.
- Suggest one sentence in THM-4867 recording this provenance.

**m11. Pre-registration record (checked; fine).**
- The predictions P1-P6 and the program were committed together in `525d6c2d2e` (author 11:59:19, committer
  11:59:26).
- `git diff 525d6c2d2e HEAD -- uh_one_arc_cut_20261104.c` is empty.
- The `.out` first appears in `addeb7d11b` (12:03:58). The file's mtime is 12:03:56.5, and it records `real 4m0.4s`,
  so the run started at about 11:59:56, after the commit.
- Git cannot exclude an earlier run.
- P1 and P2 are correctly labelled sanity checks.
- Cosmetic: the reproduction line builds `uh`, while the HYP says `uh_one_arc_cut 16`.

## Confirmed as correct

- **The theorem's proof, given Shelor.**
  - The inserted vertex beats all old vertices but its parent, so it has out-degree n - 2, i.e. in-degree 1.
  - The converse of a uniquely Hamiltonian tournament is uniquely Hamiltonian.
  - {v} is entered only by p -> v.
  - T = J({v}, v; T - v, p) matches THM-4866 (B)'s definition, with A => B and a -> b reversed.
  - hc = s_v({v}) e_p(T - v) = e_p(T - v), so e_p(T - v) = 1.
  - The arc p -> v is a cut arc, so λ(T) = 1.
  - Orders 3 and 4 are C3 and S4.
- **The equivalences.**
  - Arcs off C are never cut arcs.
  - (a) and (c) are equivalent for every strong digraph.
  - The single entering arc must lie on C, so both sides of the cut are intervals.
- **The extension lemma**, all four steps.
  - Moon 1966 gives v with T - v strong for n >= 4.
  - Each in-to-out switch inserts v into each cycle of T - v, injectively, so hc(T - v) <= hc(T) = 1.
  - hc(T) = 1 forces exactly one switch, so O is a cyclic interval.
  - Every Hamiltonian cycle reads v -> w ... u -> v, which gives the sum over O x I.
  - Completeness of the generator follows.
- **S74's program.**
  - Its canonical form is the minimum over all orderings consistent with an isomorphism-invariant equitable partition,
    with sound prefix pruning.
  - Its cycle back-tracking is sound when hc = 1.
  - `one_arc_cuts` is exact: every cut arc is the entering arc of an interval that starts at its head.
- **Every number**, which is reproduced below.

## D. Independent reproduction (commands and numbers)

Build and run from this folder:

```
gcc -O2 -fopenmp -o aud_uh aud_uh.c         && ./aud_uh all 9 16 16 20 > aud_uh_all_9_16_16_20.out   # 17 s, 12 threads
gcc -O2 -fopenmp -o aud_routes aud_routes.c && ./aud_routes route_a 7 > aud_routes_a_7.out          # 2.5 s
                                               ./aud_routes classes 9 > aud_routes_classes_9.out    # 1.7 s
python3 aud_gf.py > aud_gf.out; python3 aud_charpoly.py > aud_charpoly.out
python3 aud_badclasses.py > aud_badclasses.out; python3 aud_examples.py > aud_examples.out
gcc -O2 -fopenmp -o s74_uh ../uh_one_arc_cut_20261104.c && ./s74_uh 16 > s74_rerun_16.out           # re-run of S74's program
```

**U(n), three ways (`aud_uh_all_9_16_16_20.out`).**

| check | result |
|---|---|
| brute force over all chord orientations of the fixed n-cycle, orders 3-9 (2^27 = 134,217,728 at order 9) | 1, 1, 3, 8, 21, 55, 144 classes. 1, 4, 15, 48, 147, 440, 1296 cycle-labelled tournaments, i.e. n per class for n >= 4, so no rotational symmetry. |
| extension by every proper out-set O (not only intervals), orders 4-16 | Garey's 1, 3, ..., 121,393. The number of hc = 1 extensions with O not a cyclic interval is 0 at every order. |
| `ext1` / `ext2` counts, orders 4-16 | 6/0, 6/3, 16/14, 42/52, 110/176, 288/565, 754/1752, 1974/5301, 5168/15748, 13530/46124, 35422/133570, 92736/383243, 242786/1091184, identical to S74's `.out` |
| Shelor grammar, orders 4-16, each word checked | hc = 1 for every word; the classes are pairwise distinct (= F(2n - 6)). The class set **equals the extension set at every order 4-16 and the brute-force set at 4-9.** |
| in- and out-degree-1 vertices | present in every class at orders 4-16, and in every grammar word at orders 17-20. The inserted vertex always has in-degree 1. |
| T(n, k), orders 4-20 | rows identical to S74's at 5-16, including the full order-16 row. They equal the closed form and the generating function at every order, and all column and diagonal formulas hold at orders <= 20. |

**Route (a), labelled (`aud_routes_a_7.out`).** The table gives, for each order, the pairs (T, e) with T - e strong
and how many of those T - e are non-Hamiltonian.

| n | labelled T | strong T | pairs with T - e strong | non-Hamiltonian |
|---:|---:|---:|---:|---:|
| 4 | 64 | 24 | 72 | 24 |
| 5 | 1,024 | 544 | 4,080 | 720 |
| 6 | 32,768 | 22,320 | 292,320 | 18,480 |
| 7 | 2,097,152 | 1,677,488 | 32,981,760 | 608,160 |

The note's order-4 example is strong and non-Hamiltonian, and both of its completions are S4 (`aud_examples.out`).

**Classes, route (b), local criterion and self-converse counts (`aud_routes_classes_9.out`; own canonical form).**

| n | classes (A000568) | strong | self-converse (A002785) | uniquely Hamiltonian | 2-arc-strong | (class, arc) with T - e non-Hamiltonian |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 12 | 6 | 8 | 3 | 1 | 0 |
| 6 | 56 | 35 | 12 | 8 | 5 | 0 |
| 7 | 456 | 353 | 88 | 21 | 87 | 0 |
| 8 | 6,880 | 6,008 | 176 | 55 | **2,385** | 0 |
| 9 | 191,536 | 178,133 | 2,752 | 144 | **99,618** | **4** (in 4 classes) |

The four order-9 witnesses:
- out-rows: `34 36 33 7 15 472 159 287 95`; `18 20 17 7 456 31 175 303 111`; `10 12 9 448 15 31 183 311 119`;
  `10 12 9 464 39 15 183 311 119`;
- each has |Aut| = 9 and exactly one bad arc;
- each has a cut vertex, as Fraisse-Thomassen requires;
- each was re-checked by a permutation brute force in Python.

The local criterion holds with no exception at orders 7-9.

**Order 5 (`aud_charpoly.out`).**
- RT5 = Circ(Z_5, {1, 2}): λ^5 - 5λ^2 - 5λ - 2 = (λ - 2)(λ^4 + 2λ^3 + 4λ^2 + 3λ + 1). Its eigenvalues are 2,
  -1/2 ± 1.5388i and -1/2 ± 0.3633i; H = 15 and hc = 2.
- Q_5 (scores 3, 2, 2, 2, 1; H = 15; hc = 3): λ^5 - 4λ^2 - 3λ - 3 = (λ^2 + λ + 1)(λ^3 - λ^2 - 3).
- The uniquely Hamiltonian classes have H = 9, 9 and 11.
- In PG(2, 3), det[(0:1:0), (1:1:1), (1:2:1)] = 0, so the three points are collinear. Arcs have at most q + 1 = 4
  points.
- The self-converse proportions 1, 0.5, 0.667, 0.214, 0.193, 0.026, 0.014 at orders 3-9 are reproduced. Order 10,
  0.0009 = 8,784 / 9,733,056, matches OEIS and the S73 audit's census.

**Re-run of S74's own program (`s74_rerun_16.out`).** The program was compiled from the committed source and run
with `16`. Its output is identical to the committed `uh_one_arc_cut_20261104.out` line for line, apart from the
timings: 8 min 11 s here, run alongside other jobs, against 4 min 0 s.

**Same-day overlap.** I read `git log origin/main --since="2026-11-04 00:00" --name-only` in full (161 commits; saved
as `gitlog_since_20261104.txt`). The overlapping work is:
- klein's `a92e0b49be` (block lemma and Garey-certified generators), `74ebbdde57` (THM-4862), `a1c6c549f7` and
  `1df3d50e93` (HYP-9459) and `4b6ad434f2` (cross-citation of THM-4867): M1, m4, m5;
- the S73 lane (`3222ec5189`, and the audit in `4884563b32`): M1.

Two further klein commits are related but not overlapping: `67d95f2f65` (end-insertion and suffix chord rules on the
10,912 digit-1, i.e. uniquely Hamiltonian, classes of orders 8-13, which klein now types as consequences of Shelor's
grammar) and `2cfc7bba05` (the two-cycle law). No other lane touches cut arcs or uniquely Hamiltonian tournaments.

**Literature consulted.**
- zbMATH records: Douglas, PLMS (3) 21 (1970) 716-730; Garey, JCTB 13 (1972) 266-269; Moon, JGT 6 (1982) 303-308.
- The MaRDI record and abstract snippets for Moon 1982.
- Bang-Jensen-Havet-Yeo, arXiv:1905.11019, introduction, for Fraisse-Thomassen and Thomassen 1980.
- Crossref, for Moon, Canad. Math. Bull. 9(3) (1966) 297-301.
- The Springer and Crossref record of Fraisse-Thomassen, Graphs Combin. 3 (1987) 239-250 (abstract only).
- OEIS: A002785, A045891, A158920, A005581.
- Shelor 1977 was **not** accessible (m10).
