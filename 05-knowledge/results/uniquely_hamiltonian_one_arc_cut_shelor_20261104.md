# Every uniquely Hamiltonian tournament is a one-arc join: the classical proof (Shelor 1977), an independent check through order 16, three natural proof routes that fail, and the cut-arc distribution of the unique cycle

**Session:** opus-2026-11-04-S74. Worktree `C:/tmp/wt-s74`.

**Owner prompt.** "Start another session proving every uniquely Hamiltonian tournament is a one-arc join." Owner-supplied
inspiration:
- a "twin-prime Collatz game";
- earlier sessions on pi^(n/2)/(n/2)!;
- the tournament on 5 vertices;
- a pasted text on self-inverse tournaments, the MRB constant, the Baby Monster and Conway's M_13.

**Status.**
- **KNOWN.** The theorem itself. Shelor 1977's complete maximum-outdegree insertion construction gives every uniquely
  Hamiltonian tournament a vertex of in-degree 1, and that vertex alone is a one-arc cut. The repository already cites
  Shelor:
  - session103;
  - `operation_closed_synthesis_20261029` section II;
  - `two_entry_funnels_20261103`;
  - klein's THM-4862, but only from klein's commit at 12:13 the same day, after this session's checkpoint.

  The S73 audit pointed out that S73's HYP-9457 (b) was this classical fact (MISTAKE-705). This note records the proof
  and does not claim it.
- **Concurrent same-day work (kind-pasteur-S10, HYP-9460, pre-registered 12:28).** It announces:
  - a proof of the theorem that does not use Shelor's grammar, through local roundness of the non-critical vertices
    (THM-4868, pending its run and audit);
  - the identity #cut arcs + #non-critical vertices = n + 1;
  - as a conjecture (P7), the recurrence T(n, k) = T(n - 1, k - 1) + 2 T(n - 1, k) - T(n - 2, k - 1).

  That recurrence is the denominator of the generating function that this session's audit proved independently at
  the same time (`aud_gf.py`, 12:28; the audit's git snapshot dates from 12:17). So P7 follows from the audit's proof
  together with HYP-9460's identity (C). This covers every order, including P7's blind orders 17-18; its predicted rows
  at orders 14-18 equal the closed form of section 4 term by term.
- **PROVED.**
  - the equivalences of section 1;
  - completeness, for every n, of the extension generator of section 2. Its cyclic-interval part is klein's same-day
    block lemma (`two_branch_rearming_sigma_sweep_20261104/PROOFS.md` section 6, Lemma 5), and its |O| = 1 / |I| = 1
    case is THM-4862 (b).
- **PROVED by the independent audit** (given Shelor's grammar; machine-checked through order 20):
  - the cut-arc law k = 3 + #{b, e} - #{c};
  - the closed form of T(n, k) and its generating function (section 4, HYP-9458).
- **FINITE-EXACT.**
  - The interval-extension generation reproduces Garey's counts through order 16, and every class has vertices of
    in-degree 1 and out-degree 1. klein's Shelor-free generators already matched Garey's counts through order 14; the
    additions here are orders 15-16.
  - Three natural proof routes fail, each with explicit counterexamples.
  - The cut-arc distribution T(n, k).
- **PRE-REGISTERED** (HYP-9458): the cut-arc formulas, confirmed at orders 14-16 and then proved by the audit (section 4).
- **Inspiration sources:** checked in section 5.

Independently audited, PASS WITH CORRECTIONS, all applied (MISTAKE-706): section 8.

## 1. The theorem, its classical proof, and its equivalent forms

**Theorem (KNOWN).** Every tournament T with exactly one directed Hamiltonian cycle C has a vertex of in-degree 1, and by
converse symmetry a vertex of out-degree 1. Consequently:
1. ({v}, V - {v}) is a vertex cut crossed backward by exactly one arc, p -> v. So T is the one-arc join
   J({v}, v; T - v, p) (THM-4866 (B), session29's law). Its cycle formula gives hc(T) = 1 x e_p(T - v), so
   e_p(T - v) = 1.
2. The arc p -> v of C is a cut arc. So T has arc-strong connectivity exactly 1.

**Proof.** Orders 3 and 4 are C_3 and S_4, where vertices of in-degree 1 are visible. For n >= 5, Shelor (*On
Tournaments Having A Unique Hamiltonian Circuit*, 1977, sections 2-5) proves that every tournament with a unique
Hamiltonian cycle arises from one of order n - 1 by inserting a vertex v that beats every old vertex except its chosen
parent p, which beats v. The construction is bijective and complete, as recorded in `operation_closed_synthesis_20261029`
section II, with the five port letters a-e. That vertex v has in-degree 1. For out-degree, apply the same to the converse,
which is again uniquely Hamiltonian. QED.

**Background.**
- Douglas (Proc. London Math. Soc. 21 (1970) 716-730) characterized the class.
- Garey (JCT B 13 (1972) 266-269) counted it: F(2n - 6). Garey gives only the count, not the in-degree-1 fact.
- J. W. Moon (J. Graph Theory 6 (1982) 303-308) gave a direct proof of Garey's count.
- Moon's deletion theorem is J. W. Moon, On subtournaments of a tournament, Canad. Math. Bull. 9 (1966) 297-301.
- A local criterion (a strong tournament on n >= 6 vertices is uniquely Hamiltonian iff each strong 6-vertex
  subtournament is) appeared in a search summary dated 1994, but no reference was located. The audit's exhaustive check
  is consistent with it at orders 7-9. It is not used here.

**Equivalences (PROVED, elementary).** For a tournament with a Hamiltonian cycle C, these are equivalent:
- (a) T has a vertex set A with exactly one arc from V - A to A;
- (b) some arc of C is a cut arc;
- (c) T has arc-strong connectivity 1.

The cycle must cross back into A, so the single backward arc of a one-arc cut lies on C. Then A and V - A are
complementary intervals of C, and every chord between them points from A to V - A.

## 2. An independent check that does not use Shelor's grammar (PROVED lemma; FINITE-EXACT through order 16)

**Extension lemma.** Let T be uniquely Hamiltonian of order n >= 4.
1. By Moon (1966), T has a vertex v with T' = T - v strong.
2. Every Hamiltonian cycle of T' extends to one of T: v has in- and out-neighbours on it, so it has an in -> out
   switch. Hence T' is uniquely Hamiltonian, with cycle C'.
3. The out-neighbourhood O of v has exactly one in -> out switch along C', so O and I = V' - O are complementary cyclic
   intervals of C'.
4. hc(T' + v) = sum over w in O and u in I of P_T'(w -> u), where P counts Hamiltonian paths of T' from w to u.

So U(n) is exactly the set of extensions of U(n - 1) by a vertex whose out-neighbourhood is a proper cyclic interval,
filtered by hc = 1.

**Precursors.**
- klein's block lemma (PROOFS.md section 6, Lemma 5, commit `a92e0b49be`, 10:05) already contains the cyclic-interval
  part.
- klein's generators were already Shelor-free and matched Garey's counts through order 14.
- THM-4862 (b) is the |O| = 1 / |I| = 1 case of the cycle-count formula.
- The S73 audit generated U(10) from U(9) by the same deletion-and-insertion route.

The additions here are completeness for every n, via Moon's deletion, and the orders 15-16.

**The program.** `uh_one_arc_cut_20261104.c` generates U(n) this way, with 128-bit canonical forms and the counts P from
a subset DP of the parent. It generates by deletion and interval extension, not by Shelor's insertion, and checks each
order against Garey's F(2n - 6).

**Results through order 16** (run of HYP-9458; numbers in section 4 and the output):
- The class counts are 1, 1, 3, 8, 21, 55, 144, 377, 987, 2584, 6765, ... (Garey).
- Every class has a vertex of in-degree 1 and one of out-degree 1.
- Every class has at least two one-arc cuts and at least two cut arcs (n >= 5).

## 3. Three natural proof routes that fail (FINITE-EXACT counterexamples)

These are why the classical proof goes through Shelor's specific insertion and not a generic argument.

**(a) "A strong tournament minus one arc is Hamiltonian."** If true, it would prove the theorem at once: a cut-arc-free
T would have T - e strong and Hamiltonian for every e in C, giving a second cycle. But it is false already at order 4.
Among all labelled tournaments, the strong digraphs T - e include 24 non-Hamiltonian cases out of 72 at order 4, 720 of
4,080 at order 5, and 18,480 of 292,320 at order 6. An example at order 4:

    0 -> 3,  3 -> 1,  3 -> 2,  1 -> 0,  2 -> 0,  with the pair {1, 2} deleted.

**(b) "If T is 2-arc-strong, then T - e is Hamiltonian for every arc e."** This also suffices. It holds for every class
through order 8 (2,385 two-arc-strong classes at order 8). It fails at order 9, in exactly 4 (class, arc) pairs among
the 99,618 two-arc-strong classes, with one bad arc per class.

The classical context is Fraisse-Thomassen (Graphs Combin. 3 (1987) 239-250): every k-strong tournament has a
Hamiltonian cycle avoiding any prescribed k - 1 arcs. The four order-9 witnesses are 2-arc-strong, but each has a cut
vertex, so they are not 2-strong, as that theorem requires.

The (a) and (b) computations were run with session scripts. The audit's `aud_routes.c`, in the audit folder, reproduces
them. It also adds order 7 for (a): 608,160 non-Hamiltonian cases out of 32,981,760.

**(c) "Every non-critical vertex of a uniquely Hamiltonian tournament has in- or out-degree 1."** This would turn
Moon's deletion into an induction. It is false from order 5 on. The single-cycle extensions with |O|, |I| >= 2 number
3, 14, 52, 176, 565, 1,752, 5,301, 15,748, 46,124 at orders 5-13. These are (parent, interval) pairs, not classes. A non-critical vertex can have in- and out-degree at
least 2.

## 4. The cut arcs of the unique cycle (PROVED by the audit, given Shelor's grammar; HYP-9458 pre-registered at orders 14-16)

T(n, k) is the number of classes in U(n) whose cycle has exactly k cut arcs. The in-degree-1 and out-degree-1 vertices
give cut arcs, but the two could coincide; that k >= 2 for every n follows from the count law below. Rows for k = 2, 3, ...:

| n | T(n, 2), T(n, 3), ... |
|---:|---|
| 5 | 1, 1, 1 |
| 6 | 2, 3, 2, 1 |
| 7 | 4, 7, 6, 3, 1 |
| 8 | 8, 16, 16, 10, 4, 1 |
| 9 | 16, 36, 41, 30, 15, 5, 1 |
| 10 | 32, 80, 102, 85, 50, 21, 6, 1 |
| 11 | 64, 176, 248, 231, 155, 77, 28, 7, 1 |
| 12 | 128, 384, 592, 608, 456, 259, 112, 36, 8, 1 |
| 13 | 256, 832, 1392, 1560, 1289, 819, 406, 156, 45, 9, 1 |

**Fitted forms.**
- T(n, 2) = 2^(n-5);
- T(n, 3) = n 2^(n-7) (n >= 6);
- T(n, 4) = (n^2 + n - 8) 2^(n-10) (n >= 7);
- the top diagonals are 1, n - 4, C(n - 3, 2).

These were pre-registered at orders 14-16 before the run.

**Results.**
- P1 (Garey counts): 17,711, 46,368 and 121,393 classes.
- P2: every class has in- and out-degree-1 vertices.
- All formula predictions are CONFIRMED:

| n | T(n, 2) | T(n, 3) | T(n, 4) | top diagonals |
|---:|---:|---:|---:|---|
| 14 | 512 | 1,792 | 3,232 | 1, 10, 55 |
| 15 | 1,024 | 3,840 | 7,424 | 1, 11, 66 |
| 16 | 2,048 | 8,192 | 16,896 | 1, 12, 78 |

The full row at order 16 is 2,048, 8,192, 16,896, 23,552, 24,600, 20,256, 13,504, 7,380, 3,312, 1,210, 352, 78, 12, 1.
**The proof (supplied by the independent audit; given Shelor's grammar; HYP-9458 is now PROVED).**

*Setting.* In Shelor's grammar the parent of the new vertex x is p = pred^j(v) of the newest vertex v, where:
- j = 0 for the letters b and e;
- j = 1 for a and d;
- j = 2 for c.

*Key step.* One insertion keeps exactly those old cut arcs e whose largest one-arc cut A(e) contains p, and it adds
(p, x). Every A(e) contains the newest vertex. So:
- b and e add one cut arc;
- a and d keep the count;
- c removes one.

Hence k = 3 + #{b, e} - #{c}.

*Closed form.* With m = n - 4,

    T(n, k) = sum_(d >= 0) C(m, k - 3 + 2d) C(k - 2 + d, d),
    sum T(n, k) x^(n-4) z^k = z^2 (x + (1 - x) z) / (1 - 2x - x(1 - x) z).

At z = 1 this is Garey's 1/(1 - 3x + x^2).

*Consequences.*
- The columns are T(n, 2) = 2^(n-5), T(n, 3) = A045891(n - 4) = n 2^(n-7) and
  T(n, 4) = A158920(n - 5) = (n^2 + n - 8) 2^(n-10).
- The next diagonal is T(n, n - 4) = n(n - 4)(n - 5)/6 = A005581(n - 4). So the binomial diagonal pattern of
  HYP-9458's title stops at j = 2.
- k >= 2 for every n.

*Machine check.* The law was checked on all 9,227,464 Shelor words of orders 4-20.

*Status of the triangle.* It is not in OEIS, and there is no repository precedent. session103's strong-bridge-cut
theorem and bivariate generating function, and polarized_path_calculus FUNNELS, study related statistics on the same
grammar. klein's HYP-9459 observed related one-cycle structure. kind-pasteur-S10's HYP-9460 conjectured the equivalent
recurrence (P7) the same day (see the status list).

## 5. The owner-supplied inspiration, checked

**The twin-prime Collatz-type game.**
- The rule list ends in "the first d in a list* that divides n", and the list was not supplied, so the claimed
  equivalence ("every start reaches 1 iff there are infinitely many twin primes") cannot be checked.
- Structurally it is a FRACTRAN-style counter machine with a subtract-one escape, like the PRIMEGAME analyses of S70 and
  S71 (THM-4852, THM-4854).
- Nothing about tournaments follows from it.

**pi^(n/2)/(n/2)!.** The earlier session is session127 (`session127_gsma_norm_collapse_20261031/volume/THEORY.md`,
2026-10-31). This note's canon-only search missed it. Nothing from it is used here.

**The 5-vertex tournament and "self-inverse" (self-converse) tournaments.** Checked claims from the pasted text:

| claim | verdict |
|---|---|
| eight self-converse 5-tournaments | true (A002785) |
| n = 5 is a peak of the self-converse proportion | true only in the trivial sense of a local peak among n >= 4: 1, 0.5, 0.667, 0.214, 0.193, 0.026, 0.014, 0.0009 at orders 3-10 |
| "Regular" matrix: polynomial λ^5 - 2λ^2 - 2λ - 1, spectrum {2, ±i, golden values} | false. It is Circ(Z_5, {1, 2}) (H = 15, two Hamiltonian cycles); polynomial λ^5 - 5λ^2 - 5λ - 2 = (λ - 2)(λ^4 + 2λ^3 + 4λ^2 + 3λ + 1); eigenvalues 2, -1/2 ± 1.539i, -1/2 ± 0.363i |
| any tournament matrix | tr A = tr A^2 = 0 (no loops, no 2-cycles), so the λ^4 and λ^3 coefficients vanish, and the λ^2 coefficient is minus the number of directed 3-cycles. So λ^2(λ^3 - λ - 1) is impossible for any tournament, and the regular 5-tournament has 5 directed 3-cycles, not the 2 implied by the pasted polynomial |
| "Irregular" matrix: polynomial λ^2(λ^3 - λ - 1) | false. Scores (3, 2, 2, 2, 1), H = 15, three Hamiltonian cycles (the order-5 class Q_5 of THM-4866); polynomial λ^5 - 4λ^2 - 3λ - 3, with the pair -1/2 ± (√3/2)i |
| a pentad of PG(2, 3) with no three points collinear | false. (0:1:0), (1:1:1), (1:2:1) are collinear, and arcs in PG(2, 3) have at most 4 points |
| the MRB-constant decay law | no support. The true proportions are given in the first row |

**The inspiration that does connect.**
- Neither 5-vertex matrix is uniquely Hamiltonian. The order-5 uniquely Hamiltonian classes are the three classes with
  H = 9, 9, 11 (THM-4866).
- Converse symmetry (the "self-inverse" idea) is exactly what turns Shelor's in-degree-1 vertex into an out-degree-1
  vertex in the theorem above.

## 6. Typing

| item | type |
|---|---|
| every uniquely Hamiltonian tournament has in- and out-degree-1 vertices, hence a one-arc cut (arc-strong connectivity 1) | KNOWN (Shelor 1977; Douglas 1970; Garey 1972); already cited in the repository |
| equivalences (a)-(c) | PROVED (elementary) |
| extension lemma; generator complete for every n | PROVED (the interval part is klein's Lemma 5; THM-4862 (b) is the special case); FINITE-EXACT through order 16 (klein's generators through 14) |
| failure of routes (a), (b), (c) | FINITE-EXACT counterexamples |
| T(n, k): count law, closed form, generating function | PROVED by the audit (given Shelor's grammar); pre-registered fits confirmed at 14-16; machine-checked through order 20 |
| inspiration checks | FINITE-EXACT / CITED (A002785) |

## 7. Pre-registration record

- HYP-9458: reserved in `87166334f1`; predictions and the generator committed in `525d6c2d2e`, before the order 14-16
  run.
- Orders <= 13 were exploratory.

## 8. Audit

Independent adversarial audit: [`AUDIT.md`](../../04-computation/experiments/uh_one_arc_cut_20261104_audit_S74/AUDIT.md).
Verdict: **PASS WITH CORRECTIONS.** The theorem is correctly typed KNOWN and every proof checks.

**Reproduction.** The uniquely Hamiltonian classes were reproduced three independent ways, all giving Garey's counts
and identical class sets:
- brute force through order 9;
- extension by every out-set through order 16 (no non-interval out-set ever gives exactly one cycle);
- Shelor's grammar through order 16.

The routes, the census and the section 5 checks were also reproduced. Re-running the program gives the committed output
apart from timings.

**Proof supplied by the audit.** It proved HYP-9458 (section 4).

**Majors, applied.**
- **M1.** Credit klein's same-day block lemma, klein's Shelor-free generators, THM-4862 (b) and the S73 audit's
  deletion run.
- **M2.** Retype HYP-9458 as PROVED.

**Minors, applied.**
- The route code and Fraisse-Thomassen for (b).
- The k >= 2 gap.
- The unreferenced 1994 criterion; full Moon citations; Garey credited only for the count.
- The THM-4862 timing.
- Credit to session103, polarized_path_calculus and HYP-9459.
- session127 is the pi^(n/2) session.
- Wording and the trace argument in section 5.
- (Parent, interval) pairs in (c).
- n >= 4 in the lemma.

**Unaddressed.** Shelor's text was not accessible to the auditor (Cloudflare). The attribution rests on the
repository's records and on the grammar's computational completeness through order 16.

**Logged as MISTAKE-706.**

## 9. Reproduction

```
gcc -O2 -fopenmp -o uh 04-computation/experiments/uh_one_arc_cut_20261104.c
./uh 16
```
