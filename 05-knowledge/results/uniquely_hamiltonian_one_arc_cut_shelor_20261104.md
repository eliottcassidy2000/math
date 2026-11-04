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
  - klein's THM-4862.

  The S73 audit pointed out that S73's HYP-9457 (b) was this classical fact (MISTAKE-705). This note records the proof
  and does not claim it.
- **PROVED.**
  - the equivalences of section 1;
  - the extension lemma of section 2, which makes a generator independent of Shelor's grammar complete.
- **FINITE-EXACT.**
  - The independent generation reproduces Garey's counts through order 16, and every class has vertices of in-degree 1
    and out-degree 1.
  - Three natural proof routes fail, each with explicit counterexamples.
  - The cut-arc distribution T(n, k).
- **PRE-REGISTERED** (HYP-9458): the cut-arc formulas at orders 14-16 (section 4).
- **Inspiration sources:** checked in section 5.

Independently audited: section 8.

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
- Garey (JCT B 13 (1972) 266-269) counted it: F(2n - 6).
- Moon gave a direct proof of Garey's count.
- A 1994 local criterion states that a strong tournament on n >= 6 vertices is uniquely Hamiltonian iff each strong
  6-vertex subtournament is.

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
the 99,618 two-arc-strong classes.

**(c) "Every non-critical vertex of a uniquely Hamiltonian tournament has in- or out-degree 1."** This would turn
Moon's deletion into an induction. It is false from order 5 on. The single-cycle extensions with |O|, |I| >= 2 number
3, 14, 52, 176, 565, 1,752, 5,301, 15,748, 46,124 at orders 5-13. A non-critical vertex can have in- and out-degree at
least 2.

## 4. The cut arcs of the unique cycle (FINITE-EXACT; HYP-9458 pre-registered at orders 14-16)

T(n, k) is the number of classes in U(n) whose cycle has exactly k cut arcs. Since in- and out-degree-1 vertices exist,
k >= 2 for n >= 5. Rows for k = 2, 3, ...:

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
These closed forms are still FINITE-EXACT fits, now confirmed out of sample; a proof would come from a transfer matrix
on Shelor's port grammar (HYP-9458).

## 5. The owner-supplied inspiration, checked

**The twin-prime Collatz-type game.**
- The rule list ends in "the first d in a list* that divides n", and the list was not supplied, so the claimed
  equivalence ("every start reaches 1 iff there are infinitely many twin primes") cannot be checked.
- Structurally it is a FRACTRAN-style counter machine with a subtract-one escape, like the PRIMEGAME analyses of S70 and
  S71 (THM-4852, THM-4854).
- Nothing about tournaments follows from it.

**pi^(n/2)/(n/2)!.** A search of canon for this expression (and for the n-ball volume) found no session. A search of the
reflections timed out. Nothing is used from it.

**The 5-vertex tournament and "self-inverse" (self-converse) tournaments.** Checked claims from the pasted text:

| claim | verdict |
|---|---|
| eight self-converse 5-tournaments | true (A002785) |
| n = 5 is a peak of the self-converse proportion | true only in the trivial sense of a local peak among n >= 4: 1, 0.5, 0.667, 0.214, 0.193, 0.026, 0.014, 0.0009 at orders 3-10 |
| "Regular" matrix: polynomial λ^5 - 2λ^2 - 2λ - 1, spectrum {2, ±i, golden values} | false. It is Circ(Z_5, {1, 2}) (H = 15, two Hamiltonian cycles); polynomial λ^5 - 5λ^2 - 5λ - 2 = (λ - 2)(λ^4 + 2λ^3 + 4λ^2 + 3λ + 1); eigenvalues 2, -1/2 ± 1.539i, -1/2 ± 0.363i |
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
| extension lemma; independent generator complete | PROVED; FINITE-EXACT through order 16 (Garey's counts) |
| failure of routes (a), (b), (c) | FINITE-EXACT counterexamples |
| T(n, k) table; fitted forms | FINITE-EXACT through order 13; HYP-9458 pre-registered at 14-16 |
| inspiration checks | FINITE-EXACT / CITED (A002785) |

## 7. Pre-registration record

- HYP-9458: reserved in `87166334f1`; predictions and the generator committed in `525d6c2d2e`, before the order 14-16
  run.
- Orders <= 13 were exploratory.

## 8. Audit

(Filled after the independent audit.)

## 9. Reproduction

```
gcc -O2 -fopenmp -o uh 04-computation/experiments/uh_one_arc_cut_20261104.c
./uh 16
```
