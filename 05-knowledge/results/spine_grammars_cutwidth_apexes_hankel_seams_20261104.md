# Graded spine grammars G(k, w): Hamiltonian-path counts as an N-rational series over a cut-state automaton, exact Hankel and seam ranks, and the cutwidth censuses

**Session:** opus-2026-11-04-S79. Worktree `C:/tmp/wt-s79` (check_docs run through a `subst` drive, see housekeeping).
Canon: THM-4878. Pre-registration: `04-computation/experiments/spine_grammar_prereg_20261104.md` (commit 756414ff80,
before any census of G(0,2), G(0,3), G(1,1), G(1,2) was run). Hypothesis: HYP-9470. Independent audit: section 9.

**Owner prompt (verbatim).** "spend another session looking for recursive classifications of pieces of tournaments
structure that can compose and grow with size in order to come to the proof of the full {7,21} spectrum"

**Target.** Chan-Pak Conjecture 6.2 (TCS 1015 (2024) 114776, arXiv:2308.10214): every odd number except 7 and 21 is
the number H(T) of Hamiltonian paths of some tournament T. OPEN. THM-1370 proves that 7 and 21 never occur; Royle's
census (cited in THM-4749) verifies the conjecture below 80,557 and the strong witnesses of HYP-9029/THM-4102/4104
reach 80,265.

## 0. Summary

The session answers the prompt with one object. A **graded spine grammar** G(k, w) consists of the tournaments built
from a transitive chain v_1 < ... < v_m by reversing a set of chain arcs of **cutwidth at most w** (at most w reversed
arcs cross any cut between consecutive chain vertices) and adding **k apexes** (k <= 2) that beat arbitrary subsets of
the chain, with x_1 -> x_2. Reading the chain left to right, each vertex is a **letter** (which pending reversed arcs end
here, how many begin here, which apexes beat it), and the Hamiltonian-path count is computed by a finite automaton over
**cut states** (the fragments of the unknown path restricted to the prefix, each knowing how it will be entered and
which future vertices may not follow it). So H on G(k, w) is an N-rational series: pieces are letters, composition is
matrix multiplication, growth is word length, and the conjecture becomes a question about the range of an explicit
rational series of modest rank. The repository's main spectrum families are the smallest grades:

| grade | family (repository name) | automaton states | Hankel rank | where |
|---|---|---|---|---|
| G(1,0) | one-apex words (feedback vertex number 1) | 15 | 4 | Session 16, Session 191, HYP-9468 (binary affine state); THM-4876 (6) |
| G(0,1) | upset tournaments / Moon's interval chains | 14 | 4 | THM-4744; the all-spans-1 chain is THM-4865's tribonacci tower (also THM-4731 (6), THM-4880 note section 3) |
| G(2,0) | two-apex (two-feedback) words | 200 | 17 | Sessions 191-194; rank 17 reproduces Session 191's two-apex transfer rank; contains THM-4875's family F1 |
| G(1,1) | one apex over an interval chain | 114 | 17 | new |
| G(0,2) | cutwidth 2, contains every fa <= 2 tournament of THM-4748 | 83 | 17 | new (the bounded-span analogue is mac-mini's THM-4880, same evening) |
| G(0,3) | cutwidth 3 | 521 | 79 | new |
| G(1,2) | one apex over cutwidth 2 | 994 | 78 | new |
| G(2,1) | two apexes over an interval chain | 2544 | 78 | new |

(Hankel ranks with the standard convention H(empty) = 1, i.e. the empty cut accepting; the engine's own convention
H(empty word) = 0 adds one dimension to the apex-free grades: 5, 18, 80, audit M3. Both are exact.)

The **seam ranks** (the number of bilinear terms needed to compose two spine tournaments across a cut crossed by c
reversed arcs) are 3 for c = 1 (the Lorentzian one-arc join law h h' + R R' - e e', KNOWN), **13 for c = 2**,
**62 for c = 3** and **310 for c = 4**, the same in every ambient grammar computed. **The main finding** (section 6,
against the pre-registered prediction): the apex-free cutwidth-2 grammar G(0,2), a transitive order with reversed
arcs never more than two deep, attains **every odd number up to 8191 except 7 and 21**, with at most five reversed
arcs per value (and at most three for 217 of the 256 odd numbers below 512); its two-arc layer is THM-4748. New
PROVED items: the transfer theorem with its state bound (the method itself is KNOWN: frontier dynamic programming, DEFECT_TRANSFER D1-D2, THM-4880 (A)), the
containment dictionary, the sink-run lower bound H >= 2^r + 1, and the counting lemma that no family with a bounded
number of reversed arcs can be complete. FINITE-EXACT: the ranks, the validations and the censuses. The concurrent
bounded-**span** lane of mac-mini (THM-4880, HYP-9473, same evening, missed at first filing: section 1) reached
100,000 (span 2, misses exactly 7, 21, 47) and 5,000 (span 3, misses 7, 21) with the same kind of transfer. What the
session does **not** do: prove that any G(k, w) is onto the odd numbers other than 7 and 21; that statement, for
G(0,2) and for G(2,0), is HYP-9470 and would imply Chan-Pak.

## 1. Prior art read before the work (and what was reproduced)

- **One-apex family** (feedback vertex number 1): exact geometrically weighted binary formula, four-coordinate transfer,
  omits 13 at every order, finite-state automaton for every fixed feedback-vertex number
  (`feedback_word_hamilton_spectrum_20261020`, Session 16; Session 191 `one_apex/PROOF.md`; HYP-9468's rooted triple).
  Reproduced: G(1,0) has Hankel rank 4 and attains 560 of the 1024 odd numbers below 2048 (section 5).
- **Two-apex language**: 36-state nonnegative transfer, rational dimension 28, Hankel rank 17, witnesses for every odd
  <= 230,519 except 7, 21 (`square_sum_spectrum_overnight_20261030`, Session 191 `two_apex`), periodic-counter theorem
  and the twelve-value windows 15 * 2^k + a (Session 194 `two-apex/PROOF.md`, `WINDOW15.md`), 135-value dyadic windows
  (Session 193 `power_windows`). Reproduced: G(2,0) has Hankel rank 17 and misses exactly 7 and 21 below 2048.
- **Interval chains / upset tournaments**: Moon 1972 product formula for parts >= 2 (THM-4744 (1)); Busch-extremal =
  interval chains (Moon-Yang; THM-4744 (4), THM-4745); unique Hamiltonian cycle (THM-4866 (A), Brualdi-Li, Poet-Shader).
  Reproduced: G(0,1) attains 181 of the 500 odd numbers below 1000.
- **Two reversed arcs**: THM-4748 (five interval configurations with exact formulas; the crossing formula contains
  2^(r+2) + 2 * 3^r - 1, equal to 13 at r = 1). fa <= 2 tournaments are G(0,2) words with at most two slots.
- **One-arc joins**: the Lorentzian law H(J) = x_h y_h + x_R y_R - x_e y_e (session29 ADDENDUM; THM-4866 (B); THM-4875
  (B)); one-apex x one-apex joins attain every odd <= 2^32 except 7, 21 (THM-4875 (D), HYP-9468 CONFIRMED).
- **Towers and digits**: THM-4865 (sole-exit tower is tribonacci), THM-4873, THM-4875 (A).
- **Odd-cycle collections**: H = sum over collections of pairwise vertex-disjoint odd cycles of 2^(number of cycles)
  (THM-4876 (4), KNOWN by repository: PROP-001, THM-013/014/017), used for the sink-run bound.
- **Concurrent lanes of 2026-11-04/05** (read from `git log --since`, untruncated): Session 195's rooted grammar with
  the four-coordinate interface (h, e, R, d) closed under gadget joins across one reversed arc, sole-exit vertices and
  universal sinks (misses exactly 7, 21, 67, 173 through 20,000 with gadgets of order <= 5; two six-vertex gadgets close
  67 and 173), its endpoint-one classification and path-cover algebra; Session 192's cut sumsets, weighted-cut geometry
  and symplectic compiler; klein's digit spectra (THM-4862/4865); mac-mini's THM-4876; kind-pasteur S9/S10. The
  (h, e, R, d) interface is the c = 0 cut of G(1, *) seen from a root (section 3.4).
- **Missed at first filing, found by the audit (M1; MISTAKE-720).** (i) **mac-mini-2026-11-04-ribbon** (THM-4880,
  HYP-9473, note `ribbon_tournaments_bounded_span_transfer_20261104.md`, committed 21:40-21:48, before this session's
  first checkpoint): ribbon tournaments = a transitive order with back arcs of **span** <= w; the same profile
  (frontier) transfer with 7/27/119 states; Hankel ranks 4 (span 1) and 14 (span 2) by the same reachable/co-reachable
  pairing, an integer representation of rank 62 for span 3; "junction coupling ranks" of (u, v) -> H(u J v) with rank
  3 = the one-arc join and 2 = the near-join (the concept called a seam rank here); {1}^n = the sole-exit tower;
  strong span-2 words attain every odd <= 100,000 except 7, 21, 47 and strong span-3 words every odd <= 5,000 except
  7, 21; HYP-9473 (span-3 completeness would prove Chan-Pak). Span and cutwidth are incomparable gradings (span <= s
  implies cutwidth <= 2s - 1, cutwidth bounds no span), so the G(0,2) result is not subsumed, but an apex-free
  bounded-structure family complete to a large bound was already on the branch, with a larger range. (ii)
  `overnight_rank_return_20261029/DEFECT_TRANSFER.md` D1-D2 (Codex, 2026-10-29): the backward-arc cutwidth of a
  supplied order and the component-pairing transfer with a 2^O((k+c) log(k+c)) state bound (also
  `golden_hamilton_transfer_20261103/HANDOFF.md`, `fibonacci_parity_precision_20261103/COMPLEMENTS.md`). (iii) THM-4876
  (6) / HYP-9469 (3): the one-apex closed formula and its density one half. (iv) THM-4731 (6): the transitive
  tournament with its Hamiltonian path reversed has H = A000213(n) (the tribonacci tower).

## 2. Definitions

**Spine tournaments.** Fix an order v_1 < ... < v_m. A set D of pairs (j, i) with i < j defines the tournament with
arc v_i -> v_j for every i < j except (j, i) in D, where v_j -> v_i. Every tournament is a spine tournament for every
order of its vertices (D = the arcs against the order). The **cutwidth** of (T, order) is the maximum over t of the
number of (j, i) in D with i <= t < j.

**Letters.** Reading the chain left to right, a pending **slot** is a reversed arc whose head v_i has been placed and
whose tail has not. Slots are labelled 0, 1, ... in order of creation; the slots created at one vertex form a **head
group**. The letter at a vertex is (T, b, A): T = the set of pending slots whose tail is this vertex (at most one per
head group: two would be a double arc), b = the number of slots opened here (their head is this vertex), A = the set of
apexes that beat this vertex. A **word** of G(k, w) is a letter sequence with c_0 = 0, c_t = c_(t-1) - |T_t| + b_t <= w
and c_m = 0. The tournament T(word) has vertices v_1..v_m, x_1..x_k, chain arcs as above, x_a -> v_t iff a in A_t, and
x_1 -> x_2. G(k, w) = { T(word) }. Thus G(k, w) is the class of tournaments with a set of k vertices inducing a
transitive subtournament whose deletion leaves a spine tournament of cutwidth <= w in some order.

**Cut states.** For a Hamiltonian path P of T(word) and a prefix V_t = {v_1, ..., v_t}, the **fragments** are the
components of P restricted to V_t. A fragment F = (f_1, ..., f_r) carries a **start label**: GS if f_1 is the first
vertex of P, otherwise the pending slot (head f_1) or the apex through which P enters f_1 (its predecessor in P lies
outside V_t and beats f_1); and an **end label**: GE if f_r is the last vertex of P, otherwise OPEN with the **mask**
of pending slots whose head is f_r and of apexes beating f_r (exactly the later vertices that cannot follow f_r). The
cut state sigma_t(P) is (c_t, the head groups, the multiset of labelled fragments).

## 3. Theorems

### 3.1 The cut-state transfer (PROVED)

**Theorem A.** For each (k, w) there is a finite set Q(k, w) of cut states, nonnegative integer matrices M_ell (one per
letter ell, indexed by Q x Q), apex matrices P_1, ..., P_k, the initial vector e_0 (the empty cut) and the final vector
f (indicator of the single-fragment state with labels GS, GE) such that for every word ell_1 ... ell_m of G(k, w)

    H(T(ell_1 ... ell_m)) = e_0 M_(ell_1) ... M_(ell_m) P_1 ... P_k f.

The entry M_ell[sigma, sigma'] is the number of local choices (predecessor type, successor type) for the new vertex
that transform sigma into sigma'. In particular H restricted to G(k, w) is an N-rational series over the finite letter
alphabet, and every statistic of a prefix that a suffix can see is a linear function of the cut vector e_0 M_(ell_1)
... M_(ell_t). The number of cut states with c pending slots is at most (c + k + 2)! (c + k + 3)^(c + k) (fragments
have distinct start labels among GS, the c slots and the k apexes, so there are at most 1 + c + k of them; GE is used
at most once; the masks are disjoint subsets of the c + k pending heads). Reachable counts: G(1,0) 15, G(0,1) 14,
G(2,0) 200, G(1,1) 114, G(0,2) 83, G(0,3) 521, G(1,2) 994, G(2,1) 2544.

*Proof.* Let P be a Hamiltonian path of T(word). Placing v_(t+1), its predecessor in P is (a) nothing (v_(t+1) is the
first vertex: start label GS, allowed once), (b) the end f_r of an OPEN fragment of sigma_t with f_r -> v_(t+1), i.e.
v_(t+1) is not the tail of a slot whose head is f_r (mask disjoint from T_(t+1)), (c) a later chain vertex, necessarily
the tail of a reversed arc with head v_(t+1): one of the b_(t+1) new slots (the choice of slot is the choice of which
later vertex), or (d) an apex a in A_(t+1) not already used as a start label (an apex has one successor). Its
successor is (a') nothing (GE, once), (b') a later vertex (OPEN; the mask is the set of new slots and of apexes in
A_(t+1), the later vertices that beat v_(t+1)), or (c') the first vertex of a fragment whose start label is a slot in
T_(t+1) (v_(t+1) is that slot's tail, and the arc v_(t+1) -> head is in T(word)); (b) and (c') may not name the same
fragment (that would close a cycle). Every pair of choices produces sigma_(t+1)(P) from sigma_t(P) by merging the
named fragments through v_(t+1), terminating the slots of T_(t+1) (a fragment whose start label is a terminated slot
that was not chosen in (c') has lost its only possible predecessor, so such states are discarded) and relabelling
slots in creation order. The apex steps are the same with no new slots. Conversely a sequence of admissible choices
ending in the single fragment (GS, GE) assembles a Hamiltonian path: the fragments are joined only through arcs of
T(word) (forward chain arcs when an open end is followed by the next chain vertex, which the mask guarantees is not
reversed; reversed arcs when a tail closes the fragment its slot entered; apex arcs as recorded), and the choice
sequence is recovered from P. Hence H(T(word)) counts the accepted choice sequences, which is the matrix product. The
state bound counts labelled fragment multisets. QED

The engine `spine_grammar_census_20261104.cpp` builds Q(k, w) and the matrices on demand; `spine_cutwidth_hamilton_20261104.cpp`
is an earlier implementation with global arc identities (no apexes) used as an independent check.

### 3.2 Dictionary of the grades (PROVED, elementary)

- **G(1,0)** is exactly the one-apex class (apex x_1 last, chain transitive). Its Hankel rank 4 is the four-coordinate
  transfer of Session 16: the affine state (H, A, B) with A the binary value of the word read so far and B its last
  bit, and letter 0 adding 2(A + B) to H (HYP-9468's g(2v) = g(v) + v + (v mod 2)).
- **G(2,0)** is exactly the two-apex language with the marked arc x_1 -> x_2.
- **Connected G(0,1) words** (every internal cut crossed by one slot) are the upset tournaments: the reversed arcs form
  a path from v_m down to v_1. With all parts >= 2 they are Moon's interval chains and H is the product of 1 + 2^(p-1)
  (THM-4744 (1)). The all-parts-1 chain (letters 1, 2^(r-1), 1 in the engine's indexing: open; close-and-open; close)
  is the sole-exit tower of THM-4865: v_m beats only v_(m-1), and H = 1, 3, 5, 9, 17, 31, 57, 105, 193, 355, 653, 1201,
  2209, 4063, 7473, 13745 for r = 1..16 is the tribonacci law h_r = h_(r-1) + h_(r-2) + h_(r-3) (klein's THM-4865 (c);
  the identification of the tower with the reversed-path tournament is THM-4731 (6), H = A000213(n), and with the
  span-1 ribbon run {1}^n the ribbon note's section 3; not new, audit M1). THM-4744's parts p_i are spans (the
  interval has p_i + 1 vertices), so "all parts 1" means all spans 1.
- **fa <= 2 tournaments** (THM-4748) are the G(0,2) words with at most two slots; THM-4748's five configurations are the
  cases: disjoint (two connected words in ordered join), touching (same vertex is a tail and a head), crossing (two
  pending slots with distinct heads and distinct tails), nested, shared endpoint (one head group of size 2, or one
  letter terminating two slots).
- **One-apex joins.** THM-4875's family F1 joins two one-apex tournaments across a reversed arc between the two
  apexes; after relabelling the apexes it is a two-apex word, so F1 lies in G(2,0) (as the S78 note says), and F1
  alone has positive miss density (65 of the 512 odd numbers in [2^10, 2^11), density 0.00087 at 2^31; S78 note,
  "the bilinear family alone"). The join across a reversed arc between **chain** vertices v' -> v is a G(2,1) word
  (letters (b_i, 0) on C and (1, b'_j) on C' with one slot from v to v'), a family whose coverage has not been
  computed. The first draft inferred "G(2,1) attains every odd <= 2^32 except 7, 21" from THM-4875 (D); that is
  false (audit M2): THM-4875 (D)'s second family F2 (one-apex x arbitrary rooted tournament of order <= 7) is not a
  spine word.
- **Session 195's rooted interface** (h, e, R, d) is the cut vector of a one-apex-rooted word at an empty chain cut,
  read in the basis (paths; paths ending at the root; two-path covers starting at the root; paths of the root-deleted
  tournament); its gadget joins across one reversed arc are c = 1 seams. The c = 1 seam rank 3 below is the Lorentzian
  three-term law in that language.
- **Monotonicity of the grades.** G(k, w) is contained in G(k', w') for k <= k', w <= w'. The two gradings are
  independent axes: 13 lies in G(0,2) (THM-4748, crossing with r = 1: spine 0..4 with arcs 3 -> 0, 4 -> 1) and in
  G(2,0) (word 1, 2, 0), but in neither G(1,0) (Session 16) nor G(0,1) (section 5).

### 3.3 Ranks: minimal linear representations and seam laws (FINITE-EXACT, exact rational arithmetic)

For each grammar the reachable space R = span{ e_0 M_u } (over all words u) and the observable space O = span{ M_r P f }
(over all suffixes r) were generated by breadth-first closure with exact Gaussian elimination
(`spine_grammar_ranks_20261104.py`; the modular version `spine_grammar_ranks_mod_20261104.py`, two primes, for the
largest automata, cross-checked on the small ones). The **Hankel rank** (the rank of the pairing R x O) is the
dimension of the smallest linear representation of the series (Fliess, Schuetzenberger).

| grammar | states | reachable dim | observable dim | Hankel rank, H(empty) = 1 | engine convention H(empty word) = 0 |
|---|---|---|---|---|---|
| G(1,0) | 15 | 9 | 6 | 4 | 4 |
| G(0,1) | 14 | 8 | 6 | 4 | 5 |
| G(2,0) | 200 | 125 | 25 | 17 | 17 |
| G(1,1) | 114 | 65 | 25 | 17 | 17 |
| G(0,2) | 83 | 39 | 25 | 17 | 18 |
| G(0,3) | 521 | 211 | 120 | 79 | 80 |
| G(1,2) | 994 | 505 | 116 | 78 | 78 |
| G(2,1) | 2544 | 1549 | 116 | 78 (modular, two primes) | 78 |
| G(0,4) | 3384 | 1207 | 611 | 389 (modular, two primes) | 390 (modular) |

Convention (audit M3): the series is read on the letter alphabet (indices into the letter list of the current slot
count; the audit checked that the (T, b, A)-alphabet gives the same ranks), and the empty tournament counts one
Hamiltonian path (the empty cut is accepting), which makes the c = 0 seam the pure ordered-join product (rank 1). The
engine's own convention gives the empty word the value 0 and adds one dimension to the apex-free grades. With the
standard convention the three two-radix grades **G(0,2), G(1,1), G(2,0) all have Hankel rank 17**, and the
three-radix grades G(0,3), G(1,2), G(2,1) have 79, 78, 78 (the first draft's "+1, +1, +2" pattern was the
convention artifact). G(2,0)'s rank 17 equals the Hankel rank of the two-feedback series found by Session 191 by a different method (the
"rational dimension 28" there is the reachable space of its 36-state representation). G(1,0)'s rank 4 is Session 16's
four coordinates.

**Seam laws.** Cutting a word at a cut with c pending slots, H(l r) = sum over cut states sigma with c slots of
(e_0 M_l)[sigma] (M_r P f)[sigma]: a bilinear composition law whose minimal number of terms is the **seam rank**, the
rank of the pairing restricted to prefixes ending at a c-slot cut.

| ambient grammar | c = 0 | c = 1 | c = 2 | c = 3 |
|---|---|---|---|---|
| G(0,1) | 1 (2 with the engine convention) | 3 | | |
| G(0,2) | 1 (2) | 3 | 13 | |
| G(0,3) | 1 (2) | 3 | 13 | 62 |
| G(0,4) (3,384 states, modular) | 1 (2) | 3 | 13 | 62, and c = 4: 310 |
| G(1,0) | 4 | | | |
| G(1,1) | 4 | 13 | | |
| G(1,2) | 4 | 13 | 61 | |
| G(2,1) | 17 | 61 | | |

For k = 0 the c = 0 rank is the ordered-join product (1; 2 with the engine convention, the extra direction being the
empty suffix). The c = 1 rank 3 is the Lorentzian law (h, e, R) x (h', s', L') of session29/THM-4866 (B): KNOWN,
reproduced here as the minimal rank, and the same number as THM-4880's junction coupling rank 3 of the one-arc join
(the ribbon note's section 8 defines coupling ranks of blocks, with values 1, 2, 3, 4, 8, 9, 10, 12 for span-2 blocks
of length <= 4: the seam rank of a cut is the supremum of the coupling ranks over the junctions crossing it). The
**two-arc seam rank 13**, the **three-arc seam rank 62** and the four-arc rank 310 are new as cut invariants; 13 is
the same inside G(0,2), G(0,3) and G(0,4), so it is the rank of the two-reversed-arc composition law over all pairs
of spine tournaments of cutwidth <= 4 on the two sides (covering the crossing, nested, shared-head and shared-tail
configurations jointly); the audit confirmed 13 from below by a brute-force Hankel block (1,242 x 1,465 counts). The
number 62 also appears as the rank of THM-4880's span-3 integer representation; the coincidence is unexplained. With one apex the seam ranks rise
to 13 for c = 1 (G(1,1) and G(1,2) agree) and 61 for c = 2 (G(1,2)): the apex is a third kind of crossing. An explicit
named basis of the 13 statistics is not derived here (open task).

### 3.4 Lemmas used by the census (PROVED)

- **Monotonicity (KNOWN).** H(T) >= H(T - v): every Hamiltonian path of T - v extends to one of T by inserting v
  (Redei's insertion). Hence the prefix tournament (apexes placed) bounds every completion from below.
- **Padding.** A slot-free chain vertex beaten by every apex and followed by no chain vertex is a sink; a slot-free
  first chain vertex that beats every apex is a source. Neither changes H. (A run of slot-free vertices beaten by
  every apex that is followed by no vertex beating an apex consists of ordinal-sum cuts and is padding as well, which
  is why the engine may cut such runs; with one apex the letter A = 1 of G(2,0), beaten by x_1 and beating x_2, is also
  an ordinal-sum cut that the engine does not recognise, so some first-found two-apex witnesses are padded to the
  length cap - correct values, long words; the minimal-arc runs avoid this.)
- **Sink-run bound.** Let u_1 < ... < u_r be consecutive slot-free chain vertices beaten by every apex (letters
  (0, 0, all)), followed by a chain vertex u that beats some apex x (k >= 1), or, for k = 0, lying inside the span of a
  pending slot whose tail u comes later. Then for every odd subset S of {u_1..u_r} the cycle x -> S (increasing) -> u ->
  x (respectively u -> head -> S -> u) is an odd cycle, so T has at least 2^(r-1) odd cycles and, by the odd-cycle
  collection formula (THM-4876 (4): the empty collection and the singletons), H >= 1 + 2 * 2^(r-1) = 2^r + 1. A run
  with 2^r + 1 > X therefore cannot be continued within the target range, and a word ending in the run, or continued
  only by vertices beaten by every apex, has the count of the word without it (padding).
- **Deduplication.** Two prefixes with the same count vector have the same set of future values, and the vector
  determines the pending slot structure, so visiting a vector once loses nothing; the run length is not part of the
  key because a run cut either exceeds X or is padding with the count of an explored shorter word. A false
  deduplication (a 128-bit hash collision) could only lose values, never create them, so attainment claims are
  unaffected and the miss claims for G(1,0) and G(0,1) were reproduced by other methods (section 5, audit).
- **Closer bound.** The first future vertex u that terminates a pending slot induces, with the prefix and the
  apexes, the prefix extended by one letter (T_u, 0, A_u) with T_u nonempty and admissible; the minimum of the
  resulting counts over all such letters is therefore a lower bound for every completion (insertion lemma), and it is
  what makes the search finite without deduplication: interior patterns inside a pending span (adjacent
  transpositions, transitive modules) leave the prefix count unchanged but not the closed count.
- **Length cap.** Exact minimal-arc claims use the cap of 64 chain letters: a word with at most 4 reversed arcs and
  H <= 511 has at most 9 interior gaps r_i with 1 + 2 sum 2^(r_i - 1) <= 511 (odd-cycle collections), hence chain
  length at most 58 (audit).
- **Saturation.** With counts capped at 2^62, an entry of a state that contributes to the final count is bounded by
  that count (every configuration it counts extends to a Hamiltonian path under the same suffix), so a recorded value
  <= X is exact and a saturated final sum is correctly rejected (needed: a 64-letter word of G(2,0) overflowed in an
  earlier version and produced a false 7).

## 4. Validation of the engine

- Exhaustive: the global-identity transfer agrees with the subset dynamic programme on all 2^(C(n,2)) spine
  tournaments for n <= 7 (8 + 64 + 1,024 + 32,768 + 2,097,152 tournaments; cutwidth up to 12); random orders 9, 11,
  13: 295 + 73 + 23 tournaments, no mismatch.
- The grammar engine reproduces three known spectra exactly: G(0,1) attains 181 of 500 below 1000 with first misses
  7, 11, 13, 19, 21, 23 (per-word engine and automaton engine agree); G(1,0) attains 560 of 1024 below 2048 with first
  misses 7, 13, 21, 29, 35, 41, 43 (one-apex automaton, 12.7 million nodes to 2^21, agrees); G(2,0) misses exactly 7
  and 21 below 2048 (16,028,912 words; the two-apex witness file).
- 480 random words of G(0,2), G(0,3), G(1,1), G(1,2), G(2,1), G(2,2) (chain length <= 11, order <= 15) decoded into
  explicit tournaments agree with an independent pure-Python subset dynamic programme
  (`spine_grammar_validate_20261104.py`).
- The Hankel rank 17 of G(2,0) and the rank-3 one-arc seam are independent reproductions of Session 191 and of the
  Lorentzian law.

## 5. The one-apex class attains half of the odd numbers (FINITE-EXACT)

Enumerating the one-apex automaton (H, A, B) over all words of length <= 24 (every one-apex value <= 2^21 is attained
by a word of length <= 24, since a zero at position j contributes at least 2^(j-3)): 521,915 of the 1,048,576 odd
numbers below 2^21 are one-apex values; the misses per dyadic scale [2^j, 2^(j+1)) are 1, 2, 7, 14, 29, 57, 117, 236,
495, 998, 2035, 4077, 8180, 16382, 32923, 65646, 131959, 263502 for j = 3..20, i.e. the miss fraction tends to 1/2.
The closed formula and "about half of the odd numbers" (50,006 of the 100,000 odd numbers below 200,000, first misses
7, 13, 21, 29, 35, 41, 43, 49) are THM-4876 (6) and HYP-9469 (3) (mac-mini, the same afternoon; missed at first
filing, audit M1); the extension to 2^21 and the residue profile are the additions here: flat modulo 16 and 5
(0.501-0.504), biased modulo 7 (0.455, 0.488, 0.551, 0.581, 0.533, 0.466, 0.443 for the residues 0..6, standard
error 0.0013) and mildly modulo 3 (0.501, 0.490, 0.516) (the audit corrected the first draft's "no bias"). The class
is not sparse-missing (as 13 alone might suggest): a second radix (apex or slot) is needed for half of all odd
numbers, which is why the two-radix grades carry the conjecture.

## 6. Censuses (pre-registered; all chain lengths <= 64; values <= X; products of connected words added)

**The surprise.** The apex-free cutwidth-2 grammar G(0,2) attains **every odd number <= 8191 except 7 and 21** (4094
of 4096 with the closer-bound engine, 197,660,299 search nodes, at most five reversed arcs per value; 1022 of 1024
below 2048; 254 of 256 below 512 with the first engine, 289,792,317 nodes, 43,254,577 complete connected words; 83 cut
states, Hankel rank 18). So do G(1,1) and G(0,3) below 512 (254 of 256 each). Prediction P1 (between 10 and 60 misses below 300) is **REFUTED** in the direction
nobody in the repository had considered: the sparsest grade, a transitive order with reversed arcs never more than two
deep, is complete as far as it was computed. The decoded witnesses below 302 use at most 7 reversed arcs (24 values
with 2 arcs, 41 with 3, 46 with 4, 32 with 5, 4 with 6, 1 with 7; chain length <= 12), e.g. 13 = chain 0..5 with
3 -> 0, 5 -> 1, 5 -> 4; 35 = 0..6 with 2 -> 0, 6 -> 2, 6 -> 3; 67 = 0..6 with 2 -> 0, 3 -> 1, 4 -> 2, 6 -> 3, 6 -> 4
(Session 195's two hard numbers 67 and 173 are G(0,2) values).

**Minimal numbers of reversed arcs (exact; the no-deduplication run to 511, 163,451,097 nodes).** Every odd N <= 511
other than 7 and 21 is the count of a cutwidth-2 spine tournament with at most **five** reversed arcs, and the layers
are: 1 value with no arc (1); 8 with one arc (3, 5, 9, 17, 33, 65, 129, 257 = 1 + 2^k, THM-284); 75 with two arcs
(the THM-4748 layer: 11, 13, 15, 19, 23, 25, 27, 29, 37, 39, 41, 45, 47, 49, 51, 53, 59, 61, 69, 73, 77, 79, 81, 85,
89, 93, 95, 99, 101, 105, 109, 117, 133, 137, 141, 149, 153, 155, 159, 161, 165, 169, 177, 191, 193, 195, 197, 209,
225, 229, 241, 249, 261, 265, 281, 285, 289, 293, 297, 305, 319, 325, 329, 345, 353, 369, 383, 387, 389, 409, 419, 449,
453, 465, 489); 133 with three arcs (31, 35, 55, 57, 63, 71, 75, 83, 87, 91, 97, 103, 107, 111, 115, 119, 121, 123,
125, 131, 135, 143, 145, 147, 151, 157, 163, 167, 171, 175, 179, 181, 183, 185, 187, 199, ... , 509); 36 with four
arcs (43, 113, 127, 139, 173, 189, 211, 235, 253, 259, 271, 273, 275, 291, 309, 317, 321, 331, 339, 341, 347, 357,
359, 373, 379, 393, 397, 413, 433, 437, 443, 463, 475, 491, 501, 511); and one with five (67: chain 0..6 with 2 -> 0,
3 -> 1, 4 -> 2, 6 -> 3, 6 -> 4). The 133 three-arc witnesses fall into 15 pairwise-configuration types (each pair of
arcs disjoint / touching / crossing / nested / shared head / shared tail); the three most common are two crossings
plus a touching (26 values, e.g. 4 -> 0, 5 -> 1, 6 -> 4), crossing + shared tail + touching (22, e.g. 3 -> 0, 5 -> 1,
5 -> 3) and crossing + nested + touching (20, e.g. 4 -> 0, 6 -> 4, 7 -> 3). A closed-form census of the three-arc
cutwidth-2 layer in the style of THM-4748 is therefore the natural next theorem: with the two-arc layer it covers 217
of the 256 odd numbers below 512.

**The smallest two-radix family: one apex plus at most one reversed chain arc.** In the G(1,1) census to 2047 every
value is attained with at most one reversed chain arc: the 560 one-apex values and the 462 remaining odd numbers
below 2048 other than 7 and 21, each of which is the count of a one-apex tournament (chain transitive, apex x beating
the chain vertices with bit 1) with exactly one chain arc v_j -> v_i reversed. So the one-apex formula with a single
reversed arc is a complete candidate family in range, and the closed form for it (the one-apex cycles, the cycles
through the arc, and the interleaved pairs that bring the factor 3) is the most economical target for a numeration
proof. G(1,1) itself attains every odd number <= 2047 except 7 and 21 (1022 of 1024; 498,938,267 nodes).

| grammar | X | attained (with products) | misses besides 7, 21 | search nodes |
|---|---|---|---|---|
| G(0,1) | 2047 | 358 of 1024 | 664 | 49,184 |
| G(1,0) | 2047 | 763 of 1024 | 259 | 29,841 |
| G(2,0) | 2047 | 1022 of 1024 | 0 | 20,924,275 |
| G(0,2) | 511 | 254 of 256 | 0 | 289,792,317 (first engine); 163,451,097 without deduplication (closer bound) |
| G(1,1) | 511 | 254 of 256 | 0 | 282,279,675 (first engine) |
| **G(0,2)** | **2047** | **1022 of 1024** | **0** | 14,796,241 (closer bound; at most 5 reversed arcs per value: 1, 10, 127, 542, 328, 14 values with 0..5 arcs, upper bounds on the minima from the deduplicated run) |
| **G(0,2)** | **8191** | **4094 of 4096** | **0** | 197,660,299 (at most 5 reversed arcs per value: 1, 12, 208, 1522, 2209, 142 values with 0..5 arcs) |
| G(0,3) | 511 | 254 of 256 | 0 | 101,209,735 (521 cut states; at most 4 reversed arcs per value: 1, 8, 48, 174, 23) |
| **G(1,1)** | **2047** | **1022 of 1024** | **0** | 498,938,267 (every value with the apex and **at most one** reversed chain arc: 560 one-apex values plus 462 with one arc) |
| G(1,2) | 511 | still running at close-out | | |
| G(0,3), G(1,2) | 2047 | still running at close-out | | |
| G(0,2) | 32767 | still running at close-out | | |

(The four unfinished runs write `census_k<k>_w<w>_X<X>.txt` in the session's temporary directory
`C:/Users/Eliott/AppData/Local/Temp/s79/`; each prints its summary and an `ATTAINED_DIRECT:` line on completion and
is reproduced by `grammar k w X 64 0 301` with the committed engine.)

Outcomes of the pre-registered predictions (in time order): P1 first clause REFUTED by the 511 run (21:33), second
clause REFUTED by the 8191 run (G(0,2) misses nothing below 8192, so the miss fraction in [2^11, 2^12) is 0, not
>= 3 percent); P2 first clause CONFIRMED by the G(0,3) 511 run, second clause (a miss below 2000) NOT TESTED at
close-out (the G(0,3) 2047 run was still going; given G(0,2) it is expected to fail); P3: the 13 clause CONFIRMED at
511, the 10-percent clause CONFIRMED by the G(1,1) 2047 run (22:57; the first draft said CONFIRMED before that run
existed, audit minor 4); P5 CONFIRMED (13 is in G(0,2) and G(2,0), not in G(1,0) or G(0,1)); P4 NOT TESTED at
close-out (both G(1,2) runs still going). The minimal-grade atlas
(`spine_grammar_atlas_20261104.py`) over G(1,0), G(0,1), G(2,0) at 2047: G(0,1) misses 86 odd numbers below 300
besides 7, 21 (first 11, 13, 19, 23, 35), G(1,0) misses 45 (first 13, 29, 35, 41, 43, 49, 59, 67, 71, 77); every odd
number below 300 other than 7, 21 is attained by G(2,0) and by G(0,2); the values attained by G(2,0) only among
{G(1,0), G(0,1), G(2,0)} below 120 are 13, 35, 41, 43, 59, 67, 71, 77, 83, 109, 113.

## 7. What this gives for Chan-Pak

Chan-Pak follows from the statement "some G(k, w) is onto the odd numbers other than 7 and 21", and the data make
G(2,1) (one-apex joins, 2^32), G(2,0) (230,519) and now **G(0,2)** (511 at filing, with a few reversed arcs per value)
the candidates. G(0,2) is the most structured of the three: no apex, every vertex beaten by at most two later
vertices, an 83-state automaton of Hankel rank 18, and THM-4748's closed formulas as the two-arc layer.

**No bounded-arc family can be complete (PROVED, counting; proof corrected by the audit, which found the first
draft's "a single arc of span s gives H >= 1 + 2^(s-1) and the insertion lemma transfers this to every sub-arc" false:
the word with arcs 3 -> 0 and 3 -> 2 on the chain 0..3 has H = 3).** In a connected cutwidth-2 word with c reversed
arcs the vertices between consecutive slot events are slot-free and lie inside a pending span, so by the sink-run
bound each of the at most 2c + 1 interior gaps has length at most log_2(X - 1) once H <= X; the at most 2c event
vertices are arranged in one of C(c) ways, so at most C(c) (log_2 X + 2)^(2c + 1) values below X have witnesses with c
arcs. Hence the number of reversed arcs of any witness family for all odd N must grow at least like
(1 - o(1)) log_2 N / (2 log_2 log_2 N): a proof of Chan-Pak through G(0,2) is not a finite list of closed formulas but
a numeration theorem in which the arcs are digits (the minimal witnesses use 1 to 5 arcs below 8192, growing slowly,
as this forces).

**The two-radix heuristic (not a theorem).** In every grade the odd-cycle collections contribute powers of 2 from a
single run through an interval (THM-4748's 2^p terms, the one-apex 2^(j-i-2)) and powers of 3 exactly when two runs
can interleave on an interval: two crossing reversed arcs (THM-4748's 3^r), an apex over a reversed arc, or two apexes
(Session 194's (3^L)^k). The grades with only the binary radix are G(1,0) (half of the odd numbers, section 5) and
G(0,1) (products of 1 + 2^(p-1): density zero); every grade computed that carries both radices, G(0,2), G(1,1), G(2,0),
G(2,1), is complete in its computed range. The conjecture behind HYP-9470 (C) is therefore that the two-radix sums
produced by cutwidth-2 words are a complete numeration system for the odd numbers above 21, in the spirit of the
classical completeness of sums of distinct numbers 2^a 3^b. Each is the range question for one explicit
N-rational series: finite alphabet, explicit nonnegative matrices, Hankel rank 17 for G(2,0). Ranges of N-rational
series over a two-letter alphabet already include all squares (the series (bin w)^2 has dimension 3), so no general
decision procedure is available and the question has to be attacked with the structure of these particular matrices
(Session 194's periodic counters and windows are the first steps in that direction for G(2,0)). The grading supplies
the hierarchy in which to measure progress: the smallest grade in which each odd number first appears (section 6), the
seam ranks that bound how much a composition across c reversed arcs can transmit, and the Hankel ranks that bound the
dimension of any recursion.

## 8. Housekeeping

`agents/check_docs.py` cannot open THM-4861 (file name of 233 characters) from any worktree whose path prefix exceeds
about 25 characters on this Windows machine; the session mapped the worktree to a drive letter with `subst` and ran the
check through a wrapper that keeps `Path.resolve` on the drive letter. New theorem file names should stay short.

## 9. Audit

Independent adversarial audit (general-purpose agent, own code throughout; report
`scratchpad/s79_audit/AUDIT.md` of the session, copied to
`04-computation/experiments/spine_grammar_audit_20261104/AUDIT.md`): **PASS WITH CORRECTIONS**, 4 majors, 18 minors,
5 upgrades, all applied (MISTAKE-720). Majors: (M1) prior art missed although the day's log had been read at the
start: mac-mini's ribbon lane (THM-4880, HYP-9473; committed 21:40-21:48, before the first checkpoint), Codex's
DEFECT_TRANSFER cutwidth transfer (2026-10-29), THM-4876 (6) and THM-4731 (6); the two novelty sentences were
removed and Theorem A retyped KNOWN as a method. (M2) the inference "G(2,1) attains every odd <= 2^32" from THM-4875
(D) was false (F1 lies in G(2,0) and has positive miss density alone; F2 is not a spine word); deleted everywhere.
(M3) the apex-free Hankel ranks depended on the empty-word convention; with the standard H(empty) = 1 they are 4,
17, 79 and the c = 0 seam is 1; the "+1, +1, +2 coincidence" was an artifact; tables now show both. (M4) the
counting lemma's proof had a false step (a single-arc bound does not transfer to sub-arcs: arcs 3 -> 0, 3 -> 2 on
0..3 give H = 3); rewritten through the sink-run bound on interior gaps with the arrangement constant. Minors applied:
the state bound restated (fails only at c = k = 0; the audit's exact count); G(0,1) row 358/664; the modulo-7 bias of
the one-apex misses; P3's time order; the sink-run and padding wordings; the length-cap argument for exact minima;
the validation range (order <= 7); the alphabet and dedupe remarks; the ordinal-sum cuts behind the length-64
two-apex witnesses; THM-4880's junction ranks cited in HYP-9470 (D). Upgrades: the audit reproduced, with its own
subset-DP searches, the censuses to 511 (all grades) and to 2047 (G(0,2), G(2,0)), the one-apex and interval-chain
spectra, the 2^21 count, all listed witnesses (0 mismatches), the ranks from the dumps, and the seam ranks 3 and 13
from brute-force Hankel blocks (lower bounds meeting the automaton's upper bounds).
