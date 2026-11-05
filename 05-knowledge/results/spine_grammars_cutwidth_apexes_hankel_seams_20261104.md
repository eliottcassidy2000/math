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
| G(1,0) | one-apex words (feedback vertex number 1) | 15 | 4 | Session 16, Session 191, HYP-9468 (binary affine state) |
| G(0,1) | upset tournaments / Moon's interval chains | 14 | 5 | THM-4744; the all-parts-1 chain is THM-4865's tribonacci tower |
| G(2,0) | two-apex (two-feedback) words | 200 | 17 | Sessions 191-194; rank 17 reproduces Session 191's two-apex transfer rank |
| G(1,1) | one apex over an interval chain | 114 | 17 | new |
| G(0,2) | cutwidth 2, contains every fa <= 2 tournament of THM-4748 | 83 | 18 | new |
| G(0,3) | cutwidth 3 | 521 | 80 | new |
| G(1,2) | one apex over cutwidth 2 | 994 | 78 | new |
| G(2,1) | two apexes over an interval chain; contains the one-apex one-arc joins of THM-4875 (every odd <= 2^32 except 7, 21) | 2544 | 78 | new |

The **seam ranks** (the number of bilinear terms needed to compose two spine tournaments across a cut crossed by c
reversed arcs) are 3 for c = 1 (the Lorentzian one-arc join law h h' + R R' - e e', KNOWN), **13 for c = 2** and
**62 for c = 3**, the same in every ambient grammar computed. New PROVED items: the transfer theorem with its state
bound, the containment dictionary, and the sink-run lower bound H >= 2^r + 1. FINITE-EXACT: the ranks, the validations
and the censuses (section 6). What the session does **not** do: prove that any G(k, w) is onto the odd numbers other
than 7 and 21; that statement, for G(2,1), is HYP-9470 and would imply Chan-Pak.

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
  and symplectic compiler; klein's digit spectra (THM-4862/4865); mac-mini's THM-4876; kind-pasteur S9/S10. None of
  them grades tournaments by chain cutwidth, and none computes cut-state automata, Hankel ranks or seam ranks; the
  one-arc seam rank 3 is their Lorentzian law, and the (h, e, R, d) interface is the c = 0 cut of G(1, *) seen from a
  root (section 3.4).

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
  2209, 4063, 7473, 13745 for r = 1..16 is the tribonacci law h_r = h_(r-1) + h_(r-2) + h_(r-3) (new identification of
  the tower as a cutwidth-1 word; the law itself is klein's).
- **fa <= 2 tournaments** (THM-4748) are the G(0,2) words with at most two slots; THM-4748's five configurations are the
  cases: disjoint (two connected words in ordered join), touching (same vertex is a tail and a head), crossing (two
  pending slots with distinct heads and distinct tails), nested, shared endpoint (one head group of size 2, or one
  letter terminating two slots).
- **One-apex one-arc joins** (THM-4875 (D)) lie in G(2,1): for one-apex words b on the chain C and b' on C', the join
  C => C' with one arc v' -> v reversed is the word with letters (b_i, 0) on C and (1, b'_j) on C' and a single slot
  from v to v' (x_1 = the apex of the left word beats all of C'; the apex x_2 of the right word is beaten by all of C).
  Hence G(2,1) attains every odd number <= 2^32 except 7 and 21 (THM-4875 (D), HYP-9468 CONFIRMED), which is the
  strongest finite evidence for any single grade.
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

| grammar | states | reachable dim | observable dim | Hankel rank |
|---|---|---|---|---|
| G(1,0) | 15 | 9 | 6 | 4 |
| G(0,1) | 14 | 8 | 6 | 5 |
| G(2,0) | 200 | 125 | 25 | 17 |
| G(1,1) | 114 | 65 | 25 | 17 |
| G(0,2) | 83 | 39 | 25 | 18 |
| G(0,3) | 521 | 211 | 120 | 80 |
| G(1,2) | 994 | 505 | 116 | 78 |
| G(2,1) | 2544 | 1549 | 116 | 78 (modular, two primes) |

Unexplained coincidences to be tested on larger grades: rank G(1,1) = rank G(2,0) = 17 and rank G(1,2) = rank G(2,1)
= 78 (an apex and a slot contribute alike once one apex is present), while G(0,1) = 5, G(0,2) = 18, G(0,3) = 80 exceed
G(1,0) = 4, G(1,1) = 17, G(1,2) = 78 by 1, 1, 2. G(2,0)'s rank 17 equals the Hankel rank of the two-feedback series found by Session 191 by a different method (the
"rational dimension 28" there is the reachable space of its 36-state representation). G(1,0)'s rank 4 is Session 16's
four coordinates.

**Seam laws.** Cutting a word at a cut with c pending slots, H(l r) = sum over cut states sigma with c slots of
(e_0 M_l)[sigma] (M_r P f)[sigma]: a bilinear composition law whose minimal number of terms is the **seam rank**, the
rank of the pairing restricted to prefixes ending at a c-slot cut.

| ambient grammar | c = 0 | c = 1 | c = 2 | c = 3 |
|---|---|---|---|---|
| G(0,1) | 2 | 3 | | |
| G(0,2) | 2 | 3 | 13 | |
| G(0,3) | 2 | 3 | 13 | 62 |
| G(0,4) (3,384 states, Hankel rank 390, modular) | 2 | 3 | 13 | 62, and c = 4: 310 |
| G(1,0) | 4 | | | |
| G(1,1) | 4 | 13 | | |
| G(1,2) | 4 | 13 | 61 | |
| G(2,1) | 17 | 61 | | |

For k = 0 the c = 0 rank 2 is the ordered-join product plus the degenerate empty-suffix direction. The c = 1 rank 3 is
the Lorentzian law (h, e, R) x (h', s', L') of session29/THM-4866 (B): KNOWN, reproduced here as the minimal rank. The
**two-arc seam rank 13** and the **three-arc seam rank 62** are new; 13 is the same inside G(0,2) and G(0,3), so it is
the rank of the two-reversed-arc composition law over all pairs of spine tournaments of cutwidth <= 3 on the two sides
(covering the crossing, nested, shared-head and shared-tail configurations jointly). With one apex the seam ranks rise
to 13 for c = 1 (G(1,1) and G(1,2) agree) and 61 for c = 2 (G(1,2)): the apex is a third kind of crossing. An explicit
named basis of the 13 statistics is not derived here (open task).

### 3.4 Lemmas used by the census (PROVED)

- **Monotonicity (KNOWN).** H(T) >= H(T - v): every Hamiltonian path of T - v extends to one of T by inserting v
  (Redei's insertion). Hence the prefix tournament (apexes placed) bounds every completion from below.
- **Padding.** A chain vertex beaten by every apex and followed by no chain vertex is a sink; a first chain vertex that
  beats every apex is a source. Neither changes H.
- **Sink-run bound.** Let u_1 < ... < u_r be consecutive slot-free chain vertices beaten by every apex (letters
  (0, 0, all)), followed by a chain vertex u that beats some apex x (k >= 1), or, for k = 0, lying inside the span of a
  pending slot whose tail u comes later. Then for every odd subset S of {u_1..u_r} the cycle x -> S (increasing) -> u ->
  x (respectively u -> head -> S -> u) is an odd cycle, so T has at least 2^(r-1) odd cycles and, by the odd-cycle
  collection formula (THM-4876 (4): the empty collection and the singletons), H >= 1 + 2 * 2^(r-1) = 2^r + 1. A run
  with 2^r + 1 > X therefore cannot be continued within the target range, and a word ending in the run has the count
  of the word without it (padding).
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
495, 998, 2035, 4077, 8180, 16382, 32923, 65646, 131959, 263502 for j = 3..20, i.e. the miss fraction tends to 1/2,
with no bias modulo 16, 3, 5 or 7. The class is not sparse-missing (as 13 alone might suggest): a second radix (apex
or slot) is needed for half of all odd numbers, which is why G(2,0) and G(2,1) carry the conjecture.

## 6. Censuses (pre-registered; all chain lengths <= 64; values <= X; products of connected words added)

**The surprise.** The apex-free cutwidth-2 grammar G(0,2) attains **every odd number <= 2047 except 7 and 21** (1022
of 1024 with the closer-bound engine, 14,796,241 search nodes; 254 of 256 below 512 with the first engine, 289,792,317
nodes, 43,254,577 complete connected words; 83 cut states, Hankel rank 18). So does G(1,1) below 512 (254 of 256;
282,279,675 nodes). Prediction P1 (between 10 and 60 misses below 300) is **REFUTED** in the direction
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

| grammar | X | attained (with products) | misses besides 7, 21 | search nodes |
|---|---|---|---|---|
| G(0,1) | 2047 | 360 of 1024 | 662 | 49,184 |
| G(1,0) | 2047 | 763 of 1024 | 259 | 29,841 |
| G(2,0) | 2047 | 1022 of 1024 | 0 | 20,924,275 |
| G(0,2) | 511 | 254 of 256 | 0 | 289,792,317 (first engine); 163,451,097 without deduplication (closer bound) |
| G(1,1) | 511 | 254 of 256 | 0 | 282,279,675 (first engine) |
| **G(0,2)** | **2047** | **1022 of 1024** | **0** | 14,796,241 (closer bound; at most 5 reversed arcs per value: 1, 10, 127, 542, 328, 14 values with 0..5 arcs, upper bounds on the minima from the deduplicated run) |
| G(0,2) | 8191 | RUNNING | | |
| G(0,3), G(1,2) | 511 | RUNNING | | |
| G(1,1), G(0,3), G(1,2) | 2047 | RUNNING | | |

Outcomes of the pre-registered predictions: P1 REFUTED (G(0,2) misses nothing below 512; the [2^11, 2^12) clause was
not reachable with X <= 2047); P3 CONFIRMED (13 attained; no misses below 512); P5 CONFIRMED (13 is in G(0,2) and
G(2,0), not in G(1,0) or G(0,1)); P2 and P4 PENDING the G(0,3) and G(1,2) runs. The minimal-grade atlas
(`spine_grammar_atlas_20261104.py`) over G(1,0), G(0,1), G(2,0) at 2047: G(0,1) misses 86 odd numbers below 300
besides 7, 21 (first 11, 13, 19, 23, 35), G(1,0) misses 45 (first 13, 29, 35, 41, 43, 49, 59, 67, 71, 77); every odd
number below 300 other than 7, 21 is attained by G(2,0) and by G(0,2); the values attained by G(2,0) only among
{G(1,0), G(0,1), G(2,0)} below 120 are 13, 35, 41, 43, 59, 67, 71, 77, 83, 109, 113.

## 7. What this gives for Chan-Pak

Chan-Pak follows from the statement "some G(k, w) is onto the odd numbers other than 7 and 21", and the data make
G(2,1) (one-apex joins, 2^32), G(2,0) (230,519) and now **G(0,2)** (511 at filing, with a few reversed arcs per value)
the candidates. G(0,2) is the most structured of the three: no apex, every vertex beaten by at most two later
vertices, an 83-state automaton of Hankel rank 18, and THM-4748's closed formulas as the two-arc layer.

**No bounded-arc family can be complete (PROVED, counting).** A cutwidth-2 word with c reversed arcs is determined,
up to padding, by at most 2c integer parameters (the spans and gaps), each at most log_2 H (a single arc of span s
already gives H >= 1 + 2^(s-1), and the insertion lemma transfers this to every sub-arc), so at most (log_2 X + 1)^(2c)
values below X have witnesses with c arcs. Hence the number of reversed arcs of any witness family for all odd N must
grow at least like log N / (2 log log N): a proof of Chan-Pak through G(0,2) is not a finite list of closed formulas
but a numeration theorem in which the arcs are digits (the decoded witnesses use 2 to 7 arcs below 302, growing
slowly, as this forces).

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

PENDING.
