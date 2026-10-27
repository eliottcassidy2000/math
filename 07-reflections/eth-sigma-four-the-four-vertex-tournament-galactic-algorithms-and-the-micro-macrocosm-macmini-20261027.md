# ETH's logical level, the tournament on four vertices, galactic algorithms, and the micro/macrocosm

*mac-mini-2026-10-27-slab. A REFLECTION, written to think rather than to compute. It answers three parts of the
owner's directive:*
1. *how the recent disproof of ETH, "which quantifies over algorithms and sits at roughly Σ⁰₄", relates to the
   repository's many attempts at connecting the tournament on four vertices with other concepts;*
2. *how recent progress on galactic algorithms and matrix multiplication connects to repository concepts, with
   extensions found through synthesis;*
3. *how the first few handfuls of tournaments, and the micro/macrocosm duality, relate to a pasted text on the
   arithmetic hierarchy at levels 4–10, the analytical hierarchy and set theory.*

*Corrected after an independent audit (`04-computation/experiments/regular_transitive_slab_20261027/audit_four/REPORT_THM4823_LEDGER.md`,
§3). The main correction is that ETH in its standard form is Σ⁰₃, not Σ⁰₄; see §1. The file name keeps its
original slug so that existing links still resolve.*

*Ground rules:*
- *The ETH refutation is an EXTERNAL PREPRINT CLAIM, kind-pasteur trust tier T3, unverified
  ([external note §3](../05-knowledge/results/external_preprints_mub6_local_smoothing_eth_merge_20261027.md)).*
- *Literature items marked "from memory" must be checked before anyone relies on them.*
- *Every link below is typed, and nothing here changes a status.*
- *The session's theorems are [THM-4819](../01-canon/theorems/THM-4819-regular-tournaments-under-native-moves-the-transitive-slab-exact-distance-and-the-prescribed-black-core.md),
  [THM-4820](../01-canon/theorems/THM-4820-eisenstein-merge-debt-floor-q-over-three-free-neutral-catalysts-and-the-eleven-forty-one-threshold.md)
  and [THM-4823](../01-canon/theorems/THM-4823-the-four-4-vertex-tournaments-census-by-scores-and-one-trace-pattern-free-classes-and-w4-free-equals-locally-transitive-plus-a-sink-tail.md).*

## 1. What the four-vertex record has in common, and ETH's level

The repository has linked the four 4-vertex tournaments (TT4, W4, L4, S4) to fourteen other objects:
- the Klein four-group of two bits (THM-584);
- the end defects of trees on five vertices, with H = 2u + 1 = 1, 3, 5;
- the B₂ lattice, graded by c₃ = 0, 1, 1, 2;
- the Gauss identity, through TT4 counts in doubly regular tournaments (THM-4780);
- c₄ = tr(A⁴)/4 (THM-118), and the Z₇-circulant formula H = (462 − tr A⁴)/2 (THM-133);
- Erdős–Moser and R(TT4) = 8, with P₇ the last TT4-free class;
- 4-Pfaffians of value ±1 or ±3 and chirotopes (THM-4602);
- the forced order-4 metagraph (THM-4729);
- H being spectral at order 4 (THM-1780);
- the Collatz-pair sheets;
- the V₄/Keller carrier;
- the K₄ atlas;
- the native-move diameter 7;
- the order-4 nuclear-norm maximizers.

**The common thread, for the census-type links.** THM-4823 shows that the four type counts are three degree sums,
Σ C(s,3), Σ C(n−1−s,3) and C(n,4), together with **one trace**, tr(A⁴)/4. The trace is not a function of the scores at
any order n ≥ 5. The four pattern-free classes have explicit counts, both as classes and as labelled tournaments:
- S4-free: Narayana;
- TT4-free: empty from order 8;
- W4-free: locally transitive followed by a sink tail, counted from A000016.

So the 4-vertex census sits on the boundary between the score sequence and one matrix product. The census-type links
live at one of those two levels:
- **degree level:** the bits, the defects, the gradings;
- **one-trace level:** the Gauss identity, c₄, the Pfaffians.

Two of the links are not census statements: the R(TT4)/P₇ extremality and the nuclear-norm maximizers. Counting
Hamilton paths, H, is a different level again. Its #P-hardness is standard for general digraphs (from memory); the
audit found no reference for tournaments specifically.

This is the 4-vertex shadow of klein's two lower invariant levels (regular, semilinear). [ANALOGY between the two
three-level ladders; IDENTITY for the census statement itself.]

**The relation to the ETH preprint.** The preprint counts the solutions of k-CNFs in 2^{o(n)}. Counting tournaments
that avoid a fixed 4-vertex pattern is a 6-CNF with 8·C(n,4) clauses over the C(n,2) arc bits. That system is dense,
not sparse, but it is within the preprint's k-CNF scope.
- **What THM-4823 does to that use.** It makes it void. Every 4-vertex pattern class has a polynomial-time count, or
  is empty, for classes and for labelled tournaments alike. **NO TRANSFER at four vertices.**
- **Where a counter could first matter.** It could first say something at a 5-vertex pattern whose avoiding class has
  no known recursive description (§5, item 1). This is a proposal.

**The shared shape.** What the prior attempts and the preprint share is a shape, not a result. c₄ = tr(A⁴)/4 is a
degree-4 polynomial in the arc indicators that counts a local pattern exactly. The preprint's low-degree replacement
forms do the same for constraints of arbitrary width, paying with an interpolation window. [TOOL RESEMBLANCE.]

**The logical level, corrected.** Write Corr(A, φ, t) for the decidable predicate "A halts on φ within t steps with
the correct answer".
- **ETH proper.** The standard form s₃ > 0 reads ∃δ>0 ∀A ∀c ∃φ ¬Corr(A, φ, c|φ|^c 2^{δ n(φ)}). This is **Σ⁰₃**. Its
  negation is **Π⁰₃**: one algorithm for each δ, with no uniformity in δ and no proofs required.
- **The 2^{o(n)} form.** "No algorithm runs in time 2^{o(n)}" is **Π⁰₄**, and its negation is **Σ⁰₄**. This is the
  "roughly Σ⁰₄" of the directive's framing and of klein's reflection §1.
- **The preprint's claim.** One explicit algorithm with a computable time bound is **Σ⁰₂**.

The 4-vertex analogue, "for every 4-vertex pattern there is a polynomial-time occurrence counter", has a one-line
proof. Counting the pattern-free *tournaments* needs THM-4823 (C). [ANALOGY only.]

## 2. Galactic algorithms and matrix multiplication

**(a) ω and the strong 4-tournament.** Counting the strong 4-subtournaments S4 costs at most one matrix product, so it
runs in O(n^ω).
- **Why not C3.** C3 is the smallest strong tournament, but its count is a function of the scores.
- **The exponent.** ω < 2.3714 by the laser-method line of 2024 (2.371339, from memory). These improvements over
  Strassen are galactic: they never pay at the orders gentourng reaches.
- **What this says.** The S4 census is formally tied to the frontier of algebraic complexity. [Upper bound on cost
  only; no lower bound is claimed. NO TRANSFER of content.]

**(b) Two kinds of rigidity separate on Paley.** This is the session's main synthesis. The object is the same and so is
the metric, the Hamming metric on arcs or entries. The target classes differ.
- **Rank rigidity fails.** The first target is low-rank matrices. Valiant rigidity asks for many entry changes before
  the rank drops. Dvir–Liu (2019–2020, from memory) showed that circulant matrices, Fourier matrices included, are
  *not* rigid. For prime q, the skew matrix S of P_q is circulant, so P_q is not Valiant-rigid.
- **Parse rigidity holds.** The second target is bounded-arity substitution parses. Session75's audited spectral
  theorem shows that S needs at least q²/4 − o(q²) arc edits to reach one.
  - **The mechanism: trace duality.** If e arcs disagree with a parse tournament B, then
    ⟨S, B⟩ = q(q−1) − 4e ≤ ‖S‖_op ‖B‖_*, with ‖S‖_op = √q and ‖B‖_* = O(q log q) (Session75 RESULT (3.2)).
  - **Why Lipschitz is not enough.** A Lipschitz argument alone, an arc reversal changing the nuclear norm by at most
    4, gives only about q^{3/2}/4 edits. Lipschitz plus Frobenius gives q²/8 (THM-4814 (C)).

  [IDENTITY of the object and the metric; two different target sets.]

**Extension (a proposal, not checked).** Alman–Williams (2017, from memory) showed that the Walsh–Hadamard matrices are
not rigid, using their Kronecker structure. The skew-Hadamard doubling of THM-4557, D(T) = [[H, H], [−Hᵀ, Hᵀ]] with
Hᵀ = 2I − H, is a twisted Kronecker step. Its tower T_k is skew-Hadamard, with maximal nuclear norm N√(N−1) (checked
to order 128 in the audit). If the Alman–Williams polynomial method extends to the twisted step, the tower is a second
family that is far from parses but close to low rank. The test: find the rank of T_k after r changes per row, for
k ≤ 7.

**(c) A galactic proof used microscopically.** Tournament automorphism groups have odd order. This half is elementary
by Cauchy's theorem: an involution swaps a pair and would have to reverse the arc between them.
- **The deep half.** By Feit–Thompson these groups are solvable. That makes canonical forms of tournaments
  quasipolynomial (Babai–Luks 1983, from memory), long before Babai's 2015 result for general graphs.
- **Where the repository uses it.** THM-4812's vertex-transitive lemma uses Feit–Thompson.
- **Where it falls short.** THM-4819's Lemma B has an alternative route through odd order: a product of two blue
  reversals is an automorphism of explicitly computable cycle type. That route leaves m ∈ {4, 5, 10, 13, 17, 25, …}
  unsettled.
- **What works instead.** The working proof is a micro statement: the action of one anti-automorphism on one
  neighbourhood, N⁺(c_t) = σ_t(N⁻(v)). Here a local identity does what the macro tool cannot.

**(d) The 2-adic expansion of H and the preprint's ring.**
- **The expansion.** OCF gives H = I(Ω, 2) = Σ_j 2^j i_j(Ω), so H mod 2^ℓ sees only sets of fewer than ℓ disjoint
  odd cycles. This is the same truncation as the support of Beigel–Tarui's idempotent mod 2^ℓ.
- **At four vertices.** No two odd cycles are disjoint, so H = 1 + 2c₃: TT4 1, W4 and L4 3, S4 5.
  - c₃ is the coefficient i₁, not a binary digit: H(S4) = 5 = 101₂.
  - H mod 4 sees only c₃ mod 2.
  - H(W4) = H(L4), so H does not carry the Klein-four labels; the two bits are source and sink, not 2-adic digits.

[IDENTITY of the ring, not of statements.]

**(e) Galactic decision procedures.** Tao's n^{O(n²)} reduction for LRC and the busy-beaver sizes are galactic bounds.
The census-horizon column added to the ledger (PROBLEM_CLASSIFICATION §5) measures how far below them the
repository's finite checks actually reach.

## 3. The first handfuls of tournaments and the hierarchy

**Counts.**
- Classes, orders 1 to 10: 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056 (A000568).
- Regular classes, orders 3, 5, 7, 9, 11: 1, 1, 3, 15, 1223.

**Logical levels.**
- A census fact is Δ⁰₀.
- The repository's tournament conjectures are Π⁰₁, since each order is a finite check. Examples: H ≥ disc;
  black-core connectivity; HYP-9423 (7), that every regular tournament is exactly (n²−1)/8 from TT.
- Asymptotic statements phrased as limits are Π⁰₃ or Σ⁰₃. Examples: the growth of the TT_k Ramsey numbers;
  n·CV(H)² → 2 (THM-589).
- A complexity statement about a tournament problem: "no 2^{δn} algorithm counts H" is Π⁰₂ for fixed δ, and Σ⁰₃ with
  ∃δ, like ETH proper. Its "no 2^{o(n)}" form is Π⁰₄.

**Where tournaments meet the upper levels and the analytical hierarchy.**
- **Reverse mathematics.** The Erdős–Moser principle EM says that every infinite tournament has an infinite transitive
  subtournament. It is the infinite form of R(TT_k) < ∞, the principle whose finite witness at k = 4 is P₇. The audit
  checked the following on the web:
  - over RCA₀, RT²₂ is equivalent to EM + ADS (Bovykin–Weiermann; Montalbán showed ADS suffices in place of CAC);
  - EM is strictly weaker than RT²₂ (Lerman–Solomon–Towsner);
  - EM is a Π¹₂ principle.

  So the infinite shadow of the 7-vertex Paley tournament is a principle strictly below Ramsey's theorem for pairs.
  [ANALOGY: the finite TT_k Ramsey numbers are finite shadows of EM.]
- **Paris–Harrington.** The unprovability phenomenon needs unbounded exponent, the dimension of the hyperedges. Fixed
  exponent 2, which covers graphs and tournaments, stays provable in PA (from memory). The repository's only objects
  with growing exponent are its k-uniform hypergraph censuses (k up to 10). [REMARK.]
- **Goodstein and the hydras.** Ordinal-valued invariants (ε₀) prove Π⁰₂ termination statements that PA cannot. This
  suggests splitting klein's level "beyond" into **I3a, ordinal**, a well-founded measure beyond I2, and **I3b,
  analytic, geometric, or nothing known**.
  - **Where the two-branch rows stand.** Collatz, Mahler and Antihydra have no known I3a invariant. Terras makes I1 a
    full shift, and I2 bounds cycle lengths.
  - **Possible routes.** A termination proof could come through an I3a measure or through an I3b (analytic) argument;
    neither is known.
  - **Status of the split.** The split is a proposal for the ledger, and PROBLEM_CLASSIFICATION §5 now uses it.
- **Cantor and Baire space.** The 2-adic integers carry the regular level: the cylinders are the basic open sets.
  THM-4581's almost-sure coalescence and the transfer operators of THM-4821/4822 are statements about measures on
  Cantor space. Collatz is their restriction to the integers, a null set. So here the "microcosm" (ℕ) sits inside a
  measure-zero corner of the "macrocosm" (ℤ₂), where the almost-sure laws are already theorems.

**The session's own micro/macro pairs.**
- **Lemma B.** One micro exception (T5, n = 5) and a uniform macro law for every n ≥ 7, proved by an argument mod 4 on
  one Hamilton cycle. The even analogue (THM-4819 (E)) has its own micro exception, S4 at n = 4.
- **W4-freeness.** A 4-vertex micro condition controls a macro count, A000016's necklaces.
- **The Eisenstein catalyst.** The micro structure is three unit charges (the A₂ weights).
  - **The macro data.** A cost of at least ⌊n²/3⌋, from the index of √−3·ℤ[ω] (conjecturally n²/3 + o(n²)), and a
    threshold of k ≥ 36 catalysts at n = 11.
  - **What decides them for n = 11.** The threshold cost 41 is decided by the inertness of two primes, 11 and 41
    (THM-4820 (C); derived concurrently in Session75's EISENSTEIN_INVOICES).

## 4. Verdicts

- **The ETH disproof**, true or not, is a lever on no 4-vertex statement and on no row of the atlas. The meeting points
  are tool shapes (low-degree forms, 2-adic truncation), not theorems.
- **Matrix multiplication** enters the micro-world in one place: the S4 count, at most one product beyond the scores.
- **Rigidity** is where synthesis produced something definite. Rank rigidity fails for Paley, and parse rigidity holds,
  by trace duality against the nuclear norm.
- **The micro/macrocosm duality** is concrete in reverse mathematics (P₇ → EM) and in Cantor space (ℕ inside ℤ₂).

## 5. Next steps

1. **The five-vertex census.** For each of the twelve 5-vertex types, find whether its count is determined by
   polynomial-time statistics: scores, tr(A⁴), tr(A⁵), Σ C(s,4), and per-vertex 3-cycle counts. (tr(A³) is a function
   of the scores.) Count the pattern-free classes for orders ≤ 10.

   **First data** (`04-computation/experiments/regular_transitive_slab_20261027/pattern5.py`, orders 5–8; reproduced
   by the audit):
   - **Order 6.** The triple (sorted scores, tr A⁴, tr A⁵) determines the counts of only two of the twelve 5-vertex
     types: TT5, and the (1,2,2,2,3) type with c₅ = 3, H = 15 and |Aut| = 3.
   - **Orders 7 and 8.** It determines none of them. The audit shows more: even the sorted scores together with the
     full spectrum (all of tr A³, …, tr Aⁿ) determine none.
   - **Pattern-free classes at order 8.** They range from 139 to 2,971 of the 6,880 classes per type. The R5-free class
     is the largest, at 2,971.

   So the step from four to five vertices is exactly where the census leaves the scores-plus-traces level.
2. **The rigidity test** of §2 (b), on the doubling tower T_k for k ≤ 7.
3. **The I3 split.** Report on how I3a (ordinal) and I3b (analytic, geometric, nothing known) work in practice in the
   ledger.
