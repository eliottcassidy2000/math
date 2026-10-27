# ETH at Σ⁰₄, the tournament on four vertices, galactic algorithms, and the micro/macrocosm

*mac-mini-2026-10-27-slab. A REFLECTION, written to think rather than to compute. It answers three parts of the
owner's directive:*
1. *how the recent disproof of ETH, "which quantifies over algorithms and sits at roughly Σ⁰₄", relates to the
   repository's many attempts at connecting the tournament on four vertices with other concepts;*
2. *how recent progress on galactic algorithms and matrix multiplication connects to repository concepts, with
   extensions found through synthesis;*
3. *how the first few handfuls of tournaments, and the micro/macrocosm duality, relate to a pasted text on the
   arithmetic hierarchy at levels 4–10, the analytical hierarchy and set theory.*

*Ground rules:*
- *The ETH refutation is an EXTERNAL PREPRINT CLAIM, kind-pasteur trust tier T3, unverified
  ([external note §3](../05-knowledge/results/external_preprints_mub6_local_smoothing_eth_merge_20261027.md)).*
- *Literature items marked "from memory" must be checked before anyone relies on them.*
- *Every link below is typed, and nothing here changes a status.*
- *The session's theorems are [THM-4819](../01-canon/theorems/THM-4819-regular-tournaments-under-native-moves-the-transitive-slab-exact-distance-and-the-prescribed-black-core.md),
  [THM-4820](../01-canon/theorems/THM-4820-eisenstein-merge-debt-floor-q-over-three-free-neutral-catalysts-and-the-eleven-forty-one-threshold.md)
  and [THM-4823](../01-canon/theorems/THM-4823-the-four-4-vertex-tournaments-census-by-scores-and-one-trace-pattern-free-classes-and-w4-free-equals-locally-transitive-plus-a-sink-tail.md).*

## 1. What the four-vertex record has in common

The repository has linked the four 4-vertex tournaments (TT4, W4, L4, S4) to at least fifteen other objects:
- the Klein four-group of two bits (THM-584);
- the end defects of trees on five vertices, with H = 2u + 1 = 1, 3, 5;
- the B₂ lattice, graded by c₃ = 0, 1, 1, 2;
- the Gauss identity, through TT4 counts in doubly regular tournaments (THM-4780);
- c₄ = tr(A⁴)/4 (THM-118), and the P₇ formula H = (462 − tr A⁴)/2 (THM-133);
- Erdős–Moser and R(TT4) = 8, with P₇ the last TT4-free class;
- 4-Pfaffians of value ±1 or ±3 and chirotopes (THM-4602);
- the forced order-4 metagraph (THM-4729);
- H being spectral at order 4 (THM-1780);
- the Collatz-pair sheets;
- the V₄/Keller carrier;
- the K₄ atlas;
- the native-move diameter 7;
- the order-4 nuclear-norm maximizers.

**The common thread, now exact.** THM-4823 shows that everything the four types say about a tournament is three degree
sums, Σ C(s,3), Σ C(n−1−s,3) and C(n,4), together with **one trace**, tr(A⁴)/4. The trace is not a function of the
scores from order 5 on. The four pattern-free classes have recursive counts:
- S4-free: Narayana;
- TT4-free: empty from order 8;
- W4-free: locally transitive followed by a sink tail, counted from A000016.

So the 4-vertex tournament sits on the exact boundary between the score sequence (the tournament's Parikh vector of
out-degrees) and one matrix product. It never reaches the Hamilton-path level, where H is #P-hard. Every one of the
fifteen links lives at one of those two levels:
- **local or degree level:** the bits, the defects, the gradings;
- **one-trace level:** the Gauss identity, c₄, the Pfaffians.

This is the 4-vertex shadow of klein's two lower invariant levels (regular, semilinear). [ANALOGY between the two
three-level ladders; IDENTITY for the census statement itself.]

**The relation to the ETH preprint.** Its engine counts sparse 𝔽₂ systems of high treewidth in 2^{o(n)}. Counting
tournaments that avoid a fixed 4-vertex pattern is such a system, over the C(n,2) arc bits. THM-4823 makes that use
void: every 4-vertex pattern class is counted in polynomial time, or is empty. **NO TRANSFER at four vertices.** The
first place where a subexponential forbidden-pattern counter could say something is a 5-vertex pattern whose
avoiding class has no recursive description (§5, item 1).

**The shared shape.** What the prior attempts and the preprint share is a shape, not a result. c₄ = tr(A⁴)/4 is a
degree-4 polynomial in the arc indicators that counts a local pattern exactly. The preprint's low-degree replacement
forms do the same for constraints of arbitrary width, paying with an interpolation window. [TOOL RESEMBLANCE.]

**The logical shape.** ETH is Σ⁰₄. Its negation is Π⁰₄, a uniform family of algorithms A_δ with proofs of correctness
and running time. That is what the preprint claims to supply. The 4-vertex analogue, "for every 4-vertex pattern there
is a polynomial-time counter", is a theorem with a one-line proof. [ANALOGY only.]

## 2. Galactic algorithms and matrix multiplication

**(a) ω and the smallest strong tournament.** Counting strong 4-subtournaments costs one matrix product, so it runs in
O(n^ω). The exponent is ω < 2.3714 by the laser-method line of 2024 (from memory). These improvements over Strassen
are galactic: they never pay at the orders gentourng reaches. The census of the smallest nontrivial strong tournament
is thus formally tied to the frontier of algebraic complexity. [IDENTITY of cost; NO TRANSFER of content.]

**(b) Two kinds of rigidity separate on Paley.** This is the session's main synthesis.
- **Rank rigidity fails.** Valiant rigidity asks for many entry changes before the rank drops. Dvir–Liu (2019–2020,
  from memory) showed that circulant matrices, Fourier matrices included, are *not* rigid. For prime q, the skew matrix
  S of P_q is circulant, so P_q is not Valiant-rigid.
- **Nuclear-norm rigidity holds.** Session75's audited spectral theorem shows that S needs at least q²/4 − o(q²) arc
  edits before it admits a bounded-arity substitution parse. The mechanism is Lipschitz:
  - one arc reversal changes S by 2(E_ij − E_ji), whose nuclear norm is 4;
  - parse tournaments have nuclear norm O(n log n);
  - P_q has nuclear norm about q^{3/2}, the maximum (THM-4814 (C)).
- **The difference.** Rank is not Lipschitz. So the notion of rigidity that survives for circulants is the
  nuclear-norm one. [IDENTITY of the object; two distinct metrics.]

**Extension (a proposal, not checked).** Alman–Williams (2017, from memory) showed that the Walsh–Hadamard matrices are
not rigid, using their Kronecker structure. The skew-Hadamard doubling of THM-4557, D(T) = [[H, H], [−Hᵀ, Hᵀ]] with
Hᵀ = 2I − H, is a twisted Kronecker step. Its tower T_k is skew-Hadamard, so it has maximal nuclear norm. If the
Alman–Williams polynomial method extends to the twisted step, the tower is a second family that is far from parses but
close to low rank. The test: find the rank of T_k after r changes per row, for k ≤ 7.

**(c) A galactic proof used microscopically.** Tournament automorphism groups have odd order. This half is elementary:
an involution swaps a pair and would have to reverse the arc between them.
- **The deep half.** By Feit–Thompson these groups are solvable. That makes canonical forms of tournaments
  quasipolynomial (Babai–Luks 1983, from memory), long before Babai's 2015 result for general graphs.
- **Where the repository uses it.** THM-4812's vertex-transitive lemma uses Feit–Thompson.
- **Where it falls short.** THM-4819's Lemma B has an alternative route through odd order: a product of two blue
  reversals is an automorphism of explicitly computable cycle type. That route fails for m ∈ {5, 10, 13, 17, 25, …}.
- **What works instead.** The working proof is a micro statement: the action of one anti-automorphism on one
  neighbourhood, N⁺(c_t) = σ_t(N⁻(v)). Here a local identity does what the macro tool cannot.

**(d) The 2-adic digits of H are the preprint's ring.**
- **The digits.** OCF gives H = I(Ω, 2) = Σ_j 2^j i_j(Ω), so H mod 2^ℓ sees only sets of fewer than ℓ disjoint odd
  cycles. This is the same truncation as the support of Beigel–Tarui's idempotent mod 2^ℓ.
- **At four vertices.** No two odd cycles are disjoint, so H = 1 + 2c₃: TT4 1, W4 and L4 3, S4 5. The Klein-four
  grading is the first 2-adic digit.

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
- A complexity statement about a tournament problem, such as "no 2^{δn} algorithm counts H", would be Σ⁰₄, like ETH.

**Where tournaments meet the upper levels and the analytical hierarchy.**
- **Reverse mathematics.** The Erdős–Moser principle EM says that every infinite tournament has an infinite
  transitive subtournament. It is the infinite form of R(TT_k) < ∞, the principle whose micro witness at k = 4 is P₇.
  Over RCA₀, RT²₂ is equivalent to EM + ADS (Bovykin–Weiermann, from memory), and EM is strictly weaker than RT²₂
  (Lerman–Solomon–Towsner, from memory). So the macro shadow of the 7-vertex Paley tournament is a Π¹₂ principle
  strictly below Ramsey's theorem for pairs. [IDENTITY: the finite TT_k Ramsey numbers are the finite shadows of EM.]
- **Paris–Harrington.** The unprovability phenomenon needs unbounded exponent, the dimension of the hyperedges. Fixed
  exponent 2, which covers graphs and tournaments, stays provable in PA (from memory). The repository's only objects
  with growing exponent are its k-uniform hypergraph censuses (k up to 10). [REMARK.]
- **Goodstein and the hydras.** Ordinal-valued invariants (ε₀) prove Π⁰₂ termination statements that PA cannot. This
  suggests splitting klein's level "beyond" into **I3a, ordinal** (a well-founded measure finer than any semilinear
  one) and **I3b, nothing known**.
  - **Where the two-branch rows stand.** Collatz, Mahler and Antihydra have no known I3a invariant. Terras makes I1 a
    full shift, and I2 decides only cycles. So a termination proof for Collatz would have to supply an ordinal-type
    measure.
  - **Status of the split.** It is a proposal for the ledger, not a claim.
- **Cantor and Baire space.** The 2-adic integers carry the regular level: the cylinders are the basic open sets.
  THM-4581's almost-sure coalescence and the transfer operators of THM-4821/4822 are statements about measures on
  Cantor space. Collatz is their restriction to the integers, a null set. So here the "microcosm" (ℕ) sits inside a
  measure-zero corner of the "macrocosm" (ℤ₂), where the almost-sure laws are already theorems.

**The session's own micro/macro pairs.**
- **Lemma B.** One micro exception (T5, n = 5) and a uniform macro law for every n ≥ 7, proved by an argument mod 4 on
  one Hamilton cycle.
- **W4-freeness.** A 4-vertex micro condition controls a macro count, A000016's necklaces.
- **The Eisenstein catalyst.** The micro structure is three unit charges (the A₂ weights). The macro data are a cost of
  about n²/3, from the index of √−3·ℤ[ω], and a threshold of k ≥ 36 catalysts. Both the threshold and the cost 41 are
  decided by the inertness of two primes, 11 and 41.

## 4. Verdicts

- **The ETH disproof**, true or not, is a lever on no 4-vertex statement and on no row of the atlas. The meeting points
  are tool shapes (low-degree forms, 2-adic truncation), not theorems.
- **Matrix multiplication** enters the micro-world exactly once, in the S4 count.
- **Rigidity** is where synthesis produced something definite: rank rigidity fails for Paley, and nuclear-norm (parse)
  rigidity holds.
- **The micro/macrocosm duality** is concrete in reverse mathematics (P₇ → EM) and in Cantor space (ℕ inside ℤ₂).

## 5. Next steps

1. **The five-vertex census.** For each of the twelve 5-vertex types, find whether its count is determined by
   polynomial-time statistics: scores, tr(A³) to tr(A⁵), Σ C(s,4), and per-vertex 3-cycle counts. Count the
   pattern-free classes for orders ≤ 10. This is the first place a subexponential forbidden-pattern counter could
   matter.
2. **The rigidity test** of §2 (b), on the doubling tower T_k for k ≤ 7.
3. **The ledger split** of I3 into I3a (ordinal) and I3b (nothing known), if the atlas owners agree.
