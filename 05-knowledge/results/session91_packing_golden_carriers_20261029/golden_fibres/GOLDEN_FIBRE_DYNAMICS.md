# The golden reader is not always a dynamical quotient

Status: PROVED elementary consequence of the inherited all-period golden-reader theorem; FINITE-EXACT checks below. No external priority claim is made for automaton congruence or Markov lumpability. The new application distinguishes the operations that the arithmetic reader does and does not preserve.

## 1. Source, target and inherited theorem

Let phi satisfy phi²=phi+1. For n>=1 let W_n be binary words indexed cyclically by Z/nZ with no cyclic adjacent 11. Put

    R_n = Z[phi]/(phi^n-1),
    q(w) = sum_i w_i phi^i in R_n.

This reverses the exponent order used in some older notes; cyclic adjacency and the fibre statement are unchanged. Rotation is correspondingly multiplication by phi or its inverse, according to the chosen shift direction.

The inherited theorem is THM-4592, with its Session74 all-period correction:
`01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md` and
`05-knowledge/results/session74_eleven_mechanisms_20261027/session74-clocks/RESULT.md`.
It says q is bijective for odd n. For even n it is onto with exactly one nonsingleton fibre:

    q^{-1}(0) = {0^n, 1010..., 0101...}.

All other fibres are singletons. Its classical Fibonacci arithmetic-coding predecessors are acknowledged in the Session74 primary-literature note. We use this theorem as an inherited input, not a fresh priority claim.

The two alternating words do read as zero directly: phi²-1=phi is a unit and

    (phi²-1) sum_(j=0)^(n/2-1) phi^(2j) = phi^n-1.

The residue ring has cardinality L_n-1-(-1)^n. At n=2 it is the one-element ring; at every even n>=4 it is nonzero, so every power of phi remains a nonzero unit.

## 2. Specify the operations instead of just the cardinality

For each named position i define the totalized toggle T_i on W_n:

- remove bit i when it is1;
- add it when it is0 and both cyclic neighbours are0;
- otherwise leave the word unchanged.

Each T_i is an involution. The observable is only the current ring value q(w). In particular a blocked attempted toggle is not separately announced by a legality flag.

Call two states equivalent if every finite sequence of these named toggles leads to the same observed ring value. This is the smallest state carrier that retains q and closes under the specified future operations. It is the standard observation-equivalence definition for a finite deterministic machine.

**Theorem 1.**

1. For odd n, the reader already gives a faithful closed carrier because it is bijective.
2. For even n>=4, q alone is not a closed carrier. Its observation-equivalence refinement is the full set W_n. Every pair of distinct words is separated either immediately by q or after a single named toggle. More strongly, for any fixed position j, the pair of observations (q(w), q(T_j(w))) already reconstructs w.
3. For n=2, q is constant, so all words remain equivalent under this observer and these totalized operations. Adding a separately observable legality flag changes that conclusion.

**Proof.** Only the three zero-fibre words require attention. To separate zero from an alternating word A, choose i with A_i=0. Its neighbours are occupied, so T_i(A)=A still reads zero. But T_i(0)=e_i reads the nonzero unit phi^i. To separate the two alternating words A and B, choose i with A_i=1. T_i(A) reads -phi^i, whereas T_i(B)=B still reads zero. Thus the entire three-point fibre separates after one step. For the stronger fixed-j statement, its three output words are e_j, the alternating word with bit j zero, and the other alternating word with bit j removed. The middle output still reads zero. The first and third words are distinct nonzero words outside the zero fibre when n>=4, so the inherited singleton-fibre theorem gives distinct nonzero readings. The other fibres were already singletons. For n=2 the output ring is a singleton after every possible operation, proving the exception. QED.

An exact minimal repair for even n>=4 is therefore the ring value together with a three-way tag **only at zero**. This adds two states, not a full extra copy of the ring. A fixed-length tag for an individual zero-fibre choice needs two binary bits; it is not a claim that every state must carry two independent bits. The total cardinality is L_n rather than L_n-2.

## 3. Even randomizing the requested toggle does not remove the defect

Let P be the transition matrix obtained by choosing i uniformly and applying T_i. A partition is strongly lumpable when its transition probabilities into each partition block are independent of the representative state.

**Theorem 2.** For every even n>=4 the partition into q-fibres is not strongly lumpable for P. Specifically the transition probability into the zero fibre is

    P(0, q^{-1}(0)) = 0,
    P(A, q^{-1}(0)) = P(B, q^{-1}(0)) = 1/2.

**Proof.** Every zero-word toggle adds a nonzero unit. At an alternating word, exactly n/2 occupied bits can be removed, giving nonzero residues -phi^i; the other n/2 attempts are blocked and remain in the zero fibre. QED.

There is therefore no quotient Markov chain on R_n reproducing this chain from every initial word. A stationary conditional averaging construction would be a different contract, not a repair of strong lumpability.

For odd n the bijection permits exact transport of the whole chain, but the transitions are **guarded** additions or subtractions of phi^i. They are not automatically the translation-invariant Cayley walk. The companion operator lane analyses the resulting Lucas-cube graph and its eleven-state Paley embedding.

## 4. Why this bears on the other carriers

- The reader preserves the cyclic rotation clock at all n. Rotation invariance of the exceptional fibre does not imply closure under local toggles.
- At odd period5 the eleven words correspond exactly to R_5=F_11. The arithmetic identification preserves labels and rotation, but Euclidean packing contacts and feasible flips still require their own predicates. No physical packing optimality follows from this finite ring.
- In tournament tree compression, retaining the oriented frame makes expansion exact. Here retaining the exceptional source tag makes local updates exact. These are concrete state repairs, not an assertion that their groups or geometries are isomorphic.
- The Session90 reversal-frame family had a growing fixed-root ambiguity. Here the arithmetic fibre defect is uniformly bounded at three states, yet cannot be ignored for the given operation menu. Size of an ambiguity and whether it is observable are separate questions.
- Oddness enters through the norm/conjugate-root relation producing the three-word even-period collision. It is not the same hypothesis as odd-order tournament automorphism groups or coprime fixed-orbit descent. The shared discipline is to retain the exact quotient and action.

## 5. Reproduction

Run `python3 check_golden_fibres.py` and `python3 -O check_golden_fibres.py`.
The standard-library checker constructs Z[phi] as integer pairs and takes quotient residues using the adjugate of the multiplication-by-(phi^n-1) matrix. It independently enumerates all legal words through n=20, checks the inherited fibres and norm counts, verifies the claimed one-step observational refinement, and checks the exact randomized zero-fibre transition discrepancy.

The finite checks do not replace the inherited all-period proof. The all-order dynamic conclusions above follow from that theorem and explicit one-step witnesses, not from the finite range.
