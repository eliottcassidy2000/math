# Actual-source descent toward the tournament path spectrum

**Status: PROVED elementary extensions / FINITE-EXACT verification / OPEN global target. Same-agent audit; no independent human or agent review claimed.**

Let H(T) count directed Hamiltonian vertex words in a tournament. The target
is that its image is exactly the positive odd integers other than 7 and 21.
The two omissions are already proved; sufficiency remains **OPEN**. This is
[Chan–Pak, Conjecture 6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf),
credited there to bof and Gordon Royle. The existing repository constructors
already cover all allowed odd values through 99,999,999. This packet does not
increase that finite ceiling or claim to settle the conjecture.

## What this session adds

1. **An infinite two-atom amplification family.** For integers a,b,r >= 0 with
   a+b+r > 0, put k=9+5a+8b+r. An explicit tournament on k-a-2b-1 vertices has

   `H = 2^k - 9*5^a*19^b`.

   The inherited deficit-nine seed combines with actual strong atoms having
   (H,E)=(5,32) and (19,256), where E counts ordered two-path covers, including
   empty slots. The latter six-vertex atom was recovered from the complete
   earlier bank. These are the only nontrivial dyadic-E profile types among
   strong atoms with H<=127. SCC composition therefore classifies all dyadic-E
   parents whose individual strong components obey that cap, at every size.
   [Proof, entry guard and limitation](AMPLIFIERS.md).
2. **A complete compressed inverse with minimum construction cost.** In the
   inherited singleton/strong-four-vertex language, a binary residue recovers
   each candidate entire singleton suffix. One then removes one strong atom.
   Minimizing over actual predecessor branches gives an exact minimum-vertex
   recurrence; the algorithm is polynomial in binary input length for a fixed
   number of strong atoms. This is full-state recognition, not a proof that
   every desired scalar H enters the language. [Proof and code](DESCENT.md).
3. **An exact failure of naive least-resistance descent.** At state
   (F,q,g,h)=(17,12,0,5), the least-rank predecessor has F=7; at (49,48,0,5)
   it has F=21. At (19,14,0,5), the least-rank predecessor has the attainable
   scalar F=9 but impossible joint state (9,7,3,5). Valid alternatives exist
   in all three cases. Even a successful greedy descent can waste a vertex.
4. **A literal Jacobian-fibre comparison.** Each fixed block forgets one
   response coordinate. Retaining it gives an explicit rational linear
   automorphism with determinant H(B)^2 E(B). Its integer image still has
   nontrivial Smith guards; its actual-source subset has additional holes.
   A three-point fibre has realizable endpoints and an impossible middle.
   [Kernel, inverse, Smith form and source boundary](JACOBIAN_FIBRES.md).

The user's "path of least resistance" suggestion becomes an exact Bellman
minimum over legal histories, with admissible clock lower bounds. The
square25 carrier supplies the same methodological constraint: optimize paths
only after retaining the state on which the allowed operations are defined.
Its 33-state restricted-complement observer is not a universal spectrum
carrier. No electrical-network/tournament identification is asserted.

## Inherited frontier and the missing theorem

- [THM-1370 / h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md)
  supplies the all-size exclusions, with its repaired literature priority.
- [THM-4094 / hamiltonian-matching-deficit-and-two-prime-lane-completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md)
  reduces completeness to strong realizations of odd primes p!=7 and of 7p
  for odd primes p!=3. The carry atoms and the reduction are inherited.
- The [October 31b packet](../square_sum_spectrum_overnight_20261031b/README.md)
  supplies the four-state action, cut/interior guards and terminating inverse;
  [31h](../square_sum_spectrum_overnight_20261031h/README.md) supplies the complete
  strong-profile bank and the optimal deficit-nine source.
- Incoming [Session142](../session142_tournament_odd_spectrum_20261102/SYNTHESIS.md)
  independently clarifies why a unimodular algebraic inverse can miss the
  Boolean construction domain. Its source-guard and positive-grading results
  are credited, not rediscovered here. This packet adds the actual-state
  fibre hole and the concrete blockwise minimum-cost inverse.

The remaining task is a **uniform source-entry theorem**: for every allowed
target, find a joint response state with that F whose inverse reaches an
actual seed. Termination, an inverse over Q, a finite coverage bank and a
sparse infinite family do not imply entry. In particular, the new amplifier
family has only O((log X)^3) outputs below X; it cannot alone prove cofiniteness.

## Reproduce

From the repository root:

```text
python -B -X utf8 05-knowledge/results/spectrum_descent_jacobian_20261102/verify.py
python -B -O -X utf8 05-knowledge/results/spectrum_descent_jacobian_20261102/verify.py
```

For a single requested state, the constructor returns the minimum word and
full outgoing-neighbor bitsets, or rejects membership in this grammar:

```text
python -B 05-knowledge/results/spectrum_descent_jacobian_20261102/descent.py 17 12 0 5
```

The [receipt](verification.json) records agreement of forward Dijkstra,
one-step inversion and compressed inversion on 4,728 reachable states;
complete membership decisions on 17,836 admissible candidates; all 16 lifted
Smith forms; 80 atom-cut controls; and nine direct adjacency-DP family checks.
The all-size conclusions use the proofs, not extrapolation from those tests.
See [AUDIT.md](AUDIT.md) for sources, universes, failures and operation scope.
