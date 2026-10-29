# Proof and computation audit

Status: same-agent alternative-path review, not an independent-agent or
Lean audit. All checkers use active exceptions and work under Python -O.

## Universes and independent paths

1. Producer: every one of 65,536 labelled sign cuts of the fixed M4,M4.
   Rank modulo101 equals rational rank because every minor has absolute
   value at most16 (Hadamard). Schur Pfaffians and principal Pfaffian
   ideals give the paired invariant factors. No isomorphism quotient or
   arbitrary tournament family is inferred.
2. Audit: generates every rank-two sign cut through its unique pair of
   normalized row directions, independently giving6,272. Direct matching
   expansion gives the eight by eight Pfaffian. Gcd of six by six
   principal Pfaffians is1 for every case, proving only one nontrivial
   hyperbolic factor. The complete histogram and C3-invariant subcensus
   agree with the producer.
3. Weighted lattice producer:6,900 constructive realizations. Audit:
   all9,520 nonsingular choices d,e=1..4 and B entries=-2..2; ordinary
   determinantal ideals agree with (s,s,|t|/s,|t|/s), independently of
   the producer's Pfaffian-ideal code.
4. Bounded-weight recurrence: all d,e=1..6, L=0..3 and every possible
   paired invariant pair, including zero-count cases;109,080 directly
   enumerated nonsingular cuts. Exact packed convolution also conserves
   all104,060,401 matrices at L=50. Minimal L for H25,H27,H57 is4,4,6.
5. Parser:100,000 positive odd sources through199,999, plus360
   source/queue pairs. Every completed word retains its exact source
   cylinder, affine carry and both clocks. Audit checks the worst-letter
   rational inequalities, the four small-core paths and36 further
   source/queue pairs up to521-bit sources. The theorem uses the
   symbolic inequalities, not extrapolation from this finite sample.
6. F2 geometry: all two-dimensional isotropic subspaces, their incidence
   degrees and the six selected by a fixed split quadratic refinement.

## Errors caught before this checkpoint

- The inherited displayed M4 has Pfaffian -3, not +3. A direct eight
  by eight Pfaffian control caught the sign in the new Schur calculation.
  Prior paired groups and the product9 survive. The original note is
  corrected, and the mechanism is recorded in MISTAKES.md.
- A draft audit encoded the exponent in the general contraction window
  incorrectly. The proof inequality3^d<2^(r+d-1) was correct; its test
  is now `(3/2)^d<2^(r-1)`, and the (3,2) hostile fails as required.

## Scope gates

Integer weights, sign cuts, clock-invariant cuts and arithmetic Selmer
maps are distinct types. An unmarked paired-group isometry does not
transport specified ports. A determinant histogram does not preserve
the full matrix. A terminating auxiliary parser does not imply source
descent;3->5 and the positive3n-1 cycle5->7->5 remain hostile.

The rank30 certificate replay is inherited from August, not a new curve
discovery. Its public ordered points and coefficients match the live
curve273 JSON exactly. A rank lower bound is not an upper bound or a
proof of Sha finiteness. Current leaderboard counts are selected records,
not a population distribution. Incoming Session97 supplies stronger
independent arithmetic certificates and is integrated separately.
