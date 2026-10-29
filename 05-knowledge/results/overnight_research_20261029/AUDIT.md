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
- The port clock convention was repaired from TP=PC to T^T P=PC:
  columns in the update PEP^T act as functionals. The original full-matrix
  congruence guard was correct. MISTAKES.md records the distinction;
  dependent ports retain their embedding before reducing the guard.

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

## Second checkpoint verification

- `port_queries.py` uses direct rational inverses and determinants to
  check the resolvent, pair-probe separation, independent/dependent
  port examples and the same-response/different-group hostile.
- `fibration_height.py` checks invariant degrees8/12, the degree24
  j-map, a modular coprimality certificate, the marked section identity
  and all1,959 reduced rational parameters through height40.
- `geometry_and_family_audit.py` exhausts16 quadratic refinements and
  the720 faithful S6 actions. Its independent SymPy1.14 path computes
  the polynomial gcd over Q, and its exact specialization compares the
  family invariants with the public minimal model4239 at the named
  parameter and coordinate scale. An initial normalization bug in the
  producer's modular gcd was caught before filing: its returned constant
  now reduces modulo p, not just multiplication by the modular inverse.
- The incoming verifier for curves273/302/4239 replays13,851 checks,
  rejects12 hostile mutations and accepts the three legitimate sign
  reversals. The old weighted checker replays61 gates, including81/27
  with binary orbit totals176/512 and A000568 through20; no new longer
  sequence prefix is claimed.
- Primary external manuscripts were read for the technique review.
  Their full proof/build graphs were not audited and are not dependencies
  of our two-port, parser or height arguments.
- After integrating incoming commit4d9c3ad1ce, three additional inherited
  scripts were replayed: the exceptional-cell recurrence through level1000
  with independent inverse cells through24; the cyclic-kernel audit with
  1,099,595 gates for all odd lengths3..127; and the complete66-curve
  local-character coverage check with140,334 gates,64 binary banks and
  two ternary repairs, certifying2,005 listed point directions. These are
  replays of incoming work, not additional curve or code discoveries.
