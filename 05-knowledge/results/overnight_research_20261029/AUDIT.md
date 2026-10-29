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

## Finite-clock continuation

`finite_clock_ports.py` passes22,951 active gates normally and with -O,
with identical stdout hashes in finite_clock_receipt.json. It checks24
states: sixteen C3 star/triangle sources, a two-block source, a genuinely
nonorthogonal clock, dependent and independent-isotropic ports, three
plane clocks and a noncommutative finite group on a weighted form. The
last example is not asserted to be a tournament automorphism group.

Each actual query is evaluated by rational elimination and independently
by matching Pfaffians. Fixed-space linear equations, orbit-span rank and
the exterior-square character formula agree. The cycle formula agrees
with signed pair orbits for507 partitions through14 and fourteen named
longer cases. There are6,250 specializations of the all-word unipotent
hostile. Literal re-marking/update/sign words, terminal singular queries
and an illegal unaveraged pair probe exercise the boundaries.

The proof uses finite order and the full invariant integer coupling space.
The old Reynolds physicality hostile is explicitly addressed by retaining
invariance of the unknown and legality of the orbit sum. No new arithmetic
Selmer map, native-complement operation, or Sha bound is inferred.

## Matched ternary-clock and precision continuation

`eisenstein_clock_joins.py` passes567,405 active gates in identical normal
and optimized runs. It exhausts norm shells through4096 (14,845 labelled
points), verifies factor and divisor-character counts independently, and
checks10,800 paired-group count queries for d,e=1..6 and n=1..50 against
534,420 source/coupling pairs. That norm universe contains every possible
required shell: de+s^2 n<=36+36*50=1836. It does not impose an accidental
weight cutoff on the supposedly complete group counts.

The chart test exhausts625 arbitrary small cross matrices. Direct rational
determinants decode856 nonsingular source states, including188 singular
probe values, and900 unit re-markings preserve the chart and content.
The paired-group observation is explicitly UNMARKED: first-port orders1
and4 refute a fully marked reading of the same carrier. A single internal
block direction restores the source factors. A Fraction/integer gcd type
error in the first checker run was repaired before filing by verifying
integrality before conversion; no mathematical statement changed.

`clock_join_precision.py` passes553,697 gates normally and with -O. It
exhausts all residue triples for the fifteen moduli listed in its JSON,
checks the exact ternary equivalence relation through81, and compares
the small grid with every residue query through modulus9. At76 all438,976
states yield54,872 complete response functions. The modulo228 count is
proved via CRT, not a claimed11,852,352-state census. The new all-size
3-adic proof uses finite differences and one precise lift step; the
dyadic pattern was then only a finite census; the continuation below
now supplies its all-order proof.

The incoming two-saturation verifier initially rejected the checkout's
CRLF bytes. Their LF normalization equals the exact committed blob and
the stored source SHA256. `replay_incoming_saturation.py` instead extracts
the six pinned ORIGINAL Git blobs and runs the UNCHANGED verifier, both
normally and optimized. It passes125,459 checks, rejects six hostiles and
certifies the two index8192, quotient(C2)^13 extensions. Its source hash
guard is not bypassed or weakened. This is an inherited arithmetic replay,
not a new rank, whole Mordell-Weil group or Sha-finiteness theorem.

After incoming d87e5a8292, `overnight_rank_return_20261029/atom_precision.py`
also replayed successfully without regenerating its receipt. It verifies
the actual endpoint Smith forms, the seven-prime Hamming map and the81/27
modulo9 separator. These are credited incoming results. Our nonlinear
norm-response precision recurrence is a separate operator and does not
replace the endpoint or ordered multi-run carrier.

## Dyadic, field and fully named continuation

`dyadic_field_ports.py` passes17,803,214 active gates. Its complete
dyadic universes are all triples modulo1,2,4,8,16,32,64, with exact
coordinate-class bijections, fiber sizes and parent branching. Complete
query banks through8 check the finite-difference certificate separately.
Selected entire lift fibers through2^18 and the opposite source-parity
strata at8 and16 exercise the proof's exceptional step.

All64 F4 states,262,144 F64 states and16,777,216 F256 states are decoded
from three direct field outputs. Irreducibility and every multiplication
pair are independently checked by polynomial gcd and long division.
The C3 clock splits over the two even-degree fields, and Frobenius
exchanges its two roots; this is explicitly checked, not silently
inherited from the rational irreducible representation.
The larger ring/CRT counts are theorem evaluations, not full censuses.
`DYADIC_AND_FIELD_PORTS.md` proves the arbitrary-precision statement.

`marked_clock_ports.py` passes66,732 active gates on3,512 nonsingular
signed source states. Direct Gaussian elimination evaluates every axis
probe and both named pairing residues. There are376 singular branch
observations, correctly treated as terminal symbols. Named port orders
recover absolute source factors; a large legal cross probe recovers
their signs. Exactly2,472 tested sources reject using the primal T
as the clock on dual cokernel classes. The correct action is T^(-T).

Both checkers pass normally and with -O with identical outputs. These
are same-agent alternative-path audits, not an independent referee.
The older modular checker was rerun in both modes after its status
metadata changed from OPEN to the proved dyadic continuation.
