# Connection contracts, incoming results and a positive-charge reroute

**Status:** proved elementary bridges and explicit OPEN targets; same-agent
review. Analogies below are not equivalences between open conjectures.

## 1. Recovered arithmetic: the diagonal is a fixed-point stratum

[THM-2422, operation fibres](../../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md)
already proves the requested doubles/squares dictionary. On positive pairs,
the swap `(a,b)<->(b,a)` fixes the diagonal. Its image under addition is
`2a`, and under multiplication is `a^2`. The off-diagonal half-quadrant
is a representative space for unordered pairs, not a new tournament.

That theorem also proves why dropping the diagonal makes prime squares
look prime in the strict proper-factor fibre. The first hostile is 4.
We reuse this result rather than promote its weaker rederivation.

For `1<=a<b`, the joint map `(a,b)->(s,p)=(a+b,ab)` retains the unordered
pair exactly: its entries are the roots of `X^2-sX+p`, and
`s^2-4p=(b-a)^2`. Either operation alone forgets a different parent fibre.
On a fixed-difference family b=a+d, the image is

    4p=s^2-d^2,    s=2a+d.

Moving one step along this family sends `(s,p)` to `(s+2,p+s+1)`.
Multiplication of both entries by c sends `(s,p,d)` to `(cs,c^2 p,cd)`.
These are different operations on the same labelled parent carrier, not
an identification of addition with multiplication. THM-4057, below,
already retains the trace/norm/signed-sheet version of these coordinates.

### A precise Farey test on a difference family

**Proposition 6 (PROVED elementary bridge).** Reduce the consecutive
fractions `a/(a+d)` and `(a+1)/(a+d+1)`, where a,d are positive integers.
Put `g=gcd(a,d)` and `h=gcd(a+1,d)`. Their reduced determinant has absolute
value

    d/(gh).                                             (13)

They are Farey neighbours exactly when `d|a(a+1)`. For fixed d there are
exactly `2^omega(d)` such residue classes a modulo d, with value 1 when d=1.

Proof: before reduction the determinant is -d. Dividing the two vectors
by g and h gives (13). The integers g,h are coprime, so gh divides d.
Equality gh=d means that every full prime power in d divides either a or
a+1; these alternatives are disjoint. Equivalently d divides a(a+1).
For each prime power dividing d, the allowed residues are 0 and -1;
CRT gives the stated count. QED.

This gives an actual map from a multiplicative factorization of a
**difference** to additive neighbour events along its ray. The count uses
the radical but the event locations require the full prime powers. For
d=9 the events are a=0,8 mod 9, not all a=0,2 mod 3. At a=2,d=9 the
unreduced determinant 9 survives unchanged despite the radical test passing.
This is exactly the kind of omitted coordinate that the complete smooth
probe bank restores in the counting theorem.

Source: labelled integer pairs with fixed difference. Target: consecutive
reduced fractions and their Farey edge predicate. Preserved: determinant
after the explicitly recorded gcd division. Lost: common content g,h,
absolute height if only the yes/no edge is retained. Sidecar: a,d,g,h.
Decisive test: all a,d in the declared finite rectangle, including d=9.
No Collatz itinerary or knot type is inferred from this determinant.

## 2. Farey and torus-knot histories need their integer lift

[THM-4057, Stern–Brocot depth and rational edges](../../../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md)
already separates primitive pairs from Farey adjacency: gcd(a,b)=1
is not the determinant-one condition between two rational nodes. It also
identifies reciprocal with a mirror, retaining depth and content. The raw
rational object contains both directions on coprime pairs and neither on
noncoprime pairs, so it is not intrinsically a tournament.

[THM-2646, braid-three modular central pullback](../../../01-canon/theorems/THM-2646-braid-three-modular-central-pullback-and-full-twist-knot-fibre.md)
is a much stronger prior bridge than a numerical resemblance of Farey and
knot diagrams: `B_3` is the compatible pullback of `PSL_2(Z)` and the
integer exponent sum over `C_6`. The same exact modular image can lift to
the unknot or trefoil; one fibre contains torus knots of unbounded genus.
The missing integer is central full-twist height. It is not a knot
invariant under arbitrary Markov stabilization.

This supplies a rigorous *analogy of representation obligations*: braid
modular reduction needs central height; Collatz block sorting needs the
actual word and K; the geometric wing palette needs labelled contact
vertices. There is no declared map taking braid height to Collatz time.
The cheapest failed transfer is already explicit in each lane: identical
quotient data can have different knot type, certificate debt or boundary
extension status.

## 3. A positive-charge Collatz reroute, and what it has not proved

In the coordinate z=n-1 the exact affine odd step is

    z'=(3/2^a)z + 2^(2-a)-1.

At a=1 the constant is +1, at a=2 it is 0, and at a>=3 it is the negative
reset dividend `-(1-2^(2-a))`. For an authenticated word w set
`A=3^m/2^K` and let A_(i+1:m) be its suffix slope. Then

    z_m = A z_0 + P(w) - D(w),
    P(w)=sum_(a_i=1) A_(i+1:m),
    D(w)=sum_(a_i>=3) (1-2^(2-a_i)) A_(i+1:m).           (14)

Both sums are nonnegative. Dropping D therefore gives a source-preserving
upper bound. This is an algebraic decomposition, not an independence
assumption about the exponents.

**OPEN hypothesis (first-exit positive charge).** Whenever a positive odd
n>1 has a first strict-descent word w, its positive-charge upper bound is
already paid:

    A(n-1)+P(w)<n-1.                                   (15)

The actual descent only gives the weaker inequality with `-D(w)` present.
No implication discarding D has been proved. Nor does (15), even if proved,
assert existence of the first exit. The latter is the universal Collatz
obligation, and must not be hidden in the word selection.

**FINITE-EXACT evidence:** the companion charge probe checks all 500,000
odd sources from 3 through 1,000,001 at their first strict descent. It uses
integer numerators at the actual clock and prints the worst records for
`D(w)/(n-T_w(n))`. Every checked word satisfies (15). This is a declared
finite universe, not a random-model or uniform ROOT-current theorem.
The largest ratio is `215402584767891/281474976710656` at source63.
The initial probe's clock-scaling error is repaired and recorded in
[MISTAKES.md](../../../01-canon/MISTAKES.md); the separate suffix-sum checker
verifies both (14) and (15) throughout the same universe.

Why this may be a useful reroute: the proposed sufficient charge bound
uses only positive exponent-one emissions and actual suffix discounts.
It does not need cancellation against reset dividends. A proof would
remove one source of sign cancellation from this renewal certificate.
The original ROOT transport has additional weights and sources; no identity
with that operator has been shown, so (15) is a helper target rather than
a renamed solution of its amplitude problem.

Hostile programme: seek an integer source with a tiny first-descent margin
and a large reset dividend transported through a long exponent-one run.
Formal words alone are insufficient: enforce the exact cylinder and all
proper-prefix non-descent inequalities. Then attempt a bound using those
integer inequalities. This is more focused than assuming independent
digits or bounding the positive and negative currents separately at every
arbitrary prefix.

## 4. Incoming work: actual transfers and their boundaries

The session began at eda06857a3, after integrating eleven incoming commits.
During work it integrated the further series through 3f6fb6a438. The relevant
proofs, rather than just commit titles, were read before these connections.

* [Golden junction certificate](../scc_basins_20261020d/arithmetic_junction_certificate.md):
  21 fixed Fano probes distinguish all 105 candidate junction labels, but
  their zero constraints do not certify an arbitrary tensor's target
  support. The full target-adapted 14-probe certificate survives. Our
  three-wing gadget has the same precise warning: every two-wing projection
  is full although the joint relation excludes 384 normalized assignments.
  Neither claim identifies cubic edge-colouring with six vertex colours.
* [Golden normalization audit](../scc_basins_20261020d/arithmetic_normalization_audit.md):
  compatible four-move paths, endpoint repairs and active word inversion
  make a real termination argument. Collatz has a unique physical successor;
  arbitrary reordering is not such a path. Our sorting is only an envelope.
  The shared move is to retain the chosen context and prove every intermediate
  step, not to import the tournament's terminating potential.
* [Hollow-frame concentration](../hollow_frame_blue_extrema_20261020/concentration/HOLLOW_FRAME_CONCENTRATION.md):
  a complete positive coefficient bound gives a relative variance of order
  1/k under uniform labelled frames. This is a model-relative mean statement.
  It supplies no source-specific Collatz bound, but suggests looking for
  a positive majorant like (15) before a probabilistic interpretation.
* [Odd-blue triangle moments](../odd_blue_spectrum_hollow_cubes_20261020/ODD_BLUE_TRIANGLE_MOMENT_LAWS.md):
  the exact product `q^c product(1+q^d_j)` has a zero of exact order r at -1.
  Cancellation survives the full cube; a selected corner need not inherit
  it. This is why source selection cannot be replaced by a phase trace.
* [Complete pure-blue atom classification](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md),
  now reconciled in [THM-2454](../../../01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md):
  the rigid strong alphabet is exactly `{1,C3}`, and the pure-blue strong
  alphabet is exactly `{1,C3,regular C5}` at every order. The earlier
  finite-but-unclassified numerator has been superseded by
  `P(z)=(1+z+z^3+z^5)/(1-z^2-z^6)`. Explicit even-order prefix reroutings
  force direct violated arcs; unit-step reindexing by 2 or -2 eliminates
  the residual sign patterns, including composite orders. This proves its
  finite exception bank from actual Hamiltonian structure. Collatz cannot
  assume a finite exceptional first-exit bank by analogy.
* [Effective black-defect bound](../effective_black_defect_order_bound_20261020/PROOF.txt):
  for a strong tournament with Delta black Hamilton-path orbits,
  `n<16(Delta+2)^2(Delta+1)`. Rotation-frame control plus a removable
  power-of-two cycle produces too many distinct black orbits. This Delta
  is not the feedback slack `3f-s` of the previous packet. The useful move
  is to turn an omitted family of probes into many certified distinct
  responses, as the smooth-bank construction does; there is no established
  replacement of black path orbits by Collatz odd charges.
* [Hamilton-path relocation connectivity](../hamilton_path_relocation_connectivity_20261020/THEOREM.txt):
  all labelled Hamilton words are connected by legal single-vertex
  relocations. Nonempty insertion fibres, lifted edges and access to a
  removable vertex prove an all-order induction. This closes necessity of
  the auxiliary local forest certificate, not black metagraph connectivity.
  The Collatz analogue would need both a legal lift and universal access to
  the smaller carrier. Our block refinements preserve an existing word;
  they do not supply that missing access.
* [Regular endcaps and D3 owners](../rotation_owners_spectrum_20261020e/README.md):
  regular outer endcaps now have short authenticated black bypasses; the
  complete deficiency-three owner stratum has a degree-seven normalized
  polynomial, with external matching partners retained. This reinforces
  the need for named ports in the wing gadget. Its explicit SCC-to-two-path
  matching bijection retains block decoding for Hamiltonian weights;
  replacing two C3 blocks by one R4 changes weight9 to5 despite a matching
  count resemblance. The corrected strong-path minimum and the still-open
  odd H-spectrum completeness are kept separate from these proved rules.
* [Complete mixed-prime census through degree thirteen](../complete_mixed_prime_census_degree13_20261020/CORRECTION_NOTE.md):
  the C15(5+3+1^5) full-type defect is84, not either extrapolation88 or64.
  Relative embedding multiplicity, distinct full types and escaping sign
  witnesses contribute separately. Our lower bound deliberately counts
  response histories and wordwise rank, not distinct tournament classes;
  the same separation prevents it being mistaken for a class-count lower
  bound or for a restriction on every possible A000568 algorithm.

### Extending the sharper counter to converse-merged node counts

The incoming [exact realizability theorem](../merged_metagraph_realizability_20261020/REALIZABILITY.txt)
uses the exact converse-merged node count `K_n=(A_n+S_n)/2`. Its Burnside
formula permits for S_n at most one length-one cycle and otherwise only
cycle lengths 2c with c odd. If e is the optional fixed-point count, adding
k cycles of length 2c has multiplier

    2^[kc+2c binom(k,2)+2k R_c+ek] / ((2c)^k k!),        (16)

where `R_c=sum gcd(d,c)` over earlier half-lengths d. Thus the same odd
gcd-response recurrence computes S_n at degree horizon floor(n/2), with
a separate e=0 or 1 channel. Positive weights and the allowed half-lengths
are unchanged; the smooth-partition upper bound applies verbatim. Arithmetic
bit-size remains polynomial in n as in the ordinary counter.

**PROVED consequence:** exact merged-node counts through N have the same
upper bound `exp(sqrt(N) exp(-c sqrt(log N log log N)))` for every fixed
c<1/2, using these response recurrences. The new width lower bound applies
to their ordinary component. This does not prove the incoming realizability
problem polynomial-time, nor does it replace its coordinate/transport
certificate with a node count. The independent small partition calculation
in the control script reproduces its n=3,...,9 table.

## 5. Research move used and next questions

The [META-PATTERNS](../../../00-navigation/META-PATTERNS.md) cards used were
"Search the statement before the method", "Correct the object before
sharpening the technique", and "Separate local support from bounded-height
coverage". No new universal card is promoted from one session.

The reusable successful move was to test **all natural future probes**
before accepting a compressed state: composite smooth lengths recovered
the hidden divisor coordinates and closed the leading width scale; named
rim ports recovered a ternary colour restriction invisible to pairwise
tests. Conversely, allowing unlimited terminal clock made a uniform
two-block budget conditional on convergence and therefore unhelpful for
proving arrival by itself.

Live questions:

1. Can prefix non-descent and the exact cylinder prove (15), or produce an
   integer counterexample? Can arrival be added without circularity?
2. What is the next-order term beyond the sharp `B_N/2` interface loss?
   Can a different, non-wordwise interface actually evade that rank bound?
3. Which exact non-lattice attachment couples three wings while forcing a
   forbidden rim pattern, with every port fixed and every cross edge checked?
4. Can the divisor split in (13) organize a sparse address grammar without
   dropping content or integer clock? The current determinant law alone
   asserts no dynamical entry.
