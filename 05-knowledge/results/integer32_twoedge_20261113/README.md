# Integer 0–32 recovery and the joint-incidence metagraph presentation

2026-11-13. This packet recovers a route for every integer0 through32 and
turns several numerical connections into exact maps or explicit obstructions.
The status of each promoted statement is bound to its independent review in
`promotions.json`; finite receipts have separately declared universes. No
general Hamiltonian spectrum, all-order black-core connectivity, BSD,
Hadwiger–Nelson, abc, RH, or LRC(14) closure is asserted.

## Start with the recovered census

[TABLE.md](census/TABLE.md) gives all33 rows, with source kind and a direct
route. [RECOVERY.md](census/RECOVERY.md) explains what was actually found:
26 rows in the older role index and seven marked supplements. Dedicated
lexical censuses and quoted number-specific requests are distinguished from
synthesized lookups. The individual request history for every number was not
recovered, and a lexical frequency is not promoted to a mathematical fact.

The useful recovered chains include6/7/8 through binary augmentation and F8;
15/30 through pairings and their two Fano completions;21/22/63/64 through a
norm observation and its missing scalar phase; and23/24/25/26 through the
separately typed code/lattice routes. Their maps retain different operations.

## A faithful, checkable observer for the metagraph

The inherited Session320 theorem says that positive joint use of adjacent
edges in one Hamilton path recovers a tournament up to converse, using the
directed-support reconstruction theorem for unused arcs. Its order22 and
strong order24 hostiles therefore have a precise missing coordinate: joint
incidence. We independently reviewed those statements and their scopes in
[INHERITED_AUDIT.md](reconstruction/INHERITED_AUDIT.md).

For the full local Boolean field beta, use **two optional path slots**.
Then beta(v;u,w)=1 exactly when the arcs at v point in opposite directions.
A negative one-slot answer does not have this local meaning. For n>=3 put
M=binom(n,2), L=n binom(n-1,2). The new
[realization and repair proof](reconstruction/PROOF.md) gives:

- An exact affine binary code of length L and dimension M-1. Its image is
  characterized by XOR0 on every three-edge star and XOR1 on every vertex
  triangle; a valid field has exactly two converse realizations.
- Minimum Hamming distance2n-4 and constructive correction of n-3 arbitrary
  bit errors in O(n^4) bit operations. Radius n-2 is sharply ambiguous.
- The inherited M-1 spanning-tree queries are still sufficient without
  errors; the additional full-field coordinates supply consistency and
  error correction. These are labelled binary-query statements.

The [native presentation](native/RESULT.md) then retains legal transport:
if x records reference-arc reversals and P is a Hamilton path's edge set,

    beta = beta0 + delta x,
    native move: beta' = beta + delta chi_P,
    Hamming distance(beta,beta') = 2(n-2)^2.

Here delta is the coboundary on L(K_n). A candidate jump is legal iff one
of its two complementary cut supports is a spanning path coherently directed
by the source. Connectivity and orientation are both necessary guards. The
support is unique except at order4, where complementary Hamilton paths can
give the same jump. This part overlaps the independently incoming Session321
wedge-wiring theorem; [INTEGRATION.md](INTEGRATION.md) records the credit.

For exact weights retain `(beta, ordered path)`, then quotient by vertex
permutations. This recovers the fixed-spine word carrier and degree
`m H(T)/|Aut(T)|`, with converse-fibre size m. Quotienting path reversal too
early gives40 instead of64 darts at order5. The complete ordinary/merged
native and blue matrices were independently rebuilt through6 and agree
entrywise with the earlier census. This is an exact all-order presentation;
its general orbit enumeration remains expensive.

## The arithmetic consumer is larger than a scalar

[The corrected bridge](bridges/RESULT.md) composes the faithful wedge decoder
with `Phi(M)=M^2+M+I`, where M is row-zero binary tournament adjacency.
Incoming Session321 strengthens the target: **the full labelled Phi and its
embedded kernel are converse-invariant**. The initial claim that named
Phi-kernel vectors needed an orientation choice was false and is retracted
in [CORRECTION.md](bridges/CORRECTION.md). The actual row-zero M still retains that choice. The doubled Monsky map
on the Phi kernel is invariant at even n and changes by (0,x) at odd n.

The familiar order5 codes8/10 have equal scores and full slot polynomial
`Q=t^5+6t^3+2t`, but Phi-kernel dimensions2/0. With the retained prime labels
`(3,19,43,139,251)` and `(3,19,83,499,1163)`, respectively, these are actual
pure 2-Selmer dimensions in the specified congruent-number family, via the
cited Monsky formula. They are not a claim about the two ranks or full BSD.
Ordinal joins with transitive suffixes prove that `(n,Q,scores)` loses this
invariant at **every n>=5**; the smallest starting order is checked exhaustively.

Incoming native transport also changes Phi by the exact path-Laplacian form
`ML+LM+L^2+L`, an alternating matrix killing the all-ones vector. Its local
stencil prunes candidates but does not replace the exact wedge path guard.
The [incoming review](INTEGRATION.md) joins these constraints with the
selective strong equal-H rays while preserving their different qualifiers.

## Exact F8/Fano maps and their obstructions

The exponent action `s -> a+2^b s` becomes
`x -> alpha^a x^(2^b)` on F8*. Thus1,2,4 is literally Frobenius. In the
recovered two-coordinate blocker module, the action is

    (u,v) -> (alpha^a u^(2^b), alpha^(-a) v^(2^b)).

Restricting to the Frobenius C3 gives a Phi3 kernel `Z direct-sum Z`, where
`Z={0,alpha,alpha^2,alpha^4}`, hence an F4^2 operator module. The full C7:C3
action does not preserve it and has no linear quotient of binary dimension2
or4. This is an exact symmetry boundary, not a Selmer-group identification.

[Thirty Fano systems](bridges/FANO_ORBITS.md) yield a complete finite action
classification:

    30 = 1 + 1 + 7 + 7 + 7 + 7.

The two fixed systems are translates of `{1,2,4}` and its negative. Each
remaining system is recovered from one of four intersection types and a
stabilizer centre in F7. The tempting1+1+7+21 split is refuted by the
order-three fixed-point count. The genuine7+7+21 split belongs to35 individual
triples, a different carrier. The remaining28 are not invariant under the
larger168-element Fano automorphism action, blocking an automatic Coxeter28
identification. HN's retained pair-channel and geometric-realization guards
remain essential.

Incoming Session321 uses another carrier:128 triangle-twist masks of U21.
We independently replayed all128, their canonical seven-triangle partition,
8-to-1 majority arcs, forced Frobenius dissent ports, and12 ordinary/6 merged
classes. These six regular classes have no native edges among themselves;
this leaves paths through outside classes available. The large Hamilton DP
counts remain inherited computations, not independently rerun here.

## The norm example becomes a general recovery theorem

The [F64 example](census/NORM_CLOSURE.md) has22 norm-observation cells and
closes under multiplication alone. Adding translation exposes its lost
ternary scalar phase in the exact sequence22->63->64.

The [general theorem](census/GENERAL_CLOSURE.md) starts with scalar orbits in
F_q^d, d>=2, translates by nonzero v, and also uses an invertible linear B
with Bv outside the line F_q v. For q>=4 the exact response-refinement counts
are

    1+(q^d-1)/(q-1) -> q^d-q+3 -> q^d.

The sharp boundary is zero strict rounds for q=2 and one for q=3. If the
line-moving condition is dropped, B=identity over a prime field F_p needs
p-2 rounds; over a proper extension it leaves a nontrivial residual fibre.
Thus two rounds here have an explicit reason, not a universal relationship
to sheafification, Bott periodicity, or tournament connectivity.

## Reproduction and what remains

Run `python3 05-knowledge/results/integer32_twoedge_20261113/verify_packet.py`
to check the closed artifact manifest, source pins, reviewed promotions,
producer bindings and all paired receipts. Each proof/review lists its full
reproduction command and finite universe. Normal/-O agreement tests active
gates; the independent algorithms and separate reviews are additional evidence.

Two frozen review scripts pin the original candidate document bytes. After
status promotion, replay them with `replay_review.py bridges/check_reconstruction_audit.py`
or `replay_review.py reconstruction/check_final_bridge_audit.py` (add
`--optimized` for the second mode). The launcher reconstructs the reviewed
snapshots in a temporary directory and leaves the live files unchanged.
All mathematical producers still run directly as documented.

The coherent next question is which guarded cuts survive after a proposed
quotient, with their ordered markings and integer weights. The code image is
closed, but arbitrary cuts need not be legal native moves. Complete spectrum
classification, all-order black-core connectivity and compression to an
efficient recurrence using only earlier A000568 values remain open here.
An inherited scalar recurrence with separately computed connected inventory
coefficients already exists; this packet does not claim otherwise.
