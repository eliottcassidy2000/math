---
id: THM-4078
title: "Even-graph triangle quotient spectrum and Boolean noncommutation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT, with one separately
  marked historical FINITE-EXACT sidecar, whose cumulative first-gap
  conjecture was subsequently proved by THM-4083/4416/4427. Fourier orbit sums indexed by
  unlabeled signed switching classes form a complete joint eigenbasis for all
  multiplicity-weighted cycle operators. The triangle quotient Laplacian has
  exact gap 2(n-2), uniquely attained by the single-negative-edge switching
  orbit; its normalized gap is 12/(n(n-1)), while the exact lazy relaxation
  time is n(n-1)/6. The quotient is bipartite. Loop-deleted Boolean supports
  B_3 and B_4 fail to commute for every n>=4. The D=2 gap is proved for all
  n. The original D>=3 census remains a finite audit, while later canon
  closes all cumulative first gaps. Booleanized spectra remain OPEN.
source: codex-frontier-synthesis-creative-20260825d / even-graph spectral lane
audit: >
  PASS. The primary primal path constructs every Eulerian graph through n=6,
  performs 745,164 direct relabel gates, builds the quotient M_3 and M_4
  matrices, checks detailed balance and exact characteristic polynomials, and
  recovers the (0,1) Boolean commutator hostile for n=4,5,6. The independent
  signed-dual path gauges root edges, checks dual orbit counts 2,3,7,16,54
  through n=7, generates 9,432 exact cycle vectors, and uses integer Walsh
  transforms on 2,131,018 characters through n=8. It verifies the proved
  triangle equality orbit and records the cumulative conjecture as
  FINITE-EXACT only. Normal and optimized outputs byte-match; both scripts
  have zero assert nodes and zero floating literals.
depends_on:
  - THM-4073-even-graph-diameter-layer-exact-cycle-distance
related:
  - THM-4069-even-graph-basis-dependence-and-canonical-cycle-envelope
script: 04-computation/even_graph_triangle_quotient_spectrum_thm4078.py
output: 05-knowledge/results/even_graph_triangle_quotient_spectrum_thm4078.out
script_sha256: 85fc39ed2f15d7807aebcfcfc9d927c1b6055becd2d2470a95e34a83805cc7f0
output_sha256: 7456749c0212566ebe56fe0bad5d08b1df73ee42515864221794d064e8a3e648
independent_audit_script: 04-computation/even_graph_triangle_quotient_spectrum_thm4078_independent_audit.py
independent_audit_output: 05-knowledge/results/even_graph_triangle_quotient_spectrum_thm4078_independent_audit.out
independent_audit_script_sha256: d0ea375469c4f3f3b53679c71da506296d8de55f5fd5b3a2cec706151975e082
independent_audit_output_sha256: e9e44ff144433061c4020c5e89f3ae26948b960d2e375be1003924071b3e7c36
hash_basis: raw LF bytes
---

# THM-4078 -- the signed dual diagonalizes the even-graph quotient

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.** THM-4073
retains cycle multiplicities to recover a commuting quotient algebra. The
least-used sidecar there is the dual cut quotient. It supplies a complete
eigenbasis, an exact triangle spectral gap, the correct periodicity hostile,
and the sharp all-order boundary for Boolean noncommutation.

## 1. Complete joint quotient eigenbasis

Let

```text
Z_n={F subset E(K_n): every vertex has even degree}            (1)
```

be the binary cycle space, and let `M_k` be THM-4073's multiplicity-weighted
operator on Eulerian-graph isomorphism classes:

```text
(M_k f)([F])=sum_(C in C_(n,k)) f([F xor C]).                  (2)
```

The character group of `(1)` is

```text
Zhat_n = F_2^E(K_n) / Cut(K_n).                               (3)
```

Represent a coset by the negative-edge set `H` of a signed complete graph.
Adding a cut is vertex switching and does not change

```text
chi_H(F)=(-1)^|H intersect F|.                                (4)
```

Thus the dual objects are signed switching classes. The symmetric group
relabels them. For each `S_n`-orbit `Omega` of switching classes, define

```text
Psi_Omega(F)=sum_([H] in Omega) chi_H(F).                      (5)
```

The functions `(5)`, one per **unlabeled signed switching class**, form a
complete simultaneous eigenbasis for every `M_k`. Indeed Fourier characters
form a basis of all functions on `Z_n`, relabeling permutes them, and the
Fourier transform is `S_n`-equivariant. Therefore orbit sums form a basis of
the invariant function space, which is exactly the quotient-function space.

If

```text
N_(n,k)=|C_(n,k)|=n!/(2k(n-k)!),
c_k^-(H)=#{C in C_(n,k): |H intersect C| is odd},              (6)
```

then convolution gives

```text
boxed: M_k Psi_Omega=lambda_k(H)Psi_Omega,
lambda_k(H)=N_(n,k)-2c_k^-(H).                                (7)
```

The value is constant on `Omega`. More generally every real weighted lift
`sum_k w_k M_k` has eigenvalues

```text
sum_k w_k N_(n,k)-2 sum_k w_k c_k^-(H),                       (8)
```

one entry for every unlabeled switching class, with repeated numeric values
retained at their true multiplicities.

## 2. Exact triangle spectrum and equality orbit

Let `T(H)` be the set of negative triangles of `H`. Equation `(7)` gives the
entire quotient spectrum

```text
Spec(binomial(n,3) I-M_3)
 ={2|T(H)| : [H] an unlabeled signed switching class},         (9)
```

with one entry per switching-class orbit.

A signing with every triangle positive is a cut: fix a root and recover each
edge sign as the sum of its two root-edge signs. Hence every nontrivial
switching class has a negative triangle, and in fact

```text
|T(H)|>=n-2.                                                   (10)
```

Choose one negative triangle `abc`. For each other vertex `x`, the parities
of the four triangles in `K_{a,b,c,x}` sum to zero because every edge occurs
twice. Hence an odd number of `abx,acx,bcx` are negative, supplying at least
one new triangle for each `x`. This proves `(10)`.

Equality is rigid. If there are only `n-2` negative triangles, exactly one of
`abx,acx,bcx` is negative for every `x`, and there are no other negative
triangles. The cases `n=3,4` are immediate. For `n>=5`, four-set parity on
`abxy`, `acxy`, and `bcxy` forces every two outside vertices `x,y` to choose
the same pair, say `ab`. Therefore

```text
T(H)={triangles containing the fixed edge ab}.                 (11)
```

Every triangle of `H xor {ab}` is positive, so the preceding root argument
makes it a cut. Consequently `H` is switching-equivalent to one negative
edge. The converse is immediate.

For a positive-semidefinite Laplacian `L`, write `gap(L)` for its least
positive eigenvalue. For a stochastic operator `P`, write
`gap(P)=1-lambda_2(P)` for the **algebraic** second eigenvalue, not the
absolute spectral radius. Equations `(9)`--`(11)` prove

```text
boxed: gap(binomial(n,3) I-M_3)=2(n-2),                       (12)
boxed: gap(M_3/binomial(n,3))=12/(n(n-1)).                    (13)
```

The gap eigenspace is simple on the quotient. In the full labelled Cayley
graph its multiplicity is `binomial(n,2)` for `n>=4` and one for `n=3`.

## 3. Bipartiteness and the honest lazy gap

The all-negative signing has

```text
chi_H(F)=(-1)^|E(F)|,       lambda_3(H)=-binomial(n,3).        (14)
```

Thus both the labelled triangle Cayley graph and its isomorphism quotient are
bipartite. In particular the nonlazy algebraic gap `(13)` is not by itself a
total-variation mixing statement.

For the lazy quotient walk

```text
P_lazy=(I+M_3/binomial(n,3))/2,                               (15)
```

the exact gap and relaxation time are

```text
boxed: gap(P_lazy)=6/(n(n-1)),
boxed: t_rel=n(n-1)/6.                                       (16)
```

The walk is reversible with stationary mass proportional to the number of
labelled Eulerian graphs in each isomorphism class, by THM-4073's detailed-
balance identity.

## 4. Boolean supports fail at every possible order

Let `B_k` be the loop-deleted Boolean support of `M_k`. For every `n>=4`,

```text
(B_3 B_4)([empty],[C_3])=0.                                  (17)
```

The only `B_3` neighbor of `[empty]` is `[C_3]`, and the `B_4` loop at that
class is deleted. In the opposite order,

```text
(B_4 B_3)([empty],[C_3])=1.                                  (18)
```

Indeed the only `B_4` neighbor of `[empty]` is `[C_4]`, while the labelled
witness

```text
(12,23,34,41) xor (12,23,31)=(31,34,41)                      (19)
```

shows that `[C_4]` is a `B_3` neighbor of `[C_3]`. Isolated vertices extend
the witness to every larger `n`. Therefore

```text
boxed: B_3 B_4 != B_4 B_3 for every n>=4.                     (20)
```

This sharpens THM-4073's order-four counterexample to the full possible
range. At `n=3`, `B_4` does not exist.

## 5. Proved spectral reduction; historical cumulative conjecture and its continuation

**Current continuation:** the cumulative first-gap question below is now
proved: [THM-4083, cumulative D3/D4 gaps](THM-4083-even-graph-cumulative-d3-d4-spectral-gap.md),
[THM-4416, cumulative D5/D6 gaps](THM-4416-even-graph-cumulative-d5-d6-spectral-gap.md),
and [THM-4427, all remaining cumulative gaps](THM-4427-all-cumulative-signed-cycle-gaps-by-transposition-rigidity.md).
Their stated small-order equality exceptions remain in force. The second
weights for n>=16 are treated by
[THM-4433, signed-cycle second minima and cross-scale stability](THM-4433-signed-hamilton-second-minimum-and-cross-scale-stability.md).
These are continuations, not circular dependencies of this original
Fourier reduction. Booleanization is still outside their spectral claims.

For every `2<=D<=n-1`, the proved equations `(7)`--`(8)` reduce the natural
weighted diameter-layer Laplacian gap to

```text
2 min_([H] nontrivial) sum_(k=3)^(D+1) c_k^-(H).              (21)
```

For `D=2`, Section 2 proves the minimum in `(21)`. **Historical scope of
the original 2026-08-25 proof:** the following D>=3 extension was then
conjectural, with this theorem's census **FINITE-EXACT THROUGH `n=8` ONLY**.
That exact dual census supports

```text
min_([H] nontrivial) sum_(k=3)^(D+1)c_k^-(H)
 =sum_(k=3)^(D+1) (n-2)!/(n-k)!,                              (22)
```

the number of such cycles containing one fixed negative edge. Equality is the
single-edge switching orbit in the verified box, except at `(n,D)=(4,3)`,
where the antibalanced orbit also ties.

The suggested termwise lower bound fails, so it cannot prove `(22)`. The
all-negative signing has

```text
c_k^-(H)=N_(n,k) for odd k,        c_k^-(H)=0 for even k.     (23)
```

Thus the tempting per-length lower bound already fails completely at `k=4`;
only the cumulative prefix statement survives this hostile. The original
finite census does not prove an all-`n` claim; the continuation cited above
supplies the later proofs.

## 6. Exact audits and loss ledger

The primal audit constructs `Z_n` from anchored triangles, canonicalizes all
labelled states directly under every permutation through `n=6`, builds `M_3`
and `M_4`, and computes exact integer characteristic polynomials. The
independent dual audit gauges all root-incident signs positive, generates
simple cycles from cyclic orders, applies integer Walsh transforms through all
`2^21` characters at `n=8`, and separately recovers the dual orbit counts
through `n=7` using adjacent transpositions.

Replay from the repository root:

```bash
python3 -B 04-computation/even_graph_triangle_quotient_spectrum_thm4078.py
python3 -B -O 04-computation/even_graph_triangle_quotient_spectrum_thm4078.py
python3 -B 04-computation/even_graph_triangle_quotient_spectrum_thm4078_independent_audit.py
python3 -B -O 04-computation/even_graph_triangle_quotient_spectrum_thm4078_independent_audit.py
```

The source is THM-4073's multiplicity-weighted commuting lift. The map is the
Fourier transform followed by the switching and relabeling quotients. It
preserves the full negative-cycle vector `(c_3^-,...,c_n^-)`; Booleanization
destroys multiplicities, and `(17)`--`(20)` show that commutativity is lost.
This theorem's original proof establishes the complete spectral reduction
and its triangle specialization. The later cited continuation proves `(22)`;
no Boolean spectral equality follows from either proof.
