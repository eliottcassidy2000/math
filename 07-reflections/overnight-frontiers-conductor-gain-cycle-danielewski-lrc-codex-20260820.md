# Overnight frontier synthesis: conductor cycles, gain cycles, and a six-piece Darboux floor

**Session status (2026-08-20).**  `JC(2)` and `LRC(14)` remain open.  This
session produced and independently audited exact stopping theorems, a sharper
counterexample design box, and one reusable cycle-cohomology lemma.  It did
not construct a planar Jacobian counterexample or remove an LRC row.

## 1. Portfolio and inheritance

The live concept board was:

1. the smooth Danielewski completion of the rational Keller triple collision;
2. the constant-different nongraph normalization and its conductor nodes;
3. target-graph cubic resonance at `a=infinity`;
4. strong-dephasing gain cycles beyond the resistor shadow;
5. the fixed-prime LRC relative-response cospan.

Closest proved mechanisms were THM-3561, THM-2696/3563, THM-2473/3564,
and THM-3534.  The canonical hostiles were the rational Darboux pair
`{b,-1/c}=1`, the nonclosed target minor-Bezout form on the nongraph slice,
the conductance-identical rank-one/rank-two matrix pair, and THM-3534's pure
endpoint class.  Corrected near misses included MISTAKE-420's coefficient-
system warning, the false implication that an area-form lift is integrable,
and the temptation to read a resonant valuation tie as a factor.  The least-
used sidecars were arm valuations, conductor equalizers, and the higher-order
slow generator.

## 2. Anchor: the Danielewski near-counterexample now needs six pieces

[THM-3561](../01-canon/theorems/THM-3561-rational-keller-danielewski-polynomial-completion.md)
gives an everywhere-etale degree-three map

```text
Phi:A2 -> Y,                  Y: c^2 e=b(b+4),             (1)
```

with a literal triple collision.  If `Y` admitted polynomial functions
`P,Q` with `{P,Q}=1`, then `(P,Q)oPhi` would be a noninjective planar Keller
map.  This is the session's sharpest positive counterexample architecture:
the collision and etaleness are already present; only a polynomial Darboux
projection is missing.

[THM-3569](../01-canon/theorems/THM-3569-danielewski-two-by-three-weight-darboux-nonentry.md)
closes the complete two-by-three weight cell.  After constants are removed,
every solution must have support sizes

```text
(2,>=4),                    (>=4,2),                    or (>=3,>=3). (2)
```

The mechanism is not a bounded coefficient search.  Two complementary
weight pairs force the third weight either to be inert, recreating the
forbidden two-by-two rectangle, or to extend it arithmetically.  Arm orders
reduce the latter to three common-power ladders.  Their scalar equations
factor through one of

```text
h(hK)',       h'K+2hK',       3h'K+2hK',       h'J+2hJ', (3)
```

and `b(b+4)|h`, so none can be a unit.  The first open coefficient polygons
are now `3 x 3` and `2 x 4`, both with six pieces.

This is meaningful counterexample squeezing without pretending that a
six-piece search settles the problem.  The rational hostile `-1/c` shows
exactly what polynomiality removes: a one-pole Darboux coordinate already
exists.

## 3. Niche: a conductor cycle blocks every projection on one nongraph slice

[THM-3563](../01-canon/theorems/THM-3563-nongraph-conductor-node-cycle-keller-obstruction.md)
closes every polynomial target projection of THM-2696's sharp
constant-different normalization plane.  A hypothetical constant Jacobian
gives a polynomial primitive `phi`.  On the Laurent double curve its target
data are even, while the area term forces

```text
phi(gamma(h))-phi(gamma(-h))=16 kappa c^2 h^(-3).         (4)
```

Six normalization parameters form three source nodes.  Their branch values
give a triangle with constant edge jump `r!=0`, so descent asks for vertex
potentials satisfying

```text
B_Gamma (a,b,e)^T=(r,r,r)^T,
B_Gamma=[[-1,1,0],[0,-1,1],[1,0,-1]].                    (5)
```

The cycle vector `(1,1,1)` turns `(5)` into `3r=0`, impossible in
characteristic zero.

The proof actually excludes every regular target one-form `eta` with

```text
d(nu^* eta)=kappa du wedge dv,             kappa!=0,     (6)
```

not only decomposable forms `P dQ`.  This is the correct strengthening of
the theorem's nonclosed two-form hostile: a target two-form can pull back to
area, but no regular target primitive can do so.

### Finite node-descent lemma

Let `sigma` pair normalization branches of the same source node and `tau`
pair branches exchanged by the target involution.  The quotient graph has
vertices `S/sigma` and edges `S/tau`.  A prescribed odd jump `j(tau s)=-j(s)`
is realized by a `sigma`-invariant branch-value function exactly when its
edge cochain is a coboundary, equivalently when every signed cycle sum is
zero.  For distinct rational nodes on an affine nodal curve, Chinese
remainders make this finite condition sufficient for arbitrary node values.
It is not sufficient for the differential identity `(6)`.

If also `j(sigma s)=-j(s)`, the permutation `tau sigma` gives directed cycle
components.  A constant jump `r_C` on a component of length `m_C` descends
exactly when

```text
m_C r_C=0.                                               (7)
```

Thus adding more nodes with the same double-odd jump does not help in
characteristic zero.  Viable designs must cut the cycle and export endpoint
debt to infinity, vary the edge amplitudes so every cycle period cancels, or
retain the Kummer owner as a genuine coefficient object.

## 4. Wildcard: dephasing remembers rank one order later

The exact slow-manifold expansion in
[the gain-cycle reflection](strong-dephasing-gain-cycle-and-jacobian-wedge-rank-codex-20260820.md)
is

```text
G_lambda=lambda^(-1)G1+lambda^(-2)G2+lambda^(-3)G3+...,
G1=P K^2 P,
G2=P K^3 P,
G3=P K^4 P-2(P K^2 P)^2.                                (8)
```

Here `K=-i[H,-]` and `P` projects to populations.  `G1=2L_c` is the resistor
network in the user's original dictionary.  `G2` is the triangle sine-flux
circulation.  `G3` contains square cosine flux and is not generally a Markov
generator by itself; a unit `pi`-flux square has a negative nonedge entry.

For the bipartite Hermitian dilation of a coefficient matrix `A`, the
triangle term vanishes and

```text
(G3)_(uv)=6 |(A* A)_(uv)|^2
          -8 sum_w |A_(wu) A_(wv)|^2.                   (9)
```

Conductances plus `G3` therefore recover every column-wedge energy and
`||wedge^2 A||_F^2`.  The leading resistor shadow cannot distinguish a
rank-one matrix from the same-magnitude rank-two hostile; the square-order
correction can.  This is a precise repair of lost phase/rank information,
not a proof that dephasing solves Keller coefficient cancellation.

The literature boundary is consistent: higher-order adiabatic elimination
need not preserve complete positivity.  Hence `(8)` should be read as an
asymptotic invariant-graph generator, not termwise stochastic dynamics.

## 5. LRC anchor audit: exact static cohomology, no row movement

[THM-3534](../01-canon/theorems/THM-3534-r5-middle-response-relative-cospan-and-twisted-h1-collapse.md)
is now promoted and independently reproduced.  In its fixed finite good-
reduction bank, the rank-five defect splits `2+3`, the two response planes
agree only relative to an endpoint line, and formal twisted `H^1` has one
relative class.  Retaining the endpoint adds a second, pure-endpoint class.

This changes no LRC ledger entry: all `165` rows remain.  The `84` objects
are internal cells of one canonical non-cover row, not 84 rows.

The conductor triangle and THM-3534 use the same incidence algebra but not
the same coefficients.  THM-3563 has ordinary scalar graph cohomology; its
oddness belongs to the jump cochain.  THM-3534 has genuine chamber monodromy
on a response representation.  Installing sign transport in the conductor
triangle would solve a different problem and erase the scalar-polynomial
descent predicate.  This is the concrete MISTAKE-420 boundary.

## 6. Other exact squeezes recovered during the session

- THM-3564 restricts factor-bearing target graphs to
  `deg_a(phi)=1 mod 3`; resonance is necessary, not sufficient.
- THM-3565 classifies the reducible first-resonance family, and the proved
  THM-3568 file with slug
  `reducible-target-graph-component-euler-no-go` excludes both components
  from being `A2`.  Irreducible resonant pullbacks remain a separate gate.
- THM-3562 closes every balanced pole passport in its normalized nonsplit
  Faber chart by a common pole unit and a Lagrange sum contradiction.
- RESERVED THM-3567 is auditing whether the separated rational family
  `P=f(x),Q=y/f'(x)` has a nodal, not smooth, full-field polynomial-observable
  completion.  Until promotion, its collision/contraction packet is a proof
  candidate and is not used above.

These results sharpen architectures; none is a global planar degree
exclusion.  The current cited degree passport still leaves the leading cell
`(72,108)` as the first low-height candidate in the repository's reduced
framework.

## 7. Cheapest next experiments

1. **Danielewski `3 x 3`:** enumerate only support sums with at least two
   contributors, then solve the arm-leading equations before coefficients.
2. **Danielewski `2 x 4`:** classify the two extra-weight collision patterns;
   isolated pieces should collapse to THM-3569, while two-sided arithmetic
   extensions are the genuine survivors.
3. **Nongraph escape:** deform the node matching so its jump cochain has zero
   period, then test the stronger one-form gate `(6)` rather than minors.
4. **Gain sidecar:** add `G3` wedge energies to coefficient-fibre searches;
   do not use the sign-indefinite correction as a transition-rate matrix.
5. **LRC:** seek an actual arrival/current map into the promoted response
   carrier; another static quotient computation cannot move the ledger.

The common lesson is now exact.  Squaring amplitudes, quotienting branch
owners, or retaining only static response ranks erases the very cycle class
that later blocks construction.  The productive move is not to avoid the
quotient, but to compute the first correction or sidecar that restores its
lost cycle data.
