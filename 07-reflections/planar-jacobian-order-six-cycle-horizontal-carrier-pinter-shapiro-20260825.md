# Planar Jacobian: the order-six lift hidden in a rational three-cycle

**Status:** research synthesis after
[THM-4139](../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md)
and [THM-4146](../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md).
`JC(2)` remains **OPEN**. THM-4138 has proved the degree-`16/15` `Delta_V`
exclusion, and THM-4143 completes the old exact-`M=8` trichotomy; `M>=9`,
other cells and entry remain open. Literature scope is pinned in the
[Pintér--Shapiro source sidecar](../05-knowledge/reference/CORE-PAPERS-PINTER-SHAPIRO-2026-08-25.md).

## Inheritance pass

- **Closest proved mechanism:** THM-4139 classifies the rational cycle and
  separates its nodal normalization incidence from the horizontal `3P`;
  THM-4138 independently closes the former degree-`16/15` carrier by
  Mordell--Weil, labelled vanishing-loop and orbit-merger arguments.
- **Canonical hostile:** the same source roots `{-6,6,0}` enter three
  different maps: a Möbius root permutation, the `q=0` nodal normalization,
  and THM-4134's polynomial horizontal section. Their target data disagree.
- **Corrected near miss:** THM-4139 refutes equality of the normalization
  points with horizontal `3P`; THM-4146 retains only the weaker equality of
  source root divisors. A `2:3` formula, root set or cycle-type match alone
  transports no Keller predicate.
- **Closest central-sign hostile:** THM-3742 already separates a signed Pell
  `C_14` lift from its projective `C_7` quotient and restores the sign by a
  stereographic sidecar. The least-used adjacent tools are THM-3465's finite-
  character rigidity and THM-4057's rational/Pythagorean carrier.

## Live concept board

| Concept | Exact carrier | Load-bearing sidecar |
|---|---|---|
| closed `Delta_V` carrier | THM-4138 actual BC section/loops | `a,rho,q`, pole pair, target image, meridians |
| rational quadratic cycle | `x^2-29/16` | exact denominator and height |
| projective/linear lift | `PGL_2` order `3`, `SL_2` order `6` | central sign / defining section |
| `3:4:5` tree root | signed leg/hypotenuse cycle | scale and ordered roles |
| Petersen-minor flow | nonzero `F_2^2` palette | incidence and vertex conservation |
| double-well passages | Poisson sector count | parity, defects, uniform tail |
| `63=2^6-1` | derivative degree / matrix scalar | prime-power and sign data |

The productive collision was not “all of these contain a three.” It was the
same typed quotient appearing three ways:

```text
six-state lift --forget central sign--> three-state projective palette,
ordered passages --forget exact count--> even/odd spectral character,
defining cubic --forget line character--> invariant three-line divisor.
```

## Exact gains

### 1. Scalar period six was separated from the genuine six-lift

The entire rational preperiodic graph is now classified. There is no rational
six-cycle. The exact-period-six dynatomic instead has irreducible degree `54`,
eighteen real roots and three real cycles. The canonical rational length-six
object is the lift

```text
A=[[1,-13],[1,3]],       A^3=-64I,
B=A/4 in SL_2(Q),        B^3=-I,        B^6=I.
```

This prevents two common errors: calling the projective three-cycle a hidden
scalar six-cycle, and mistaking `63` for the exact-period-six point count.

### 2. The cycle is a complete Eisenstein hexagon

The six integral lifts are exactly the solutions of

```text
X^2+2XY+13Y^2=48.
```

In `(X+Y,2Y)` coordinates this is an Eisenstein norm. Antipodal quotienting
produces the three slopes `-7,5,-1`. The homogenized cycle cubic satisfies

```text
L(Bz)=-L(z).
```

Thus the divisor is projectively invariant while its defining section is
anti-invariant. This is a literal, executable defining-section character in
the comparison model, not an established normal or meridian character on the
planar boundary carrier.

### 3. `3:4:5` and `29` are forced

For a right triangle `a^2+b^2=h^2`, requiring

```text
-(a+b) -> h -> -(b-a) -> -(a+b)
```

under `(y^2-D)/b` forces `h=3a-b`, hence `4a=3b` and

```text
(a,b,h,D)=(3k,4k,5k,29k^2).
```

The primitive recurrence is precisely `(y^2-29)/4`. Separately, the cycle is
the unique rational quadratic three-cycle whose points form an arithmetic
progression. These are two independent characterizations.

### 4. The THM-4134 collision is real but has a uniform fibre obstruction

Writing `t=4rho/(3a)` turns the horizontal section into

```text
U=a/2+t^2,       V=-t(t^2+3a/4).
```

At `a=-48`, `V=-t(t-6)(t+6)`, and `t=y+1` identifies its zero polynomial
with the cycle cubic. The same `-48` is the characteristic discriminant of
the centered integer order-six matrix. Explicitly, if
`S=[[1,1],[0,1]]`, then `A_t=SAS^-1=[[2,-12],[1,2]]`, and its projective map
`t |-> 2(t-6)/(t+2)` cycles `-6 -> 6 -> 0`. This action is introduced only on
the root divisor; it is not an automorphism or monodromy action of the
horizontal surface or BC cover.

The common-fibre transfer fails uniformly, not accidentally. If
`r^2=-3a/4`, then

```text
t=0:       q=a^3/2,
t=+-r:     q=5a^3/64.
```

Their difference is `27a^3/64`. The outer roots collapse to one full
`(U,V,q)` point; only the `t`/`rho` sheet sign separates them. Therefore this
is a section-polynomial/discriminant collision, not a three-point orbit on one
elliptic fibre. The negative result is useful: it names the first failed
implication and the missing coordinates.

Incoming THM-4139 sharpens the picture. Its zero-fibre normalization sends
the same parameters to `(-12,216),(-12,-216),(-48,0)` at `q=0`, whereas the
horizontal section sends them to the two values above. The Möbius map cycles
only the parameters. Thus there is no conflict between the two results: they
prove source-divisor coincidence and target-map separation, respectively.

### 5. `63` has three typed roles and one hostile

- `deg((g^6)')=63` by the exact chain-rule product of degrees
  `1+2+4+8+16+32`.
- The three integer pair sums multiply to `64`, and `A^3=-64I`.
- `63=3^2*7` has no primitive prime divisor at exponent six; its prime
  support is inherited from exponents two and three, although the prime-power
  and sign sidecars still produce order six modulo `63`.
- Scalar iteration modulo `63` has three cycles, all of length three. Hence
  `63` belongs to derivative-degree and modular order-six/lift structure, not
  to the affine cycle multiplier or a scalar six-cycle.

## What the two incoming papers contributed

### Pintér: finite orbit representatives with full conservation retained

[Pintér v2](https://arxiv.org/abs/2607.22267v2) proves the proper-Petersen-
minor four-flow theorem through a strong architecture: reduce a minimal
obstruction to a rigid core, quotient the exceptional extensions by genuine
automorphisms, and close each representative with an explicit minor model.

The exact local bridge is the nonzero `F_2^2` palette. Modulo two, the standard
order-six matrix cycles its three nonzero vectors. In a cubic graph,
conservation forces every vertex to see all three once, so the palette can be
identified equivariantly with the normalized rational cycle
`{-2,1,-1/2}`. The map is useful only with the cubic incidence/conservation
sidecar. It gives no graph-minor-to-Jacobian implication.

### Shapiro: fixed-complexity closure is not a global closure

[Shapiro v1](https://arxiv.org/abs/2608.23342) proves a Poisson limit for
double-well passage counts by extracting a rank-one kernel and controlling
three defects, then proving a summable bound uniform in the sector number.
The parity of passage number is the exact character separating the two low
spectral sectors.

The planar lesson is methodological and sharp: closing every fixed
monodromy-word length, jet depth or character weight does not close their
union. The analogue of Shapiro's uniform tail would have to be an algebraic
degree/height/ramification majorant. No such majorant is imported.

## Consequences for the planar frontier

### Positive control

The order-six matrix is a determinant-one linear polynomial automorphism.
After real conjugacy it is a sixty-degree rotation, so THM-3465 closes the
real dagger-paired pure-character Keller sector. This calibrates the symmetry
lane but does not supply an entry theorem. Complex cyclic-equivariant maps
can contain nonlinear triangular automorphisms.

### Retrospective obstruction address

The three-point projective divisor alone could not have closed the old
`Delta_V` branch. THM-4138 succeeds because it retains the actual BC section,
labelled meridians and orbit-merger budget. Any future use of this comparison
route in a new cell must retain

```text
(root divisor, model defining-section character, a, t/rho sheet, q fibre,
 actual BC normalization and target map, labelled meridians when constructed).
```

The first coordinate has projective period three and the comparison model's
defining cubic detects the central sign; the fibre and sheet data refute the
common-fibre inference. This is not a normal-bundle character on the BC cover
and is not a new ingredient in the already proved THM-4138 exclusion.

## Generated next tasks

1. **Anchor -- exact-`M=8` closure extraction for `M>=9`.** Isolate which
   THM-4130/4138/4141/4143 inputs are uniform in `M` and which depend on the
   exact-`M=8` packet. The first deliverable is a typed obstruction matrix for
   one smallest `M=9` cell, not an analogy with the closed wall.

2. **Anchor -- three-map separation lemma.** Abstract the source divisor,
   normalization map and horizontal-section map into a commuting-or-failing
   diagram. State necessary and sufficient sidecars for a root permutation to
   descend to a target or cover action; use THM-4139/4146 as the hostile.

3. **Niche -- decorated finite passport certificates.** Following Pintér's
   method, quotient an actual open-cell attachment bank by decorated-dual-graph
   automorphisms, not bare graph isomorphism. Retain pole orders, local degrees,
   intersection numbers, root labels, actual monodromy and any computed normal
   character; keep the model defining-section character as comparison only.

4. **Niche -- separate equivariant calibration sector.** Declare source and
   target actions `B_s,B_t` and the law `F(B_s z)=B_t F(z)` before diagonalizing
   over `Q(sqrt(-3))`. Decompose both components into residues modulo six and
   write every character component of the Poisson bracket. Use THM-2230 to
   remove only genuine target-shear directions in a fixed response fibre.
   Without an entry map this calibrates equivariant Keller maps; it is not
   progress on an open Jacobian cell.

5. **Wildcard -- labelled monodromy transfer operator.** For one new `M>=9`
   cell, first construct the complete labelled meridian generating set. Only
   then form a lazy operator and compute its exact stationary projector and
   defect spectrum. Retain the ordered Hurwitz word and commutator cycle as a
   sidecar; averaging alone destroys the obstruction THM-4130 uses.

6. **Wildcard -- uniform-complexity gate.** Search for an algebraic analogue
   of Shapiro's uniform tail: a bound that controls all word lengths or all
   response weights simultaneously from degree, pole width and ramification.
   Fixed-depth eliminations do not compose without it.

7. **Wildcard -- decorated two-fibre classification.** For
   `V=-t(t^2+3a/4)`, classify affine root-divisor identifications while
   retaining the singleton `q=a^3/2` root, two reduced outer source points in
   the `q=5a^3/64` fibre and both `rho` sheets. Compute the decorated stabilizer
   and test compatibility with possible defining-section characters. A common
   `q` fibre is already
   impossible for every `a!=0`, independent of relabelling or trace shape.

## Stopping reasons

- No actual Keller map or cover action realizes the order-six comparison
  symmetry on the horizontal carrier.
- The three roots do not share one `q` fibre; any claimed elliptic orbit would
  be false.
- Pintér's minor operations do not automatically model resolution moves, and
  a bare dual graph loses Keller data.
- Shapiro's positivity, spectral gap and semiclassical parameter have no
  established algebraic counterpart here.
- THM-4138 and the exact-`M=8` trichotomy are already closed by independent
  proofs. The correct promotion here is the narrower THM-4146 extension and
  new `M>=9` tasks, not a new `JC(2)` consequence.
