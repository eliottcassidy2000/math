# Planar Jacobian: owner descent, cyclic projectors, and the torsion-envelope frontier

**Current status (2026-08-26).** This session produced two canonical results:
[THM-4248](../01-canon/theorems/THM-4248-weight-eleven-z-zero-owner-descent-planar-jacobian-exclusion.md)
and
[THM-4249](../01-canon/theorems/THM-4249-w0-cyclic-projector-missing-eigenline-attachment-squeeze.md).
They are local advances inside the inherited reduced `(2,3)` seam. They do
not prove seam entry or `JC(2)`.

## Portfolio and inheritance

**Anchor — exact `M=11`.** The closest mechanism was THM-4222/4232's
primitive-CM degree-zero specialization. The hostile was the hidden elliptic
side face from THM-4218. The corrected near miss was the incomplete-support
failure recorded in MISTAKE-487/522. The least-used sidecar was the affine
relation

```text
K=2848/45-(7/6)Delta,
```

which makes the first-surviving-owner split exhaustive and prevents a false
`Delta=K=0` stratum.

**Niche — exact `M=12`, `W=0`.** The closest mechanism was THM-4230's finite
marked-ratio reduction, with THM-4247's concurrent degree-twelve involution
exclusion as an independent same-geometry control. The hostile was THM-4241's index-four
visible-hidden glue, which invalidates a rational-eigenspace sieve. The
corrected near miss was MISTAKE-521: rational isogeny data is not the integral
Hom lattice. The least-used sidecar was the integral `O[C_12]` action.

**Wildcard — incidence geometry rather than another shell enumeration.** The
torsion condition turns the residual into a finite bipartite incidence graph
between map orbits and ratio orbits. This is not a tournament: the observable
is symmetric incidence, has many ties/nonedges, and carries no intrinsic
orientation. Forcing it into a tournament would discard the annihilator ideal,
which is exactly the coordinate that makes the workload small.

## Live concept board

The board ended with six concepts:

1. first surviving Newton owner;
2. positive-genus simplicity versus elliptic targets;
3. integral cyclic projectors;
4. low hidden-shell cyclic determinants;
5. CM annihilator ideals and target torsion;
6. attachment-orbit/root-choice symmetry.

The important interactions were:

- The owner split does more than classify faces: it selects the base exponent
  and hence the exact regular model. The relation between `Delta` and `K`
  closes the apparent missing corner.
- A rational projector would lose the index-four glue. Multiplying the
  spectral factors inside `O[T]` instead produces integral maps to `a_u u`,
  `d v`, and `2ell`; the unavoidable factors of two become useful torsion.
- The hidden determinant is not merely a lattice index. At degrees `12` and
  `24`, it converts collapse of one cyclic vector and its translate into
  simultaneous two-torsion for the explicit basis maps `f,g`.
- Once the `u` coordinate vanishes, the `v` coordinate turns marked ratios
  into a quotient of CM torsion: `R=-1/(X^3+1)`. The target unit quotient is
  therefore the native representation, not an after-the-fact symmetry cut.
- Retaining the actual ideal `(d)` after taking the global torsion union
  changes the next problem from `308*(55+35)` naive pairs to exactly `1,512`
  map-ratio incidences.
- The common ratio `1/3` exposes a second integral denominator: reduction
  modulo the ramified Eisenstein prime `pi` turns the cyclic determinant into
  the anisotropic form `A^2+B^2` over `F_3`. Divisibility of that determinant
  would divide the hidden projector by `pi`, landing at a degree forbidden by
  the hidden Gram. This prime-local obstruction interacts with every one of
  the 132 degree-42 map orbits at once.

## Anchor result: why the `Z=0` wall closes

The complete support admits exactly five owner strata, according to the
first surviving pure-`p` owner and whether the `K` tail is present. Two
independent exact hull engines check all `26,624` optional-support and
collision patterns, including the fixed residual point `(2,0,1)`.

Every stratum has one genus-five main component

```text
C: 1-A S P^5-B S^3P^4=0
```

and only rational side components. Its cyclic `11`-cover has primitive CM
type `{4,5,8,9,10}`, so the cited primitive-CM theorem makes `J(C)` simple.
Dimension then forbids a nonzero map to the good elliptic target. Resolving
the generic rational map adds only rational exceptional curves, so all
special components are constant and degree conservation contradicts a
positive-degree generic Keller morphism.

The equality/failure boundary is geometric, not bookkeeping:

- `A=0` replaces the main face;
- `B=0` creates a genus-four/elliptic hostile; and
- `A+B=0` coalesces eleven transverse contacts into a double edge.

Together with THM-4232, this closes exact `M=11` for arbitrary `U,Z` whenever
`A*B*(A+B)!=0`. Exactly those three walls remain.

## Niche result: why the missing eigenline is forced

For `m=a_u u+b f+c g+d h`, the integral cyclic projectors give

```text
P_u(m)=a_u u,       P_v(m)=d v,       P_L(m)=2ell,
q(m)=4N(a_u)+N(d)+3K,                 q(2ell)=12K.
```

Attachment collapse is the vanishing of all adjacent difference classes, so
it survives every polynomial in the source action. The visible maps force
`N(a_u)>=7` and, when `d!=0`, `N(d)>=3`. The only nonzero-`a_u` degree
profiles then have `K=1` or `2`.

Those two hidden shells collapse to the same hostile. Exact cyclic
determinants force `[2]f(Q_0)=[2]g(Q_0)=O`. The explicit `Y`-numerators are,
up to units,

```text
y(t^2+rho^3),             y(1+rho^3t^2),
```

where `rho^4-2rho^3-2rho+1=0`. The gate makes both points affine and `y`
nonzero, so simultaneous two-torsion would give `rho^6=1`. The resultant is
`-108`, hence impossible. Therefore `a_u=0` for every degree-`34/42`
candidate.

THM-4247, promoted concurrently, proves the degree-twelve half of this
statement by an orthogonal method: its hidden denominator has degree at most
three and factors as `t(t^2-1)`, whose reciprocal intersection lies only on
the excluded `Z=0` wall. Agreement between that pole-divisor route and the
cyclic-determinant/resultant route is a genuine hostile control.

The exact residual is:

```text
degree 34: 4224 vectors / 176 source-target orbits,
degree 42: 3168 vectors / 132 source-target orbits.
```

The remaining `d` coordinates imply `[d(omega^2-1)]v(Q_0)=O`. Ideal nesting
then gives torsion unions of `336` and `216` points. After target-unit
quotient and removal of the three gate-boundary orbits, the raw envelopes
contain exactly `55` and `35` marked ratios.

The common degree-42 `E[3]` orbit gives `R=1/3` for all 132 residual map
orbits. Write the hidden projector as `H=A f+B g`, with

```text
q(H)=12K,       K in {5,7,10,11,13},
Delta=A^2+omega B^2.
```

Collapse forces both `H(Q_0)` and `(TH)(Q_0)` to vanish, hence
`[Delta]f(Q_0)=0`. If the norm-three prime `pi=omega^2-1` divided `Delta`,
then the anisotropy of `A^2+B^2` over `F_3` would force `pi|A,B`. The divided
hidden vector would have degree `4K`, impossible because all hidden degrees
are divisible by six and no listed `K` is divisible by three. Therefore
`3` does not divide `N(Delta)`.

At the exact good specialization to `F_397`, the point `f(Q_0)=(340,181)`
has additive order 18. Any characteristic-zero collapse specializes to
`[Delta]f(Q_0)=0`; multiplying by the conjugate would make `N(Delta)` kill
an order-18 point, contradicting the preceding prime-local obstruction. Two
independent exact implementations reconstruct all 3,168 vectors and 132
orbits and retain direct per-orbit nonvanishing as a hostile control. Thus
`1/3` is excluded, the final ratio envelopes have cardinalities `55/34`, and
the exact residual frontier has `864+648=1,512` incidences.

## Typed connection ledger

| source | target | map | preserved predicate | destroyed information | needed sidecar | cheapest decisive test |
|---|---|---|---|---|---|---|
| full Newton support | owner strata | lower-hull normal fan | complete face ownership | lower coefficients inside a face | fixed residual point and `K(Delta)` | exhaustive collision census |
| main face | elliptic response | Jacobian/CM specialization | existence of nonconstant component map | rational paths and graph cycles | proper regular model | simplicity plus resolved degree conservation |
| full Hom lattice | three eigensectors | integral `O[T]` polynomials | attachment-difference vanishing | mixed cancellation between sectors | index-four glue basis | exact projector multiplication |
| hidden shell | two-torsion hostile | cyclic determinant/adjugate | simultaneous attachment collapse | individual target values up to units | explicit `f,g` formulas | quartic-versus-sixth-root resultant |
| residual `d` coordinate | marked ratio | `P=v(Q_0)`, then `X^3` | necessary collapse condition | hidden coordinates of the map | CM ideal and gate orbits | kernel union/intersection census |
| common `E[3]` ratio | all degree-42 map orbits | ramified-prime determinant descent | collapse of `H,TH` and integrality | behavior at other target-torsion primes | hidden-degree divisibility and good reduction | order-18 point at `q=397` |
| ratio envelope | exact attachment question | map-ratio incidence graph | every possible witness | algebraic equality of mixed maps | explicit normalized full basis | good-prime pole-ideal evaluation |

## Hostile probes and stopping reasons

- **Owner continuity failed.** The `K=0` face must be deleted and changes the
  packet; it cannot be recovered by a limit argument. The repaired atlas
  treats it as a separate stratum.
- **Rational spectral projection failed.** THM-4241's mixed degree-four map
  proves the direct-sum lattice has index four. Integral group-ring
  polynomials repair the loss.
- **Degree/gcd alone failed.** Raw shells contain tens of thousands of
  candidates and every residual orbit has size 24. Fibre degree only creates
  the first norm bounds.
- **The torsion envelope is not an equivalence.** A point in `E[d*pi]` makes
  the visible projector compatible, but hidden-visible cancellation can
  still fail. The ratio `1/3` is the first exact false positive; no member of
  `S_34` or `S_42` has been produced.
- **A finite-field zero would prove nothing in the needed direction.** The
  successful certificate instead finds nonvanishing at a prime where the
  algebraic constants, node radicals, source, target, and map denominators
  all have good reduction. Characteristic-zero collapse would survive that
  specialization, so the contrapositive is sound.
- **A direct symbolic addition prototype was too slow.** This is an
  engineering stopping reason after the mathematical universe became finite,
  not evidence that the incidences survive.

## Sharp next fronts

1. **Evaluate the 1,512 incidences.** Use the explicit basis `[u,f,g,h]` over
   a good finite-field embedding, represent each mixed map by pole ideals or
   reduced function-field coordinates, and test equality on one canonical
   twelve-node orbit. A characteristic-zero degree bound must accompany any
   modular nonvanishing claim.
2. **Prime-stratify the remaining ratios.** For each target torsion orbit,
   intersect the prime factors of its annihilator with those of the cyclic
   determinant `Delta=A^2+omega B^2`. Test whether a residue-field norm form
   forces divisibility of both hidden coefficients and hence another
   forbidden divided degree. The cheapest next lane is a ratio shared by the
   largest number of actual `(d)` ideals, so one proof removes many
   incidences rather than one map-ratio pair.
3. **Attack the three `M=11` walls separately.** `A=0` asks for the new
   genus-five CM type, `B=0` for the genus-four elliptic quotient, and
   `A+B=0` for a nontransverse stable model. They are different geometries,
   not one coefficient-deletion problem.
4. **Lift away from `W=0`.** The projector proof uses the enlarged `C_12`
   action and will not persist verbatim. The transferable object is the
   annihilator-ideal sidecar: search other hidden-Hom points for a smaller
   symmetry algebra that still produces integral projectors.
5. **Keep entry separate.** Even closure of all exact-weight cells inside the
   seam would not prove that every planar Keller counterexample enters it.

No meta-pattern was promoted from this single session. The promising reusable
move is “replace rational eigenspace projection by an integral group-ring
annihilator and keep its denominator as torsion,” but it needs evidence from a
second problem before entering `META-PATTERNS.md`.
