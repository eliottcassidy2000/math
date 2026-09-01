# The `D=0` square face: critical-value collapse and the missing source-row invoice

**NON-CANONICAL SESSION SYNTHESIS.** Proved claims below route to THM-4299;
the source-row discussion is a research invoice, not a seam-entry theorem.

## Trigger and inheritance

The incoming frontier had closed exact weight twelve on `U*Z*D!=0` and made
the three wall coordinates visible in source-normal `G` rows, but it left
`D=0` and actual row existence open. The inheritance pass selected:

- THM-4297's all-component good-differential extinction as the closest
  mechanism;
- `(U,W,Z)=(1,-2,1)` as the hostile where both `D` and `Lambda` vanish;
- MISTAKE-531 as the warning against transporting simple-root response data;
- the unspent sidecar `r^2(ar+b)`, the actual double-edge-root coordinate.

The meta-pattern chosen was **divide exceptional multiplicity before judging
a wall**. Here it means dividing by `r-1` only after proving it is a unit, then
retaining the critical value rather than the raw quadratic coefficients.

## Anchor / Niche / Wildcard

**Anchor.** Determine whether the open `D=0`, `UZ!=0` wall supports an actual
degree-carrying special component. THM-4299 proves that the off-corner stratum
`D=0`, `U*Z*Lambda!=0` is empty of nonautomorphic Keller candidates.

**Niche.** Audit whether THM-4298's visible rows already exist in the fixed
THM-3992 gauge. They do not. THM-3997 reaches the first positive diagonals,
THM-4007 reaches `G[t^4]`, and the proposed next lane THM-4020 is explicitly
superseded/unproved. The absolute eliminations through `G[t^5]`, `G[t^6]`,
`G[t^7]`, and `G[t^8]` remain an invoice. In particular, existing proved row
equations force or forbid none of `U=0`, `Z=0`, or `D=0`.

The useful scratch identity is that a new symplectic row cancels from the same
order of the target invariant. If `theta_n` is its remaining tangent
polynomial and `q=gamma*x^2+lambda`, its effect on the next invariant row is

```text
delta G_(n+1)=q' theta_n-q theta_n'/(n+1).                (S1)
```

This identifies a codimension-one horizontal compatibility at each step, but
does not supply the missing absolute row. At rows `6,7,8`, the weight-twelve
channels mix with lower weights `6,8,10`; the diagonal visibility transform
has deliberately forgotten that horizontal sidecar. The cheapest honest next
test is therefore the absolute `G[t^5]` elimination, not another inversion of
the already-proved triangular observer.

**Wildcard.** Treat `D=0` as rank one before doing any tail calculation. Over
`C`,

```text
(U,W,Z)=(a^2,2ab,b^2),
1-UP^6-WS^2P^5-ZS^4P^4
 =(1-aP^3-bS^2P^2)(1+aP^3+bS^2P^2).                     (S2)
```

This wildcard overtook the anchor. It replaced one genus-seven curve by two
`j=0` elliptics and exposed a single order-six contact. The lower rows then
collapsed, by parameterized Morse reduction, to the one-variable critical
series `kappa(t)`. That was the decisive compression because it retained the
first splitter rather than merely observing the vanished discriminant.

## Concept-board update

The live board was

```text
square face | deck character | A11 critical value |
good differential | degree conservation | cubic corner.
```

- The square face changed the deck lane: `tau` exchanges the two elliptics,
  while `tau^2` has the wrong target character for a nonconstant component
  map.
- It changed the differential lane: order nine is checked separately on two
  component generic points, not imported from a smooth genus-seven fibre.
- The order-six contact changed the tail lane: every lower coefficient enters
  only through `kappa(t)`, so the repeated-critical ladder is one-dimensional.
- The critical series changed the genus lane: `m=ord(kappa)` gives the exact
  tail genus `floor((11-m)/2)` and the exact positive form order `9s+5beta`.
- Degree conservation converts those local statements into the actual
  exclusion only after the full component inventory is named.
- At `Lambda=0`, `r-1` ceases to be a unit and the entire Morse connection
  breaks at its source. The first face is cubic, not a limiting quadratic.

## Connection and loss ledger

The exact map is

```text
(r,t,z; all lower rows)
  -> (w,t,z; kappa(t))
  -> first Newton face
  -> normalized tail
  -> vertical order of phi^*eta_0.                       (S3)
```

It preserves divisorial valuations, tail genus, and vanishing of the good
target differential. It destroys the individual lower-row labels and their
global source-normal provenance. `kappa(t)` restores precisely the local
information needed for extinction; it does not restore seam entry. This is
why THM-4299 can close an already-entered wall stratum while the niche invoice
through `G[t^8]` remains fully open.

## Positive and hostile controls

- Positive tail: `(U,W,Z)=(1,2,1)` with one weight-eleven perturbation has
  `kappa=t-4t^2+...`. It produces a genuine genus-five tail, but its good-form
  order is `104`, so it is Keller-constant.
- Failure boundary: `(1,-2,1)` has `D=Lambda=0`. The three branches meet
  pairwise to order six and the first local face is
  `q^3-qz^12-t^12/2` up to units. Division by `r-1` is illegal.
- Endpoint controls: `U=0` or `Z=0` makes the main face split into rational
  factors but also deletes a Newton-hull endpoint. A rational central fibre
  is not a component-inventory proof.

## Status and next tests

- **PROVED:** THM-4299 closes `D=0`, `U*Z*Lambda!=0` inside the inherited
  exact-`M=12` seam.
- **OPEN:** the cubic corner `D=Lambda=0`, the endpoint walls `U=0` and
  `Z=0`, absolute source-normal rows through `G[t^8]`, seam entry, and
  `JC(2)`.
- **Cheapest anchor test:** normalize the cubic face at `(1,-2,1)` and compute
  the first face on which the good differential can have order zero.
- **Cheapest niche test:** derive the absolute `G[t^5]` compatibility in the
  fixed gauge, retaining every horizontal lower-weight coordinate.

The stopping certificate is typed: the quadratic wall has been exhausted;
the surviving corner loses exactly the unit that enabled the Morse quotient.
