# Independent hostile audit: the seven cells carry a genuine mixed finite spectrum

**Status: FINITE-EXACT INDEPENDENT AUDIT ACCEPTS THE SCOPED `7 x 13`
SOURCE-TIME PACKAGE.  LRC(14) remains open.**  The audit reconstructed the
source tensor from the hash-pinned THM-2594 interval engine and rebuilt all
thirteen endpoint residue matrices from the independent THM-3514 atom engine.
The submitted seven-cell probe was neither imported nor read until the
clean-room semantic surface had been fixed.

## Verdict

Retaining `cell_ell(y)` before the first-collision pair integral gives a
literal four-coordinate finite tensor

```text
M(omega,nu,ell,c),
```

not a post-hoc split of the previously marginalized owner profile.  It has
the reported 362 atom pairs and 5150 nonzero entries out of 138411.  Every
one of the thirteen `(1,0,t)` endpoint residue matrices contracts to a
`7 x 13` table with all 91 Fourier coefficients nonzero.  For the fixed
class `(1,0,6)`, output double-centering removes exactly the DC and two axis
sectors and leaves all 72 mixed modes nonzero.

This is a finite source-time simple-kernel residue pushforward.  The endpoint
`AX/BY` factors are still preintegrated scalars, so it is not an integrand-level
physical current, a THM-2512 bridge, an exact-address `C(a;X,m)`, a row
exclusion, or LRC(14).

## Cell construction before marginalization

The seven half-open THM-2594 cells were rebuilt directly:

```text
cell_ell = {y : ||y-ell/7|| < 1/14},   ell in F_7.
```

They are pairwise disjoint, tile the full outer base, and each has exact grid
mass `42548128262640`.  The source endpoint set was split into all 39 guard
atoms before `P^2`; each atom was transferred, restricted by `Q`, and folded
through `P_(13^5)` separately.  Summing the atom profiles restores `P^2 e`,
`U`, and `V` exactly.

For each aligned atom pair and common offset `c`, the product profile was
formed once and integrated over each of the seven cells.  The seven masses
sum to the unsplit product mass entry-by-entry.  Thus summing `ell` recovers
all `39 x 39 x 13` source entries, whose digest is

```text
96337f74f2599044870e90878b8bcfb2ce4a04dc970c88b28f23a4f237b3ea53.
```

The audit stored the last two axes as `[c][ell]`, giving digest
`5dadb2ec...af06`; the candidate used `[ell][c]`.  Explicitly transposing
those axes gives the candidate digest exactly:

```text
39d7a0b4e5b2d8b85631d682ed1967091e44dc41e17b33a77e7184d3dc93e0cf.
```

The initial digest difference was therefore serialization, not mathematics.

## Complete `F_7 x F_13` closure

The independent endpoint reconstruction reproduces the complete pair-bank
digest

```text
c28119c8b54f47e5b7a46f1508fbba604b0e3997eaadb05b03ad28edd9aed468.
```

For each residue `t`, contract this matrix with `M` and transform in the cell
and common-offset coordinates.  Every residue has support census

```text
(total, DC, F_7 axis, F_13 axis, mixed) = (91,1,6,12,72).
```

The coordinate and spectrum banks agree exactly with the candidate:

```text
989dafc220a6d09aeacfce4af0e9a4fe13eedacc79fa66032ea39bc107fd8efb
5f173227c5e203309f61bdfd9d47cc64a3b49ae8f14abd0f7bfc469eda278533.
```

At `(1,0,6)`, the endpoint total remains
`225010624370142818572`, and the first mixed coefficient is
`218019411785559321795`.  The cell marginal also recovers the previously
audited thirteen owner coefficients exactly.

## The two ANOVA operations must not be conflated

An initially plausible stronger conjecture failed cleanly: endpoint-pair
ANOVA components do not route into separate output Fourier sectors.  In
fact, after contraction, each of the flat, left-only, right-only, and
doubly-centred endpoint-pair components has support `91/91` for every
residue.  In particular the fixed class's endpoint-pair interaction has
shape `(91,1,6,12,72)`.

Output ANOVA is different.  Double-center the already contracted `7 x 13`
fixed-class table.  Its transform has shape

```text
(72,0,0,0,72),
```

so all and only the mixed modes remain.  This distinction explains why the
endpoint-pair interaction can be full while output centering is purely mixed;
the two centering operators act on different index sets and do not commute
with the contraction in a sector-separating way.

The coordinate erasures provide exact hostile controls.  Replacing all cell
rows by their average gives `(13,1,0,12,0)`; replacing all owner columns by
their average gives `(7,1,6,0,0)`.  Either operation kills every mixed mode,
while preserving the surviving axis.  Thus both retained coordinates are
load-bearing for the mixed block.

## Typing boundary

The positive result says that a finite pair function of aligned source-time
guard atoms has maximal spectral capacity after a lawful one-base cell
refinement.  It does not move the endpoint integrations themselves onto that
base.  Consequently:

- `(1,0,6)` denotes the THM-2334 residue class, not an isolated exact address;
- the coefficients `AX(omega)` and `BY(nu)` remain endpoint integrals already
  evaluated elsewhere;
- resemblance to the `F_7 x F_13` shape in THM-2512 supplies no bridge theorem;
  and
- the `U_full`/`U_clock` two-transplant temporal boundary remains untouched.

The next obstruction is therefore not spectral closure.  It is an
integrand-level temporal/Fubini identification plus exact-address isolation.

## Reproducibility

Normal and optimized runs are byte-identical.  The audit semantic digest is
`a0534b8d995956c126ce204117a9488b222ad25dd8b483b269e37740ffe13ccb`;
the pinned script LF digest is
`c9a07252d1573cca47ce5a161d8dade2821078fc68f0f8aa03d0fd92763dae5a`.
