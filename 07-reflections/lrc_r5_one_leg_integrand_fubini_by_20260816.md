# The literal BY leg obeys exact product Fubini but does not move to the ancestry stalk

**Status: FINITE-EXACT VERIFIED IN NORMAL AND OPTIMIZED MODES; external
hostile audit pending.**  This note closes the cheapest one-leg equality
test requested by
`lrc-r5-source-aligned-relation-residue-7x13-spectrum-codex-20260816.md`.
It proves an exact external-product Fubini identity with every declared guard
and retained source coordinate.  It does **not** construct a same-stalk
temporal transport, a THM-2449 response table, a THM-2512 bridge, a grouped
exact-address coefficient, a physical current, a row exclusion, or an
LRC(14) result.  LRC(14) remains open.

## Inheritance and the one-leg choice

The closest positive input is the source-aligned tensor

```text
M(omega,nu,ell,c),       ell in F_7, c in F_13,
```

which places the two guard-atom labels on one THM-2471 first-collision
ancestry base before integrating that base.  The closest corrected boundary
is THM-2538: separately integrated endpoint banks do not determine their
mixed transportation coupling, while a genuine source/arrival product must
be formed on one ancestry base before integration.  MISTAKE-281 and
MISTAKE-283 forbid reading matching labels, nonzero Fourier factors, or a
common row as that missing fibre product.  MISTAKE-313 adds the factor-level
warning: the full list of guards and clock factors, rather than an output
hash alone, determines the carrier.

Of the two preintegrated endpoint factors, `BY` is the lawful cheapest leg.
It is a direct endpoint sum on the actual half-open intervals.  `AX` also
contains the transformed `x/Q` overlap sweep and therefore remains frozen in
this test.

## Exact source, target, and map

For every endpoint target pair `(alpha,beta)` and right atom `nu`, let
`I_(alpha,beta,nu)` be the exact half-open interval list produced after
removing only the central guard and partitioning by the 39 endpoint atoms.
In the split embedding, an interval `I=[L,R)` contributes the literal BY
endpoint numerator

```text
D_B(I)
 =root^(w_B 169 L)-root^(w_B 169 R),

w_B=13+U_full[TARGET_B]=742599.                              (1)
```

The direct interval route and the inherited optimized endpoint route give

```text
BY_(alpha,beta)(nu)
 =sum_(I in I_(alpha,beta,nu)) D_B(I)                        (2)
```

for every one of the `13^2*39=6,591` atom entries.  Equation (2) was checked
by iterating all `7,108,460` split interval pieces; it was not inferred from
the old table hash.

Write

```text
G_(omega,nu)(tau)
 =safe(omega,tau)safe(nu,tau),
```

with the sheet shifts and chamber danger sets of the audited endpoint
factorization, and leave `AX_(alpha,beta)(omega)` preintegrated.  The old
endpoint matrix is

```text
E_t(omega,nu)
 =1/13^3 sum_(alpha,beta,tau)
    zeta_13^(beta-alpha-t tau)
    G_(omega,nu)(tau)
    AX_(alpha,beta)(omega) BY_(alpha,beta)(nu).              (3)
```

Replacing (2) inside (3) gives the one-leg object

```text
E^lit_t(omega,nu)
 =1/13^3 sum_(alpha,beta,tau)
    sum_(I in I_(alpha,beta,nu))
    zeta_13^(beta-alpha-t tau)
    G_(omega,nu)(tau)
    AX_(alpha,beta)(omega) D_B(I).                           (4)
```

All sums are finite, so distributivity gives exact Fubini equality

```text
E^lit_t(omega,nu)=E_t(omega,nu)                              (5)
```

without a convergence premise.  Computationally, (5) agrees in all
`13*39^2=19,773` entries and reproduces the pinned pair-bank digest

```text
c28119c8b54f47e5b7a46f1508fbba604b0e3997eaadb05b03ad28edd9aed468.
```

The target map then remains exactly

```text
A_t(ell,c)
 =1/DEN sum_(omega,nu) M(omega,nu,ell,c)E^lit_t(omega,nu).   (6)
```

All `13*7*13=1,183` values in (6) equal the prior coordinate bank.  Its
digest is

```text
989dafc220a6d09aeacfce4af0e9a4fe13eedacc79fa66032ea39bc107fd8efb,
```

and its complete spectrum digest is

```text
5f173227c5e203309f61bdfd9d47cc64a3b49ae8f14abd0f7bfc469eda278533.
```

Every refined residue `t` still has spectrum shape
`(91,1,6,12,72)`.  At the fixed class `(1,0,6)`, the mixed witness remains

```text
Ahat_6(1,1)=218019411785559321795 mod 572252886246508880869.  (7)
```

Thus the previous nonvanishing was not an artifact of storing BY as one
scalar per atom.

## Connection contract and exact type boundary

| field | exact content |
|---|---|
| source | the common-ancestry cell tensor `M(omega,nu,ell,c)` |
| target | the same refined endpoint output `A_t(ell,c)` |
| map | expand BY into tagged endpoint intervals, integrate the BY factor, and contract the resulting guarded endpoint matrix against `M` |
| preserved predicate | both guard atoms, every `tau` guard value, the source inverse-branch alignment, `ell`, `c`, and all thirteen refined residues |
| destroyed data | the location of the BY endpoint variable after (2), the AX endpoint/Q variable already hidden in AX, endpoint chronology, and every coupling between either endpoint variable and the ancestry point |
| needed sidecar | an actual measurable cospan or natural-extension section placing the endpoint variable on the inherited ancestry base, followed by the AX analogue and the grouped exact-address projector |

The variables remain typed as follows.

- `y`, the THM-2471 ancestry variable, is integrated inside `M` while
  `(ell,c)` are retained.
- `v`, the literal BY endpoint variable represented by `[L,R)`, is integrated
  on a separate endpoint factor.
- `x` and its `Q` overlap are still separately preintegrated inside AX.
- `alpha,beta,tau,t` are finite Fourier/address labels, not circle points or
  chronological maps.

Therefore the answer to the narrow temporal question is **no**: this is not
a lawful one-leg temporal transport onto the common ancestry stalk.  The
positive statement is an exact, lawful **iterated/product-space Fubini
lift**.  One may realize (6) on the external product of the source and BY
endpoint spaces, but no diagonal, section, or chronological cospan identifies
`v` with `y`.  Fubini changes integration order; it does not manufacture that
missing map.

## Positive and hostile controls

The positive controls reproduce the interval census, all 6,591 BY entries,
all 19,773 guarded endpoint entries, the frozen source tensor, all 1,183
retained coordinates, both prior output hashes, every residue support shape,
and (7).  The literal endpoint primitive is coded independently of
`fast_endpoint_sum`; both routes agree exactly.  Normal and `python -O`
replays are byte-identical, and all truth-bearing gates use explicit
`require` checks.

Four hostiles localize what equality and spectrum do not say.

1. Deleting only the right-leg guard changes `8,112/19,773` endpoint entries.
   The fixed relation nevertheless retains spectrum shape
   `(91,1,6,12,72)`.
2. Replacing the BY character endpoint differences by raw interval measures
   changes `8,788/19,773` entries.  It too retains all 91 fixed-relation
   Fourier modes.
3. Averaging away `ell` leaves only `(13,1,0,12,0)`.
4. Averaging away `c` leaves only `(7,1,6,0,0)`.

The first two controls sharpen the interpretation of the earlier spectral
closure: full support is a capacity property and survives severe typing
errors.  It cannot certify the guards, the character integrand, or temporal
lawfulness.  The latter two controls reconfirm that the retained common-base
coordinates jointly carry the mixed spectrum.

## Why the major frontiers remain open

- **THM-2449:** its table is a nonnegative rational integral of the complete
  lawful source/deep/word packet at its typed clock.  Equations (3)--(6) are
  split-field endpoint-character contractions on separate factors, not that
  response density.
- **THM-2512:** the shared `7 x 13` shape and all 72 mixed modes are only the
  algebraic target shape.  There is no rational lawful THM-2449 input table,
  no identified ANOVA defect `d_A`, and no Boolean owner/arrival/deep
  toothpick on one ancestry fibre.
- **`U_full`:** the endpoint leg uses the fixed `U_full` word, but no theorem
  identifies its endpoint tuple, clock, or exact relation address with the
  source-side `K=2` ancestry construction or with `U_clock`.
- **Current:** neither endpoint variable has been transported to `y`; no
  chronological arrival predicate, relative phase on one physical atom, or
  grouped coefficient `C(a;X,m)` is present.
- **LRC(14):** no live row is excluded, so the exact ledger is unchanged and
  LRC(14) remains open.

Expanding AX next would prove a useful three-factor Fubini identity and audit
the transformed `x/Q` sweep.  It would still leave three separate integration
variables.  The decisive next object is a typed endpoint-to-ancestry cospan
or diagonal support theorem, with the endpoint guard and exact-address
sidecars retained; only after that object exists does a two-leg literal
expansion become a temporal bridge test rather than another external-product
identity.

## Reproduction

```text
python -B 04-computation/lrc_r5_one_leg_integrand_fubini_by_probe_20260816.py
python -B -O 04-computation/lrc_r5_one_leg_integrand_fubini_by_probe_20260816.py
```

Both commands reproduce
`05-knowledge/results/lrc_r5_one_leg_integrand_fubini_by_probe_20260816.out`
byte-for-byte.  The semantic SHA-256 is

```text
6b87dc762517773254e3388dffd0bed9a5101685730dd3b93ebf4c639f224e49.
```
