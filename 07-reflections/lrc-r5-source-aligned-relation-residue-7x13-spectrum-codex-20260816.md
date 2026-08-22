# The source-aligned fixed relation has full 7 x 13 mixed spectrum

**Status: FINITE-EXACT CANDIDATE; independent audit pending.**  This note
records a positive source-time simple-kernel spectrum and its exact typing
boundary.  It is not a physical relation current, a THM-2512 bridge theorem,
an exact-address coefficient, a scalar-row exclusion, or an LRC(14) result.

## Inheritance pass

The closest proved mechanism is THM-2594: its seven half-open cells live on
one common ancestry base and support a nonzero theta-slaved contraction.  The
closest target mechanism is THM-3479/THM-3514: the fixed `U_full` endpoint
bank has a complete refined `(1,0,t)` residue transform and full endpoint
`K4 x F_13` spectrum.  The canonical hostile is the raw source-to-arrival
guard transition: at `K=2` its `39 x 39` matrix has full support but rank
three, with all thirteen source-sheet rows identical inside each chamber, so
every nontrivial source-sheet character is erased.  The least-used relevant
sidecar is THM-2471's retained inverse branch `n in {0,...,168}`; its high
base-13 digit recovers the source sheet exactly.

The preceding source-aligned probe used that sidecar to realize

```text
e_omega=e P_omega,
f_omega^src=1_Q P^2(e P_omega),
```

before the two packets meet on one first-collision base.  It retained the
common torsor offset `c in F_13` but integrated out THM-2594's cell label.
The present computation performs the missing refinement before that final
marginal.

## The retained tensor and its gates

For endpoint/source-time guard atoms `omega,nu`, seven-cell label
`ell in F_7`, and common offset `c in F_13`, define the integer numerator

```text
M(omega,nu,ell,c)
```

by inserting `cell_ell(y)` into the exact product profile on the common base.
The seven cells are the half-open partition

```text
cell_ell={y: ||y-ell/7||<1/14}.
```

Every cell has exact grid mass `42548128262640`.  The tensor has
`362/1521` supported atom pairs and `5150/138411` supported
`(omega,nu,ell,c)` entries.  Summing over `ell` recovers every one of the
previous `39 x 39 x 13` common-offset entries, not only its total or support:

```text
sha256(M)=39d7a0b4e5b2d8b85631d682ed1967091e44dc41e17b33a77e7184d3dc93e0cf.
```

This is the important order of operations.  The cell is inserted before
integration; it is not reconstructed from an already marginalized owner
profile.

## Complete residue-bank result

Let `E_t(omega,nu)` be the actual preintegrated endpoint atom-pair matrix for
the refined THM-3479 class `(1,0,t)`.  Contract on the source tensor:

```text
A_t(ell,c)
 =sum_(omega,nu) M(omega,nu,ell,c) E_t(omega,nu) / DEN.
```

All arithmetic is reduced through the certified split embedding

```text
p=572252886246508880869,
zeta_7=353818603057943120846,
zeta_13=505438565698892403012.
```

The denominator is a unit modulo `p`.  For every `t in F_13`, the complete
two-dimensional transform

```text
Ahat_t(h,k)
 =sum_(ell,c) A_t(ell,c) zeta_7^(-h ell) zeta_13^(-k c)
```

has support

```text
(total, DC, F_7 axis, F_13 axis, mixed)=(91,1,6,12,72).
```

Thus all thirteen refined residues have all `91/91` coefficients nonzero,
including all `13*72=936` mixed coefficients.  The `h=0` rows recover all
thirteen previously computed owner profiles exactly.  This makes the result
strictly stronger than an accidental nonzero scalar at the fixed relation
class.

The fixed all-unit relation is `(1,0,6)`.  Its endpoint aggregate is

```text
225010624370142818572 mod p,
```

and its first mixed witness is

```text
Ahat_6(1,1)=218019411785559321795 mod p.
```

After doubly centering the `7 x 13` table itself, both axes vanish and every
one of the `6*12=72` mixed modes remains nonzero.  Independently doubly
centering the `39 x 39` endpoint pair matrix before source contraction still
leaves all `91/91` output modes nonzero.  Hence neither endpoint main effects
nor output-table main effects account for the mixed signal.

## The hostiles localize both load-bearing coordinates

Two exact coordinate-erasure controls separate axis support from coupling.
Replacing every cell row by the cell average leaves exactly

```text
(13,1,0,12,0),
```

while replacing every owner column by the owner average leaves

```text
(7,1,6,0,0).
```

In both controls all 72 mixed modes vanish.  The positive mixed block is
therefore carried jointly by the retained base cell and the retained common
owner offset.  This does not prove a causal uniqueness statement about any
deeper THM-2594 slaving coordinate; it only identifies the two coordinates
load-bearing in this tensor.

## Connection contract and exact loss boundary

| field | exact content |
|---|---|
| source | source-time Boolean tensor `M(omega,nu,ell,c)` on one THM-2471 base |
| target | refined endpoint residue matrix `E_t(omega,nu)` |
| map | finite contraction over the two aligned guard atoms |
| preserved | source sheet via inverse branch, seven-cell label, common offset, all thirteen refined residues |
| destroyed | endpoint chronology, exact relation address, the separate endpoint integration variables, characteristic-zero order/sign |
| needed sidecar | an integrand-level realization of `AX` and `BY` on the same temporal stalk, plus the grouped exact-address projector |
| cheapest next test | replace one preintegrated endpoint atom table by its interval integrand while retaining `(ell,c)` and verify Fubini/type equality before attempting both legs |

The shape now matches the `F_7 x F_13` mixed-frequency object used by
THM-2512, and the fixed relation has the strongest possible nonvanishing
inside that finite shape.  But `E_t` is still formed from separately
preintegrated `AX` and `BY` values.  The computation therefore proves a
residue pushforward simple kernel, not a THM-2449 response table or a lawful
physical current.  It also remains an unrestricted class sum rather than
the grouped exact-address coefficient `C(a;X,m)`.

This shifts the frontier cleanly: spectral capacity is no longer the first
missing gate for the source-aligned fixed relation.  The missing gate is a
temporal/Fubini theorem that transports the endpoint integrands themselves,
followed by the `U_full`/`U_clock` same-tuple reconciliation.  LRC(14)
remains open.

## Reproduction

```text
python -B 04-computation/lrc_r5_source_aligned_relation_residue_7x13_spectrum_probe_20260816.py
python -B -O 04-computation/lrc_r5_source_aligned_relation_residue_7x13_spectrum_probe_20260816.py
```

The coordinate-bank, spectrum-bank, and semantic SHA-256 values are

```text
989dafc220a6d09aeacfce4af0e9a4fe13eedacc79fa66032ea39bc107fd8efb
5f173227c5e203309f61bdfd9d47cc64a3b49ae8f14abd0f7bfc469eda278533
cd55336bb1dfe5f37f020c242c4bca5b7c6be339ec57e95d69e10bbe68d9dbaa.
```
