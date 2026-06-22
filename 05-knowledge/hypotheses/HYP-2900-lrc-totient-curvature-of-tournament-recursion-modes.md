---
id: HYP-2900
status: COMPUTATIONAL SIGNAL / proof-route synthesis
source: codex-2026-06-22-S113
tags: [lrc14, euler-totient, coprime-density, multiplicative-functions, tournament-recursion, exact-period]
depends_on:
  - THM-442
  - THM-550
  - THM-554
  - HYP-2628
  - HYP-2629
  - HYP-2630
  - HYP-2886
related:
  - HYP-2899
  - HYP-2898
  - HYP-2897
  - HYP-2896
  - HYP-2895
  - HYP-2890
  - HYP-2883
  - HYP-2877
  - OPEN-Q-108
results:
  - 04-computation/lrc_tournament_totient_recursion_modes_codex_s113.py
  - 05-knowledge/results/lrc_tournament_totient_recursion_modes_codex_s113.out
---

# HYP-2900: totient curvature is the multiplicative defect exposed by tournament recursion modes

The three tournament recursion modes are cell-address recurrences:

```text
full:      A+B+C-D-E-F+G
even half: A+B-C
odd half:  A+B-C+D-E-F+G
```

They are exact for the full and half tiling cell counts, but they are not
recurrences for Euler-totient packet counts, coprime density, Hamiltonian paths,
or LRC witness counts.  HYP-2900's claim is that this failure is the useful
object: the signed residual is an **Euler-factor curvature** measuring which
CRT / coprime-density labels must be retained before scalarizing.

## Exact audit

Script:

- `04-computation/lrc_tournament_totient_recursion_modes_codex_s113.py`
- output: `05-knowledge/results/lrc_tournament_totient_recursion_modes_codex_s113.out`

The script verifies:

```text
n = sum_{d|n} phi(d)      for n<=120
```

and checks the three cell recurrences through `n<=80`.  All cell recurrences
are exact on their proper carriers.  In contrast, `phi` has nonzero residual in
every tested sample for all three recursion modes.

## Slot-size correction

The three formulas use different subtournament sizes.

Full mode at size `n`:

```text
A,B,C: n-1
D,E,F: n-2
G:     n-3
```

Even half mode at even size `n`:

```text
A,B: n-1
C:   n-2
```

Odd half mode at odd size `n`:

```text
A,B: n-1
C,D: n-2  (equal cardinality, different geometric slots)
E,F: n-3
G:   n-4
```

This matters at the LRC14 boundary:

```text
full n=14:
  A,B,C size 13; D,E,F size 12; G size 11

even_half n=14:
  A,B size 13; C size 12

odd_half n=15:
  A,B size 14; C,D size 13; E,F size 12; G size 11
```

So size `14` is not produced by the odd half recurrence; it is an input carrier
for the odd `15` recurrence.  The LRC14 composite `2*7` sits exactly where the
even half formula compares prime `13` against `12=2^2*3`.

## LRC14 boundary residuals

The audit records exact rational residuals for `rho(n)=phi(n)/n`:

```text
full n=14      residual rho = -2252/1001
even_half n=14 residual rho = -296/273
full n=15      residual rho = 766/455
odd_half n=15  residual rho = -218/385
```

The prime-exponent curvature at `n=14` is:

```text
full n=14:
  {2: 7, 3: 3, 7: 1, 11: -1, 13: -3}

even_half n=14:
  {2: 3, 3: 1, 7: 1, 13: -2}
```

Interpretation: the additive slot recurrence wants to reuse prime-size `13`
and the old `2^2*3` carrier.  The actual exact-period packet at size `14`
introduces the `7`-seam and a different dyadic exponent.  That is the same
reason HYP-2628/HYP-2630 insist on raw exact-period `phi` packets before the
squarefree/coimage/character projection.

## Relation to multiplicative functions

For `rho(n)=phi(n)/n`, the audit finds zero multiplicativity failures for
coprime pairs `a,b<=40`:

```text
rho(ab)=rho(a)rho(b) when gcd(a,b)=1.
```

The largest non-coprime defects are dyadic:

```text
rho(32*32)=1/2,  rho(32)rho(32)=1/4,  defect=1/4.
```

This separates two mechanisms:

1. `phi` / `rho` are multiplicative over CRT factors.
2. The tournament recurrences are additive boundary operators over slot sizes.

The bridge is not a scalar recurrence.  It is a labelled residual carrying the
prime exponents and CRT overlaps created by the slot boundary.

## Incoming character/totient synthesis

After this audit, incoming KPS S31q and mac-mini S44 made the same bridge from
the character side:

- KPS S31q identifies the sign words as arithmetic characters: `+++---+` is
  the Mobius / inclusion-exclusion channel, `++-+--+` is the Legendre `chi_7`
  channel on the seven subset slots, and the `++-` / omega package is the
  Eisenstein `chi_3` channel.
- mac-mini S44 identifies resonance killing as totient-weighted: a runner
  divisible by `b` kills all `phi(b)` primitive Farey points of denominator
  `b`, with `sum_{b<=14} phi(b)=64`, and the analytic floor is the familiar
  coprime-density / zeta(2) object.
- KPS S31r adds the parity-stratified operator composition: Mobius applies on
  all sizes, Eisenstein is the even/pronic order-2 mode `A+B-C`, and Legendre
  is the odd/square order-4 mode.  For LRC14 this reads
  `14=2*7 = Eisenstein(even fold 14 -> 7) o Legendre(odd apex 7)` on the
  Mobius floor.

S113 adds the slot-size correction to that picture.  The character signs say
which channel is being sampled; the slot atlas says which integer sizes and
prime exponents enter that sample.  At `n=14`, the even-half channel compares
two size-13 prime slots with a size-12 dyadic/triadic slot while the actual
exact-period carrier is `14=2*7`.  Hence the residual is not arbitrary
noise: it is the character-labelled Euler curvature caused by asking an
additive slot operator to act on multiplicative exact-period packets.
S31r explains why this exact seam is the right one to inspect: the even
Eisenstein fold exposes the apex prime `7`, and the residual records the
Euler-factor cost of passing from the even fold into the odd Legendre channel.

The strengthened proof target is therefore:

```text
exact-period phi packets
-> Mobius / chi_7 / chi_3 character decomposition
-> parity composition 2q -> q (Eisenstein fold, then Legendre apex)
-> Euler-factor slot curvature with correct subtournament sizes
-> CRT/chi7/coimage signed-current packets
-> HYP-2890 Fejer residual leak and scalar cap only at the end.
```

## Proof-route implication for LRC14

HYP-2900 refines HYP-2899's product-Mobius packet ledger and reinforces the
current LRC14 architecture:

```text
exact-period phi packets
-> Euler-factor curvature residual
-> slot-labelled recursion atlas
-> half-tiling / mirror fixed-line address
-> scalar cap, floor, or Fejer estimate only after labels are retained.
```

This dovetails with:

- HYP-2628/HYP-2629: exact-period `phi` packets refine squarefree mask mass.
- HYP-2630/HYP-2632/HYP-2883: copy mass alone is not enough; retain `chi_7`,
  affine lane, and signed-current labels.
- HYP-2898: scalar additive energy is too coarse; labelled Fejer/residual data
  is the durable route.
- HYP-2897/HYP-2896/HYP-2895: exact-tiler boundary and one-tail disproof
  branches route away, leaving the moderate resonant support-six middle.

The next theorem should not be "totient obeys the tournament recurrence."
It should be:

```text
the Euler-curvature residual of the three slot recurrences factors through
exact-period CRT packets, and after deleting divisibility-killed denominators
the remaining labelled curvature is controlled by the HYP-2883/HYP-2890
signed-current / Fejer residual machinery.
```

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
exact_period_phi_packets
> euler_curvature_residual
> three_recursion_slot_atlas
> half_tiling_fixed_line
> scalar_cell_count
> raw_subtournament_size
> raw_runner_vertices
```

Pairwise observable: number of proof-relevant labels preserved, overlap, and
declaration order.  Fingerprints:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
Hamiltonian_path_count=1
```

Assumption challenged: the `A..G` formulas should directly transfer to
coprime-density or LRC witness counts.  They do not.  They are address
operators; their multiplicative residual is the useful packet data.
