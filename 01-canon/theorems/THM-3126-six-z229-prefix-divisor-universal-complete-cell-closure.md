---
id: THM-3126
title: "Six z229 prefix divisor-universal complete-cell closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Six explicit projected z1=229 three-drift
  prefixes cannot be extended by any fourth tail: every nonaligned
  denominator projection exceeds the sharp translated danger-band capacity,
  while an aligned extension lies in THM-2928's closed four-aligned branch.
  This statement is independent of candidate THM-3113's screen exhaustion.
  It does not promote THM-3113, close z1=228, lower the proved k=3 cap 229,
  classify physical covers outside these prefixes, or prove LRC(14).
source: root/frontier-synthesis-2026-08-02
audit: >
  An independent hostile audit rebuilt every closed safe cell from literal
  interval intersection with the strict danger arcs, reduced the resulting
  Boolean group-ring carriers modulo every divisor, and reproduced all six
  sizes, all 792 support inequalities, the 744 density closures, the 48
  full-support boundary pairs, the ten exceptional moduli, and the sharp
  unit slack.  It separately checked the complete-cell projection typing,
  the d=1 aligned handoff, the d=2 cardinality hostile, and the translated
  d=28 boundary.  Fresh normal and optimized runs are byte-identical after
  LF normalization to the stored transcript and the declared hashes match.
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
related:
  - THM-3111-projected-k3-z230-exact-screen-and-compressed-complete-cell-descent
  - THM-3113-projected-k3-z229-terminal-and-z228-screen-double-layer-descent
  - MISTAKE-334
script: 04-computation/lrc14_j7_k3_z229_six_prefix_divisor_universal_support.py
output: 05-knowledge/results/lrc14_j7_k3_z229_six_prefix_divisor_universal_support.out
script_sha256: ba8ef1ef8d93a5d51f5933ab9753c62b6012f24661e108455b4976fc827ce7fe
output_sha256: 1a820b6b356662587c563b84b3a56a9bbd24b70fbbe98b87d17bd8ecbc0bacf1
semantic_sha256: 9f58b66c5ed2930023e1023efeb9fb9d873a856afe895ff2768a31889bd9500a
hash_basis: LF-normalized bytes
---

# THM-3126 -- six z229 prefix divisor-universal complete-cell closure

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**

## 1. Exact scoped statement

For a six-body set `E`, its body ruler `L`, and two fixed low drift labels
`B`, let `J(E,B)` be the complete closed `1/L` cells safe from every label in

```text
E union {229} union B.                                  (1)
```

The danger combs are strict-open, so the weak endpoint inequalities used to
construct `J(E,B)` make it an inner carrier.  The six carriers are

| `E` | `L` | `B` | `|J(E,B)|` |
|---|---:|---|---:|
| `(1,4,5,9,11,14)` | `194040` | `(234,243)` | `34442` |
| `(1,5,6,8,9,14)` | `35280` | `(234,246)` | `6812` |
| `(1,5,9,11,12,14)` | `194040` | `(234,243)` | `35738` |
| `(1,8,10,11,12,14)` | `129360` | `(260,312)` | `25358` |
| `(2,5,8,9,11,14)` | `388080` | `(234,246)` | `73130` |
| `(2,8,10,11,12,14)` | `129360` | `(260,364)` | `25798` |

For every divisor `d>1` of the corresponding ruler, reduction of cell
addresses satisfies

```text
|pi_d J(E,B)| > ceil(d/7).                              (2)
```

Consequently no prefix `(1)` can be completed by any fourth tail and the
remaining aligned tails.  This closes the six explicit prefix families, not
the full projected `z1=229` layer.

## 2. Divisor-support proof

Every residue modulo `d` has exactly `L/d` lifts.  With
`rho=|J(E,B)|/L`, pigeonhole therefore gives

```text
|pi_d J(E,B)| >= ceil(rho d).                           (3)
```

The minimum density among the six carriers is

```text
rho_0=17221/97020,
rho_0-1/7=3361/97020.                                   (4)
```

For `d>=25`,

```text
(rho_0-1/7)d-6/7
 >= (rho_0-1/7)25-6/7
 =173/19404>0.                                          (5)
```

Since `ceil(d/7)-d/7<=6/7`, equations `(3)`--`(5)` prove `(2)` on the
uniform tail.  Direct carrier-specific density closes `744` of the complete
`792` carrier/divisor pairs.

Only `48` small-modulus pairs tie at the coarse bound.  Their moduli are

```text
2,3,4,5,8,9,10,11,15,22.                               (6)
```

Exact reduction in two independent aggregations shows that every one of
these `48` projections is the full `Z/dZ`.  Thus `(2)` holds in every case;
the minimum support-minus-capacity slack is one.

## 3. Completion consequence

Let `z` be a proposed fourth tail.  If its endpoint-grid denominator is
`d>1`, then at each normalized local coordinate its bad cell addresses form
one translated strict-open cyclic band.  Proved THM-2984 gives the sharp
capacity

```text
kappa(d)=ceil(d/7).                                     (7)
```

Equation `(2)` leaves one complete prefix-safe cell outside that band at
every local coordinate.  Hence the lossless projected residual of THM-2941
is the full circle:

```text
P_(E,{229} union B union {z})=T.                        (8)
```

Literal completion would imply `P subset U_A` by THM-2941 `(25g)`, whereas
THM-1166 gives `mu(U_A)<=36/91<1` for the three aligned multipliers.  This
contradicts `(8)`.  If `d=1`, then `z` is aligned and the row lies in the
four-aligned/three-drift branch already closed uniformly by THM-2928.

This proof is independent of the sign, primitive unit, and height of `z`.

## 4. Hostile boundary and correction discipline

At `d=2`, the coarse lower bound equals the translated capacity:

```text
ceil(2|J|/L)=1=kappa(2).                                (9)
```

A same-cardinality set supported only on even cells realizes support one,
while every actual carrier above has support two.  Therefore total carrier
mass cannot replace the residue-support sidecar.  The verifier constructs
this hostile explicitly.

MISTAKE-334 is also load-bearing.  At `d=28`, the centered count is
`beta(28)=3`, but an arbitrarily translated strict band has capacity
`kappa(28)=4`.  Every gate here uses `(7)`, never the centered proxy.

## 5. Unconditional result versus candidate consumer

The definitions, all `792` divisor checks, and the six prefix closures are
independent of THM-3113.  Candidate THM-3113 is needed only for the stronger
composition claims that:

1. its `z1=229` screen is exhaustive;
2. all screen residuals land on precisely these six prefixes; and
3. the complete `z1=228` layer has no residual state.

Until that screen receives its required immutable independent replay,
THM-3126 does **not** prove the proposed ledger subtraction
`374263-48=374215` or the proposed cap `z1<=227`.  The proved projected
`k=3` cap remains `229` by THM-3111.

## 6. Exact verification

The companion reconstructs every complete cell directly from integer weak
endpoint inequalities.  For every divisor it independently computes the
support as a residue set and as a byte histogram, requiring equality before
using the cardinality.  It freezes all six carrier hashes, the semantic
digest, the `744+48` partition, the same-size parity hostile, and the
MISTAKE-334 regression.  There are no truth-bearing Python `assert`s.

Normal and optimized commands are

```bash
python3 04-computation/lrc14_j7_k3_z229_six_prefix_divisor_universal_support.py --output /tmp/thm3126.out
python3 -O 04-computation/lrc14_j7_k3_z229_six_prefix_divisor_universal_support.py --output /tmp/thm3126.opt.out
cmp /tmp/thm3126.out /tmp/thm3126.opt.out
```

Both modes are byte-identical to the stored transcript.  The source, output,
and semantic LF-normalized SHA-256 values are those in the frontmatter.

This is a proved finite-exact closure of the six displayed prefix families.
The independent audit does not enlarge the scope stated in Section 5.
