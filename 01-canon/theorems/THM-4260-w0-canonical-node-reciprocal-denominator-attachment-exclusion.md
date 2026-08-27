---
id: THM-4260
title: "W=0 canonical-node reciprocal-denominator attachment exclusion"
status: >
  PROVED RELATIVE TO THM-4247/4249/4253 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. Intrinsically on the already-specialized W=0
  gate-interior fibre U*Z*(U+Z)!=0, none of the final 176 degree-34 or 104
  degree-42 source-target map classes collapses the twelve attachments.
  Thus all 864+648=1,512 live class-ratio incidences are impossible and
  S_34=S_42=empty on this fibre. No W-adic, neighbourhood, seam-entry,
  M=12, JC(2), or DC(2) conclusion follows.
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-4247-w0-involution-degree-twelve-attachment-exclusion
  - THM-4249-w0-cyclic-projector-missing-eigenline-attachment-squeeze
  - THM-4253-w0-degree-forty-two-norm-three-profile-exclusion
related:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4241-w0-hidden-cm-saturation-and-visible-hidden-index-four-gluing
  - THM-4251-w0-hidden-degree-twenty-four-attachment-exclusion
  - THM-4255-specialization-kernel-and-transverse-hasse-jet-repair
  - THM-4258-w0-three-sample-attachment-recurrence-and-two-torsion-sidecar
script: 04-computation/jc23_w0_canonical_node_attachment_exclusion_thm4260.py
output: 05-knowledge/results/jc23_w0_canonical_node_attachment_exclusion_thm4260.out
script_sha256: cff97177dc253d63b09ed2ab9bfd3a3da74512a86c01fa400350083b63ef780a
output_sha256: 829313e434bf60e4e4712d744e940cb577262ea08fefb5621e86b30790275d07
proof_only_output_sha256: 27e32aa329c4ed23ce859284109a9ecf8d6ee7bf207e975a5741639c4ba4e2e8
independent_script: 04-computation/jc23_w0_canonical_node_attachment_exclusion_independent_audit_thm4260.py
independent_output: 05-knowledge/results/jc23_w0_canonical_node_attachment_exclusion_independent_audit_thm4260.out
independent_script_sha256: 238af711b734a0e054a3da6d8a836f5f2996e34550d5e507111718d5eff2314f
independent_output_sha256: 481594bbddfba2451e3f7780916ffd6dd04ca7c394ae3b56311ceb5c3d1731f8
hash_basis: raw LF bytes
audit: >
  PASS/ACCEPT. The frozen compiler reconstructs the post-THM-4253 class and
  incidence universe, audits the node-orbit gauge, and performs every
  reciprocal-denominator calculation in exact GF(q)(t) arithmetic. The
  proof embedding q=397 attains the sharp characteristic-zero degree bound
  and has gcd t^2-1 for all 280 classes; q=577 independently reproduces all
  ten profile results. A hostile clean-room referee rebuilt the lattice from
  the standard-library THM-4249 audit and used a separate native polynomial
  engine. Ordinary, optimized, and fixed-hash-seed replays agree.
---

# THM-4260 -- `W=0` canonical-node reciprocal-denominator exclusion

**PROVED RELATIVE TO THM-4247/4249/4253 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THIS IS A THEOREM ON THE ALREADY-SPECIALIZED `W=0`
FIBRE ONLY.**

## 1. Statement and exact inherited frontier

Retain the exact-weight-twelve notation and normalization of the cited
theorems:

```text
C_0: x^6+y^4=1,                 E_0: Y^2=X^3+1,
O=Z[omega],                     omega^2+omega+1=0,
Q_j=tau^j Q_0,                  0<=j<12,
W=0,                            U*Z*(U+Z)!=0.          (1)
```

For `e in {34,42}`, let `S_e` be THM-4230/4249's finite set of marked
ratios `U/Z` for which some degree-`e` map sends all twelve attachments to
one target point. THM-4249 makes every remaining candidate have zero
`u`-eigenline coordinate. THM-4253 removes the degree-`42` norm-three
profile. The exact frontier is therefore

```text
degree 34: 176 source-target classes, 864 class-ratio incidences;
degree 42: 104 source-target classes, 648 class-ratio incidences. (2)
```

Every class has `24` raw representatives under target units and the source
`T` action. Its retained `(N(d),K)` profiles are

| degree | `(N(d),K): number of classes` | total |
|---:|:---|---:|
| `34` | `(4,10):36, (7,9):52, (13,7):32, (16,6):24, (19,5):24, (25,3):8` | `176` |
| `42` | `(9,11):24, (12,10):36, (21,7):32, (27,5):12` | `104` |

> **Theorem.** No class in `(2)` collapses the twelve attachments at a
> gate-interior node. Consequently
>
> ```text
> S_34=S_42=empty                                      (3)
> ```
>
> on the already-specialized `W=0` fibre `(1)`. Thus the inherited
> nonautomorphic planar Keller candidate is excluded on this fibre.

The native finite carrier is the bipartite graph

```text
source-target map class -- admissible CM-torsion unit orbit.       (4)
```

It has `280` class vertices, `82` untyped ratio vertices, and `1,512` edges.
The degree-`34/42` ratio envelopes have `55/34` vertices; seven occur in
both typed degrees. There is no intrinsic oriented pairwise observable, so
no tournament, orientation gauge, or tie convention is introduced.

## 2. The necessary hidden projection

Use THM-4241's integral basis `[u,f,g,h]`. Since the `u` coordinate is zero,
write a candidate as

```text
m=b f+c g+d h.                                        (5)
```

For the source involution `iota:(x,y)->(x,-y)`, define the THM-4247
anti-invariant projection

```text
ell=(1-iota)m=m-m composed_with iota
   =(2b+omega^2 d)f+(2c+d)g.                          (6)
```

If `m(Q_j)` is independent of `j`, then `iota Q_j=Q_(j+6)` gives

```text
ell(Q_j)=O                    for every j.             (7)
```

Only this necessary direction is used; the invariant component discarded
by `(6)` is never reconstructed. In the exact residual lattice,

```text
q(ell)=12K.                                           (8)
```

Put `A_ell=X_ell/x`, and let `d_ell(t)` be the **reduced coefficient
denominator of `A_ell=X_ell/x`**. The coefficient-form calculation of
THM-4230/4247 gives

```text
d_ell(t)=t D_ell(t^2).                                (9)
```

Thus `d_ell` is odd. The pole ledger counts a nonzero finite denominator
root with cost at least six, with the origin correction of two, and gives

```text
6 deg(d_ell)-2 <= 2q(ell)=24K.                        (10)
```

Integrality alone gives `deg(d_ell)<=4K`; oddness sharpens this to the
characteristic-zero bound

```text
deg(d_ell)<=4K-1.                                     (11)
```

Precomposition by `T` sends the attachment coordinate `t` to `-1/t`.
Consequently `(7)` requires a common gate-interior root of `d_ell` and its
reciprocal

```text
d_ell^*(t)=t^deg(d_ell) d_ell(-1/t).                  (12)
```

## 3. The `q=397` proof embedding and characteristic-zero lift

The frozen exact compiler reconstructs all `280` classes before performing
the denominator test. Grouping by `(degree,K)`, its complete proof universe
is

```text
(34,3):8, (34,5):24, (34,6):24, (34,7):32,
(34,9):52, (34,10):36,
(42,5):12, (42,7):32, (42,10):36, (42,11):24.         (13)
```

At the good embedding

```text
(q,zeta_12,r,s)=(397,157,161,27),                     (14)
```

exact rational-function arithmetic proves, class by class,

```text
deg(d_ell mod 397)=4K-1,
ord_t(d_ell mod 397)=1,
gcd(d_ell(t),d_ell^*(t)) mod 397=t^2-1.               (15)
```

The first equality is not merely a finite-field statistic: together with
the upper bound `(11)`, it proves that reduction has lost no denominator
degree. The second says that reduction has not raised the order of the
distinguished factor `t`.

For completeness, let `O_p` be the local DVR of the coefficient field at the
good prime above `397`, and put `n=4K-1`. The equality
`deg(d_ell mod 397)=n` makes the top coefficient of `d_ell` a local unit.
The equality `ord_t(d_ell mod 397)=1` makes the top coefficient of
`d_ell^*=t^n d_ell(-1/t)` a local unit as well. Thus the two polynomials can
be normalized **separately** to monic polynomials in `O_p[t]`. Any monic
common divisor over the fraction field is integral over `O_p` (Gauss lemma
plus integral closure of a DVR), and its reduction remains monic of the same
degree. Hence a characteristic-zero common factor of degree greater than
two would reduce to a common factor of degree greater than two in `(15)`, a
contradiction.

There is already a degree-two common factor. The points `t=+1,-1` are the
`iota`-fixed fibres, so `(6)` vanishes there and supplies `t^2-1`. Therefore

```text
gcd_K(t)(d_ell,d_ell^*)=t^2-1                         (16)
```

up to a nonzero scalar for every class. These two roots are exactly the
geometric wall:

```text
Z/U=((t^2-1)/(2t))^2=0.                               (17)
```

They lie outside the gate `U*Z*(U+Z)!=0`. Equations `(7),(12),(16),(17)`
exclude every class in `(13)` before any of its torsion-ratio neighbours is
consulted, proving `(3)`.

As a hostile field control, the independent embedding

```text
(q,zeta_12,r,s)=(577,57,224,25)                       (18)
```

reproduces the sharp degree and the monic gcd `t^2-1` in all ten profiles.
The proof does not average fields or infer a characteristic-zero zero from
a modular zero: `q=397`, `(11)`, and the monic-DVR argument are already
logically sufficient; `q=577` tests the specialization mechanism.

## 4. Canonical-node gauge and loss ledger

For fixed `U,Z`, the `24` radical choices split into two cyclic `C_12` node
orbits. The compiler chooses one. The source automorphism

```text
rho:(x,y)->(-x,-y)                                    (19)
```

exchanges the two. On the final `a_u=0` lattice, direct calculation for all
`280` representatives gives the classwise identity

```text
rho=(-omega^2)T^2.                                    (20)
```

The right side is a source action followed by a target unit. It therefore
fixes every source-target class while carrying the chosen node orbit to the
other. Auditing one node orbit loses no class and no collapse possibility.

Quotienting the raw maps loses a representative but retains degree, orbit
size `24`, the integral coefficient vector, hidden profile, and all incident
ratio tokens. Quotienting torsion points by target units retains their full
canonical orbit token; in particular the seven cross-degree overlaps are
stored once in the untyped union but remain typed on the edges. The frozen
ledger hashes are

```text
class ledger:  d95cc8b06d88e1b18159db224e79d156d9fe35d4067b9676f9300fa17933a7a2
edge ledger:   9642bfaf19497e7ed02f04093e1c970f8fa1f317c4ea8f0fecaf2895a841003e
ratio ledger:  69f35a4c5f80111f9d58df75e35f7ef57a6d249aec317ff44fe73f9c84d37b83
two-field denominator audit:
               63208069e8c76ba6693ea13e3b28806d4a3f4397115d9721b399347803b8bfa1 (21)
```

## 5. Orthogonal observer and specialization firewalls

THM-4258 is a related, not load-bearing, observer theorem. Its cubic
attachment recurrence reduces a direct equality table from `18,144` to
`4,536` group-value rows, but explicitly excludes no incidence. The present
theorem takes an orthogonal route: the hidden denominator obstruction kills
each map class before direct incidence evaluation. No concrete normalization
or two-torsion halving from THM-4258 is used.

THM-4251 is related rather than a dependency. It is a positive prototype for
a reciprocal-denominator obstruction, but none of its degree-`24` shell
calculations proves `(13)--(16)`.

THM-4255 supplies the mandatory specialization firewall. This theorem begins
**intrinsically on the already-specialized `W=0` fibre**. It does not infer a
coefficientwise statement before setting `W=0`, and the finite-field step is
arithmetic good reduction of degree-preserved monic polynomials, not
evaluation along a formal graph in the `W` direction. In particular, no

- `W`-adic strict transform,
- transverse or conormal Hasse jet,
- neighbourhood exclusion,
- entry into the whole exact `M=12` gate, or
- `JC(2)` or `DC(2)` conclusion

follows from `(3)`. Any transfer off `W=0` must retain the `W`-order and the
leading transverse jet as separate data.

## 6. Scope and reproduction

The theorem proves `S_34=S_42=empty` and closes the inherited candidate on
the gate-interior `W=0` fibre. It does not classify the hidden-Hom locus away
from that fibre, cross `U=0`, `Z=0`, or `U+Z=0`, close the full exact-weight-
twelve gate, prove seam entry, or prove `JC(2)`/`DC(2)`.

The compiler pins the exact upstream git revision and raw dependency hashes
before rebuilding the lattice and torsion carrier. From the repository root:

```bash
PYTHONHASHSEED=0 python3 -B \
  04-computation/jc23_w0_canonical_node_attachment_exclusion_thm4260.py \
  --repo . --proof-only

PYTHONHASHSEED=0 python3 -B \
  04-computation/jc23_w0_canonical_node_attachment_exclusion_thm4260.py \
  --repo .

PYTHONHASHSEED=0 python3 -B \
  04-computation/jc23_w0_canonical_node_attachment_exclusion_independent_audit_thm4260.py \
  --repo .
```

The proof-only stdout SHA-256 is
`27e32aa329c4ed23ce859284109a9ecf8d6ee7bf207e975a5741639c4ba4e2e8`.
The proof-plus-hostile stdout SHA-256 is
`829313e434bf60e4e4712d744e940cb577262ea08fefb5621e86b30790275d07`,
identical to the frozen output. Normal, optimized, and fixed-seed runs
byte-match. **QED.**
