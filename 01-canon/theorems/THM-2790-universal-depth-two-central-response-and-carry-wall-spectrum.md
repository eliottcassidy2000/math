---
id: THM-2790
title: "Universal depth-two central response and carry-wall spectrum"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  nonzero depth-two central
  direction of the canonical THM-2625 endpoint current, all 169 forward
  differences and all 156 nontrivial cycle characters survive in both
  certified finite fields.  This gives 28,392 spatial responses and 26,208
  modes, including the full 2,184/2,016 carry-wall censuses.  It is an
  unnormalized coefficient-side quotient theorem, not physical high-sheet
  descent, a positive magnitude floor, a line bundle, a row exclusion, or
  LRC(14).
source: root/twin-transgression-central-response-2026-07-28
depends_on:
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
related:
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
  - THM-2791-full-arm-orbit-transfer-and-lower-central-chord
script: 04-computation/lrc14_universal_central_response_thm2790.py
output: 05-knowledge/results/lrc14_universal_central_response_thm2790.out
independent_script: 04-computation/lrc14_universal_central_response_independent_audit_thm2790.py
independent_output: 05-knowledge/results/lrc14_universal_central_response_independent_audit_thm2790.out
script_sha256: f2501a8b9b6a3b8d8d3f30b0bc5267f9b3f161f7fe71f6333d0b5fb36e261b82
output_sha256: 720b118e98c25d95195ac824647f534a14cb964fc967624ecbf68f4b178c4f9e
independent_script_sha256: c61470e08927f67896b6a07c6a3a346a68e48e57ee0a162aa06f3d2196643b17
independent_output_sha256: d6a2e6e532ac14e9a3c9bb6ea343d9879df0822ca606cf766b27b81d62be1787
hash_basis: LF-normalized bytes
---

# THM-2790 -- universal depth-two central response and carry-wall spectrum

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2625 constructs one canonical coefficient-side endpoint current with full
joint support, and THM-2779 supplies the normalized Heisenberg endpoint-deck
action.  This theorem computes the entire first central response of that
current.  THM-2788 identifies this action as the common depth-two commutator
shadow of the carry-suppressed and physical modular extensions.  The result is
therefore quotient-central, not a physical high-sheet descent theorem.

## 1. Current and central response

Let `P_L,Q_R` be the separate endpoint factors of the canonical THM-2625
marked current.  For every nonzero

```text
s in V=F_13^2
```

define the **unnormalized** fixed-increment endpoint current

```text
J_s(R)=P_(R+s) Q_R.                                             (1)
```

The fully normalized THM-2625 marked current is
`C0 J_s(R)/169^2`, where the inherited scalar `C0` is nonzero.  Thus
all support and Fourier-nonvanishing statements below are unchanged by this
normalization, while the finite-field companion lawfully works with the
cyclotomic-integer numerator `(1)`.

The central endpoint-deck generator associated to \(s\) is

```text
Z_s:R -> R+s.                                                   (2)
```

The intrinsic first response is therefore

```text
((Z_s-1)J_s)(R)
 =P_(R+2s)Q_(R+s)-P_(R+s)Q_R.                                  (3)
```

The probe reconstructs `P,Q` in both certified THM-2625 finite-field
specializations and evaluates (1)--(3) for every

```text
168 nonzero s * 169 origins R = 28,392 cells per field.         (4)
```

## 2. Maximal spatial response

In both fields:

```text
support(J)=28,392/28,392,
support((Z-1)J)=28,392/28,392.                                 (5)
```

For every individual direction \(s\):

```text
support_R((Z_s-1)J_s)=169/169.                                 (6)
```

Translation by \(s\) partitions the endpoint plane into thirteen cycles,
equivalently the determinant fibres

```text
Delta=det(s,R) in F_13.
```

Indeed, after choosing any normalized transverse vector `t` with
`det(s,t)=1`, every origin is uniquely `R=w s+v t` and
`det(s,R)=v`.  Thus `w` is the central-cycle coordinate and `v` is the
two-digit low coordinate from THM-2779.

All thirteen cycles are active for every \(s\).  In particular, the
two-digit carry wall

```text
Delta=12
```

has

```text
2,184/2,184 nonzero responses,
13/13 for every s.                                              (7)
```

This is a genuinely new product-current noncancellation statement beyond
THM-2625 full vertex support and THM-2779's separate signed-edge/Hadamard
gate.  It is not implied by either earlier support statement.  The canonical
current is nowhere constant along the central direction.

A sharp abstract hostile isolates that distinction.  On one `C_13` cycle,
put

```text
P_j=zeta_13^j,                 Q_j=zeta_13^(-j).
```

For odd order, every adjacent sum and difference of each separate factor is
nonzero, so the THM-2779 edge gates pass.  Nevertheless
`P_(j+1)Q_j=zeta_13` is constant: its central response and all twelve
nontrivial central modes vanish.  The exact computation below is therefore
load-bearing rather than a formal corollary of the separate-factor gates.

## 3. Every nontrivial central character survives

Fix one determinant cycle and enumerate it as

```text
R_j=R_0+j s,                    j in F_13.
```

For \(k=1,\ldots,12\), put

```text
Jhat_s,Delta(k)=sum_j J_s(R_j) zeta_13^(-k j).                  (8)
```

The exact response identity is

```text
widehat((Z_s-1)J_s)(k)
 =(zeta_13^k-1) Jhat_s,Delta(k).                               (9)
```

The probe checks (9) on every cycle in both fields.  More strongly,

```text
Jhat_s,Delta(k) !=0
```

for every

```text
s!=0,          Delta in F_13,          k in F_13^*.             (10)
```

Thus each field has

```text
168*13*12=26,208/26,208 nonzero central modes,
156/156 for every s.                                           (11)
```

On the carry-wall cycle alone:

```text
168*12=2,016/2,016 nonzero modes,
12/12 for every s.                                             (12)
```

The two finite-field maps are certified cyclotomic ring homomorphisms from
the THM-2625 construction.  A nonzero image proves that the corresponding
characteristic-zero cyclotomic integer is nonzero.  Thus (5)--(12) are
rigorous coefficient-side nonvanishing statements, not numerical
approximations.

## 4. Deterministic digests

Digest serialization is:

```text
J and response:
  s in lexicographic order, then R in lexicographic order;

central modes:
  s in lexicographic order;
  determinant cycles ordered by their lexicographically least point;
  character k=1,...,12;

each field element:
  one unsigned eight-byte big-endian integer.
```

Changing the chosen origin on a cycle multiplies a nontrivial mode by a
thirteenth root of unity.  Hence the displayed mode digests and phases use the
stated lexicographic origin convention; only their nonvanishing is
origin-independent.

For

```text
p=352341050142921841:
  J        fb738f7aca8dc6c08e0a6178e0fd9603f60deffec4bb26e16c7b9f4abae33d3f
  response 4bf698c20fe9099ade588112c44f8d0506a4a23b7261a2db3d2c744c240476da
  modes    29916cc2fb6a66c939cda174be1623e57643e7d4a0707702e2b4f3a66189b4c2

p=956354278959359281:
  J        af3c9c215bf56b9ecb80f0f6d83ac3d746aed7cf0cf65100d94ce3af0fe5fe71
  response 7635f9584e36aa3227432125749393f35169ad7bfba1bc3220087625db44cb72
  modes    e07cd3af43dd3969c2d4afb5706debc178375e04e76440f0a624d1aebe840ff4
```

## 5. Gauge audit

Three gauge statements must be separated.

1. **THM-2625 representative gauge.**  The separate factors `P` and `Q`
   already descend through the representative gauge.  Therefore (1)--(12)
   are well-defined on its endpoint coefficient plane.
2. **THM-2779 frame torsor.**  Changing the normalized frame changes the
   central direction \(s\).  Since (5)--(12) hold for every nonzero \(s\),
   the support predicate and exact support lower bounds are
   frame-choice invariant.  The individual coefficient values, phases,
   and digests are not thereby identified across frame gauges; no
   `H_13`-equivariance of the current has been proved.
3. **THM-2779 decoder torsor.**  The THM-2625 current does not use a target
   decoder choice.  Its central response therefore does not select a point
   of the decoder torsor or build the missing decoder/frame naturality
   square.

There is also an origin gauge inside a fixed determinant cycle.  If
`det(s,t)=1` and `t` is replaced by `t+c s`, then the chosen origin on the
`Delta` cycle advances by `c Delta` central steps.  Formula (8) is therefore
multiplied by `zeta_13^(k c Delta)`.  Its value and digest change, while its
nonvanishing does not.

The exact lower bound supplied here is a support lower bound:

```text
all 28,392 spatial responses and all 26,208 nontrivial modes survive.
```

Finite-field nonvanishing does not provide a positive Archimedean magnitude
floor for a chosen complex embedding.

## 6. Physical and holotopy boundary

On the depth-two carrier `Omega=Z/169`, both the carry-suppressed
Heisenberg action and the physical modular action have the same central
quotient generator

```text
Z_1=O^13 mod <Z_2>.
```

Equations (5)--(12) prove complete charge under that quotient center.
They do not choose between the two extensions, because those extensions
have the same alternating commutator and differ in their \(p\)-power
Bockstein coordinate.

On the full physical depth-six carrier, `Z_1` has order \(13^5\) and

```text
[X,Z_1]=Z_2.
```

THM-2782 supplies a semantic physical `Z_1` segment with nonzero local
thirteen-bin polyphase modes, but its coefficient changes from nonzero to
zero after the `j=13` shift `Z_2`.  It therefore does not descend to the
depth-two central cycle.

Consequently, the THM-2625 result here and THM-2782 provide **parallel
nonzero spectra**, not a line-bundle or charge-matching map.  A fibrewise
Fourier ratio can be written after arbitrary choices because both sides
are nonzero, but it does not commute with:

```text
decoder/frame gauge,
the physical p-power map,
the Z_2,...,Z_5 lower-central transitions,
semantic ancestry.
```

The nonzero-to-zero THM-2782 transition is not invertible descent data.
The data are therefore better modeled by a filtered or constructible
coefficient system with a vanishing boundary than by a line bundle with
invertible transitions.  No sheaf or gluing functor is constructed here.

## 7. Exact companions and audit status

Run:

```bash
python 04-computation/lrc14_universal_central_response_thm2790.py
python -O 04-computation/lrc14_universal_central_response_thm2790.py
python 04-computation/lrc14_universal_central_response_independent_audit_thm2790.py
python -O 04-computation/lrc14_universal_central_response_independent_audit_thm2790.py
```

The primary companion reconstructs both certified THM-2625 finite-field
endpoint factors, gates the exact primes and deterministic digests, checks
all spatial and Fourier support counts, and includes the constant-product
edge-gate hostile.  Normal and optimized runs byte-match

```text
05-knowledge/results/lrc14_universal_central_response_thm2790.out
```

after LF normalization.

The independent companion does not use the primary orbit partition or
analysis routine.  For each nonzero `s`, it chooses a canonical transverse
`t_s` with `det(s,t_s)=1`, parametrizes every determinant fibre directly
as

```text
R(w,Delta)=w s+Delta t_s,
```

and independently rechecks the plane partition, forward difference,
Fourier-multiplier sign, all support counts, and both carry-wall counts in
both fields.  Its normal and optimized runs byte-match

```text
05-knowledge/results/lrc14_universal_central_response_independent_audit_thm2790.out.
```

The independent proof audit also rederived the representative/frame gauges,
the characteristic-zero specialization implication, the THM-2782/2788
quotient boundary, and the unnormalized-current typing.  The current script
hashes are recorded in the frontmatter.

No Boolean-carrier identification, same-ancestry wing-current attachment,
physical high-sheet descent, positive Archimedean amplitude floor, row
exclusion, or LRC(14) decrement follows.

QED.
