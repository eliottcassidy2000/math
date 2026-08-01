---
id: THM-3035
title: "Level-two Farey anharmonic quartic S3 orbit diamond"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.
  The anharmonic S3 action on P1(F_p) has an exact boundary, harmonic,
  equianharmonic, and generic orbit atlas, with structural p=2 and p=3
  degenerations.  Its co-occurrence cover retains the stabilizer phase lost
  by the coarse orbit.  The PSL/PGL reflection sign and an exact THM-2056
  Farey hostile show why the finite frame does not determine an integral gate.
source: codex-modular-farey-anharmonic-orbit-diamond-2026-08-01
depends_on:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-2996-prime-modular-affine-defect-trichotomy-and-spherical-quartic-uniqueness
related:
  - THM-2998-quartic-star-complement-cross-wall-real-locus-and-c3-x0-14-quotient
  - THM-3034-ordered-quartic-cross-wall-x1-14-and-diamond-quotient
script: 04-computation/modular_farey_anharmonic_quartic_orbit_diamond_thm3035.py
output: 05-knowledge/results/modular_farey_anharmonic_quartic_orbit_diamond_thm3035.out
script_sha256: af88b3644cbea786e39de2509ab01c2ff55407f3769007690f009b200f0e3bd2
output_sha256: f5cf56eb89f8e25c2be7fdc823edb2d77c9c747f162f31dd129dd474e853c352
hash_basis: LF-normalized bytes
---

# THM-3035 -- the Farey / anharmonic / quartic `S_3` orbit diamond

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

## 1. Inheritance

[THM-2056](THM-2056-kelvin-polar-farey-defect-certificate.md) proves that a
Farey child defect contains an integral Gram cross term; mod-two directions
alone do not determine its sign.  [THM-2632](THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar.md)
identifies a mod-two Farey frame with an ordering of the three nonzero
directions of `V_4=F_2^2`.  The six orders form a regular `S_3`-set.
[THM-2984](THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate.md)
isolates the six anharmonic transforms, while
[THM-2996](THM-2996-prime-modular-affine-defect-trichotomy-and-spherical-quartic-uniqueness.md)
proves the exceptional quartic identity

```text
AGL_2(F_2)=S_4,                         S_4/V_4=S_3.       (1)
```

Finally, [THM-3034](THM-3034-ordered-quartic-cross-wall-x1-14-and-diamond-quotient.md)
identifies two ordered quartic sign sheets with a common cyclic order-three
phase and an odd sheet exchange.

The missing object is an exact orbit-and-stabilizer dictionary that retains
which binary and ternary faces co-occur, together with the orientation and
integral coordinates destroyed by its coarse projections.

## 2. The prime anharmonic orbit theorem

For a prime `p`, let

```text
H_p=<s,c> <= PGL_2(F_p),
s(r)=1/r,                         c(r)=-1/(1+r).          (2)
```

The three points `{infinity,0,-1}` are permuted faithfully, so `H_p=S_3`.
For every prime `p>=5`, the action on `P^1(F_p)` has the following named
orbits:

```text
boundary:       {infinity,0,-1}                 =S3/C2,
harmonic:       {1,-2,-1/2}                     =S3/C2,
equianharmonic: {r:r^2+r+1=0}                   =S3/C3,
generic:        every remaining orbit           =S3/1.    (3)
```

The equianharmonic orbit is present exactly when `p=1 (mod 3)`.  Therefore
the number of generic regular orbits is

```text
(p-7)/6  if p=1 (mod 3),
(p-5)/6  if p=2 (mod 3).                              (4)
```

The two small primes are not harmless missing cases:

```text
p=2: P^1(F_2) is the boundary orbit; no transverse point remains.
p=3: boundary={infinity,0,-1}, plus {1}=S3/S3.          (5)
```

At `p=3`, the harmonic and equianharmonic conditions collapse at the single
`S_3`-fixed point `1`.  Thus prime two is an exhaustion and prime three is a
fixed-point collapse.

For the LRC-relevant prime `13`, the atlas is

```text
boundary        {0,12,infinity},
harmonic        {1,6,11},
generic         {2,4,5,7,8,10},
equianharmonic  {3,9}.                                  (6)
```

Their permutation characters on `(identity, transposition, 3-cycle)` are

```text
regular six:   (6,0,0),
natural three: (3,1,0),
parity two:    (2,0,2),
trivial one:   (1,1,1).                                 (7)
```

This gives the exact representation ladder

```text
six Farey/V4 orders                    regular S3,
three nonzero V4 directions/matchings  S3/C2,
two ordered quartic sign sheets        S3/C3,
prime-three collapsed point            S3/S3.            (8)
```

It is a dictionary of `S_3`-sets, not a transfer of a physical observable.

## 3. The co-occurrence cover

For a chosen ratio `r`, define

```text
C_r={(sigma,sigma.r):sigma in S3}.                       (9)
```

Projection to the first coordinate is the regular six-state Farey/`V_4`
frame.  Projection to the second coordinate is the anharmonic orbit of `r`,
and every second fibre has size `|Stab(r)|`.  The exact multiplicities are

```text
generic 1,       harmonic/boundary 2,
equianharmonic 3, prime-three singleton 6.              (10)
```

Thus keeping only the ratio orbit forgets exactly the stabilizer phase.
Boundary and harmonic packets are isomorphic three-state `S_3`-sets: a coarse
three-state quotient cannot distinguish a trivial-leg boundary from a
transverse harmonic packet without a labelled sidecar.

## 4. The `PSL/PGL` reflection-sign boundary

The order-three transform and anharmonic reflection have matrices

```text
c=[[0,-1],[1,1]],        det(c)=1,
s=[[0, 1],[1,0]],        det(s)=-1.                     (11)
```

For odd `p`, projective determinant square class is well-defined, and

```text
H_p intersect PSL_2(F_p)=H_p   if p=1 (mod 4),
                         =C_3   if p=3 (mod 4).          (12)
```

The actual order-two generator inherited from
`PSL_2(Z)=C_2*C_3` acts by

```text
alpha(r)=-1/r,                                            (13)
```

not `s(r)=1/r`.  They alias only in characteristic two.  The minimal prime
three hostile is

```text
s(1)=1,                         alpha(1)=-1=2.            (14)
```

So the anharmonic reflection fixes the collapsed transverse point while the
modular involution sends it into the boundary orbit.  The common `C_3` is
robust; the binary reflection requires a determinant-square or orientation
sidecar.

This also prevents an action-level misidentification with THM-3034.  Its
`C_3` acts freely by translation on each genus-one sheet; the ratio `C_3`
fixes the equianharmonic `S_3/C_3` points.  The same abstract group has
opposite fixed-point geometry.

## 5. Exact Farey hostile: the frame is not the gate

Use the THM-2056 defect

```text
F_w(d)=||d||^2-91 w.d
```

and take

```text
w=(1,0),       u=(1,0),       v0=(88,1),       v1=(90,1).  (15)
```

Both `(u,v0)` and `(u,v1)` are acute determinant-one Farey frames, and their
ordered reductions modulo two agree.  Their endpoint signs also agree:

```text
F(u)=-90,             F(v0)=-263,             F(v1)=-89.
```

But the two children lie on opposite sides:

```text
F(u+v0)=-177,                         F(u+v1)=1.          (16)
```

Hence

```text
same regular S3/V4 frame + same endpoint gate bits
  does not determine the Farey child gate.                              (17)
```

The lost coordinate is the integral Gram/carry term in THM-2056.  This is
why the orbit diamond organizes the two tower grammars but does not close an
LRC or Keller inequality by itself.

## 6. Connection contract

```text
source:    an ordered level-two Farey frame and a transverse ratio;
map:       (sigma,r) -> (sigma,sigma.r);
target:    the regular/natural/parity/trivial S3 orbit ladder;
preserved: C2/C3 stabilizer type and labelled co-occurrence before projection;
destroyed: integral Gram term, determinant-square sign, affine origin,
           physical owner, and quartic sheet phase after coarse projection;
sidecars:  the THM-2056 Gram/carry coordinate and a PSL/PGL orientation bit;
hostiles:  (14), boundary-versus-harmonic aliasing, and (16).             (18)
```

No LRC row, Keller chart, physical current, or ledger decrement follows.

## 7. Exact companion

Run

```text
python 04-computation/modular_farey_anharmonic_quartic_orbit_diamond_thm3035.py
python -O 04-computation/modular_farey_anharmonic_quartic_orbit_diamond_thm3035.py
```

Both modes byte-match

```text
05-knowledge/results/modular_farey_anharmonic_quartic_orbit_diamond_thm3035.out.
```

The companion checks every orbit and stabilizer for primes
`2,3,5,7,11,13,17,19`, all four character rows, the co-occurrence fibres,
the `PSL/PGL` split, all six mod-two Farey frames, and the exact determinant,
acute-angle, endpoint, and child values in `(15)--(16)`.  It contains no
truth-bearing Python `assert`.
