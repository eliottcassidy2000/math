---
id: THM-3035
title: "Level-two Farey anharmonic quartic S3 orbit diamond"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.
  The anharmonic S3 action on P1(F_p) splits into boundary, harmonic,
  equianharmonic, and generic orbits, with exact stabilizer fibres and bad-prime
  degenerations.  The PSL/PGL reflection sign and a THM-2056 Farey hostile show
  precisely why this orbit dictionary does not determine an integral child
  gate.  No Keller, LRC(14), or physical-current consequence is asserted.
source: codex-modular-farey-anharmonic-orbit-diamond-2026-08-01
depends_on:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-2996-prime-modular-affine-defect-trichotomy-and-spherical-quartic-uniqueness
  - THM-2998-quartic-star-complement-cross-wall-real-locus-and-c3-x0-14-quotient
related:
  - THM-3034-ordered-quartic-cross-wall-x1-14-and-diamond-quotient
script: 04-computation/farey_anharmonic_s3_orbit_diamond_thm3035.py
output: 05-knowledge/results/farey_anharmonic_s3_orbit_diamond_thm3035.out
script_sha256: a2e7c8e79340023499353ef17419dd1ffa095b17fed5bf1b8a837c9f6eee6672
output_sha256: f5cf56eb89f8e25c2be7fdc823edb2d77c9c747f162f31dd129dd474e853c352
hash_basis: LF-normalized bytes
---

# THM-3035 -- level-two Farey / anharmonic / quartic `S_3` orbit diamond

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

## 1. Inheritance

- THM-2056 proves that the Farey child defect contains the exact Gram cross
  term `2u.v`; mod-two direction data alone is not a gate certificate.
- THM-2632 identifies a mod-two Farey frame with an ordering of the three
  nonzero directions of `V_4=F_2^2`; the six frames form the regular `S_3`
  set.
- THM-2984 gives the anharmonic transforms

  ```text
  s(r)=1/r,                  c(r)=-1/(1+r),
  ```

  and their abstract `S_3` action.
- THM-2996 identifies the unique quartic affine coincidence
  `AGL_2(F_2)=S_4`, hence `S_4/V_4=S_3`, while retaining the bad-prime
  boundary.
- THM-2998 supplies two ordered sign sheets with a free common `C_3` action.

The missing exact statement is not another abstract group isomorphism.  It is
an orbit-and-stabilizer dictionary together with the sign coordinate that is
lost when the modular involution is confused with the anharmonic reflection.

## 2. The orbit theorem

For a prime `p`, let

```text
H_p=<s,c> <= PGL_2(F_p),
s(r)=1/r,                   c(r)=-1/(1+r).
```

The three points `{infinity,0,-1}` are permuted faithfully, so
`H_p~=S_3`.
For every prime `p>=5`, its action on `P^1(F_p)` has the following named
orbits.

```text
boundary:       {infinity,0,-1}                  = S3/C2,
harmonic:       {1,-2,-1/2}                      = S3/C2,
equianharmonic: {r:r^2+r+1=0}                    = S3/C3,
generic:        every remaining orbit            = S3/1.       (1)
```

The equianharmonic orbit is present exactly when `p=1 mod 3`.  Consequently
the number of generic regular orbits is

```text
(p-7)/6  if p=1 mod3,
(p-5)/6  if p=2 mod3.                              (2)
```

The two bad primes are structural.

```text
p=2: P^1(F2) is the boundary orbit; no transverse point remains.
p=3: boundary={infinity,0,-1}, plus the singleton {1}=S3/S3.   (3)
```

At `p=3`, harmonic and equianharmonic conditions collapse at `1`.  This is
the exact prime-three anomaly complementary to the prime-two exhaustion.

### Proof of the classification

The formulas give

```text
s^2=c^3=1,                 scs=c^(-1).
```

Their action on `{infinity,0,-1}` is faithful, so `H_p` has order six.
For `p>=5`, the fixed locus of `<c>` is exactly
`r^2+r+1=0`.  Its points are the nontrivial cube roots of unity, so the
two-point orbit exists exactly for `p=1 mod3`.  The three transpositions
and their fixed pairs are

```text
r -> 1/r:          {-1,1},
r -> -1-r:         {infinity,-1/2},
r -> -r/(1+r):     {0,-2}.
```

Their fixed loci therefore comprise the disjoint boundary and harmonic
triples.  These exhaust all points with nontrivial stabilizer because every
proper nontrivial subgroup of `S_3` has order two or three.  Removing the
named orbits from the `p+1` points and dividing by six gives `(2)`.
Direct substitution gives the two degenerations in `(3)`.

At `p=13`, the four orbits are

```text
boundary        {0,12,infinity},
harmonic        {1,6,11},
generic         {2,4,5,7,8,10},
equianharmonic  {3,9}.                               (4)
```

Their permutation characters on `(identity, transposition, 3-cycle)` are

```text
regular six:  (6,0,0),
natural three:(3,1,0),
parity two:   (2,0,2),
trivial one:  (1,1,1).                               (5)
```

Thus the recurring repo carriers fit one exact representation ladder:

```text
six Farey/V4 orders                    regular S3,
three nonzero V4 directions/matchings  S3/C2,
two ordered quartic sign sheets        S3/C3,
prime-three collapsed point            S3/S3.          (6)
```

This is a dictionary of `S_3`-sets, not an automatic transfer of a physical
observable.

## 3. The co-occurrence frame

For a chosen ratio `r`, define

```text
C_r={(sigma,sigma.r):sigma in S3}.                   (7)
```

Projection to the first coordinate is the regular six-state Farey/V4 frame.
Projection to the second is the anharmonic orbit of `r`, and every second
fibre has size `|Stab(r)|`.  The multiplicities are therefore

```text
generic 1,       harmonic/boundary 2,
equianharmonic 3, prime-three singleton 6.            (8)
```

This is the precise free-factor co-occurrence carrier: keeping only the
second projection forgets exactly the stabilizer phase.  In particular,
boundary and harmonic packets are isomorphic three-state `S_3`-sets.  No
three-state quotient can distinguish a trivial-leg boundary from a genuinely
transverse harmonic packet without a labelled sidecar.

## 4. The `PSL/PGL` reflection-sign boundary

The order-three transform has matrix

```text
c = [[0,-1],[1,1]],             det(c)=1,
```

but the anharmonic reflection has

```text
s = [[0,1],[1,0]],              det(s)=-1.             (9)
```

For odd `p`, projective determinant square class is well-defined, and

```text
H_p intersect PSL_2(F_p) = H_p   if p=1 mod4,
                         = C_3   if p=3 mod4.           (10)
```

Indeed, rescaling a projective matrix changes its determinant by a square.
The even words in `s` have determinant-square class `1`, while every odd
word has class `-1`; and `-1` is a square precisely for `p=1 mod4`.

The actual order-two generator inherited from `PSL_2(Z)=C_2*C_3` is

```text
alpha(r)=-1/r,                                          (11)
```

not `s(r)=1/r`.  They alias only in characteristic two.  At `p=3`, for the
minimal hostile,

```text
s(1)=1,                  alpha(1)=-1=2,
```

so `s` fixes the transverse singleton while `alpha` sends it to the boundary.
The common `C_3` is robust; the binary reflection requires an orientation or
determinant-square sidecar.

This also explains the correct relation to THM-2998: its `C_3` is a free
translation action on each genus-one sheet, whereas the equianharmonic
`S_3/C_3` points are fixed by the ratio `C_3`.  The same abstract group has
opposite fixed-point geometry.  They must not be identified as actions.

## 5. Exact THM-2056 hostile: the frame is not the gate

Let

```text
F_w(d)=||d||^2-91 w.d,
w=(1,0),       u=(1,0),       v0=(88,1),       v1=(90,1).
```

Both pairs `(u,v0)` and `(u,v1)` are acute determinant-one Farey frames and
have the same ordered reduction modulo two.  Their endpoint signs also agree:

```text
F(u)=-90,        F(v0)=-263,        F(v1)=-89.
```

But their mediants lie on opposite sides:

```text
F(u+v0)=-177,                F(u+v1)=1.             (12)
```

The first failed implication is therefore exact:

```text
same regular S3/V4 frame + same endpoint gate bits
  does not determine the Farey child gate.
```

The lost coordinate is the integral Gram/carry term in THM-2056.  This is
the sharp reason the orbit diamond organizes the two tower grammars but does
not close an LRC or Keller inequality by itself.

## 6. Connection contract

```text
source:    ordered mod-two Farey frame and a transverse ratio;
map:       (sigma,r) -> (sigma,sigma.r);
target:    regular/natural/parity/trivial S3 orbit ladder;
preserved: C2/C3 stabilizer type and labelled co-occurrence before projection;
destroyed: integral Gram term, determinant-square sign, affine origin,
           quartic sheet phase after coarse projection;
sidecars:  THM-2056 Gram/carry coordinate and a PSL/PGL orientation bit;
hostiles:  p=3 alpha-versus-s, boundary-versus-harmonic alias, and (12).
```

The canonical companion

```text
python 04-computation/farey_anharmonic_s3_orbit_diamond_thm3035.py
python -O 04-computation/farey_anharmonic_s3_orbit_diamond_thm3035.py
```

checks every orbit and stabilizer for primes `2,3,5,7,11,13,17,19`, the
character and co-occurrence tables, the `PSL/PGL` split, all six mod-two
Farey frames, and the exact hostile `(12)`.

The candidate LF-normalized hashes are

```text
script  a2e7c8e79340023499353ef17419dd1ffa095b17fed5bf1b8a837c9f6eee6672
output  f5cf56eb89f8e25c2be7fdc823edb2d77c9c747f162f31dd129dd474e853c352
```

The finite companion uses explicit `require` checks rather than truth-bearing
Python `assert` statements.  It is evidence for the finite controls; the
all-prime quantifiers are proved symbolically above.  Promotion awaits an
independent proof and hostile audit.
