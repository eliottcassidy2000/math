---
id: THM-3191
title: "Factorial-block exterior Clifford law and global carry Smith profile"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The normalized exterior transfer of one fixed factorial prime block is an
  extension-independent rank-two operator D satisfying D^3=s^2 Delta D.
  Through K blocks, the exact carry thickness and p-free factorial unit
  multiply D^K, giving the global-across-blocks p-primary state Smith type
  (1,p^H,p^H), exterior type (p^H,p^H,p^(2H)), and an explicit squared-scale
  return over Z_p.  No Smith claim at primes different from p is made.
audit: >
  The pure-integer companion checks 300 one-block parameter systems, 38,336
  nonzero quotient vectors in the complementary chart atlas, 456 carried
  local blocks, and 392 global fixed-parameter products.  Normal and optimized
  replay agree with the stored transcript.  An independent immutable audit
  rederived the extension cancellation, cubic law, two-chart cover, carried
  unit layer, p-primary determinantal divisors, and squared missing-wedge
  return, and separately replayed both modes against the stored output.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
depends_on:
  - THM-3185-iterated-factorial-frobenius-descent-and-witt-carry-reset-hierarchy
  - THM-3188-quadratic-character-pre-reset-holonomy-and-exterior-flag-rigidity
related:
  - THM-3186-full-exterior-continuant-path-convolution-and-cancellation-wall
  - THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary
script: 04-computation/factorial_block_exterior_clifford_carry_thm3191.py
output: 05-knowledge/results/factorial_block_exterior_clifford_carry_thm3191.out
script_sha256: b95da466c338c396ccd4e5e805a30c9070bd89189f1bed83ec7ea77e93fd023c
output_sha256: 77127ab63593b3264f982e2c52893013479d7cdc4a6b0a08c5c8751bfac26a58
hash_basis: LF-normalized bytes
---

# THM-3191 -- factorial-block exterior Clifford law and global carry Smith profile

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3185 gives the thickness of each factorial reset, while THM-3188 gives
the exact quadratic-character holonomy of the invertible pre-reset tail.  The
two pieces compose more rigidly than their separate Smith profiles suggest:
the unknown boundary-extension entries cancel, leaving one rank-two operator
with a cubic Clifford law.  Iterating it closes the global carry atlas.

## 1. Fixed prime-block notation

Fix an odd prime `p` and one pair of residues

```text
s!=0,                    Delta=1-4sv!=0 in F_p,
chi=Delta^((p-1)/2),     a=-chi.                            (1)
```

Let `T_n` be the exact factorial Gauss--Manin transfer of THM-3185, with
fixed integral parameters `d==s (mod p)` and `v`.  For `ell>=1`, define

```text
B_ell=T_(ell p-1)...T_((ell-1)p),
h_ell=v_p(ell p)=1+v_p(ell),
u_ell=ell/p^(v_p(ell)) in F_p^*.                            (2)
```

Thus `B_ell` contains one length-`p-1` invertible tail and the reset at
`N=ell p`.

## 2. The extension-free one-block operator

In the ordered exterior basis

```text
w_01=M wedge X,       w_02=M wedge D,       w_12=X wedge D, (3)
```

THM-3185 gives the reset layer

```text
C_0=
[0   2s^2            -s       ]
[0   -s               2sv     ]
[0   -s(2s+1)         s(2v+1)].                            (4)
```

THM-3188 writes the pre-reset exterior holonomy as

```text
A_tilde=
[1   a(beta+gamma)    -a beta]
[0        a               0  ]
[0       -a               a  ].                            (5)
```

The first column of `C_0` is zero.  Therefore both unknown extension entries
`beta,gamma` disappear from the product.  Entrywise multiplication gives

```text
D=C_0 A_tilde=a s E,                                       (6)

E=
[0    2s+1             -1  ]
[0   -(1+2v)            2v ]
[0   -2(s+v+1)         2v+1].                              (7)
```

At an ordinary block, `D` is exactly `p^(-1)Lambda^2B_1 mod p`.  More
generally the carried local layer is

```text
p^(-h_ell)Lambda^2B_ell==u_ell D.                          (8)
```

Indeed every pre-reset tail has the same reduction `(5)`, while division of
the reset matrix by `p^h_ell` instead of by the exact integer `ell p`
contributes precisely its unit part `u_ell`.

## 3. Cubic law and two-chart atlas

The lower-right `2x2` block of `E` is

```text
F=[-(1+2v)       2v;
   -2(s+v+1)     2v+1].                                    (9)
```

It has

```text
tr F=0,                       det F=4sv-1=-Delta,           (10)
```

so Cayley--Hamilton gives `F^2=Delta I`.  Since `E` has block form
`[0,*;0,F]`, this implies

```text
E^3=Delta E,                    D^3=s^2 Delta D.            (11)
```

Both `D` and `D^2` have rank two.  For every positive power they have the
same right kernel `span(w_01)` and left kernel `span(1,-1,1)`.  Explicitly,

```text
D^(2r+1)=(s^2 Delta)^r D,
D^(2r)  =(s^2 Delta)^(r-1)D^2             (r>=1).          (12)
```

There is also a minimal adaptive projective atlas.  The image of `D` lies in

```text
Y_01-Y_02+Y_12=0.                                            (13)
```

Its last two coordinates are `as F(Y_02,Y_12)^t`.  Since `F` is invertible,
they cannot vanish simultaneously on a nonzero quotient vector.  Hence the
two charts

```text
Y_02!=0,                         Y_12!=0                    (14)
```

cover `P(im D)`.  A selected visible coordinate may vanish, but its
complementary coordinate is then forced nonzero.  This is an exact
Pluecker-chart wall/transition statement for the factorial exterior layer;
no identification with a PRS-selected coordinate is asserted.

## 4. Global carry holotopy

Let

```text
P_K=B_K...B_1=T_(Kp-1)...T_0,
H_K=sum_(ell=1)^K h_ell=v_p((Kp)!)=K+v_p(K!),
U_K=prod_(ell=1)^K u_ell
   ==K!/p^(v_p(K!))                              (mod p).  (15)
```

Compound multiplicativity and `(8)` give the exact first global layer

```text
p^(-H_K)Lambda^2P_K==U_K D^K                    (mod p).   (16)
```

Thus the carry depth lives in `H_K`, the p-free factorial phase lives in
`U_K`, and the projective exterior state alternates only between `D` and
`D^2` according to the parity of `K`.  Equations `(12)` and `(16)` are a
closed two-state discrete holotopy rather than an unbounded family of new
tail flags.

The state determinant has valuation `2H_K`, and THM-3185's rank-one block
reduction gives one unit invariant factor.  Since the right side of `(16)`
has rank two, its two smallest exterior determinantal divisors have valuation
`H_K`.  Therefore, over `Z_p` (equivalently over the DVR `Z_(p)` for these
valuation statements),

```text
Smith_p(P_K)=(1,p^H_K,p^H_K),
Smith_p(Lambda^2P_K)=(p^H_K,p^H_K,p^(2H_K)).                (17)
```

This global-across-blocks p-primary Smith statement is stronger than merely
summing `p`-adic determinant lengths: it proves that the two transverse
`p`-thicknesses remain exactly equal after every carried block.  Other prime
divisors coming from the integral parameters, `Delta`, or intermediate
determinantal divisors are uncontrolled and are not part of `(17)`.

## 5. Global missing-wedge return

Each pre-reset tail preserves the integral line `w_01`, and the reset at
`N=ell p` sends that line to itself first at squared thickness.  Multiplying
the exact local returns gives

```text
p^(-2H_K)Lambda^2P_K(w_01)
  ==U_K^2(-Delta)^K w_01                         (mod p).   (18)
```

The coefficient is a unit.  Thus the third exterior invariant in `(17)` is
not only forced by determinant length; its oriented generator and exact first
nonzero residue are explicit.

## 6. Scope, boundaries, and connections

The fixed-parameter hypothesis is load-bearing.  The same `s,v,Delta` occur
in every block because one factorial system is being iterated.  Products of
heterogeneous operators with varying residues need not commute and are not
governed by `(11)--(16)`.  The discriminant wall `Delta=0` is excluded and is
already the sharp homogeneous/exterior collapse of THM-3188.

The two-chart atlas `(14)` says that vanishing of one selected exterior
coordinate is a chart wall, not exterior death.  To use it in THM-3186's
continuant or in a PRS Newton argument one still needs an explicit map from
the selected path coordinate to one of these two charts, plus an adaptive
permission to switch charts.  Visible path sums can cancel, and neither
Smith data nor `(11)` chooses the required projection.

THM-2624 is a useful hostile boundary for this last step: two signed charts
can give exact tomography without furnishing a physical carrier transport.
Here the two charts share one exterior image plane, but no proved map sends a
PRS-selected coordinate into that plane with a lawful chart switch.  The
atlas is therefore reconstruction data, not yet a selector theorem.

The odd/even alternation in `(12)` and the p-free unit cocycle in `(15)` are
formal analogues of central parity and carry towers elsewhere in the
repository.  No owner, root, current, or LRC object map is known.  No
`NC(2)`, `GMC(2)`, arbitrary-support, or `LRC(14)` conclusion is claimed.

## 7. Exact evidence

Run

```text
python 04-computation/factorial_block_exterior_clifford_carry_thm3191.py
python -O 04-computation/factorial_block_exterior_clifford_carry_thm3191.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
integer arithmetic only.  It checks all 300 off-wall parameter systems over
`p=3,5,7,11,13`; every one of 38,336 nonzero quotient vectors in their two
projective charts; 456 ordinary/carried local blocks through `ell=p+2`; and
392 exact global products through seven blocks for `p=3,5,7`.  The global
checks derive the `p`-primary state and exterior determinantal-divisor
valuations and the squared return from the full integer transfer product;
they make no claim about invariant factors at other primes.  There is no
floating point, random sampling, imported executable, or assertion-sensitive
test.

**QED.**
