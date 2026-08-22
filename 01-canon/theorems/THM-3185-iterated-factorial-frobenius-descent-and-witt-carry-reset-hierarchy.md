---
id: THM-3185
title: "Iterated factorial Frobenius descent and Witt-carry reset hierarchy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The factorial Gauss--Manin state descends through arbitrarily many prime
  blocks: z_(Kp+a)^[d]=s^K z_a^[s] modulo p when d=s modulo p.  At the
  ell-th reset, h=1+v_p(ell) is exactly the thickness of both transverse
  directions in the x-weighted lattice.  Scalar, exterior, gauge, and
  discriminant-wall layers are explicit.
audit: >
  The pure-integer companion checks 720 complete multiblock state descents,
  744 full rank-one block-propagator identities, the resonant-pair corollary,
  reset determinantal divisors at ordinary and carried block numbers,
  exact compound-matrix reconstruction and both normalized exterior layers,
  and discriminant walls of orders one and two.  Normal and optimized replay
  agree with the stored transcript.  Two independent immutable audits
  rederived the state and matrix descents, every Smith/gauge/wall formula,
  both exterior layers, and the global ledger.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
depends_on:
  - THM-3182-factorial-gauss-manin-rank-one-reset-and-two-transverse-smith-bands
related:
  - THM-3148-fixed-offset-frobenius-endpoint-resultant-classification
  - THM-3178-squarefree-resultant-tangent-cone-and-first-witt-norm
script: 04-computation/factorial_multiblock_witt_carry_reset_thm3185.py
output: 05-knowledge/results/factorial_multiblock_witt_carry_reset_thm3185.out
script_sha256: 7b30fe78ca9d93bd843a6d0b110396650723c13c2aa2578e87424cf259e3b733
output_sha256: 1c0b963c5e4eb15f213829feb67d2bb6968a264a356f429c08560ca68bcab407
hash_basis: LF-normalized bytes
---

# THM-3185 -- iterated factorial Frobenius descent and Witt-carry reset hierarchy

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3182 isolates the first prime reset.  The same mechanism repeats at every
prime block, but not with constant thickness: the block number contributes
its own `p`-adic valuation.  This produces a literal carry hierarchy in the
integral Gauss--Manin lattice.

## 1. Arbitrary-block state descent

Use THM-3182's exact notation

```text
q_d(x)=d-x+v x^2,
M_n^[d]=int_0^infinity q_d(x)^n e^(-x) dx,
X_n^[d]=int_0^infinity x q_d(x)^n e^(-x) dx,
z_n^[d]=(M_n^[d],X_n^[d],d^n)^t.                            (1)
```

Let `p` be prime, let `d==s (mod p)` with

```text
1<=s<p,                                                     (2)
```

and write

```text
n=Kp+a,                    K>=0,             0<=a<p.        (3)
```

Then, coefficientwise in `v`,

```text
z_(Kp+a)^[d](v)==s^K z_a^[s](v)                  (mod p).  (4)
```

### Proof

Modulo `p`, `q_d=q_s` and

```text
q_s^(Kp+a)=((q_s)^p)^K q_s^a.                               (5)
```

Every nonconstant term of `(q_s)^p` has `x`-degree at least `p`.
The factorial functional sends `x^j` to `j!`, so every contribution using
one such term vanishes modulo `p`; only the constant `(s^p)^K==s^K`
survives.  This proves the `M` coordinate.  For `X`, the monomial value is
`(j+1)!`, and every new degree is at least `p`, so the same argument applies.
The third coordinate is Fermat:

```text
d^(Kp+a)==s^(K+a)==s^K s^a.                                (6)
```

Equivalently, THM-3182's transfer satisfies
`T_(jp+r)^[d]==T_r^[s]`; each rank-one reset contributes one additional
factor `s`.  This gives a second derivation of `(4)`.

The matrix statement is stronger and makes the reset mechanism explicit.  If

```text
B_s=T_(p-1)^[s] ... T_0^[s],
e_D=(0,0,1)^t,                 1=(1,1,1)^t,
```

then

```text
B_s==s 1 e_D^t,                                               (6a)
T_(Kp+a-1)^[d] ... T_0^[d]
   ==s^K z_a^[s] e_D^t             (K>=1, mod p).            (6b)
```

Indeed the last transfer in a block is `s 1 e_D^t`, while the third row of
every preceding transfer is `s e_D^t`; hence
`e_D^t T_(p-2)...T_0=s^(p-1)e_D^t=e_D^t`.  Also
`T_(a-1)...T_0 1=z_a^[s]`.  Thus a complete block erases both transverse
input coordinates and remembers only the `D` coordinate, multiplied by `s`.
The unit case `s=1` still satisfies `(4)` and `(6)`; only the pair indexing in
the next section degenerates.  The nonunit case `s=0` is a separate boundary
and is not asserted here.

## 2. Arbitrary-quotient resonant-pair corollary

Suppose

```text
d=Kp+s,               K>=1,               2<=s<p.          (7)
```

Taking `a=s-2,s-1` in `(4)` gives

```text
M_(d-2)^[d]==s^K M_(s-2)^[s],
M_(d-1)^[d]==s^K M_(s-1)^[s]                    (mod p).   (8)
```

If `p>2(s-1)`, the two small rows retain degrees `s-2,s-1` modulo `p`.
Thus `(8)` is the complete height-zero residual pair

```text
s^K F_(s-2,s),                  s^K F_(s-1,s).              (9)
```

In the notation of THM-3148, if `p` does not divide `delta_s`, `(9)` has no
common affine root.  Hence the large resonant pair has no common slope-zero
factor over `Q_p`.  Unlike THM-3148's first-block formulation, this conclusion
allows every quotient `K`.

This is only a zero-face result.  Repeated resets create additional positive
valuation layers; `(8)` alone does not exclude their common factors.

## 3. Reset thickness at every block

Let

```text
N=ell*p,                  ell>=1,
h=v_p(N)=1+v_p(ell),
Delta=1-4dv.                                                  (10)
```

At the block boundary `n=N-1`, THM-3182's `x`-weighted transfer is

```text
G_(N-1)=
[-N                 2vN                    d       ]
[-N(2N+1+2d)        N(1+2v(2N+1))         d(2N+1)]
[0                   0                      d       ],        (11)
```

with

```text
det G_(N-1)=-dN^2 Delta.                                    (12)
```

Work over `Z_p`, or an unramified `p`-adic DVR with valuation normalized by
`v_p(p)=1`, and assume first that `d Delta` is a unit.  The entry `d` is a
unit, every `2x2` minor is divisible by `N`, and the minor on rows `(1,3)` and
columns `(1,3)` is exactly `-Nd`.  Therefore

```text
Smith(G_(N-1))=(1,p^h,p^h)                                  (13)
```

up to units.  The carry depth `v_p(ell)` thickens both transverse directions
by exactly that amount.

## 4. Scalar companion and the adjacent lattice modification

For the scalar state

```text
Y_n=(M_n,M_(n-1),d^n)^t,                                   (14)
```

the reset matrix is

```text
S_(N-1)=
[2N(2N-1)v    (N-1)N Delta    d-N]
[1             0                0  ]
[0             0                d  ].                       (15)
```

Its `2x2` minor on rows `(1,2)` and columns `(1,3)` is `N-d`, a unit, while

```text
det S_(N-1)=-(N-1)N Delta d.                                (16)
```

Hence

```text
Smith(S_(N-1))=(1,1,p^h).                                  (17)
```

For odd `p` with `v` a unit, the exact gauge from THM-3182 has determinant

```text
det P_n=n Delta/(2v).                                       (18)
```

Thus `P_(N-1)` is a unit gauge and `P_N` has index `p^h`.  The determinant
length identity is

```text
2h=h+h-0.                                                   (19)
```

So the second thick direction in `(13)` is the adjacent output-lattice
modification.  It is real in the specified `X` lattice, but it is not a
rational-conjugacy invariant.

## 5. Exterior-square carry layer

In the basis `(M wedge X,M wedge D,X wedge D)`, put
`C_N=Lambda^2 G_(N-1)`.  Exact minors give

```text
C_N=
[-N^2 Delta   2Nd^2                -Nd                    ]
[0             -Nd                   2Nvd                  ]
[0             -Nd(2N+1+2d)         Nd(1+2v(2N+1))].       (20)
```

Every entry is divisible by the exact integer `N`.  After division by `N`
and reduction modulo `p`, with `dbar=s`,

```text
N^(-1)C_N ==
[0  2s^2       -s       ]
[0  -s          2sv     ]
[0  -s(2s+1)    s(2v+1)]                         (mod p).   (21)
```

This universal first layer is independent of `ell`.  It has rank two, right
kernel `M wedge X`, left kernel `(1,-1,1)`, and all three nonzero `2x2` minors
equal to

```text
-s^2(1-4sv).                                                (22)
```

The missing direction returns at the exact squared scale:

```text
N^(-2)C_N(M wedge X)==-Delta(M wedge X)           (mod p). (23)
```

Consequently `(13)` induces

```text
Smith(Lambda^2 G_(N-1))=(p^h,p^h,p^(2h)).                  (24)
```

Normalizing by `N`, rather than merely by `p^h`, is load-bearing in
`(21)--(23)`: the quotient `N/p^h` is a residue-field unit which otherwise
rephases the displayed layer.

## 6. Complete discriminant-wall valuation

More generally retain `d` as a `p`-adic unit, suppose `Delta!=0`, and put

```text
t=v_p(Delta)>=0.                                            (25)
```

The same exact minor and determinant calculations give

```text
Smith(G_(N-1))=(1,p^h,p^(h+t)),
Smith(S_(N-1))=(1,1,p^(h+t)),
Smith(Lambda^2G_(N-1))=(p^h,p^(h+t),p^(2h+t)).              (26)
```

For the gauge in `(18)`, the input and output indices have valuations

```text
t,                         h+t,                              (27)
```

and the determinant identity becomes

```text
2h+t=(h+t)+(h+t)-t.                                        (28)
```

Thus a discriminant wall thickens only the final invariant direction after
the carry thickness has been accounted for.  If `Delta=0` exactly, the final
invariant factor is zero and no finite Smith exponent is asserted.

## 7. Global carry ledger and block monodromy

The local thicknesses telescope to a closed global formula.  Assume `d` and
`Delta` are units and compose the transfers from `n=1` through `n=Kp-1`.
For the `x`-weighted lattice, THM-3182's transfer determinant gives total
length

```text
L_X=2 v_p((Kp)!)=2(K+v_p(K!)).                              (29)
```

The last equality is Legendre's formula.  Equivalently,

```text
v_p((Kp)!)=sum_(ell=1)^K [1+v_p(ell)],                      (30)
```

so `(29)` is exactly twice the sum of the reset thicknesses in `(10)`.

For the scalar companion, `(16)` gives

```text
L_scalar=v_p((Kp-1)!)+v_p((Kp)!)
        =2v_p((Kp)!)-v_p(Kp).                               (31)
```

The endpoint gauge contributes

```text
v_p(det P_(Kp))-v_p(det P_1)=v_p(Kp),                       (32)
```

Here the gauge comparison assumes, as in `(18)`, that `p` is odd and `v` is
a unit.  The determinant ledgers `(29)` and `(31)` themselves do not need
that gauge hypothesis.  Under the gauge hypothesis, `(31)` and `(32)` give

```text
L_X=L_scalar+v_p(Kp).                                       (33)
```

This is the global form of the local length identity `(19)`: adjacent scalar
singular transfers and output-lattice modifications telescope to the doubled
weighted carry ledger.

There is an independent unit-colour coordinate.  The multiplier in `(4)` is

```text
s^K in F_p^*.                                               (34)
```

If `s` is primitive, it traverses all `p-1` unit colours with period `p-1`,
while `1+v_p(ell)` records the independent `p`-adic odometer/carry depth.  At
`p=13` this is an exact `C_12` unit-colour times `13`-adic carry skeleton.
This is a precise formal parallel to the colour/carry towers used on the LRC
side of the repository, but there is no object map, owner map, or reduction
between the two problems.

## 8. Scope and next bridge

The theorem gives an exact multiblock/carry atlas for the specified integral
states.  It does not determine the Newton polygon of the moment coordinate:
the non-unimodular scalar/weighted gauge and the oriented finite-tail
projection are still load-bearing.  In particular, it does not prove the
generic `floor(s/2)` Euclidean-depth staircase.

Equation `(9)` suggests a new arithmetic sieve: any prime divisor of `d-s`
with `p>2(s-1)` and `p` not dividing `delta_s` kills the slope-zero face.
Positive slopes, common residual selectors, and carried discriminant walls
must still be handled before this becomes an all-degree result.  No `NC(2)`,
`GMC(2)`, arbitrary-support, or `LRC(14)` conclusion is claimed.

## 9. Exact evidence

Run

```text
python 04-computation/factorial_multiblock_witt_carry_reset_thm3185.py
python -O 04-computation/factorial_multiblock_witt_carry_reset_thm3185.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
integer arithmetic only.  It reconstructs direct `M,X,D` rows for 720
multiblock tests, verifies 744 complete rank-one block propagators, checks the
resonant-pair specialization, and verifies ordinary, carried, and twice-carried
reset thicknesses.  It reconstructs the exact compound transfer from `(11)`,
checks its divisibility by `N`, derives `(21)`, and separately derives the
`N^2`-normalized missing-column return `(23)`.  It also checks both finite
discriminant-wall orders from the actual weighted, scalar, and exterior
determinantal divisors, the global Legendre/telescoping ledger, and primitive
block-colour cycles.  There is no floating point, random sampling, imported
executable, or assertion-sensitive test.

**QED.**
