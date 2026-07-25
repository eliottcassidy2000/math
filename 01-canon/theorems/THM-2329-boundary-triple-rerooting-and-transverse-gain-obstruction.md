---
id: THM-2329
title: "Boundary-triple rerooting and the transverse-gain obstruction"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. Re-rooting
  the nonzero THM-2327 word/deepest-comb/bare Fourier triangle at its
  three legs gives the same nonzero monomial with projective shallow
  input labels [1:0], [0:1], and [1:-1]. Under THM-2321's canonical
  label map these are the pure-a axis, pure-b axis, and distinguished
  mixed gain -1. They are exactly the S_3 orbit on which one of the
  three character legs is trivial. Each of the other eleven directions
  in P^1(F_13) requires all three shallow characters to be nonzero, so
  no same-root c_3 edge can realize it by rerooting. This is a
  projective-label and coefficient-incidence theorem, not target
  polarization or relation-address incidence. No scalar profile is
  excluded and LRC(14) remains open.
source: codex-2026-07-25-boundary-triple-rerooting
depends_on:
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2321-prescribed-root-character-bispectrum-slice-positivity
  - THM-2327-two-colour-marked-unit-c3-triangle
related:
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2325-prescribed-target-gain-full-lattice-91-unit-needle-bank
script: 04-computation/lrc14_boundary_triple_rerooting_thm2329.py
output: 05-knowledge/results/lrc14_boundary_triple_rerooting_thm2329.out
script_sha256: f126dfd7e16576cb21453b529e727995d0e78d12882f2c116d7927b55fe73504
output_sha256: 1e7150fc0d04720f78c89b4a6b080853ac0985110d0405673678bb79162e40ff
hash_basis: working-tree bytes (LF)
---

# THM-2329 -- boundary-triple rerooting and transverse gain

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2327 produces an actual nonzero heterogeneous Fourier triangle

```text
X+B=Y,

(1_(E_Q))_hat(X)
(1_(D_(c_3)))_hat(B)
conjugate((1_E)_hat(Y))
 !=0,                                               (1)
```

where

```text
B=m c_3,
gcd(m,91)=1,

nu_13(X)=nu_13(Y)=lambda<nu_13(c_3),                (2)
```

and `X,Y` have the same prescribed nonzero shallow root character. The
first reading of (1) is an axis: the deep leg is trivial modulo thirteen
at the shallow scale. But an additive triangle has no preferred output.
Re-rooting it at the deep leg produces the distinguished mixed direction
`[1:-1]`.

This small operation identifies exactly how far the new analytic incidence
reaches in the target-gain corolla:

```text
same coefficientwise triangle
  -> two pure-axis readings
  -> one mixed gain -1 reading
  -> eleven genuinely transverse directions still missing.           (3)
```

## 1. The projective boundary triple

Let `p` be an odd prime. For a nonzero zero-sum character triple

```text
r+s+t=0                         in F_p,             (4)
```

choose two legs as inputs and the negative of the third as output. The
input pair has a projective label in `P^1(F_p)`.

One of the three character legs in (4) is zero exactly for the three
directions

```text
B_p={
  [1:0],
  [0:1],
  [1:-1]
}.                                                   (5)
```

Indeed, the first two points are `s=0` and `r=0`; the third is
`t=-(r+s)=0`. Conversely, each equality gives the displayed point.

Permuting the three legs acts on the projective input direction by the
usual anharmonic `S_3` action. Thus (5) is one complete `S_3` orbit: it
is the orbit of a triangle with one trivial character leg.

The complement has

```text
|P^1(F_p) minus B_p|=(p+1)-3=p-2                  (6)
```

directions. A direction lies in this complement exactly when

```text
r!=0,
s!=0,
r+s!=0.                                             (7)
```

These are the **transverse directions**: all three character legs are
nontrivial. At `p=13`, there are exactly

```text
13-2=11                                             (8)
```

of them.

## 2. Exact Fourier rerooting

Use the Fourier convention

```text
f_hat(n)=integral f(x)exp(-2*pi*i*n*x)dx.
```

For real functions `f,d,e` and integers `x,b,y` with `x+b=y`, define

```text
B_(f,d;e)(x,b)
 =f_hat(x)d_hat(b)conjugate(e_hat(y)).              (9)
```

Re-root at the old first input. Since `b+(-y)=-x`, realness gives

```text
B_(d,e;f)(b,-y)
 =d_hat(b)e_hat(-y)conjugate(f_hat(-x))
 =B_(f,d;e)(x,b).                                  (10)
```

Re-root at the old output. Since

```text
x+(-y)=-b,
e_hat(-y)=conjugate(e_hat(y)),
d_hat(-b)=conjugate(d_hat(b)),
```

one also has

```text
B_(f,e;d)(x,-y)
 =f_hat(x)e_hat(-y)conjugate(d_hat(-b))
 =B_(f,d;e)(x,b).                                  (11)
```

Thus all three expressions are literally the same complex monomial, not
three sums whose phases might cancel.

Apply this to (1):

```text
f=1_(E_Q),
d=1_(D_(c_3)),
e=1_E,

x=X, b=B, y=Y.                                     (12)
```

Every factor is real-valued. Equations (1), (10), and (11) prove that all
three rerootings are nonzero.

## 3. Their exact shallow-character labels

Divide the frequencies by `13^lambda` and reduce modulo thirteen. Put

```text
kappa=X/13^lambda mod 13.                           (13)
```

Equation (2) gives

```text
kappa!=0,
B/13^lambda=0 mod 13,
Y/13^lambda=kappa mod 13.                          (14)
```

The three rerootings (9)--(11) therefore have input-character pairs

```text
(X,B):       (kappa,0)       -> [1:0],
(B,-Y):      (0,-kappa)      -> [0:1],
(X,-Y):      (kappa,-kappa)  -> [1:-1].             (15)
```

They exhaust `B_13`. No choice of multiplier `m`, gauge, owner grade, or
root character changes this projective set: the load-bearing fact is
that a `c_3` shift is trivial at every shallower thirteen-adic scale.
Here “trivial” always means trivial **shallow character after division by
`13^lambda` and reduction modulo thirteen**. The ordinary frequencies
`B` and `-B` are nonzero, and their nonzero deepest-comb coefficients are
load-bearing in (1).
All six ordered rerootings give these three directions twice. On a
general slope `g`, the six slopes are

```text
g, 1/g, -(1+g), -1/(1+g), -g/(1+g), -(1+g)/g,
```

with the usual limiting interpretation at `0,infinity,-1`. For the
trivial-leg triangle this orbit collapses to exactly `B_13`.

## 4. The target-gain label map

THM-2321 identifies the root-character projective line with THM-2315's
two-target quotient by

```text
Phi([r:s])=[r e_a+s e_b].                           (16)
```

Consequently (15) maps to

```text
[1:0]   -> pure a-axis,
[0:1]   -> pure b-axis,
[1:-1]  -> mixed gain -1.                          (17)
```

The function-role typing is:

```text
[1:0]:   word + deep  -> bare,
[0:1]:   deep + bare  -> word,
[1:-1]:  word + bare  -> deep.                     (17a)
```

The third point is an honest mixed label, not a collapsed axis. It is
also fixed by target reversal:

```text
(-1)^(-1)=-1.                                      (18)
```

Thus the mixed direction forced by a trivial-output rerooting is exactly
one of THM-2315's unavoidable tournament ties. The other inversion-fixed
gain `+1` is already transverse: its character triple is proportional to

```text
(1,1,-2),
```

whose three entries are nonzero in `F_13`.

Equations (15)--(17) are an exact projective-label incidence statement for
the nonzero coefficient (1). They do **not** say that the analytic legs
are separately polarized by target `a` and target `b`. The function legs
remain word source, deepest danger comb, and bare source. Rerooting changes
which function occupies each input/output role; only a direction compatible
with the actual word `Q` has even support-level target meaning.
In particular, none of (17a) is THM-2321's homogeneous word current
`M_r M_s conjugate(M_(r+s))`.

## 5. The sharp transverse obstruction

Let `g in F_13^*` be a mixed target gain. A representative input pair is

```text
[1:g],
```

and the output character is `1+g`. Therefore:

```text
g=-1:
  the output character is trivial and THM-2327 supplies the boundary
  rerooting (15);

g!= -1:
  the two inputs and the output are all nontrivial.                  (19)
```

It follows that none of the other eleven projective directions can be
obtained by rerooting any same-root `c_3` edge. Replacing one such edge by
another, increasing its multiplier, or adding more vertices in the same
root-character graph never changes (14).

The missing analytic object is now exact:

```text
a coefficientwise nonzero triangle

  f_hat(A) g_hat(B) conjugate(h_hat(A+B))

with all three shallow root characters nonzero,
retaining the positive word and a lawful target polarization.        (20)
```

THM-2321 proves positive aggregate currents in every transverse
projective label. THM-2325 proves enormous exact relation-address banks
in every corresponding target-vector fibre. What remains is to make one
address and one target-polarized coefficientwise triangle coincide.

## 6. Tournament, Fano, and information loss

The faithful finite carrier here is

```text
(zero-sum character triple; chosen output; projective input label;
 exact three function legs; complex monomial).                       (21)
```

A tournament head selector loses two essential facts:

- the additive triangle can be rerooted, so no output orientation is
  intrinsic before the function roles are retained;
- the only mixed boundary point reached in (17) is reversal-fixed and
  must be a tie in every target-swap-equivariant tournament shadow.

The `S_3` above is the anharmonic subgroup of `PGL_2(F_13)` that permutes
the three boundary points. It is not the full projective group and is not
a proved physical target symmetry. Only the input swap is automatically
the target swap from THM-2315.

Likewise, the three points in (5) are not a Fano line supplying sevenfold
incidence. They are the trivial-leg orbit in a projective **line**. The
other eleven points are not generated by repeating this boundary orbit;
they require the transverse condition (7). Any Fano or `chi_7` probe that
forgets which character leg is zero therefore erases the exact obstruction.

No scalar profile is excluded. The scope is exactly the positive
shallow-owner word strata on which THM-2327 applies; the repeated-first
and alternative resonance branches remain outside that theorem. LRC(14)
remains open.

## 7. Exact verification

The companion enumerates all `14` projective directions over `F_13`, all
`36` nonzero-root rerootings, and all `168` nonzero character pairs. It
checks that precisely `132` vector pairs, representing precisely `11`
projective directions, have all three legs nonzero. It independently
checks the three rerooted monomials in exact Gaussian-rational arithmetic
and verifies shallow triviality/root preservation on `244440` graded
`c_3` shifts.

Reproduce with

```bash
python3 04-computation/lrc14_boundary_triple_rerooting_thm2329.py
python3 -O 04-computation/lrc14_boundary_triple_rerooting_thm2329.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_boundary_triple_rerooting_thm2329.out
```

byte-for-byte after LF normalization.
