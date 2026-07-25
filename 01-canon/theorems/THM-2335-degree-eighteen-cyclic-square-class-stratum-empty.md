---
id: THM-2335
title: "Degree-eighteen cyclic square-class stratum is empty"
status: >
  PROVED + VERIFIED-EXACT. For the structured quartic P and sextic Q of
  the degree-eighteen Jacobian spectral cone, 4P(y)^3+49Q(y)^2 is a
  square in C[y] if and only if B=C=D=W=0. Consequently the projective
  cyclic Hurwitz stratum H=1 isolated by THM-2332 is empty, including
  every support, not merely the residual support-three/four locus.
  This does not eliminate the H_2 S_5^2 or H_4 S_4^2 strata, and JC(2)
  remains open.
source: codex-2026-07-25-cyclic-square-class-elimination
depends_on:
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
related:
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2314-degree-eighteen-bd-ratio-bank-closure
  - THM-2316-degree-eighteen-cd-ratio-bank-closure
  - THM-2320-degree-eighteen-dw-ratio-bank-closure
  - THM-2324-degree-eighteen-bc-rational-ratio-closure
  - THM-2328-degree-eighteen-bw-ratio-bank-closure
script: 04-computation/jc2_degree18_cyclic_square_class_thm2335.py
output: 05-knowledge/results/jc2_degree18_cyclic_square_class_thm2335.out
script_sha256: 06edb041dda9640a3b495a90c2ca7d02298b79b34b82bcefc7525ed937aed1cd
output_sha256: c108c4a29d86b12baff3923d9e256f484588835103af5dad02833b1aa038130e
hash_basis: working-tree bytes (LF)
---

# THM-2335 -- the cyclic square-class stratum is empty

**PROVED + VERIFIED-EXACT.**

Use the structured polynomials from THM-2332:

```text
P
 =245y^4+1890By^2-24300B^2+122472D,

Q
 =539y^6+11340By^4+183708Cy^3
  +(72900B^2-367416D)y^2
  +(2361960BC+2480058W)y,

F=4P^3+49Q^2.                                      (1)
```

Then

```text
F is a square in C[y]

       iff

(B,C,D,W)=(0,0,0,0).                               (2)
```

In particular, the weighted-projective parameter cone has no cyclic
square-class point. The `H=1` branch in THM-2332 is empty.

## 1. The square root is forced from infinity

The leading coefficient of `F` is the parameter-independent number

```text
73060029=69*1029^2.                                (3)
```

Suppose `F=S^2`. Choose `r in C` with `r^2=69` and change the sign of
`S` if necessary. Then `S=rT`, where `T` has degree six and leading
coefficient `1029`.

Comparing coefficients successively from `y^11` down through `y^6`
uniquely forces

```text
T
 =1029y^6
 +(317520/23)B y^4
 +(1571724/23)C y^3
 -(20412/529)(1825B^2-12558D)y^2
 +(91854/529)(8060BC+5313W)y
 -(36741600/12167)(65B^3-69BD-3105C^2).            (4)
```

There is no `y^5` term. The recursion is legitimate because at every
stage the new coefficient is multiplied by the same nonzero number
`2*69*1029`.

After (4), `F-69T^2` has degree at most five. Its six coefficients are
nonzero rational constants times the following primitive equations:

```text
E5 =
 1560B^2C+805BW-9016CD,

E4 =
 617375B^4-6609510B^2D+7172550BC^2
 +22495725CW+21595896D^2,

E3 =
 233080B^3C+192395B^2W-1302168BCD
 -637560C^3-1088682DW,

E2 =
 1199915000B^5-12521098800B^3D+23960664000B^2C^2
 +41362139700BCW+34086322368BD^2
 -65507551200C^2D+25352682075W^2,

E1 =
 (8060BC+5313W)(65B^3-69BD-3105C^2),

E0 =
 -105225921875B^6+1530475458750B^4D
 +445024125000B^3C^2-7669002612600B^2D^2
 -472410225000BC^2D-10629230062500C^4
 +12875106064968D^3.                               (5)
```

Thus (5) is necessary and sufficient for `F` to be a square. The proof
below needs only `E5,E4,E3,E2,E1`; `E0` is retained both for the exact
equivalence and as an independent checksum.

## 2. The face B=0 collapses to the apex

Set `B=0`. Then

```text
E5=-9016CD,                                        (6)
```

so `C=0` or `D=0`.

If `C=0`, equation `E4=0` gives

```text
21595896D^2=0,
```

and then `E2=0` gives

```text
25352682075W^2=0.
```

If instead `D=0`, equation `E4=0` gives

```text
22495725CW=0.
```

The subcase `C=0` is already finished. If `C!=0`, then `W=0`, while
`E3=0` becomes `-637560C^3=0`, a contradiction. Hence

```text
B=0 and F square  =>  C=D=W=0.                    (7)
```

## 3. Normalize the chart B!=0

The parameters have weights

```text
wt(B,C,D,W)=(2,3,4,5).                             (8)
```

For every nonzero `lambda`, putting

```text
(B',C',D',W')=(lambda^2 B,lambda^3 C,
               lambda^4 D,lambda^5 W)
```

gives

```text
F(lambda y;B',C',D',W')=lambda^12 F(y;B,C,D,W).
                                                               (9)
```

Therefore squareness and support are invariant under this action. Over
`C`, choose `lambda` so that `B'=1`.

With `B=1`, equation `E5=0` uniquely gives

```text
W=8C(1127D-195)/805.                              (10)
```

After (10), the next two useful residual equations become

```text
E4 =
 251952120C^2D-36421650C^2
 +21595896D^2-6609510D+617375,                    (11)

E3 =
 -(8/5)C(398475C^2+7620774D^2
         -1851500D+87350).                         (12)
```

### 3.1 The subface C=0

If `C=0`, equation (10) gives `W=0`. Equations `E4=E2=0` are then the
two quadratics

```text
21595896D^2-6609510D+617375=0,

34086322368D^2-12521098800D+1199915000=0.          (13)
```

Their exact resultant is

```text
14658253353277217122876416000000 !=0.              (14)
```

They have no common complex root.

### 3.2 The subchart C!=0

After (10), `E1=0` splits as

```text
(4/5)C(74382D-2795)
   * (65-69D-3105C^2)=0.                           (15)
```

Since `C!=0`, there are two branches.

On the first branch,

```text
D=2795/74382.                                      (16)
```

Equation (11) and the parenthesized factor of (12) demand,
respectively,

```text
C^2= 14418035/972766179,

C^2=-28025/391314.                                 (17)
```

The difference between these two numbers is

```text
20348640145/235409415318 !=0.                      (18)
```

On the second branch,

```text
C^2=(65-69D)/3105.                                 (19)
```

Substituting (19) into (11) and the parenthesized factor of (12) gives,
up to nonzero constants,

```text
f(D)=3199392D^2-105156D-29015,

g(D)=22862322D^2-5581065D+287075.                  (20)
```

Their exact resultant is

```text
Res_D(f,g)=-466513778248576599586500 !=0.           (21)
```

So this branch also has no complex point. Equations (14), (18), and
(21) exhaust the normalized `B!=0` chart.

## 4. Converse and frontier effect

At the affine apex,

```text
F(y;0,0,0,0)=69(1029y^6)^2,                       (22)
```

so the converse in (2) holds. The apex is absent from the
weighted-projective cone. Therefore the cyclic branch

```text
4P^3+49Q^2=S_6^2                                  (23)
```

is empty on every support.

Combining this theorem with THM-2332 changes the residual
degree-eighteen Jacobian classification from three Hurwitz strata to
exactly two:

```text
mixed branch:
  4P^3+49Q^2=H_2 S_5^2,       H_2 squarefree;

four-simple branch:
  4P^3+49Q^2=H_4 S_4^2,       H_4 squarefree.      (24)
```

The discarded branch was the alternating-monodromy case: square
discriminant would force the connected cubic cover's monodromy into
`A_3`, hence make it cyclic with signature `(3,3)`. Equations (5)--(21)
show that the structured Jacobian coefficient cone never realizes this
apparently admissible permutation-theoretic object. This is the precise
place where the coefficient geometry is stronger than the abstract
tournament-like sheet relation.

No point of either branch in (24) is excluded here. The Faber,
one-form, flux, deck, other Newton-edge, `JC(2)`, and `DC(2)` problems
remain open.

## 5. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_cyclic_square_class_thm2335.py
python3 -O 04-computation/jc2_degree18_cyclic_square_class_thm2335.py
```

Both transcripts are byte-identical to the stored output. The companion
reconstructs (4) recursively rather than assuming it, verifies all six
coefficient identities in (5), checks every face split and normalized
substitution, computes both exact resultants and the first-branch gap,
and includes the affine-apex positive control and a rejecting `B`-axis
hostile control. No executable check uses Python `assert`.

Independent audit is pending. QED.
