# Gauge-transverse half-step paths behind the THM-2825 collar

**Status: PROVED abstract algebra + FINITE-EXACT application, scratch
candidate.**  This note uses the hash-pinned primary and hostile-audit
outputs for proposed THM-2825 and performs a new all-cell path
decomposition.  Until THM-2825 is promoted, this remains outside the proved
dependency graph.  Even after promotion it gives no row exclusion and no
LRC(14) conclusion.

## 1. Inheritance and the missing type

The closest positive mechanism is the THM-2825 `+2h` collar: all `587`
labelled right-cofiber pieces have a unique common mate two half-steps to the
right with the same delayed semantic content.  The hostile side is just as
important:

- every right piece has empty source carrier support;
- every `+2h` mate has source and target carrier support at delta zero;
- native factor masks change;
- only `74/587` source endpoint profiles agree, although all `587/587`
  target endpoint profiles agree;
- local ancestry chambers agree in all `587` pairs.

Thus the collar should not be asked to be a carrier-preserving endomorphism.
The missing type is an **off-diagonal morphism from cofiber degree to common
degree**.  This is the least expensive way to make the absent source carrier
part of the statement rather than an unrecorded failure.

The incoming unnumbered half-step/Hasse note asked specifically for either a
collision-free even-half-step common collar or a decorated odd collar.  The
new result supplies the first object on the free labelled interval module,
and sharply records why it still does not descend to a physical carrier
action.

## 2. The exact path-quiver decomposition

Fix one of the `193` nonempty labelled cells `c=(e,s,t)`.  Write `R_c` for
its right pieces and `M_c` for its common pieces.  Every piece has the same
length and, within the cell, the same positive weight.  Let

```text
h = T/(2*13^5).
```

Put a directed edge `x -> x+h` whenever `x` is in `R_c union M_c` and
`x+h` is in `M_c`.  Across all cells this directed graph is a disjoint union
of finite paths:

```text
685 common path components,
587 preceded by one right-cofiber root,
98 preceded by no bank atom,
all 685 terminating outside the bank.
```

The `587` cofiber-rooted paths are collision-free and cover

```text
54,754 / 63,308
```

common atoms.  Their common-path length histogram is

```text
  2:60,   14:9,   17:16,  21:7,   28:67,  40:9,
 43:9,    54:37,  66:9,   69:16,  73:42,  80:37,
 92:9,    95:96, 118:9,  121:7,  144:36, 145:7,
158:14, 184:14, 210:21, 236:21, 262:21, 288:14.
```

The remaining `8,554` common atoms form `98` carrier-only paths:

```text
length 13:14 paths,   39:14,   65:14,
       91:14,        117:14, 143:28.                 (1)
```

In particular every carrier-only length is `13m` with `m` odd.  This
prime-power regularity was not part of the nearest-collar census.

## 3. Candidate lemma: a transverse/tangent factorization

Let `K` be a field of characteristic not two (in particular `Q` or `F_13`)
and let

```text
V = K<all labelled pieces R_c union M_c>.
```

This is the free coefficient module on the labelled interval
decomposition; it is not the original physical allocation module.  On its
distinguished equal-weight basis define

```text
N e_x = e_(x+h)  if x+h lies in M_c,
        0         otherwise,

P e_x = e_x       for x in M_c,
        0          for x in R_c.
```

Decompose

```text
d = P N (1-P),                  a = P N P.             (2)
```

Here `d` is the `587`-edge carrier-transverse collar and `a` is the
`62,623`-edge carrier-tangent common shift.  Exact path incidence gives

```text
N=d+a,       d^2=0,       d a=0,
rank(d)=587, rank(a)=62,623.                            (3)
```

The semantic grading is

```text
J e_x = +e_x  on live content,
        -e_x  on zero content.
```

The complete THM-2825 parity relation proves alternation on every
cofiber-rooted component, and the new companion directly checks all `8,554`
atoms of the carrier-only components.  Therefore

```text
Jd=-dJ,                   Ja=-aJ.                      (4)
```

The desired even collar is not a new primitive translation.  It factors as

```text
S := a d = P N^2 (1-P).                                (5)
```

Consequently

```text
rank(S)=587,             S^2=0,
JS=SJ,                   [P,S]=S.                      (6)
```

Equation `(6)` is the exact promised type: `S` preserves semantic content
while crossing from carrier-absent to carrier-present degree.  If
`Gamma=2P-1`, then `{Gamma,S}=0`.

On the normalized real basis, `S` is a positive partial isometry.  If
`Q_R` and `Q_2` project onto the right roots and their `+2h` images, then

```text
S* S=Q_R,              S S*=Q_2.
```

Hence

```text
F=S+S*
```

is a positive-cone-preserving permutation and a non-idempotent involution on
`Q_R V direct_sum Q_2 V`; it commutes with `J` and anticommutes with
`Gamma`.  This is a lawful operator on the free labelled coefficient
module, not merely a numerical support coincidence.

There are two multiplications here and they should not be conflated.  Under
operator composition, `S` is nilpotent and `F` is an involution.  The
entrywise support mask of `S` is still an idempotent Schur multiplier.
Because `N`, `a^(k-1)d`, `S`, and `F` are partial permutation matrices,
their Hilbert-space operator norms are one; their matching-support Schur
multipliers are contractive as well.  This gives an exact dimension-free
norm-one positive control for the signed-blocky/Schur lane, but only on this
selected labelled bank.  It does not establish the uniform wall
decomposition sought by HYP-9046.

### 3.1 Cross-agent synthesis: the absent root completes `M_2` to `M_3`

The independent Schur/scale audit proves that the two collar images

```text
M_1={r+h:r in R},                M_2={r+2h:r in R}
```

are disjoint and that their common-side linking algebra is
`M_2 tensor I_587`.  Adjoining the right root rather than discarding it
upgrades this canonically.  Let

```text
W_0 : H_R -> H_R    be the identity,
W_1 : H_R -> H_M1   be d,
W_2 : H_R -> H_M2   be S.
```

The three ranges are orthogonal and

```text
W_i* W_j = delta_(i,j) I_R.
```

Therefore

```text
E_ij=W_i W_j*,             0<=i,j<=2
```

satisfy the matrix-unit law

```text
E_ij E_kl = delta_(j,k) E_il.                          (6a)
```

They generate an exact

```text
M_3 tensor I_587
```

on the `1,761`-dimensional collar triple
`H_R direct_sum H_M1 direct_sum H_M2`.

If `B` is the compression of the half-step incidence operator to this
triple, then

```text
B=E_10+E_21,           B^2=E_20=S,          B^3=0.    (6b)
```

Thus semantic preservation at `+2h` is literally the square of the
three-state nilpotent ladder.  Relative to each right root, the semantic
grading is

```text
diag(1,-1,1),
```

while common-carrier degree is

```text
diag(-1,1,1).
```

This is the cleanest algebraic role for the absent source carrier: it is the
third summand of a linking algebra, not a failed common atom.  It also
strictly extends the common-side Pauli `M_2` found by the Schur audit.
Nothing in this `M_3` construction repairs the source endpoint or native
factor mismatch; its scope is still the labelled coefficient module.

The two gradings make the same fact a punctured `V_4` torsor.  Relative to a
root, the occupied bidegrees are

```text
R  : (semantic,carrier)=(0,0),
M1 : (1,1),
M2 : (0,1),
```

and the missing corner is `(1,0)`, an opposite-semantic cofiber state.  The
three arrows have the three nonzero `V_4` degrees

```text
deg(d)=(1,1),        deg(a)=(1,0),        deg(ad)=(0,1).
```

Their addition law is exactly the factorization `ad=S`.  This triangle
cannot be completed inside the present right bank: its semantic populations
are `573` live and `14` dead, so no semantic-flipping permutation of `R`
exists even before endpoint and factor decorations are imposed.  A genuine
fourth corner must therefore enlarge the cofiber object or use signed
multiplicity; it cannot be obtained by relabelling the current `587` roots.

The construction is gauge-covariant.  For any diagonal change of labelled
basis `G` commuting with `P` and `J`, replace `N` by `G N G^-1`; all of
`(2)--(6)` are conjugated identities.  Because every component is a path,
all nonzero edge phases can be gauged to one.  Thus the ranks and gradings
are intrinsic, but there is no loop holonomy or new gauge-invariant scalar
phase.

## 4. Exact nilpotence and the odometer boundary

The rooted path forest has sharp nilpotence

```text
N^289=0,                 a^288=0.                      (7)
```

More generally, `a^(k-1)d` is translation by `kh` on precisely the roots
whose paths have length at least `k`.  Its complete rank profile is encoded
by the histogram in Section 2; the compact plateaus are printed by the
companion.

At the odometer scale,

```text
rank(a^13 d)=527,        nullity on R = 60.            (8)
```

This exactly reproduces the independent direct `+14h` census:

```text
527 land in common,      60 land outside.
```

So the equality `14h=7(2h)` is genuinely realized by seven local even
steps on the `527` long paths.  It fails on the `60` length-two paths for a
typed boundary reason, not because direct and iterated translations
disagree.

## 5. The rootless remainder is a signed prime-power unit

Let one carrier-only path in `(1)` have length `13m`, where
`m in {1,3,5,7,9,11}`.  Root it at its left endpoint and reduce its
half-step word in

```text
F_13[C_13] = F_13[X]/(X^13-1).
```

All `98` paths start live and alternate.  Put

```text
N_13 = sum_(r=0)^12 X^r,
A_13 = sum_(r=0)^12 (-1)^r X^r.
```

If `R` is the raw mask and `L` the live mask, exact residue counting gives

```text
R = m N_13,
2L-R = A_13,
(1+X)A_13=2.                                          (9)
```

Thus the raw rootless block is a nonunit norm element, but both its
semantic contrast and its positive live mask are units:

```text
(2L-R)^(-1) = (1+X)/2,
L^(-1)      = (1+X)(1-m N_13).                        (10)
```

The companion verifies `(9)--(10)` for all `98` components.  This gives a
uniform two-tap decoder for the signed contrast and a closed inverse for
the live mask.  It also identifies the mechanism: semantic alternation
turns a rootless augmentation-zero `C_13` norm tower into a full-spectrum
unit.

The inverse in `(10)` is signed, as required by the positive-convolution
boundary.  No claim is made that it is a positive physical transport.

## 6. Physical sidecars and the exact stopping line

The hash-pinned physical and independent outputs give the following complete
typing of `S`.

```text
right native-factor holes:
  E3       319,       E3+c2     37,
  c2       217,       q1        14;

+2h image:
  all six native factors present in all 587;

carrier:
  right  = source empty, target delta0,
  image  = source delta0, target delta0;

endpoint profile:
  target agrees with delta0 translation in 587/587,
  source agrees in only 74/587;

ancestry chamber:
  agrees in 587/587.
```

Therefore `S` is a target-endpoint-preserving, ancestry-preserving,
semantic-preserving labelled coefficient morphism which deliberately
crosses the source-carrier and native-factor boundary.  It is **not** a
lawful endomorphism of the original physical carrier.  Calling it a
physical transport would erase exactly the coordinate it was designed to
expose.

The interval translation does descend after forgetting repeated cell labels
at the level of source/image interval multiplicities:

```text
587 labelled roots collapse to 37 physical interval/weight triples,
and their 587 images collapse to the translated 37 triples.

root multiplicities:
  7:24, 9:7, 42:1, 49:1, 56:2, 72:1, 81:1.
```

This still does not descend the semantic `q`-pair, owner, factor, or
source-endpoint decorations, so the labelled direct sum remains
load-bearing.

The target-label convolution axis is also independent of this half-step
operator.  The earlier `t=12` equivariance hostile remains valid: `S` does
not turn physical half-step translation into target convolution.

## 7. Cheapest next test

The new operator pays the former “absent source carrier” objection by
changing its type.  The next decisive test is now narrower:

1. enrich the `587` arrows of `S` by the **source endpoint defect class**
   (distances/counts `0:74, 9:187, 10:245, 81:81`);
2. ask whether the defect is a coboundary along each of the `587` rooted
   paths, with zero boundary on the target endpoint;
3. if so, conjugate `N` by that path potential and test whether one of the
   `527` length-at-least-fourteen paths gives a decorated `a^13d` cospan at
   the odometer scale;
4. use the `60` length-two paths as the hostile boundary.

Path acyclicity makes this the cheapest possible cohomological test: a
one-dimensional defect always has a formal path primitive, so the only
substantive question is whether that primitive is a lawful endpoint/address
gauge.  A negative answer should name the first factor or endpoint coordinate
that prevents it.

## 8. Reproduction

```text
python3 04-computation/lrc14_nearest_half_step_common_right_collar_path_operator_thm2825.py
python3 -O 04-computation/lrc14_nearest_half_step_common_right_collar_path_operator_thm2825.py
```

The final normal, optimized, and stored scratch transcripts agree byte for
byte.  The script contains no Python `assert` node.

```text
script_sha256 =
  7f0780e70161cbfafc4499d02a7d1c5aa40366b6dfa9935b00dc652e3d54c8e0

output_sha256 =
  6cff846944c59d8e243e74afe8829046feb288f0cf23da1c998037dcf70411b2
```
