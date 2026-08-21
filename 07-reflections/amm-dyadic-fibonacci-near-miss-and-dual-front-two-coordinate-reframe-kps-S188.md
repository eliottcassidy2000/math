# AMM dyadic Fibonacci near-miss and the dual/front two-coordinate reframe

Status: **EXACT OBSERVATION + REFUTED LAW + OPEN REFRAME** (2026-08-21).

## 1. Inheritance pass

- THM-3330 and its archived all-`R` transient companion give the exact Rule-A
  conjugacy, front/stall dynamics, and adjacent death/closure thresholds.
- THM-3331 gives the superblock continuation obstruction in its stated rule
  class.
- THM-3588 gives the first positive truncated-adjoint Farkas wall at `R=512`.
- THM-3597 extends that dual wall exactly through `R=1024` and replaces the
  repeated-cut computation by one backward sweep.

These are three different observables: policy survival threshold, front
capture/death, and earliest forced departure from a doomed prefix.  None is a
synonym for feasibility of the full AMM entry polytope.

## 2. A seductive Fibonacci law, killed by the archived hostile scale

The exact first-surviving Rule-A offsets at dyadic epochs begin

```text
R:       256  512  1024  2048  4096  8192
D0*:       1    5    15    38    89   192.
```

Writing `R=256*2^m`, the first five values fit perfectly

```text
T_m=2T_(m-1)+F_(m+3)=9*2^m-F_(m+6),                  (1)
```

which gives `1,5,15,38,89`.  The fit is not merely numerical decoration:
the degree word itself is defined by a Fibonacci--Lucas comparison, so a
Fibonacci correction has a plausible source.  Nevertheless `(1)` predicts

```text
T_5=9*32-F_11=199,
```

while the exact adjacent pair in the older transient atlas is

```text
R=8192: D0=191 dies at row 2045; D0=192 closes.       (2)
```

Thus `(1)` is **REFUTED**, with minimal displayed witness `192!=199`.  The
same archive places `D0*(16384)` in `[401,416]`, whereas `(1)` predicts 432,
so the discrepancy grows rather than oscillating by a floor error.

The lesson is specific.  Fibonacci--Lucas arithmetic determines the moving
degree word, but the survival threshold also depends on the nonlinear
high-cell clamp and its evolving front.  A recurrence fitted only to the
deadline word has forgotten the state coordinate that first becomes visible
at `R=8192`.

## 3. The two-coordinate proof diagram

For a failing Rule-A trace, retain both

```text
j(R,D) = first fatal top row,
q(R,D) = last row by which every surviving alternative must differ. (3)
```

The first is a primal policy observable.  The second is a dual
proof-obligation horizon.  At the last failing offsets currently covered by
the adjoint atlas,

| `R,D` | `j` | `q` | `q/R` |
|---:|---:|---:|---:|
| `256,0` | 61 | 24 | 0.093750 |
| `512,4` | 121 | 44 | 0.085938 |
| `1024,14` | 250 | 94 | 0.091797 |
| `2048,37` | 508 | 195 | 0.095215 |

The normalized dual horizon is roughly stable near nine percent across these
four scales even though the greedy offset threshold is moving toward a
positive linear slack.  The fourth point sharpens a possible factorization:

```text
j/R = 0.248046875 ~ 1/4,
q/j = 0.3838582677... ~ phi^-2,
q/R = 0.09521484375 ~ phi^-2/4 = 0.0954915028....       (4)
```

The appearance of `phi^-2` is not numerology-free evidence: THM-3027 proves
that the same number is the universal binding mass fraction in the
archimedean capacity problem.  But the four-point sequence is nonmonotone
(`R=512` is a visible hostile point), so this remains a pre-registered target,
not a limit claim.  The decisive next datum is the terminal failed offset at
`R=4096`.  Conceptually the right object is a limiting dual profile
`B(theta;epsilon)` whose first sign change is `theta=q/R`.

The superblock theorem supplies a third coordinate: the row at which a
same-sign front pair enters an absorbing death regime.  A lawful AMM no-go
would require an overlap statement of the form

```text
every departure before q enters a superblock basin before escape.          (5)
```

No current theorem proves `(4)`.  The adjoint wall constrains *when* a repair
must occur; the front theorem constrains *what happens after a specified
state enters its basin*.  Their state spaces and allowed continuation classes
must be matched explicitly before composition.

## 4. The adjoint is an exact top-distance walk

The original certificate routine recomputed the same transposed Pascal cone
for every cut.  If `lambda_i` is the level propagated from the fatal top cell,
then all cuts arise in one reverse pass:

```text
B_i=B_(i+1)+<lambda_i,u_i^A-L_i>,
lambda_(i-1)=truncate(K_i^T lambda_i).                 (6)
```

This reduces the cut atlas from repeated backward cones to one cone per
failed offset.  THM-3597 pins `(6)` against the legacy coefficient-cancellation
routine at both boundary cuts of the hardest offset in each epoch, and
THM-3601 does so again at `R=2048`.

There is a further exact compression.  Put

```text
delta_i=d_i-d_(i-1) in {0,1},
r=d_i-1-t,
M_i(z)=sum_r lambda_i(d_i-1-r) z^r.                    (7)
```

The fatal top cell initializes `M_(j-1)=1`.  Directly changing coordinates in
the Pascal transpose gives

```text
delta_i=0: M_(i-1)(z)=(1+z)M_i(z),
delta_i=1: M_(i-1)(z)=Pi_(>=0)[z^(-1)(1+z)^2 M_i(z)].  (8)
```

Thus the dual multiplier is a one-dimensional, time-inhomogeneous boundary
walk in *distance below the current top*.  A flat degree step permits `r` to
stay or increase by one; a rising degree step permits moves `-1,0,+1` with
weights `1,2,1`, while deleting the illegal move below `r=0`.  This proves
that the dual cone and the primal high-cell front already live in the same
coordinate.  The missing overlap theorem is therefore not a mysterious map
between two state spaces: it is a sign/escape statement for the charges

```text
h_i(r)=u_i^A(d_i-1-r)-L_i(d_i-1-r)                    (9)
```

sampled by the walk `(8)`.

This is the most concrete new analytic target.  Derive a scaling limit for
the Sturmian `delta_i` word and the charged walk `(8)--(9)`.  If its zero is
`theta=phi^-2/4`, the coincidence in `(4)` has a mechanism; if not, the
`R=4096` exact wall should expose the missing state variable.

The completed exact target is

```text
R=2048, D0=37: j=508, first negative cut=196, q=195.    (10)
```

## 5. Cross-frontier pattern

This episode repeats the exact failure mode caught in THM-3592: quotienting
or compressing by a symmetry of the easy representation can destroy a
downstream side condition.  There it was simultaneous support reversal versus
signed regularity floors.  Here it is Fibonacci degree-word arithmetic versus
stateful high-cell clamp dynamics.  The reusable move is to carry one extra
sidecar coordinate through the compression and demand the first hostile scale
before naming a law.
