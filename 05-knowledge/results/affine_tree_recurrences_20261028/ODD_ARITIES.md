# Incoming odd-arity clocks, the base-four tilt, and a variance separator

Status: **PROVED** synthesis of the incoming periodic-clock compiler and
the affine recurrences in this packet. No identification of probability
laws on different sources is made.

## 1. Inherit the stronger clock classification

The concurrent
[Session95 residue-clock theorem](../session95_niche_recurrences_20261029/valuation/RESIDUE_CLOCK_BRANCHING.md)
already realizes **every odd arity** `d>=3` critically, not just Mersenne
arities. Let `h=ord_d(2)`, and choose the unique subset `S` of `{1,...,h}`
such that

\[
 \sum_{j\in S}2^{h-j}=(2^h-1)/d.
\]

A node has `d` children when `a modulo h` lies in `S` (residues use
`1,...,h`). Thus `p=P(branch)=1/d`. The new step is to propagate its
**ordered affine carry**, rather than derive the already-known size law
again. As before, parse the same chronological word depth-first; changing
the drawing does not reorder the arithmetic execution.

## 2. An exact binary-to-base-four transformation

Let `epsilon_a` be the periodic indicator of this branching event. Then

\[
 {1\over d}=\sum_{a\ge1}\epsilon_a2^{-a},\qquad
 q=\sum_{a\ge1}3\epsilon_a4^{-a}
   ={3\sum_{j\in S}4^{h-j}\over4^h-1}.               \tag{26}
\]

So the binary digits of `1/d`, read as base-four digits `0` and `3`, are
exactly the branching probability under the slope-tilted law. This is an
explicit map to the base-four Cantor digit set, not the numerical square
of `1/d`: squaring would introduce cross terms between digit positions.

It preserves the entire selector and its least period. Both digit
expansions are unambiguous for these periodic, non-dyadic inputs. Forgetting
the selector and retaining only a tree size still loses the actual
valuation residues and source carry.

There are useful rigorous bounds:

\[
 1/d^2<q<3/d^2,\qquad dq<3/4.                         \tag{27}
\]

For the first pair, write `x_a=epsilon_a 2^{-a}`. The upper bound follows
from `sum x_a^2 < (sum x_a)^2`. For the lower bound use
`sum_{b>a}x_b<=2^{-a}` to obtain
`(sum x_a)^2<=3 sum x_a^2=q`; equality would require a full eventual
tail of ones, giving a dyadic input. The upper equality would require
only one nonzero digit. Both are excluded here. Finally `epsilon_1=0`
since `1/d<1/2`, and `3/2^a<=3/4` on every selected letter, strictly
smaller on infinitely many selected letters. Hence `q<(3/4)p`.

The slope-tilted tree is therefore subcritical for every such odd arity.
Equations (10)–(12) remain valid with `p=1/d` and the `q` in (26):

\[
 M=1-1/d+(t/d)M^d,\quad L=1-q+tqL^d,
 \quad B={L/3\over1-(t/d)\sum_{i=0}^{d-1}L^iM^{d-1-i}}. \tag{28}
\]

The proof uses only the probabilities of the two letter types and
`beta_a=lambda_a/3`; it does not require a singleton selector. Thus
`E Lambda=1`, `EB=infinity`, and

\[
 E_*J={1\over1-dq},\quad E_*D={4/3\over1-dq},\qquad
 \lim_{t\uparrow1}\sqrt{1-t}\,B(t)
        ={1\over3}\sqrt{{2d\over d-1}}.              \tag{29}
\]

For the last identity `1-M~sqrt(2(1-t)/(d(d-1)))` while
`1-L=O(1-t)` by `dq<1`. The leading carry singularity is independent of
the exact tilted probability `q`, though the finite coefficients are not.

## 3. Two arithmetic families and a five-child control

For `d=2^r-1`, the selector is `{r}` and `q=3/[d(d+2)]`, recovering
MOMENTS.md. For `d=2^r+1`, the least period is `2r`, the selector is
`{r+1,...,2r}`, and

\[
                  q={1\over4^r+1}.                  \tag{30}
\]

Indeed `(2^{2r}-1)/(2^r+1)=2^r-1` is the final block of `r` binary ones;
the base-four evaluation is a geometric sum. The period is `2r`: a
proper divisor of `2r` is at most `r`, and `2^h-1<2^r+1` for `h<=r`.
No prime-value assumption occurs in either family.

| `d` | `h,S` | `q` | Tilted `E_*J` |
|---|---|---:|---:|
| 3 | `2,{2}` | `1/5` | `5/2` |
| 5 | `4,{3,4}` | `1/17` | `17/12` |
| 7 | `3,{3}` | `1/21` | `3/2` |
| 9 | `6,{4,5,6}` | `1/65` | `65/56` |
| 15 | `4,{4}` | `1/85` | `17/14` |
| 17 | `8,{5,6,7,8}` | `1/257` | `257/240` |

**Variance separator.** Critical ternary stopping has `E Lambda^2=infinity`.
For the critical five-child selector, the squared-slope letter masses are

\[
 b=9\sum_{a\bmod4\in\{3,4\}}8^{-a}=9/455,\qquad
 a=9/7-b=576/455.
\]

Consequently its squared-slope moment is the minimal nonnegative root of
`Y=a+bY^5`. At `Y=2` the right side is `864/455<2`; iteration from zero
is bounded. Therefore

\[
       \boxed{E\Lambda^2\le2\quad\hbox{but}\quad EB=\infty}. \tag{31}
\]

The carry obstruction is not just another name for infinite slope
variance. Changing the clock selector can remove the latter without
removing the former. These are auxiliary first-extinction rules over
the same source law, not alternative Collatz dynamics.

## 4. Recurrence memory and the incoming port/flow work

The incoming [pole-memory theorem](../session95_niche_recurrences_20261029/jets/POLE_RECURRENCE_MEMORY.md)
proves that a known support bound `B` makes `2B` moments sufficient for
an unknown finite weighted source, with a sharp uniform boundary. Apply
that general Vandermonde argument to the finite carry distribution at
fixed `(J,D)`: a finite jet becomes faithful once an appropriate support
bound is retained. The bound need not stay small as `D` increases. This
strengthens RECURRENCES.md's known-support formulation; it does not
justify a depth-independent fixed jet bank.

The incoming [Schur boundary recurrence](../session95_niche_recurrences_20261029/flows/SCHUR_BOUNDARY_RECURRENCES.md)
similarly retains joint input/output partitions and closes under its
specified complete interchild attachments. Its scalar circulation
series alone has a degree-12 hostile. Our port-indexed affine product
retains joint endpoints for precisely the same compositional reason,
but neither result identifies flow degree, valuation clock, moment order,
or the physical port menu with one another.
