# Critical and subcritical affine stopping grammars

Status: **PROVED**, self-contained for the probability model specified
below. Native integer orbit consequences are explicitly separated.

## 1. A family of stopping rules

Let `a_1,a_2,...` be independent valuations with `P(a=j)=2^{-j}`, `j>=1`.
This is exactly the finite-word law of normalized odd Haar measure for
Collatz, as inherited from
`session94_recursive_ports_crt_clocks_20261029/dynamics/VALUATION_TREE_CLOCKS.md`.
Fix integers `r>=2` and `1<=d<=2^r-1`. A letter has `d` ordered children
when `r|a`, and no children otherwise. The chronological first-extinction
word is defined by queue increments `d*1[r|a]-1`, starting at one.

Put

\[
 p={1\over2^r-1},\quad q={3\over4^r-1},\quad
 \mu=dp,\quad \nu=dq.
\]

The original tree has mean offspring `mu<=1`; its offspring law is
nondegenerate. Its extinction probability is one: the least fixed point
of `s=1-p+ps^d` in `[0,1]` is one (also immediate for `d=1`). The tilted
law below has mean `nu<=3/(2^r+1)<1`.

Let `m` count branching vertices. A completed tree has `J=dm+1` vertices
and `(d-1)m+1` leaves. Define nonnegative series

\[
 M(t)=E[t^m],\quad L(t)=E[t^m\Lambda],\quad
 B(t)=E[t^m B_w],\qquad \Lambda=3^J2^{-D},\ B_w=C_w2^{-D}.
\]

`B(t)` denotes a series of normalized carries, distinct from the
integer-carry moments `F_r` in RECURRENCES.md.

## 2. Three coupled recurrences, including the ordered carry

For `0<=t<1`, and also as formal series,

\[
 M=1-p+tpM^d,\qquad L=1-q+tqL^d,                      \tag{10}
\]

\[
 \boxed{B={L/3\over1-tp\sum_{i=0}^{d-1}L^iM^{d-1-i}}}. \tag{11}
\]

To prove the second recurrence in (10), note that the letter tilt is

\[
 2^{-a}(3/2^a)=3\,4^{-a},\quad \sum_{a>=1}3\,4^{-a}=1,
 \quad\sum_{r|a}3\,4^{-a}=q.
\]

For (11), apply the root map first and the `d` child maps in order.
Its root intercept contributes `beta_0` times the product of all child
slopes. Child `i` contributes its intercept times every later child slope;
earlier children still contribute their probability masses. Their sum is
`B*sum L^i M^(d-1-i)`. Also `beta_0=1/2^a=lambda_0/3`, so the root's total
contribution is `[(1-q)+tqL^d]/3=L/3`. This proves (11), retaining the
chronological asymmetry rather than replacing all children by one mean.

The denominator is positive for `t<1`: `M,L<=1`, so its subtracted term
is at most `t*mu<1`. Formal coefficients can be computed without radicals:
with `S=sum L^i M^(d-1-i)`,

\[
 M_0=1-p,\quad L_0=1-q,\quad
 M_m=p[t^{m-1}]M^d,\quad L_m=q[t^{m-1}]L^d,
\]
\[
 B_m=L_m/3+p\sum_{j=0}^{m-1}S_jB_{m-1-j}.             \tag{12}
\]

The first two series also have closed coefficients

\[
 M_m={1\over dm+1}{dm+1\choose m}(1-p)^{(d-1)m+1}p^m,
\]
\[
 L_m={1\over dm+1}{dm+1\choose m}(1-q)^{(d-1)m+1}q^m. \tag{13}
\]

For completeness the common counting factor follows by the cyclic lemma:
among the `dm+1` cyclically marked positions of a word with `m` increments
`d-1` and the remaining increments `-1` (sum `-1`), exactly one starting
position has nonnegative proper partial sums. There is no periodic-word
exception because a proper period would divide the total sum `-1`.
The resulting good words are the ordered tree traversals.

## 3. A mean-one slope does not pay the carry

**Theorem.** For every `(d,r)` above,

\[
 E\Lambda=1.
\]

If `mu<1`, then

\[
 EJ={1\over1-\mu},\quad ED={2\over1-\mu},\quad
 E(D+J)={3\over1-\mu},\quad EB_w={1\over3(1-\mu)}.     \tag{14}
\]

If `mu=1`, all four expectations in (14) are infinite. The slope-weighted
probability law `P*(w)=Lambda_w P(w)` has instead

\[
 E_*J={1\over1-\nu},\quad E_*D={4/3\over1-\nu}.        \tag{15}
\]

Proof. A finite word's tilted probability is the product of `3*4^{-a}`.
The tilted tree is subcritical, so these prefix-free finite words have
total tilted mass one. This proves `E Lambda=1`, or `L(1)=1`. The original
law similarly gives `M(1)=1`. Letting `t` increase to one in (11), using
monotone convergence for its nonnegative coefficients, proves the carry
claim. For the clocks, the root contributes mean `1` to `J` and `2` to
`D`; its mean number of independent children is `mu`. Summing the expected
generation contributions gives (14), with divergence at criticality.
Under the tilted letter law the root valuation mean is `4/3`, and the
mean number of children is `nu`, giving (15).

| Rule `(d,r)` | Original mean offspring | `EJ` | `ED` | `EB_w` |
|---|---:|---:|---:|---:|
| parity ternary `(3,2)` | `1` | infinite | infinite | infinite |
| valuation-multiple-of-three ternary `(3,3)` | `3/7` | `7/4` | `7/2` | `7/12` |
| critical seven-child `(7,3)` | `1` | infinite | infinite | infinite |
| critical fifteen-child `(15,4)` | `1` | infinite | infinite | infinite |

For the parity ternary rule, (15) gives `E_*J=5/2`, `E_*D=10/3`, and
`E_*(D+J)=35/6`. The `5/2` here is an explicitly tilted tree-size mean;
no identity with an earlier parity-amplitude process is asserted.

These are source **measure** identities. They do not assert independence
along a fixed integer orbit. For example the critical parity rule never
stops on the valuation sequence `(2,2,2,...)` at the already-converged
source `1`. Conversely the first-extinction word `(2,1,1,1)` carries
`41` to `107`. Thus tree extinction is neither necessary nor sufficient
for descent.

## 4. Exact critical singularity and the Mersenne arities

For `(d,r)=(3,2)`,

\[
 M=2/3+tM^3/3,\quad L=4/5+tL^3/5,\quad
 B={L\over3-t(L^2+LM+M^2)}.                            \tag{16}
\]

In particular `E[Lambda|m]=(6/5)(108/125)^m`, by (13). Conditional slopes
are less than one for `m>=2`, but averaging over all sizes gives exactly
one. Discarding the small-size terms would change the original source law.

Let `delta=1-t`. Expanding the first two equations at one gives

\[
 1-M\sim\sqrt{\delta/3},\qquad 1-L\sim\delta/2,
 \qquad \boxed{\lim_{t\uparrow1}\sqrt{1-t}\,B(t)=1/\sqrt3}. \tag{17}
\]

To justify the expansions, both roots increase to one; the linear term
in the equation for `M` cancels and its quadratic coefficient is nonzero.
For `L` the derivative coefficient `1-3/5=2/5` is nonzero. Substitution
in (16) then proves the limit. It is an Abelian limit for the generating
function; no coefficient asymptotic or fixed-source tail bound is claimed.

More generally take the **critical Mersenne arity** `d=2^r-1`, `r>=2`.
Here `p=1/d`, `q=3/[d(d+2)]`, and the same expansion proves

\[
 1-M\sim\sqrt{{2(1-t)\over d(d-1)}},\quad
 1-L\sim {3(1-t)\over d(d-1)},\quad
 \lim_{t\uparrow1}\sqrt{1-t}\,B(t)
       ={1\over3}\sqrt{{2d\over d-1}}.                \tag{18}
\]

This is an exact link between `3,7,15,31,...`: the dyadic probability of
a valuation divisible by `r` is exactly `1/(2^r-1)`. Primality of the arity
plays no role. It does not identify these trees with tournament towers.

## 5. An exact slope-moment boundary for the ternary parity rule

For real `s>=0`, set

\[
 x=2^{-(s+1)},\quad a_s={3^s x\over1-x^2},\quad
 b_s={3^s x^2\over1-x^2},\quad
 \Delta(s)={27\over4}a_s^2b_s.
\]

Then `T_s(t)=E[t^m Lambda^s]` is the minimal nonnegative solution of
`T_s=a_s+t b_s T_s^3`, possibly infinite at `t=1`. Therefore

\[
 E\Lambda^s<\infty\quad\Longleftrightarrow\quad\Delta(s)\le1. \tag{19}
\]

Indeed the cubic fixed-point equation has a positive solution exactly
when its maximum `y-b_s y^3` is at least `a_s`. Iteration from zero yields
the smaller solution; at the double-root boundary it is still finite.
If no root exists the increasing finite-tree sums diverge.

Explicitly

\[
 \Delta(s)={27\over4}{27^s2^{-4(s+1)}\over(1-4^{-(s+1)})^3},
 \quad\Delta(0)=1,\quad\Delta(1)=108/125,\quad
 \Delta(2)=432/343.
\]

The logarithmic derivative is
`log(27/16)-3 log(4)*4^{-(s+1)}/(1-4^{-(s+1)})`, strictly increasing,
negative at zero and tending to a positive number. Consequently there
is exactly one further root `s_*` of `Delta(s)=1`, with `1<s_*<2`.
Finite moments occur exactly on `[0,s_*]`. In particular `E Lambda=1`
but `E Lambda^2=infinity`.

## 6. What this changes in the research target

The critical rule has two distinct debts: a heavy slope tail already
visible in (19), and an infinite mean additive carry visible only after
retaining the extension coordinate. A subcritical rule removes the latter
infinite-mean obstruction, with explicit cost (14), while retaining mean
slope one. It therefore creates a better finite-cost input language, but
does not create a descending certificate by itself.

Three precise **OPEN** helper questions are:

1. For the `(3,3)` rule, can one bound a carry-sensitive deficit conditional
   on the original integer source, rather than under Haar averaging?
2. Can a finite port-state refinement distinguish the expanding first-stop
   words and send them to a paid second controller, with a subcritical
   weighted return operator and no source resampling?
3. Which finite carry jets, together with source residue and occupied-port
   sidecars, suffice for such an inequality under every legal attachment?

These ask for a source-preserving descent mechanism. Universal Collatz,
uniform ROOT bounds, six-colour forcing and all-order native tournament
wiring remain open.
