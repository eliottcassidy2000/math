# Universal integer entry and five exact expanding return types

**PROVED:** the cheaper ternary queue terminates for every positive odd integer,
with explicit bounds on both native clocks. **PROVED with a finite exact
certificate:** its return increases the source precisely on five stated residue
classes. **OPEN:** convergence of repeated returns. This is not a Collatz proof.

We inherit the `(d,r)=(3,3)` grammar and affine/source conventions from
[MOMENTS](../affine_tree_recurrences_20261028/MOMENTS.md) and
[REFINEMENT](../affine_tree_recurrences_20261028/REFINEMENT.md). No independent
source is sampled at a failed return. All numbers below are positive odd integers
unless a different domain is explicitly stated.

## 1. An integer certificate, independent of Haar measure

Set `n_0=n`, `Q_0=1`, and, while `Q>0`, perform the actual odd Collatz step

\[
 a=v_2(3n+1),\qquad n'=(3n+1)/2^a,\qquad
 Q'=Q-1+3\mathbf1_{3\mid a}.
\]

Let `J` be the first time `Q=0`, `D` the sum of the actual valuations, and
`R(n)=n_J`. These quantities will now be proved to exist, not assumed.

For `n>=7`, take `W(n,Q)=n^3 4^Q`. If `3` does not divide `a`, then `a>=1`
and `Q'-Q=-1`; otherwise `a>=3` and `Q'-Q=2`. In either case

\[
 \frac{W(n',Q')}{W(n,Q)}
 \le \frac{(3+1/n)^3}{32}
 \le \rho:=\frac{1331}{1372}<1.                 \tag{1}
\]

An odd Collatz image is never divisible by three. Thus a trajectory starting
at `n>=7` subsequently visits only `1,5`, or values at least seven. At `1`
the ratio in (1) is `1/4`; at `5` it is `1/500`. Hence (1) holds all the way
to the actual endpoint, including the last step. If it never terminated, its
positive integer potential would decrease forever, a contradiction. The
sources `1,3,5` each terminate in one step, with endpoints `1,5,1`.

There is also a certificate for **every initial queue size**: replace `n^3`
by `f(n)`, where `f(3)=32` and `f(n)=n^3` otherwise. At `3 -> 5` the ratio
is `125/128`; no step has target three. Thus `f(n)4^Q` contracts by at most
`125/128` for every positive odd state. Absorbing terminal states can be assigned
potential zero if a conventional ranking function is wanted.

For `K=max(1,ceil(log_2 n))`, the exact integer inequality
`2*1331^23 < 1372^23` gives

\[
 \boxed{J\le69K,\qquad D\le139K,\qquad J+D\le208K.}       \tag{2}
\]

Indeed, a nonterminal state has `W>=4`, whereas after `69K` steps (1) gives
`W<4n^3/2^(3K)<=4`. Also `3n_i+1<=4n_i`, so telescoping gives
`2^D n_J<=4^J n_0` and `D<=2J+floor(log_2 n_0)`. These are bounds on the
odd-step, division, and usual unaccelerated integer clocks respectively.
For even initial integers, first remove their finite power of two; that
initial division cost is additional and explicitly known.

This eliminates exceptional **positive integer non-entry sources for this
auxiliary grammar**. It does not eliminate exceptions to Collatz convergence.
In particular, `3 -> 5` remains a terminating but expanding return.

## 2. Every non-descending block has an extremely small alphabet

Suppose `n>=7` and `R(n)>=n`. Its trajectory cannot visit `1` or `5`, since
from those values it stays at or reaches one. Thus every intermediate source
is at least seven. Write `m` for the number of branching valuations. Queue
balance forces `J=3m+1`, and

\[
 D\ge5m+1,
 \qquad \frac{R(n)}n\le\frac{(22/7)^{3m+1}}{2^D}.          \tag{3}
\]

If `D>=5m+2`, the right side is at most `(11/14)rho^m<1`. Therefore equality
holds in the minimum valuation budget: **all leaves have valuation one and
all branches valuation three**. Moreover (3) is then at most
`(11/7)rho^m`; the exact inequality `11*1331^15<7*1372^15` forces `m<=14`.
Equivalently, any block with at least 46 odd steps strictly descends.

## 3. Closing the finite exceptional window

For the `{1,3}` tree words with `m` branches, `J=3m+1`, `D=5m+1` are fixed.
Write their affine action as `(3^J n+C)/2^D`. The maximum possible numerator
`C_m` is computed by the following finite recurrence. Start `C_0=1`. For every
ordered triple `i+j+k=m-1`, start `(c,d)=(1,3)` for the root and append its
three children in order, updating

\[
 (c,d)\longmapsto
 \bigl(3^{3u+1}c+2^d C_u,\ d+5u+1\bigr),\quad u=i,j,k.
                                                               \tag{4}
\]

Take the largest resulting `c`. This is exact by induction: child lengths and
valuation sums are fixed, and every child numerator has a positive coefficient.
For `3<=m<=14`, the slope is below one, so non-descent requires

\[
 n\le\frac{C_m}{2^{5m+1}-3^{3m+1}}.
\]

The largest of these twelve exact rationals is

\[
 \frac{362205083235814619852970625}{2032926274040285529221}
 \quad\text{with floor }178169.                              \tag{5}
\]

The accompanying certificate checks **every odd source from 1 to 178169**
using two implementations (bit valuation versus literal repeated division).
None has a non-descending block with `m>=3`. The recurrence supplies the
unbounded reduction; the exhaustive finite check closes the remaining window.
Independent enumeration of all tree words through `m=8` also checks (4).

For `m=0,1,2`, there are exactly five words. All five have slope greater than
one and positive carry, so every positive source realizing them strictly grows:

| Valuation word | `J` | `D` | `C` | Exact original-source guard |
|---|---:|---:|---:|---|
| `1` | 1 | 1 | 1 | `n=3 mod 4` |
| `3111` | 4 | 6 | 179 | `n=125 mod 128` |
| `3311111` | 7 | 11 | 16177 | `n=1357 mod 4096` |
| `3131111` | 7 | 11 | 12289 | `n=2269 mod 4096` |
| `3113111` | 7 | 11 | 9697 | `n=2877 mod 4096` |

The guard for any exact valuation word is
`n=(2^D-C)*3^(-J) mod 2^(D+1)`. It includes the last oddness condition;
reducing the modulus to `2^D` would merge the terminal valuation alternatives.
Directly check `n=1,3,5` to finish all positive odd sources. Consequently

\[
 \boxed{R(n)>n\iff n\text{ lies in the five classes above};
 \quad R(n)=n\iff n=1.}                                    \tag{6}
\]

This is a computer-assisted universal classification, not a conjecture based
on the larger test through one million. No source or endpoint is quotiented
away when selecting a word.

## 4. Why resetting the queue is still an unpaid step

`R` is a total acceleration of positive odd Collatz: each application consumes
between one and finitely many native odd steps, with the exact division clock.
The original orbit reaches one iff the iterated return map does. If a block
reaches one internally, subsequent valuations are two and drain its queue.

The certificate (1) does not telescope across restarted blocks for free:
resetting `Q=0` to `Q=1` multiplies `W` by four. A typical failure `3 -> 5`
makes the issue explicit. The remaining target is a source-sensitive payment
for these repeated resets and for intervening descending blocks. Section 5 of
[RENEWAL](RENEWAL.md) states the surviving pointwise obligation.

No claim of literature priority or independent-agent audit is made. The proof
has two exact finite implementations and symbolic inequalities; formalization
or an external audit would be a further check.
