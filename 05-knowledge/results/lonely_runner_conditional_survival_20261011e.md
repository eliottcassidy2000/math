# Conditional survival for lonely runners: products, sheets and closed walls

**PROVED:** elementary source-bound conditional identities, finite CRT gluing,
the complete source-clock criterion and closed-wall recovery below.
**FINITE-EXACT:** the declared executable controls.
**OPEN:** arbitrary LRC(14). No independence heuristic, positive average or
failed clock grid is promoted to a universal loneliness statement.

## 1. Inheritance and the coordinate that products omit

The anchor is a common-time survival certificate. The niche is exact
conditional sheet counting; the wildcard is recovery of an isolated witness
from nearby positive-measure sets. The concept board is **radius / common
time / conditioned residue profile / positioned interval / closed wall**.

- [LEM-042, pair-overlap law](../../01-canon/theorems/LEM-042-pair-overlap-law.md)
  supplies the correct all-wrap pair formula. The two-branch formula in
  [THM-594](../../01-canon/theorems/THM-594-pair-overlap-law-mirsky-newman-floor.md)
  Part B is explicitly superseded: at radius $1/14$, $(50,51)$ has overlap
  $73/3570$, not $137/17850$ (MISTAKE-327).
- [THM-2064, multi-tail sheet capacity](../../01-canon/theorems/THM-2064-multitail-sheet-capacity-and-dyadic-seam.md)
  already keeps all tails on one common fibre, uses strict danger and proves
  a sufficient capacity criterion. Its dyadic seam remains a residual.
- [THM-1162, corrected whole-safe-set conditioning](../../01-canon/theorems/THM-1162-whole-safe-set-conditioning.md)
  retains positioned components, not a single invented bad run.
  [THM-3910](../../01-canon/theorems/THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response.md)
  and [THM-4112](../../01-canon/theorems/THM-4112-antipodal-component-ancestry-chain-and-scale-separated-lrc-families.md)
  already use component counts and ancestry for scoped scale-separated rules.
- MISTAKE-235 rejects the identification of the zero-frequency term
  $(6/7)^{13}$ with an LRC existence theorem. A tight AP has zero strict
  survival measure and still has weak witnesses.

The incoming [context/prime/ROOT note](context_prime_root_20261011.md)
suggests a useful discipline: choose a context that retains the marked
source and its native map. Here the retained source is the speed tuple and
the actual clock; a residue label without that clock does not identify time.
The weighted gluing below is an arithmetic counterpart of the earlier
[boundary-context count](boundary_context_weights_20261011d.md).
These are elementary reformulations and a reusable verifier, not claims
of novelty for CRT, conditioning, or the tight AP.

## 2. What $5/6$ and $6/7$ measure

Fix distinct positive integer relative speeds $V=(v_1,\ldots,v_r)$ and
$0<\delta<1/2$. On one time period with normalized Lebesgue measure, set

\[
 A_v(\delta)=\{t:\|vt\|\ge\delta\},\qquad
 S_i=\bigcap_{j=1}^i A_{v_j}(\delta),\qquad \mu_i=|S_i|.
\]

Every marginal has measure $s=1-2\delta$. The same formula holds for strict
safety, because one runner has only finitely many boundary times. Existence
is a different predicate and does depend on those boundaries.

| Total runners, including the stationary reference | Relative speeds | Usual radius | One-runner safe fraction |
|---|---:|---:|---:|
| 6 | 5 | $1/6$ | $2/3$ |
| 7 | 6 | $1/7$ | $5/7$ |
| 12 | 11 | $1/12$ | $5/6$ |
| 14 | 13 | $1/14$ | $6/7$ |

Independent uniform positions $\Theta_1,\ldots,\Theta_r$ on an $r$-torus
give survival exactly $s^r$. This is a different sample space from
$(v_1t,\ldots,v_rt)$ for one uniform time $t$.
For the latter, the exact chain rule is

\[
 \mu_r=\prod_{i=1}^r \rho_i,\qquad
 \rho_i=\frac{|S_{i-1}\cap A_{v_i}|}{|S_{i-1}|},\quad \mu_0=1, \tag{1}
\]

as long as preceding denominators are positive. If a factor becomes zero,
all later masses are zero; later conditional probabilities under this
measure are undefined. The closed survivor set can nevertheless be nonempty.
Ordering changes the factors, not the final event.

Two signs of dependence already occur at $\delta=1/14$:

\[
 |A_1\cap A_2|=11/14>(6/7)^2,
\]
\[
 |A_{13}\cap A_{15}|=(6/7)^2-6/(49\cdot13\cdot15)<(6/7)^2.
\]

For arbitrary radius, direct arc counting gives the finite all-wrap law
for $a\le b$, $g=\gcd(a,b)$:

\[
 |D_a\cap D_b|=
 g\sum_{j\in\mathbb Z}\min\left(
 \left(\frac{\delta(a+b)-|j|g}{ab}\right)_+,\frac{2\delta}{b}\right).
 \tag{2}
\]

There are finitely many nonzero summands. Each summand is the overlap of
two intervals centred at the corresponding difference of comb centres;
every difference has multiplicity $g$. At $\delta=1/14$ this is LEM-042.
Since $|A_a\cap A_b|=1-4\delta+|D_a\cap D_b|$, its credits have both signs
relative to the independent value. Full-period credits are not automatically
conditional credits on an inherited window.

## 3. Exact positioned conditional updates

Let $S$ be the complete closed survivor set of a fixed finite source. Write
its positive-length components as $[a_j,b_j]$, $1\le j\le C$; singleton
components are retained separately and contribute zero mass. Let $F$ be the
continuous periodic primitive of $1_{[\delta,1-\delta]}-s$, normalized by
$F(0)=0$. For $z\in[0,1]$ it is

\[
 F(z)=\begin{cases}
 -sz,&0\le z\le\delta,\\
 -s\delta+2\delta(z-\delta),&\delta\le z\le1-\delta,\\
 s(1-z),&1-\delta\le z\le1.
 \end{cases}
\]

Its range is $[-\delta s,\delta s]$. Integration on every component gives
the exact boundary bill

\[
 |S\cap A_v|=s|S|+
 \frac1v\sum_j\{F(vb_j)-F(va_j)\}.                       \tag{3}
\]

In particular

\[
 |S\cap A_v|\ge s|S|-\frac{2\delta sC}{v}.                \tag{4}
\]

If $|S|>0$ and $v|S|>2\delta C$, the appended runner leaves positive
measure, hence a witness. Iteration is valid when this inequality is checked
against each **actual preceding** survivor. This is a sufficient adaptive
criterion, not a universal speed-order theorem. The exact signed bill (3)
is stronger than its coarse bound and keeps endpoint positions that mass
and component count alone omit. The earlier component and ancestry results
are not replaced by this simple estimate.

## 4. Conditional CRT products on one native time grid

Let $A\subseteq\mathbb Z/Q$ be the already surviving residues and
$B\subseteq\mathbb Z/q$ the next allowed residues. Put
$d=\gcd(Q,q)$ and $L=\operatorname{lcm}(Q,q)$. The labelled profiles are

\[
 a_c=|\{a\in A:a=c\pmod d\}|,\qquad
 b_c=|\{b\in B:b=c\pmod d\}|.
\]

The generalized CRT gives exactly one residue modulo $L$ for each compatible
pair. Therefore the complete lifted set has cardinality

\[
          Z=\sum_{c\bmod d}a_cb_c.                        \tag{5}
\]

Every old residue has $q/d$ lifts before the new test. For $A\ne\varnothing$
the exact conditional fraction is

\[
          \rho=\frac{\sum_c a_cb_c}{|A|q/d}.               \tag{6}
\]

Survival is equivalent to some shared residue having $a_cb_c>0$.
A lower bound $b_c\ge\eta q/d$ on every occupied old class gives
$\rho\ge\eta$. When $d=1$, equation (6) is $|B|/q$: this is an actual
independence theorem on a common finite clock. Without the profile, equal
marginal counts can have opposite outcomes: take $Q=6,q=2,B=\{0\}$;
$A=\{0,2,4\}$ gives $\rho=1$, while $A=\{1,3,5\}$ gives $\rho=0$.

For the original speeds sampled at $t=k/N$, speed $v$ has native period
$q=N/\gcd(N,v)$ and multiplier $u=v/\gcd(N,v)$, with

\[
 B=\{k\bmod q:\|uk/q\|\ge\delta\}.
\]

The implementation glues these actual sets. The final period can be less
than $N$; a reported residue still means time $k/N$, **not** $k/q$ or the
current accumulated-period quotient. Repetition to modulus $N$ is uniform.
For speed two on clock four, the safe residue one has period two but means
time $1/4$; substituting time $1/2$ makes that runner unsafe.

As a constructive native family, choose pairwise coprime $q_i\ge2$,
$N=\prod q_i$ and $v_i=N/q_i$. Then the common-time survival probability is

\[
 \prod_i\frac{q_i-2\lceil q_i\delta\rceil+1}{q_i}.        \tag{7}
\]

Each factor must be positive for a witness; this is automatic for
$\delta\le1/3$. Thus this supplies an infinite cofactor family, by ordinary
CRT, rather than a random-phase surrogate. The finite factors differ from
$1-2\delta$ by rounding and closed-boundary effects. At radius $1/14$,
$(q_i)=(2,17,19)$ has survival $112/323$ and an actual witness $3/646$.
With shared factors, even the native product fails: periods $(6,10)$ give
$23/30$ survival at that radius, whereas $(5/6)(9/10)=3/4$.

## 5. Closed boundary, complete clock, and a quantitative limiting bridge

For the tight AP $V=\{1,\ldots,k-1\}$, at $\delta=1/k$ the complete weak
safe set is

\[
             \{a/k:1\le a<k,\ \gcd(a,k)=1\}.              \tag{8}
\]

Indeed the $k$ points $0,t,\ldots,(k-1)t$ are pairwise at distance at least
$1/k$. Their circular consecutive gaps are therefore all at least $1/k$;
their sum is one, so every gap equals $1/k$. The points form the regular
$k$-grid, and $t$ must generate it. Conversely every unit numerator works.
This proves (8) for $k\ge2$; the executable radius interface uses $k\ge3$.

At $k=12$ and $14$, Haar survival is zero despite respectively four and six
valid times. Under uniform counting on the natural clocks $12$ and $14$,
the same sets have probability $1/3$ and $3/7$. The last actual mass factor
in (1) is zero, rather than $5/6$ or $6/7$.

More generally let $\delta=p/d$ in lowest terms and $V\ne\varnothing$.
Every endpoint of every one-runner safe interval lies on the clock

\[
                     N=d\operatorname{lcm}_{v\in V}v.     \tag{9}
\]

The intersection is a finite union of closed intervals, possibly singleton
intervals; a nonempty one has an endpoint from these native intervals.
Consequently weak loneliness for this fixed source is **equivalent** to
positive survivor count on clock (9). This finite decision rule can be
expensive and makes no uniform small-clock claim. A different failed clock
is only an uncertified grid: speeds $(1,2)$ at radius $1/3$ have witnesses
$1/3,2/3$, while clock two has no safe point.

There is also an exact way to keep a boundary witness in a measure statement.
For a fixed nonempty finite source put $M=\max V$. Then

\[
 \exists t\ \min_v\|vt\|\ge\delta
 \quad\Longleftrightarrow\quad
 |S(\delta-\varepsilon)|\ge\frac{2\varepsilon}{M}
 \text{ for every }0<\varepsilon<\delta.                 \tag{10}
\]

Forward: distance to an integer is one-Lipschitz, so a circle arc of radius
$\varepsilon/M$ about a weak witness remains safe at the relaxed radius.
Reverse: the asserted positive masses imply nonempty compact relaxed safe
sets. Along any $\varepsilon_j\downarrow0$ they are nested; compactness of
the circle gives a point safe at radius $\delta$. It is enough to check
such a sequence. The lower bound is stronger than needed for this reverse
direction. It is an exact conditional target, not an independently supplied
lower bound for arbitrary thirteen-speed inputs.

## 6. Certificate contract and bounded evidence

[Script](../../04-computation/experiments/lonely_runner_conditional_survival_20261011e.py)
and [saved output](lonely_runner_conditional_survival_20261011e.out).

`survivor` binds the ordered distinct positive integer speeds and exact
rational radius to the complete closed interval set. `audit_survivor`
recomputes it before a consumer uses mass, an endpoint bill or a witness.
`atom_engine` independently checks every native endpoint and open atom.
`conditional_lift` retains both labelled shared-residue profiles and solves
the actual CRT pairs; `clock_survival` binds the resulting classes to the
original clock. Empty conditional support returns an undefined later factor,
not an invented probability. Exact bool/float aliases and forged certificates
are rejected.

Reproduce with:

```text
python -B -X utf8 04-computation/experiments/lonely_runner_conditional_survival_20261011e.py
python -B -O -X utf8 04-computation/experiments/lonely_runner_conditional_survival_20261011e.py
```

The universe is 704 continuous cases: all subsets of $1,\ldots,10$ of size
at most three, at radii $1/6,1/7,1/12,1/14$; 930 all-wrap pair comparisons;
all tight APs from 3 to 20 total runners; every pair of subsets of moduli
1 through 5 (3,844 CRT cases); 288 native clock cases; and 123 complete-clock
cases, including empty and isolated survival. The protected $(50,51)$
correction, both signs of pair dependence, wrong-grid and equal-count profile
hostiles are explicit. This is finite validation of proved formulas, not
evidence of an unproved universal product law.

Normal and optimized runs reproduce **13,842 checks** and the saved output.
LF-normalized SHA-256:
`cb011c1ac0c7a39c1f6eaa2796da88b21474fca475855fecfdb66cd204fd199c`.
