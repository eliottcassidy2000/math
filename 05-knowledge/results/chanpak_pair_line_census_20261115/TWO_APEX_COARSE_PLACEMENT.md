# Two apexes give a physical coarse placement rule

**Status:** proved all-order normalized placement; finite-exact census
of the declared strong subclass through order 12. No exact target
selector, density-one claim, or Chan–Pak proof is asserted.

## 1. The literal source class and exact transfer

Use a transitive spine of size $m$ and two named exceptional vertices
$A,B$, fixing $A\to B$. A bit in each apex's word says whether it beats
that spine vertex. Require $A$ to beat the first and lose to the last
spine vertex; $B$ beats at least one spine vertex. Then the $A$-spine
part is strong, can reach $B$ via $A$, and $B$ can return to it.
Consequently all these sources are strong.

Use the inherited odd-cycle packing expansion. Every cycle uses at
least one apex, so the contribution types are:

1. No cycle.
2. One cycle using just one apex: it has one increasing spine run.
3. One cycle using both apexes: it has two increasing spine runs,
   possibly an empty $A$-to-$B$ run since $A\to B$.
4. Two disjoint cycles: they contain $A$ and $B$ separately, one
   nonempty spine run apiece.

For types 3 and 4, assign each spine vertex to neither run, to the
first run, or to the second; an unused choice contributes $z$.
Each run tracks its empty flag, parity, and eligibility of its
current last vertex. The first selected vertex must be beaten by
the departure apex, and the last selected vertex must beat the
arrival apex. This gives a 25-state exact two-run transfer in time
linear in $m$ per evaluation. For a shared cycle require odd total
spine selection size; for separate cycles require a positive even
selection size in each. The extra factors are $2z$ for a shared
cycle and $4z^2$ for a pair. These rules follow directly from the
cycle powers in the slot polynomial.

At order 12, the entire declared strong subclass has 261,888
presentations, 215,841 distinct $(H,Q)$, and 178,579 different
triangle-reader outputs $F=H+Q$ between 64,491 and 2,734,957.
The longest consecutive-odd run has length 15, and the run
starting at the minimum has length 1. This is only the stated
two-apex subclass, not the complete tournament census. At odd
source orders the reader remains constrained to multiples of 3.
The C++ transfer passed 220 direct independent subset-DP controls.

## 2. A Bernoulli rule for arbitrary coarse targets

Uniformly for every real $x\in[0,6]$ and all sufficiently large $m$,
a source in the declared class can be constructed with
$$
 \left|\frac{F}{4^m}-x\right|
   = O\left((3/4)^{m/2}\right).
$$
Moreover over all two-apex sources, $F\le6\cdot4^m+O(3^m)$.
Thus the limiting normalized range is exactly $[0,6]$, not just
dense in some subband.

**Proof.** At $z=2$, normalize the two-run assignment weight by
$4^m$. It is exactly a product random assignment with unused,
first-run, second-run probabilities $(1/2,1/4,1/4)$.
The remaining coefficients are 4 for the shared cycle and
16 for the two separate cycles. The empty and one-apex-cycle
terms in $Q$ are $O(3^m)$, and so is $H$: at $z=1$ the
two-run assignments have weight 1 and total $3^m$.
The shared-cycle assignment's selected-total parity is odd with
probability $1/2$. The probability that both run sizes are even
is $1/4+2^{-m-1}$, by the parity characters of one trial.
Start/end requirements only lower these probabilities, yielding
the proposed upper maximum after adding the lower-growth terms.

Now set $r=\lfloor (m-1)/2\rfloor$. On the first $r$ spine positions
set the $A$ bits to all 1 and choose the $B$ bits below. Set both
words to all 0 on the last $r$ positions. In the middle (of length
$L=m-2r$), use arbitrary bits and ensure at least one $B$ bit is 1
to keep the source strong. On the event that the prefix and suffix
both contain both run labels, all start/end constraints depend
only on the first selected $B$-run position: its $B$ bit must be
1. This is the same test for shared and split cycles. Missing
one of the four boundary run labels has probability at most
$4(3/4)^r$. The probability associated with the prefix test,
ignoring missing labels, is
$$
  \sum_{i=0}^{r-1}\frac{b_i}{4}(3/4)^i .
$$
Choose the bits so this is within $(3/4)^r$ of $x/6$:
the two maps $t\mapsto (3/4)t$ and $t\mapsto1/4+(3/4)t$ cover
$[0,1]$ by overlapping images. Invert one containing the current
target at each step; the unresolved remainder after $r$ steps is
at most $(3/4)^r$.

Use the stopping times $S,T$ for seeing both labels from the front
and the back, on their respective disjoint $r$-position guards. On
success, reveal only the first $S$ assignments and last $T$; all
eligibility and the prefix test are then resolved. The middle still
consists of independent assignments, and has positive length.
Its odd total parity has probability exactly $1/2$.
The split parity correction, conditioned on these revealed boundaries,
is at most $2^{-(m-S-T)-1}$. Averaging it costs only
$$
 O\left(2^{-m}
  \Big[\sum_{s=1}^r 2^s\Pr(S=s)\Big]
  \Big[\sum_{t=1}^r 2^t\Pr(T=t)\Big]\right)
 =O((3/4)^m),
$$
since $\Pr(S=s)$ and $\Pr(T=s)$ are at most $2(3/4)^{s-1}$
(the probability one label has not appeared after $s-1$ trials).
The split parity estimate follows from the two single-run characters
with mean $2^{-(m-S-T)}$ and the joint character with mean zero.
Thus the sum of the leading contributions to $F/4^m$ equals $x$
up to $O((3/4)^r)$; the lower-growth contributions are
$O((3/4)^m)$. Since $r=\lfloor(m-1)/2\rfloor$, the estimate follows. ∎

This yields unconditional approximate *absolute* placement: for
large prescribed $N$, choose the least even $m$ with $N\le6\cdot4^m$,
so $N/4^m\in(6/16,6]$. The construction gives an even-order source
with $|F-N|=O(N^{1-(1/2)\log_4(4/3)})$. The exponent is below 1,
but the error is still much larger than an integer unit. The finite
census makes the distinction visible. The next exact problem is
to turn the middle run assignments (which were used only as a
parity buffer here) into a lower-digit correction without consuming
the coarse placement. This recasts the desired return as a
three-choice carry problem with two independent run end tests.

## Reproduction

Run from this directory:
```
c++ -std=c++17 -O3 two_apex_band.cpp -o /tmp/pairline_two_apex
python verify_two_apex.py /tmp/pairline_two_apex
/tmp/pairline_two_apex 10
```
The maximum spine parameter 10 gives source order 12. The command
line cap is intended for this census, not arbitrary orders with
the fixed-width integer implementation.
