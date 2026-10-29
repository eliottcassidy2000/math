> Historical source snapshot. Original: https://github.com/eliottcassidy2000/math/blob/201fa96604d51171909318cc482220e00df54c38/05-knowledge/results/atomic_tail_reset_20261029/ATOMIC.md . Relative links resolved to that pinned source. Claims retain their original status; no full74 search rerun in Session106.

# Atomic survivor bounds and the price of a failed return

**Status: PROVED scoped inequalities; FINITE-EXACT minima through sixty
returns; CONDITIONAL moment bounds; OPEN unbounded escape and reset payment.**
Same-agent proof review and two independent search algorithms, not an
independent human/agent proof audit or a literature-priority claim.

## 1. Inherited language, actual sources and two clocks

For positive odd n, perform the native odd Collatz step
`n -> (3n+1)/2^a`, with a the exact valuation. The cheaper ternary queue
starts at Q=1 and updates `Q -> Q-1+3*[3 divides a]`; stop at Q=0.
Universal termination of this auxiliary queue is inherited. Its five
strictly growing complete return words are:

|letter|native valuations|odd clock J|division clock D|carry C|source guard|
|---|---|---:|---:|---:|---|
|0|1|1|1|1|3 mod 4|
|1|3,1,1,1|4|6|179|125 mod 128|
|2|3,3,1,1,1,1,1|7|11|16177|1357 mod 4096|
|3|3,1,3,1,1,1,1|7|11|12289|2269 mod 4096|
|4|3,1,1,3,1,1,1|7|11|9697|2877 mod 4096|

The actual endpoint is `(3^J n+C)/2^D`. Distinct carries are distinct
branches even when their clocks agree. Appending v after u gives
`(J,D,C)=(Ju+Jv,Du+Dv,3^Jv Cu+2^Du Cv)`. The exact source guard is
`n=(2^D-C)3^(-J) mod 2^(D+1)`; the final oddness bit is necessary.

Let E_k be the positive odd sources whose first k queue returns all grow;
E_0 is all positive odds. Write s_k=min E_k and m_k=(s_k+1)/2. The
Mersenne source `2^(k+1)-1` proves E_k nonempty for every finite k.

Put m=(n+1)/2 and give this **original source** probability

    w(m)=1/(2*4^floor(log2 m)).

These weights sum to one; denote the resulting measure by nu. They are
not reassigned to current endpoints. A residue class m=r mod 2^L has mass
`4^(-L)+[r>0]w(r)`, with 0<=r<2^L. Thus the prior packet proves

    nu(E_k)=q^k+A_k,     q=1049603/4194304,

where A_k sums the named canonical source atoms of the completed-word
cylinders. A_k tends to zero iff the infinite intersection of E_k has
no positive integer. This is only the growing-return exceptional set,
not the full set of possible divergent Collatz orbits.

Sources: [atomic mass and reset boundary](https://github.com/eliottcassidy2000/math/blob/201fa96604d51171909318cc482220e00df54c38/05-knowledge/results/atomic_descent_wedge_20261029/DESCENT.md),
[local integer rank](https://github.com/eliottcassidy2000/math/blob/201fa96604d51171909318cc482220e00df54c38/05-knowledge/results/atomic_descent_wedge_20261029/INCOMING.md), and
[exact exceptional-address recurrence](https://github.com/eliottcassidy2000/math/blob/201fa96604d51171909318cc482220e00df54c38/05-knowledge/results/session97_overnight_arithmetic_carriers_20261029/controllers/EXCEPTIONAL_CANTOR_RECURRENCE.md).

## 2. New all-source inequality from occupied cells

The inherited inverse system has maps

    g_i(y)=(2^Di y-Ci)/3^Ji

on the odd 2-adics and disjoint images. Its compact attractor K meets N_L
normalized level-L cells, where

    N_0,...,N_10 = 1,2,3,4,6,9,13,18,24,31,40,
    N_L=N_(L-1)+N_(L-6)+3N_(L-11)  (L>=11).

Its dimension s and growth parameter lambda=2^s satisfy

    lambda^11-lambda^10-lambda^5-3=0,
    0.441<s<0.442,   N_L <= (17/8)lambda^L.

Here s is a dimension; s_k is the least integer survivor.

**Projection lemma.** For L<=k, the occupied level-L cells of E_k and K
coincide. Each k-letter cylinder has division depth at least k, and every
such cylinder admits an infinite 2-adic continuation. Changing its tail
cannot change its first L digits. Conversely every K address has such a
prefix. Positive integers realize every finite cylinder; no assertion is
made that a positive integer realizes an infinite address.

Set ell_k=floor(log2 m_k), so ell_k<=k. Then

    nu(E_k) <= B(k,ell_k)
       := N_k/4^k + (1/2) sum_(j=ell_k)^(k-1) N_(j+1)/4^j.       (1)

**Proof.** Cover E_k by its N_k occupied level-k cells in source-index
coordinates. Beyond the canonical representative each cell has weight
exactly 4^(-k), giving the first term. A surviving representative below
2^k must be at least m_k. In shell `[2^j,2^(j+1))`, every representative
has weight `1/(2*4^j)`, and its residue at level j+1 is distinct. The
projection lemma allows at most N_(j+1) such representatives. Summing
the shells proves (1). This argument covers **all positive integer
sources**, without an upper truncation on the input integer.

Since E_k decrease, replace B(k,ell_k) by
`min_{0<=h<=k} B(h,ell_h)` whenever that is sharper. In particular,
`A_k <= min B - q^k`. The subtracted term is the cylinder tail, not a
Haar-measure estimate for the whole survivor mass.

With R=lambda/4<1, (1) gives

    nu(E_k) <= (17/8)R^k + (17lambda/16) R^ell_k/(1-R).           (2)

Equivalently the second term is O(m_k^(-(2-s))), with an explicit constant
from the displayed floor. The generic unrestricted atomic tail gives
only O(1/m_k). The stronger exponent 2-s, approximately 1.558829, comes
from the arithmetic cell sparsity. It does not supply a lower bound on m_k.

## 3. Exact finite search, including the branch that breaks the shortcut

[least.cpp](https://github.com/eliottcassidy2000/math/blob/201fa96604d51171909318cc482220e00df54c38/05-knowledge/results/atomic_tail_reset_20261029/least.cpp) searches all five-branch prefix cylinders, keeping
the least original source r and its actual endpoint h. If a prefix has
odd clock J and division clock D, appending a branch of source guard r_i
and depth d requires the unique lift

    t=((r_i-h)/2)*3^(-J) mod 2^d,
    r'=r+2^(D+1)t,
    h'=(3^j*(h+2*3^J*t)+c)/2^d.                                (3)

Canonical source lifts never decrease. A valid Mersenne witness supplies
the initial incumbent; a branch with r' at least the incumbent cannot
improve it. There is no restriction to a sampled integer range or to two
of the five letters. This proves exhaustiveness of the branch-and-bound
search at the declared finite horizons.

The implementation is limited to k<=60. Sources are below 2^61. A prefix
at original source r has `h+1 <= (3/2)^depth*(r+1)`, hence h<2*3^60;
each tentative numerator is below `2*3^67+16177 < 2^128`. Nonzero source
lifts are computed only when their modulus is below the incumbent;
then J<=D<=59. The other branches need only small modular inverses. Thus
the stated unsigned integer widths are justified, not inferred from a
successful run. The original and native division clocks are both retained.

[audit_minimum.py](https://github.com/eliottcassidy2000/math/blob/201fa96604d51171909318cc482220e00df54c38/05-knowledge/results/atomic_tail_reset_20261029/audit_minimum.py) is an independent implementation:
it composes the total affine carry in arbitrary-precision integers,
solves the complete source congruence, and switches to literal native
execution when a cylinder contains at most one possible smaller source.
It imports no producer code, authenticates each supplied upper bound
before using it, and checks the entire smaller-source universe.

|k|exact s_k|independent affine audit|
|---:|---:|---|
|40|406360948735|PASS|
|50|716331825496063|PASS|
|55|59501570680684543|PASS|
|60|1090945607561904127|PASS|

The first forty results agree with the inherited receipt. At k=55 the
minimizer uses letter 4 at position 32, then letter 1 at positions 35 and
36 (positions start at one). Previous minimizers through 54 use only
letters 0 and 1. Thus the earlier observed two-letter pattern cannot be
promoted to a safe search reduction. This is a new hostile example, not
a correction to the earlier explicitly finite observation.

At k=60, J=63 and D=65; the endpoint is
33844837030466085603018630125. Here ell_60=58, and (1) gives exactly

    nu(E_60) <= 1556803587 / 2^120 < 1.17121*10^(-27).            (4)

Every positive odd source below s_60 leaves this growing-return language
within sixty returns. This is not a claim that it reaches one within
sixty returns, or a new global Collatz verification threshold.

## 4. A quantitative target that would actually pay an excursion envelope

Each of the five growing blocks satisfies

    R(n)+1 <= (3/2)(n+1).

For letter 0 this is equality. For the other four, substitution of the
guard's least positive source proves the affine inequality, whose slope
is already smaller than 3/2; increasing n then preserves it. After G
growing retries the peak at their endpoints is at most `(3/2)^G(n0+1)`.
The first non-growing return cannot increase that peak.

**Conditional moment statement.** If independently

    m_k >= c*alpha^k   for all k, with c>0,
    alpha^(2-s)>3/2,                                          (5)

then A_k tends to zero and the original-source-weighted expectation
`E_nu[(3/2)^G]` is finite. Indeed (2) bounds `Pr(G>=k)` by a sum of two
exponentials; the tail-sum formula

    E[beta^G]=1+(beta-1)sum_(k>=1) beta^(k-1)Pr(G>=k)

converges for beta=3/2. The critical per-return alpha is about 1.297071.
This is a conditional estimate for one first-exit excursion. It is not
universal repayment below the **original** source, and restarting the
prior at each endpoint is not permitted.

The exact finite window scan includes every k with k+r<=60 and widths
r=1,...,12. Plateaus disprove a fixed positive per-return growth factor:
even width three has a ratio of one at k=32. Width seven is the first
tested width whose worst ratio exceeds the above payment threshold:

    min_(0<=k<=53) (s_(k+7)+1)/(s_k+1)
       =7468447/1048576 > 7.                                  (6)

**OPEN candidate, not a theorem:** `s_(k+7)+1 >= 7(s_k+1)` for every k.
It was selected after the declared finite window scan, not predicted
independently. A proof would imply (5), because
`7^(2-s) > 7^(3/2) > (3/2)^7`; the last comparison squares to the exact
integer inequality `7^3*2^14>3^14`. Failure beyond sixty would be useful
evidence for variable windows, not grounds to omit that branch.

A constructive sufficient target is stronger and more informative: for
every n in E_(k+7), produce a legal marked predecessor y in E_k with
`y+1 <= (n+1)/7`. This forces (6) at every k. The map must retain the
original source label, native clock relation, and the future continuation
guard. A smaller integer with an unrelated itinerary does not suffice.

The recovered [THM-4555 collision/deletion theorem](https://github.com/eliottcassidy2000/math/blob/201fa96604d51171909318cc482220e00df54c38/01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md)
already gives smaller equal-clock sources for certain trailing-one
prefixes. The missing condition for the proposed target is that the
new prefix remains in the required growing language with the required
number of completed returns. The reset-2 obstruction and the sporadic
125 collision are inherited; they are not newly discovered here.

### Out-of-sample hostile search through return67

After selecting the candidate from horizons0..60, run the full five-letter
search for horizons h=61,...,67 with exclusive source bound
`7*(s_(h-7)+1)-1`. A witness below that bound would refute the next
seven-return inequality. [window_hostile.cpp](https://github.com/eliottcassidy2000/math/blob/201fa96604d51171909318cc482220e00df54c38/05-knowledge/results/atomic_tail_reset_20261029/window_hostile.cpp) finds
none in all seven additional windows; [window_checks.json](https://github.com/eliottcassidy2000/math/blob/201fa96604d51171909318cc482220e00df54c38/05-knowledge/results/atomic_tail_reset_20261029/window_checks.json)
retains every bound and node count. A positive control with the horizon55
minimum just inside the bound recovers the actual letter4 witness.

These are lower bounds, not newly computed exact minima beyond60. In
particular,

    s_67 >= 7636619252933328895.

The implementation has horizon at most67 and supplied source bounds
at most2^63. The actual prefix endpoint is below
`2^63*(3/2)^67`, so every tentative numerator is below2^115; all
source lifts are guarded before shifting. This gives an independent
width argument for this enlarged finite universe.

With the resulting lower bound ell_67>=61, (1) gives

    nu(E_67) <= 1000514718963 / 2^134
              < 4.595*10^(-29).                              (7)

Thus A_67 is bounded above by the same rational minus q^67. All earlier
finite bounds remain valid. The additional tests support the candidate
but leave its unbounded quantifier OPEN. The arbitrary-precision
[bounded affine auditor](https://github.com/eliottcassidy2000/math/blob/201fa96604d51171909318cc482220e00df54c38/05-knowledge/results/atomic_tail_reset_20261029/audit_bound.py) separately checks the full
horizon67 exclusion; its result is stored in `bound67_audit.json`.

## 5. What remains unpaid

Even universal finite G would merely give a local non-growing return
after an arbitrarily large growing excursion. It need not finish below
n0. The source-preserving signed coordinate in
[WEDGE_GEOMETRY.md](https://github.com/eliottcassidy2000/math/blob/201fa96604d51171909318cc482220e00df54c38/05-knowledge/results/atomic_tail_reset_20261029/WEDGE_GEOMETRY.md) measures exactly that difference.
To obtain a natural-number descent certificate, either prove a paid
return below n0 with an independently bounded native clock, or construct
a source-preserving rank that absorbs and pays every reset. The inherited
local queue rank alone does neither.

The weakest next target is unbounded s_k, equivalent to A_k->0. The
quantitative next target is a uniform version of (5) or (6). The final
Collatz target additionally needs the reset/source comparison. These
three statements remain distinct and OPEN.

