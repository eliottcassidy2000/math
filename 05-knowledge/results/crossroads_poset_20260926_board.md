# Beyond 27, finite-poset boundaries, and a carrier that retains each integer

**Current synthesis, 2026-09-26. PROVED / FINITE-EXACT / CITED submitted
claims / OPEN Collatz.** Two elementary results are independently audited:
[THM-4502, 41-tail growth family](../../01-canon/theorems/THM-4502-collatz-41-tail-growth-family.md)
and [THM-4503, parity posets and arithmetic height selection](../../01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md).
The attached global poset results remain claims under review here. Their
internal arguments were examined in detail without certifying every cited
input or either asserted quantitative rate. No proof of Collatz or of the
full 1/3–2/3 conjecture is claimed.

## 1. What the session changed

1. An explicit convergent family containing 27 has arbitrarily long initial
   growth, a ternary phase recursion, Theta(log X) occurrence, and a
   zero-dimensional dyadic closure. It extends the old `4n+1` comb.
2. Expanding fixed-count parity words are exactly extensions of a two-chain
   poset. Exact carries recover their arithmetic starting residues. Swapping
   incomparable events changes the starting integer.
3. Six events involving 27 and 31 show that a height cut can destroy the
   poset law, convex support, and an XYZ covariance inequality. An exact
   total-variation formula restores uniformity over many residue periods.
4. A full-support generating function retains every fixed positive integer.
   Its missing decay estimate is explicit; strict decrease at every step
   is already false. This is an exact reduction, not a decay proof.
5. Actual segments refute the literal endpoint-free HYP-9161. An endpoint-
   corrected or whole-injective-orbit version remains open. A separate
   peak-count proof gap was repaired without changing its asymptotic bound.

Full lanes: [family and bridge](crossroads_poset_20260926_bridge.md),
[fixed-integer audit](crossroads_poset_20260926_integer.md),
[finite-poset barrier audit](crossroads_poset_20260926_barrier.md), and
[Kahn–Saks width audit](crossroads_poset_20260926_width.md).
The [reproduction audit](crossroads_poset_20260926_audit.md) records universes,
hostiles, independent checks, and unresolved dependencies.

## 2. Inheritance, portfolio, and the board after the experiments

The anchor is the user's fixed-integer Collatz problem. The orthogonal niche
is conditioning a uniform extension law by arithmetic height. The wildcard
is the finite Fibonacci boundary, including an exact 233 connection.
LRC(14) remains the repository's broader open anchor; no transfer to it is
claimed without a concrete object map.

Closest mechanisms: THM-4501's actual orbit-family realization and THM-4500's
auxiliary Bellman contraction plus a separate global policy construction.
Canonical hostiles: the all-odd -1 address, positive cycles of nearby affine
maps, and the new 27/31 conditioning example. Corrected near miss: equating
finite-word freedom with freedom to change the word of one fixed source.
The neglected sidecar is the law on sources, alongside exact height and carry.

| Live concept | What the new object adds | Surviving obligation |
|---|---|---|
| recursive 27 families | exact convergence, arbitrary growth, logarithmic frequency | distinguish forall k exists n from exists n forall k |
| poset extensions | exact width-two prefix-slope encoding | retain source law; many words do not imply large width |
| affine carry | exact source decoding and swap displacement | preserve this coordinate under every rearrangement |
| finite boundaries | attained extrema survive where averages vanish | later time is allowed; a changed source needs justification |
| atomic weighted mass | positive mass for each integer, uniform tail | decay at fixed weight parameter without normalization |
| landing transport | linear average on truncated actual segments | retain endpoints or prove an entire-orbit inequality |

The bridge and lane notes specify source, target, map, preserved predicate,
destroyed data, sidecar, and cheapest hostile test for each connection.

## 3. A family and the appropriate meaning of fractal recursion

There are two useful meanings of "the next 27." Under the property that
each new start takes longer than every smaller start to first fall below
itself, the exhaustive census through one million gives:

| start | shortcut first-descent steps |
|---:|---:|
| 2 | 1 |
| 3 | 4 |
| 7 | 7 |
| 27 | 59 |
| 703 | 81 |
| 10087 | 105 |
| 35655 | 135 |
| 270271 | 164 |
| 362343 | 165 |
| 381727 | 173 |
| 626331 | 176 |

Thus **703 is the next 27 under this precise record criterion**. None of
the later records visits 27 or41, so this family is not defined by sharing
their tail. This is a finite census, not a proved asymptotic occurrence law.
The records determine the plateaus of `m_k=min S_k`: for example m_k=27
at depths7..58 and703 at depths59..80. For fixed0<z<1,
`z^(m_k)<=Z_k(z)<=z^(m_k)/(1-z)`. Hence the least surviving source, not
just the number of surviving words, controls the atomic mass. The census
determines m_k through depth175, and only m_176>1000000 afterwards.

A different family permits an exact infinite theorem:

Let r_k be the least solution of `41*2^r=-1 mod 3^k`, of period
`p_k=2*3^(k-1)`. Then

    n_(k,t)=2^k(41*2^(r_k+p_k t)+1)/3^k-1

has k prescribed odd shortcut steps to `41*2^(r_k+p_k t)`, then the known
tail. Least-phase examples begin 27,291,12439,2173995791,1449330527,966220351.
Their exact shortcut stopping times are 70,75,82,101,102,103. These times
count the odd step and its first halving together. They are not claimed
to be record holders. For fixed k, the starts follow an integer affine
recursion; its first two rows are `4n+1` and `64n+35`.

Below X the union has Theta(log X) elements. Its dyadic closure consists
of the family, `(2/3)^k-1`, and -1, with only O(J) residue classes modulo
`2^J`. Hence upper box and Hausdorff dimension are zero. Exact recursive
arithmetic need not have positive fractal dimension. The larger positive-
slope language governed by THM-4495's entropy is a different set. All
members of the new family converge; its limit points do not convert this
into a claim about every positive integer.

## 4. The papers: actual contributions and audit boundaries

The [September 25 author announcement](https://igorpak.wordpress.com/2026/09/25/how-to-make-small-improvements-and-break-barriers/)
identifies the new manuscripts and their intended improvement. It confirms
provenance of the claims, not independent acceptance. Exact attachment
versions, page counts, SHA256 hashes, and visually checked pages are pinned
in the two audit notes.

The Aires–Chan–Pak–Panova barrier paper claims an absolute positive
improvement over `C=(5-sqrt(5))/10`, not the full 1/3 bound. It partitions
finite posets into bounded range, large width, and bounded-width/large-range
cases. In the first regime it combines a local-event probability floor,
near-equality rigidity of a three-element configuration, and a triple
selected at maximal expected-rank displacement. The visually checked local
improvement is `(D+1)^(-4(D+1))/4096`; 4096 is a denominator. Its enormous
global epsilon numeral is not certified here: the arithmetic derivation
is omitted in the manuscript.

The Aires Kahn–Saks manuscript claims balance tending to 1/2 as width tends
to infinity. Its structural alternative gives either nearly uniform orders
of a selected set or nearly uniform insertion positions into a selected
chain. Both alternatives matter. The final comparisons are coupled inside
one finite random extension. The internal proof audit found no concrete
gap. The cited selection, variance, logconcavity and gap results and the
asserted final quantitative rate have not all been independently re-proved.

The common useful concept is a joint law with boundary control, not merely
a balance constant. No external manuscript theorem is a dependency of
THM-4502 or THM-4503. The lane notes distinguish submitted claims, checks
of internal derivations, finite controls, and elementary proved additions.

## 5. The exact bridge and its transfer boundary

For a odd events and b even events, take two chains `A1<...<Aa`,
`B1<...<Bb`, with `A_(m_j)<B_j`, `m_j=min{i:3^i>2^(i+j)}`. In the feasible
endpoint range, extensions are exactly words whose prefix multipliers
exceed one. Its width is at most two. Large-width machinery does not become
applicable merely because it has many extensions.

The sidecar is `T_w(n)=(3^a n+C_w)/2^ell` and
`r_w=-C_w*3^(-a) mod 2^ell`. A legal swap 10 to 01 at j changes the source
residue at binary digit j. Balancing different extensions thus balances
different inputs. A fixed-input argument must retain or compensate this.

At ell=6,a=5 the four least residues are 27,39,47,31. Cutting at31 selects
the two extreme insertion positions. Their common comparisons still allow
all four, so no poset on these six labels has exactly the selected extension
law. Its order-statistic XYZ covariance changes from -5/32 to +1/4; its
support is nonconvex. This refutes a transfer, not a poset theorem.

The positive repair is exact. With m selected classes modulo N, X=qN+s,
and t selected classes in the final partial period,

    TV(height-selected law, uniform law)=t(m-t)/(m(qm+t)) <= 1/(4q)

for q>=1. Any comparison balance loses at most TV. The poset correspondence
therefore transfers when many full periods fit below X. It does not
automatically survive ell comparable to log_2 X, or ell tending to infinity
for one fixed source. This arithmetic boundary is now explicit.

## 6. Replacing the supposed computability obstacle

Every iterate of one integer is already computable. The missing fact is
that the search for descent terminates, not an algorithm for another step.
A finite program can produce unbounded integers. Positive cycles of 3n-1
and 5n+1 are direct hostile controls against arguments using only determinism,
positivity or computability.

The integer note examines the assumptions individually: real versus well-
founded ranks, boundedness versus cycle classification, pointwise versus
averaged decay, finite binary support versus a computable infinite address,
fixed versus moving height, and the quantifiers required for first descent.
Compatible least residues satisfy `r_(k+1)=r_k+e_k2^k`. They represent a
nonnegative integer iff they stabilize, equivalently iff
`limsup r_k/2^k<1/2`. The threshold is sharp. For already fixed n this
stabilization is automatic; proving it again gives no descent information.
An infinite bad branch must be proved unable to stabilize at a positive n.

A common horizon for all n is unnecessary and false. The target is
`forall n>1 exists k(n):T^k(n)<n`, equivalently one common finite horizon
for each bounded set of starts. Finite bounded-height compactness is valid;
an unbounded moving height does not provide it.

The exact no-descent carrier is the residue plus

    H_w=min_(j:3^(o_j)<2^j) floor(C_j/(2^j-3^(o_j))).

The surviving starts are exactly the positive members of the cylinder at
most H_w. This retains the additive correction instead of replacing actual
no descent by a slope approximation.

The strongest replacement here is a summable atomic law:

    S_k={n>=2:T^j(n)>=n for 1<=j<=k},
    Z_k(z)=sum_(n in S_k) z^n,  0<z<1.

At finite k this is an explicitly computable rational function: sum the
finite or infinite geometric series of the capped cylinders. The tail
above H is uniformly at most `z^(H+1)/(1-z)`. Consequently

    Collatz holds iff lim_k Z_k(z)=0 for one fixed z in (0,1).

Each integer has positive mass, so this limit would exclude every individual
counterexample. The exact map to the old density is

    lim_(z increases to1) (1-z)Z_k(z)=W_k/2^k.

That normalization erases each fixed atom. The replacement retains the
missing information, but leaves a hard explicit decay problem. It is an
equivalent reduction, not a proved approach to its required bound.

Hostile tests already rule out the easiest contraction: `Z_2=Z_3=z^3/(1-z^4)`.
At height27 the whole surviving mass is `z^27` at depths7 through58, before
first descent at59. A viable potential must tolerate delayed progress.
An exponential bound would imply first-descent time O(n), stronger than
needed; it must not be smuggled in as an assumption.

## 7. Returning to 233, Fibonacci structure, and boundary terms

In the finite Fibonacci poset on N labels, `i<j` iff j-i>=2, extensions
are matchings of the N-vertex path: disjoint adjacent swaps. There are
F_(N+1) of them. Thus N=12 gives **233** extensions. Their matching
indicators are 11-bit strings with no adjacent ones, exactly the legal
Zeckendorf strings for 0 through232 using weights F2,...,F12. This is a
representation map preserving a local exclusion rule. No such map to the
earlier tournament with 233 Hamiltonian paths is established; that
tournament sequence is not Fibonacci.

The boundary calculation is also exact. The probability of swapping edge
j is `p_j=F_j F_(N-j)/F_(N+1)`; expected-rank displacement is
`d_j=p_j-p_(j-1)`. On the symmetric interval N=2m+1,

    |d_i|=F_(2|i|)/F_(2m+2).

All normalized finite moments vanish, while the maximum tends to
`(3-sqrt(5))/2` and total absolute displacement tends to sqrt(5)-1.
Fixed central comparisons approach the BFT constant, while the extremal
boundary comparison remains better balanced. The paper deliberately uses
an attained finite maximum before the limit can erase it.

This gives a precise proof-design analogy for the earlier Bernoulli B1
idea: endpoint corrections can survive while bulk symmetry cancels.
No identity connecting Bernoulli numbers to Collatz carries is proved here.
The actionable point is to retain finite endpoint or atomic information
until its desired consequence is established. For Collatz a witness moving
to later time is allowed; changing the starting integer is the unsupported
step. These are different quantifier issues.

## 8. Corrections and strongest next targets

HYP-9161 quantified over arbitrary positive orbit segments without an
endpoint allowance. A critical hover followed by logarithmically many
halvings gives Omega(L) average multiplicity despite theta L=O(log L).
The m=96 example has99 dippers and4 landing positions. Its whole length
is O(k), so an O(k) boundary term survives. The hypothesis and current
routes now preserve the correction. The record-peak estimate survives
a suffix proof with the stronger -1-bit barrier. No established thinness
theorem is refuted.

The next targets, ordered by what is now concrete:

1. **Atomic-mass decay.** Find an adaptive block inequality or an auxiliary
   arithmetic potential controlling Z_k at fixed z. It must pass the 27
   plateau and the certified family without assuming termination.
2. **Weighted extension inequalities.** Retain C_w,r_w and height weights;
   derive a defect sharper than TV when q=0. The six-event covariance
   reversal is the mandatory first test.
3. **Fixed-source certificate posets.** Randomize the order of verifying
   arithmetic facts about one trace, or schedules of commuting carry
   operations. This preserves the source, but balance concerns certificate
   order, not automatically descent. Data-independent circuit posets
   cannot distinguish starts; arithmetic labels must enter an implication.
4. **Endpoint-stable landing transport.** Prove or refute a bound with
   `+C'k`, or over an entire injective orbit. Success could improve the
   logarithmic thinness factor, but would not itself exclude divergence.

The exact objects replace an ambiguous computability requirement with
specific missing inequalities and reproducible hostile cases. The open
problem remains the pointwise dynamical implication, not finite encoding.

## 9. Incoming work integrated without changing its probability space

During close-out, [THM-4504, Moran function and families of27](../../01-canon/theorems/THM-4504-families-of-27-moran-function-of-the-inverse-tree.md)
arrived from a concurrent session. Its exact fair-word martingale and
finite integer-window comparison are useful additional carriers. A uniform
parity word is still not the law of one fixed source. The backward-tree
growth, path-record count, delay constant and finite branch fractions must
retain their stated model or empirical scope. They are not imported here
as actual asymptotic orbit theorems.

The independent integration audit repaired three precise failures: the
finite hitting conditional is undefined when the hit event is empty;
`sup M_j>=W` need not give a finite hit on every address, so the open set
must be defined using an actual finite hit; and the glide count's factor
`2^(k-m)` changes the height exponent to `1-(1-h)m/k`. Fair-word almost-sure
equivalence does not justify a pointwise identification of the hitting
sets. These repairs preserve the valid limiting probability identities.
See [the Moran audit](crossroads_poset_20260926_moran_audit.md).

The incoming landing audit independently repaired the peak-denominator
gap and the earlier climb cartoon. We retained those corrections, added
the stronger suffix proof and the actual average-multiplicity refutation,
and repaired the negative-b scope of its one-bit band lemma. The ordinary
positive Collatz bound survives; general signed maps require high windows
or a separately charged finite core. No chronological log was used as a
replacement for these current proof statements.

There is a useful exact synthesis beyond correcting the limit: the greedy
word in the supremum counterexample is the **same ceiling-critical word**
used in our finite hover/drop construction. Choose odd when M_j<2 and even
otherwise, starting M_0=1. Then `1<=M_j<3` and the odd count is
`o_j=ceil(j log_3 2)`. This word is computable by integer comparisons.

**PROVED corollary of THM-4476:** its unique compatible 2-adic start is
irrational over Q, although it is computable. Suppose instead z=p/q with
q>0 odd. Multiplication by q preserves parity, so
`q T_1^j(z)=T_q^j(p)`, an actual integer orbit for the map with odd constant
b=q. Each odd step adds `1/(3M_j)` to the normalized carry h, hence
`o_j/9<=h_j<=o_j/3`. Its scaled orbit is
`x_j=M_j(p+q h_j)`, eventually positive and Theta(j), with
`|x_j|<=3|p|+qj`. The irrational odd frequency excludes repetition or a
visit to zero, since either would make parity eventually periodic with
rational frequency. Thus there are linearly many distinct integer orbit
values below X, contradicting THM-4476 for the fixed odd b=q. This rules
out every rational element of Z_2, in particular every ordinary integer.
Consequently the computable least nonnegative representatives of this
branch tend to infinity. Root and flow independently audited this corollary.

Every finite hover is arithmetically realizable; this particular infinite
hover is rigorously excluded by an orbit-coupled growth bound. That is a
concrete successful version of the requested fixed-integer test. It does
not extend automatically to all infinite bad words, whose multipliers
need not stay bounded. The corollary is inherited from thinness, not a
new proof of the full conjecture.
