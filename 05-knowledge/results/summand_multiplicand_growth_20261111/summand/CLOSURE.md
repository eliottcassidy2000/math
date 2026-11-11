# Strict closure is a mixed-seed ideal, with exact finite birth updates

Date: 2026-11-11. **PROVED:** arbitrary finite-seed closure, the joint
sum/product theorem, finite Apéry-state updates, logarithmic synchronous
generation depth, and eventual stabilization of the top three values with
an exact finite compiler. **FINITE-EXACT:** small hostile and positive
controls in `check.py`. No claim about shortest addition circuits or the
minimum stabilization time is made.

## Inheritance and the unused coordinate

The oldest recovered definition is the May3
`07-reflections/summand-graph-fermat-zeckendorf.md`: a strict witness is a
paired-parent cospan `{a,b}→a+b` with distinct positive VALUES a,b. The
faithful carrier is the pair, rather than ordinary reachability in its
two-arrow shadow. THM-2422,
`01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md`,
separates the fibre, shadow, and repeated generation; its corrected dyadic
law supersedes the old Fibonacci-depth guess. The May31 product-sum lane
`07-reflections/natural-mode-graphs-sum-product-s365.md` retains the same
distinctness convention for multiplication, with factors at least2.

The immediate proved input is S259,
`05-knowledge/results/session259_arithmetic_graph_growth_20261111/summand/RESEARCH.md`,
§2: for two seeds a<b, the closure is `{a,b}∪(a+b+⟨a,b⟩)`. Its §7 leaves
classification with three or more seeds open. The following ideal theorem
resolves that finite-seed classification. Section 5 also resolves existence
of S259's eventual top-frontier regime for arbitrary finite seed sets,
while leaving its minimum onset time and detailed transient open.

Board: paired parents / seed-support type / semigroup ideal / finite residue
state / synchronous depth. Closest mechanism: the two-seed translated
semigroup. Hostile: one seed alone cannot double itself. Corrected near miss:
the ordinary generated monoid includes values not yet constructible under
strict gates. Least-used sidecar: which original seed types occur in an
expression, before evaluating the expression as an integer.

## 1. The arbitrary-seed theorem

Let S be a finite nonempty set of positive integers, H(S) its additive
monoid including0, and C(S) its closure under sums of distinct available
values. Set

    I(S)=union_(a,b in S, a<b) [a+b+H(S)].

Then exactly

    C(S)=S union I(S).                                      (1)

Equivalently, besides the seeds themselves, the attainable integers are
those having **some** expression as a nonnegative integer combination of
seed values in which at least two distinct seed types occur. The expression
need not be unique; repeated uses of a seed are allowed.

*Proof.* The right side contains S and is closed under every strict sum.
Two distinct seeds produce one of the generating sums; a sum involving a
member of I remains in I because I+H⊆I. Thus C is contained in the right side.

Conversely, take any mixed expression. Start by adding its two largest
distinct seed values a<b. Their sum exceeds every seed value occurring in
the expression. Add the remaining seed copies one at a time, in any order.
The running value always exceeds the next seed, so every gate is legal.
The resulting value is the expression's value. Every member of I has such
a mixed expression. QED.

This argument also explains an incremental collision repair when an old
genealogy must be retained. If a generated x has a legal split x=u+v,
then the forbidden x+x can be replaced by

    (x+u)+v=2x.

Both gates are legal, since u<x and v<x+u. Reordering the original mixed
expression proves(1) without this extra gate; preserving an existing
expression can require it. These are distinct optimization objectives.

For S={4,7,11}, the existing seed11 is also4+7. An expression that tries
11+11 can instead use(11+7)+4. Its numerical coincidence with a seed does
not invalidate the mixed-support classification.

### Joint sum/product version

More generally let M be a commutative monoid, and let its finite seed set S
have a positive additive weight w: products satisfy w(xy)=w(x)+w(y), and
w(s)>0 for every seed. Let M(S) be the generated monoid and use the operation
only on distinct available elements. Then

    Cl_strict(S)=S union union_(a≠b in S) ab M(S).            (2)

The same proof works: begin with two distinct seed types having the largest
weights in the expression. Their product has weight greater than every
remaining seed, so equality at a later gate is impossible. Cancellation is
not needed for this argument. Distinct elements may have equal weights.
The identity also holds for an infinite seed set: each individual finite
construction or mixed expression uses only a finite subset of the seeds.

For ordinary multiplication of integer seeds≥2, take w=log, or simply note
that the product of the two largest seed values exceeds all remaining
seeds. Thus(2) is the exact strict multiplicative closure theorem as well.
It does not replace integer factorization by an unrestricted real logarithm:
the actual monoid elements and their seed expressions remain attached.
For seeds1 and2 in addition, exponentiation `a↦p^a` gives seeds p,p² in
multiplication, and both strict closures are the corresponding full ray.
The parent lane develops this actual prime-power embedding further.

## 2. One new seed has an exact update on the whole infinite closure

Adjoin a positive integer c not already in the seed set S. It may already
be in C(S), or may be a currently missing value. Write H⁺=H\{0}. Then

    H(S∪{c})=H+N0 c,
    I(S∪{c})=[I union (c+H⁺)]+N0 c.                       (3)

*Proof.* A mixed expression either already uses two different old seed
types, placing its old contribution in I, or uses c and at least one old
seed type, placing it in c+H⁺. Further copies of c are unrestricted.
Conversely each expression in the right side is mixed. QED.

H⁺ is the ordinary monoid, **not** the old strict closure. This distinction
is useful rather than problematic. From the single old seed2, the monoid
contains4 even though4 is not strictly constructible. Adjoining3 makes
3+4=7 attainable via(3+2)+2, avoiding the inaccessible old subexpression4.

If c already belongs to C(S) but is not a seed, it belongs to I. Then H
and I in(3) are unchanged, so the update adds no attainable values. If c is
missing, (3) gives every newly unlocked value without an infinite forward
search. For example adjoining4 to{2,3} fills4 and6; it cannot create1.

## 3. A fixed-dimensional exact compiler for successive seed births

Fix any seed m retained throughout subsequent births. For residues r mod m,
record the two minimum vectors

    h_r=min{n∈H:n≡r mod m},
    i_r=min{n∈I:n≡r mod m},                              (4)

using∞ for an empty class. These are Apéry-type minima. Because m∈H and
I+H⊆I, a residue contains precisely the arithmetic progression starting
at its minimum. Thus(4), together with the finite exceptional seed set S,
reconstructs the entire infinite closure:

    n∈C(S) iff n∈S or n≥i_(n mod m).                     (5)

No finite target cutoff is hidden in this statement. The state has2m
integer/∞ entries, and m need not change when a smaller new seed arrives.

Put p_0=m and p_r=h_r for r≠0; these are the positive-monoid minima.
For a new c, let ell=m/gcd(m,c). Equations(3) yield the finite updates

    h′_r=min_(0≤j<ell) [jc+h_(r−jc)],
    i′_r=min_(0≤j<ell) {jc+i_(r−jc),
                        (j+1)c+p_(r−(j+1)c)}.             (6)

All residue subscripts are modulo m. Terms with∞ remain∞. Taking another
ell copies of c preserves the residue and increases the value by ell c,
a positive multiple of m; it therefore cannot improve the minimum. This
proves that the finite range in(6) is exact.

Initialization from one seed m is `h_0=0`, all other h_r=∞, and every i_r=∞.
Adding any finite sequence of distinct seed values and applying(6) gives
the exact closure at each stage. This is a concrete one-node growth law,
with memory independent of the largest target one later asks about.

When gcd(S)=1 and |S|≥2, every i_r is finite. For a direct reason, the
residues of H form a submonoid of the finite cyclic group, hence a subgroup;
gcd1 makes it the whole group. Translating by one legal seed pair then
supplies every residue to I. Thus C is cofinite. For general gcd d,
`C(S)=d C(S/d)` and the same statement holds within the positive multiples
of d; other residues remain absent forever.

There is also an exact hole count in the gcd1 case. Let
e=|{s∈S:s<i_(s mod m)}|. Then

    #positive holes = (sum_r i_r)/m −(m−1)/2 −1−e.         (7)

In residue0 there are i_0/m−1 missing positive values below the ideal;
in residue r>0 there are(i_r−r)/m. Sum and restore the e exceptional seeds.
For the final hole, start at i_r−m in every residue and step backwards by m
while the value is positive and belongs to S. The largest positive value
remaining is the final hole; if none remains the conductor is1. This
extracts the exact conductor and genus from a finite state, for any finite
seed set, without claiming a simple scalar formula in all seed values.

## 4. Universal logarithmic synchronous depth

Now normalize gcd(S)=1, with |S|≥2. Let S_0=S and let S_(t+1) adjoin all
strict sums of members of S_t simultaneously. Write d_S(n) for the first
generation containing n. Let c be the exact conductor from§3,
m0=min S, M0=max S, and

    t0=max(0,floor(2c/m0)−1).

Every mixed expression for n has at most floor(n/m0) seed occurrences.
The sequential construction in§1 has one fewer gates, hence depth no larger
than that number. Seeds already have depth0. Therefore

    [c,2c] ⊆ S_t0.                                      (8)

**Interval-growth lemma.** If `[c,M]⊆S_t` and M≥2c, then

    [c,2M−1] ⊆ S_(t+1).                                 (9)

For each new z∈[M+1,2M−1], choose
`a=max(c,z−M)`, `b=z−a`. Then c≤a<b≤M. This supplies an actual legal
parent pair for every new target; the earlier interval is retained.
Iterating(9) gives

    [c,(2c−1)2^(t−t0)+1] ⊆ S_t,       t≥t0.              (10)

The upper frontier of ANY strict state satisfies
`max S_(t+1)≤2 max S_t−1`, since the largest distinct parents differ by at
least1. Hence, for n≥max(c,2),

    max{0,ceil[log2((n−1)/(M0−1))]} ≤ d_S(n)
       ≤ t0+max{0,ceil[log2((n−1)/(2c−1))]}.              (11)

In particular

    d_S(n)=log2 n+O_S(1)                                  (12)

as n grows through attainable values. For gcd d>1 apply this to S/d and
n/d. This supplies a bounded additive generation-time loss relative to the
doubling frontier for every finite seed set with at least two values.

The theorem does not say that every value up to the actual moving maximum
has already filled. Formula(10) supplies a uniform expanding interior with
a fixed positive fraction of the maximal allowed scale. It complements
S259's exact and potentially very long parity transient at the top edge.
No fixed-point iteration of a finite truncation is being mistaken for the
full infinite closure.

## 5. A consecutive block catches the moving maximum

**Block-chase lemma.** Let a finite positive state X contain the three
consecutive values A,A+1,A+2, and put M=max X. Its synchronous update is
`X′=X union {x+y:x<y in X}`. If B=A+2<M, let g be the gap from M to its
second-largest value. Then

    M′=2M−g,
    [M+A,M+B] ⊆ X′,
    δ′=M′−(M+B)=δ−g,       δ=M−B.                    (13)

Indeed all three old block values are strictly below M, so adding M to
each is legal. Since B is an old value below M, `1≤g≤M−B=δ`. Thus δ′ is
nonnegative and strictly smaller than δ. After at most δ synchronous
updates this translated three-block reaches the maximum. If it already
ends at M, the three pairs among M−2,M−1,M give the new top three values
`2M−3,2M−2,2M−1`. Therefore the top-three property persists, and

    M_(t+1)=2M_t−1                                      (14)

forever after that point.

The positivity of A matters: the new translated block lies strictly above
the old maximum. No growth rate or density heuristic is used here; the
tracked block has an actual three-pair construction at every step.

**Universal stabilization theorem.** Every finite seed set with at least
two distinct positive values eventually has its top three values spaced
by `d=gcd(S)`. From that point onward

    M_(t+1)=2M_t−d,
    M_(T+j)=(M_T−d)2^j+d.                              (15)

To prove this, first divide all seeds by d and work with gcd 1. Section 3
gives an exact conductor c. Let m0=min S, M0=max S, and

    T0=max(0,floor((c+2)/m0)−1).

Section 1 constructs each of c,c+1,c+2 by generation T0. Apply (13) to
this block. The bound `M_T0≤(M0−1)2^T0+1` gives the explicit finite bound

    T ≤ T0+M_T0−c−2
      ≤ T0+(M0−1)2^T0−c−1.                            (16)

This bound can be very loose. It establishes a finite transition for
every seed set, including arbitrarily long parity transients; it does not
identify the earliest transition or claim that all lower holes have filled.
For d>1, the top three mean `M−2d,M−d,M`; values outside dZ never occur.

### The exact shrinking-window compiler

Once any three-block `[A,A+2]` is known in a state, one need not retain
the growing state to execute the chase. Retain only

    M=max X,  W=M−A,  D={M−x:x∈X, 0≤M−x≤W}.           (17)

Thus D contains 0,W−2,W−1,W. If W=2, the top-three regime is already
reached. Otherwise let `g=min(D\{0})`. Since W−2 is positive and in D,
`1≤g≤W−2`. The exact update to the moving window is

    M′=2M−g,   W′=W−g,
    D′={u+v−g:u<v in D, u+v≤W}.                        (18)

*Proof of exactness.* The new window has lower edge
`M′−W′=2M−W=M+A>M`, so it contains no retained old vertices. A legal sum
of old values with deficits u,v has new deficit u+v−g. It belongs to the
new window exactly when u+v≤W. Both old deficits must then be at most W,
so every needed parent is already in D. No discarded vertex can enter
the new window. The pairs of deficit 0 with W−2,W−1,W show that D′ contains
its required three-block deficits W′−2,W′−1,W′. QED.

W decreases by at least 1 until it reaches 2. Thus the compiler uses a
shrinking subset of the initial W+1 possible deficits and at most W−2
steps. It returns an exact generation with a closed frontier recurrence;
it may return later than the first such generation because a different
three-block could reach the top sooner. Its loss boundary is explicit:
lower vertices and their genealogy are discarded, but they can never
affect the retained frontier during this chase.

## Reproduction and transfer contract

`python3 -B check.py` tests all nonempty seed subsets of[1,8] against literal
numeric-birth closure, checks all successive Apéry updates, and probes(9)
with explicit parent pairs. Selected synchronous states test(11) without
constructing exponentially huge late states. Every finite state in[1,12]
that contains a three-block tests(18) against its literal full update;
selected seed states then compare a complete compiler chase with literal
generations. `--write` regenerates the scoped receipt;
`python3 -O -B check.py` retains every explicit check.

Source→target: seed expressions→mixed-support ideal. Preserved: exact
attainability. Lost: genealogy, minimal circuit size, and synchronous depth;
§4 restores a depth bound, not an optimal individual genealogy.

Source→target: infinite ideal→residue minima. Preserved: complete membership
and exact one-seed updates. Lost: representation counts and alternative
parents. Sidecar: the seed exceptions and chosen modulus m. Hostile:
single-seed{2}, whose monoid contains4 while its strict closure does not.

Source→target: synchronous state plus a consecutive three-block→shrinking
upper deficit window. Preserved: exact maximum and every retained frontier
value until the block catches the maximum, then the exact maximum forever.
Lost: lower holes, representation counts, and the minimum possible onset.
Cheapest hostile: apply (18) at W=2 without stopping; that would incorrectly
use a translated block containing the forbidden pair M+M.
