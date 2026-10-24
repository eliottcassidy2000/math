# Source-preserving recursive certificates and unbounded-depth return cylinders

**PROVED elementary guarded-word theorems / FINITE-EXACT bank comparisons
and policy tests / OPEN universal coverage and Collatz.**

There is an explicit, disjoint dyadic cylinder for every k>=0 on which the
first descent below the original source takes exactly2k+2 odd Collatz steps.
Each certificate uses at most three parameterized macro nodes. The cylinders
for k>=3 are wholly outside the preceding171-row/65-cylinder core bank.
Their membership test carries an unbounded integer counter but is directly
computable for every fixed input. This is a positive replacement for a
fixed-depth rule; it is not a proof that every orbit enters or is discharged
by these cylinders.

Write U(n)=oddpart(3n+1) on positive odd integers, and
tau(n)=min{j>=1:U^j(n)<n}, allowing infinity.

## 1. Inheritance and live coordinates

The closest proved mechanism is the affine precision lift, including its
clipped last division, in [the swap-lift note](reset_20260926_swaplift.md).
Its finite bank verifies171 odd coefficients q<=341 and yields65 complete
dyadic source classes. The canonical hostile here is27->41->31: the child
41 has descended but the original27 has not. The corrected near miss is
equating bounded proof-node count with bounded orbit duration. A macro
parameter can encode an arbitrarily long proved word.

The least-used sidecar is the original source retained through every call,
together with the exact affine carry. The older
[barrier atlas](collatz_procgen_20260922_barrier_atlas.md) already records
the danger of relaxing a rewrite system until it forgets the actual path;
no external theorem from that atlas is a dependency here.

The board has six coordinates: original threshold, actual current integer,
legal repeated word, total affine map, final valuation budget, and typed
child calls. The new source cylinders strengthen the final budget. The
nested27 family below shows why repeated-word compression alone does not
pay that budget. The companion [entry note](entry_20260927_entry.md) handles
the separate logical distinction between local-descent entry and entry to
a convergence-certified set.

## 2. Main all-height return-cylinder theorem

**PROVED.** For k>=0 put m=k+1. Let t=t(k) be the least integer t>=1 with

    2^t(8^m-5) > 9^m-5.                              (R1)

Let beta_k be the unique odd residue in[1,2^t) satisfying

    beta_k*9^m = 5 mod2^t.                            (R2)

For every positive integer b=beta_k mod2^t, the odd positive integer

    n=b*8^m-5                                        (R3)

has

    tau(n)=2k+2.                                     (R4)

The defining source cylinder is exactly

    n = beta_k*2^(3k+3)-5 mod2^(3k+3+t(k)).            (R5)

The least positive representative in(R5) is at least3, so there is no
small-source exception. The choice in(R1) exists for every k. It is a
sufficient uniform-in-b valuation budget; no assertion is made that it is
the smallest possible budget after exploiting the actual beta_k.

**Proof.** More generally, a positive n=a*8^k-5 with positive even a has the actual
word(1,2)^k, with endpoints

    U^(2j)(n)=a*9^j*8^(k-j)-5,
    U^(2j+1)(n)=(3a/2)*9^j*8^(k-j)-7.                (R6)

For1<=j<=k the even endpoints exceed n. For0<=j<k the odd endpoints
exceed n because their difference is at least(a/2)*8^k-2>0.
The valuation claims follow by substituting in3x+1; this is an actual
word, not just an affine formula with possibly illegal divisions.

In(R3), a=8b, so after2k steps the actual integer is

    y=8b*9^k-5.

The next valuation is1 and the next integer is12b*9^k-7, again greater
than n. Its successor is

    w=oddpart(b*9^m-5).                              (R7)

Here b*9^m-5 is positive, and(R2) gives valuation at least t. Thus

    w <= (b*9^m-5)/2^t.

Let D=2^t*8^m-9^m. Equation(R1) says D>5(2^t-1)>=0. Since b>=1,

    2^t*n-(b*9^m-5)=bD-5(2^t-1)>0.                 (R8)

Consequently w<n. Every preceding iterate was above n, giving(R4).
The integrality and oddness of beta_k follow from invertibility of9^m
modulo2^t and oddness of5. This completes the proof.

The true terminal valuation is

    v2(3U(y)+1)=2+v2(b*9^m-5).

Thus the final two steps consume at least t+3 binary divisions in total.
That distinction prevents the factor-of-two near miss described below.

## 3. Failure boundary, density, and fixed-input computability

The first budgets are t=1 for k=0..4 and t=2 for k=5..10. These give,
respectively, the five complete valuation classes

    v2(n+5)=3k+3, k=0..4,

and six narrower cylinders

    n=2^(3k+3)-5 mod2^(3k+5), k=5..10.

The next band is t=3, beta=5 for k=11..16. These are sufficient bands
obtained from(R1), not an empirical pattern extrapolated to all k.

**Exact boundary controls.** Keeping only the preceding budget can fail:

| k | b | actual v2(b*9^(k+1)-5) | source n | endpoint after2k+2 steps |
|---:|---:|---:|---:|---:|
|5|3|1|786427|797159|
|11|1|2|68719476731|70607384119|

In both cases every preceding step is also above the original source.
These refute the corresponding weakened selected-return guarantees;
they do not exclude a later descent.

**Corrected near miss, never a proved dependency.** The provisional
implication `y=3 mod16 -> U^2(y)<=(9y+5)/32` is false:19->29->11, while
the proposed upper bound is176/32. The universally valid denominator is16.
The missing bit is b mod4: b=1 mod4 guarantees the additional division
and repairs the /32 statement. The general repair is the exact budget
v2(b*9^m-5)>=t in(R2). The theorem above does not use the false bound.

**PROVED density statement.** Distinct k cylinders are disjoint, since
each has exactly v2(n+5)=3k+3. Put K_k=3k+3+t(k). Their union has natural
density

    delta=sum_(k>=0) 2^(-K_k).                        (R9)

To justify existence, the entire tail k>=L lies in n=-5 mod2^(3L+3).
Its upper density therefore tends to zero as L grows; finite unions
approximate the union from inside. Also

    (9^m-5)/(8^m-5) > (9/8)^m,

so(R1) gives2^(-K_k)<9^(-k-1). Consequently the omitted density after
the first L cylinders is less than1/(8*9^L).

This gives an exact arithmetic explanation for the density scale: the
final valuation budget compensates for the repeated factor9/8 expansion.
The common denominators are dyadic; the exponential density bound is
governed by9. Neither fact grants independence along a fixed orbit.

**FINITE-EXACT comparison, all-height consequence.** Independent
regeneration of the old171 rows and exact prefix comparisons show that
its65-cylinder union misses the entire cylinder n=-5 mod4096. All new
cylinders k>=3 lie in that missed cylinder. The first three new cylinders
are wholly contained in the old bank. Their overlap has density73/1024;
there are no partial overlaps, including for untested larger k.

For the first20 cylinders, the density added to the old bank is exactly

    2553380107527241 / 2^64.

The augmented density is

    6990313556829423891 / 2^65
      = approximately0.18947282861673337

among all positive integers (twice this among odd integers). The infinite
augmentation differs by less than1/(8*9^20). These numbers concern a set
with certified first descent, not a set already proved to converge to1.

**Fixed-input computability.** Given n, compute s=v2(n+5). Membership
requires s>=3 and s divisible by3; then k=s/3-1 is uniquely determined.
Compute t(k), b=(n+5)/2^s, and test(R2). No infinite search or choice of a
2-adic input replaces this fixed integer. The counter is unbounded as
n varies. A constant number of macro nodes is not a constant-bit or
constant-time certificate.

## 4. A sound recursive certificate grammar

An evaluated proof carries the original source N and a transfer

    x=(M*N+C)/D,

with the current actual odd integer x. A legal word with transfer
(p,d,c) composes by

    (M,D,C) -> (pM,dD,pC+Dc).                         (R10)

All divisibility guards are checked on the actual x. The source N is
immutable. At a terminal node, descent is discharged exactly when

    (D-M)*N > C.                                    (R11)

This terminal comparison is an accounting identity, not by itself a new
estimate. The new estimate is(R8), which pays it on the entire family.

The implemented primitive rules are:

* `Step`: one actual U step, with its exact valuation.
* `Repeat1(l)`: guard l>=1 and v2(x+1)>=l+1; transfer
  (3^l,2^l,3^l-2^l), encoding l valuation1 steps.
* `Repeat12(k)`: guard k>=1 and v2(x+5)>=3k+1; transfer
  (9^k,8^k,5(9^k-8^k)), encoding2k steps.
* `ToFixed(c)`: c positive odd and(3x+1)/c a power of2 greater than1;
  this certifies the actual one-step endpoint c, without claiming c<N.
* `Core(q)`: the actual x must lie in that bank row's source cylinder.
  The actual lifted endpoint and affine word are evaluated; x is not
  replaced by the smaller core3q.
* `Call(x,body)`: the child source must equal the current actual integer.
  Its finite child proof discharges descent below x; the parent retains N
  and subsequently checks its own terminal comparison.

Soundness follows by induction on the finite proof tree: each leaf is a
legal actual word, composition preserves(R10), and every discharged
conclusion checks(R11) for the correct source. There are no circular
proof references or assumed eventual descents. Parameters and proof-tree
sizes may be unbounded. A literal known finite descent is representable
using Step nodes; that elementary completeness fact supplies no method
for producing a proof for an arbitrary input.

For(R3), the proof has `Repeat12(k), Step, Step` (omit the empty repeat
when k=0), so at most three nodes certify the unbounded duration2k+2.
This is the concrete positive separation between orbit horizon and
proof-node count.

Hostile typed-call controls include27->41->31, where the valid child
41->31 fails the parent target27, and a forbidden call replacing actual27
with smaller9. A positive nested control is7->11 followed by the child
11->17->13->5, whose endpoint is below both original thresholds.

## 5. The fixed family containing27 still carries unpaid growth

**PROVED.** For n_k=4*8^k-5, k>=1, put l=4+v2(k). Then

    tau(n_k)>2k+l.                                  (R12)

After the first(1,2)^k block, the actual integer is4*9^k-5. The elementary
valuation identity v2(9^k-1)=3+v2(k) gives v2(4*9^k-4)=5+v2(k).
Hence the next l valuations are all1. Every one of those steps increases
the already larger value. Their endpoint is

    e_k=3^l(9^k-1)/2^(l-2)-1.                        (R13)

Thus a second compressed word can increase unpaid growth; no change of
macro counter or arrival at a locally descending residue pays the source.

For odd k, l=4 and

    e_k=(81*9^k-85)/4,
    a_k=v2(243*9^k-251)-2.

The very next step pays the original source precisely when

    (243*9^k-251)/2^(a_k+2) < 4*8^k-5.               (R14)

At k=3 this holds:2043 reaches14741 after10 steps, then691 at step11.
At k=1 the endpoint is161, then121, still above27.

**PROVED sparse-gate restriction.** If odd k<j both satisfy(R14), then

    j-k > (243/32)*(9/8)^k.                          (R15)

Indeed(R14) implies2^a_k>(243/16)(9/8)^k; the strict lower bound follows
by cross multiplication from1215*9^k>1004*8^k. Subtract the two divisibility
congruences243*9^k=251 mod2^(a_k+2). Since
v2(9^(j-k)-1)=3+v2(j-k), the gap is divisible by
2^(min(a_k,a_j)-1). Both necessary payment bounds then give(R15).
This restriction applies to this specified next-step gate, not to all
possible later certificates. In particular it does not prove that k=3
is the only successful exponent. Exact testing of odd k<=1999 finds only3.
For k=1 mod4, a_k=2 by reduction modulo32, and(R14) always fails.

## 6. A distinct completed family containing27

The [exponent-completion note](entry_20260927_families.md) constructs
integers k>=1,h>=0 with47*4^h=-7 mod3^(2k+1), and

    a=2(47*4^h+7)/3^(2k+1),  n=a*8^k-5.

Here `Repeat12(k), ToFixed(47)` is an actual two-node word. All members
other than27 are larger than47, and this word certifies exact first
descent2k+1. The special parameter k=1,h=0 gives27->...->47; the two-node
word must be rejected as a descent certificate for27. The fixed actual
47->23 suffix, of34 odd steps, repairs this case and gives its first
descent at37 steps. The test stores that suffix explicitly: three outer
nodes include34 child Step nodes, so its stored tree has37 nodes, not3.

This construction also admits a genuine parameter recursion:
U^2N(k,h)=N(k-1,h), followed by the fixed terminal endpoint47. The
parameter decreases while actual integers grow. Soundness comes from the
legal word and the final original-source comparison, rather than an
assumption that each recursive integer is smaller. This is a model for
an allowed replacement of an unsuccessful numeric rank. It does not
establish closure or coverage for the fixed-a=4 family of Section5.

## 7. Exact experiment and remaining target

Reproduce from the repository root:

    python 04-computation/experiments/entry_20260927_recursive.py
    python -O 04-computation/experiments/entry_20260927_recursive.py

The [retained output](entry_20260927_recursive.out) records:

* 1154 legal repeat controls and11122 rejected guards, over every odd
 source3..2047 and repeat lengths1..6, compared with literal U words;
* 800 new-cylinder controls, k=0..99 and eight coefficient lifts each,
 checking every earlier iterate against the original source;
* the two exact weakened-budget hostiles, and all-height bank disjointness
 checked against the entire containing residue n=-5 mod4096;
* 128 nested27-family controls;1000 odd next-gate exponents and8128
 pairwise valuation-congruence checks;
* independent terminal-family controls, including the27 target rejection;
* a baseline source-preserving policy on4095 odd sources3..8191, the first
128 fixed27-family sources, and128 Mersenne sources2^(2k)-1.

The baseline policy uses the old bank, maximal repeat words, and actual
one-step returns; it was not changed to prioritize the new cylinders.
All tested certificates close. Their largest outer macro counts are20,
63,196 in the three respective universes; largest selected endpoint times
are51,346,627. Policy endpoints can skip an earlier first descent inside
a macro, so those policy times are not automatically first-descent times.
Only the individually checked first-descent claims above carry that label.

**OPEN executable target.** The useful next question is whether a larger
guarded-word grammar admits a source-preserving reduction for every residual
input, with a well-founded parameter independent of assumed Collatz
termination. The existing proof establishes this for explicit infinitely
many cylinders and for an independently completed27-containing family.
It does not prove that all residual inputs enter those families, that a
finite node bound works universally, or that a local descent after earlier
growth repays the original source.
