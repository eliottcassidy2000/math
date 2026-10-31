# A polynomial-bit realizer for the closed interval conductor

2026-10-31. **PROVED WITH INDEPENDENT PROOF/IMPLEMENTATION AUDIT.** No unrestricted LRC claim.

## Input and result

Input: a rational closed interval I=[a,b], a<b; a rational radius0<r<1/2;
positive integer speeds v_1,...,v_k; rational phase shifts phi_i. Let

    B_i={t∈I: ||v_i t+phi_i||<r},  W=b-a.

Suppose each μ(B_i)<W and Σ_i μ(B_i)<=W. There is an exact algorithm returning
rational t∈I with ||v_i t+phi_i||>=r for every i. Its arithmetic-operation
count is O(k(1+log^+(W Q²V²))), where V=max v_i and Q is the lcm of the
positive denominators of r and all phi_i. All intermediate bit lengths are
polynomial in the total binary input length, so this is polynomial bit time.
No enumeration of O(W Σv_i) periodic teeth is used.

The hypotheses are decidable by exact rational arithmetic. Session122's
interval conductor implies them for k>=2: each bad measure is <=W/k<W.
Thus that all-large-tail existence theorem has a polynomial-bit rational
realizer, including its equality boundary. A core-safe interval remains
safe because the algorithm never leaves it. For k=1 at the exact single-gap
boundary use the separate elementary first-safe-endpoint formula.

## Exact measure and component oracles

Put h=2r and F_h(x)=h floor(x)+min(frac(x),h), valid also for negative x.
Then

    μ(B_i∩[l,u])=[F_h(v_i u+phi_i+r)-F_h(v_i l+phi_i+r)]/v_i.

This is an exact periodic antiderivative: shifting by r makes the bad arc
(0,h) modulo1. Open endpoints have zero measure.

If t is bad for speed i, n=floor(v_i t+phi_i+1/2) is its unique nearest
integer and its open bad component is

    ((n-r-phi_i)/v_i, (n+r-phi_i)/v_i).

Both oracles use a bounded number of rational operations and floor calls.
They do not enumerate any intervening periods.

## Equality repair by one actual overlap

First test a. If safe, return it. Write D=W-Σμ(B_i)>=0.
If D>0, use the density f=1-Σ1_(B_i).

If D=0, choose any i bad at a, and let e be the right endpoint of its bad
component. Since μ(B_i)<W and a is in that open component, a<e<b.
Test e. If safe, return it. Otherwise choose j bad at e. Necessarily j!=i.
The two bad components have an actual open overlap O contained in I, with
positive rational length. Compute it by intersecting their endpoints with I.
Then use

    f=1-Σ1_(B_i)+1_O.

Pointwise f<=1_(I\union B_i): on O at least two bad indicators are1;
outside O at least one bad indicator makes f<=0. Its integral is
D+|O|>0. Thus equality is converted either to an immediate exact boundary
witness or to a strictly positive, exactly integrable witness-search density.
No marginal witness is substituted for a common time.

## Positive-density bisection

Maintain J⊆I with ∫_J f>0. Test both endpoints for safety; return any safe
endpoint. Otherwise split J into equal closed halves at its rational midpoint.
The integral is additive (the shared point has measure zero). If the left
integral is positive, keep it; otherwise the right integral is positive.
Compute each integral from the k exact measure oracles and the interval O.
No numerical approximation or guessed margin is used.

Every bad-set event is t=(n±r-phi_i)/v_i and, in reduced form, has denominator
at most QV. Hence two distinct event points differ by at least1/(Q²V²).
The overlap interval's interior endpoints are bad-set events; clipping at a
or b introduces no new interior event of I. Continue until

    length(J)<1/(Q²V²).

Then J has at most one distinct bad-set event in its interior. On each side
of that event every bad indicator is constant. Positive ∫_J f implies
positive-measure safe times, because f<=1_safe. Therefore one of the two
open subintervals (or all J if no event) is safe. Its adjacent endpoint of J
is safe too: an open bad set containing that endpoint would extend into the
safe subinterval. Thus an endpoint safety test must succeed. A safe isolated
point alone would not imply positive integral; such equality cases are
handled by the earlier direct component-endpoint step.

The number of halvings is at most the least nonnegative d with
W/2^d<1/(Q²V²), which is O(1+log^+(WQ²V²)). Each stage makes O(k) oracle and
safety calls. Dyadic refinement adds at most d bits to the input endpoint
encodings. For example, if D_I is a common denominator of the input interval endpoints,
D_I·2^d·Q·product_i(v_i) is a common denominator for all interval-integral
values encountered at depth d. Its bit length is O(N+d), where N is the
total binary input length; numerators also have polynomial bit length.
Rational products, sums, lcm computation and floors therefore have polynomial
bit cost in N and d. The claim is polynomial bit time,
not a strongly polynomial arithmetic algorithm or a small universal constant.

## What is actually realized

Source: a complete rational core-safe interval with the proved measure
budget. Target: one rational common-time witness, checked on every original
speed. The map retains the shared clock, phases, strict bad/closed safe
convention and input interval. It forgets the rest of the safe set and makes
no promise of continuous dependence on input or maximum loneliness.

The bad-measure budget is sufficient, not necessary. An input can have common
safe times while violating this budget; the algorithm rejects that domain
rather than claiming it unsafe. Exact interval geometry is more informative
than nonemptiness alone. No claim of a realizability topos is made from this
algorithm by itself.

## Corollary: a finite certificate can be compiled into a selector

Consider a fixed, already verified finite attachment certificate with:

* rational complete core-safe covers;
* finite small-speed branches covering every integer below each stated cutoff;
* well-founded inherited-family edges that retain the original constraints;
* finite leaf witnesses or an exact one-tail first-arrival operation; and
* all-large leaves discharged by the Session122 conductor.

Replace each all-large existential leaf by the algorithm above. The resulting
program returns an exact common-time witness for every input in the certified
family and runs in polynomial time in the binary encoding of its arbitrary
attached speeds (with the finite certificate fixed).

Indeed all branch speeds on a finite small-speed path belong to a fixed finite
universe, so every intermediate cover and cutoff is fixed finite rational data.
An inherited edge is well-founded by hypothesis. The program takes finitely
many fixed-table decisions and eventually reaches either a fixed witness,
a one-tail rational first-arrival calculation, or the polynomial-bit conductor
realizer. Constraint promotion or padding only strengthens the safe set, so
its witness remains valid for the original row. Sorting a fixed number of
input speeds and testing membership/comparisons has polynomial bit cost.

This applies to the already certified Session122 eight-core/five-tail family,
with its inherited nine-core/four-tail and earlier finite families. It adds
constructive witness extraction, not a larger family or a new radius. The
terminal realizer and a complete direct dispatcher for this bounded-core family
are implemented and tested here; see [the family selector](FAMILY_SELECTOR.md).
The dispatcher uses universal completion to avoid enumerating inherited
certificate tables. The generic fixed-dimension ILP method is a separate
baseline and is not used in this proof.
