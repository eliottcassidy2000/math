# Two hub orders force a doubled ternary coefficient

The Wave 6 [two-hub frame proposal](https://github.com/eliottcassidy2000/math/blob/f66f56b/05-knowledge/results/chanpak_waves_20261113/README.md#user-content-wave-6-the-digit-structure-of-two-hub-counts)
counts Hamilton paths by extremal spine data and suggests compiling
positive smooth antichains into ordinary path counts. The frames
are indeed smooth monomials. They cannot be independently turned on.

## 1. Exact involution and the unavoidable second copy

Fix a transitive spine and hubs x,y. A path through x then y has
increasing spine runs S0,x,S1,y,S2. In a nonempty-extrema frame let

    a=max S0, b=min S1, c=max S1, d=min S2.

The four fixed join arcs are a→x, x→b, c→y, y→d. Suppose a>d
and c>b, with the four extrema distinct. The swapped order y then
x has valid join arcs for transformed frame (c,d,a,b). An ordinary
spine vertex i outside those four anchors sees

    [i<a]+[b<i<c]+[i>d]
      = [i<a]+[i>b]+[i<c]+[i>d]−1
      = [i<c]+[d<i<a]+[i>b]

available runs in both frames. Each run sorts increasing and the
anchors guarantee the joins, so one may transport the k-th legal
run choice at each vertex to the k-th legal run choice in the
other frame. This gives an explicit fixed-point-free involution
on the corresponding paths across the two hub orders. Every frame
with a triple-covered vertex falls in this case: its left anchors
b,d lie before that vertex and its right anchors a,c lie after.
Thus every triple-containing monomial 2^α3^β comes with an exactly
equal partner. In particular literal distinct-term antichains that
require a single 3-bearing frame are not the unmarked sum of such
an ordinary two-hub host. No conclusion about alternative
representations of the same numerical N follows from pairing alone.

## 2. Extra routes sharply strengthen the physical lower bound

Suppose any frame has β≥1 triply covered sites. Restrict to only
those sites, the four anchors, and the two hubs. In this induced
core the four locked anchors allow two complete ternary chambers:
3^β paths in each hub order. Neither hub is deletable in these
paths, because deleting x joins a backward pair of spine vertices,
and deleting y does the same.

There are further, disjoint families coming from one-hub insertion.
When using x alone, force the bridge a→x→b. The β triply covered
spine sites can be assigned freely to the two increasing runs.
There are also e_x other anchors strictly between b and a: the
left anchor d if b<d, and the right anchor c if c<a. Similarly the
bridge c→y→d gives e_y. Consequently e_x+e_y=2. The one-hub
Hamilton paths with these bridges number 2^(β+e_x) and 2^(β+e_y).

Into each such Hamilton path insert the other hub using any fixed
legal insertion rule. Deleting that other hub recovers the original
word, hence each lift is injective. Their lifted families are also
disjoint. If both hubs in a common path were nonadjacent, each
backward bridge would already appear next to its hub, so the
opposite deletion would create an illegal backward spine arc.
If the hubs were adjacent, the shared external predecessor and
successor would have to be both (a,b) and (c,d), impossible. Neither
lift enters either of the locked ternary chambers, since at least
one hub is deletable. Finally, inserting both hubs successively into
the increasing spine path gives a further path in none of these
families, because its spine order is globally increasing.

Every Hamilton path of the induced core can be extended
injectively back to the original tournament by ordinary insertion
of the omitted vertices. Therefore

    H ≥ 2·3^β + (2^e_x+2^e_y) 2^β + 1
      ≥ 2·3^β + 4·2^β + 1.

When e_x,e_y are 0 and 2 the binary term improves to 5·2^β.
## 3. The next debt scale is unavoidable

Number the β triple sites in spine order, and look at site j with
l=j−1, r=β−j. It makes one new forward-order frame for each hub.
For hub x, if the site beats x, replace the right anchor a with
the site; if x beats it, replace the left anchor b. For y do the
same at c or d. The replacing arc is therefore guaranteed by the
actual tournament. An incoming-anchor replacement gives 3 choices
for all l earlier triple sites and 2 for all r later sites; an
outgoing-anchor replacement interchanges those arities. The old
anchor displaced by the site is allowed in at least one run (end
run 2 for incoming replacements, end run 0 for outgoing ones).
All other fixed joins stay legal.

All these alternate frames are distinct: only one of the four
extrema is changed, and to a specifically recorded site. Their
swapped counterparts are distinct too by the involution. They
retain backward joins across both hubs, so none belongs to the
deletable lifts or the increasing-spine path of §2. Hence

    H ≥ 2·3^β + 4 Σ_(j=1..β) min(3^(j−1)2^(β−j),
                              2^(j−1)3^(β−j)) + 4·2^β + 1.

The middle debt term is 8·6^k−8·4^k for β=2k and
20·6^k−16·4^k for β=2k+1. This proves that the excess above
2·3^β is at least a constant multiple of (√6)^β. At β=2
the inequality gives 51; at β=3 it gives 143. This is
a quantitative obstruction to unmarked compilation, not merely
a budget for how many vertices a triply covered frame needs.

### The 34/45/69 control, for any antichain representation

For a positive 2–3 antichain sum the lowest 3-exponent is v3 of
the integer: every 3-exponent appears at most once, so its lowest
coefficient cannot be canceled by higher ones. Similarly its
smallest 2-exponent is v2. For 45 a largest 3-exponent at most 2
would give only one term 9·2^a, impossible. For 69 such an exponent
cap allows b=1 alone, or b=1,2; the latter with odd total forces
23=2^a+3, impossible. Thus *every* antichain for 45 or 69 has a
β≥3 term, which the H≥143 bound excludes from an ordinary
two-hub realization at those respective targets. (It does not
exclude the targets themselves from the spectrum.)

For 34 a cap β≤1 would allow a one-term power or a pair
2^a+3·2^c with a>c. In the pair v2(34)=1 forces c=1 and
17=2^(a−1)+3, impossible. Hence every antichain for 34 has β≥2.
The exact minimal induced core at β=2 is 105 (finite complete
check of all forced-arc orientations on the 8 named vertices).
Indeed one may order left anchors b,d by swapping the hubs, and
there are two orders for the two right anchors: then the other
four hub/anchor bits, the hub edge, and every hub/site bit are free.
Literal path counts over these cases give the minimum 105.
For β≥3 use 143 or larger above. Thus even the shifted
proposal H=1+2S with S an arbitrary 2–3 antichain fails for
H=69=1+2·34, as a literal selection of those frames.

There are also infinite families where the *raw* antichain
hypothesis necessarily fails while the unmarked counts themselves
are attained. A sum for 3^m can only be its single term (its smallest
b is m, and any larger b overshoots). In a sum for 5·3^m, some b is
at least m+1, since 5 is not a power of two. In sums for 13·3^m or
17·3^m, some b is at least m+2: after dividing by 3^m, if one used
only b=0,1 an odd two-term sum would have to be 2^a+3 (neither case
fits). But respectively 2·3^m, 2·3^(m+1), or 2·3^(m+2) already
exceeds these four target families. This applies when the invoked
b is positive, as required for the frame lower bound (the
3^m family uses m≥1). All these values themselves are tournament
path counts by ordinally joining copies of a directed triangle
to any of the small known 5,13,17 sources.

## 4. Actual strong family, asymptotically attaining the factor two

Take x→y and the following spine labels, where the bits record
whether x and y beat the corresponding spine vertex:

    10, 01, (00)^p, (11)^q, 01, 10,     p,q≥1.

The four outer labels make the required anchor bridges for both
hub orders. Hence there is a frame with β=p+q. The tournament is
strong: each hub lies on a backward-return cycle over a spine
interval, and the two intervals overlap and cover the full spine.

Its unmarked absolute count is

    H(p,q) = 2·3^(p+q)
           + 6(2^p 3^q + 3^p 2^q)
           + 14·2^(p+q) + 2(3^p+3^q) − 1.

Here is a direct transfer proof. Process spine vertices in increasing
order. For the fixed hub order x,y, track the last status in S0
(0 absent, 1 beats x, 2 loses to x), the last status in S1
(0 absent, 1 beats y, 2 loses to y), and whether S2 has started.
On reading label (u,v):

- add it to S0 freely, setting its status to u+1;
- add it to S1 when S1 has started or u=1, setting its status to v+1;
- add it to S2 when S2 has started or v=1, marking S2 present.

Reject states with last status 2 for either S0 or S1. Also reject
empty S1 if the direct arc of the chosen hub order is not legal.
The other order swaps the bits and the direct-arc criterion.

For either repeated letter 00 or 11, the transfer on these 18
states satisfies t(t−1)(t−2)(t−3)=0. In fact the diagonal
weight at a state is the number of actions that leave its flags
unchanged, and each state-changing action increases this number
by one: the action then becomes a self-loop, with all other
actions' availability unchanged (except possibly more actions
at starts, which only increase it further). Thus successive
factors t−j eliminate states in increasing diagonal-weight
order. Each positive block length is therefore a combination of
1^length,2^length,3^length, and the two-block answer is a linear
combination of i^p j^q for i,j in {1,2,3}. The exact 3×3 table
at lengths 1,2,3 is

    [157,369,913; 369,853,2085; 913,2085,5053].

Inversion of the two power matrices gives the displayed closed
formula; both the annihilator and the values are replayed.

For p=q=k it specializes to the literal infinite ray

    H = 2·9^k + 12·6^k + 14·4^k + 4·3^k − 1.

Taking p,q about equally large makes H(p,q)/3^(p+q) tend to 2.
Together with the lower bound this proves that the minimum ordinary
H among all two-hub tournaments with a β-site ternary frame is

    2·3^β + Θ((√6)^β).

Neither the factor two nor the correction's exponential scale can
be removed by a better choice of the two hub bit words. It does not
make the remaining debt tunable at a prescribed arbitrary target.
