# An exact excursion forest with ordered binary carries

**PROVED elementary identities / FINITE-EXACT forest implementation /
OPEN descent inequality, 2026-09-26.** No Collatz conclusion or literature-
priority claim. This continues the [crossing board](crossroads_crossing_20260926_board.md),
especially the finite-feature obstruction
[THM-4507](../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md).

## 1. Inheritance and the object

Anchor: a forest retaining one actual integer's excursion arithmetic.
Niche: Bernoulli/Fermat congruence projectors and ordered digit carries.
Wildcard: local tiling curvature, boundary closure, and signed differences.

Closest proved mechanism: the exact crossing identity and affine word map.
Canonical hostile: a positive shadow of an expanding rational cycle.
Corrected near miss: closed height excursions need not decrease the integer.
Least-used sidecar: the positions of all binary carries, not their total.
The concept board is the forest, affine curvature, ordered carries,
constructible residue projectors, tiling boundaries, and difference signs.
Method cards: retain local incidence through counting; controlled forgetting
requires a sidecar; attack a bound before extending it.

Use shortcut Collatz T(n)=(3n+1)/2 for odd n and n/2 for even n, and
height h(n)=floor(log_2 n). Every height increment is -1,0,or1. Stack-match
each downcrossing with the latest unmatched upcrossing. A matched interval
I=[i,j] has h(n_i)=h(n_j), with all strict interior heights larger.
Such intervals are disjoint or nested, forming an ordered forest.

Retain at every node:

- endpoints i,j and exact integers n_i,n_j;
- all chronological children, flat steps and boundary steps;
- the parity word w, length L, odd count a and integer affine carry C;
- each odd step's valuation and full ordered binary carry stream.

Unmatched downcrossings and still-open nodes are explicit boundary data.
A finite prefix is never completed by pretending its open excursions close.
This is an exact computable object for every finite orbit prefix, regardless
of whether the orbit eventually terminates.

For each node and every formally substituted x,

    F_w(x)=(3^a x+C)/2^L = A_w x+B_w.

Here C depends on the order of the word. For earlier block w and later v,

    C_(wv)=3^(a_v) C_w + 2^(L_w) C_v.

Consequently nested composition is exact. Deleting flat edges or forgetting
which child comes first is not an exact quotient. The inherited
[affine blueprint](collatz_blueprint_20260921_affine.md) already treats
the affine group and its noncommutativity; the forest adds a chronological
decomposition of an actual source, not a new abstract affine group.

## 2. A rational curvature and its composition rule

Every matched node contains an odd step, so B_w>0. Define

    kappa(w)=(A_w-1)/B_w=(3^a-2^L)/C.

Then

    F_w(n)-n = B_w(1+kappa(w)n).

Thus its exact descent certificate is kappa(w)<-1/n. Merely returning to
the same band gives no sign: the primitive excursion11->17->26->13 has
w=110, F_w(n)=(9n+5)/8 and kappa=1/5. The same word gives27->31.
If kappa<0, the fixed-point threshold is r=-1/kappa=C/(2^L-3^a): descent
occurs exactly above r, equality at r, growth below it. Integer admissibility
is a separate requirement.

For affine maps F(x)=A x+B and G(x)=C x+D with A,C>0 and B,D>0,

    kappa(F o G) = [B kappa(F)+A D kappa(G)]/[B+A D].

This is a positive weighted average. It preserves the minimum and maximum
of the constituent curvatures, but does not prove which side of the
source-dependent threshold -1/n contains the average. Pure even blocks
have B=0 and remain exact matrix atoms rather than receiving a fictitious
finite curvature.

This law is a rigorous point of comparison with angle deficit in a tiling:
there is a local scalar and a gluing law, but gluing does not discard the
boundary/source data. The laws are different (weighted averaging versus
additive angle sum); no map equating them is claimed.

## 3. The intrinsic tournament becomes an order, and sorting changes source

The pairwise observable most directly relevant to endpoint size is the
effect of swapping two affine blocks:

    (F o G)(x)-(G o F)(x)
       =(A-1)D-(C-1)B
       =BD [kappa(F)-kappa(G)].

It is independent of x. Distinct curvatures therefore give a transitive
tournament, while equal curvatures are genuine ties and the maps commute.
The vertices are arithmetic blocks; the observable is endpoint difference;
the orientation gauge is fixed by composition order; ties remain ties.
Redei's existence theorem adds no new arithmetic assertion to this order.

More seriously, the formally better order is generally unavailable to
one fixed integer. A length-L parity word with a odd letters has its unique
source residue

    n = -3^(-a) C_w mod 2^L.

All 2^L words give distinct residues. This follows inductively: each parity
step fixes one new binary input bit; on either branch the odd numerator
slope is invertible modulo the remaining power of two. Equivalently, an
integer has exactly one L-step parity prefix, and there are exactly 2^L
source classes.

Hence if wv and vw are different binary words, their source residues are
different. No fixed integer realizes both orderings. For example, the
formal blocks10 and110 yield words10110 and11010, requiring sources25 and11
modulo32 respectively. A comparison theorem can organize the order in
which facts are verified; it cannot reorder the actual arithmetic steps.

Preserved by the affine comparator: formal endpoint ordering.
Destroyed by sorting: source residue, chronology and height legality.
Needed sidecar: the original residue/carry address. The two-word example
is the cheapest decisive hostile test.

## 4. Three carry states, an unbounded path

The [Bernoulli/carry lane](forest_20260926_bernoulli.md) derives a stronger
sidecar than a fixed list of valuations. If n has binary digits b_s, define

    c_0=1,
    c_(s+1)=floor((c_s+3b_s)/2),
    c_s=floor((3(n mod 2^s)+1)/2^s) for s>=1.

All c_s lie in{0,1,2}. From each state, the b=0 and b=1 successors are
different, so the ordered carry path recovers every input bit. A
nonnegative integer corresponds exactly to a path eventually staying at0.
A three-state graph with an unbounded ordered path is therefore lossless;
the graph alone is not the integer.

The total K(n)=sum_(s>=1)c_s satisfies

    K(n)=1+3 s_2(n)-s_2(3n+1),

where s_2 is binary population count. It is unbounded inside every fixed
finite polynomial-valuation fibre selected by the construction in that
lane. This explicitly meets the escape condition from THM-4507.

But the total loses order. The growing transition51->77 has equal
(s_2,K)=(4,9), and the scalable family

    n_L=15*2^L+51,
    U(n_L)=45*2^(L-1)+77, L>=9,

has equal(s_2,K)=(8,17) at arbitrarily large endpoints. Thus
c log n+f(s_2(n),K(n)), c>0, cannot be globally nonincreasing. This
particular one-edge hostile does not by itself exclude adding an arbitrary
bounded correction; its height ratio is bounded. Retain that scope.

There is a sharper hostile to a natural way of weighting depth. Let
`D_n(z)=sum_s b_s z^s`. For the same padded family, direct binary expansion
gives the polynomial identity

    D_(U(n_L))(z)-D_(n_L)(z)
      =z^(L-1)(z-1)(z^4-1)+z(z-1)(z^4-z^2+1).

Both terms are positive for every real z>1, and the first tends to infinity
with L. Thus no positive multiple of D_n(z), plus a nonnegative logarithmic
height term and a bounded correction, is an every-step rank. The same
witness defeats any nonempty finite nonnegative combination of such
evaluations with at least one positive coefficient. This does not exclude
signed combinations, nonlinear functionals, or selected-return rules.
At z=1 the difference is zero, explaining the digit-total collision.

For odd iterates n_j and S_j=s_2(n_j),

    S_(j+1)=3 S_j+1-K(n_j).

The positive discounted carry budget is exact:

    S_0=sum_(j<k) [K(n_j)-1]/3^(j+1) + S_k/3^k,
    S_0=sum_(j>=0) [K(n_j)-1]/3^(j+1).

For odd positive n, c_1=2,c_2>=1, hence K>=3. Also n_j<=2^j n_0 gives
S_j<=log_2(n_0)+j+1, so the remainder vanishes independently of Collatz.
Thus this positive budget is not a termination proof. A useful forest
estimate must do more than restate its exponentially discounted tail.

## 5. A polynomial gluing law retaining carry positions

Let D_n(z)=sum b_s(n)z^s and C_n(z)=sum_(s>=1)c_s(n)z^s. The exact
one-step coefficient relation is

    z D_(Tn)(z)=3D_n(z)+1-(2-z)C_n(z)/z  for odd n,
    z D_(Tn)(z)=D_n(z)                    for even n.

For a chronological interval of length L with a odd steps, define

    B_w(z)=sum_(t odd-source) 3^(odd steps after t) z^t,
    A_I(z)=sum_(t odd-source) 3^(odd steps after t) z^t C_(n_t)(z)/z.

Both are polynomials with nonnegative integer coefficients. Direct
composition gives

    z^L D_y(z)=3^a D_x(z)+B_w(z)-(2-z)A_I(z).             (F1)

B_w depends only on the parity word; A_I retains the ordered carries of
the actual source. For adjacent earlier interval I and later interval J,

    A_(IJ)(z)=3^(a_J) A_I(z)+z^(L_I) A_J(z),
    B_(IJ)(z)=3^(a_J) B_I(z)+z^(L_I) B_J(z).              (F2)

Thus the same exact gluing law works across an arbitrarily deep forest.
These positive polynomials are a concrete boundary/bulk carrier, not an
assumed drift estimate. They retain information unavailable in finitely
many valuations or one scalar total.

At z=1, (F1) is the weighted digit-count budget. At z=2, B_w(2) is the
ordinary affine carry C_w, while the entire carry-position term vanishes:

    2^L y=3^a x+C_w.

This is the key obstruction to a direct positivity argument: positivity
of A_I does not by itself give a height gain at the endpoint z=2 where
the ordinary integer is read. Differentiating retains it:

    L*2^(L-1)*y+2^L D'_y(2)
      =3^a D'_x(2)+B'_w(2)+A_I(2).                      (F3)

The remaining problem is to bound this positional information across
resets in a way that constrains the affine endpoint, rather than merely
expressing it in another form. Equations(F1)--(F3) are exact independent
of termination. They provide a precise candidate object for that estimate.

## 6. An eight-tile space-time realization

The carry sidecar has an exact local two-dimensional realization. At one
shortcut step set p=n mod2, m=1+2p, c_0=p and let b_s be its input bits.
Use tiles with data(p,b,c,q,c') satisfying

    m b+c=q+2c',   q in{0,1}.

There are six tiles in odd mode p=1 (b in{0,1},c in{0,1,2}), and two
in even mode p=0,c=c'=0. Left/right edge labels(p,c) make p constant
along a row and propagate its carry. The top bit is b and bottom bit q.
The left boundary additionally imposes b_0=p,c_0=p; hence q_0=0.
Discard q_0 and shift the remaining output one place to obtain T(n).

Place digit s of time t at horizontal position s+t. Then bottom bit q_s
matches input bit s-1 at time t+1, and the discarded zero lies outside
the next row's diagonal left boundary. Finite binary support at the right
of each row is retained. Thus these eight tile types exactly implement
shortcut Collatz on a diagonal half-plane with a fixed finite binary top
word padded by infinitely many zero bits.
This is a coloured-square local constraint construction, not a claim that
these eight logical types are regular-polygon vertex figures.

Its virtue is explicit: the scalar source is replaced by a spatially
ordered boundary word, which is unbounded across integers and exact at
each forest node. Its limitation is equally explicit: a finite tile set
can support arbitrarily deep space-time diagrams. Showing that every
positive finite top boundary eventually reaches the1<->2 pattern remains
Collatz itself. The construction provides a carrier on which to seek a
boundary obstruction; it does not supply that obstruction.

## 7. The next precise target

The forest is now an implemented carrier rather than a proposed picture.
It retains an unbounded carry word at each event and exact composition
across nesting. A possible proof route is a source-preserving rule selecting
completed blocks for which the rational certificate kappa<-1/n follows
from a bound on their ordered carry patterns. The rule must be effectively
defined before knowing whether the orbit terminates.

Three separate obligations remain OPEN:

1. Define an ordered-carry boundary quantity that survives nesting, yet
   is not bounded on every fixed finite-feature fibre.
2. Control its creation at resets, including the expanding rational
   shadow families of THM-4507 and the infinite-bank obstruction of
   [THM-4483](../../01-canon/theorems/THM-4483-forced-charges-two-place-ranks.md).
3. Prove a decreasing return is found after finite work for each fixed
   source, without assuming the desired stopping time or replacing that
   source with a residue-population average.

The finite carry graph, curvature average, and positive budget discharge
none of the three alone. The tiling boundary and signed-difference lanes
give explicit examples of what can be lost when gluing or taking an
absolute value, making these requirements testable rather than metaphorical.

## 8. Reproduction and boundaries

Run `python -B 04-computation/experiments/forest_20260926_excursions.py`.
The [script](../../04-computation/experiments/forest_20260926_excursions.py)
checks all parity words through length10, all3249 pairs of positive-carry
word maps of lengths<=5, and510 source orbits through termination or a
4000-step cap. There are21,123 steps and6624 closed nodes in this finite
universe;1606 nodes return higher and5018 lower. Six truncated27-prefixes
also test open-boundary handling. No finite census is a convergence theorem.

The exact eight-tile update is checked at every visited integer, including
the source mode, discarded zero, right boundary and reconstructed output.
The weighted forest identity has105 exact source/horizon/evaluation
controls, at z=1,3/2,2; its polynomial proof is direct composition above.
The independent affine replay checks actual source parity and the unwrapped
integer identity. Each forest node is independently recomposed from its
children and remaining single-step atoms. Every integer's carry stream is
compared with the population-count formula. All checks use explicit
exceptions and stay active with Python -O. Retained
[output](forest_20260926_excursions.out) records the consequence objects.
