# Even all marked ratios together cannot fill absolute counts with a bounded two-port boundary

This obstruction goes beyond the single 24:13 grammar. It allows
the common marked ratio to be chosen separately for every target,
with no bound on its denominator. What remains fixed is the
transitive two-apex spine, the use of the same line for both
address vectors, and a bounded boundary supporting two disjoint
odd cycles.

## Count the zero determinants without fixing a ratio

Let a,b be the forward binary addresses of the two rows of an
n-site tape. Write C=2^n−1 and r(t)=C−rev_n(t), the address of its
reversed complement. The same-line condition is

    a r(b) = b r(a).

For 0<a,b<C this is the assertion that
F_n(a)=F_n(b), where F_n(t)=r(t)/t. At the endpoints the only
additional solutions are (0,0) and (C,C).

In each nonempty fiber F_n(t)=p/q>0, with p,q coprime, there
are at most 2^ceil(n/2) addresses. It suffices to prove this
when p≤q: the involution r takes the fiber p/q bijectively to q/p.
Let m=floor(n/2), and write t=2^m x+y, 0≤y<2^m, using the
fixed widths for reversal. For fixed x the equation

    p t = q(C−rev_n(t))

requires that p y + q·2^(n−m) rev_m(y) take a particular
constant value. This is injective in y: a nonzero difference
in the reversed term has magnitude at least q·2^(n−m),
larger than p(2^m−1). There are only 2^ceil(n/2) high parts x.
Thus summing squared fiber sizes, or choosing a and counting
its possible partners b, bounds the number of zero-determinant
tapes across **all ratios** by O(2^(3n/2)).

## Put the addresses on the Hamiltonian-count scale

Fix a constant k. Suppose both apices beat distinct sites among
the first k, and are beaten by distinct sites among the last k.
Narrow boundary patterns such as the first two and last two
digits of the 24:13 language satisfy this with k=2. For n≥2k
choose these two entrance and two exit sites, and assign each
of the n−2k other interior sites to the A cycle, to the B cycle,
or to neither. Retain assignments with each chosen subset even.
There are (3^(n−2k)+2+(-1)^(n−2k))/4 assignments, by the two
parity filters. The assignments yield distinct positive packs
of two directed odd cycles, contributing four each to H. Hence

    H ≥ 3^(n−2k).

The short cases n<2k contribute only a fixed finite number
of counts. Any counted value H≤X therefore has
n≤log_3 X+2k. Combine the all-ratio tape bound with summing
this geometric series and allowing both mutual orientations.
The number of distinct count values at most X is at most

    O_k(X^((3/2) log_3 2)),

with exponent <1 because 2^3=8<9=3^2. In particular not
even an infinite palette of target-chosen ratios can cover
a cofinite arithmetic lane in this bounded-boundary model.

This is not a bound on all two-apex tournaments. The determinant
need not vanish if the interface can read and compensate its
actual nonzero value. The boundary may also have no room for
the two distinct cycles or may grow with the target. A more
promising descent wheel should program a changing signed
determinant, or use a different notion of compatibility altogether,
rather than just add more zero-determinant rays.

## A finite list of nonzero determinants is not enough either

The same exponent holds for each specified signed determinant D,
uniformly in D and n. For a fixed a with r(a)≤a and a>0, the
equation is

    r(a)b + a rev_n(b) = aC − D.

For fixed high half of b the reversed term again dominates
the linear low-half term, so there is at most one choice of its
low half. The case a=0 fixes b directly if n≥1. If r(a)>a,
apply reversed complement simultaneously to both addresses.
This interchanges a and r(a) and negates D, giving the same
bound on that part by injection. Thus at most O(2^(3n/2))
tapes can realize any prescribed D.

For a menu of K determinant values, the upper bound at count
scale X, under the same two-port boundary, becomes
O_k(K X^((3/2) log_3 2)). In particular a fixed or polylogarithmic
menu is still sparse. To cover Ω(X) distinct values through X
requires a menu of size Ω(X^(1−(3/2)log_3 2)) in this model,
even if selecting a determinant and an address costs nothing.
This is a necessary capacity, not a sufficient construction.
