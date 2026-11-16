# The 67/69 seam: a real contrast, its selector, and its absolute limit

The [verified 67/69 packet](https://github.com/eliottcassidy2000/math/blob/9fb47c1a516353dbce1154a1d59d923841ac73a1/05-knowledge/results/chanpak_pair_67_69_20261114/README.md)
uses outgoing rows

    A = (2,12,17,37,75,23,47),  B = (2,12,17,37,43,71,31).

These give H=67,69, and have the same labelled outdegree tuple
(1,2,2,3,4,4,5). Relabel B by the two swaps (2 3)(4 5): the new B
differs from A by only the arc between 2 and 3. The numerical contrast
is therefore a real edge move, not a projected pairing of unattached
numbers. It does not say that the prime 67 or the factorization
69=3·23 generated the move.

## One precise finite Vitali analogy

The [Vitali recovery](../vitali_atom_recovery_20261011.md) and its
[selector certificate](../vitali_selector_certificate_20261005.md)
separate two unrelated issues: selecting one point of a quotient fibre,
and computing a function that actually descends to that quotient.
The test is exact: the function descends iff it is constant on each
fibre. The older equal-lambda witnesses 109/111 fail that test even
for labelled triangle-pair counts. The 67/69 witnesses fail it for
the labelled score tuple, not for the finer triangle-pair observation.

The sharper 67/69 analogue uses the two independently inflatable
marked slots 0,1. Their projected path-word tensors agree at every
entry except the (1,1) entry, which rises by 2. Thus for **every**
pair of transitive slot sizes the eight higher-visit entries, even
if observed exactly, do not determine which of two adjacent odd
counts occurs. If T,U occupy the two slots, the missing quantity
is exactly **2H(T)H(U)**, not always 2. A selector that wishes to
recover the exact count must keep the arm bit and the relevant
module factors. This is a finite marked-fibre statement; the
classical Vitali nonmeasurability paradox is not invoked.

Classical Vitali interval selection makes a different point. Suppose
a finite collection of physical intervals each certifies every odd
integer it contains. Greedily choose a longest remaining interval
and discard every interval that intersects it. Each discarded
interval lies inside the concentric triple of the picked interval:
its radius is at most the chosen radius, and it meets the chosen
interval. Hence the selected disjoint physical intervals carry at
least one third of the union's total real length. For infinite
families with bounded windows one can first take finite subfamilies.
This gives a packing/density bound, not pointwise coverage.
The proof assumes the **whole intervals** are physically covered,
not merely that the endpoints of the intervals are counts.

No density conclusion alone gives cofiniteness: even the set of all
positive odds except the odd numbers 2^j+1 (j≥2) has relative density
one and arbitrarily large adjacent pairs in every residue class
modulo every fixed even modulus. The omitted sequence is still
infinite. This is the precise limit of a Vitali-style step.

## A signed contrast atlas with a closed-form inverse

The [old third-module law](https://github.com/eliottcassidy2000/math/blob/9fb47c1a516353dbce1154a1d59d923841ac73a1/05-knowledge/results/chanpak_pair_67_69_20261114/THIRD_MODULE.md)
lets nonempty tournaments X,Y,Z replace old slots 0,1,6. Let
h_T=H(T), and c_T be the number of ordered two-slot path covers with
both slots nonempty. Its signed path-count difference is

    Δ = H(B')-H(A') = 2 h_X h_Y h_Z - (2 h_X+c_X)(h_Y+c_Y)c_Z.

For transitive slots of sizes n,m,r (all positive), normalise d=Δ/2.
Since h=1 and c=2^s-2 for a slot of size s, the exact readout is

    d = 1 - 2^n (2^m-1)(2^(r-1)-1).

At r=1 this reproduces +2; at r=2,n=m=1 it gives -2,
the physical pair 135/133. The new value is 133=7·19, illustrating
a 7p lane without supplying all such lanes. For r≥2 the contrast
is negative. Thus growing three physical modules can store large
signed displacements, but its register is sparse.

To decide whether a given negative d belongs to this atlas, compute
n=v2(1-d) (must be at least 1) and odd Q=(1-d)/2^n.
If Q+1 is a power of 2, the unordered Mersenne exponents are
{1,log2(Q+1)}. Otherwise any two exponents x,y must both be at
least two, and

    2^(x+y-1) < (2^x-1)(2^y-1) < 2^(x+y).

Hence L=ceil(log2 Q) must be x+y, and 2^L+1-Q must be 2^x+2^y.
Its binary bits determine {x,y}, also when x=y and the sum is a
single power of 2. Check the exponent sum and each lower bound.
The two cases cannot collide: with both exponents at least two a
Mersenne product Q=2^L-1 would force 2^x+2^y=2. This is an exact
inverse up to interchange of m=x and r-1=y, not a density argument.

It follows that the number of distinct negative half-gaps of
magnitude at most X is

    (log_2 X)^3/12 + O((log_2 X)^2).

Indeed log2(1-d) lies between n+x+y−2 and n+x+y, so the logarithmic
boundary differs by only O(1) layers. The unordered triples
(n,{x,y}) with positive coordinates and sum at most L number
L^3/12+O(L^2); the diagonal x=y accounts for at most O(L^2).
The formula explains why arbitrarily growing dilation is weaker
than control of a growing interval of **all** signed contrasts.

## What kind of conjectural implication survives these tests?

A genuine implication from Collatz would need a tournament lift
for the **specified predecessor** (2^a q−1)/3, with q's witness
converted to a compatible marked state; arithmetic termination
would then be a valid induction rank once paths through excluded
counts have literal bypasses. A genuine Goldbach implication would
need a physical inverse for the **specific selected prime pair**,
not merely the existence of two prime labels in a modulus class.
The 67/69 tensor, even with three independently growing slots,
supplies neither lift. These are missing structural lemmas, not
consequences of the named arithmetic conjectures alone. Likewise
the small prime 67 being adjacent to 69 does not establish any
prime-pair distribution statement.
