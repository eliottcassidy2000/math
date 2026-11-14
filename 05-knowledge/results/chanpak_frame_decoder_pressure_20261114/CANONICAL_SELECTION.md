# What a natural representative can and cannot prove

A first answer is **minimum order, then least adjacency code after
canonical labeling**. That really is one isomorphism class for each
attained value, and there is a finite exact decision procedure for it.
It is not, on its own, a proof that the domain is every odd value
except 7 and 21.

## Global normalized selector

Every tournament is the ordinal sum of its strong components in a
unique order; the counts multiply. Every minimum-order witness for
(N>1) has no singleton components, since singleton deletion
preserves the count. Each surviving component has count (h_i\ge3).
[Busch's strong lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
gives (h_i\ge C5^{n_i/3}) with an absolute constant (C>0).
Since (prod h_i=N) and the number of components is at most
(log_3 N), this gives (O(\log N)) vertices for a
minimum-order witness of any attained (N>1).
The case (N=1) uses the singleton.

One can therefore define (\mu(h)) as the least strong order
for count (h\ge3), with infinity for failure, and compute
[
 m(N)=\min\{\mu(N),\min_{1<d<N,\ d\mid N}(\mu(d)+m(N/d))\}.
]
The equality follows by peeling off the first strong component,
and conversely by ordinal joining the minimizing witnesses.
Compare all minimum-cost concatenations by the chosen canonical
adjacency code (one can use a tuple of the canonical codes of the
ordered components; their order is isomorphism invariant). This
is an exact global representative rule **on the attained domain**,
including actual ordinary isomorphism rather than just equivalence
under converse. The analytic bound makes its search finite, not
efficient in the bit length of (N). The first-terrace decoder
replaces the strong-component search with a structural inverse
only inside its restricted domain.

The inherited prime and 7-prime reduction from Sessions 323–324
offers a different representative rule that intentionally does not
minimize order: fix one canonical strong witness for each primitive
target, and use arithmetic factorization to assemble the rest.
For (v_7(N)=0) use prime factors in increasing order.
For (v_7(N)\ge2), write the 7-power as a product of fixed
49 and 343 witnesses, using one 343 iff the exponent is odd,
and then append the other prime factors.
For (v_7(N)=1), pair the 7 with the smallest prime divisor
(p\ge5) of the cofactor; if no such prime exists, use 63 and
append the leftover factors of 3. The last case fails only for
7 and 21. This rule becomes unconditional as soon as one proves
strong witnesses for all odd primes except 7 and all
(7p) for primes (p\ge5, p\ne7). It has no branching ambiguity
and strictly decreases the remaining target; the difficult step is
producing the strong witnesses, not choosing how to assemble them.

## Exact obstruction to choosing a smaller scalar parent

The first occurrence of 63 is the two converse classes at order
eight, by the [complete count-63 census](https://github.com/eliottcassidy2000/math/blob/1b91742ce/01-canon/theorems/THM-344-n8-h63-complete-omega.md)
(and the earlier order census). Busch's bound is 75 at order nine,
and all larger strong orders have larger bounds. A decomposable
witness with two nontrivial factors would require one factor
divisible by 7, necessarily 7 or 21; both are impossible.
So any count-63 witness is one of those two cores with
optional singleton padding.

Use the ascending spine plus hubs (x\to y), with word
`100110` specifying the two arcs from the hubs to each spine site.
The eight deletion counts are
[
 (11,37,37,33,33,25,25,1).
]
Taking converse preserves this deletion-count multiset.
The nontrivial values occur by order six, from the two-hub
words `110`, `1200`, `1120`, and `1121` respectively
(assigning 37 to `1200`, 33 to `1120`, 25 to `1121`);
count 1 occurs at order one. The literal replay checks every
assertion about these path counts.

Thus *no* minimum-order witness of 63 can be reached by a
count-increasing single insertion from a minimum-order witness
for its smaller scalar count. Padding does not help: remove
a padding singleton and the count stays 63; remove a core
vertex and the redundant order-seven environment remains.
This failure is independent of the canonical tie-break.

A prospective global selector should choose after keeping several
possible endpoint states, not force its source to be the scalar
representative. The clean induction obligation is:
for every unresolved primitive target and some permitted endpoint
state, exhibit an exact inverse move to a smaller target in a
state whose availability is guaranteed. Then a canonical rule
chooses the least verified branch. Preserving the state is what
makes the recursion close. The witness selection machinery here
does not prove that universal inverse move.
