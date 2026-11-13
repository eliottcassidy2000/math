# Chan–Pak is a decidable membership predicate and a universal sentence

The question is whether every positive odd integer except 7 and 21 is
the number of Hamiltonian paths of a tournament. Do not confuse a
decision problem over inputs, a truth statement about one explicitly
defined set, and independence from an axiom system.

## A self-contained effective cutoff

In a tournament the strong components form an ordinal chain. A
Hamilton path traverses them contiguously in order, so the count is
the product of their counts. Singleton components contribute 1 and
may be deleted. Every nontrivial strong component has a Hamilton
cycle (Camion's theorem), whose rotations give at least its order
many Hamilton paths. Therefore if h>1 occurs, it occurs on at most
h vertices: for integers a_i>=3 the sum of component sizes is at
most sum a_i, which is at most their product for two or more
components (and the single-component case is immediate). For h=1
use a singleton. This elementary cutoff suffices for every logical
conclusion below.

## A sharper cutoff, and the complexity in bit length

[Busch's strong lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
is H(S)>=5^((n−1)/3) for a strong tournament of order n>=3.
It gives the tidy compressed integer bound
```
3^|T| <= H(T)^3
```
for every nonempty singleton-free tournament T. For a strong
triangle this is equality. For n>=4, 5^(n−1)>=3^n (check 125>=81
then multiply by 5 versus 3). Apply the inequality per component
and multiply. The bound is sharp for the order of arbitrary
compressed witnesses, as an ordinal sum of k triangles has
order 3k and H=3^k. This does not assert it is sharp for the
**minimum** witness size at a given count.

Thus a binary L-bit target needs only O(L) vertices and O(L²)
orientation bits. Exhaustive labelled enumeration with subset
DP gives a decision algorithm in 2^O(L²) time. Exact H of a
guessed tournament is a #P quantity. In the usual complexity
model the decision problem is also in NP with one PP query:
use the sign of the GapP function 1−(H−h)^2 to test equality.
The integer h is itself representable as a #P count by accepting
exactly the binary guesses below it. These upper bounds are
not hardness or efficiency claims.

For a prime target, the nontrivial part of any witness must be a
single strong component. The same is true at a target 7p with p
prime, since 7 itself is excluded. This explains why primality
and 7p inputs are useful to the counterexample search; it does not
prove they occur. The logical analysis does not require factoring.

## The logical consequence

The membership predicate D(h) is decidable by bounded, primitive
recursive operations. With standard arithmetization (or the
equivalent absence of a counterexample search trace), Chan–Pak is
a Pi^0_1 sentence, not merely the superficial forall-h exists-T
Pi^0_2 presentation. There is a fixed machine that enumerates
admissible odd targets, decides each completely, and halts at
the first failure. A false conjecture has a finite verifiable
counterexample: the target and a finite exhaustive nonmembership
computation. In particular PA can prove the negation if the
conjecture is false. Independence from PA, if proved, would
therefore imply that the conjecture is true in the standard integers.

There can be no computable many-one reduction of general halting
to **ordinary tournament count membership**, because that would
decide halting. Reducing one fixed sentence to another machine's
behavior is possible and is useful for comparing proof strategies,
but it does not transfer undecidability of the whole program class
to that one machine. See the [Collatz-style compiler](FRACTRAN_ORBIT.md)
and the separate uniform-family results.
