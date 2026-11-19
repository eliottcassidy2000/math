# A directed p-gon cannot be completed neutrally; use a compressed core

**Proved all orders:** the directed odd polygon already forces a new
Hamilton path once just its distance-two chords are oriented. **Exact
constructions:** five selected large prime counts, including a 46-vertex
witness with 10,000,000,000,037 paths. No proof of prime-count infinitude
or prescribed-prime coverage is claimed.

## The distance-two obstruction

Let n>=5 be odd, and orient the polygon 0→1→...→n−1→0. It has exactly
n Hamilton paths before adding chords. For every i modulo n give its
distance-two chord a sign: + for i→i+2 and − for i+2→i.
These chord pairs are distinct for n>=5.

* If all signs agree, the step +2 (or −2) forms another directed
  Hamilton cycle, because gcd(2,n)=1. Breaking it gives n additional
  paths, using none of the original cycle edges.
* Otherwise some consecutive signs are +,−, say at i,i+1. Then
  the following is a new spanning path using only specified arcs:

      i+4, i+5, ..., i (clockwise), i+2, i+3, i+1.

  The first chord is i→i+2; the final chord is i+3→i+1.
  Its vertex list contains every vertex once. Different +,−
  transitions give different terminal vertices, hence different paths.

Thus even these distance-two arcs cannot be filled without creating
extra paths, regardless of the other chords. If there are t +,−
transitions in the mixed case, a full tournament has at least
n+2 ceil(t/2) paths by Rédei oddness. In the uniform case it has at
least 2n+1. See Busch for a classical source of the oddness and
insertion facts: https://doi.org/10.37236/1141 .

Even adding more vertices cannot rescue a p-gon: every path on its
induced subtournament extends by deterministic single-vertex
insertion, injectively since deletion recovers the old order.
This does not rule out a differently structured partial tournament
having p paths and some neutral completion; it rules out a
directed p-cycle as its source of p choices.

## Prime counts require broken rotational symmetry

A prime-count witness has just one nontrivial strong component
(counts multiply across components, and components of count one
are singletons). For prime p>=5 this component has at most
1+3 log_5(p)<p vertices, by Busch's strong-order bound. A tournament
automorphism acts freely on ordered Hamilton paths; its group order
thus divides p. If it were nontrivial it would have order p and
would act faithfully on the strong component, since every singleton
component is individually fixed by their linear order. A nontrivial
permutation of order p needs at least p vertices, a contradiction.
So prime-count tournaments are rigid for p>=5. Cyclic indexing
of the paths is fine; a physical rotation symmetry of period p
is not.

## Compressed digit witnesses

Use the inherited one-apex join F1 from prime_witness_frontier/
targeted/PROOF.md. For positive addresses b,e, encode their binary
integers b−1,e−1 as two apex words on transitive spines. Place the
converse of the second block on the right; all cross arcs point
right to left except a bridge from the left apex to the right apex.
Let F(u) be the Hamilton count of the one-apex block for address u,
with F(0)=1, F(2z+1)=F(z), and F(2z)=F(z)+z+[z odd] for 2z>0.

The exact F1 count is be+F(b−1)F(e−1)−[b,e both odd].
The first term counts the independent ordered covers when the
bridge is used; the other terms count the paths avoiding it
(two whole blocks minus the forbidden apex-to-apex seam).
There is no need to assume extra polygon chords are neutral.

| Prescribed prime p | b | e | Order |
|---:|---:|---:|---:|
| 1,000,000,007 | 3 | 191311747 | 32 |
| 1,000,000,009 | 7 | 86246185 | 32 |
| 1,000,000,000,039 | 1613 | 583862840 | 43 |
| 1,000,000,000,061 | 7 | 88583919141 | 42 |
| 10,000,000,000,037 | 19 | 404043168698 | 46 |

These are finite exact constructions, not minimum-order claims.
The supplied code prints the outgoing row masks, verifies primality
by trial division, rechecks F independently by summing over possible
neighbors of the apex, checks 256 small joins by literal subset DP,
and tests every chord sign word at n=5,7,9,11. Large counts use the
proved separator formula, not subset enumeration on 46 vertices.

## What would prove an infinite prime supply?

One full reduced arithmetic progression of realized counts would
suffice by Dirichlet. Alternatively, a proof of o(X/log X) missed
odd counts below X would force all but o(X/log X) primes to occur
by the prime number theorem. A power-saving miss bound would do it.
But relative density one among odds, full congruence coverage,
multiplicative closure and any finite verified prefix do not suffice:
take the odd composites except 21, together with 1 and the finite
verified odd primes other than 7. That set has all those properties
and no further prime values. It is a logical countermodel, not a
candidate tournament spectrum.

The productive replacement for the p-gon is therefore a compressed
binary-address construction with an exact range theorem, or a
sufficiently sharp estimate on how many targets its joins miss.
The example search certifies some prescribed primes, not a
universal or unbounded prime selector. Chan–Pak's original statement:
https://arxiv.org/abs/2308.10214 .

Run python3 polygon_prime_replay.py.
