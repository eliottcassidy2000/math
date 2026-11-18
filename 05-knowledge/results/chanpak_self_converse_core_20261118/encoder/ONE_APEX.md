# All strong self-converse tournaments with a transitive deletion

A transitive deletion is a vertex contained in every directed triangle:
the remaining tournament is transitive exactly when it has no directed
triangle. In a strong self-converse tournament with such a vertex, there
are just two descriptions (with allowable overlap at the triangle).

* If there are at least two possible deletion vertices a,b, order all
  other vertices transitively. Since deleting a or b gives a transitive
  tournament, each of a,b has an insertion cut in this order. Say a->b.
  If the cut of a is no later than that of b, the entire tournament is
  transitive. Otherwise the entire interval between the two cuts forms
  the spine of a nearly transitive tournament: put b, the interval,
  and a in increasing order and reverse only a,b. Everything before
  the earlier cut dominates this central interval and everything after
  the later cut is dominated by it. Strongness precludes either
  outside tail. Thus the strong tournament is exactly the transitive
  order with its extreme arc reversed, with H=2^(n-2)+1.
* If there is exactly one possible deletion vertex v, every reversing
  isomorphism fixes v. Its action on the remaining transitive tournament
  must reverse its unique order. This forces the length-2k out-tag word
  to equal its reverse complement. Strongness forces first tag 1 and
  last tag 0 (otherwise the first or last spine vertex is a source or
  sink in the full tournament). Hence the word is precisely u followed
  by the reverse-complement of u, for a positive binary integer u
  of width k. Conversely these tags make a directed spanning cycle
  and exhibit the desired anti-involution.

Using the [one-apex address law](https://github.com/eliottcassidy2000/math/blob/395537b8edc50a89a42bb7742e61c20828620ab3/05-knowledge/results/chanpak_exponent_dimension_20261117/apex/PROOF.md),
the second alternative reads

    H = u²+2F(u)−1+[u odd].

Thus the whole strong single-deletion symmetric sector has only
Theta(sqrt X) distinct counts through X. For the upper bound H>u²
and the extreme-flip family contributes only O(log X).
For the lower bound, induction in the address recursion gives
F(u)<=u+1; if u is odd, F(u)=F((u−1)/2)<=(u+1)/2. Hence each
diagonal count lies in (u²,(u+1)²], and is distinct. For example
u=1,2,3,4,5 gives 3,9,11,25,31, respectively. This is a structural
reason that one deletion interface, however efficiently computed,
cannot cover the odd counts by itself. A second deletion vertex
isn't a second independent arithmetic knob here; it actually puts
the strong tournament back in the extreme-flip family.
