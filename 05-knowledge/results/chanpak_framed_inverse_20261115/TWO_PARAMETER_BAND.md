# A second scale of same-graph adjacent pairs

The coordinated placement in the
[complement-star inverse](COMPLEMENT_STAR_INVERSE.md) survives a
second independently changing length. This is stronger than padding
the final tournaments by vertices that leave the count unchanged.

Let e encode the nonterminal mask of the rooted remainder as before.
Before attaching the single outlet, append t>=0 new unmarked
vertices to its transitive spine, immediately before the root.
Thus each original marked position is t positions farther from
the final root. The resulting even path response is E_t=2^t e.
The ordinary count of the remainder is 1+p(2^(t+1)e)
by the same ordered-word formula; the recurrence for p gives
the physical target
```
N_t(e)=1+2^(t+2)e−J(e).
```
Here J is the pair energy on the original set bits and does not
change just because a zero was appended at the low end.
The continuation shape is still a complete graph minus a star.

Put A=2^(t+2), and separate lower and higher scales by taking
e=x+2^k y with k>=4. The exact response is
```
N_t(e)=2^k y(A−rho(x)) + 1 + A*x − J(x) − J(y).
```
Choose
```
(x,y) = (6,8A−5)       or       (5,8A−3).
```
Their slopes coincide: 8(A−rho(6))=8A−3 and
8(A−rho(5))=8A−5. They also have equal bit lengths and equal
popcount. To calculate their constants, J(6)=J(5)=2 and
```
J(8A−5)=7A−6,     J(8A−3)=6A−4.
```
For these last identities start with A=4: the energies are
22 and 20. Doubling A inserts a new highest bit of weight 8A
in each y, with other bits unchanged. The interaction added
is the old bit reversal plus one for the adjacent old
highest bit. Since the y words have all bits set except
position 2 or 1, their reversals are respectively
7A−1 and 6A−1. This gives the recurrence and proves both
formulas at every t.

Thus for every k>=4,t>=0 the **literal strong adjacent pair**
in a single unlabelled graph shape is
```
(8A−5)(8A−3) 2^(k−3) + 3−A,
(8A−5)(8A−3) 2^(k−3) + 5−A.
```
The full tournaments have order k+2t+8. Their complement
charts are the complete graph on k+2t+7 gaps with t+6
edges removed from one star; the distinguished gap has
degree k+t, so the ordinary undirected Hamilton count
is positive and the same on both charts. The same
isomorphism assertion holds before the final outlet.
At t=0 this recovers the first adjacent pair 1565,1567.
The count is not unchanged by varying t.

A target-first decoder for these bands tries powers A=2^(t+2)
up to sqrt(N): the expression bounds A² by N. For each A,
divide N+A−3 or N+A−5 by (8A−5)(8A−3) and require an
exact power of two of exponent at least one. Returning the
corresponding x,y,k gives the physical directions and
the actual rooted remainder, without enumerating graphs.
This is an inverse on a larger explicit bank, still at most
O(log² X) distinct outputs through X. It does not fill the
gaps between bands.
