# One locally complete five-clock is far sparser than the full terrace

Let \(\mathcal V\) be the set of distinct ordinary path counts in
the family whose corridors are repetitions of \(4,2,5,3,7\) with
arbitrary positive widths (as in [the local proof](ODD_MONODROMY.md)).
Then
\[
 \sum_{N\in\mathcal V}N^{-203/500}<\infty,\qquad
 \#\{N\in\mathcal V:N\le X\}=o(X^{0.406}).
\]
This is not an upper bound for the full terrace, for which the
different exponent 0.952 was proved, or for arbitrary tournaments.
The ordered-word Dirichlet sum already diverges at 0.401, which
is **not** a lower bound on distinct-value growth.

In column state \(v=(a,b)\) divide by its readout \(a+b\).
Set \(p=b/(a+b)\in[2/3,1]\), which is invariant because \(b'/a'\ge2\).
At seam \(q\) with width \(r\ge1\), the physical readout is
multiplied by
\[
 G_q(p,r)=e_q+(f_q-e_q)p + 2^r\bigl(f_q+(g_q-f_q)p\bigr).
\]
It is increasing in \(p\). Thus it is bounded below by
\(G_q(2/3,r)\) and above by \(G_q(1,r)\), uniformly in
length and order. With
\(s=203/500\), the replay gives the strict inequality
\[
 \prod_{q\in(4,2,5,3,7)}
       \sum_{r\ge1}G_q(2/3,r)^{-s}<0.99.
\]
Each series is bounded using exact rational ceilings for the
first 32 terms and a geometric tail, since
\(G_q(2/3,r)\ge (f_q+2g_q)2^r/3\). This factors the cost of one
clock period. The initial-width sum
\(\sum_{r_0\ge1}(1+2^{r_0})^{-s}\) is finite.
Summing over all periods gives a convergent geometric series,
even when ordered word repetitions and all coincidences are counted
with multiplicity. The little-\(o\) count follows by splitting the
sum at an arbitrary fixed \(N_0\): the finite initial counts
are \(o(X^s)\), and the rest contribute at most \(X^s\) times
a tail tending to zero with \(N_0\).

Conversely, at \(s=401/1000\), the replay certifies
\[
 \prod_{q\in(4,2,5,3,7)}
       \sum_{r=1}^{32}G_q(1,r)^{-s}>1.01 .
\]
Using only these widths and any fixed initial width gives exponential
growth of the word sums at this exponent. This is a pressure bracket
for descriptions, not for distinct integers, and does not use
a heuristic model of collisions.

The five active lengths sum to 21. This is an actual cost in the
word: every clock adds at least \(21+5(1+1)=31\) vertices, even
before making large congruence representatives. The arithmetic
coincidence does not by itself prove any additional forbidden count;
the growth theorem instead uses the exact readout gain phase by phase.

## Contrast: a recoverable sublanguage and two-source obstruction

In the same family there are \(\Omega(X^{49/250})\) distinct
physical totals. Restrict widths to \(1\le r\le31\), except for
\(q=3\), after which require \(2\le r\le32\).
The initial state is fixed at \((1,2)\).
The entry \(a\) stays odd; \(b\) stays even.
For each active \(q\) the outgoing \(b'=2^r(f_q a+g_q b)\)
satisfies
\[
 v_2(b')=r+v_2(f_q),\qquad
 (v_2(f_4),v_2(f_2),v_2(f_5),v_2(f_3),v_2(f_7))=(0,0,2,1,2).
\]
For \(q=4,2\) this is odd plus even inside the bracket.
At \(q=5\), \(g_5=24\) and \(b\) is even, so its term has
larger valuation than the \(f_5a\) term.
At \(q=3\), \(b\) has valuation at least 3 from the preceding
\(q=5\). At \(q=7\), the preceding \(q=3\) leaves valuation
at least 3 by the width restriction.
Thus the outgoing two-vector determines its last width through
the indicated valuation. Undo this width and invert the
nonsingular rational kernel to recover the previous state.
This is an injective code for words made of full cycles; the
phase on reversing is known and the seed is the only state in
the code of scalar 3, since every step grows strictly.

At \(\sigma=49/125\), sum the weights \(c^{-\sigma}\) over all
allowed tuples of widths of one complete cycle, where
\(c=\prod_{q}G_q(1,r_q)\) is an upper bound on its gain.
The exact replay certifies this sum exceeds 1.02.
Draw independently a cycle at each step with probability
proportional to this weight. Stop when the product of the cycle
costs first exceeds \(Y\). It is between \(Y\) and \(C Y\), for
a fixed maximum allowed cost \(C\), and each stopped word has
probability at most \(Y^{-\sigma}\). This finite probability-one
set therefore contains at least \(Y^\sigma\) different words,
and by the inverse code they give different two-vectors of
size \(O(Y)\).

To expose the entire vector in an ordinary count, append one
further complete physical cycle with its first four widths
equal to 1 and its final width-power equal to a common
power of two \(P\) of size \(O(Y)\).
Before multiplication by \(D_P\), its two coordinates
\(A(a,b),B(a,b)\) are an invertible linear image over
\(\mathbb Q\) of the old state, both positive and \(O(Y)\).
Choose \(P\) greater than every \(A(a,b)\) in the stopped set.
Then the literal scalar \(A(a,b)+P B(a,b)\) encodes both
coordinates by quotient and remainder, and is \(O(Y^2)\).
Taking \(Y\) proportional to \(\sqrt X\) gives the lower
bound. The closing width can exceed the alphabet for the
injection proof, but remains a legal clock width.

Finally, suppose a proposed absolute-target compiler selects
two clock counts \(x,y\), allows only a fixed finite menu of
output choices per pair, and has output \(N\ge c\max(x,y)\)
for some fixed \(c>0\). Then its distinct-value inverse-power
series converges at \(2s=0.812\), since
\[
 \max(x,y)^{-2s}\le x^{-s}y^{-s}
\]
and the clock inverse-power series at \(s=0.406\) converges.
A finite number of such compilers and subsequent arbitrary
products of their nonunit outputs still give convergence at
0.812 by the multiset Euler product. This excludes a universal
two-reservoir program of this specified kind. The cost condition
holds in particular when both selected source tournaments remain
induced subtournaments, by insertion monotonicity; the
finite-menu condition is an independent restriction. A signed
correction inspired by MRB escapes the **numerical** inequality
only if the operation can actually undermine that retained-source
condition or retain more free parameters. The theorem makes
no assertion that three sources suffice.
