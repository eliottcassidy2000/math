# Bounded marked repairs have an exact dyadic inverse

Fix a nonempty source A, write h_A=H(A) and q_A=P_A(2) for
the count of two ordered optional paths. For a nonempty host B
with marked set O, write h_B=H(B), M_B for ordered optional
covers whose second path has marked start, and W_B for those
with both marked first-path end and marked second-path start.
The physical closer (B ⇒ A, then add an apex beating precisely
O) has
~~~text
 N=(h_B-W_B)h_A + M_B q_A.
~~~
This is the marked-closer identity from the pair-line session.
Here W is the actual correction against the baseline
(h_B h_A + M_B q_A), not an arbitrary difference of outputs.

## 1. Every pair of marks imposes a positive repair

For each unordered pair of marks x,y, let W_{xy} count covers
using x,y as the tested endpoints in either order. Contract
x,y to a vertex z in a new tournament on one fewer vertex.
Keep all arcs not incident with z, and let z lose precisely
to the common in-neighbors of x,y. In a Hamilton path of the
contracted tournament, split at z. If the path has a successor
after z, at least one of x,y beats that successor. Use one
such vertex (choose by fixed order if both) for the start of
the second path and the other for the end of the first. The
predecessor, if present, beats both. If there is no successor,
either choice is valid. Both paths in the resulting cover are
nonempty, and contracting the two tested endpoints back to z
recovers the original path. Hence (W_{xy} ≥ H(B/xy) ≥ 1).

Apply the closer to a singleton source with exactly x,y
marked. Its output is the odd number (h_B-W_{xy}+2M),
so (W_{xy}) is even. Consequently, for (k=|O|),
~~~text
 W_B = sum over pairs {x,y} in O of W_{xy} ≥ k(k-1).
~~~
In particular W=0 iff there are at most one mark. This rules
out the idea of choosing arbitrary slopes by adding many
mark sites while keeping a zero repair. The contraction
argument also gives a graph-dependent bound, even when both
marks lie in one strongly connected component.

## 2. Width budget from an SCC cut

Write (B=P⇒C⇒S), where C spans all the SCCs from the
first through the last marked component, inclusive. P,S are
unmarked; when empty, their h and q are 1. Then exactly
~~~text
 h_B=h_P h_C h_S, M_B=h_P M_C q_S,
 W_B=h_P W_C h_S .
~~~
These follow by splitting a path at the SCC cut: a marked
start in C requires all vertices of P to lie in slot 1; both
marked endpoints in C additionally force all of S into slot 2.

If the first and last marked components F,L are distinct,
write (C=F⇒U⇒L). Pick a 2-slot cover of F whose second
slot starts marked, and one of L whose first slot ends
marked. Such covers exist: a marked vertex alone in the
specified slot and any spanning path on the other vertices
suffices. Put any optional cover of U in between and concatenate
by slots. All these covers survive, so (W_C ≥ q_U ≥ 2^{|U|}).
No bound on the individual SCCs in U is needed here.

Suppose the actual coefficient (W_B ≤ R), there are at
least two marks, and the two boundary marked SCCs each have
at most d vertices (including the case of a shared SCC).
Then
~~~text
 |C| ≤ 2d+floor(log_2 R),
 h_P h_S W_C ≤ R, k(k-1) ≤ R.
~~~
Thus bounded repair forces a bounded interval of marked SCCs;
it does not merely constrain the number of marks.

## 3. Target inverse, with literal reconstruction

Under those bounds and for the fixed source A, the possible
positive targets are a finite union of explicit rays
~~~text
 u+v2^t, t ≥ 0,
~~~
where (with (S_0) the suffix with singleton SCCs stripped)
~~~text
 u=h_P(h_C-W_C)h_{S_0}h_A,
 v=h_P M_C q_{S_0}q_A>0.
~~~
To decide and reconstruct for a prescribed N, enumerate marked
C from the width bound, whose endpoints contain marks, and
reduced prefixes (P_0) and suffixes (S_0) subject to the
product budget. Check whether ((N-u)/v) is a positive
**integral power of two**, and if so append that many binary
clock bits, i.e. t singleton SCCs, to the suffix. Singleton
prefix sites have no effect. There is no suffix bound that
needs to depend on the target. Each enumerated ray gives a
literal witness, not only a modular condition. To see that
the reduced tails are finite, a nontrivial strong tournament
of order n has a spanning directed cycle (Camion), hence at
least n Hamilton paths by cutting that cycle. All those SCCs
have n≥3. The product of their h counts is at most R, so there
are finitely many possible components and finitely many
possible lists. Their total order is at most that product.
For completeness, take a longest directed cycle in a strong
tournament. Any vertex outside it that both beats and loses to
cycle vertices can be inserted at a consecutive change of
orientation. Otherwise the outside vertices partition into
those beaten by the whole cycle and those beating the whole
cycle. Strong connectivity forces an arc from the first part
to the second; its two endpoints can be inserted together into
the cycle. This contradiction proves the claimed spanning
cycle fact without a separate citation.

The stronger Moon–Busch lower bound quoted in the earlier
pair-line session shrinks the enumeration to O(log R) vertices.

If N is prime, (h_P=1), because it divides N and the divided
bracket is bigger than one: (W_C ≤ M_C), (q_{S_0} ≥ h_{S_0}),
and (q_A ≥ 2h_A). For more general N only divisors of N
need to be tried for h_P. Any fixed finite source bank and fixed
R,d give O(log X) outputs up to X. Thus a strategy with two
or more marks, bounded absolute repair, a fixed source bank,
and bounded marked SCC boundaries cannot place a baseline at
every target. To defeat the bound one must grow at least one of
those budgets, use the one-mark case, or change the construction.

### Completely explicit two-mark inversion for transitive hosts

A simple nontrivial infinite example is two marks at positions
whose distance is d≥1, followed by t transitive singleton sites;
leading unmarked sites do not matter. Then
(M=(1+2^d)2^t), (W=2^{d-1}+[d=1]), and (h_B=1).
The exact formula is
~~~text
 N=(1-2^{d-1}-[d=1])h_A +(1+2^d)q_A2^t.
~~~
With the singleton source, d=2, t=3 gives exactly N=79.
For any prescribed repair cap W≤R, enumerate at most
(1+floor(log_2 R)) distances and apply the quotient
test. As in the general theorem, matching a residue or an
approximate dyadic ratio is insufficient.

## 4. An all-order limitation of transitive bounded repair

In a transitive host with at least two marks, the extremal
pair alone contributes (2^{d-1}+[d=1]) for distance d.
After stripping unmarked leading and trailing vertices,
the host is therefore a binary word of length at most
(2+floor(log_2 R)), with first and last bits 1.
For a fixed source A this gives **O(R log(2+X))** possible
outputs at most X with actual repair W≤R, even as R grows
with X (the constant depends on A). Thus allowing a
sublinear repair (Rle X^{1-epsilon}) cannot make a
fixed transitive-host source bank cover a positive fraction
of targets. The claim is not about small differences from
an independently chosen baseline, only the physical W
coefficient in the closer formula.

A target-first search should therefore vary the source profile
and the marked core together instead of first insisting on
arbitrarily small W. The exact ray test provides a way to close
the bounded branch rigorously before escalating either budget.
