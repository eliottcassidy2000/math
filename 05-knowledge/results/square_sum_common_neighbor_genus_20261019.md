# The common-neighbor genus jump for prescribed AP petals

**Status: PROVED genus formula; the rational-point finiteness corollary uses Faltings.** This distinguishes many neighbors of one petal from one simultaneous neighbor of many petals.

Let \(r\ge1\), and let \(X_1,\ldots,X_r\) be pairwise distinct rational spreads with \(0<X_i<1/2\). Fix arbitrary signs \(\epsilon_i,\eta_i\in\{\pm1\}\). A candidate new spread \(W\) must satisfy the AP-petal conditions and all prescribed contacts:
\[
1+W=a^2,\quad1-W=b^2,\qquad
1+\epsilon_iX_i+\eta_iW=c_i^2\quad(1\le i\le r).
\tag{1}
\]
The smooth projective curve defined by this multiquadratic function field is geometrically connected and has genus
\[
\boxed{g_r=1+2^r(r-1).}
\tag{2}
\]
In particular \(g_1=1\), \(g_2=5\), and \(g_3=17\).

## All-sign branch proof

The two AP linear forms vanish at \(-1,1\). The i-th contact form vanishes at
\[
\beta_i=-\eta_i(1+\epsilon_iX_i).
\]
Since \(1+\epsilon_iX_i\in(1/2,3/2)\), the sign of \(\beta_i\) is \(-\eta_i\). It cannot equal \(\pm1\), since that would require \(X_i=0\). If \(\beta_i=\beta_j\), their signs first force \(\eta_i=\eta_j\), after which
\(\epsilon_iX_i=\epsilon_jX_j\). Positivity forces \(\epsilon_i=\epsilon_j\), and then \(X_i=X_j\), contradicting the hypothesis. Thus all \(r+2\) finite branch values are distinct for every sign choice.

Over an algebraic closure, the \(r+2\) square classes of these linear forms are independent: each has an odd valuation at its own finite zero and even valuation at every other zero. Adjoining their square roots therefore gives a degree \(d=2^{r+2}\) connected cover of \(\mathbb P^1_W\).

At each finite zero, precisely one generator ramifies, so the inertia group has order two. At infinity all the linear forms have odd valuation, but their ratios are units with square roots in the completed algebraic-closure local field; hence there is again only one order-two inertia group. There is no other ramification. Thus there are \(r+3\) branch values, each contributing \(d/2\) to Riemann–Hurwitz:
\[
2g_r-2=-2d+(r+3)d/2=2^{r+1}(r-1),
\]
which is (2).

## Relation to rational petals and bridges

The conic \(a^2+b^2=2\) is rationally parametrized by
\[
a=\frac{t^2+2t-1}{t^2+1},\qquad
b=\frac{t^2-2t-1}{t^2+1},\qquad
W=\frac{4t(t^2-1)}{(t^2+1)^2}.
\]
Thus (1), after restricting to \(0<W<1/2\) and excluding the finitely many degenerate encodings, is precisely the simultaneous rational AP-petal contact condition. Different square-root sign choices may encode the same petal, which does not affect the genus calculation or finiteness conclusion.

For \(r\ge2\), Faltings's theorem gives only finitely many rational points on the smooth projective curve (1), hence only finitely many simultaneous rational neighbor petals. This remains true if the port signs are allowed to vary: there are only finitely many sign patterns.

For a one-petal insertion between two distinct existing petals, choose the two child-port signs opposite. The resulting bridge curve has genus five. There is therefore no nonconstant morphism from a genus-one curve parametrizing infinitely many such bridges, and only finitely many rational bridges exist for those fixed endpoints.

**Boundary.** Finiteness is not nonexistence. This theorem does not rule out one suitable bridge for every endpoint pair, does not give an effective bound or a method to find every bridge, and does not exclude insertion of two or more petals at once. It does not alter the inherited uniform-count theorem, which constructs arbitrarily many neighbors of one fixed parent and is the \(r=1\), genus-one situation.

It concerns the auxiliary AP-petal architecture. It is not an obstruction to arbitrary square-sum paths or cycles on sparse label sets, which have elementary polynomial-numerical constructions, and it has no bearing on Gerbicz's established ordinary consecutive-interval existence theorem. Prescribed-edge and common-center restrictions are separate questions.

The external finiteness theorem is G. Faltings, *Endlichkeitssätze für abelsche Varietäten über Zahlkörpern*, Invent. Math. **73** (1983), 349–366, Satz 7, p. 365, DOI [10.1007/BF01388432](https://doi.org/10.1007/BF01388432).
