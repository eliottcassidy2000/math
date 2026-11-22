# A path-conflict grammar for the full scalar and one-arc isomorphism series

**PROVED.** There is an exact recursion in objects other than tournaments:
partial assignments forced by entire directed Hamilton paths. A weighted
independence-set transform gives the polynomial that counts isomorphism
classes by ordinary Hamiltonian-path count at each size. Coloring the path
atoms by which side of a distinguished flip they belong to gives the
joint polynomial for literal one-arc pairs. The two polynomials provide
different minimum-order data.

## 1. Work in path atoms, not graph classes

Fix \(n\) named vertices and a permutation \(g\). If \(g\) has an even
cycle it fixes no tournament, since a half-rotation of that cycle swaps
a pair of vertices and would reverse its arc. Otherwise \(g\) has odd
order \(r\). On each orbit of unordered vertex pairs, an invariant
orientation amounts to one Boolean variable: propagate the selected
orientation by \(g\), with no conflicting return because \(r\) is odd.
Write \(m\) for the number of variables.

A directed spanning path candidate is a listing of all vertices,
whether or not compatible with a chosen tournament. The action of \(g\)
on these candidates has free orbits of size \(r\), since a listing fixed
by a power of \(g\) would fix every vertex. A path orbit is an **atom**:
its consecutive pairs prescribe values of some Boolean variables,
by requiring every path in that orbit. Drop any atom prescribing
opposite values to the same variable. The retained atoms are partial
assignments. Two atoms conflict if they prescribe opposite values
somewhere. Pairwise disjointness of conflicts is sufficient for
simultaneous realization, since it is just agreement on their shared
coordinates. This is an ordinary conflict graph on path orbits; neither
its vertices nor its faces are tournaments.

A compatible set \(\mathcal I\) of atoms prescribes \(q(\mathcal I)\)
variables and has exactly \(2^{m-q(\mathcal I)}\) invariant completions.
For the set of available atoms in one completion, each complete orbit
contributes exactly \(r\) Hamiltonian paths. Expanding \(u^r=1+(u^r-1)\)
for each available atom, then reversing the two sums, proves
\[
 F_g(u):=\sum_{T:gT=T} u^{H(T)}
 =\sum_{\mathcal I\ \mathrm{independent}}
       2^{m-q(\mathcal I)}(u^r-1)^{|\mathcal I|}.
\]
Burnside applied separately to the invariant statistic \(H\) gives the
desired orbit polynomial
\[
   F_n(u)=\sum_{\lambda\vdash n\,:\,\text{parts odd}}
              \frac{F_{g_\lambda}(u)}{\prod_{j\ge1} j^{m_j}m_j!},
\]
with \(m_j\) the number of \(j\)-cycles in \(\lambda\).
The requested size series for \(h\) is \(\sum_n z^n[u^h]F_n(u)\).
For that coefficient one can omit all \(g\) with \(r\nmid h\).
This recovers the effect of stabilizers without needing to canonicalize
any tournament in the construction.

For a direct recursive evaluator let the state contain the accumulated
partial assignment \(D\), and drop or keep one atom \(L\). The drop
child retains \(D\), the keep child changes it to \(D\cup L\) and
multiplies by \(u^r-1\). A contradictory keep branch is zero.
Terminate at \(2^{m-|D|}\); \(|D|\) counts assigned orbit variables.
Two atom groups with the same partial assignment can be coalesced:
\(b\) identical atoms become one branch of weight \(u^{rb}-1\).
The recursion is exact, not asserted efficient at large \(n\).
It differs from searching over all orientations one-by-one, since
it contracts partial assignments forced by path families as a unit.

## 2. Color atoms to count literal flips

Fix one ordered pair of named vertices \(s,t\). The first tournament
has \(s\to t\), and the other has \(t\to s\). Both share every other
arc. For simultaneous isomorphism the acting subgroup is all
permutations of the other \(n-2\) vertices, fixing \(s,t\) individually:
the tail and head are intrinsically determined from the first
tournament and the changed pair. Thus no factor two or binomial
multiplicity is needed.

For odd-order \(g\) in this subgroup, remove the marked pair's Boolean
variable. A candidate path orbit either uses the fixed marked arc in
one direction on every member, or avoids it on every member. Give
its atom weight \(u^r-1\), \(v^r-1\), or \((uv)^r-1\), respectively,
according as it belongs to the first side only, the second only, or
both. All other literals, conflicts, and completion weights are as
above. The elementary same expansion proves
\[
 F^{\rm flip}_g(u,v)
 = \sum_{\mathcal I\ \mathrm{independent}}
       2^{m-q(\mathcal I)}\prod_{L\in\mathcal I} w(L),
\]
where \(m\) now excludes the marked pair. Burnside over \(S_{n-2}\)
gives an exact class polynomial at size \(n\) for counts
\((H(T),H(T'))\). In this application the coalescing rule multiplies
the monomials \(1+w(L)\) of atoms with the same assignment before
subtracting one. For \((a,b)\) one can omit every order
\(r\nmid\gcd(a,b)\). In particular for a literal unit seam, **only
the identity term survives**. Equivalently every simultaneous
automorphism acts freely on the paths of each endpoint, so its order
divides the two coprime counts. Thus all literal unit pairs have
trivial simultaneous automorphism group. This does not imply either
unmarked endpoint is rigid individually.

Combining this polynomial with [the common-cut grammar](PRIMITIVE-SEAMS.md)
gives the primitive polynomial without guessing which strong cores
supply it. If the desired scalar target has one first-born class,
that places no analogous restriction on which of its arcs lead to
the neighboring count. The quantity for a perfect-number straddle is
the **joint** coefficient of \(u^{E_p-1}v^{E_p+1}\), not a product
of marginal coefficients.

The replay compares both path-atom transforms against direct orbit
enumeration through order four. This is a non-tournament representation
of the problem, not a proof that the joint coefficient is nonzero
at every exponent.
