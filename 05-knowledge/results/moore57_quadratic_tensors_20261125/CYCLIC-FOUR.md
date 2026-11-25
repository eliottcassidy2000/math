# A single sheet cycle permits at most three fixed stars

**PROVED.** In a rooted Moore permutation chart with \(q\) sheets per branch, suppose a symmetry fixes some branch heads and acts by a \(q\)-cycle on one fixed branch. For even \(q\) it fixes at most **three** heads, and the triangle transport on any fixed triple is precisely the half-turn. For odd \(q\) it fixes at most **two**. The argument needs only the complete matchings incident to the fixed stars and their saturated pair boards, so it applies during construction before completing any exterior-exterior matchings.

In particular no permutation of the outside branch labels and no change from translations to affine or nonlinear multipliers can develop four individually fixed whole stars from one source site at degree 57. The original three-star Schur prototype attains the bound for a partial chart: translation by one fixes all three central stars and its center triangle acts by 28. Its other fixed heads do not have completed stars. Even cycling the three original central heads cannot evade the obstruction if the source of the proposed fourth is fixed and has a 56-cycle: the cube of the symmetry fixes all four heads and still has a 56-cycle there. [Three-star prototype](three_stars.py) · [Complete branch and pair requirements](../moore57_thirteen_20261125/LATIN.md)

## The phase of a fixed triangle

More generally suppose the generator acts on one fixed fiber in \(r\) cycles, all of length \(t\), where \(q=rt\). All fixed fibers have the same cycle structure because their matchings intertwine the actions. Synchronize two fixed fibers \(U,V\) by their direct matching. Their labels can now be written \((i,z)\) with \(i\in\{1,\ldots,r\}\) and \(z\in\mathbb Z/t\); the generator adds one to \(z\). Transport through a third fixed branch \(B\), returning to the common label set, is some permutation
\[
             P(i,z)=(\pi(i),z+d_i).
\]
Its total shift \(\Theta(P)=\sum_i d_i\pmod t\) is independent of how the cycle origins and cycle indices are chosen, since changing origins adds a difference telescoping over \(\pi\).

**Fixed-triple identity.**
\[
                \Theta(P)=r^2\sum_{j=0}^{t-1}j\pmod t.             \tag{1}
\]
To prove it, pair each noncentral exterior site with its two neighbors in \(U,V\). Saturation of the \(U,V\) board identifies these sites bijectively with all ordered label pairs except the direct pairs \((a,a)\) and the pairs \((a,P(a))\) whose route goes through \(B\). Every orbit of pairs has length \(t\). For each ordered pair of cycle indices there is exactly one pair orbit at each relative phase \(j\in\mathbb Z/t\), measured as second phase minus first. The direct pairs remove zero phase. The \(P\)-pairs remove total phase \(\Theta(P)\).

Pick one source site per cycle in \(B\). The union of their exterior neighborhoods takes one representative of every remaining pair orbit. Indeed the unique matched neighbor in \(B\) of the sites in any one orbit runs exactly once through a complete source cycle. For each source representative, let its two direct neighbors in \(U,V\) have common-coordinate labels \(a,b\). The projections of its exterior neighbors onto \(U,V\) each enumerate *precisely the complement of the same two labels* \(a,b\). The direct route to either central branch and the route via the other account for those two holes. Therefore the sum of second phases minus first phases is zero already at each source. Summing over source representatives proves (1). The two-hole identity for a single board was established in [the local Latin obstruction](LOCAL-LATIN.md); accounting for entire symmetry orbits gives the fixed-triple restriction.

If \(r=1\), then \(P\) is translation by \(\Theta(P)\). Its shift must be nonzero by the direct route. The sum of all elements in \(\mathbb Z/q\) is zero for odd \(q\) and \(q/2\) for even \(q\). A fourth fixed branch would be a second intermediate with the same shift, violating uniqueness. For the original proposed one-representative-per-offset construction the contradiction is also visible numerically: with \(D=\mathbb Z/56\setminus\{0,28\}\), the chosen original offsets sum to \(\sum_{h\in D\setminus\{b\}}h=-b\) modulo 56, while the prescribed holes require 28. Formula (1) explains why permuting the original offsets by a full-cycle action cannot repair this.

The same orbit count works when a finite group acts regularly on each fixed sheet. In left regular coordinates the invariant relative displacement of a pair \((u,v)\) is \(u^{-1}v\). Sum its class in the abelianization over pair orbits. The equal projected label lists still cancel. The displacement \(d\) of a fixed third branch hence obeys \(\bar d=\sum_{g\in G}\bar g\) in \(G_{\rm ab}\). In particular an arbitrary regular **abelian** symmetry permits at most three fixed branches; having three requires a unique involution. For a nonabelian group the condition constrains a coset of the commutator subgroup.

## Smaller translations on the Schur prototype

Suppose translation by \(r=56/t\) acts on the central sheets, with \(t\mid56\), and the prospective fixed branch \(b\) retains its three original incident central maps. Other exterior branch colors and the proposed new edges are free subject to equivariance and the matching equations. A translation of the slot indices by \(d\) has total shift \(d\bmod t\): use cycle representatives \(0,\ldots,r-1\) and sum \(\lfloor(i+d)/r\rfloor\) over them. The direct transports via the prototype branch \(b\) have shifts \(-b\) and \(-s(b)\). Both must agree with the right side of (1). Testing this necessary filter gives the following.

| Order \(t\) of the translation | Required total shift modulo \(t\) | Source branches passing |
|---:|---:|---|
| 56 | 28 | none |
| 28 | 0 | none |
| 14 | 0 | none |
| 8 | 4 | 12, 44 |
| 7 | 0 | none |
| 4 | 0 | 8, 12, 20, 36, 44, 48 |
| 2 | 0 | 2, 8, 12, 14, 16, 20, 36, 40, 42, 44, 48, 54 |

Thus the largest surviving orbit length is **eight**, with seven source orbits, at \((b,s(b))=(12,20)\) or \((44,36)\). Exploiting this choice would still require an equivariant recoloring of the other exterior branches and mutually compatible placements for the seven source representatives. Leaving every original translation row untouched fails the mod-56 sum even for one site. The table is a filter on modified equivariant constructions, not a fourth-star witness.

To reproduce from this directory:

~~~python
from three_stars import A, B, SLOTS, make_chart
maps = make_chart(())
for t in (56, 28, 14, 8, 7, 4, 2):
    r = 56 // t
    phase = (r*r*t*(t-1)//2) % t
    assert 28 % t == phase
    keep = [b for b in SLOTS if maps[b,A][0] % t == maps[b,B][0] % t == phase]
    print(t, phase, keep)
~~~
