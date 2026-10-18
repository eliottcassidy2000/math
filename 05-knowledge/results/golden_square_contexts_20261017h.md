# Golden square contexts, a signature obstruction, and a controlled field escape

2026-10-17, clock_duality_oct15. **PROVED, scoped, self-audited:** the
signature obstruction and the explicit quadratic escape. **FINITE-EXACT:**
four complete induced graph extensions, with positive boundary witnesses
or constructive extension rules. **OPEN:** six-colour forcing and the
chromatic number of the plane. Field degree alone proves no colouring bound.

## 1. Fix the inherited source before changing the operation

Use [shared golden contacts](shared_golden_contact_20261017g.md) verbatim:
F=Q(zeta15,j,w), j^2=-11, w^3+2w^2-4w-4=0 with1<w<2, and[F:Q]=48.
Put alpha=zeta5, d0=1+alpha, u=(w+3-j)/(w+3+j), Z=d0+u,
N=Z conjugate(Z), D=1-4/N, s^2=D, and

    v=Z(-1+s)/2,       E=F(s),       [E:Q]=96.

The intended embedding has37/15<N<59/20<4, so s is imaginary and
conjugation sends s to -s. The relative automorphism sigma fixes F and
also sends s to -s, but fixes the coefficients instead of conjugating them.
These are different operations. The original source is G=A+uC+vC, with
5400 vertices,30856 edges, and the same eight ordered literal ports P.

Closest mechanism: the complete unit-direction classification in K+uK+vK
and the inherited570 boundary-extension witnesses. Canonical hostile:
all those patterns survive even though the old relative six-core is large.
Corrected near miss: a new field coordinate or a new mixed edge need not
alter the marked boundary relation. Least-used operation: quadratic field
conjugation fixing **every old port**, rather than another independent copy.
This is the geometric niche in the [signed outer-packet portfolio](outer_packet_golden_clocks_20261017h.md).

## 2. Three operations tested on the complete source

The exact new induced graphs are:

| source | vertices | edges | common vertices of the two copies | conclusion |
|---|---:|---:|---:|---|
|G union(G+v^2)|10800|67112|0|Only5400 matching edges join the copies.|
|G union vG|10790|61697|10|Every edge already belongs to one of the copies.|
|G union sigma(G)|10260|59442|540|Every edge already belongs to one of the copies; there is a fold onto G.|

All three have chromatic number exactly4: the receipts give a four-colouring,
and the old spindle remains. All570 original six-colour port patterns extend.
The first and third conclusions are stronger than a570-pattern census:

* For G union(G+v^2), extend any k-colouring of G, k>=2, by adding1 modulo k
  on the translated copy. The matching edges stay proper.
* For G union sigma(G), the intersection is exactly A+uC, fixed pointwise
  by sigma. Identity on G and sigma on its second copy define an edge-
  preserving fold onto G. Thus every colouring of G extends, with every
  original vertex fixed.

The rotated source has only the finite six-colour relation certified here;
no analogous general retraction is asserted. Its two copies overlap in
vC, whose15 edges explain30856*2-15=61697. The conjugate copies share2270
edges, explaining30856*2-2270=59442. The full pair sieves are independent
of these edge-count predictions: they retain68673,63245,60908 candidates,
then decide respectively109,209,132 distinct characteristic-zero norms.

The v-squared identity does leave the old linear coefficient space:

    v^2=-Zv-Z/conjugate(Z).

Nevertheless the selected square layer is a Cartesian matching context.
The missing interaction is now explicit rather than inferred from degree.

## 3. Why a point cannot be joined to its quadratic mate in E

**G1 (PROVED).** No nonzero element x of the anti-invariant subspace F*s
has x conjugate(x) in Q_(>0). In particular no x in E is at unit distance
from sigma(x).

Proof. At the intended complex embedding, x=b*s gives
x conjugate(x)=-D b conjugate(b)>0 for b!=0. Consider the F-embedding tau
that fixes zeta15,w and sends j to -j. It sends u to conjugate(u), so it
changes N to |d0+conjugate(u)|^2. The same rational intervals used in the
inherited real-branch certificate give

    |d0|^2>13/5,    Re(d0)>13/10,    Im(d0)>9/10,
    Re(u)>1/3,      -Im(u)>9/10.

Consequently

    tau(N)>13/5+1+2[(13/10)(1/3)+(9/10)^2]
          =913/150>4,                                    (1)

and tau(D)>0. On F this embedding commutes with ordinary complex
conjugation; hence tau(b conjugate(b))=|tau(b)|^2>0. Thus

    tau(-D b conjugate(b))<0,

which is impossible for a positive rational. For a point x=a+b*s,
x-sigma(x)=2b*s, proving the distance claim. This excludes an entire
matching mechanism throughout the field, not only in the finite graphs.
It does not exclude unit edges joining different sigma-orbits, nor prove
that the whole field is four-colourable. The orbit quotient has no loops;
its other edges can still be difficult.

## 4. Change the exact obstructing signature, with the source retained

Set R=4-N. Adjoin a square root h of R that is positive at the intended
embedding, and define

    e=(v-sigma(v))/h=Z*s/h,
    e^2=-Z/conjugate(Z)=:r,       e conjugate(e)=1.          (2)

This is a purpose-built new unit direction. Extend sigma by sigma(h)=h;
then sigma(e)=-e. Hence e/2 and -e/2 are genuine unit-distance mates.

**G2 (PROVED degree and escape).** E'=F(s,h)=E(e) has degree192 over Q.
It preserves every old coordinate and every old unit edge, while the
matching prohibited by G1 is now possible.

Two exact degree-one places of F at31 certify the degree:

| (zeta15,j,w) mod31 | N | D | R | quadratic characters(D,R) |
|---|---:|---:|---:|---|
|(7,19,25)|12|11|23|(-1,-1)|
|(10,19,25)|28|23|7|(-1,+1)|

All defining polynomials vanish, all denominators are units, and both rows
are direct modular evaluations. The square classes of D and R are
independent, so adjoining their roots has relative degree4 over F.
Equivalently the square class of r is that of -N=R/D. This proves that e
is a genuine new quadratic direction over E, rather than a renamed old
coordinate. The real interval gives R>0 at the selected embedding; at tau,
R<0. The positivity premise of G1 is therefore no longer available in
the enlarged coefficient field. There is no conflict with that theorem.
This degree192 field must not be identified with a different closed carrier
merely because their degrees coincide.

### An authenticated attachment to every old port

Let a=P_5=(3-j sqrt(5))/8, so |a|=|a-3/4|=1. The seven-point walk

    L_+=(0,1,1-a,1/4,5/4,5/4-a,1/2)

has six unit steps. Put L=L_+ union(-L_+), with its shared0. It has13
vertices. Its induced unit graph is exactly the odd13-cycle formed by the
two six-step walks and the edge from1/2 to-1/2. Multiplying by e supplies
an explicit unit path from each old port p to p+e/2, and its sigma-image
ends at p-e/2. Their joining edge cannot be preserved by the old sigma-fold.
Thus the added edge is attached by authenticated walks, not left isolated.

The tested source retains the entire old graph:

    G'=G union(P+eL).

**FINITE-EXACT:** G' has5496 vertices and31068 edges,212 of them new. Every
new pair is tested in E', with the original ports unchanged. The new part
is exactly the Cartesian product of the port graph and the13-cycle, glued
along P. Its four-colouring and all570 boundary extensions are constructive:
choose a three-colouring c_L of the cycle with c_L(0)=0, and use

    c(p+e*l)=c(p)+c_L(l) mod k,       k>=3.                (3)

Every colouring of the **whole old G** with k>=3 colours therefore extends.
The new graph is exactly four-chromatic, and six-colour forcing is still
absent. The new exterior also peels completely by degree-at-most-five
removal, another certificate that it cannot kill an old six-colouring.

This is a controlled reroute with a recorded boundary: it defeats the
whole-field anti-invariant-unit obstruction but remains Cartesian on the
selected source. The next geometric target is a unit edge joining different
port and e-coefficient fibres, together with a nonempty relative six-core;
even that necessary structure will still need an authenticated restrictive
boundary relation. Adding more disjoint or Cartesian copies cannot do it.

## 5. Relation to the Collatz lane and reproducibility

The common move is to expose the failed matching before changing generators.
For the outer Collatz word, terminal columns are isolated but signed values
may cancel across columns. For the old golden field, sigma-mate unit edges
are forbidden by a real-embedding sign; adjoining h changes precisely that
condition. The new finite geometry still has an explicit colour extension,
just as the positive Collatz kernel still lacks the required quantitative
signed estimate. These are two applications of a research method, not a
map transporting a colouring theorem into Collatz convergence.

The [manifest](../../04-computation/experiments/outer_golden_20261017h/manifest.json)
records scripts, exact receipts, and hashes. Rung, fold, and Cartesian
extension witnesses are stored as reversible recipes referring to the
inherited graph; the rotated case stores each complete colouring. The
separate verifier reconstructs all570 assignments for every graph and
checks every edge and literal port. Good-prime reduction only rejects
candidates; retained new differences receive exact characteristic-zero
norm tests. Normal and optimized runs are compared. No numerical search
failure is used as a non-colourability proof.
