# Minimum-order flips at discriminant eleven and an inverse for equal-height deletions

## Two prescribed midpoints, each at its first possible order

**FINITE-EXACT.** The CM integer \(32^3+744=33{,}512\) and the nearest integer \(33{,}506\) to \(e^{\pi\sqrt{11}}\) both occur as midpoints of single-arc reversals on eleven vertices. The two rounding conventions produce different inputs. For each row, build an apex over the ten-vertex deletion. Bit \(j\) of outgoing row \(i\) means \(i\to j\), and the apex beats the vertices in the mask. Toggling the displayed bit changes the lower graph into the upper.

| Midpoint | Lower and upper counts | Deletion rows | Mask | Bit |
|---:|---:|---|---:|---:|
| 33,506 | 33,505 and 33,507 | 236,189,136,272,933,76,158,296,103,495 |625|8|
| 33,512 | 33,511 and 33,513 | 712,821,817,374,129,465,918,270,529,184 |797|6|

**PROVED using the cited finite maximum.** Eleven is the minimum order for every one of these four counts. The maximum at order ten is 15,745, and adjoining a universal sink embeds each smaller spectrum in the order-ten spectrum. Moreover every displayed endpoint has at least two nonisomorphic witnesses at its minimum order. Reversing all arcs preserves the count by path reversal and replaces each score \(s\) by \(10-s\). Both graphs in the first row have sorted scores \((3,3,3,4,5,5,5,6,6,7,8)\). Both in the second have \((2,3,4,4,5,5,6,6,6,7,7)\). Each list differs from its reflected list. Direct counts, arc orientations and the score comparisons take a fraction of a second to reproduce with `python3 -B minimum_eleven.py`. [Gordon Royle's extremal tournament counts](https://mathoverflow.net/questions/415359/maximum-number-of-hamilton-paths-in-a-tournament-on-n-vertices)

**FINITE-EXACT, inherited.** If \(A_h(z)=\sum_n a_n(h)z^n\) counts tournament isomorphism classes by order at fixed height \(h\), the four leading coefficients at \(z^{11}\) are the following. All these minimum-order classes are strong. Summing the rows at each height in the saved order-11 census over the cyclic-triangle levels gives the indicated totals and zero for the nonstrong classes. The new edges above identify particular neighboring pairs within these four large sets. [Exact order-11 census](../path_spectrum_two_angles_metagraph_20261106/data/census_11.txt.gz) · [Least-order and multiplicity table](../path_spectrum_two_angles_metagraph_20261106/data/least_order_11.txt) · [THM-4889, Hamiltonian-path spectrum by order and by value](../../../01-canon/theorems/THM-4889-the-hamiltonian-path-spectrum-by-order-and-by-value-exact-censuses-to-order-11-resolved-by-cyclic-triangles-least-orders-strong-windows-level-floors-and-the-zone-law.md)

| \(h\) | 33,505 | 33,507 | 33,511 | 33,513 |
|---|---:|---:|---:|---:|
| \(a_{11}(h)\) | 37,863 | 38,660 | 38,257 | 38,645 |

Converse acts as an involution on each of these four sets. Its two-element orbits contribute an even number of classes. Therefore the odd multiplicities at 33,505, 33,511 and 33,513 force an odd number of self-converse minimum classes at each of those heights. The displayed representatives at those heights also give explicit two-element orbits. At 33,507 the number of self-converse classes is even.

## The geometry of a prescribed change

**PROVED.** Take two physical marked-apex charts related by reversing one arc and having the same deletion height \(h\). This includes equal-height internal edits as well as apex edits. Specify positive left and right addresses \(x,y\), their parities, and a signed change \(g\). In one chamber the old height is
\[
 H_0=K+qxy+cy+tx+h(f(x)+f(y)),
\]
where the fixed parities are absorbed in \(K\), the coefficients \(q,c,t\) are nonnegative, and \(q\ge2h\). The binary recurrence is \(f(0)=0,\ f(2w+1)=f(w),\ f(2w)=f(w)+w+(w\bmod2)\). Thus \(0\le f(w)\le w\). Subtracting the charts cancels both binary reads. With \(a=\Delta q,b=\Delta c,c'=\Delta t,d=g-\Delta K\), the exact change equation is \(axy+by+c'x=d\). [Signed physical flank law](../chanpak_decimal_core_selection_20261117/FLANKS.md)

* If \(a\ne0\), then \((ax+b)(ay+c')=ad+bc'\). A nonzero product gives a finite list of signed divisors. Recover the two addresses by division by \(a\), requiring integrality, positivity and both parities. A zero product gives coordinate lines, which are single-variable digit readers.
* If \(a=0\) and \(bc'\ne0\), solve \(by+c'x=d\) with the two parities by the Euclidean algorithm. Coefficients of the same sign give a finite interval of solutions. Opposite signs give one increasing ray or the empty set. The other linear cases give coordinate lines, the empty set, or the entire chamber.

On an increasing ray write \((x,y)=(x_0+i u,y_0+i v)\) with positive increments. Set \(P(i)=K+qxy+cy+tx\). The height belongs to \([P(i),P(i)+h(x+y)]\). Consecutive intervals are strictly disjoint because their baselines increase by at least \(q(uy+vx+uv)>h(x+y)\). Given an exact old height \(M\), one integer square root selects the largest \(i\ge0\) with \(P(i)\le M\), if one exists. Only this candidate needs the digit read. This extends the spoke inverse to a classification by the *change in the two-slot count* across arbitrary equal-height edited deletions. [Preceding spoke inverse](../chanpak_perfect_spoke_inverse_20261122/README.md)

Consequently, for any fixed finite catalogue, the part arising from finite sets and increasing oblique rays represents only \(O(\sqrt X)\) midpoints up to \(X\). More precisely, each ray contributes \(\sqrt{X/(quv)}+O(1)\) parameters, since its exact midpoint is \(quv i^2+O(i+1)\) and is strictly increasing. Count ray parameters with their labels for the additive leading coefficient. Different rays may represent the same midpoint. This leaves coordinate lines, whole chambers and changing core families as distinct routes to denser images.

## The specified nine-vertex switching orbit

**FINITE-EXACT.** In all 23,180 oriented equal-deletion-height charts of the specified switching orbit, examine both changes \(g=\pm2\), four address parities, and arbitrary **positive** binary addresses. The finite branches give 214 labelled points, all from internal edits. Their largest midpoint is 198,689,928. The oblique branches give 91,008 labelled rays, including 2,272 from internal edits and 88,736 from spoke edits. For this entire part of the atlas the exact nearest midpoint offsets are as follows. Negative entries lie below the target. The variable address lengths have no cutoff. Coordinate-line and full-chamber edits are outside this census. [Switching orbit and axis analysis](NINE.md)

| Target center | Nearest midpoint minus target |
|---:|---:|
|33,506|48|
|33,512|42|
|885,480|−262|
|884,736,744|21,552|
|147,197,952,744|−28,090|
|262,537,412,640,768,744|−142,493,648|
|\(2^{60}(2^{61}-1)\)|453,678,458,321,142,734|

For nearest-midpoint recovery the intervals give an ordered sequence. If the one square-root candidate is below the target, compare its successor. If it is above, compare its predecessor when present. A target below the first baseline needs just the first midpoint. The reproduction command also prints each nearest physical configuration. The two eleven-vertex graphs above supply literal crossings where this whole oblique family skips over the target.

~~~sh
c++ -O3 -std=c++17 nine_switch.cpp -o /tmp/heegner-nine-charts
/tmp/heegner-nine-charts > /tmp/heegner-nine-bank
python3 -B oblique_inverse.py /tmp/heegner-nine-bank
~~~
