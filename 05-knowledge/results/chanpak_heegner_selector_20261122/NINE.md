# Exact displacement between two unit flips on the nine-vertex switching orbit

**PROVED with a finite-exact template cover.** For every prescribed positive integer \(k\) there are unit flips on the same labelled vertex set whose four Hamilton counts are
\[
                 h,\quad h+2,\quad h+2k,\quad h+2k+2 .
\]
At most four other arc reversals connect the lower hosts. The same marked unit flip survives at every intermediate vertex of the resulting cube of orientations. The construction uses \(O(\log(2+k))\) vertices and \(h=O(k+1)\). For a connecting **single** flank arc in the construction below, the possible half-displacements are exactly the positive integers coprime to 15. The integer prescribed here is the *displacement*; the absolute height \(h\) is supplied by the construction.

For \(2k=640320^3+744\), the [target decoder](nine_ribbon.py) uses 68 vertices and three flank arcs. It puts literal unit flips at the even midpoints 119,509,003,775,753,594,498 and 119,771,541,188,394,363,242. Their difference is precisely the nearest integer to Ramanujan's constant. This is a distance between midpoints.

## Why these particular nine vertices matter

Use the inherited switching fixture with edge bits \(\mathtt{0x7103d59dc}\) in lexicographic order on \(0,\ldots,8\); one means the smaller label beats the larger. Its 256 cuts and 36 individual arc choices give 66 unit contacts covering every even midpoint residue modulo 24. The hostile case for modular selection is that these contacts occupy only 29 distinct absolute midpoints. The extra coordinates below are an apex, its two binary flanks, and the binary carry of the changing address. [Earlier switching census and convention](../chanpak_mod24_moonshine_ports_20261122/README.md)

**FINITE-EXACT classification with unbounded flanks.** Choose any cut, any physical arc (allowing the apex as endpoint), and any apex. Hold one exterior address at any nonnegative integer, allow the other address to vary over either entire parity class, and require the signed arc difference to equal a fixed \(+2\) or \(-2\) throughout. Of all \(256\cdot36\cdot9=82,944\) markings, only the following two succeed. Here 134 and 158 are masks on the eight low vertex labels. Both have base count 371 and give 369 after the edit.

| Cut | Apex | Mark | Free address | Opposite address |
|---:|---:|---|---|---:|
| 134 | 7 | \(\{1,2\}\) | left | 0 |
| 158 | 2 | \(\{6,7\}\) | right | 0 |

Both parities work, with exactly the same midpoint function
\[
                C(x)=370+388x+15f(x),
\]
where \(f(0)=0,\ f(2v+1)=f(v),\ f(2v)=f(v)+v+[v\ {m odd}]\). For example the first upper chart has profile \((h,q,c,t,e,r,K)=(15,864,766,388,8,0,371)\), and the edit gives \((15,864,751,388,7,0,369)\). Thus its signed upper-minus-lower difference at arbitrary nonnegative addresses \((x,y)\) is \(2+15y+[y\ {\rm odd}]\). [Literal row constructor and two-flank formula](../chanpak_decimal_kernel_20261122/fermat_receiver.py)

For completeness, with \(x\) free and \(y\) fixed, let \(\delta\) denote the change of one chart coordinate. Constancy on either parity of \(x\) forces \(\delta h=0\): at \(x=0,2,4,6\) the even-branch reads are \(0,2,4,4\), and at \(x=1,3,5,7\) the odd-branch reads are \(0,0,2,0\). Neither is affine. The remaining conditions are precisely
\[
 \delta q\,y+\delta t=0,\qquad
 \delta c\,y+\delta e[y\ {m odd}]+\delta r[x\ {m odd}]+\delta K=\pm2.
\]
Exchanging left and right gives the reflected equations. The native census computes tagged profiles by subset path counts and outputs all 23,180 equal-deletion-height markings. On each parity the two displayed linear equations permit zero, one, or all fixed addresses. The exact classifier checks all three possibilities. No nonempty two-dimensional parity face survives. These quantifiers concern a single fixed arc and apex in this specified orbit.

## Locally universal and globally sparse

**PROVED.** \(C(x)\) attains every even target modulo every modulus. Its absolute image has upper density at most \(1/194\) among the even integers, since \(f(x)\ge0\) and a midpoint at most \(Z\) requires \(x\le(Z-370)/388\). On a positive absolute target \(E\), the necessary exact interval and phase are
\[
  (E-370)/403\le x\le(E-370)/388,\qquad x\equiv5-E/2\pmod{15},
\]
using \(f(x)\le x\) and evenness of \(f\). Testing the retained addresses in \(C\) suffices.

Here is an explicit joint modular lift. Let \(M=2^a m\) with \(m\) odd. Choose \(L\ge\max(1,a)\) with \(2^L\equiv1\pmod m\). The binary word of \(L\) ones followed by \(L\) zeros represents \(B=(2^L-1)2^L\equiv0\pmod M\). It has \(f(B)=(2^L-1)^2+1\), a unit at every prime dividing \(m\), and congruent to two modulo \(2^a\). Repeating this word keeps the address at zero modulo \(M\) and its parity even; each copy adds the same \(f(B)\) modulo \(M\). This follows from the [binary split identity](../chanpak_decimal_kernel_20261122/REVERSED-PHASE.md), displayed in the next section. Append any suffix of prescribed address residue. Its read is added to the freely adjustable register. Thus \(x\bmod M\) and an arbitrary even \(f(x)\bmod M\) can be prescribed jointly. The assertion follows from \(\gcd(388,30)=2\). The construction specializes the earlier [general local-word controller](../chanpak_decimal_kernel_20261122/LOCAL-WORDS.md) to a single repeating block.

## The high address corrects an absolute displacement

Adding one to bit \(j\) that was zero, with \(u=\lfloor x/2^{j+1}\rfloor\) and lower \(j\) bits fixed, changes \(C\) by \(388\cdot2^j\) modulo 30. Here and below a change of any same-width flank bit is one physical apex arc. As \(-194\equiv1\pmod{15}\), the unsigned half-displacement runs among \(\pm2^j\pmod{15}\), the eight units. Conversely, pick a sign of \(2k\) and a phase of \(j\) so \(D=\pm2k\equiv388\cdot2^j\pmod{30}\). Choose \(j\ge1\) in that phase large enough and set all lower \(j\) bits to zero. At the addresses \(2^{j+1}u\) and \(2^{j+1}u+2^j\) with positive even \(u\), the signed difference is \(403\cdot2^j-15u\), by the binary recurrence. The division \(u=(403\cdot2^j-D)/15\) is exact and even. This proves the one-bit classification.

For all widths, write the two addresses \(2^Lu+a\) and \(2^Lu+b\) with the same positive high word. With \(R_L\) denoting reversal of precisely \(L\) bits, the split identity is
\[
 f(2^Lu+a)=f(u)+f(a)+(2^L-1-R_L(a))u+[a<2^{L-1}][u\ {m odd}].
\]
Set \(d=R_L(b)-R_L(a)>0\) and \(u=2v+e\). The first midpoint minus the second, divided by two, is **exactly** \(15d v+I\), where
\[
 2I=15de+388(a-b)+15\bigl(f(a)-f(b)-e([b<2^{L-1}]-[a<2^{L-1}])\bigr).
\]
The two addresses differ in as many flank arcs as the Hamming distance of \(a,b\); every mixture of those bits is still in the same unrestricted unit family.

The certificate uses lengths 1 through 8, at most four changed bits, \(d\in\{1,3,5,15\}\), and \(e\in\{0,1\}\). Each pair yields all signed half-displacements congruent to \(I\pmod{15d}\) starting at \(I+15d[e=0]\). Keep the least starting value for each progression. Since \(15d\mid225\), a scan of 225 residues decides the union above a finite threshold: every signed half-displacement above 3030 occurs in the positive orientation. For each of 1 through 3030 at least one orientation occurs. The [replay](nine_ribbon.py) computes these bounds from the suffixes, chooses the template from the target, and recovers \(v\) by exact division. The modulus and the number of changed bits are independent of \(k\). Because the retained templates are finite, \(u=O(k+1)\) and both claimed growth bounds follow. Minimal connecting distance for nonunit gaps remains open.

## What the five cut shells preserve

For a nine-vertex tournament \(T\), let \(S_s\) be the sum of \(H(T^A)\) over all \(s\)-element vertex switches for \(0\le s\le4\). Each cut occurs once. Switching \(s\) vertices flips \(w_s=s(9-s)\) of the 36 arcs. Its five possible distances are \(0,8,14,18,20\). Writing \(t=4-s\) gives the exact reflection
\[
       w_s=20-t(t+1),\qquad t^2+t+p=p+20-w_s .
\]
In particular, for \(p=41\) the five shell values of the quadratic in increasing \(t\) are 41, 43, 47, 53 and 61. The five constants \(p=3,5,11,17,41\) give \((p-1)/2=1,2,5,8,20\) and \(4p-1=11,19,43,67,163\). The integers 8 and 20 really are the minimum nonempty and maximum cut sizes in this gauge. This reflection preserves the number of reversed arcs, and discards which arcs they are.

**PROVED.** The five shell means \(S_s/\binom9s\) lie on a cubic of degree at most three in \(w_s\). In particular every tournament of order nine satisfies the exact equation
\[
                 9H(T)=7S_1-5S_2+3S_3-S_4.
\]
By the [switching forest expansion, THM-5313](../../../01-canon/theorems/THM-5313-switching-derivatives-of-hamiltonian-paths-and-cycles.md), a switch character has even support consisting of the leaves of disjoint even-length paths. On nine vertices the possible support sizes are \(0,2,4,6\). The mean of a support-\(2j\) character on the \(s\)-shell is the coefficient of \(z^{2j}\) in \((1-z)^s(1+z)^{9-s}\), divided by \(\binom9{2j}\). It is a polynomial of degree at most \(2j\) in \(s\), unchanged under \(s\leftrightarrow9-s\), hence a polynomial of degree at most \(j\) in \(w_s\). For the displayed identity, each of the nine weight-eight characters is absent. Their sum on a size-\(s\) cut is \((-1)^s(9-2s)\), so their summed orthogonality relations give the integer weights shown.

There is a local integer recovery rule as well. Fix one gauge sign and index the 256 cuts by \(u\in\{0,1\}^8\). Multiplying the observations by \((-1)^{|u|}\), the nine missing-character checks become the constant and the eight first-coordinate moment equations (take differences and divide by two over the integers). Any three missing true observations can be recovered over the rationals. For four, recovery from these checks holds precisely when their four bit vectors do not form an affine rectangle. Indeed the check columns are \((1,u)\), and no three distinct binary vectors are collinear. In a dependence of four distinct vectors every coefficient is nonzero, their sum is zero, and every nonconstant coordinate picks a pair of opposite coefficients. To distinguish all four columns needs two different pairings, forcing two equal positive and two equal negative coefficients up to sign and scale. This says exactly that two opposite corner sums coincide. This is a recovery criterion using these linear checks. Actual tournament constraints can carry additional information.

**FINITE-EXACT modular selection.** In the fixed switching fixture the 256 cuts realize 125 distinct heights from 59 to 3255. The least modulus distinguishing all 125 *already attainable* values is 699. From \(\{3,5,11,17,41\}\), three primes suffice, for example 3, 11, 41, and no pair does. From \(\{11,19,43,67,163\}\) two suffice, for example 11 and 163, and no singleton does. The five moduli \(1,2,5,8,20\) jointly see only residues modulo 40, insufficient here. An injective residue lookup selects the unique candidate height in this finite image, whose equality with the prescribed integer is then checked directly. The open Chan–Pak assertion concerns all odd heights across all orders.

Run from this directory:

~~~sh
c++ -std=c++17 -O3 -Wall -Wextra -Werror nine_switch.cpp -o /tmp/nine_switch
python3 -B nine_decoder.py --native /tmp/nine_switch
python3 -B nine_ribbon.py 36 --rows
python3 -B nine_ribbon.py 131268706320384372
~~~
