# Two odd labels for every primitive right triangle

**PROVED.** Every primitive right triangle carries exactly two positive odd integer labels in the construction below. Multiplication by three preserves its absolute leg difference; the other two neighboring branches replace that difference by the old leg sum, with opposite signs. The first nonunit difference is 7. Its entire set of labels is
\[
                 3^j\{7,11,13,17\},\qquad j\ge0.
\]
Exactly 7 and 21 are absent Hamilton path counts in that sector. More generally we determine the tournament spectrum on all 594 infinite tripling rays with primitive leg difference at most 1,200.

## The branch rule and the forgotten digit

Put \(k=(h-1)/2\) for an odd positive label \(h\). The maps \(3h-2,3h,3h+2\) append the ordinary ternary digits \(0,1,2\) to \(k\). For Euclid parameters \((m,n)\) use the same digit to apply
\[
 B_0(m,n)=(m+2n,n),\quad B_1(m,n)=(2m+n,m),\quad B_2(m,n)=(2m-n,m).
\]
For \(m>n>0\) coprime of opposite parity, the odd leg, even leg and hypotenuse are \((a,b,c)=(m^2-n^2,2mn,m^2+n^2)\). These are the Berggren branches on this choice of Euclid parameters. For example, the three children of \((2,1)\) give \((15,8,17),(21,20,29),(5,12,13)\) in our digit order. At slopes \(m/n>3,\ 2<m/n<3,\ 1<m/n<2\) the unique parent is respectively \((m-2n,n),(n,m-2n),(n,2n-m)\). Each inverse preserves the primitive mixed-parity condition by unimodularity and decreases the coordinate sum. At slope two the only primitive pair is \((2,1)\); slope three is excluded by parity. [Berggren descent, Keith Conrad](https://kconrad.math.uconn.edu/blurbs/linmultialg/descentPythag.pdf) · [Earlier odd-root chart and its physical boundary](../chanpak_graceful_ternary_20261115/BERGGREN.md)

Start at the cusp \((1,0)\). Digit 0 fixes it, while digits 1 and 2 both send it to \((2,1)\). Reading \(k\) from most to least significant digit now gives a definite triangle for every \(h>1\); the label 1 denotes only the cusp. The inverse just proved recovers all digits after the first. If its length is \(d\) and its ternary value is \(v\), the two labels are precisely
\[
               2(3^d+v)+1,\qquad 2(2\cdot3^d+v)+1 .
\]
Thus the leading digit is the entire lost information. For example, 7 and 13 name the same triangle \((15,8,17)\); 21 and 39 both name \((65,72,97)\). In each case one twin is a Hamilton count and the other is excluded. The earlier [physical ternary code](../chanpak_berggren_square_bridge_20261117/root/PROOF.md) assigns a distinct realized count to each encoded word at fixed depth. Here the labels are the prescribed odds themselves, so the lost digit matters to tournament realizability.

**PROVED, affine classification.** Consider an integral affine map \(L(h)=Ah+b\), \(A>0\), taking odd labels to odd labels. It maps each pair of twin labels to a pair of twins for all sufficiently large inputs if and only if
\[
                 A=3^r,\qquad b\in2\mathbb Z,\qquad |b|\le3^r-1.
\]
The original twin gap is \(2\cdot3^d\); the new gap forces \(A=3^r\) even if the requirement is made only for one pair. On ternary coordinates the map becomes \(k\mapsto3^r k+C\) with \(C=(3^r+b-1)/2\). At depth \(d\), input \(k=3^d+v\) must map into the smaller leading-digit cylinder \([3^{d+r},2\cdot3^{d+r}-1]\). Taking \(v=0\) and \(v=3^d-1\) gives \(0\le C<3^r\). Conversely every such \(C\) appends its padded word of length \(r\), preserving twinhood. In particular the only ordinal-product multipliers on Hamilton counts which act on the triangle quotient are powers of three, all realized by ordinal sums with directed triangles. Every other positive integral multiplier separates every twin pair.

## Why the first layer has difference seven

Write \(\delta=a-b=(m-n)^2-2n^2\). Direct substitution gives the three signed child differences, in digit order,
\[
                   a+b,\qquad -\delta,\qquad -(a+b).
\]
At the root these are \(7,1,-7\). The two outer branches create the same new magnitude. On the middle branch a triangle simply moves along its existing difference ray. Since every nonroot leg sum exceeds 7, the only beginnings of difference-seven rays are the two outer children of the root. Their scalar twin labels are \((7,13)\) and \((11,17)\). Tripling is the middle branch, establishing the claimed full list of labels. Literal tournaments with counts 11,13,17 and 63 appear in the replay. Ordinal multiplication by \(H(C_3)=3\) covers every later label on these four rays after 63 on the first. The two global exceptions are inherited. [THM-1370, exclusion of 7 and 21](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) · [Earlier sevenfold ray](../chanpak_moonshine_prime_blocks_20261118/SEVENFOLD_GUARDED_LIFT.md)

The odd leg \(a\) is prime precisely when \(m-n=1\) and \(a=m+n\) is prime. Then \(b=(a^2-1)/2\), and \(a+b=2a+1\) forces \((a-3)(a+1)=0\). Hence the root is the unique primitive triangle with prime odd leg \(p\) whose leg sum is \(2p+1\). For \(p=3\) this is genuinely the origin of the new gap 7. The value \(21=2\cdot11-1\) is a scalar tripling of 7 on that gap. It cannot itself be a primitive leg sum or difference, as the next theorem shows.

## Complete classification of difference rays

**PROVED.** A positive integer \(g\) is a primitive absolute leg difference exactly when it is odd and every prime divisor is \(1\) or \(7\pmod8\). For each admissible \(g\) there are exactly \(2^{\omega(g)}\) geometric rays and \(2^{1+\omega(g)}\) scalar tripling rays, including \(\omega(1)=0\).

Use \(z=(m-n)+n\sqrt2=s+n\sqrt2\) in \(\mathbb Z[\sqrt2]\). Its norm is the signed leg difference, and the middle branch multiplies it by the unit \(1+\sqrt2\). Since \(\gcd(s,n)=1\), an inert prime cannot divide its norm: if either coefficient vanished modulo that prime, the other would have to vanish too; otherwise 2 would be a square. An odd prime splits precisely when \((2/p)=1\), equivalently \(p\equiv1,7\pmod8\). The prime 2 cannot divide the norm since \(s\) is odd.

For completeness the ring is norm-Euclidean: approximating two rational coefficients by nearest integers gives errors \(|u|,|v|\le1/2\), hence \(|u^2-2v^2|\le1/2\). For a split \(p\), the two evaluation ideals sending \(\sqrt2\) to the two distinct roots mod \(p\) are distinct, relatively prime and have index \(p\). Their product is \((p)\). Unique factorization of ideals now gives a choice of one of the two primes to receive the *entire* exponent of each \(p\mid g\), since using both would make \(p\) divide both coordinates of \(z\). Every such choice gives a principal ideal of index \(g\), with a primitive generator. Conversely the index of \((z)\) is the absolute value of the determinant of multiplication by \(z\), hence equals its absolute norm.

Normalize the generator to be positive in the real embedding with \(\sqrt2>0\). All positive units are powers of \(1+\sqrt2\). Indeed, scaling a unit into \([1,1+\sqrt2)\) makes its conjugate of absolute value at most one. In that range its integer coefficient of \(\sqrt2\) can only be 0 or 1; in the latter case its integer rational coefficient would be 0 (norm \(-2\)) or force an endpoint outside the range. Every sufficiently high power of the positive unit times a given positive generator has both coefficients positive: the conjugate shrinks to zero in absolute value. These positive pairs give one forward ray since \((s,n)\mapsto(s+2n,s+n)\) preserves positivity and backward iteration eventually leaves the cone. Norm odd gives \(s\) odd and hence mixed parity in \((m,n)\). This proves both the count and the absence of extra rays.

Consequently all three primes \(2,3,11\) from the odd-square bracket fail as primitive leg differences. For the [proper-divisor identity](DIVISOR-TYPES.md), a balanced type \(p,p^3,p^2qr\) can have respectively 4, 4, or 16 scalar rays, provided every constituent prime splits. A nonsplit or even constituent gives none. The squarefree support selects the number of rays; it cannot distinguish the two sheets or guarantee an individual Hamilton count.

## A sharp bound converts a height window to infinite sectors

At depth \(k\ge1\), consider the new rays born there. Their minimum gap is
\[
                            2(k+1)^2-1.
\]
The parent has depth \(k-1\), with \(M=m+n\ge3+2(k-1)=2k+1\), because the three branches increase \(M\) by \(2n,2m,2(m-n)\), respectively. Its leg sum \(M^2-2n^2\) is at least \((M^2+2M-1)/2\), since \(n\le(M-1)/2\). Equality at the claimed minimum comes from the parent reached by \(k-1\) digits 2: it has \((m,n)=(k+1,k)\). In particular a bank realizing every allowed odd target through \(2\cdot3^L-1\), with \(L\ge4\), realizes all allowed labels in every gap \(g<2(L+1)^2-1\). All starting labels there have depth at most \(L-1\) and are at most \(2\cdot3^L-1\); multiplication by 3 extends their rays, with 63 handling the two holes.

**FINITE-EXACT input, PROVED infinite extension.** For \(g\le1200\) there are 124 distinct gaps and 594 scalar seeds. Of these, 569 lie in the established finite tournament window through \(10^9\), seventeen are ordinal products of counts in that window, and 7 is the single excluded starting seed. The replay supplies literal tournaments for each of the other seven starting counts:

| Gap | Seed counts outside the finite window and its factor test |
|---:|---|
| 721 | 1,549,681,951 and 2,324,522,929 |
| 799 | 6,973,568,801 |
| 881 | 20,920,706,405 |
| 1,057 | 188,286,357,653 |
| 1,151 | 376,572,715,303 and 564,859,072,961 |

The window input uses the finite-exact root-change bank, whose 29 missing word counts are all at most 2647, together with the published small complete interval apart from 7 and 21. [HYP-9494(D), finite range](../../hypotheses/HYP-9494-root-change-carmichael-hierarchy-two-root-change-coverage-c1-depth-three-adic-factor.md) · [Root-change census, §§5–6](../root_change_carmichael_chamber_20261107/PROOF.md) · [Chan and Pak, §6.5, Royle's small interval](https://doi.org/10.1016/j.tcs.2024.114776)

To enumerate every ray beginning at positive depth, the script enumerates all primitive parents with \(m^2+2mn-n^2\le1200\), using \(m\le\lfloor\sqrt{1200}\rfloor\). It takes both outer children and both choices of leading digit. Each new large witness is a small recorded core with a single apex over transitive flanks. The check counts paths directly on its adjacency rows by splitting at the apex; it also checks that 11,13,17,63 have the claimed ordinary tournament counts. From this directory run `python3 -B berggren.py`. Unbounded gaps remain open.

## How many ray seeds must an absolute inverse address?

**PROVED.** Let \(r(g)\) count the scalar tripling rays at leg difference \(g\). For \(\operatorname{Re}(s)>1\), the classification gives
\[
 \sum_{g\ge1}\frac{r(g)}{g^s}
 =2\prod_{p\equiv1,7\ (8)}\frac{1+p^{-s}}{1-p^{-s}}
 =\frac{2\zeta(s)L(s,\chi_8)}{\zeta(2s)(1+2^{-s})},
\]
where \(\chi_8\) is 1 on 1,7 and \(-1\) on 3,5 modulo eight, and vanishes on even integers. The Euler product follows by summing the two choices at each positive prime exponent. In the quotient, the factors at the inert primes cancel and the displayed correction cancels the factor at 2.

There is an independent elementary count by area:
\[
 \sum_{g\le D}r(g)
 =\frac{4\sqrt2\log(1+\sqrt2)}{\pi^2}D+O(\sqrt D\log(2D)).
\]
Besides the two root labels, count four starting labels per coprime mixed-parity parent with \(m>n>0\) and \(m^2+2mn-n^2\le D\). Setting \(u=n/m\), the region area is \(D\int_0^1(1+2u-u^2)^{-1}du/2=D\log(1+\sqrt2)/(2\sqrt2)\). Its boundary length is \(O(\sqrt D)\). Half the lattice points have opposite parity, with a boundary error \(O(\sqrt D+1)\) by counting the two permitted residue classes modulo two. On these pairs use Möbius inversion only over odd common divisors \(d\le\sqrt D\). The scaled boundary errors total \(O(\sqrt D\log(2D))\); the area terms have density \(\tfrac12\prod_{p>2}(1-p^{-2})=4/\pi^2\), with the truncated series contributing \(O(\sqrt D)\). The resulting linear main term counts distinct scalar seeds. It does not require or imply a density for exceptional Hamilton path counts.
