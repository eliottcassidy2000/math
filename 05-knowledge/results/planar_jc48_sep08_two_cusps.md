# A literal sextic with two co-projected ordinary cusps

**Status: PROVED GEOMETRY AND THREE-MERIDIAN BOUND + FINITE-EXACT +
INDEPENDENTLY AUDITED.** The [independent audit](planar_jc48_sep08_two_cusps_audit.md)
accepts the exact geometry, positive-meridian bound, and conditional Euler
ledger. The full complement group and whole-Keller exclusion remain OPEN.
This is one parameter value, not a classification of the surrounding family.

## 1. Object and inheritance

The bounded target is
\[
 U=t^4-2t^2,\qquad V=t^6-\frac32t^4+\frac13t^3-t.
\tag{1}
\]
Its image is a birationally normalized rational sextic with two ordinary
\((2,3)\) cusps, four ordinary nodes, and one \((2,9)\) branch at
infinity. The two cusps have the same U-coordinate but distinct images.

The closest current mechanism is the
[higher-cusp rational braid bundle](planar_jc48_sep07_higher_braid.md):
monic vertical quartics and actual positive meridian relations can
outperform an infinity-only obstruction. Its one-cusp sheet-count
consumer is not transferred to this two-cusp object. A more distant
proved sidecar is
[THM-3844, two-cusp polynomial branch quadratic-resolvent design gate](../../01-canon/theorems/THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate.md).
That theorem concerns a different quartic, with normalization
\((t^2(2t-3),t^3(3t-4))\), and proves that a global three-class and its
monogenic cubic completion still do not supply a Keller source. No class
group or completion theorem is transported from it here.

The four-concept board is: full collision branches, actual cusp access
paths, monic-quartic meridians, and the distinction between branch
geometry and an admissible source. The canonical hostile is the
monogenic completion in THM-3844. The corrected near miss is losing the
\(s+t=0\) collision branch or assuming separate projection values for
the cusps. The least-used sidecar is the marked access-path information
needed to combine local cusp relations globally.

## 2. Exact finite singularity inventory

The derivatives satisfy
\[
 U'=4t(t^2-1),\quad V'=(t^2-1)(6t^3+1),\quad
 \gcd(U',V')=t^2-1.
\tag{2}
\]
The only critical parameters are \(t=1,-1\), with images
\[
 (-1,-7/6),\qquad(-1,1/6).
\]
At these parameters \(U''=8\), and
\(U''V'''-V''U'''\) equals 352 and -416, respectively. Thus both
critical branches are ordinary cusps. Also
\[
 \gcd(U+1,V+7/6)=(t-1)^2,\quad
 \gcd(U+1,V-1/6)=(t+1)^2,
\tag{3}
\]
so neither cusp image has an extra parameter preimage.

Write \(N=(U(s)-U(t))/(s-t)\), \(M=(V(s)-V(t))/(s-t)\), and
\(p=s+t,q=st\). The complete equations are
\[
 N=p(p^2-2q-2),\qquad
 6M=6p^5-24p^3q-9p^3+2p^2+18pq^2+18pq-2q-6.
\tag{4}
\]
There are two branches of the first equation, and both must be kept:

* If \(p=0\), then \(M=-(q+3)/3\), so \(q=-3\).
  This gives the pair \(t=\pm\sqrt3\), whose common image is
  \((3,27/2)\).
* If \(q=(p^2-2)/2\), then
  \[
  M=-\frac{(p-2)(p+2)(3p^3-2)}{12}.
  \tag{5}
  \]
  The factors \(p=\pm2\) are the two diagonal cusp parameters.
  The three roots of \(3p^3-2\) give three further unordered pairs.
  Their discriminant \((s-t)^2=4-p^2\) is nonzero, and none has
  \(p=0\). The intersection \(p=0,q=-1\) of the two branches
  of N is not a zero of M.

Thus the collision list contains exactly four off-diagonal pairs. The
full ordered-pair resultant retains both cusp multiplicities:
\[
 \operatorname{Res}_s(N,M)=
 -\frac{(t-1)^2(t+1)^2(t^2-3)}{27}
 (18t^6-54t^4+6t^3+54t^2-18t-17).
\tag{6}
\]

For \(T=U'(s)V'(t)-V'(s)U'(t)\), the exact lexicographic Gröbner
basis of \((N,M,T)\), in variables \((s,t)\), is
\[
 \{s+t^3-2t,\ (t^2-1)^2\}.
\tag{7}
\]
Its set-theoretic support is just the two diagonal cusp pairs. Hence all
four off-diagonal intersections are transverse. There are no triple
images: for a putative image \((A,B)\),
\[
 R=\operatorname{rem}_t(V-B,U-A)
   =t^3/3+(A+1)t^2-t+A/2-B.
\tag{8}
\]
Its leading coefficient is the fixed nonzero number \(1/3\), and the
coefficients of \(\operatorname{rem}_t(U-A,R)\) generate the unit
ideal in \(\mathbb Q[A,B]\). Three distinct common roots would force
the cubic R to divide \(U-A\), a contradiction. Equations (2)--(8)
therefore exhaust all finite singularities and show that the four
pairs have distinct ordinary-node images.

The parametrization is finite, since t satisfies the monic equation
\(t^4-2t^2-U=0\) over the coordinate ring of the image. The finite
collision list makes it generically injective and hence birational;
its normalization is \(\mathbb A^1\).

## 3. Infinity and the actual vertical equation

The degree-six homogenization has no basepoint and is birational onto
its image, proving that the image is a sextic. At the unique infinity
parameter use
\[
 X=U(1/z)/V(1/z),\qquad Z=1/V(1/z).
\]
Then \(X\sim z^2\), \(Z\sim z^6\), and
\[
 [z^7](Z-X^3)=0,\quad
 \lim_{z\to0}(Z-X^3)/X^4=3,\quad
 [z^9](Z-X^3-3X^4)=2/3.
\tag{9}
\]
The infinity branch is therefore \((2,9)\), with line contact six.
The genus check is \(1+1+4+4=10\), accounting for the two affine
cusps, infinity cusp, and four nodes.

The source stores the complete monic quartic
\[
 F(u,v)=\operatorname{Res}_t(U-u,V-v),\qquad\deg_vF=4,
\tag{10}
\]
recomputes that resultant, and checks \(F(U,V)=0\). Birationality and
\([\mathbb C(t):\mathbb C(U)]=4\) ensure this is the actual minimal
equation. Its discriminant is
\[
 -\frac{256}{531441}u(u-3)^2(u+1)^6
 (324u^3+972u^2+1080u+289)^2.
\tag{11}
\]
The cubic node-projection factor has discriminant -59592250800 and
resultant 868900175 against \(u(u-3)(u+1)\). Thus the four nodes
have separate projection values, none equal to a cusp or smooth
critical value. The fibres at the special values are
\[
\begin{aligned}
 F(-1,v)&=(6v-1)^2(6v+7)^2/1296,\\
 F(0,v)&=v^2(9v^2-36v+34)/9,\\
 F(3,v)&=(2v-27)^2(36v^2+180v+289)/144.
\end{aligned}\tag{12}
\]
In particular the two cusp targets co-project at \(u=-1\), but are
distinct points in that fibre. The value \(u=0\) has one smooth
vertical fold, at t=0, because \(U''(0)=-4\) and \(V'(0)=-1\).

## 4. What this pays for monodromy, and what remains

The monic equation supplies the same outside-root section and fibre
surjection \(F_4\twoheadrightarrow\pi_1(\mathbb C^2\setminus C)\)
as in the audited braid suppliers. A small loop around \(u=0\),
accessed through the regular base, is a conjugate of one positive half
twist: its only collision is the smooth fold in (12). After an adapted
positive-meridian basis change, its fixedness identifies two of the
four meridians. Consequently **the group is generated by at most three
positive meridians**. The independent audit proves this analytic consequence
from the actual fold and the positive adapted basis; no numerical braid
word is needed.

The local braid at \(u=-1\) is, in a basis adapted to disjoint cusp
neighborhoods, a product of two disjoint cubes of positive half twists.
The four node values give conjugates of squares of positive half
twists. The unresolved datum is their simultaneous access-path gauge
relative to the three-meridian basis. The present source does not claim
those conjugators, a full group presentation, two-meridian generation,
or a transitive representation.

The previous one-cusp Euler ledger cannot simply be reused. Suppose C is
the whole irreducible nonproperness curve of a nonautomorphic polynomial
Keller map. For two unibranch cusps, actual retained count a on its whole
smooth stratum, generic covering degree d,
actual cusp counts \(n_1,n_2\), and deleted-complement node overlaps
\(\omega_p\), the same stratification instead gives
\[
 1=-a+n_1+n_2+\sum_p\omega_p.
\tag{13}
\]
Indeed the curve's smooth stratum has Euler characteristic
\(1-2N-2=-1-2N\), its complement has Euler characteristic N, and
each node fibre contributes \(2a-d+\omega_p\). Hence the two-meridian
one-cusp exclusion is not paid by the geometry in this file. The
cheapest next test is a certified common system of loops for (11), or
a global source obstruction respecting (13).

Two parameter hostiles delimit this literal control. In the ambient
family \(V_\lambda=t^6-3t^4/2+\lambda(t^3/3-t)\), the cusp jet
determinants are \(32(12-\lambda)\) and \(-32(12+\lambda)\).
The ordinary-cusp test fails at \(\lambda=\pm12\). At \(\lambda=0\)
both U and V are even and the parametrization loses birationality.
These exact controls do not classify any other parameter.

## 5. Reproduction

The explicit universe is the single curve (1), its full pair equations
including the zero-sum branch, its exact tangent/triple ideals, its
infinity jets and full discriminant, plus the three named parameter
hostiles. No floating-point root tracking or permutation census occurs.

Source: [planar_jc48_sep08_two_cusps.py](../../04-computation/planar_jc48_sep08_two_cusps.py).
Output: [planar_jc48_sep08_two_cusps.out](planar_jc48_sep08_two_cusps.out).

```bash
python3 04-computation/planar_jc48_sep08_two_cusps.py
python3 -O 04-computation/planar_jc48_sep08_two_cusps.py
```

Both runs pass the same 42 always-active exact gates and reproduce the
505-byte output byte for byte. Frozen SHA-256 pins:

* Source: `d3c1bffe31c3147557337066dcce1ee4dc97d19cd2f77832244fb1b0317f4af0`.
* Output: `c1981a35fe9293387f60c356c3d4bad70ef58e25ef65ad1c35c7779c122831a5`.

The independent analytic and source audit passes with both fresh replays
matching the frozen bytes. The full global group and any Keller realization
remain OPEN. Source and output are frozen.
