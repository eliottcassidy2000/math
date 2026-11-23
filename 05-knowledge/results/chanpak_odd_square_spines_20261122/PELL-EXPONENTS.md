# The primes 17 and 29 select the sharp Pell costs at every odd exponent

**PROVED with finite-exact phase inputs.** Take every odd \(e\ge3\) and target \(N_e=2^{e-1}(2^e-1)-1\). Require a primitive positive normal \((a,b)\) to give minimum projected support simultaneously on Farey spines four and eight. The following rule gives the negative defect \(D=2a^2-b^2\) nearest zero and the least positive defect. Exponent classes with incompatible spines have no such normal at any slope.

| \(e\bmod28\) | Further condition | Signed extremal defects | Number of odd classes modulo 112 |
|---|---|---|---:|
| 7, 11, 13, 23 | All | No common minimizing direction | 16 |
| 3, 9, 17, 25, 27 | All | \((-62,62)\) | 20 |
| 1, 5, 15, 19, 21 | \(e\equiv1\bmod8\) | \((-2,2)\) | 6 |
| 1, 5, 15, 19, 21 | \(e\not\equiv1\bmod8\) | \((-2,62)\) | 14 |

The corresponding one-sided lower limits of \(b^2|a/b-1/\sqrt2|\) are obtained by dividing the defect magnitudes by \(2\sqrt2\). At exponent 31 both costs are \(31/\sqrt2\). At exponent 61 only the upper cost is \(31/\sqrt2\). A common optimal direction concerns two modular projections. These conditions do not assert or exclude the physical tournament flip.

## Why the prime 29 decides compatibility

At 29 both parity correction coefficients vanish on each spine. With the notation of the [projected hyperbola](SQUARE-PRIMES.md), its equation is \(uv=C\) with \(C=q(N_e-k)+ct\bmod29\). For \(C\ne0\) the smallest image of a non-axis direction has \((29-1)/2=14\) residues and occurs exactly when \(\chi_{29}(b/a)=-\chi_{29}(C)\). An axis uses 28 residues. For \(C=0\) the equation is a union of two axes and every nonzero normal attains all 29 residues. This local factor is shared by all four physical parity masks, so it factors out of the size of their union. The following small table uses \(2^{28}=1\bmod29\). A star means all directions minimize.

| \(e\bmod28\) | Required character at spine 4 | Required character at spine 8 |
|---|:---:|:---:|
| 1, 5, 15, 19, 21 | \(-\) | \(-\) |
| 3, 27 | \(+\) | \(+\) |
| 9, 25 | \(*\) | \(+\) |
| 17 | \(+\) | \(*\) |
| 7, 11, 23 | \(+\) | \(-\) |
| 13 | \(-\) | \(+\) |

For all 40 compatible classes, the necessary marginal tests at 2 and 5 still force \(a\) odd, \(10\mid b\) and \(5\nmid a\). The residues \(D\bmod16\in\{2,14\}\) and \(D\bmod5\in\{2,3\}\) leave only \(\pm2,\pm18\) at magnitude below 62. Either sign of 18 forces \(3\mid a,b\). Thus only the two unit Pell orbits need exhaustion below the exhibited magnitude.

The unique nonnegative reduced seeds for \(D=2,-2\) are respectively \((1,0),(1,2)\). Indeed \(T^{-1}(a,b)=(3a-2b,3b-4a)\) stays nonnegative for \(4a/3\le b\le3a/2\), and there it decreases \(a+b\) strictly for positive inputs. For \(D=2\) a reduced seed has \(b<4a/3\), hence \(a^2<9\); for \(D=-2\) it has \(b>3a/2\), hence \(a^2<8\). Checking these finite bounds and the axes gives the two seeds. Under \(T\) their projective orbits at 29 are the same set of slopes \(0,11,27,2,18\). Every nonzero one is a nonsquare. A common positive character at 29 therefore rules out both signs of 2. When the common character is negative, the additional obstruction at 17 described below distinguishes the signs.

For sufficiency the replay works with the original four-parity union on each spine. For each of the 40 classes it explicitly finds an optimal normal on the first viable orbit in each of \((1,2),(1,8)\) below and \((1,0),(7,6)\) above. The second seeds have defects \(-62,+62\). One full period of 480 projective steps suffices. The script prints each selected seed and exponent of \(T\), checks positivity and primitivity of the actual integer vector, and checks it against both joint projections. Iterating \(T^{480}\) supplies the infinite sequence for the associated lower limit. Thus the conditions at 257 and the shared parity choices are retained even though the final classification uses only \(e\bmod28\) and \(e\bmod8\).

## The thirteen phases selected by the cycle direction

**PROVED with finite-exact phase inputs.** Let \(e\ge3\) be odd, \(N_e=2^{e-1}(2^e-1)-1\), and restrict \(e\) to the thirteen classes in which \((7,10)\) minimizes the projected residue support on both Farey spines. Among all primitive positive normals \((a,b)\) simultaneously minimizing both projections, the least magnitude of a negative defect \(D=2a^2-b^2\) is always 2. The least positive defect is 2 when \(e\equiv1\bmod8\) and 62 otherwise. [Classification of the thirteen classes](SQUARE-PRIMES.md)

| \(e\bmod112\) | Closest negative defect | Least positive defect | Positive witness |
|---|---:|---:|---|
| 1, 33, 49 | \(-2\) | 2 | \((99,140)\) |
| 43, 75 | \(-2\) | 62 | \((191,270)\) |
| 5, 21, 29, 61, 71, 77, 85, 103 | \(-2\) | 62 | \((293465818991,415023341310)\) |

In particular the one-sided constants for approximation of \(1/\sqrt2\) are both \(1/\sqrt2\) in the first row. In the other rows the lower-side constant is \(1/\sqrt2\) and the upper-side constant is \(31/\sqrt2\). The appearance of 31 is governed by the exponent modulo eight within this bank. It occurs at ten classes, including \(e=61\), and does not require \(2^e-1\) to be prime.

## Why the two sides split

For each exponent, at 2 and 5 every simultaneous minimizer has normal \([1:0]\). Thus \(a\) is odd, \(10\mid b\), and \(5\nmid a\). The only possibilities for \(0<D<62\) are 2 and 18, by residues \(2,14\bmod16\) and \(2,3\bmod5\). If \(D=18\), reducing modulo three forces both coordinates to be multiples of three. Below zero the witness \((7,10)\) gives \(D=-2\), the smallest possible magnitude because \(D\) is even.

The only possible upper improvement is now the primitive norm-one equation \(a^2-2(b/2)^2=1\). By the [Pell reduction](PELL-DIRECTIONS.md), all its positive solutions are in the \(T\)-orbit of \((1,0)\), where \(T(a,b)=(3a+2b,4a+3b)\). Its projective slopes \(b/a\bmod17\) are \(0,7,\infty,10\).

At 17 the marginal minimizing directions for spine four are all 16 finite nonzero slopes when \(e=1\bmod8\), and only the eight nonzero squares otherwise. The source of this change is visible directly in the two parities of the right address. For \(D_j=q(N_e-k-e_{\rm tag}j)+ct\bmod17\), where \(e_{\rm tag}=1740\) is the coefficient of the right parity in the spine-four chart, the discriminant test is that \((qS+ac+bt)^2-4ab D_j\) is a square or zero. The left parity coefficient vanishes modulo 17. Counting the union over \(j=0,1\) gives the following sizes. The axis case covers either vanishing weight.

| \(e\bmod8\) | \((D_0,D_1)\) | \(b/a\) a nonzero square | \(b/a\) a nonsquare | Axis |
|---:|---|---:|---:|---:|
| 1 | (10,4) | 13 | 13 | 16 |
| 3 | (3,14) | 12 | 13 | 16 |
| 5 | (12,6) | 12 | 13 | 16 |
| 7 | (11,5) | 12 | 13 | 16 |

For fixed nonzero weights, rescaling the discriminant variable identifies all directions with the same quadratic character of \(ab\), equal to that of \(b/a\). Checking \(b/a=1,3\) and either axis on the 17-element field therefore proves the local table. Full joint minimization in the replay confirms that the surviving 17-directions are exactly the marginals indicated. Since 7 and 10 are nonsquares modulo 17, the norm-one orbit is disjoint from the allowed slopes for the ten asymmetric exponent classes.

Checking the original joint projections at all five primes gives an upper witness \(T^3(1,0)=(99,140)\) in the other three classes. For the asymmetric classes \(T^2(7,6)\) and \(T^{14}(7,6)\) give the two listed witnesses of primitive norm 31. All checked seeds are primitive and \(T\) is unimodular. The projective orders at \(2,5,17,29,257\) divide 480, so successive powers of \(T^{480}\) from each witness remain optimal. They converge in slope to \(1/\sqrt2\) with unbounded denominator. The identity \(b^2\lvert a/b-1/\sqrt2\rvert=\lvert D\rvert/[2(a/b+1/\sqrt2)]\) proves both one-sided lower limits.

This classification answers an approximation question about the search directions. It leaves the absolute digital height and the existence of the exponent-61 flip open. Run the [period and optimization replay](pell_exponents.py) from this directory with the command below.

    python3 -B pell_exponents.py
