# Why a single growing transitive module cannot follow either perfect-number neighbor

Let \(H(w)\) count ordinary Hamilton paths after replacing a vertex of a fixed finite tournament by a transitive tournament of size \(w\ge1\). For integers \(p\ge1\) write \(E_p=2^{p-1}(2^p-1)\). This formula gives an even perfect number whenever \(2^p-1\) is prime. Fix a shift \(d\ge1\) and either sign throughout.

**PROVED.** The widths at which \(H(w)=E_{w+d}\pm1\) contain no four-term arithmetic progression of positive common difference. They cannot even contain three consecutive integers starting at \(w\ge2\), or three consecutive odd or even widths starting at \(w\ge3\). In particular this holds separately for either orientation of any one-arc edit outside the growing module.

### The integral readings at consecutive widths

Write \(C_j\) for the number of directed symbolic spanning paths that visit the growing module \(j\) times. Every fixed vertex appears once. The \(j-1\) intervening excursions are nonempty, have the required initial and terminal contacts, and have disjoint vertex sets. Permuting them acts freely. Hence \((j-1)!\mid C_j\). Filling each run in its unique transitive order gives
\[
 H(w)=\sum_{j\ge1} C_j\,j!\,{w\brace j}=\sum_{\ell\ge1} A_\ell\,\ell^w,\qquad
 A_\ell=\sum_{j\ge\ell}(-1)^{j-\ell}\binom j\ell C_j.
\]
In particular \((\ell-1)!\mid A_\ell\). These finite sums and the excursion action are the [homogeneous-reader factorial theorem](../chanpak_carry_admission_20261115/RICH_CORE_ADMISSION.md).

For any \(k\ge1\) let \(s(k,i)\) denote the signed Stirling numbers of the first kind. For every \(m\ge0\),
\[
              k!(k-1)!\ \mid\ \sum_{i=1}^k s(k,i)H(m+i).
\]
Indeed the multiplier of \(A_\ell\ell^m\) on the right is \(k!\binom{\ell}{k}\), and every nonzero term has \(\ell\ge k\). At \(k=4\) this gives
\[
144\mid H(m+4)-6H(m+3)+11H(m+2)-6H(m+1).
\]
The target function is \(F(w)=A4^w-B2^w\pm1\), where \(A=2^{2d-1}\), \(B=2^{d-1}\). The same expression for \(F\) is \(24 A4^m\), which has just one factor of three. This already excludes four consecutive matches at any width. It also pinpoints the coefficient failure for the proposed infinite ray. Uniqueness of finite power sums forces \(C_j=0\) for \(j>4\) and \(C_4=A\), contrary to \(6\mid C_4\). Uniqueness applies on every infinite subsequence by dividing by the largest base with a nonzero difference.

The modulus \(k!(k-1)!\) is sharp over strong sources. Use \(k-1\) pairs \(a_i,b_i\) beside a port \(v\), with \(v\to a_i\to b_i\to v\) and \(b_j\to a_i\) for \(j\ne i\). Orient the two internal groups arbitrarily. Every \(k\)-visit symbolic path must use all \(2(k-1)\) fixed vertices as \(k-1\) two-vertex excursions, giving precisely the \((k-1)!\) permutations of the matched pairs. Larger visit counts are impossible. Thus its \(k\)-point expression at \(m=0\) equals \(k!(k-1)!\). In fact the integer span of all sample vectors \((H(1),\ldots,H(k))\) is exactly the lattice given by these \(k\) consecutive-width divisibilities at \(m=0\). Stirling inversion turns the normalized left sides into \(C_j/(j-1)!\). The displayed sources for \(j=1,\ldots,k\) give a triangular matrix with unit diagonal in these coordinates. This is a statement about integral linear combinations of sources. [Physical triangular probes for the factorial interface](../chanpak_symmetry_interfaces_20261115/OBSERVABLES.md)

### Three tests at the same parity

For a positive step \(s\), put \(\nu=(1+v_3(s))(1+[\,2\mid s\,])\). Whenever \(a\ge\nu+1\),
\[
 3^{\nu+1}\mid H(a+2s)-(1+2^s)H(a+s)+2^sH(a).
\]
The multiplier of \(A_\ell\) is \(\ell^a(\ell^s-1)(\ell^s-2^s)\). Bases 1 and 2 disappear. A multiple of three contributes \(3^{\nu+1}\) from \(\ell^a\). At every remaining base \(\ell\ge4\) there is a three in \(A_\ell\). For odd \(s\), choose whichever of \(1,2\) agrees with \(\ell\) modulo three. The corresponding difference of \(s\)-th powers is divisible by \(3^{1+v_3(s)}\) by lifting the exponent. For even \(s\) both differences have this divisibility, by taking the \((s/2)\)-th powers of squares. At \(\ell=4\) the product has valuation exactly \(\nu\), since \(4-1=3\), and for even \(s\) also \(v_3(2^s-1)=1+v_3(s)\). Thus \(F\) fails the displayed test.

At step one the modulus is 9. At step two the concrete law is
\[
                  27\mid H(a+4)-5H(a+2)+4H(a)\qquad(a\ge3).
\]
The perfect-number reading gives \(180 A4^a\), divisible by 9 exactly. A source proposed for the upper counts at exponents \(57,59,61\), for example, fails if its first width is at least three.

For every \(s\) the exponent \(\nu+1\) of the universal three-adic divisor is sharp. Use the three-pair source just described. Its coefficient \(A_4=6\) has valuation one, and it has no higher bases. Whenever \(a\ge\nu+2\), the base-three contribution vanishes modulo \(3^{\nu+2}\), so the three-term expression has valuation exactly \(\nu+1\).

For arbitrary \(a\ge1\), add the factor \((\ell^s-3^s)\) to the multiplier. The resulting linear combination uses the four widths \(a,a+s,a+2s,a+3s\). Base three disappears. At a larger multiple of three, \(\ell^a(\ell^s-3^s)\) contributes at least \(a+s\ge\nu\) powers of three, and \(3\mid A_\ell\). The proof for unit bases is unchanged, and the extra factor at base four is a unit. Hence no four terms of an arithmetic progression can agree with \(F\), even near width one. This statement concerns a common source over several exponents. A hit at one isolated exponent remains possible.

### Which supports can carry synchronized growth

**PROVED.** Let any positive finite number of distinct marked vertices be replaced by transitive tournaments with respective widths \(w_i=n+d_i\), where the offsets are fixed integers. At sufficiently large \(n\) all widths are positive. If the ordinary height agrees with \(E_{n+d}\pm1\) at infinitely many integers \(n\), write \(C_{4i}\) for the number of symbolic paths with four visits to \(i\) and one to all other ports. Write \(C_{22ij}\) for two visits to \(i,j\) and one elsewhere. Then
\[
 \sum_i C_{4i}4^{d_i}+\sum_{i<j}C_{22ij}2^{d_i+d_j}=2^{2d-1},
 \qquad
 \sum_{i<j}(-1)^{d_i+d_j}C_{22ij}\equiv2\pmod3 .
\]
The congruence takes place in \(\mathbb Z[1/2]\) reduced modulo three. At least two ports must be revisited in the same spanning path and the signed count of these joint returns must survive modulo three. Both orientations of a putative physical seam have to pass with the same right-hand side.

First the marked vertices must induce a transitive tournament. Otherwise, in any spanning path in the original source, insert repeated copies of a directed cycle of marks at an occurrence of one of its vertices, returning to that vertex each time. A fixed resulting symbolic word with visit vector \(\mathbf r\) has \(\prod_i r_i!\,{n+d_i\brace r_i}\) fillings. Its exponential growth base \(\prod_i r_i\) can exceed four, contradicting agreement along an infinite sequence. With acyclic marks every segment between fixed vertices visits the marks in increasing order, giving finite symbolic support. All top growth contributions are positive. The maximum \(\prod_i r_i\) must therefore be four, attained only at a single quadruple or two doubles. For any port \(i\) each of its excursions contains a fixed vertex since the induced marked tournament is acyclic. These fixed vertices distinguish excursions even when other port labels repeat. Thus \((r_i-1)!\) divides the number \(C_{\mathbf r}\) at a fixed visit vector. In particular \(6\mid C_{4i}\), proving the congruence after equating the leading coefficients. This reduces the two-port search to configurations with genuine joint returns. [Prior cap-two source census](../chanpak_absolute_liouville_leech_20261122/MODULE-OFFSETS.md)

### A four-width test without a run cap

For acyclic marked ports the balanced-return test has a finite version. Choose \(m\) such that all \(m+d_i\ge1\), and form
\[
 L=\bigl(H(m+4)-6H(m+3)+11H(m+2)-6H(m+1)\bigr)/24 .
\]
For \(D\) a subset of the ports let \(C_D\) count symbolic paths visiting precisely the ports of \(D\) twice and visiting each other port once. Then
\[
 L\equiv\sum_D C_D\sum_{S\subseteq D}
 (-1)^{\sum_{i\in S}(m+d_i)}\,\left\lfloor |S|/2\right\rfloor \pmod3 .
\]
The formula holds with arbitrarily large finite visit bounds. Expanding each surjection count into powers gives a product base \(b\) multiplied by \(\binom b4\) after the division by 24, so each term is integral. A port visited at least four times forces a factor of three in \(C_{\mathbf r}\). For a three-visit port the choices of bases 1 and 2 carry \(\binom31=\binom32=3\), and base 3 carries \(3^{m+d_i}\). The surviving two-visit factors read \(2^{m+d_i}2^j-2\) at \(n=m+j\). Modulo three, picking their exponential terms on \(S\) gives the sign displayed, and \(-2=1\bmod3\) elsewhere. The remaining factor is \(\binom{2^{|S|}}4\), congruent to \(\lfloor|S|/2\rfloor\bmod3\). Indeed \(2^j\bmod9\) cycles through \(1,2,4,8,7,5\) and extracting the coefficient of \(x^4\) in \((1+x)^b\bmod3\) gives \((b\bmod3)\lfloor b/3\rfloor\bmod3\). Thus four consecutive perfect matches already require joint returns with a nonzero signed weight, before any asymptotic hypothesis.

For a sharp distinction, take the ordinal sum of two directed triangles and widen one vertex in each to widths \(n,n+1\). A single triangle has three one-visit paths and one two-visit path, so the combined height is \((2^n+1)(2^{n+1}+1)=2\cdot4^n+3\cdot2^n+1\). The unique double-return choice gives the required leading coefficient and congruence at \(d=1\). The absolute discrepancy from \(E_{n+1}+1\) is \(4\cdot2^n\), entirely in the next coefficient.

### Where the two returning ports must lie

**PROVED.** In a synchronized upper family valid at infinitely many exponents, all growing ports capable of a return belong to the same strong component of the fixed source. Indeed the component heights multiply. If two separate components support returns, their leading bases are integers at least two, and their product across the source is four. Thus each of the two has maximum visit product two and its height is a nonconstant affine polynomial in \(X=2^n\). The remaining components have constant height. The discriminant of the upper polynomial \(2^{2d-1}X^2-2^{d-1}X+1\) is \(-7\cdot2^{2d-2}\), making this factorization impossible. This applies the [upper irreducibility theorem](../chanpak_cech_perfect_seams_20261122/NORM-AND-PORTS.md) to the visit bound established above.

For the lower family the same split forces the two component heights to equal \(2^{n+d}+1\) and \(2^{n+d-1}-1\), in either order. Each component affine has integer constant term, because each one-visit filling is one and the constant term of a two-visit filling is \(-2\), even with shifted widths. The product of their constant terms and all fixed heights is \(-1\). Thus the fixed heights equal one and the component constants are \(1,-1\). Rational factorization of \(E_{n+d}-1=(2^{n+d}+1)(2^{n+d-1}-1)\) with positive leading coefficients gives the stated normalization. In particular the required joint returns on the upper side of a physical seam have to occur inside one strong component.
