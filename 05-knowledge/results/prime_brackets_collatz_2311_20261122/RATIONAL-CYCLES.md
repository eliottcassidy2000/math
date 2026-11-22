# Independent residues at 2, 3 and 11 on signed rational Collatz cycles

**PROVED.** Fix either real sign, a residue modulo \(2^b\), a unit residue modulo \(3^s\), and a ball \(z+11^u\mathbb Z_{11}\) in \(\mathbb Q_{11}\). There are genuine rational shortcut-Collatz cycles of that sign, of arbitrarily large **minimal** period, whose initial point meets all three prescriptions. All points have odd reduced denominator. One can prescribe the exact 11-adic valuation of the **unreduced clock denominator** as any positive integer at least the required reduced denominator valuation. This is a statement about rational periodic points, not new positive or negative *integer* cycles.

Use the shortcut branches \(C_0(x)=x/2\) and \(C_1(x)=(3x+1)/2\) on \(\mathbb Z_2\), according to the parity of \(x\). A chronological binary word with length \(n\), weight \(k\), and numerator \(B\) has affine map \(x\mapsto (3^k x+B)/2^n\), with \(B>0\) when \(k>0\). Its fixedpoint \(x=B/(2^n-3^k)\) is an actual rational cycle. To see that its branches are correct, the affine inverses \(y\mapsto2y\), \(y\mapsto(2y-1)/3\) send \(\mathbb Z_2\) into the required even or odd class, respectively. Their backwards composition contracts by \(2^{-n}\) in 2-adic norm and has that rational fixedpoint, whose orbit therefore reads the given word. Cyclic shifts give fixedpoints with the same denominator sign, so the whole orbit has the stated real sign.

## Clock resonance and the base-3 exception at 11

The element 2 has order 10 modulo 11 and \(3\equiv2^8\). The first resonance is consequently \(n\equiv8k\bmod10\). Here \(n\) is even. Since \(2^{10}-1=93\cdot11\) while \(3^5-1=2\cdot11^2\), define
\[
 \theta=\frac{2\log(3^5)}{\log(2^{10})}\in11\mathbb Z_{11},
 \qquad \theta/11\equiv3\pmod{11}.
\]
For resonant \((n,k)\),
\[
 v_{11}(2^n-3^k)=1+v_{11}(n-\theta k).
\]
Indeed, the ratio \(2^n/3^k\equiv1\bmod11\), and taking its fifth power preserves the valuation of the difference from 1. The fifth power is \((2^{10})^{n/2}/(3^5)^k\). On \(1+11\mathbb Z_{11}\) the logarithm preserves this valuation (its linear term has strictly lower valuation than every higher term), and its value is \(\frac12\log(2^{10})(n-\theta k)\). The two log valuations and the normalized congruence for \(\theta\) follow by reducing the first logarithm term at the first nonzero digit. Thus resonance gives valuation one when \(11\nmid n\), at least two exactly when \(11\mid n\), and at least three exactly when \(n/11\equiv3k\bmod11\). For nonresonant words the valuation is zero.

## Thue swaps write a chosen 11-adic translation without changing the clock

For positive sign use base clock \((n_0,k_0)=(6,2)\), and start with the two words \(001010,001100\); for negative sign use \((4,3)\) and words \(0111,1011\). The first slope is \(a=9/64<1\), the second \(a=27/16>1\). They each satisfy \(v_{11}(a-1)=1\). The two maps at depth zero have the same slope and translation difference of valuation zero: changing adjacent \(01\) to \(10\) has affine translation difference \(\pm1/4\) before multiplication by the fixed later slope, which is an 11-adic unit.

From \(A_0,B_0\) form \(A_{j+1}=A_jB_j,\ B_{j+1}=B_jA_j\) as chronological concatenations. They consume the same clock, \(2^j(n_0,k_0)\), at depth \(j\). Their translation difference \(\Delta_j\) satisfies
\[
 v_{11}(\Delta_j)=j.
\]
Indeed \(\Delta_{j+1}=\pm(a^{2^j}-1)\Delta_j\) up to the order chosen for differences, and \(a^{2^j}-1\) has valuation one by LTE. This gives exact 11-adic digits even though the word is written in binary.

To prescribe translation modulo \(11^d\), install ten choices between \(A_j,B_j\) at each depth \(0\le j<d\). Their translations add affinely; the multiplier on each individual difference is the product of slopes of all subsequent words. It is independent of the choices and is \(1\bmod11\) for subsequent matching-clock blocks. Therefore the ten differences at level \(j\), divided by \(11^j\), have one common nonzero residue modulo 11. Toggling \(0,\ldots,10\) of these choices selects the next digit, without changing any previously selected digit. All choices have the same complete binary length and weight, consuming \(10(2^d-1)\) base clocks. Any fixed extra prefix or suffix of matching base clocks preserves this argument; in particular a common tail multiplies these differences by the same 11-adic unit.

## Independent outer guards, and primitive periods

**At 2:** given any residue modulo \(2^b\), take its first \(b\) shortcut parity bits, and force them at the start of the word. That parity prefix uniquely determines the initial residue: two arguments with identical first \(b\) bits have 2-adically integral images after those steps and their difference has been multiplied by \(3^k/2^b\). Therefore their difference is divisible by \(2^b\). Pad the prefix by a fixed number of ones and zeros to equal a whole number of base clocks. This is possible since both \(k_0\) and \(n_0-k_0\) are positive.

**At 3:** any rational cycle with at least one odd step has all cycle points in \(\mathbb Z_3^\times\). Its denominator \(2^n-3^k\) is a 3-adic unit; an odd step lands at a unit and even steps preserve units. More strongly, the fixed tail \(1^s0^e\) forces the initial point (the terminal point of the periodic word) to satisfy
\[
 x\equiv-2^{-e}\pmod{3^s},
\]
because \(C_1^s(y)+1=(3/2)^s(y+1)\) for every \(y\in\mathbb Z_3\). Since 2 generates the units modulo \(3^s\) (its order is \(2\cdot3^{s-1}\) by LTE), choose \(e\) to match any prescribed unit.

**Primitive marker:** reserve a final fixed block of \(m\) base clocks. Write it as
\[
 0\,1^{\,k_0m-s}\,0^{\,(n_0-k_0)m-e-1}\,1^s0^e.
\]
Choose \(m\) large enough that the middle zero run is nonempty and its first run of ones is strictly longer than the prefix block's full length plus \(s+2n_0\). Its total length and weight are \(mn_0,mk_0\). In the remaining body (base words and Thue swaps of them), any one-run has length at most \(2n_0\): each individual base word contains a zero. Thus the long marker is the unique longest cyclic one-run, even if the terminal \(1^s\) touches the initial prefix. A word that is a proper power would have at least two such runs. This forces the word to be primitive, hence the actual rational orbit has minimal period its entire length.

Finally choose total base-clock length \(L\) with \(v_{11}(L)=r-1\), for any desired \(r\ge1\), as large as desired after accounting for the prescribed blocks. Fill missing clocks with a fixed base word, before the final block. By LTE \(v_{11}(1-a^L)=r\), which is also the valuation of the unreduced clock denominator since powers of 2 are 11-units. Take \(r\) large enough that \((1-a^L)z\in\mathbb Z_{11}\); this only requires \(r\ge -v_{11}(z)\) when \(z\ne0\). With \(d\ge\max(0,r+u)\), program the total translation to be \((1-a^L)z\bmod11^d\) using the middle bank; \(d\) is chosen before \(L\), since it depends only on \(r,u\). Then the fixedpoint lies in the prescribed 11-adic ball after division by \(1-a^L\). All three controls coexist at one exact clock and the taper does not change the sign determined by \(a<1\) or \(a>1\).

The structural distinction is exact: 2 reads a finite parity guard, 3 permits a contracting terminal branch with unit-only range, and 11 permits translation steering at a resonant clock. Local compatibility imposes no automatic global integrality: primes outside these places can survive in the reduced denominator.

In particular, every integer seed prime to 3 can be mimicked to arbitrary simultaneous precisions at these three places by signed primitive rational cycles. Even exact agreement with a finite parity prefix, together with all fixed local precisions at 3 and 11, is not an integer-cycle certificate.

Run `python3 rational-lift-check.py` for the eight exact primitive rational cycles and the short-clock valuation controls.
