# A 66-vertex residue lift from the actual exponent-31 edge

**PROVED.** One insertion of six bits into the 60-vertex exponent-31 edge yields a strong single-arc seam whose center matches \(E_{61}=2^{60}(2^{61}-1)\) simultaneously modulo \(4,81,5,7,13\). The same insertion with length \(m=6+12t\) realizes arbitrarily high joint precision at \(3,5,7,13\) for every \(E_r\) with \(r\equiv1\pmod{12}\). Its centers have binary valuation exactly two for \(m>0\). This is a small witness for the original six-label residue test. Absolute exponent-61 placement and finer binary precision require more information. [Exact exponent-31 witnesses](../chanpak_e31_physical_20261122/README.md) · [Independent unequal-flank decoder](README.md)

The absolute search now has two explicit bracketing edges and a finite inverse for both the receiver and its address. The exact lower construction has a two-class neutral component with no upper exit at exponent 61. [Receiver and component results](../chanpak_e61_fibonacci_target_20261122/README.md)

Use the third of the three explicit witnesses. Its deletion data read \(h=1749,c=53334\), so at fixed odd address parity its center has form \(K+cy+hp(y)\). The original right address is \(y_0=42773862359093=16u+5\) with \(u=2673366397443\) odd. Append \(m\) ones to this higher prefix and retain the four-bit suffix 0101:
\[
                  y_m=16((u+1)2^m-1)+5.
\]
The high prefix keeps its one-apex reading \(p(u)\) and its odd parity. The padded reversal of the suffix is 10. Thus the actual center, including \(m=0\), is
\[
 Z_m=E_{31}+S(2^m-1),\qquad
 S=(u+1)(16c+5h)=2304679764206100516.
\]
The single arc in the deletion and its two endpoint corrections are unchanged by the insertion. The starting apex-deletion configuration is strong after adjoining the original suffix. Inserting more letters strictly between the first and last vertices of that transitive suffix preserves strongness: any new vertex is reached via the first old suffix vertex and reaches the last old one.

At \(m=6\), \(E_{61}-Z_6\) has valuations 4, 1, 1, 1 at the four odd primes. Since \(v_3(S)=1\) and \(S\) is a unit at \(5,7,13\), the odd image for \(m=6+12t\) consists of all lifts of the phase of \(Z_6\) modulo \(27,5,7,13\). This uses the exact period argument from the main note. Every \(r\equiv1\pmod{12}\) shares these base residues: \(2^r-2^{61}\) is divisible by 4095 and \(2^r+2^{61}-1\) is divisible by 3, so the factored difference of centers vanishes modulo \(3\cdot4095\). At 2 one has \(v_2(S)=2\) and \(v_2(E_{31})=30\), giving the stated binary limit. The conclusion regarding agreement modulo 4 applies to these exponents for \(r\ge3\).

**FINITE-EXACT, prescribed insertion family.** In each of the three stated exponent-31 right addresses, choose a suffix width \(1\le s\le\operatorname{bitlength}(y_0)\) and insert a positive run of zeros or ones immediately to its left. The last bit is retained. Across all these choices the ray above is the unique affine ray reaching \(E_{61}\pmod{4095}\). The same run can be cut at widths 4, 5, or 6; all three descriptions give identical addresses. Thus this particular insertion family cannot repair its binary obstruction by moving the cut.

For the finite test write \(y_0=BU+v\), \(B=2^s\), \(R=\operatorname{rev}_s(v)\). The full center after a positive-width insertion is \(E_{31}+\delta+S(2^m-1)\). For ones the coefficient and boundary adjustment are \((U+1)(cB+h(B-1-R))\) and \(h[U\text{ even}][2v<B]\). For zeros they are \(U(cB+h(B-R))\) and \(h[U\text{ odd}][2v\ge B]\). These follow by applying the suffix recurrence before and after the edit, using \(p(2^m U)=p(U)+(2^m-1)U+[U\text{ odd}]\) for \(m>0\). Enumerating the twelve values of \(m\bmod12\) completes the finite test. Run `python3 -B check_61.py` in this directory to reproduce it.
