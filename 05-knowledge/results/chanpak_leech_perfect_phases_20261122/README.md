# Independent binary and odd corrections on a physical unequal-flank seam

**PROVED.** Let \(N\) be any odd target. Fix \(a\ge1\) and nonnegative precisions \(b_3,b_5,b_7,b_{13}\). The same specified arc of a six-vertex core has strong realizations with counts \(L,L+2\) satisfying
\[
 L\equiv N\pmod{2^a3^{b_3}5^{b_5}7^{b_7}13^{b_{13}}}.
\]
The left address programs the binary digits while keeping the starting odd residues fixed. The length of a high block of ones on the right then programs all four odd lifts. The left word has length \(O((a+1)^2)\), uniformly in the odd precisions. The right word has length \(O(a+1+\prod_\ell\ell^{\max(0,b_\ell-v_\ell(4095))})\). The sizes increase with precision. Equality at the absolute target requires additional height control.

The earlier one-apex result gives arbitrary fixed-modulus phases with nearly tied flanks. Here both the arc and its six-vertex deletion core stay fixed, the two arms control different parts of the modulus, and an actual perfect-number seam starts a second, sharper lift. The hostile example is exponent 19 transported toward exponent 31: odd accuracy can grow indefinitely on that single ray, but its eighth binary digit always fails. [One-apex phases](../chanpak_apex_flip_local_global_20261122/README.md) · [Unequal-flank formulas and cube-root boundary](../chanpak_gap_star_ports_20261122/README.md) · [Physical centers at perfect numbers](../chanpak_gap_star_ports_20261122/PERFECT-TARGETS.md)

A different core gives an absolute edge at exponent 31. Its 60-vertex witness also seeds a lift toward exponent 61. [Absolute exponent-31 edge](../chanpak_e31_physical_20261122/README.md) · [Transfer to exponent 61](EXPONENT-61.md)

## 1. An odd wheel with an explicit inverse

Write \(p(w)=H_{\rm apex}(w)-1\), with \(p(0)=0\), \(p(2u+1)=p(u)\), \(p(2u)=p(u)+u+[u\text{ odd}]\). Use the two even addresses \(x,y>0\) on the six-core of the cited proof. The center \(Z\) of its actual single-arc pair is
\[
 Z=112+162x+29p(x)+(129+320x)y+29p(y).
\]
Choose \(B=2^s\), an even suffix \(0\le r<B\), and \(y=B(2^m-1)+r\), \(m\ge1\). With \(c=129+320x\), the suffix recurrence gives \(Z=I+Q2^m\), where
\[
 Q=cB+29(B-1-\operatorname{rev}_s(r)),\qquad
 I=112+162x+29p(x)+cr+29\big(p(r)+[2r<B]\big)-Q.
\]
Indeed the leading \(2^m-1\) has \(p=0\) and is odd. Reversal uses all \(s\) digits including padded zeros.

**FINITE-EXACT base covering.** At \(x=2\), enumerate \(1\le s\le15\), every even \(r<2^s\), and \(0\le m_0<12\). Retain only \(\gcd(Q,4095)=1\). The residues \(I+Q2^{m_0}\) cover every class modulo \(4095=2^{12}-1=3^2\cdot5\cdot7\cdot13\). The numbers covered cumulatively by suffix width are 12, 23, 23, 33, 80, 173, 301, 593, 1052, 1791, 2693, 3544, 3977, 4087, 4095. The script implements precisely this finite enumeration. A base residue with \(m_0=0\) is used at positive multiples of 12.

Put \(v_\ell=v_\ell(4095)\) for the four odd primes. Given a target \(C\) matching its seed modulo 4095, take \(m=m_0+12t\). Every residue in this initial class modulo \(\ell^{b_\ell}\) is attained, with a unique \(t\bmod \ell^{\max(0,b_\ell-v_\ell)}\). Indeed \(2^{12}=1+4095\), and raising \(1+4095\) to an integer \(u>0\) gives difference from one of valuation \(v_\ell+v_\ell(u)\), by binomial expansion. The unit multiplier \(Q2^{m_0}\) preserves that valuation. The period thus equals the number of lifts, proving surjectivity. The four periods are pairwise coprime and the Chinese remainder theorem supplies a common \(t\). Increasing it by their product lets us require \(m\ge\max(1,a)\) within the stated right-width bound.

## 2. High binary digits that preserve the odd seed

For an even address \(x\), write its occupied binary positions as positive exponents \(d_i\). The recurrence gives
\[
 p(x)=\sum_i2^{d_i}-\sum_{d_i>d_j}\big(2^{d_i-d_j-1}+[d_i-d_j=1]\big).
\]
A given mark contributes \(2^{d_i}\) in isolation; each lower mark removes its indicated zero-position contribution in the recurrence. Start at \(x=2\). We add groups of twelve high marks. In decreasing order their exponents have residues modulo 12 equal, up to a common translation, to
\[
                       (8,3,6,9,7,0,4,2,5,1,10,11).
\]
They include each residue exactly once. The sum of the corresponding powers of two vanishes modulo 4095. Their internal pair cost does too: for this ordered list \(d\),
\(\sum_{i<j}2^{(d_i-d_j-1)\bmod12}=12285=3\cdot4095\). Choose their actual successive gaps at least \(a+1\) and at least two. Any specified relative residue is then possible within an additional eleven positions. Between two groups the modular cross contribution factors through the sum of powers on the higher group, hence also vanishes. The same argument handles cross contributions to the original low address. Every gap used below is at least two.

Place the first group above position \(a+1\), at distance at least \(a+1\) from the old mark. Between neighboring groups choose the gap from the highest mark of the lower group to the lowest of the upper group to be \(j+1\), with \(1\le j<a\). A common translation of the latter group's twelve residues accommodates this gap. Modulo \(2^a\) the added marks leave \(x\) fixed. All pair costs vanish except the one across each such junction, which contributes \(2^j\). Thus \(p(x)\) changes by minus the sum of the requested powers modulo \(2^a\), and both \(x\) and \(p(x)\) stay fixed modulo 4095. No unintended cross terms survive since every other cross-gap contains an internal gap of at least \(a+1\). At most \(a\) groups suffice by binary expansion. This proves the word-length estimate.

Given \(C=N+1\), choose an odd seed in §1. Compute its initial \(I_0\). It is even, since \(I_0+Q2^m\) for \(m>0\) is the center of neighboring odd counts. Put junctions at the 1 digits of \(29^{-1}(I_0-C)\bmod2^a\); this representative is even. The altered \(I\) now equals \(C\) modulo \(2^a\) and retains its original value modulo 4095. The same is true for \(Q\) modulo 4095, so it is still a unit there. Choose \(m\) by §1 to complete the proof. The lower and upper tournaments differ only at the certified internal arc for *all* these choices.

For \(C=E_{31}=2^{30}(2^{31}-1)\) and precisions \((a,b_3,b_5,b_7,b_{13})=(32,3,2,2,2)\), the algorithm gives a left address of 6965 bits, \(s=11,r=630,m=3204\). Its center matches \(E_{31}\) modulo 24,007,471,320,268,800. Its product term \(320xy\) exceeds \(E_{31}\). This illustrates exactly what the modular inverse retains and what it costs in absolute height.

## 3. Transporting the physical exponent-19 seam

The earlier exponent-19 lower witness has deletion read 143, tagged-start mass 173, lower apex count 291, and right address \(y_0=553528478\). Write \(y_0=2^{17}\cdot4223+11422\). Insert \(m\) ones after the high prefix, so that \(y_m=2^{17}(4224\cdot2^m-1)+11422\). The new high part is obtained by appending ones to 4223, leaving its \(p\) value and odd parity intact even at \(m=0\). The one-sided center is \(292+173y_m+143p(y_m)\), with endpoints joined by the original physical arc. The suffix formula now gives exactly
\[
                   Z_m=E_{19}+S(2^m-1),\qquad S=137468394624 .
\]
Here \(v_2(S)=7,\ v_3(S)=1\) and \(\gcd(S,5\cdot7\cdot13)=1\). Consequently for \(m=12t\), the odd image to any precisions is precisely the class of \(E_{19}\) at starting depths \((3,1,1,1)\), truncated at the requested precisions, by the same lifting argument as §1.

Every \(E_u\) for \(u\equiv7\pmod{12}\) lies in this odd class. To see the extra power of three, set \(X=2^u,Y=2^{19}\). The difference is \((X-Y)(X+Y-1)/2\); the first factor vanishes modulo 4095 and the second modulo 3. In particular arbitrarily high joint precision at all four odd primes transfers the known exponent-19 edge to the residues of exponent 31. For any \(m>0\) the binary valuation of \(Z_m\) is exactly seven. At each such exponent \(u\ge8\) the physical ray also matches \(E_u\) modulo 128 but fails modulo 256 if \(u\ge9\). Changing the *other* flank in §2 removes this obstruction.

## 4. Four exponent phases and the missing sign

For each prime exponent \(u>3\), its residue class modulo 12 determines \(E_u\bmod4095\). The rows for \(u\equiv1,5,7,11\) respectively give \(1,496,4033,3583\). The same table applies to all positive exponents in these classes, by the period of 2 modulo 4095. To determine which high-precision tests distinguish exponents inside a row, write \(z_u=2^{u+1}-1\). Then \(8E_u=z_u^2-1\).

Let \(v=v_\ell(4095)\) and fix an exponent \(u\) in the row. As the exponent changes by multiples of 12, \(z\) visits every lift of \(z_u\bmod\ell^v\). For any precision \(b>2w\), \(w=v_\ell(z_u)\), equality \(E=E_u\pmod{\ell^b}\) holds exactly at \(z\equiv z_u\) or \(-z_u\pmod{\ell^{b-w}}\). Indeed equality of squares forces the same valuation, and after dividing out \(\ell^{2w}\) at most one of the two factors is divisible by the odd prime. The reflected solution persists inside the chosen exponent row at arbitrarily high precision exactly when \(w\ge v\). The orders of 2 modulo \(9,5,7,13\) are \(6,4,3,12\). Thus the persistent reflected primes in the four rows are respectively
\[
             \varnothing,\qquad \{3,7\},\qquad \{5\},\qquad \{3,5,7,13\}.
\]
At exponents 19 and 31 the valuations of \(E_{31}-E_{19}\) at \(3,5,7,13\) are exactly \(3,2,1,1\); the first distinguishing powers are \(81,125,49,169\). Only the 5-adic reflected exponent branch lives in their common row. The order and lifting argument is a different observable from the earlier test for divisibility of a half-perfect center. [Wieferich width periods](../chanpak_perfect_wieferich_seams_20261122/PROOF.md)

## 5. Testing a two-coordinate gluing rule at the Leech holes

The six labels are the Coxeter numbers \(h\) of the pure-A Niemeier root systems \(A_{h-1}^{k}\), \(k(h-1)=24\), having **even** component count \(k\). The other pure-A types have \(h=9,25\) and odd \(k\). The Niemeier–Leech correspondence supplies the associated deep-hole type; the Leech lattice itself has no roots. [Earlier six-label dictionary](../chanpak_mrb_leech_clock_20261114/MRB_LEECH_SOURCES.md) · [Niemeier lattice catalog](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/) · [Möller and Scheithauer](https://annals.math.princeton.edu/2023/197-1/p04)

**Theorem.** For \(u\ge3\), place \(E_u-1,E_u+1\) in the first two discriminant coordinates of each such root lattice and zero in the others. Allow a choice of component signs. Actual Niemeier glue accepts this vector precisely for \(h=13\) and \(u\equiv11\pmod{12}\).

The discriminant group of \(A_{h-1}\) is cyclic of order \(h\). In the sum-zero model the first fundamental weight has norm \((h-1)/h\), and its multiples generate this group. A coset \(j\) has minimum norm \(j(h-j)/h\) for \(0\le j<h\): one minimizes the sum of squares of coordinates with the same fractional part and sum zero by taking the two nearest possible coordinate values. Glue elements must have total norm zero modulo \(2\mathbb Z\). At \(h=2,4\) the two labels are odd, and their squared sum is 2 modulo 8; it fails even isotropy at both labels. At \(h=3,7\) the norm condition would require \(E_u^2+1=0\pmod h\), impossible.

At \(h=5\), the same norm condition holds exactly if \(u\equiv3\pmod4\). Then \(E_u=-2\pmod5\), giving absolute coset labels 2 and 1. A glueword supported on just these two components would produce a new norm-two vector of norm \(6/5+4/5\) outside the root lattice, contradicting that the full root system is \(A_4^6\). At \(h=13\), the root discriminant group is \((\mathbb Z/13)^2\), so the index of an even unimodular overlattice is 13 and its glue is one isotropic line. Its equation is \(y=\pm5x\); changing a component sign exchanges the two choices. Among powers \(2^u\bmod13\), the norm test passes just at \(u\equiv11\pmod{12}\), where \(E_u=8\pmod{13}\) and the two labels give \((7,9)\), on the line \(y=5x\). The elementary residue assertions at 5 and 13 can be checked in a period of length twelve.

The two possible split primes have a common explanation. For \(X=2^u\), the function \(E=X(X-1)/2\) has its only odd-characteristic critical point at \(2X-1=0\). There \(64(E^2+1)=65\), so exactly the odd primes 5 and 13 are ramified zeros of this norm polynomial. Indeed its derivative is \(E(2X-1)\), and at a zero \(E\ne0\). On each passing residue above, \(z=2^{u+1}-1\) is divisible by the respective prime. The identity
\[
                     64(E_u^2+1)=z^4-2z^2+65
\]
shows that its valuation there is exactly one. Even the passing norm residue cannot lift to its square modulus. This is a concrete test of the proposed two-coordinate map; the map forgets the marked arc and its two flanks. In particular the exponent-19 seam exists in tournaments despite failing every actual gluing test in this six-type dictionary. [Signed six-label tournament constructions](../chanpak_six_signed_height_pairs_20261118/README.md) · [THM-4758, deep-hole theta identities](../../../01-canon/theorems/THM-4758-leech-deep-hole-swap-theorem-and-the-deep-hole-theta-identity.md)

Run `python3 -B check.py` in this directory for the finite seed covering, physical phase controls, and the exponent-31 decoder. To prescribe a different even center and precisions use `--center C --two-depth A --odd-depths B3 B5 B7 B13`. The returned width can be used in compressed form when an explicit tournament would have very large order.
