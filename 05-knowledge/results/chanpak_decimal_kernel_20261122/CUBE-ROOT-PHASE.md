# Cube-root bracketing in a prescribed congruence class

**PROVED.** Fix any one of the eight canonical Farey charts, a positive integer \(m\), and an odd phase \(a\). For every sufficiently large real target \(T\) there are two positive flank pairs whose lower counts satisfy
\[
       H_-<T<H_+,\qquad H_-\equiv H_+\equiv a\pmod m,
       \qquad H_+-H_-=O_{m,\mathrm{chart}}(T^{1/3}).
\]
Each lower tournament is strong and reverses its designated physical arc to a strong tournament with two more paths. This controls the congruence and the absolute error in the *same* construction, including when \(m\) contains powers of all 15 supersingular primes. The phase can be chosen from the target when \(T\) is an odd integer. The threshold and constant depend on \(m\).

The inherited cube-root construction uses complementary low digits to cancel the main fluctuation in the two digital readings. An independent local inverse could change the sum of the two addresses, costing order \(\sqrt T\) in height. Here the [complementary controller](COMPLEMENTARY-CONTROLLER.md) changes the residues while keeping the sum of the flanks *exactly fixed*. Its common length bounds the step between successive controlled points. [Original cube-root unit-gap theorem](../chanpak_gap_star_ports_20261122/README.md) · [One-word local inverse](LOCAL-WORDS.md)

Write \(L=k+qxy+t x+c y+r[x\text{ odd}]+e[y\text{ odd}]+h(f(x)+f(y))\). In these charts \(h\) is odd, \(\gcd(t,h)=1\), and \(q,c,t,r\) are even. Write the requested modulus as \(2^v M\), \(M\) odd. The controller supplies a fixed length \(\ell\) depending on that modulus. Put \(W=2^\ell\), and \(g=\gcd(h,M)\). For a sufficiently divisible dyadic factor in \(S+1\) and
\[
                         t S\equiv a-k-r\pmod g,
\]
every pair of high words with \(U+V=(S+1)/W-1\) has complementary \(\ell\)-bit continuations \(x=WU+i,\ y=WV+(W-1-i)\) for which \(L\equiv a\pmod m\). In particular, the local adjustment can take place in *every* interval of \(W\) consecutive left addresses starting at a multiple of \(W\) on this same line.

### The absolute scale on the controlled line

Set \(z=\lfloor\sqrt{T/q}\rfloor\). Choose a power of two \(B\asymp z^{1/3}\); for large enough \(T\) it is a multiple of any prescribed fixed dyadic factor. Take \(S=Bv'-1\) with the displayed congruence, choosing the least such \(S\ge 2z+2B\). As \(tB\) is a unit modulo \(g\), \(2z+2B\le S<2z+(g+2)B\). Put \(w=S/2\), a half-integer. Set \(x=w+d, y=w-d\). For any half-integer \(d=O_m(\sqrt{zB})\), write \(x=BU'+j,y=BV'+(B-1-j)\). Then \(U'+V'=v'-1\). The concatenation identity gives
\[
 f(x)+f(y)=\tfrac12(v'-1)(B-1)+O(z/B+B+|d|+1).
\]
Indeed the two suffix reversals sum to \(B-1\); their weighted high-word contribution differs from the stated mean by at most \(B|U'-V'|/2=O(|d|+B)\). The remaining terms are at most \(U'+V'+2B+2\).

Consequently the actual height is within \(O_m(z^{2/3})\) of
\[
 G(d)=q(w^2-d^2)+(t+c)w+(t-c)d+k+\tfrac h2(v'-1)(B-1).
\]
Choose one controlled point per \(W\)-block, starting with the first full block at \(x\ge w\). Its \(d=O(W)\), giving \(G(d)-T\ge 2qzB-O_m(z)>0\) for large \(z\). End with the first full block at \(x\ge w+C\sqrt{zB}\), where \(C^2>2(g+3)\). Because \(w^2-z^2\le(g+2)zB+O_m(B^2)\), the last actual height is below \(T\) by order \(zB\). All addresses are positive. Consecutive selected points have \(|\Delta d|<2W\), so \(|\Delta G|=O_m(\sqrt{zB})\). The approximation then gives \(|\Delta L|=O_m(z^{2/3})\) for consecutive controlled points. Use the crossings at \(T\) (or at \(T\pm1/2\) for an integer target) to obtain strict brackets. Equivalently take the first height at most \(T\) and the first strictly below \(T\); if an exact equality intervenes, each nearest strict endpoint is within one step of that equality. Since \(z^{2/3}=O(T^{1/3})\), this gives the result.

An exact prescribed hit would follow at an individual target if a phase-matched witness were certified at distance *strictly smaller* than its modulus. The theorem supplies no such uniform inequality as the modulus grows. In particular it does not yield an exact edge at exponent 61. [Absolute exponent-61 search](FAREY-HYPERBOLA.md)

**FINITE-EXACT at exponent 61.** On spine 4, a search with the target phase modulo 510 (the primes 2, 3, 5 and 17) gives an absolute lower bracket at offset −16,163,825,630,511,139,740 and an upper bracket at offset +10,343,962,259,434,468,860 from the lower target. The corresponding addresses are (1076664420585741, 1076468036218610) and (1077021891657411, 1076110748975420). The same search also gives strict brackets for the full repeated-prime modulus 3600. Here the implementation looks at short intervals around the positive quadratic root and checks each integer height and phase directly. The smaller unconstrained offsets on this chart are in the [order-118 certificate](FAREY-BOUNDARIES.md). Run `python3 -B cube_root_phase.py` from this directory. [Exact finite search](cube_root_phase.py)

### Allowing the precision to grow with the target

**PROVED.** Fix one canonical chart. For any sequence of targets tending to infinity and moduli \(m=o(\sqrt{\log T})\), both brackets can be taken in any prescribed odd phase modulo \(m\) and within \(T^{1/3+o(1)}\) of \(T\). In particular one may require powers of all 15 supersingular primes with the minimum exponent tending to infinity sufficiently slowly.

Indeed, the complementary controller has \(\ell=O(m^2)\): its odd block length can be chosen even and at most \(2\varphi(M)+4\) by Euler's theorem, and its terminal controller has length \(O(\log^2(2m))\). Thus \(2^{a+1}W=T^{o(1)}\), eventually less than \(B\asymp T^{1/6}\), where \(2^a\) is the dyadic part of \(m\). The progression step \(gB\) is at most \(hB\) for this fixed chart. The constant in the two-address approximation and the choice of the terminal displacement can therefore be independent of \(m\). For a step of at most \(2W\), the smooth height varies by \(O_{\rm chart}(W\sqrt{zB}+W^2)\); the nonsmooth error remains \(O_{\rm chart}(z^{2/3})\). Combining the bounds proves the growing-precision claim.
