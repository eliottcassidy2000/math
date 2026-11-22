# A single apex flip at prescribed gap, phase, and scale

**PROVED.** Fix positive integers \(d,M\). At every sufficiently large **even** order, each residue of an odd integer modulo \(M\) occurs as the lower endpoint \(L\) of an edge \(L\to L+2d\) in the arc-flip graph of strong tournaments. Both endpoints have the same transitive spine and differ at one apex contact. More precisely, for each permitted residue and every sufficiently large real \(X\), such a lower endpoint lies in
\[
[X, X+O_{M,d}(\sqrt X)].
\]
The exponent is sharp for this one-contact family, even with \(M=1\). At sufficiently large **odd** order there are exactly \(4d\) positively oriented marked address pairs in the whole normalized family. For \(d=1\), all four lower counts converge to \(-1\) in the 2-adic integers.

This gives a congruence and spacing theorem for the actual one-arc rungs proposed in [THM-5271, flip potential and seam ladder](../../../01-canon/theorems/THM-5271-the-count-is-a-potential-on-the-flip-metagraph-the-seam-ladder-implies-chan-pak-the-microcosm-3-1-6-1-21-the-square-sum-cycle-at-32-and-the-order-10-near-join-reach.md). Chan–Pak and the full ladder are **OPEN**.

## Inherited interface

On a transitive spine with increasing arcs, a word records which vertices the apex beats (bit 1). In a contact flip \(U0V\to U1V\) put \(u=\operatorname{bin}(U)\), and let \(w\) be the binary value of \(V\) after reversing and complementing its bits. Strongness of both endpoints is equivalent to the entire words starting in 1 and ending in 0. Necessity follows from the source or sink at an incorrect end, and sufficiency from the directed cycle through the apex and the whole spine. Thus for a strong pair both flanks are nonempty and use the standard positive addresses without leading padding. The common Hamiltonian path consisting of the apex followed by the spine survives the flip.

Write \(F(0)=1\), \(F(2q+1)=F(q)\) and \(F(2q)=F(q)+q+[q\text{ odd}]\) for \(q>0\) in the last rule. The [one-site inverse and capacity theorem](../chanpak_mixed_defects_20261115/SINGLE_SITE.md), using [THM-5070, port-chain transfer](../../../01-canon/theorems/THM-5070-port-chain-transfer-for-odd-cycle-packing-counts-and-ternary-coverage-to-a-million.md), gives
\[
L=2uw+u+[u\text{ odd}]+F(u)+F(w)-1,\qquad
H(U1V)-L=2(\lceil w/2\rceil-\lceil u/2\rceil).
\tag{1}
\]
The negative orientation is obtained by converse and reversal of the spine. The cited theorem already gives four inverse candidates per target and \(\Theta_d(\sqrt X)\) distinct lower endpoints up to \(X\). The [small-component theorem](../chanpak_physical_routes_20261115/ONE_APEX.md) also bounds how many *successive* small contact flips can be used. Our construction selects distinct physical edges.

The previous [all-modulus one-apex theorem](../one_apex_all_residues_20261116/README.md) controls one scalar. The [fixed-order dyadic pair decoder](../chanpak_exponent_dimension_20261117/apex/DYADIC.md) controls a two-contact gap. Here one common prefix steers two dependent addresses while keeping one contact, their gap, and their order fixed. The additional information is which side of a power of two each address occupies.

## A binary switch with a prescribed additive order

Set \(Q=\mathrm{lcm}(4,M)=2^a m\) with \(m\) odd and \(a\ge2\). Choose \(\ell\ge a\) with \(2^\ell\equiv1\pmod m\) (vacuous at \(m=1\)). Starting from a positive prefix \(q\equiv-1\pmod {2^a}\), append either of two equal-width blocks
\[
A=1^{2\ell},\qquad B=0^\ell1^\ell.
\]
Both leave \(q\bmod m\) fixed and end at \(-1\bmod 2^a\). The all-ones block leaves \(F\) fixed. The zero run in \(B\) contributes
\[
(2^\ell-1)q+1
\]
to \(F\), since exactly its first step has an odd prefix. This increment has a fixed residue \(D\pmod Q\), given by \(D=2\pmod {2^a}\) and \(D=1\pmod m\). In particular, \(2D\) generates the multiples of four modulo \(Q\). The recurrence for \(F\) also shows that, when fixed bits are appended, the extra increment to \(F\bmod Q\) only reads the old prefix modulo \(Q\). Thus changing the number of switched blocks translates **both** future reads by the same amount.

## Every phase at every sufficiently large even order

Take an arbitrary positive integer prefix \(y\), append \(\ell\) ones, and then append \(J=Q/4-1\) blocks, each \(A\) or \(B\). Call the resulting integer \(q\). Put \(s_0=8d+1\), choose \(K\ge3\) with \(2^K>s_0\), and set
\[
s=2^K q+s_0,\qquad u=2(s-d),\qquad w\in\{2s-1,2s\}.
\]
All three addresses are the same prefix \(q\) with fixed \((K+1)\)-bit suffixes. Both choices of \(w\) give gap \(2d\) in (1). If \(j\) blocks are switched, \(u,w\bmod Q\) do not change, and the lower endpoint increases by \(2jD\bmod Q\). For each branch this reaches one entire coset modulo four.

The two branches give the two odd cosets. Indeed let \(s=2t+1\); since \(s=1\pmod8\), \(t=0\pmod4\). The lower endpoint for \(w=2s\) minus the one for \(w=2s-1\) equals
\[
2u+F(2s)-F(2s-1)=2u+t+2-[t\text{ odd}]\equiv2\pmod4.
\]
Here \(u\) is even. Since any permitted residue modulo \(M\) lifts to an odd residue modulo \(Q\), this supplies the requested endpoint by choosing one branch and \(j\) between 0 and \(J\). An inverse requires only a linear congruence for \(j\). With \(b\) the bit length of \(y\), both addresses have length \(b+P\) for \(P=\ell+2\ell J+K+1\). The order is \(2(b+P)+2\). Every \(b\ge1\) works.

## The real scale of this local completeness

Let \(C=2^P\). The construction works for **every** prefix \(y\), placing \(u\) in \([Cy,C(y+1))\). For a positive fixed gap, (1) gives \(2d-1\le w-u\le 2d+1\). The recurrence gives \(1\le F(x)\le x+1\): at an odd argument the bound strengthens to \((x+1)/2\), which controls the parity bonus at its even child. Consequently every lower endpoint satisfies \(2u^2\le L\le 2u^2+O_d(u+1)\). Choose \(y=\lceil\sqrt{X/2}/C\rceil\). The selected residue then has a lower endpoint at or above \(X\), within \(O_{M,d}(\sqrt X)\) of it. In particular each permitted residue supplies \(\Theta_{M,d}(\sqrt X)\) distinct lower counts through \(X\): a lower count allows only \(O_d(1)\) choices of \(u\) by the same quadratic bounds, whereas the prefix intervals are disjoint.

Conversely, the inherited capacity estimate bounds even the entire one-site family in \([X,2X]\) by \(O_d(\sqrt X)\) points. Pigeonholing gives empty subintervals of length \(\Omega_d(\sqrt X)\). This proves sharpness of the exponent and identifies the absolute accuracy still needed for the Chan–Pak ladder.

## Odd-order addresses and their 2-adic limit

For positive gap write \(u=2r-\epsilon,\ w=2(r+d)-\eta\) with independent bits \(\epsilon,\eta\). At large odd order the widths are unequal. Since \(w-u\le2d+1\), their only possible widths are \(k,k+1\), at order \(2k+3\). Indeed an unbounded larger width forces \(u\) to infinity and hence \(w<2u\). Put \(T=2^k, p=T/2\). Then \(u<T\le w\). At \(r=p-d\) there are two choices (\(w=T\)); at \(p-d<r<p\) there are four each; at \(r=p\) there are two (\(u=T-1\)). All these addresses have the stipulated widths when \(k\) is sufficiently large. This gives exactly \(4d\) oriented marked address pairs, before identifications of isomorphic tournaments.

For each such pair \((u,w)=(T-i,T+j)\), the offsets depend only on \(d\) and the branch. Choose a fixed suffix width \(h\) with \(2^h>\max(i,j)\). The addresses have high words \(2^{k-h}-1\) and \(2^{k-h}\), followed by fixed suffixes. Appending a fixed suffix adds to \(F(q)\) an expression affine in \(q\) plus a fixed multiple of \([q\text{ odd}]\), by the recurrence. The high-word reads are 1 and \(2^{k-h}+1\) for \(k>h\). Hence (1) is a fixed integer quadratic in \(2^{k-h}\) on each branch. Each lower endpoint converges 2-adically to its branch's integer constant. By comparison, for even orders every odd 2-adic integer is a subsequential limit: apply the phase construction with \(M=2^a\) at increasing orders.

For \(d=1\) and \(k\ge3\), the full list of four *centers* \(L+1\) is
\[
2T^2-15T/4,\quad 2T^2-3T/2,\quad
2T^2+3T/2,\quad 2T^2+15T/4.
\]
In the same order the offsets \((u-T,w-T)\) are \((-3,0),(-2,0),(-1,1),(-1,2)\). Substitute \(F(T-3)=T/4+1\), \(F(T-2)=T/2+1\), \(F(T-1)=1\), \(F(T)=T+1\), \(F(T+1)=T/2+1\), and \(F(T+2)=3T/4+3\) from the recurrence. At order nine the centers are 98, 116, 140 and 158. All four lower counts are \(-1\bmod 2^{k-2}\) at order \(2k+3\). Thus a constraint on physical width selects a single 2-adic limit from a family which, without that constraint, realizes every odd phase.

## The forced factor at an odd-order unit seam

For every one of the four odd-order branches with \(d=1\), the center is \(2T^2\equiv2\pmod3\). Hence the lower endpoint is \(1\pmod3\), and the upper endpoint is divisible by 3. More precisely, on each individual branch the natural density of width indices \(k\) for which the upper has exact 3-adic valuation \(j\) equals \(2/3^j\), for every \(j\ge1\).

To see this, write the upper as \(G(T)=2T^2+cT+1\), where \(c\) is one of \(\pm3/2,\pm15/4\). All denominators are units modulo every power of 3. Modulo 3 both nonzero arguments solve \(G(T)=0\), and \(G'(T)=4T+c\) is nonzero there. Moreover \(2^{2\cdot3^{j-1}}=1+3^j b_j\) with \(3\nmid b_j\), by cubing from \(4=1+3\). Above any solution in \(k\bmod 2\cdot3^{j-1}\), its three lifts to the period \(2\cdot3^j\) change \(T\bmod 3^{j+1}\) by the three multiples of \(3^j b_j T\). The derivative gives precisely one lift solving the next congruence. Starting at \(j=1\) yields two solution classes in every period \(2\cdot3^{j-1}\). Subtracting the density for divisibility by \(3^{j+1}\) from that for \(3^j\) proves the claim.

At even orders the phase theorem gives a different prediction. For any fixed finite set of primes \(S\) and any half-gap \(d\), both endpoint counts can simultaneously be chosen coprime to \(\prod_{p\in S}p\) with the stated square-root spacing. At odd \(p\) choose any lower residue outside \(\{0,-2d\}\); at 2 choose 1. The Chinese remainder theorem supplies a single allowed phase. This removes any prescribed finite list of prime divisors in this construction. Primality of both endpoints remains open.

The map to \((d,L\bmod M)\) forgets both absolute height and the edge's location in the flip graph. The prefix interval restores height to square-root precision. The side of the binary width boundary restores the order distinction. The native single-arc edges here need no identification with the path-reversal or tiling-line metagraph. The applicable [research pattern](../../../00-navigation/META-PATTERNS.md#separate-local-support-from-bounded-height-coverage) is separation of local support from bounded-height coverage.

## Reproduction

Run `python3 replay.py` in this directory. Literal subset path counts check every normalized marked flip through order 9. The program also decodes every permitted residue for the moduli and gaps in its printed universe, checks the odd-order boundary independently by width enumeration, and reconstructs a requested address with `python3 replay.py --decode MODULUS RESIDUE HALF_GAP PREFIX`. The all-order statements follow from the proofs above.
