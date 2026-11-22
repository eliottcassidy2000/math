# An exact unequal-flank selector and a periodic-carry obstruction

**PROVED:** the absolute one-arc target test below is necessary and sufficient for a specified pair of marked star cores. It is not a uniform hitting theorem. For the exponent-31 even perfect number, mere imbalance does not repair the error in the tested persistent small-core bank. A separate periodic-carry theorem excludes using a fixed repeating binary correction to extend such a repair across exponents.

## Factor completion after eliminating the full binary read

Use the [two-star-arm identity](../chanpak_gap_star_ports_20261122/README.md). Fix a labelled deletion core, a tag set and a core arc to flip. Let the two marked core profiles have Hamilton counts \(h,h'\), optional two-path reads \(Q,Q'\), and unsigned end and start reads \(R,R'\) and \(S,S'\). For each choice of parities \(x\bmod2,y\bmod2\), absorb the two signed parity terms into the constants \(K,K'\). The exact counts are
\[
L=Qxy+Rx+Sy+K+h[p(x)+p(y)],\qquad
L'=Q'xy+R'x+S'y+K'+h'[p(x)+p(y)].
\]
Here \(p(0)=0\), \(p(2u+1)=p(u)\), and \(p(2u)=p(u)+u+[u\text{ odd}]\) for \(u\ge1\). In particular \(0\le p(w)\le w\); for the even update use that \(p(u)\le u-1\) for positive odd \(u\), directly from the odd recursion. The two positive binary addresses guarantee a strong host and the two hosts differ exactly at the chosen core arc.

To realize the prescribed **absolute** pair \(N,N+2\), set
\[
A=h'Q-hQ',\quad B=h'S-hS',\quad C=h'R-hR',\quad
T=(h'-h)N-2h-h'K+hK'.
\]
Eliminating the **entire** term \(p(x)+p(y)\) gives \(Axy+By+Cx=T\). Thus for \(A\ne0\)
\[
                 (Ax+B)(Ay+C)=AT+BC. \tag{1}
\]
If the right side is nonzero, enumerate its signed divisor pairs, reject wrong residues modulo \(A\), recover both positive addresses with the required parities, and test \(L=N\) with the exact binary recurrence. No second count check is needed mathematically: equation (1) guarantees \(h'L-hL'=h'N-h(N+2)\), and \(h>0\). The signed divisor enumeration matters even though addresses are positive. If the right side vanishes, one factor is zero, reducing the problem to the inherited fixed-left or fixed-right inverse. If \(A=0\), solve the linear equation \(By+Cx=T\); if it is identically zero, the carry remains an ordinary two-address test. In every case the fallback search is finite because \(Qxy\le N-K\) and \(x,y\ge1\). Exchange the roles of the two profiles to test either orientation of the physical arc. This is not a polynomial-time factoring guarantee.

If \(Q\ne Q'\), there is an independent and useful finite cap on the smaller address:
\[
 \min(x,y)\ \le\ 
 \frac{2|h'-h|+|R'-R|+|S'-S|+|K'-K-2|}{|Q'-Q|}.
\]
Divide the exact gap equation by \(\max(x,y)\), use \(p(x)+p(y)\le x+y\), and round down. Enumerate the capped address in either position; the eliminated linear equation gives the other address directly except on the explicitly handled zero-coefficient wall. Thus in this regime no factorization of a giant integer is needed. Also, if \(|AT+BC|\) is prime, (1) can never give two addresses with both \(|Ax+B|,|Ay+C|>1\). This is an absolute obstruction to such aspect ratios even if every previously tested congruence passes.

**Nonpersistent literal control.** Use lower core \((10,4,9,18,7)\), reverse only \(0\!-\!1\), and take tag mask 2. The profiles \((h,Q,(e,o),(e',o'),k)\) are \((13,104,(5,9),(27,38),23)\) and \((11,100,(9,11),(32,40),29)\). With \(N=1221\) and even parities, \(A=-156,B=-106,C=-221,T=-2344\). The signed factors \(-730,-533\) give \(x=4,y=2\). The binary recurrence gives \(1221,1223\), and literal Hamilton DP and the single marked arc both check. The formulas cannot be treated as equal-slope rails: both \(h\) and \(Q\) change.

## Exponent 31: what was tested and what was not

For \(E_{31}=2305843008139952128\), set \(N=E_{31}-1\) and require one physical pair with heights \(N,N+2\). In the supplied order \(4,\dots,8\) core classes, every core arc reversal and nonempty proper tag mask was considered in the following restricted tests:

- **Persistent unequal-flank charts:** equal \(h,Q,R,S\) on both sides and adjusted intercepts differing by two. There are 132 distinct scalar charts, with parity split \(40,26,26,40\) in order \(00,01,10,11\). The bank is exactly closed under exchanging flanks, \(Q\ge320\), and the smallest adjusted intercept is 111. Exact suffix inversion of the larger address after enumerating the smaller positive address up to 85,000,000 found no straddling pair. This cutoff exceeds \(\lfloor\sqrt{(N-111)/320}\rfloor\) and is complete for this two-positive-address bank.
- **Fixed-left walls:** equal \(h\), but unequal slopes may meet at a single \(x\) via \(S+Qx=S'+Q'x\); an unequal intercept slope can likewise select one \(x\), separately by parities. The deduplicated bank has 766 exact scalar charts (including empty fixed flanks). Their one-variable exact inverses give no straddle.
- The exact linear tests when \(h=h',Q=Q'\) and both unsigned slope differences are nonzero also found no straddle in those orders. Some other nonpersistent regimes were explored, but this packet does **not** claim a complete census of the remaining factorization cases.

The existing [star-arm packet](../chanpak_gap_star_ports_20261122/STAR-ARMS.md) already realizes \(E_{31}+1\) as an ordinary count using factors; the lower neighbor also has the known signed-rail scalar realization. Neither identifies endpoints of one arc. The previously successful \(17,19\) case uses different literal sources for different offsets. The current [nontransitive-gate construction](../chanpak_e31_physical_20261122/README.md) does produce a literal exponent-31 straddle on 60 and 62 vertices. It changes the small-core setup by growing a marked nontransitive receiver, so it is not a counterexample to the scoped negative tests here. The \(\{4,2,5,3,7,13\}\) [signed-height portfolio](../chanpak_six_signed_height_pairs_20261118/README.md) supplies exact mantissa/height selectors and multiplicative companions, not a same-arc selector. Likewise the earlier [Golay/Leech and mod-24 comparison](../chanpak_mod24_moonshine_ports_20261122/README.md) supplies residue observations, not the exact signed divisor and carry agreement of (1). Neither arithmetic analogy warrants transferring scalar attainment to a marked edge.

## Fixed repeated blocks cannot repair across exponents

Let \(B=b_1\cdots b_k\) be a fixed block with binary value \(v\), \(q=2^k\), and \(a=q-1-\operatorname{rev}_k(v)\). A fixed nonempty prefix, then \(n\) repetitions of \(B\), then a fixed suffix gives an address with an eventually exact form \(\alpha q^n+\beta\). The binary surcharge is eventually
\[
             \gamma_1 q^n+\gamma(B)n+\gamma_0,\qquad
 \gamma(B)=p(v)+b_k(1-b_1)-\frac{av}{q-1}.
\]
Indeed appending \(B\) increments \(p\) by \(aw+p(v)+(1-b_1)[w\text{ odd}]\) on current address \(w\); after one repeat the latter parity is \(b_k\). Summing \(w'=qw+v\) gives the formula. A fixed suffix adds only an exponential term and a constant.

**Strict drift.** The exact inversion expansion
\[
 p(v)=\sum_{i<j}b_i(1-b_j)2^{j-i-1}+\#\{i<k:b_i=1,b_{i+1}=0\}
\]
follows from the appending recursion. On the bi-infinite periodic extension \(b_i\), expand \(av/(q-1)\) geometrically, subtract the finite weighted-inversion terms, and obtain
\[
 \gamma(B)=\#\{\text{cyclic }10\}-\sum_{i:b_i=1}\sum_{d\ge1}(1-b_{i-d})2^{-d-1},
\]
where \(i\) ranges over one period. In each cyclic run of \(r\) ones, its \(j\)-th one (\(j=0,\ldots,r-1\)) has \(j\) previous ones so the inner sum is at most \(2^{-j-1}\). Thus the run pays less than one to the double sum, while contributing exactly one cyclic \(10\). Therefore \(\gamma(B)>0\) for every nonconstant block; for all-zero or all-one blocks it is zero. In particular reflected-complementary periodic tails still cannot cancel their drifts.

Apply a fixed star core to any two addresses built as fixed prefixes, repeated fixed blocks, and fixed suffixes, permitting a constant flank. Along an exponent progression \(p=p_0+d n\), with fixed repeat rates, all address, product and unsigned terms are combinations of constants and powers of 2 raised to linear functions of \(n\); fixed parities give fixed signed terms. In contrast the binary contribution has linear coefficient \(h\) times the sum of the two nonnegative repeat drifts (weighted by their repeat rates). The perfect target \(E_p+c\) has no linear term. Distinct growing exponentials dominate successively, then the nonzero linear term dominates the constant. Thus if either repeated block is nonconstant and grows, exact equality cannot hold for infinitely many \(n\).

This includes an attempt to absorb an odd factor of \(Q\) by fixed rational affine addresses \(x=A_0 2^p+B_0\), \(y=C_0 2^p+D_0\) on a progression on which the addresses are integral: a positive rational non-dyadic leading coefficient has a nonconstant eventually periodic binary period, and a fixed additive correction changes only a bounded suffix (a nonconstant period prevents an arbitrarily long carry or borrow). Truly dyadic coefficients can avoid the periodic drift, but if both flanks scale as \(2^p\), matching \(Q A_0 C_0=1/2\) with positive dyadic \(A_0,C_0\) forces \(Q\) to be a power of two. None of the 132 persistent charts in the supplied order-8-and-below bank has that property. General, nonperiodically target-dependent choices remain open.

## Cross-exponent width rigidity in the new nontransitive-gate family

The [successful exponent-31 packet](../chanpak_e31_physical_20261122/README.md) also gives an exact unrestricted pair of addresses \(v\ge1,y\ge0\), with lower count
\[
L(v,y)=25[(1+p(v))(1+p(y))+(v+1)(y+1)-[v,y\text{ both even}]]
       +28[Q_v(y+1)-m_v]-10 .
\]
Here \(V(v)\) consists of its transitive binary spine and a sole tagged apex; \(Q_v\) is the two-slot count and \(m_v\) is the weighted path-ending-at-apex read. This is a different family from a fixed core with two transitive flanks, so the preceding periodic-carry argument does not itself apply.

**Width law.** Write \(t\) for the number of trailing ones in the binary word of \(v\), and \(b\) for the remaining bit length (zero if the word is all ones). Then uniformly
\[
Q_v\asymp 2^t3^b,\qquad
L(v,y)\asymp 2^t3^b(y+1).
\]
Indeed the trailing ones correspond to ordinal sinks, each multiplying the two-slot count by two. The remaining nonempty word starts in one and ends in zero. For this word, select the slot containing the apex in two ways. Assign its first spine vertex after the apex in that slot, and its last spine vertex before the apex. Each interior vertex may freely choose the before path, after path, or other slot; transitivity and the two fixed contacts guarantee legality. This gives \(2\cdot3^{b-2}\) covers; the ambient three-choice description gives upper bound \(2\cdot3^b\). With no remaining spine \(Q=2\). To bound the height, \(2m_v\le Q_v\) because choosing the slot of the selected apex-ending path gives two disjoint subclasses of all covers. The direct formula then gives a positive constant times \(Q_v(y+1)\) as a lower bound. Conversely \(p(v)\le v,\ p(y)\le y\) and \(v+1\le 2^{\operatorname{bitlength}(v)}\le Q_v\) imply the stated upper bound.

Consequently any exact \(E_p+O(1)\) realization in this two-word gate satisfies
\[
b\log_2 3+t+\log_2(y+1)=2p+O(1).
\]
In particular there cannot be infinitely many exact hits under fixed positive rational linear active width \(b=\alpha p+O(1)\), rational linear padding \(t=\beta p+O(1)\), and rational linear suffix width \(\log_2(y+1)=\delta p+O(1)\): that would force \(\alpha\log_2 3+\beta+\delta=2\), impossible since distinct powers of 2 and 3 never agree. The successful isolated exponents and bounded active receivers are not excluded. This identifies the specific cost of genuinely growing the nontransitive gate at perfect-number scale: its natural width is ternary, not dyadic.
