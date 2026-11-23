# When a modular arc target becomes an exact integer target

A finite modulus suffices once the lower tournament count is fixed **exactly** and the same physical arc is used at that modulus. Two instances are different in size. For one flip from 199 the target residues for 201 and 701 both become exact at ten binary digits, and both precisions are sharp. For the exponent-61 perfect number 276 binary digits suffice uniformly over all tournaments at the correct lower height. In the strong lower family, 274 suffice. The large bounds are sufficient precisions, with the existence of a matching arc still open.

## A bound after cutting out the arc endpoints

Write \(h=H(T)\) and \(\Delta=H(T')-H(T)\) for one arc reversal. For an \(n\)-vertex host with \(n\ge2\),
\[
                 1-h\ \le\ \Delta\ \le\ 2^{n-4}(h+3).
\]
Indeed each newly gained full path uses the new arc, say \(v\to u\). Remove its two endpoints and retain its ordered prefix and suffix. This injects the gained paths into ordered two-path covers of \(T-\{u,v\}\), allowing an empty slot. The two-slot bound \(Q(X)\le2^{|X|-2}(H(X)+3)\) follows by weighting the odd-cycle packings of \(X\): each nonempty packing uses at least three vertices per cycle, giving it at most one quarter of its ordinary weight after division by \(2^{|X|}\). Also \(H(T-\{u,v\})\le h\) by two path insertions. Positivity gives the lower bound. [Two-slot bound including the empty term](../chanpak_crossing_scales_20261115/FIXED_COUNT_TIGHTNESS.md) · [THM-1370, insertion monotonicity](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md)

If the edit is restricted to a consecutive interval of components, with common exterior multiplier \(g\), apply this bound to the interval and multiply by \(g\). With active order \(n\), the upper bound becomes \(2^{n-4}(h+3g)\). Multiplication of the exterior holds even if the edit joins the affected components.

**PROVED at exponent 61.** Put \(h=2^{60}(2^{61}-1)-1\). Then \(h+3<2^{121}\) and \(5^{156}\le h^3<5^{157}\). A strong tournament of height \(a\) has order at most \(1+3\log_5 a\) by Busch's bound (the formula also bounds a singleton at \(a=1\)). Hence a strong host of height \(h\) has order at most 157. The preceding bound gives \(\Delta<2^{274}\). Since \(1-h>2-2^{274}\), the residue \(\Delta\equiv2\pmod{2^{274}}\) implies \(\Delta=2\) for any such host. In an exact lower host already known to have order 93 the same argument needs only 210 bits. [Strong lower witness on 93 vertices](../chanpak_e61_compact_volchkov_20261122/STRONG-93.md)

For the all-order statement, assume only \(H(T)=h\) and \(\Delta\equiv2\pmod{2^{276}}\). If the edit lies inside one component, take that component for the active interval. Its order is at most 157 and \(h+3g\le4h\), giving \(\Delta<2^{276}\). For a cross-component edit the valuation of its gain is one. The interval between the endpoint components is empty or consists of a single vertex. Indeed a nonempty middle has its two-slot mass as a factor of the gain, and four divides that mass at order at least two. To see the latter, pair each subset and its complement: every summand is odd and there are an even number of complementary pairs. [Exact cross-component formula and valuation budget](../chanpak_kaprekar_decimal_20261123/FIRST-JUMP-SQUARES.md)

If the endpoint components have heights \(a,b\), then \(ab=h/g\). The whole active order is therefore at most \(3+3\log_5(h/g)\). For \(g=1\) this gives \(n\le159\), hence \(\Delta\le 2^{155}(h+3)<2^{276}\). If \(g>1\), then \(g\ge3\) and \(3\log_5 g>2\). We get \(n\le157\) and again use \(h+3g\le4h\). The same lower bound on negative changes excludes the negative congruence aliases. This proves the all-order criterion. A verifiable solution can use any method to produce the required exact lower host and modular arc. Independent solutions at a succession of moduli do not supply such a host at the specified precision. [Earlier correct lower height with a modular exit](../chanpak_decimal_kernel_20261122/BACKCHAIN-LOWER.md)

## Sharp precisions and a coherent two-edit path at 199

**FINITE-EXACT, all orders.** There are 347 distinct one-flip readings opposite in mod-four type to the starting height 199. Among them the complete list of square-sum neighbors, now without any restriction to the 72 prescribed labels, is
\[
  57,125,201,285,377,585,701,825,1097,1245,2505,3897,4157,6857 .
\]
The exact finite reduction is especially small for this prime. There is one nontrivial strong component, with 39 possible isotypes. Take every internal edit. If a singleton is an endpoint of an external edit and the other endpoint is in the core, test all core ports and either zero or one intervening singleton. With a nonempty middle, two or more intervening vertices would force gain divisible by four. If both endpoints are singleton and the core is outside, a single intervening singleton gives height \(3\cdot199\). If the core is inside, its two-slot mass is divisible by four. This exhausts both directions including losses from internal edits. [Complete core generator and port formula](PHYSICAL.md)

| Target after one flip from exactly 199 | Least binary precision, all moves | Other actual readings still possible one bit earlier | Least precision also requiring a square sum |
|---|---:|---|---:|
|201|10|713, 3785, 6857|10 (6857 survives at 9)|
|701|10|189, 2237|8 (4157 survives at 7)|

For example, an actual one-arc neighbor congruent to 701 modulo 256 is exactly 701 if its sum with 199 is square. This makes the role of the extra observable quantifiable. The two incident edges of the 72-cycle nonetheless require distinct cores.

A different pair of moves through 57 does compose. Start at the nine-vertex tournament with outgoing masks
\[
                   (510,508,128,228,172,388,436,256,28).
\]
Reverse the pair \(\{5,6\}\) and then \(\{0,3\}\), indexed from zero. The three exact counts are \(43,57,199\). Their adjacent sums are \(10^2,16^2\). The first numerical relation, read backwards, is the accelerated positive Collatz step \(57\to43\); the second square reflection provides a genuine continuation by a different operation. In the opposite edit order the counts are \(43,141,199\). The two arcs are disjoint on the labelled host. This extends the first contact in the [parallel signed-square decoder](../chanpak_signed_square36_20261123/README.md).

Run `python3 -B precision.py` in this directory, with the same `labelg` and C++ prerequisites as the all-order arc census.
