# Two readings, three exit phases: curvature and physical controllers

**PROVED scoped / FINITE-EXACT; absolute inverse outside named banks OPEN.** This extends the inverse [cofactor test](../chanpak_cofactor_contradiction_20261123/README.md) and its [uniform persistence obstruction](../chanpak_cofactor_contradiction_20261123/UNIFORM-GUARDS.md). Work with the same literal [factor-seven exterior and two marked reversals](../chanpak_ternary_buckets_20261122/PHYSICAL-FORKS.md).

Write \(D=H(M)\), \(Q=Q_2(M)\) for a nonempty filler, \(h=7(D+40Q)\), and \(C=3h=21D+840Q\). Then the actual center and exit counts are \(C,C\pm2D\). The observations have different jobs: \(D\) sets the exit gap, and \(Q\) changes the common center independently whenever a controller can vary it at fixed \(D\).

## 1. The real geometry and the modular geometry are different

Ordinal attachment of \(V\) is the exact action \((D,Q)\mapsto(D_VD,Q_VQ)\). Each slot restricts to a slot in each ordinal block, and concatenating the restrictions is unique. In coordinates \(G=2D,C=21D+840Q\) the quadratic form \(G(C-21G/2)\) is proportional to \(DQ\); every attachment preserves it up to a positive scalar. Up to scaling this is the **only nondegenerate common conformal quadratic form**: even a single transitive vertex acts on \(D^2,DQ,Q^2\) with distinct eigenvalues \(1,2,4\); the two square forms are degenerate. After determinant normalization every nonempty ordinal attachment is a genuine Lorentz boost, of rapidity \(\frac12\log(Q_V/D_V)\) up to choice of axis. Since \(Q_V\ge(|V|+1)D_V\) by splitting each Hamilton path into two ordered slots in \(|V|+1\) ways, the projective ratio
\[
R=(C-21G/2)/G=420Q/D
\]
strictly grows under every nonempty attachment. In the \(\log R\) coordinate this is an exact positive translation. A modular inverse may be synthesized by repeating blocks; it is **never** an inverse that returns the exact real ratio.

A second geometry comes from writing the marked exit as \(T:(C,G)\mapsto(C+G,G)\), a parabolic shear. On projective lines adjoin the **read-chart change** \(S:(C,G)\mapsto(-G,C)\). Matrices obey \(S^2=(ST)^3=-I\), and \(T^m=I\bmod m\). Thus their projective action modulo \(m\) is a quotient of the orientation-preserving \((2,3,m)\) triangle group. This is the precise use of the [previous triangle comparison](../hidden_compositions_5_6_7_eleven_20261113/README.md): \(m=5,6,7\) yields the spherical, Euclidean, and hyperbolic source groups, respectively, from excess \(1/m-1/6\). The comparison \(m=7,8,9\) instead stays hyperbolic, with deficits \(1/42,1/24,1/18\). The finite projective action does not make the source group finite. This gives a second explanation for why looking only at periodic residues cannot detect real projective drift.

These chart generators act transitively on primitive projective lines modulo \(m\) (primitive pairs taken up to multiplication by a unit). Indeed, given a primitive column \((x,y)\), use upper shears so that \(x+ky\) is a unit. For each prime divisor of \(m\) at most one choice of \(k\) is excluded (none when \(y=0\bmod p\)); combine choices by CRT. The lower shear, conjugated from \(T\) by \(S\), then kills \(y\). **Only the one-step shears out of the marked center are supplied by the fork.** The formal swap \(S\) is not an arc reversal or a proof that an exit is a new fork center. The earlier [supports \(\{5,6,7\}\), \(\{7,8,9\}\)](../session64_lca_reverse_valuation_20261026/session64-repairs/next11/README.md) are marked endpoint neighborhoods with different repair failures, not evidence of literal continued shears. The present analogue of those support repairs is to keep the two coordinates rather than just a center phase.

The native physical reductions at those moduli are much smaller than the abstract projective action. Exactly, for a filler of order \(n\):
\[
\begin{array}{c|c}
m &(C-2D,C,C+2D)\bmod m\\\hline
5 &(-D,D,3D)\\
6 &(D,3D,5D)\\
7 &(-2D,0,2D)\\
8 &(3D,5D,7D)\\
9 &(B-2D,B,B+2D),\quad B=6D\ (n\ {\rm even}),\quad B=0\ (n\ {\rm odd})
\end{array}
\]
For 9 use \(Q_2\equiv(-1)^nD\bmod3\), by the parity of the slot polynomial. At 7 the center is forced to zero while the exits need not be zero. At 5, 6, 7 and 8 a pure change of \(Q\) cannot change any fork residue. Thus crossing the triangle-group curvature threshold and having two effective physical modular actuators are different phenomena.

## 2. Constructive modular inverse with a variable block

Let \(N_k\) be the tournament obtained by adding a single last-to-first reversal to the transitive tournament of order \(k+2\). The relevant pair is
\[
(D_k,Q_k)=(2^k+1,\ 4\cdot 2^k+2\cdot3^k-2).
\]
For the first reading, a path either does not use the reversed arc (one path) or chooses an arbitrary subset of the interior vertices before that arc. For the second, two empty-allowed ordered slots may put both ends in different slots (two choices, each interior vertex chooses a slot); or in the same slot (two choices). In the same slot there is an increasing version for each binary distribution of the interior vertices, and the reversed arc version requires that any interior vertex lie before the terminal of that arc, after its initial, or in the other slot. Exclude the all-other-slot assignment from the increasing version because the reversed arc forbids the direct increasing adjacency of the two endpoints. Summing gives the displayed formula.

For an odd modulus \(m\) prime to 6, set \(e=\varphi(m)\). Then \(N_e\) has \((D,Q)\equiv(2,4)\bmod m\); a singleton has \((1,2)\). Starting from a seed \((D_0,Q_0)\), form an ordinal sum with \(a\) copies of \(N_e\) and \(b\) singleton pads. The result has
\[
(D,Q)\equiv(D_0\,2^a,\ Q_0\,2^{2a+b})\pmod m.
\]
When 2 generates the units and both seed reads are units this gives every target unit pair by two discrete logs. Take the least nonnegative exponents (add \(e\) singletons if otherwise empty). The order cost for an empty seed is at most \((e-1)(e+3)\) for a nonempty recipe, or \(e\) for the empty one replaced by padding. To set a specified center \(c\) and gap \(g\) when \((m,840)=1\), solve \(D=g/2,\ Q=(c-21D)/840\) in the ring first; the physical exits then have residues \(c\pm g\).

For each \(p\) below 2 generates the units and there are also axis seeds in the *same* family. The exponents in the middle column give \(D=0,Q\ne0\bmod p\), and those on the right give \(D\ne0,Q=0\bmod p\).
\[
\begin{array}{r|rr}
p &D=0&Q=0\\\hline
11&5&7\\13&6&8\\19&9&13\\29&14&9\\37&18&7\\59&29&31\\61&30&32
\end{array}
\]
Scaling each axis seed with the same controller gives every point on that axis. Joining the two seeds gives the origin. This covers the **entire residue plane**, including any \(c,g\) modulo these primes. Nothing in this modular construction requires divisibility to persist on an exit. It gives no bound making a particular lift equal to a prescribed absolute height.

## 3. An exact integer inverse on an enlarged bank

Unlike the modular construction, the following solves its own exact target problem. Observed strong or nonempty block profiles are:
\[
\begin{array}{r|r|l}
D& \log_2 Q & \text{outgoing rows}\\\hline
5&5&(2,4,9,3)\\
19&8&(38,4,8,19,39,14)\\
217&12&(70,60,24,17,33,77,158,63)\\
223&12&(230,28,8,17,37,78,158,62)\\
599&13&(178,100,137,19,38,76,157,58)
\end{array}
\]
Rows use the lowest bit for vertex 0. The last three were found in the anchored order-eight chart. Ordinally combine copies of these blocks and singleton padding: then \(Q=2^e\) and \(D\) is a product of their five displayed bases, with \(e\) at least the weighted sum of the base exponents. For an absolute prescribed cofactor \(z\), enumerate \(e\) with \(40\cdot2^e<z\), compute \(D=z-40\cdot2^e\), and test by repeated division by the fixed bases. They are pairwise coprime, so greedy division determines every multiplicity uniquely; the weighted remainder is the required number of singletons. The identical procedure for a prescribed **absolute exit or center** \(z\) enumerates \(840\cdot2^e<z\) and sets \(D=(z-840\cdot2^e)/(21+2k)\) for \(k=-1,0,1\) respectively. Reject a noninteger or nonpositive \(D\), and reject an empty recipe. This performs \(O(\log z)\) outer trials, without an integer factorization oracle.

An example at \(D=223,e=12\) gives the cofactor \(164063\), exact parent height \(1148441\), and the exits \(3444877,3445769\) around center \(3445323\). The two exits are actual single-arc neighbors of this center. Additional exact prime cofactors in the enlarged bank are 164057, 655583, and 5243479, verified by trial division and the exact inverse. The fixed bank's output count through \(X\) is only \(O((\log X)^6)\), bounding separately the five multiplicities and \(e\); the algorithm cannot cover all targets. Also \(D=1\) is required for a unit fork and reduces the bank to pure singleton padding. The center phase controller and this exact bank therefore solve different problems.

See also [the exact obstruction to uniform multiplication by seven](UNIFORM-AMPLIFIER.md).

See [fixed-order, one-arc control](PHYSICAL-TOGGLES.md) for marked generators up to 100, and [prime-power ordinal control](ORDINAL-REAGENTS.md) for a larger certified range. Replay the arithmetic and targeted inverse with `python3 fork-phase.py`.

**Renewal boundary.** [Matching depth and the exact observation lattice](REENTRY-DEPTH.md) predict when two readings cease to determine the next fork. [Restricted renewal](RESTRICTED-RENEWAL.md) restores an exact two-read decoder on two wrapper towers.
