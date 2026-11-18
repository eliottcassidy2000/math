# Three, five, twenty-seven: from modularity history to a testable tournament reduction

The useful analogy is not that the integers recur. It is that a coarse invariant can stay fixed while the behavior at one local place changes, and a theorem for a restricted class can suffice if all remaining cases can be brought into it. The elliptic and tournament invariants are **different objects**.

## 1. The historical distinction and a literal Weierstrass bridge

[Wiles's semistable modularity theorem](https://wstein.org/edu/Spring2003/21n/papers/Wiles%2CModular_Elliptic_Curves_and_Fermats_Last_Theorem.pdf), with the Taylor–Wiles repair and Ribet's level-lowering implication, already settled FLT. The 3–5 switch is inside that proof: when the residual 3-representation cannot start the argument, replace the curve by one with the same mod-5 representation but whose mod-3 representation can start it. The later results concern **extension of modularity**: Diamond allowed semistability just at 3 and 5 (see the early account in [Conrad](https://math.stanford.edu/~conrad/papers/ramdef.pdf)); [Conrad–Diamond–Taylor](https://math.stanford.edu/~conrad/papers/cdtmaster.pdf) treated conductor not divisible by 27; [Breuil–Conrad–Diamond–Taylor](https://math.stanford.edu/~conrad/papers/tswfinal.pdf) finished all rational elliptic curves. These are not successive prerequisites to FLT. The elliptic conductor's exponent at 3 is neither the discriminant valuation nor the valuation of a tournament count.

There is a genuine 27 here, not just a numerical association. [Li–Long–Tu, §2.4 and Remark 4.2](https://emis.de/ft/24022) identify the Fermat cubic U³+V³=1 with X₀(27). For s=U+V, d=U−V, the identity s³+3sd²=4 gives the birational transformation

    X=12/s, Y=36d/s,       Y²=X³−432.

Writing X=4x,Y=8y+4 gives y²+y=x³−7, the [level-27 optimal curve](https://www.lmfdb.org/EllipticCurve/Q/27/a/3). The database records conductor exponent 3 at 3, discriminant exponent 9 there, and no other bad prime. The [Weierstrass equation](https://dlmf.nist.gov/23.3) gives X=℘(z), Y=℘′(z)/2 for Eisenstein-lattice invariants g₂=0,g₃=1728; if ω³=1 then ℘(ωz)=ω℘(z), ℘′(ωz)=℘′(z). In characteristic 0 the Y projection forgets a three-sheet rotation; its branch values are the two roots of Y²+432 and infinity. The integral model reduced at 3 becomes the cusp w²=x³ after y=1+w. On its normalization x=t²,w=t³ the projected coordinate is a cube, inseparable in characteristic 3. At 5 the same integral cubic is smooth (partial derivatives and the equation have no common zero). This is a precise example of different local behavior at 3 and 5, not an identification of graph sheets with ramification.

## 2. Archive matches and a corrected symmetry boundary

- [Pythagorean Frey / Belyi THM-4839](../../../01-canon/theorems/THM-4839-unit-cayley-graphs-toroidal-maps-petersen-klein-four-belyi-pythagorean-frey.md) links Euclid to the Klein four cover X(4)→X(2). It is a marked torsion construction, not an unrestricted graph-to-elliptic map. Its Gaussian/Eisenstein unit Cayley graphs suggest a literal orientation comparison: multiplication by i on F5 is a 4-cycle plus one fixed point, so its square swaps endpoints of an edge; it cannot preserve any tournament. Multiplication by 2 on F7 is the Eisenstein order-three action, with 1+3+3 vertices and seven size-three edge orbits. Exactly 128 labeled orientations are possible, all Hamilton counts divisible by 3 because the rotation is free on total path orders. The replay sees counts 9,45,99,129,135,171,189. Requiring translations as well leaves just the two Paley orientations: nonzero differences must use {1,2,4} or its negative. Both have H=189=7·27. The factor 27 does not follow from rotation alone.
- The [Legendre pencil normalizer](../legendre_pencil_normalizer_20261011.md) explicitly shows that matching j forgets the relative differential; the [square-sum elliptic pencil](../square_sum_symmetric_elliptic_contact_pencil_20261019.md) gets consecutive elliptic contact but not the continuing interior port. These are close analogues for the difference between a count match and a reusable profile.
- [THM-2517](../../../01-canon/theorems/THM-2517-cubic-anchored-spectrum-flt3-and-three-time-boolean-lift.md) already uses FLT(3) directly to expose mixed modes by a cubic measurement, but this does not furnish a common physical clock. The [V4/C3 reflection](../../../07-reflections/modular-v4-c3-two-prime-point-direction-decoder-codex-20260801.md) contrasts Walsh Smith factors 1,2,2,4 and a 27 Fourier determinant norm; inverse arithmetic demands an extra coordinate. Neither is a conductor theorem.
- A caveat to the older THM-1810 order-four analogy: an equianharmonic four-point configuration admits an order-three permutation with cycle type 3+1. The compatible tournaments are precisely a cyclic triangle with an ordinal singleton (two converse classes, H=3). The strongly connected order-four H=5 type has no order-three symmetry, as its score classes are two pairs; indeed its automorphism group is trivial. Full A4 cannot preserve any tournament orientation, since an involution interchanges endpoints of an edge. Likewise the smooth equianharmonic model above is not the singular cubic y²=x³.

## 3. The exact local observer on the tournament side

For a core T write C_k for its ordered **nonempty** path-cover counts, and define

    h = H(T)=C_1,   raw = gcd C_k,   g = gcd (k−1)! C_k.

[Session293](../session293_cycle_returns_20261113/control/MODULE_IDEAL.md) proves that g is the **exact gcd over all uniform exteriors**. Excursions can be permuted; their factorial is an arithmetic tariff independent of automorphisms of T. It is not an elliptic conductor. The local feature analogous to keeping a second coordinate is the polynomial f(x)=Σ C_k binom(x,k). Odd-cycle packings express it as an integral monic parity polynomial (each packed odd cycle contracts to one atom), f(0)=0, f(1)=h, f(2)=2h+C₂. Newton inversion with factorial weights gives g=gcd_j (j−1)! f(j). This latter weighted basis change is already in the [27 reset ladder](../eight_nine_twentyseven_20261112/content/REPAIR_LADDER.md).

**Jet refinement of [Session294's prime-cube observer](../session294_context_27_20261113/README.md).** For odd p, p² or p³ divides g iff the corresponding precision holds at f(1),…,f(p). Let R=(f mod p)/(x^p−x)^2 when p²|g; this quotient is a polynomial. If p³|g, then for a=1,…,p−1,

    f(a+pt) ≡ p² R(a) t(t+1) mod p³;  f(pt)≡0 mod p³.

In particular the raw content reaches p³ iff (x^p−x)^3 divides f mod p, and a nonzero R(a) tells exactly which local addresses leak. At p=3 only one t mod 3 leaks in a given unit class; at p=5 there are three potential leak addresses. This is an exact higher-digit observer, not an assertion that the formal polynomial is physically realized.

Proof. Necessity tests the first p weighted evaluations, whose factorial weights are units. Conversely the p² tests and parity give f(a),f(a−p)≡0 mod p² for each nonzero a, and f(0),f(p)≡0. Taylor expansion kills the first derivative mod p on all residue classes; thus (x^p−x)^2 divides the reduced polynomial, and all evaluations are divisible by p². The indices j>p carry at least one additional p in (j−1)!, so the extra first-p precision suffices for p³. At the p³ precision, use f(a−p)=0 mod p³ in Taylor to infer p f′(a)≡p² f″(a)/2, and then expand f(a+pt). The reduced quadratic Taylor coefficient is R(a) since (x^p−x)′=−1 in F_p. At zero, parity makes either the linear or quadratic term vanish identically, and the test at p forces the other term to vanish mod p³. The same reasoning gives R(0)=0. The claim about the raw gcd follows from unweighted Newton inversion and the fact t(t+1) is nonzero for some t when p is odd.

For the actual H513 core, f reduces mod 3 to x²(x³−x)², so R(±1)=1. That explains the older residual f=18 on x≡±4 mod 9, and it is the missing bit between raw 9 and uniform-context 27. The monic formal polynomial x∏_(j=1)^p(x²−j²) shows the degree 2p+1 bound for this leak cannot be improved using just parity/precision: R=x, and f(p+1) exhibits valuation 2. The old physical survey rules out the formal order-seven bound at p=3 and witnesses the leak first at eight. At p=5 this formal bound is eleven, not a tournament example.

## 4. Literal same-count switches at 5, 27, and 7

The contrast between closed count and open response is directly observable:

| Closed h | Orders | Switch in uniform-context divisor g | Mechanism |
|---:|---|---|---|
| 15 | 5→5 | 15→3 | two distinct first-order models; C₂=90 vs 78, so the factor 5 disappears |
| 513 | 8→8 | 27→3 | inherited native switch; raw 9→3, so an excursion factor also disappears |
| 4151=7·593 | 10→10 | 7→1 | replace the cyclic seven-sheet witness by a different model, C₂=79954 vs 84802 |
| 16079=7·2297 | 11 | g=1 | new primitive model already at the first possible order, versus the paired presentation at 12 |

The two 4151 witnesses are both first-order since the exhaustively certified order-nine maximum is 3357. For 16079 use the [audited](../first_creator_generation_20261108/audit/AUDIT.md) order-ten maximum 15745. In the 4151 row the original C7 action ensures divisibility by 7 in *all* uniform exteriors; the alternate C₂ is nonzero mod 7, so a two-slot read detects the difference. The explicit new masks and full cover checks are printed by replay.py. The result preserves the scalar count, **not** the Paley sheet action, the original wrapper, or its continuation ports. These are not new unattained counts. They are a physical normalization at the same first order and a specific test for a subsequent wrapper.

A general reason the target-prime case is comparatively tame at this observer: for a strong H=p witness of a prime p≥5, Busch implies n<p (the minimum for n≥5 exceeds n), and g divides both H and (n−1)!n!. Hence g=1. For a strong H=7p witness with p≥11, n<p by the stronger M(n)>7n for n≥11 (verify bases 11,12,13 and multiply by five every three). Thus g is 1 or 7, and the only place to clear is 7. In contrast 3|H always implies 3|g: parity makes the reduced f vanish at zero and ±1. A reusable count-preserving normalization at 7 for all 7p witnesses is still unproved.

## 5. A reduction of the target theorem, not merely a numerical parallel

One can restrict Chan–Pak to **odd targets coprime to 15**, excluding 7, without losing any cases. Strip 3^a5^b by ordinal multiplication (C3 gives 3 and the nontransitive four-vertex tournament gives 5). The remaining unit m is covered by the restricted theorem unless m=7. If b≥1 in that exception, use the known 35 witness and reattach the remaining factors. If b=0 and a≥2, use the known 63 witness and reattach 3^(a−2). The two remaining cases are exactly 7 and 21. This is a rigorous “restricted class suffices” bridge. It also means covering all odd N with 27∤N (except the holes) would suffice, but is a weaker restriction than needed. Neither reformulation proves the required coprime core coverage.

The next test is now concrete: find a change of source at a prescribed 7p that keeps the needed continuation ports **as well as** the scalar count, while clearing the residual 7 factor in the uniform reader. At levels 3, 5, 7, and 27, merely matching H cannot certify that.

Reproduce with `python replay.py | diff -u replay.out -`; `g++ -std=c++17 -O3 seek.cpp -o /tmp/seek` reproduces the deterministic searches for `/tmp/seek 10 4151` and `/tmp/seek 11 16079`.
