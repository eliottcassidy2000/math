# Session35 arithmetic: 57, Leyland norms, and intrinsic sign fibres

Status: PROVED elementary statements below; FINITE-EXACT bounded receipt; no prime-distribution or major open-problem claim. No theorem IDs reserved. Worktree read at 1b91742ce55a. All work written outside the repository.

## Inheritance and concept board

Closest mechanisms: the base-2 Leyland divisibility sieve in `05-knowledge/results/special_primes_curling_genus_zero_langlands_20261015.md`, §3; actual balanced deletion in `05-knowledge/results/zero_heavy_ports_20261020g/microcosm.md`, §§2–4; the marked regular actions in `05-knowledge/results/metagraph_fibre_network_20261020k/README.md`, §3. Sparse-checkout missing files were read with `git show HEAD:path`, not treated as absent research.

Least-used relevant sidecar: **THM-3370, `01-canon/theorems/THM-3370-berggren-two-cube-biquadratic-norm-collision.md`**, especially its Q(sqrt(-2)) norm in (6). That theorem's Berggren parameter and cube representation are NOT transported here; its primitive quadratic-norm splitting mechanism is.

Hostiles: repeated-prime modulus9 has three roots of zero; nonsquare target3 has no roots modulo57; four roots alone do not specify a Klein-four operation; a Legendre character modulo a Blum composite becomes symmetric, not tournament-oriented. Corrected near miss: this note initially considered target3 for a two-root fibre, but the exact test rejects it; target9 is the repaired example.

Board: Leyland values / norm splitting / square-root sign torsors / literal deletion bits / dimension typing / composite Paley obstruction. Anchor: arithmetic structure of57. Niche: recovered quadratic-norm carrier. Wildcard: canonical principal-root marking gives an operation-level 4/2/1 bridge.

## 1. A norm filter for Leyland-Blum values

Let L(x,y)=x^y+y^x with integers x,y>1. A Blum integer means a product of two DISTINCT primes each congruent3 modulo4.

**Proposition A (necessary filters, PROVED).** If L(x,y) is Blum, x,y are coprime and of opposite parity. Write e for the even argument and o for the odd argument. Then e is not a square, and every prime p dividing L satisfies the Legendre condition (-e/p)=1.

Proof: equal parities make L even. If a prime divides both arguments, its square divides both terms, hence L is not squarefree. For coprime opposite-parity arguments,

L=e*(e^((o-1)/2))^2+(o^(e/2))^2.

Every prime p|L is coprime to e and o, so division by the first square gives (-e/p)=1. If e is a square, this implies (-1/p)=1, excluding either Blum prime. All statements are necessary only. For instance L(2,3)=17 passes coprimality/opposite parity/nonsquare-base tests but is prime.

Writing e=d*s^2 with d squarefree gives the literal quadratic norm

L = Norm_Q(sqrt(-d))/Q (o^(e/2) + s^o*d^((o-1)/2)*sqrt(-d)).

This is a genuine common carrier with THM-3370, not a shared numerical label. The map preserves value and prime splitting; it forgets the exponential constraints on the two norm coordinates. Those constraints must be retained for a converse. It imports no Pythagorean spine or cube decomposition.

For e=2*s^2, each Blum prime must be3 modulo8: (-e/p)=(-2/p)=1 permits p=1,3 modulo8, and p=3 modulo4 leaves3. The latter residue evaluation is the usual supplementary quadratic-reciprocity formula.

### Base2 refinement

For odd y>=3, N_y=2^y+y^2 is1 modulo8 and equals Norm(y+2^((y-1)/2)*sqrt(-2)).

If 3 does not divide y, then 3 divides N_y. In this subfamily:

**N_y is Blum iff N_y/3 is prime.**

Indeed N_y>=57, and N_y/3 is3 modulo8, hence is not3 and is a Blum prime when prime. Conversely a Blum N_y divisible3 must have this prime quotient. Modulo9, the odd residue classes y modulo18 give respectively

1,3,5,7,9,11,13,15,17 -> 3,8,3,6,8,0,0,8,6.

Consequently y=11 or13 modulo18 is excluded by a repeated factor3. In the divisible-by3 subfamily the surviving necessary classes are1,5,7,17 modulo18, and none is sufficient: y=19 gives3*179*977. For y divisible3 there is no asserted classification.

The first instance is57=2^5+5^2=3*19. Its norm representation is5^2+2*4^2; in Z[sqrt(-2)],

(1+sqrt(-2))*(1+3sqrt(-2))=-5+4sqrt(-2).

Thus its two Blum factors are themselves norms3 and19. Also Z[sqrt(-2)]/(5+4sqrt(-2)) is isomorphic to Z/57 by sqrt(-2)->13: 13^2+2=3*57, 5+4*13=57, and both quotient rings have cardinality57 (multiplication by5+4sqrt(-2) has determinant57).

## 2. Intrinsic square-root fibres, not arbitrary four-point labels

Let N=product p_i be odd and squarefree. Under CRT, S_N:a<-r^2 is the product of the local squaring maps. A target a has roots iff each local target is zero or a quadratic residue. If exactly s local targets are nonzero squares, its fibre has2^s points.

The operations are intrinsic: independently negate each local CRT root coordinate. These commuting involutions act transitively; a zero coordinate makes its sign flip trivial. The effective action is the regular F2^s action. This proves the fibre law and identifies the operations, not just its cardinality.

For Blum primes p_i=3 modulo4, a further canonical origin exists: b_i=a_i^((p_i+1)/4) is the unique root in QR(F_p_i), or zero if a_i=0. Its square is a_i by Euler's criterion; it lies in QR because a_i does. The other root is not QR since -1 is not QR. Therefore every root has unique signs relative to b_i on the active coordinates. No arbitrary ordering of four roots is used. The principal-root fact is classical: Menezes–van Oorschot–Vanstone, *Handbook of Applied Cryptography*, §§2.157–2.160, [primary chapter](https://cacr.uwaterloo.ca/hac/about/chap2.pdf). The zero-coordinate extension and the deletion comparison are the local arguments here, not attributed as new classical results.

For57 the CRT idempotents are19 (1 mod3,0 mod19) and39 (0 mod3,1 mod19). Its roots of1 are1,20,37,56; roots of9 are3,54; root of0 is0. For target1, ordered bits(00,10,01,11) map explicitly to(1,20,37,56). For target9 the unique principal root is54 (its mod19 coordinate16 is a square), and the other root is3. There are9 targets with four roots,10 with two,1 with one, and37 with none. Thus the square image has20 points and weighted source count4*9+2*10+1=57.

### Explicit deletion-fibre normal form

The inherited path decoder P(I,K) and central deletion retain the two ordered subset coordinates and forget the newest finite bit in each; infinity is fixed. Fix a parent path and a square target a modulo pq whose zero-coordinate pattern matches the parent's infinity pattern, with p assigned to the first coordinate and q to the second. Map a lift with newest bits epsilon_i to the unique CRT root whose coordinates are (-1)^epsilon_i*b_i; infinity coordinates map to0.

This map is a bijection and intertwines **the actual deleted-bit toggles** with **local root negations**. Its inverse reads the sign of the local root relative to the unique principal root, then applies the inherited literal path decoder. It is canonical relative to the specified coordinate assignment, target, and parent; there is no claimed canonical choice of target for every parent. Changing the target changes the map.

Preserved: the two labelled involutions, their commutation, stabilizer caused by boundary coordinates, and the effective regular action. Lost: parent path identity if not retained; path ports; tournament arcs; blue/black status; metagraph adjacency; arithmetic multiplication between different fibres. In particular, bit toggles are path-set operations, not established black moves. This is a common marked normal-form theorem, not a natural equivalence of entire metagraphs.

For primes1 modulo4 both signs can be QR, e.g. +/-2 modulo17, so the principal-QR-root marking fails. For nonsquarefree modulus9, zero has three roots: squarefree scope is necessary. Non-Blum squarefree semiprimes still have4/2/1 counts, but not this canonical QR marking.

## 3. Dimension arithmetic and the orientation boundary

For s active sign coordinates, the root fibre has2^s points; its functions over C form a2^s-dimensional representation of F2^s. The Walsh characters indexed by subsets of active coordinates give2^s one-dimensional eigenspaces. This is an exact representation dimension, while s is the torsor's binary-coordinate count. Its finite reduced algebraic-set dimension is0. These three meanings of dimension must not be conflated.

There is no field F57 because finite fields have prime-power cardinality. Z/57 is the product ring F3 x F19, not a vector space over a common finite field. An unrelated57-dimensional complex representation can of course exist, but cardinality57 does not identify it with this ring.

For Blum N=pq, Jacobi(-1/N)=(-1)*(-1)=+1. Thus Jacobi(y-x/N)=Jacobi(x-y/N); the composite Jacobi sign is symmetric on units and zero on nonunits. It cannot orient a Paley tournament. Using only the local Legendre character produces ties whenever that local difference vanishes. A lexicographic product of the two local Paley tournaments resolves those ties, but choosing priority is additional structure, not intrinsic to the scalar Jacobi character.

## 4. Exact experiment and stopping boundaries

Run `python3 session35-arithmetic/check.py` from the workspace root. Standard-library-only trial division with all primes through10^6 gives exact factorizations for every value<=10^12. The universe is all2<=x<y with L(x,y)<=10^12: y<=39 follows from2^40>10^12, and larger x cannot occur since x<y. There are111 pairs and12 Blum values. The receipt lists all factors. Equal arguments are excluded in the universe and cannot be Blum anyway because their value is even.

Independent test path: brute-force root bucket enumeration is compared with CRT-generated sign orbits for every residue of every product of distinct Blum primes below50, checking3800 nonempty fibres. It verifies local involutions, commutation, canonical QR origins and boundary multiplicity. Mod9 sieve tested for odd3<=y<2000. Hostiles are hard assertions. No probabilistic primality calls occur.

Loops completed: parity/coprimality candidate -> prime17 hostile -> necessary-only norm filter; base2 divisible3 candidate -> y11 repeated-prime and y19 three-prime hostiles -> exact quotient-primality iff; four-point-count analogy -> operation test -> canonical principal-root marking; composite Paley guess -> even character obstruction -> priority-marked lexicographic survivor.

Open: infinitude/density of Leyland-Blum values, distribution beyond the finite bound, and any metagraph connectivity consequence. No claimed LRC, Collatz, HN, Jacobian or other prize closure.

Independent audit: `session35-lie-arithmetic-audit` independently reproduced the111-pair/12-Blum census by direct trial division and checked literal path/deletion plus compatible mod57 CRT fibre maps for k=2..5. All substantive claims passed; source pin refined to HAC §§2.157–2.160.
