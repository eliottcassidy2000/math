# The one-factor-of-seven axis and an exact growing step of fourteen

**Proved:** a literal minimal-order one-arc switch 77→91; a synchronized two-contact switch giving exact gap 14 under every common continuation; and a selector placing such a pair within O_M(√x) of any sufficiently large prescribed quotient x, in any prescribed compatible residue class modulo a fixed M. A different port architecture has a necessary sevenfold scaling constraint and still misses known targets. No general prime placement or Chan–Pak proof is claimed.

## 1. What remains after deleting composite quotients

For odd H with v₇(H)=1, write H=7q with q odd and 7∤q. Consecutive odd q give a step 14 in H. The full lane is 7 mod 14, with the multiples of 49 removed. Requiring prime q removes the unit q=1 and the composite q; apart from q=3 (the hole 21), a prime pair q,q+2 must start at q≡−1 mod 3 and avoid q≡0,−2 mod 7. Thus the shortest step is 14, but even the first local sieve rejects some 14-steps.

Deleting composite quotients removes genuine graph constructions, not just candidates. Ordinal product multiplies counts. If q has any prime factor p other than 3 or 7, a witness for 7p and a witness for q/p suffice for 7q. Pure powers q=3^j are already covered for j≥2 by 63 and ordinal triangles; q=1,3 are exactly the known holes. Therefore, assuming the ordinary prime axis, the only new generators on the v₇=1 axis are 7p for odd primes p≠3,7. This is not a claim that a witness for a composite quotient can be divided by a factor.

More precisely a tournament of count 7q has exactly one strongly connected component whose count is divisible by 7. It has count 7d, d|q, with a 7-free ordinal prefix and suffix. If F_r(z) is the size enumerator of all classes at count r (allow the empty tournament at r=1), and S_r(z) that of strong classes, then

    F_(7q)(z) = Σ_(d|q) S_(7d)(z) Σ_(ab=q/d) F_a(z) F_b(z).

Each witness has a unique ordered decomposition of this form. If the whole class is self-converse, the unique marked component must be self-converse and the prefix and suffix are converses. Hence the self-converse size series is Σ_(a²|q) S^SC_(7q/a²)(z) F_a(z²). In particular when q is squarefree, all nonunit mass in a self-converse witness must be carried in the marked strong component.

## 2. A single arc realizes the first twin-prime gap

The following order-seven tournaments (outgoing bit masks) differ only on the arc between vertices 3 and 4, zero-based:

    (74,32,3,6,15,29,62)   H=77=7·11,
    (74,32,3,22,7,29,62)   H=91=7·13.

The exhaustive labelled replay below checks all tournaments through order seven and all their single arc flips. No one-arc +14 between strict v₇=1 counts occurs below order seven. At order seven there are 50,400 upward labelled flips in that channel, of which 20,160 have prime quotients at both ends. These are labelled flips, not an isomorphism-class count.

The two goals “one arc” and “repeatable” can also coexist at unbounded order, but prime access needs care. Start with rows (32,13,25,113,99,198,135,31), flip 2↔6, and replace vertex 0 by a transitive m-block. It has only one possible exit vertex in either core, so it is visited in at most two runs. The one-run read is unchanged by substitution, and the two-run read is multiplied by 2^m−2. The m=1,2 counts give exactly 246·2^m−95 and 246·2^m−81. Both are strict v₇=1 on m≡2,8,11,14,20 mod 21 (reduce 2^m mod 49). But the upper count is always divisible by 3, so the example misses the prime-pair problem despite preserving a literal single arc and the gap at all scales.

## 3. A common continuation keeping the full +14

Take the same nine-vertex deletion D=(262,268,280,49,35,7,319,383,56). Add a marked root with outgoing mask 183 or 159. That exchanges exactly two incidences: it stops beating vertex 5 and starts beating vertex 3. Define g=H(D). The marked E and O count the ordered covers by a path starting at the root followed by an optional second path, classified by the even or odd number of vertices in the first path. Literal path DP gives

    (H,E,O,g) = (1947,2698,2748,175) or (1961,2698,2748,175).

Why not make just one recontact? On a fixed deletion, E+O is the singleton contribution g plus a sum of positive terms indexed by root out-neighbors, according to the second vertex of the first path. Changing one root arc strictly changes E+O, so one recontact cannot preserve the whole continuation state. The raw pair above exchanges vertices with equal marked continuation counts.

Let N append a vertex that beats only the root, and S append an ordinary sink. The fixed-root binary extension rules read

    N: (H,E,O,g) → (H+2E,E+O−g,E+O,g),
    S: (H,E,O,g) → (H,E+O,E+O,g).

After one common N the two states are (7343,5271,5446,175) and (7357,5271,5446,175). Thus the prime quotients are 1049 and 1051. Both states have all four coordinates divisible by 7; common words preserve this, retain identical E,O,g, and hold the gap exactly 14. Both root masks beat vertex 7, an ordinal source in D, so the root reaches all deletion vertices. After N all deletion vertices reach the root via the new vertex, so both graphs are strong. An ending N restores strength after any appended ordinary sinks: they reach the root via the last vertex and are reachable from the root directly.

The inherited normalized binary count law (O=E+g) reads the continuation word with N=1, S=0 as its low-to-high binary address z≥1. It gives the lower quotient

    Q(z)=−482+1506(z+1)+25 F(z),
    F(0)=1, F(2t+1)=F(t), F(2t)=F(t)+t+(t mod 2).

The other read is exactly Q(z)+2. For the binary law, peel the first N or first block S^kN. A lone N with remaining address t gives z=2t+1 and doubles A=2E. A block S^kN with k≥1 gives z=2^k(2t+1), adds 2^k(A+g) to H, and sends A to 2^(k+1)(A+g)−2g. Induct using F(2^k(2t+1))=F(t)+2^k+(2^(k+1)−2)t. For instance iterating only N gives 753·2^t−457 and two more, t≥1 including the normalized starting pair. Those pure-N rails both have v₇=1 exactly for t≡0 or 1 mod 3. The general words below let the target choose much more than this pure phase.

## 4. Choose the exact phase and place the physical pair near the target

Here is a joint estimate, not just separate real and modular statements. For any fixed M ≥1 and any residue c modulo M compatible with oddness, every sufficiently large x admits a pair of strong tournaments with counts 7Q and 7(Q+2), where Q≡c (mod M) and |Q−x|≤C_M√x.

First, for even u, use the low-digit operation T_k(u)=2(2^k(u+1)−1), k≥1. Set v=u+1. The identities are v′=2^(k+1)v−1 and F(u′)=F(u)+2^k v. If M=2^a m with m odd, toggling k=1,…,a−1 programs the odd residue of F(u) modulo 2^a. All later k are at least max(1,a). Modulo m a sufficiently large k with 2^(k+1)≡1 makes v′=v−1; a k with 2^(k+1)≡2 fixes v=1 and adds 1 to F there. Move to v=1, adjust F, then move to the desired v while precompensating the known F drift. This makes v and F independently programmable modulo m. The result holds for m=1 by omitting the odd step.

For a finite menu of suffixes, take any L≥max(1,a), B=1506·2^L. Since gcd(1506,25)=1, the preceding argument gives r≥0 realizing every compatible class for −482+B(r+1)+25F(r) modulo M. Only finitely many B residues occur, so choose one suffix per compatible class and residue and pad all suffixes to a common bit width k. This width depends on M only. To choose near x, take L with √x≤B<2√x, and pick q, a nonnegative multiple of 2M, by rounding x/(B·2^k) down to such a multiple. Insert the chosen k low bits r below q and put u=2^k q+r. The binary recurrence gives

    F(u) ≡ F(q)+F(r)−1 (mod M).

Indeed deleting the low bits contributes an affine term in q plus a term in its parity; both vanish modulo M by the choice of q. Set z=2^L(u+1)−1. The final L low ones do not change F, so Q(z)=−482+B(u+1)+25F(u). Choose the suffix class by subtracting 25(F(q)−1) from the desired phase. This preserves the parity requirement. The error from rounding q is at most 2M B·2^k. The remaining error is bounded by B·2^k+25·2^k(q+1)+482, using F(u)≤u+1. All are O_M(√x).

In particular choose M divisible by 7 and avoid c≡0,−2 mod 7 to guarantee strict v₇=1 at both outputs. Adding any fixed finite prime sieve to M lets both quotients avoid every sieving prime simultaneously (choose c≡2 mod 3). This does not guarantee primes or exact absolute placement. The explicit compiler uses O(a²+m(a+m)) suffix bits, as the odd part needs at most three runs of m moves and each exponent is O(a+m). This even permits M≤c_ε√log x with error O(x^(1/2+ε)), by choosing c_ε small enough. It is nowhere near the modulus exceeding the absolute error.

## 5. A forced factor of seven in the movable Farey port chart

In the [earlier port chart](../least_hole_near_twins_20261117/README.md), N=2(2R+u)(S+v)−uv and the port-change gap is δ=uS−vR+2uv. Put K=uv and W=4vR−K+2δ. Eliminating S gives W²=13K²+(4N−12δ)K+4δ². For N=7p and δ=14, modulo 7 this is W²≡−K², so 7 divides both. If 7 divides just u or just v, the gap identity forces its matching cover R or S divisible by 7. If both u and v are divisible by 7, N≡4RS mod 7 forces 7 to divide at least one cover coordinate, but not both (then δ is divisible by 49). Hence **exactly one** of (u,R),(v,S) is coordinatewise divisible by 7.

Dividing that profile by 7 algebraically sends the pair to the +2 port equations at p. It does not necessarily give a physical normalized profile. It does give the same necessary conic bounds because each ratio P₂/H is unchanged and still at least 2. Thus the earlier integer conic exclusions for 10,007, 30,011 and 1,000,000,007 also exclude a +14 switch at their sevenfold targets in this specific port chart, at every order. The synchronized rooted construction above is not constrained to that chart.

Reproduce with `python replay.py`, then `g++ -O3 scan.cpp -o scan && ./scan`. The binary recurrence and seed data are checked by subset DP; the simultaneous residue/error result is the all-order argument above.
