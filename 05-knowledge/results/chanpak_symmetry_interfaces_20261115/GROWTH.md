# Slot-polynomial growth and why two all-size wrappers do not cover targets

This note separates two uses of recursion. One can grow a full count with explicit functional operators, or look for a new tournament whose full profile can be used in place of an earlier one. Knowing the ordinary count alone is not enough for the second operation.

## 1. Exact cost of appending a vertex

Write P_T(k) for the optional ordered slot polynomial. If T'=T+v is any one-vertex extension, the [odd-cycle packing identity](https://arxiv.org/abs/2307.05569) gives the polynomial equation
```
P_T'(k) = k P_T(k) + 2 k sum_{C directed odd cycle through v} P_{T'−V(C)}(k).
```
The first term counts packings in which v is a singleton. In each other packing exactly one selected cycle contains v; remove it and its coefficient 2k. Cycles here have length at least three, up to directed rotation.

In particular for every integer k≥1 there is at least k-fold growth, and equality holds if and only if v belongs to a singleton strong component. Indeed if v lies on no directed triangle, all its in-neighbors beat all its out-neighbors, ruling out any cycle through v and therefore any nonsingleton component. Conversely a triangle gives a strictly positive summand at every k≥1 by Rédei path existence on its complement (including the empty complement). Taking k=1 recovers strict count growth at every addition of a vertex lying in a nontrivial strong component. Taking k=2 gives the stronger necessary condition Q(T')≥2Q(T) even if H is allowed to change.

There is also an elementary colored-path injection: select one of the k slots, insert v by the path insertion rule in that slot, and keep all others unchanged. The slot containing v and deletion recover the old cover. This proves the lower bound independent of symmetric functions, while the cycle argument identifies precisely when it is strict.

## 2. Neutral vertex padding cannot synchronize the two binary hosts

Use the matched binary states from [PORTS2.md](PORTS2.md): for a nonempty reservoir W, the strong sources F_4(W) and G(W) both have H=h+16q, for h=P_W(1), q=P_W(2). But F_4(W) has one more vertex. Writing s=P_W(3), their second readings are
```
P_F4(2) = 64q + 162s − 2h,
P_G(2)  = 44q + 98s − 22h.
```
Thus P_F4(2)−P_G(2)=20q+64s+20h>0, but
```
P_F4(2)−2P_G(2) = −24q−34s+42h < 0.
```
The inequality uses q≥2h and s≥3h by putting a full path in just one of the slots. So there is no one-vertex extension of G(W) meeting the F_4(W) slot state, even before asking that it leave H unchanged. Simply padding the shorter of the matched scalar sources cannot repair the continuation discrepancy. A change of the reservoir or of the actual reader is necessary.

## 3. Complete small-modulus behavior coexists with a power-sparse grammar

Let A(T)=C3[T,1,1] and B(T)=C3[T,TT2,1]. Start from the singleton, and apply these two wrappers alone in any order. Every such tournament is physical and strong. The [operator formulas](COMPOSITION.md) and [finite residue replay](verify_modular.py) show that this grammar hits all odd classes modulo 6, 10, 14 and 18 at every sufficiently large *order* (respectively n≥12,12,13,26).

Joint residues are more restrictive: modulo 30 at fixed order this same grammar ultimately misses one odd class, namely 25 at even order and 5 at odd order, despite full separate marginals modulo 3 and 5. The [norm invariant](COMPOSITION.md) proves the exclusions; the finite state replay proves eventual attainment of the other classes.

Nevertheless its distinct exact counts up to X are at most O(X^θ), θ=log₂ ρ≈0.406, where ρ is the positive root of ρ³=ρ+1. Proof: A and B add two and three vertices respectively. Their counts are H(T)+P_T(2) and 3P_T(2)−H(T), both at least P_T(2). Every assignment of T's vertices to two slots can be filled by directed paths, so P_T(2)≥2^|T|. A wrapper output of count ≤X must therefore have order ≤log₂ X+3. The number of wrapper words with total order cost n satisfies the recurrence with steps 2,3; the elementary induction bound O(ρ^n) and geometric summation prove the claim. No independence or distinctness assumption about their outputs is required. By the Euler-product argument for reciprocally summable banks in the [previous session](../chanpak_dyadic_collision_ports_20261115/AMPLIFICATION.md), arbitrary ordinal products of this particular wrapper bank still have density zero.

This is a physical example of why residue reachability, even with exact control for every sufficiently large order, does not replace absolute target placement. The promising operation is a selective port or a new conditional reader; iterating these two symmetric wrappers without changing the interface cannot close Chan–Pak.

## 4. Even changing both transitive sleeves at every step is sparse

Enlarge the grammar to arbitrary wrappers `C3[T,TT_a,TT_b]` with a,b≥1 at every step; start at a singleton. This contains the one-singleton infinite palette of [COMPOSITION.md](COMPOSITION.md), but no restriction from that palette's mod-15 argument is assumed here. The number of distinct attained counts up to X in the enlarged grammar is
```
O(X^log₂ ρ),     ρ=2cos(π/7)=1.801937...<2 .
```
In particular it still has density zero. This gives a quantitative reason to introduce a nontransitive sleeve, a changing marked reader, or a different quotient rather than merely enlarging the two transitive sleeves.

First, an output of order n has H≥2^n/12. Write m=|T|. When m≥2, the two-nonempty-slot count obeys `c_2(T)≥P_T(2)/3≥2^m/3`: cut each full Hamiltonian path in its m−1 interior positions (giving c_2≥(m−1)H), then use P_T(2)=2H+c_2 and the two-color assignment lower bound. With a,b≥2, count just the cyclic schedules with two runs in all three modules. There are three such schedules and c_2(TT_a)=2^a−2≥2^(a−1), so the count is at least 2^(n−2). With one sleeve singleton and the other of size b≥2, a cyclic word with two T and two TT_b runs and one singleton gives at least 2^(m+b)/6=2^n/12. If both sleeves are singleton, H=H(T)+P_T(2)≥2^m. If m=1, the literal cyclic schedules give H=(2^a−1)(2^b−1)+2≥2^(n−3). This proves the stated size bound in all cases.

Second, one need only count unordered pairs of sleeve sizes at each application. Reversing the quotient and all modules gives a converse with the two sleeves exchanged and T replaced by its converse. Path reversal in each slot preserves the full slot polynomial, and it is sufficient for substitution by the schedule formula. Consequently the two sleeve orders give the same full slot state, including all subsequent outputs. A size increment s=a+b has at most floor(s/2) choices. Let x solve
```
sum_(s≥2) floor(s/2) x^s = x²/[(1−x)²(1+x)] = 1,  0<x<1.
```
The number of recipes of total vertex cost j is O(ρ^j) for ρ=1/x, by the positive-step recurrence and induction. The equation gives ρ³−ρ²−2ρ+1=0, whose root here is 2cos(π/7): put z=exp(iπ/7) and expand z^3−z²+z−1+z^(−1)−z^(−2)+z^(−3)=0 in z+z^(−1). The size bound gives n≤log₂ X+O(1), completing the estimate. The 7 here comes from the three-term recurrence for unoriented sleeve choices; it is not an explanation of the forbidden Hamiltonian count 7.

The appearance of 7 has a precise walk interpretation. Let M be the adjacency matrix of a three-vertex path with a loop at its last vertex, i.e. rows (0,1,0), (1,0,1), (0,1,1). This is the reflection-folded path on six vertices. The exact recipe generating function is
```
1/(1 − x²/[(1−x)(1−x²)])
   = 1 + x²/(1−x−2x²+x³)
   = 1 + sum_(j≥0) (M^j)_(1,3) x^j.
```
The last identity is the (1,3) entry of (I−xM)^(-1) by cofactors. Thus the exponential recipe growth is the Perron growth of walks on the folded path, not an assumption about random tournaments. On the tournament side, for a fixed child T and size increment a+b, the first non-leading coefficient of the slot polynomial already distinguishes the unordered sleeve sizes: the number of directed triangles increases by |T|ab (exactly the triples using all three modules), and the coefficient of k^(n−2) is twice the triangle count by odd-cycle packing. This does not identify the child of an arbitrary output or establish distinctness across different recipes.
