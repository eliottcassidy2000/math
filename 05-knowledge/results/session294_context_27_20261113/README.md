# Session 294: the first 3-adic context jump, and its amplification

**Proved in notes / finite-exact census. Chan–Pak remains open.**

Session 293 showed that a uniform module's universal context gcd is
I=gcd_k((k−1)! F_k), where F_k=k! pc_k counts ordered path covers.
The extra factorial comes from freely permuting the k−1 exterior
excursions. A triangular family of physical exterior probes establishes
exactness. This session asks what the numbers 8, 9 and 27 mean structurally.

## A short polynomial and an odd-prime observer

For a nonempty n-vertex tournament define f(x) by summing
2^(number of cycles) x^(n−sum(length−1)) over disjoint packings of
directed odd cycles, including the empty packing. Apply the inherited
odd-cycle Hamilton count identity (THM-002) in every block of a
vertex partition into k nonempty blocks. After contracting the occupied
cycles to atoms, the number of such partitions is the Stirling number
S(j,k), where j is the exponent of x. Therefore pc_k is the Stirling
transform of f and F_k=Δ^k f(0). This refines the existing path-cover
and endpoint-jet perspective (THM-1975, THM-3121, THM-3134). It also gives
f(1)=H and f multiplies under ordinal sums.

**For any odd prime p and integral parity polynomial f with f(0)=0,
p³|I iff f(1),...,f(p) are all zero modulo p³.** The factorial weights
for k≤p are units, and triangular finite-difference inversion proves
necessity. For sufficiency, f(a−p)=±f(p−a) for 1≤a<p; f(0)=f(p)=0.
Taylor modulo p² forces f′(a)=0 modulo p for each residue a. Hence all
integer evaluations vanish modulo p². Terms k>p gain the missing factor
from the excursion factorial; terms k≤p already vanish modulo p³.

## First birth at eight

For n≤6 there is no abstract monic parity polynomial of degree n
with f(0)=0 passing f(1),f(2),f(3) modulo 27 (substitute all lower
coefficients modulo 27). For n=7 put f=x^7+A x^5+B x³+C x.
The conditions reduce to A≡B≡4 mod9 and 1+A+B+C≡0 mod27.
Here A=2·triangles, B=2·5-cycles+4·disjoint triangle pairs, and
C=2·spanning 7-cycles. A complete direct replay of 2²¹ labelled
orientations gives 266,000 passing A, 5,600 passing A and B, zero
passing all three. An independent path-cover DP on the inherited
Session24 unlabelled catalog corroborates this.

For n=8 put f=x^8+A x^6+B x⁴+C x². The conditions are f(1)≡0 mod27,
A+2B≡3 mod9, C≡0 mod3; they follow from f(2)−4f(1) and f(3).
Path-cover DP on the complete Session24 catalog gives exactly 14
classes, all raw gcd 9 and actual context gcd 27:

| (A,B,C) | H | Classes | Automorphism orders |
|---|---:|---:|---|
| (22,76,36) | 135 | 2 | 1,1 |
| (34,214,210) | 459 | 8 | all 1 |
| (34,232,246) | 513 | 2 | 1,1 |
| (40,292,288) | 621 | 2 | 3,9 |

The four different polynomials give exactly four full uniform-context
profiles: Stirling inversion and the exterior probes are both
triangular. But all four have exactly the same residual modulo 27:
f(x)=18 on x≡±4 mod9 and zero otherwise. For multiples of 3, use
C≡0 mod3. For units write ±x=1+3j. Taylor and f(1)=f(−2)=0 mod27,
f″(1)≡2 mod3 give f(x)≡9j(j+1). This shows exactly where the unscaled
cover evaluator leaks a ternary digit and why the first such leak
(at Δ⁴) is invisible in actual uniform contexts (3! supplies a factor).

## Infinite residue complementarity

For an ordinal sum of r exceptional blocks of any of the four types
and s directed triangles, r≥1, the exact valuations are
v₃(raw)=min(2r+2s,3r+s), v₃(I)=min(2r+2s+1,3r+s), and v₃(H)=3r+s.
The product polynomial has valuation at least 2r+2s on the residual
pair x≡±4 mod9, with equality at 4; it has valuation at least 3r+s
elsewhere, with equality at 1. This uses x³+2x divisible by 3
everywhere and by 9 on the residual pair. Raw differences have
the same gcd as integer evaluations. If the residual branch strictly
wins, the fourth difference attains it (the earlier evaluations have
higher valuation), and all k≥4 gain at least one. Otherwise H itself
attains the minimum. Thus with s=0, the hidden gap between scalar
and full context valuations is r−1, arbitrarily large even when all
blocks are asymmetric. Three ordinary directed triangles, at order
nine, instead give raw and context gcd both 27. Two ways to buy a
factor of 27, one through independent internal factors and one through
a missing factor provided only by exterior excursions.

For the H=513 core with rows (56,117,129,22,100,76,141,59), cyclic
substitution of three copies yields H=163567386770762025 (v₃=9)
and raw/context gcd both 3⁸, by the THM-3121 run-word DP. The ordinal
sum also has v₃(H)=9 but raw/context gcds 3⁶ and 3⁷. Independent
local excursion factorials need not multiply: separately they divide
the same exterior coefficient, so the lcm is always justified but
the product need not be.

## Square-sum and Collatz: exact bridges with different carriers

For bidirected square-sum G8=P4⊔P2⊔P2 and G9=P4⊔P3⊔P2 the same
uniform-context argument applies, but the exact raw/context gcds are
48 and 96 at both sizes. In any bidirected linear forest of c
components the first nonzero coefficient pc_c=2^(number of
nonsingleton components) is an odd-prime unit. Consequently the odd
prime gcd valuations are v_p(c!) and v_p(c!(c−1)!). The induced graph
on {1,2,4,6,8,9,10,11} has edges 1–8, 6–10 and no others, hence
pc(z)=z⁶(2+z)² and exact gcds 2880, 345600: the same 3-adic jump
9→27 at order eight, but caused by absence of low-run terms rather
than divisibility of nonzero terms. Nine-scaling gives c disjoint
square-sum edges using labels {2·9^j,14·9^j : 0≤j<c}; cross-layer
sums have nonsquare unit 2 mod3. Multiples of 9 induce a rescaled
square-sum graph, while the nonzero ±1 classes among labels divided
by 3 give a bipartite sector inside the multiples of 3.
Abstract uniform contexts are not literal square-sum attachments:
labels 1,8 have no distinct common neighbor (factor their difference
of squares, 7). This is the same capacity issue as the older
square_sum_common_center_capacity_20261019 result.

Session53 (session53_g27_collatz_cooper_20261024) provides a useful
opposite regime: G27 attaches to 9,22, giving three endpoint extensions
and 32 joins from 1,116 two-covers. Selecting specific endpoints can
break the uniform excursion action. The same session proves the
long-descent Collatz cylinder 27+2⁵⁹t→23+3³⁷t. For any fixed shortcut
parity word of length k with o odd steps, x₀+2^k t→y₀+3^o t; reducing
the free t modulo 3^s loses exactly min(s,o) ternary digits, with
fibres all of size 3^min(s,o). The short prefix 27→41→62→31→47
gives 27+16t→47+27t. Compare the branch-guarded 3-adic tournament
recursion of THM-4754, not a claim about arbitrary Collatz iterates.

A separate exact clock: 8=2³=−1+9, ord₉(8)=2 and ord₂₇(8)=6.
By parity f(8^j) modulo 27 reads f(1) and 9 f′(1) only. All six
clock reads vanish for the order-eight types, yet f(4)=18.
THM-4553's eight-point PSL(2,7) → seven-point Paley → six-point
torus carries genuine point-stabilizer group actions, in contrast to
a Collatz branch cylinder or an exterior excursion action. These
are three different reasons a factor of 3 appears.

## What this suggests for Chan–Pak

A uniform module can transmit a prime-power restriction without an
automorphism. Splitting by endpoint incidence or changing the core
is necessary to evade such an obstruction when it excludes the
target. A next experiment is to use a bounded number of selected
endpoint colors to erase the growing restriction in the ordinal
families while keeping a renewable interface. Content 1 alone is a
gcd statement, not positive representation by a single physical
context, so this is not a coverage proof.

Arithmetic and order-seven proof independently reviewed. The order-eight
class claim uses the completeness of the inherited Session24 catalog.
