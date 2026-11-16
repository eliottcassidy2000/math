# Seven, eight, nine: signed quotients and a raw-gap obstruction

This note does **not** identify a modular-curve genus with a Hamilton
path count, or prove Collatz convergence. It identifies places where
sign is retained or lost and gives independent exact small controls.

## The seven-point arrow table and the genuine missing sign

Use the cyclic Fano convention on residues modulo seven: for each j, the
product of e_j and e_(j+1) is e_(j+3), with cyclic permutations giving the
same sign and reversed order giving the negative sign. These seven
oriented triples are associative lines; the algebra is the real
octonions. Retain the archive's distinction between this single Fano
plane (seven triples), the Paley tournament's whole directed triangle
set (fourteen), its arcs (twenty-one), and Hamilton count (189).
See the [Paley-Fano audit](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_mod6_20260922_paley_fano_octonion_design.md)
and [Sévennec's primary paper](https://cml.centre-mersenne.org/item/10.5802/cml.9.pdf).

On forgetting signs, the imaginary axes have the Fano index set
F₂³ minus 0. Every unsigned Fano linear permutation has exactly eight
signed lifts preserving the fixed multiplication table; the kernel
consists of the eight sign characters of F₂³, since invariance forces
s(x+y)=s(x)s(y). The entire signed group is the order-1344 nonsplit
extension described by [Koca, Koc and Koca](https://arxiv.org/abs/1612.06107).
The companion exact enumerator verifies the lift count without needing
to infer it from the order. The purely unsigned, multiplication-preserving
subgroup is exactly the order-21 Paley affine group: the signs of e_x e_y
literally orient the pair (x,y), so a pure permutation must preserve
the tournament, and the exact check confirms that only the 21 affine
maps with square slope do. In particular an index-eight quotient does
not make eight extra unsigned arrow-preserving transformations.

A local proof of the nonsplitting obstruction is convenient here.
Fix a nonzero c in F₂³. The maps
`T_l(x)=x+l(x)c` for the four linear functionals with l(c)=0 form a
Klein four-group. Any order-two signed lift of a nonidentity T_l must
negate the fixed e_c: choose x with l(x)=1, so T_l swaps x and x+c.
For a lift of order two, the two signs on the swapped pair multiply
to +1, but e_x e_(x+c)=±e_c changes sign when the pair is swapped.
The product of lifts of two distinct nonidentity elements therefore
fixes e_c. Under a hypothetical homomorphic section, that product
would be an involutive lift of the third nonidentity element, and
would have to negate e_c. Contradiction. Thus even the local
four-group action cannot be lifted coherently without signs.

## Nine path orbits and the one distinguished sign

The Paley affine group acts regularly on oriented arcs (unique
square slope and translation take any arc to 0→1). Its action
on Hamilton paths is free; an ordered path names all vertices.
Consequently the nine paths starting 0→1 are the exact orbit chart
for the 189 paths. The table below uses left-nested multiplication of
all seven imaginary units; full products are scalar because the seven
indices XOR to zero.

| Normalized path | Left product | Reverse-negate, then renormalize |
|---|---:|---|
| 0123456 | −1 | 0123456 |
| 0124563 | +1 | 0135246 |
| 0124635 | +1 | 0134526 |
| 0126345 | +1 | 0126345 |
| 0134526 | +1 | 0124635 |
| 0135246 | +1 | 0124563 |
| 0152346 | +1 | 0152463 |
| 0152463 | +1 | 0152346 |
| 0152634 | +1 | 0152634 |

The table is exhaustive by testing the 5! tails after 0→1; the
replay independently checks all 7! possible orders as well.
A full path's left-product sign is permutation parity times
the reference product of 0123456 (which is −1):
right multiplication by distinct imaginary units anticommutes
by alternativity. An octonion sign-gauge change multiplies
all full products by the same product of seven gauge signs, so
it preserves their relative split. The one negative orbit in
this convention is precisely the arithmetic progressions of
square step; its complement has eight orbits. The reversal-negation
involution fixes three orbit representatives and exchanges three
pairs. It preserves relative product sign (both reverse and
negate have odd permutation sign separately).

Thus three different normalized counts are genuinely encoded:
the unsigned total **9**, the positive part **8**, and the signed
difference **7**. Undoing the normalization multiplies by 21:
there are 168 positive, 21 negative, and signed difference 147.
This is an exact occurrence of 7/8/9 around the same object as
the modular group of order 168. It uses orbit counts or signed
counts, not the ordinary Hamilton count; it cannot infer either
hole in Chan–Pak. In fact no larger group of vertex permutations
can preserve the 168 positive paths in this fixed tournament:
their set of first arcs is the full 21-arc set (already the
orbit of any one positive path supplies it). So preserving
positive paths forces preservation of all oriented arcs, and
the group is exactly Paley's order-21 group, **not** the modular
group of order 168. The eight sign characters act trivially on
unsigned paths and multiply their signs by +1 (each nontrivial
character flips exactly four imaginary axes). This rules out
the obvious ways of turning the 168 equality into a compatible
regular action on paths in the fixed orientation. The Paley and modular connection at order 168
already appears in [the point-stripping ladder](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md).

## Catalan boundary: prime raw gaps in the 7 mod 24 slice

A Collatz odd valuation word w=(k_1,…,k_L) has k_i≥1 and
K=Σk_i. Its affine denominator before reduction at an integer
fixed point is Δ=2^K−3^L. Recover the precise carry distinction
from [the Catalan signed-cycle audit](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/catalan_elliptic_20260921_catalan.md).

**Theorem.** If |Δ| is a prime ≡7 mod 24, then (K,L)=(4,2) and
Δ=+7. Neither absolute gap 9 nor 21 is possible. Moreover no
primitive cycle at signed additive parameter ±1 belongs to one
of these three raw gap sizes.

For the prime claim K=1,2 gives only gaps 1 or 5. For K≥3
the negative case Δ≡−7 mod 8 is impossible because 3^L is
1 or 3 mod 8. Thus Δ≡7 mod 8 and ≡1 mod 3, forcing L and K even
respectively. Set K=2k,L=2s. Then the gap is
`(2^k−3^s)(2^k+3^s)`. Primality forces 2^k−3^s=1.
The only solution is (k,s)=(2,1): k is even by reduction modulo 3, and the two
coprime positive odd factors (2^(k/2)−1),(2^(k/2)+1) multiply
to 3^s, forcing the first to be one. The gaps 9,21 are impossible
because Δ≡2^K mod 3 is nonzero.

At (K,L)=(4,2) the ordered words (1,3),(2,2),(3,1) have carries
3+2^k1 = 5,7,11. Only the middle has carry divisible by seven.
Its candidate is n=b for additive parameter b=±1 on signed
integers, and its native odd step has valuation two and fixes n.
So that length-two word is repetition, not a primitive cycle.
The argument involves a special raw-gap slice, not a global
Collatz exclusion for sign ±1.

The same factorization applies to any positive Δ≡7 mod24
(K≥3). If it is a Carmichael number, write its displayed factor
pair u<v. If v were prime, Korselt would force v−1|uv−1, hence
v−1|u−1, impossible. If the smaller u is prime, then
u−1|v−1 and so u−1|v−u=2·3^s. This is a usable necessary
factor test on Carmichael candidates in that slice. It does not
turn λ(D)|D−1 into a Collatz orbit guard. See the archive's
[Carmichael successor blocks](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/carmichael_successor_blocks_20261014b.md),
which keeps the prime-factor certificate and the predecessor
lcm (e.g. 561 has slack (561−1)/λ(561)=7). This is unrelated
to claiming Carmichael's *totient* singleton conjecture from
the existence of a duplicate arithmetic representation.

## Bott is a stable periodicity, not an arrow-period rule

The relevant archives are the [Fano lattice carrier](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/creation_fano_20260925.md)
and [marked periodicity](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_bott_marked_periodicity_20261007d.md).
Seven signed octonion operators give the negative-square
Clifford representation; adding an eighth on the doubled
space is a different construction. Classical eight-step
Clifford stabilization multiplies the matrix size by **sixteen**,
while Bott's homotopy theorem is about stable groups.
See [Bott's primary paper](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/bott4.pdf)
and [Atiyah–Bott–Shapiro](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/abs.pdf).
Keeping an explicit marked matrix corner can preserve all affine
carry data. Keeping only a homotopy class cannot: the archive
exhibits both expanding and contracting words in the same
homotopy class. Likewise a signed lift of the Fano action cannot
be replaced by its unsigned projective data as a group action.
No numerical period or pointwise convergence follows from
matching eights.

Replay: `python check_signs.py` (signed multiplications, all paths,
all monomial lifts, local four-group, finite raw-gap controls).
