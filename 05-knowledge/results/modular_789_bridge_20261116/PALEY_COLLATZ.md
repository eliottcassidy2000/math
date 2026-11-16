# The signed mod-seven Collatz branch action

Fix b=+1 or −1 and use accelerated odd letters on their native valuation
charts. A letter of exponent k≥1 induces on F_7 the affine map
`f_k(x)=2^(−k)(3x+b)`. This is a **marked** map: the exponent must
come from the original integer. But no residue x mod7 excludes a
single letter: the equation 3n+b≡2^k mod2^(k+1) selects odd n with
exact valuation k, and CRT gives infinitely many positive such n
in each residue modulo seven. This local realization does not make
the unmarked reduction into a deterministic map. More strongly,
every finite exponent word occurs on a full 2-adic source cylinder:
working backwards, if the required next residue is odd r modulo
2^s, impose 3n+b≡2^k r modulo 2^(s+k). This uniquely selects
an odd preceding residue, makes its k native, and gives the
required quotient residue. Start with odd terminal residue mod2
and iterate. CRT with modulo seven now lets any such finite chart
start at any of the seven marked residues. This does not provide
one source for every requested **infinite** word.

In the archived Paley convention arrows have square difference in
{1,2,4}. Since 2 is square and 3 is not, every f_k reverses every
arrow, for both signs of b. It also exchanges the two Fano line
planes dev{0,1,3} and dev{0,1,5}: all square-affine maps preserve
each, and negation exchanges them. Thus a word with L accelerated odd letters
preserves the marked plane and arrows iff L is even. Here the Fano
line table before and after the switch must be distinguished; this
is not a single octonion automorphism on a fixed signed table.

As permutations the letters generate precisely AGL(1,7), of
order 42, including its arrow-preserving half of order 21.
Indeed f_k f_(k+1)^(-1) is doubling, so the group also contains
u(x)=3x+b. Doubling conjugates u to 3x+2b; composing that
with u^(-1) gives translation by b, which generates all translations.
Finally 3 has order six modulo seven, hence all slopes occur.
These are permutation-group operations on the seven labels, not
claims about simultaneous inverse itineraries on the positive integers.

For a word of L letters of total exponent K, the slope is
3^L 2^(-K). It is one exactly when L is even and K≡L/2 mod3,
using 3²≡2 and the order three of 2. Its mod-seven fixed-point
gate is then complete: if the slope is not one, the marked word
has exactly one algebraic fixed residue; if the slope is one,
there are no fixed residues unless its translation vanishes, in
which case all seven residues are fixed. For two letters and
K≡1 mod3, the translation vanishes precisely when both exponents
are 2 mod3: before the denominator the relevant carry is
b(3+2^k1). This directly interprets the little K=4, L=2
raw-gap test in SIGNS.md as the only zero-translation chart
among (1,3), (2,2), (3,1). Again the mod-seven fixed residue
is not a guarantee of a native integer cycle.

Prior material: [the marked Paley bridge](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_paley_bridge_20261001.md)
distinguishes the trivial-cycle Frobenius from global Collatz symmetry;
[the audited two-plane Paley result](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_mod6_20260922_paley_fano_octonion_design.md)
proves the Fano plane decomposition and its exchange under negation.
The interface retained here is the valuation-marked affine letter,
not the unlabelled integer dynamics.
