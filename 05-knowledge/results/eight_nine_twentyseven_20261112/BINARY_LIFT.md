# A physical mod-eight selector and the role of multiplication by nine

2026-11-12. **PROVED:** every binary precision is reachable in the stated
physical family, with an exact digit-lifting algorithm and common order.
This is congruence coverage, not equality to every prescribed odd integer.
It depends on the literal selector in [constructors/PROOF.md](constructors/PROOF.md).
Scalar attainability at every modulus was already proved in
[Session137, width_phase2](../session137_overnight_proof_frontiers_20261102/width_phase2/PROOF.md)
and inherited by Session139. The contribution here is this explicit
fixed-core, retained-exterior recipe, its unique addresses, and its
common-order realization; it is not a new general modular-existence claim.

## 1. The four configurations retain their high carry

Fix any uniform exterior E around the eight-vertex core C. The two disjoint
arc switches in the constructor proof give four tournaments of the same
order. Let

    h=H(E[C]),   A=H(E[TT4]),   B=H(E[TT2]).

Their counts, with switches a,b in {0,1}, are exactly

    N_ab=h+2a*A+4b*B-8ab*(A+2B).                      (1)

H of any tournament is odd, so A and B are odd. Consequently the four
values give each odd residue modulo8 exactly once: b toggles the residue
by4, while a moves to the other class modulo4. The mixed term is a genuine
negative response of size8 times a positive context count. It cannot be
discarded when requesting an exact integer target.

For empty E the four values are513,515,517,495, with residues1,3,5,7.
All four configurations keep the same eight vertices. No quotient by
isomorphism or replacement of a missing arc by a formal weight is used.

## 2. Nine moves through every higher binary digit

For s>=3, the powers of9 form precisely the subgroup

    {u mod2^s : u=1 mod8},

and its order is2^(s-3). To see this directly, repeated squaring gives

    v_2(9^(2^j)-1)=j+3       for every j>=0.           (2)

The base case is9-1=8. If v_2(y-1)>=3, then y+1 has valuation1,
so squaring raises that valuation by exactly one. For arbitrary positive
exponent t, factoring out its odd part gives v_2(9^t-1)=3+v_2(t).
This proves the stated order; every power is1 mod8, and there are exactly
2^(s-3) such residues, so the subgroup is complete.

It follows that the indexed family

    N_ab * 9^t,
    a,b in {0,1},  0<=t<2^(s-3),                     (3)

contains each odd residue modulo2^s exactly once. First choose the unique
N_ab with the desired residue modulo8. Its quotient from the desired
unit modulo2^s belongs to1+8Z, which has a unique power-of9 address.
The division here is modular inversion of an odd number, not division
of a graph or an assertion that the target integer factors.

These are literal tournament counts. Let D be the ordinal join of two
directed triangles, so |D|=6 and H(D)=9. Append t copies of D by ordinal
join. If all outputs must have the same order, add neutral transitive
vertices until each has order

    |E[C]|+6*(2^(s-3)-1).                            (4)

Thus every odd binary residue is attained at this common order while the
original exterior and one of its four core configurations remain induced.
The bound is explicit, not claimed efficient or minimal.

### A terminating digit algorithm

After choosing N_ab, let q be the required ratio, q=1 mod8. Start t=0.
For j=0,...,s-4 compare9^t and q modulo2^(j+4). If they disagree, add2^j
to t. Equation(2) says that multiplication by9^(2^j) toggles exactly the
next unresolved bit while preserving all earlier ones. At the end,
9^t=q mod2^s and0<=t<2^(s-3). This uses s-3 binary decisions, even though
the resulting literal graph order can be much larger.

The precise role of8 and9 here is therefore a residue selector followed
by a generator of the subgroup above its fixed low bits. The construction
does not turn the mod-eight high-carry term in(1) into an arbitrary digit.

## 3. A second exact 8/9 split appears in the primitive-return guard

The constructor's infinite rooted return reduces coprimality to excluding
`2^k=28 mod73`. Here

    73=8^2+8+1=Phi_9(2)=8*9+1.

The element2 has order9 modulo73: its first nine powers are distinct and
2^9=1. The power map chi(x)=x^9 on the nonzero residues has this nine-element
subgroup as its kernel and eight distinct phase values. Explicitly,

    28^9=10 mod73,   10^4=-1 mod73,   10^8=1 mod73.

Thus28 has a nontrivial phase whereas every2^k has phase1. This gives the
same all-k exclusion as the literal nine-state table. In fact28 has order72:
28^8=16 has order9, and28^9=10 has order8.

This is an actual order-nine orbit with eight cosets in a finite field.
It is distinct from the four odd residue classes modulo8 in section1;
we do not identify those two sets merely because the number8 appears.
Both supply exact finite address maps used by physical constructors.

## 4. Why this universal selector is specifically binary

No operation restricted to the substituted core can guarantee all residues
modulo an odd prime p for **every** uniform exterior. Choose as exterior
any cyclic tournament on p vertices, all pointing into the core. Its
translation group of order p acts freely on Hamilton paths, so its count
is divisible by p. The ordinal join therefore has count divisible by p
for every replacement core, regardless of the number of configurations.
For p=3 the exterior can simply be a directed triangle.

Thus every modulus with an odd prime factor has such an exterior obstruction
to universal residue coverage. Powers of2 are different here for a proved
reason: every tournament Hamilton count is odd, so both response weights
in(1) remain units in the relevant binary calculation. This is a boundary
of the universal-in-exterior claim, not a restriction on all constructors
that are allowed to change the exterior as well.

## Scope and reproduction

`binary_lift.py` implements the residue selector and bit algorithm, then
constructs the ordinal joins and neutral padding. The controls exhaust
all odd residues through precision2^10 and directly count the eight
common-order outputs at precision2^4. Larger conclusions follow from
(1)–(4), not from the finite test. A printed graph realizes its emitted
count modulo the requested modulus; exact-target coverage remains open.
