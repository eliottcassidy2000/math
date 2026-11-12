# Contraction decks: exact module-switch counts and a residue/interval separation

**PROVED deductions from the inherited substitution kernel; finite controls
below. No Chan–Pak completeness claim.** This is the constructive continuation
of the one-count-code lane. The new organizing object is the *contraction deck*
of a physical module-switch family, not another universal classification claim.

## 1. Independent local cyclicizations have an exact positive expansion

Let B be any tournament containing k pairwise disjoint three-vertex modules
M_1,...,M_k, each inducing TT_3. Every exterior vertex sees each module
uniformly. For E⊆[k], let T_E replace the modules indexed by E by directed
triangles C3, leaving all other arcs fixed. This is a same-size physical
construction; it can be performed by reversing the long arc in each TT3.

For S⊆[k], let B/S contract each M_i, i∈S, to one vertex, retaining the
uniform external arcs and leaving other modules unchanged. Put a_S=H(B/S),
so every a_S is a positive odd integer, and a_empty=H(B).

**Exact module-switch law.**

    H(T_E) = sum_(S⊆E) 2^|S| a_S.                       (A1)

**Proof.** Ordered path-cover profiles of TT3 and C3 are respectively
(1,6,6) and (3,6,6). Their difference is2 times the singleton profile,
zero-extended. The inherited uniform substitution contraction is multilinear
in each factor profile. Expanding independently at every changed module,
choosing a singleton profile at sites S is precisely the Hamilton count of
B/S, multiplied by2^|S|. This proves (A1) for arbitrary exterior quotient
and modules, not just a transitive exterior. QED.

The path-cover version also holds coordinatewise, with a_S replaced by the
ordered path-cover profile of B/S, with zero extension as appropriate. This
last statement concerns a linear profile identity across differing orders;
it does not identify the physical graphs or their isomorphism classes.

## 2. Fixed-precision truncation is an exact recursive correction rule

For every r>=1,

    H(T_E) ≡ sum_(S⊆E, |S|<r) 2^|S| a_S  (mod 2^r).     (A2)

The degree-j contraction deck is needed only modulo2^(r-j). Thus fixed r
uses at most sum_(j=0)^(r-1) binom(k,j) deck entries, O(k^(r-1)) for r>=2.
The omitted terms vanish *exactly*, not because they are small in absolute
value. Computing each deck entry can still be expensive. This is a reduction
in the number of response inputs for fixed low-bit precision, not a general
polynomial-time Hamiltonian counter or ordinary-class enumerator.

It is an actual correction hierarchy compatible with the family grammar:
coarse data a_empty are supplemented by singleton contractions, pair
contractions, and so on. All degrees restore the exact integer. The single
physical code in PROOF.md can encode the profile at any node, but gives no
free shortcut for obtaining these distinct contraction nodes.

## 3. Complete mod8 criterion

Assume k>=1. Write `a_{i} ≡ 1+2b_i (mod4)` with b_i∈{0,1}. For E of size s,
(A2) and oddness of every pair contraction give

    H(T_E)-H(B) ≡ 2s² + 4 sum_(i∈E)b_i  (mod8).          (A3)

Therefore:

- if all b_i=0, the output residues are exactly H(B)+{0,2};
- if all b_i=1, the output residues are exactly H(B)+{0,6};
- if both types occur, all four odd residues modulo8 occur.

For the last statement choose E empty, a type0 singleton, a type1 singleton,
and that mixed pair. Their offsets are0,2,6,4. The first two statements
follow because s² is0 modulo4 for even s and1 for odd s, and because empty
and singleton choices occur. For k=0 there is only the base residue.

The contraction coordinate cannot be dropped. A tempting but FALSE shortcut
would replace every a_i by H(B) modulo4 and conclude universal
`H(T_E) ≡ 3^|E| H(B) mod8`. The smallest stated hostile is
`B=C3[TT3,1,1]`, where H(B)=9 and the changed count is15, whereas3*9=27.
The contracted count is3, not9 modulo4. This failed guess was rejected before
publication; (A3) retains the missing coordinate.

## 4. An explicit mixed-type physical cube

Take

    B = TT3 -> C3[TT3,TT3,TT3].

Choose M_1 as the source TT3 and M_2 as one module inside the cyclic core.
The four literal twelve-vertex outputs are

    E={}      : 2721,
    E={2}     : 2823,
    E={1}     : 8163,
    E={1,2}   : 8469.

Here a_1=2721, a_2=51, a_12=51. The single contractions differ modulo4,
so the four residues are1,7,3,5 modulo8. The values are direct consequences
of (A1), independently checked by vertex-subset path counting. This is a
small explicit residue controller, not a new global residue-surjectivity
result: Session281 already supplies a much stronger all-modulus theorem.

## 5. Why independent switches do not make an unbroken interval

Every a_S>=1, so

    H(T_[k])-H(B) >= 3^k-1.                              (A4)

There are at most2^k distinct outputs. If k>=2, then

    (H(T_[k])-H(B))/2 +1 >= (3^k+1)/2 > 2^k.

Thus **every independent TT3/C3 switch cube with at least two modules misses
some odd integer between its own minimum and maximum**, at every core size
and with arbitrary exterior structure. All coefficients in (A1) are positive,
so the empty and full subsets really are its extrema. Relabelling or merging
isomorphism classes can only reduce the number of output values further.
The theorem does not bound the length of subintervals inside the cube, or the
union of cubes with different bases, or compensated non-module edits.

This yields the precise constructive bridge obligation: the family grammar
can compose and classify these switches and can realize all low residues in
the mixed case, but a full Chan–Pak selector must either (i) supply different
absolute bases whose intervals overlap, (ii) compensate the positive higher
interactions, or (iii) leave this independent uniform-module switch language.
The receiver-preserving compensated gates of Sessions280–282 are candidates
for (ii)/(iii); (A1) alone does not supply their missing renewal/anchoring step.

## Reproduction and scope

`python3 switches.py` and `python3 -O switches.py` check the subset expansion,
full low-bit truncations, mod8 classification, and literal output counts on
small complete switch universes. These checks do not turn marked assignments
into an A000568 isomorphism census. Generic multilinearity is inherited;
(A1)–(A4) are the explicit special-family deductions and boundaries here.
