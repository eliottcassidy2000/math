# A four-coordinate semilinear census beyond p-group symmetry

2026-10-20. **PROVED family theorem; FINITE-EXACT full-automorphism census at
q=25.** Independently reconstructed and audited from the staged computations.
No historical priority claim is made.

## 1. Results and exact scope

Let p be an odd prime with p≡3 or 5 mod8, let s≥1 be odd, and put

    q=p^(2s),   r=(q−1)/8,   H=F_q⁺ ⋊ C_r,

where C_r is the subgroup of order r in F_q^*, acting by multiplication.
Then r and |H|=qr are odd. The number τ_T(H) of ordinary tournament isomorphism
types whose automorphism group contains a conjugate of this permutation group
H equals the number τ_E(H) of such **sign-even** graph types, and

    τ_T(H)=τ_E(H)=1.

The sole graph type is the empty graph. For q>9 these are genuinely non-p-group
symmetries, since r>1 and gcd(p,r)=1. Examples include group orders 75 at q=25,
1815 at q=121, 3549 at q=169, and 66339 at q=729. The q=9 boundary has r=1 and
is a p-group control, not a new mixed-prime example.

Here sign-even means every graph automorphism reverses an even number of edges
relative to a reference orientation. This character does not depend on the
reference orientation. It is not even degree. In fact every graph in the
fixed cube below has even degree.

For q=25, the full normalizer and two independent complete rooted-isomorphism
searches give a sharper exact census: the 16 invariant tournaments form one
global type with automorphism group H of order75; the 16 invariant graphs
form six global types, only the empty one sign-even.

This is a family of embedded-subgroup symmetry classes, not a census of all
odd non-p subgroups in S_q. There is no assertion for arbitrary extensions,
for cyclic rings Z/p^a, or for all tournaments and graphs of order q.

## 2. Inheritance and the changed coordinate

The [prior session map](../prime_power_bijections_session_20261020/README.md)
proposed F25⁺⋊C3 as a small test beyond the ambient-Sylow base route. Its
translation subgroup has order25, while a Sylow5 subgroup of S25 has order
5^6; the ambient-Sylow theorem therefore does not directly apply to this base.

Closest mechanism: the odd fixed-cube signed Fourier representation in
[PRONORMAL_QUOTIENT_THEOREM.md](../pronormal_quotient_bijections_20261020/PRONORMAL_QUOTIENT_THEOREM.md).
Least-used sidecar: the distinction between multiplier-relative and fully
sign-even graphs in
[GENERATOR_FREE_NATURALITY.md](../generator_free_circulants_20261020/GENERATOR_FREE_NATURALITY.md).
Hostile: the adjacent two-orbit graph at q=25, which is multiplier-even but has
an odd-sign Frobenius automorphism.

The operative concept board is: the odd fixed cube; signed multiplier action;
semilinear source fusion; automorphism-sign escape; and full versus relative
isomorphism. The family theorem uses a concrete semilinear subgroup, not an
unproved pronormality or full-normalizer assumption.

## 3. The four-coordinate signed cube

Since s is odd,

    v₂(q−1)=v₂(p²−1)=3.

Indeed (p²)^s−1=(p²−1)(1+p²+⋯+p^(2s−2)), with an odd second factor, and
p≡3,5 mod8 gives v₂(p²−1)=3. Thus r is odd.

Choose a primitive element g of F_q^*. Put U=⟨g^8⟩, of order r, and
C_j=g^jU for j modulo8. The eight directed difference classes are paired by
negation: −C_j=C_(j+4), since −1=g^(4r) and r is odd.

Translation invariance makes every invariant object a Cayley object. U
invariance makes every connection set a union of the C_j. Therefore:

- An H-fixed tournament chooses exactly one of C_j,C_(j+4) for each
  j=0,1,2,3. There are exactly16 choices.
- An H-fixed graph chooses any of the four unions C_j∪C_(j+4). There are
  exactly16 choices. A mask b means its selected indices are those j with
  the j-th bit of b equal to1.
- Each unordered-pair orbit has size qr: every vertex has2r differences in
  one union, and division by2 removes the two ordered presentations. This
  orbit size is odd.

Multiplication by g^k and field Frobenius induce on the directed classes

    j ↦ k+p^e j mod8,    e=0,1.

Only p modulo8 matters. The resulting action on the directed classes is
C8⋊C2, of order16. It is an induced action: Frobenius on the full field may
have order2s rather than2, and no claim of a faithful order16 subgroup of the
full permutation group is needed.

Orient every pair orbit using differences in C_j, for j=0,1,2,3. If the
semilinear map stabilizes a graph with selected coordinate set S⊆{0,1,2,3},
its edge-reversal sign is

    (−1)^#{j∈S : (k+p^e j mod8)∈{4,5,6,7}}.             (1)

To justify (1), each selected pair orbit is either mapped with its reference
orientation preserved or reversed. A reversal contributes its entire odd
size qr. Orbit permutation itself contributes no sign to this edge-reversal
character. This last point is essential: the character is not the sign of
the permutation of the unordered edges.

## 4. Proof of the one-type family theorem

**Tournament side.** A skew selector in Z/8 cannot be fixed by any nontrivial
translation. Otherwise its stabilizer, a nontrivial subgroup of C8, would
contain translation by4, which takes every skew selector to its complement.
Hence the16 source objects form two scalar orbits of size8.

Take the selector I={0,1,2,3}. Frobenius maps it to

- {0,1,3,6} when p≡3 mod8;
- {0,2,5,7} when p≡5 mod8.

Neither set is a translate of four consecutive residues. Thus Frobenius
interchanges the two scalar orbits, and all16 tournaments form one semilinear
orbit. They consequently form one ordinary isomorphism type as well.

**Graph side.** Up to cyclic rotation, the six graph masks are
0,1,3,5,7,15. The following stabilizers have negative sign by (1):

| Nonzero mask | Selected coordinates | k | e | Number of reversed coordinate orbits |
|---|---|---:|---:|---:|
| 1 | {0} | 4 | 0 | 1 |
| 3, p≡3 mod8 | {0,1} | 1 | 1 | 1 |
| 3, p≡5 mod8 | {0,1} | 0 | 1 | 1 |
| 5 | {0,2} | 2 | 0 | 1 |
| 7 | {0,1,2} | 4 | 0 | 3 |
| 15 | {0,1,2,3} | 1 | 0 | 1 |

For example, when p≡3 mod8 the map j↦1+3j sends {0,1} to {1,4}; modulo4
it stabilizes the graph, but exactly one selected orientation reverses.
Every nonempty invariant graph is therefore globally sign-odd. The empty
graph is globally sign-even regardless of its automorphism group. This proves
τ_T(H)=τ_E(H)=1. No computation at unbounded q or classification of full
Cayley-graph automorphism groups is used. ∎

This resolves the named q=25 frontier and simultaneously supplies an infinite
family: fixing p=3 and allowing odd s already gives unbounded q. It does not
repair the unrestricted odd-group conjecture, which remains false by the
inherited C15 example on10 vertices.

## 5. Full q=25 census and exact normalizer

Use F25=F5[t]/(t²−2), labelled a+5b for a+bt. The polynomial t²−2 is
irreducible over F5, and g=2+t, label7, has order24. Thus U=⟨g^8⟩ has
order3 and H has order75.

**Normalizer proof.** The translation group P=F25⁺ is the unique Sylow5
subgroup of H, since its number divides3 and is1 modulo5. It is therefore
characteristic in H, so any normalizer of H lies in N_(S25)(P)=AGL(2,5).
The multiplier of order3 acts irreducibly over F5, since3 does not divide4.
Its centralizer in GL(2,5) is F25^*, and conjugation on its C3 subgroup is
either identity or inversion. Frobenius realizes inversion. Consequently

    N_(S25)(H)=AΓL(1,25),       |N|=1200.

The ordinary affine-multiplier subgroup L=AGL(1,25) has order600.

The exhaustive computations yield:

| Graph mask | Masks in its N/global orbit | Degree | Full automorphism order | Globally sign-even? |
|---|---|---:|---:|---|
| 0 | 0 | 0 | 25! | yes |
| 1 | 1,2,4,8 | 6 | 300 | no |
| 3 | 3,6,9,12 | 12 | 300 | no |
| 5 | 5,10 | 12 | 600 | no |
| 7 | 7,11,13,14 | 18 | 300 | no |
| 15 | 15 | 24 | 25! | no |

There is one tournament N/global orbit, containing every mask0 through15.
Its full automorphism order is75, exactly H. The two L-orbits have masks

    {0,1,3,7,8,12,14,15}  and  {2,4,5,6,9,10,11,13}.

For every nonempty, noncomplete graph representative and the tournament
representative, both complete rooted searches return exactly the rooted
normalizer automorphisms. Translations show that multiplying the rooted count
by25 gives the full automorphism order. For empty and complete graphs the
full group is S25 directly; the empty sign is trivial, and a transposition
makes the complete-graph sign nontrivial.

Degrees separate all graph orbits except masks3 and5. Those two are separated
by the common-neighbor histogram relative to vertex0:

    mask3: 4,5,6,7 common neighbors, each occurring6 times;
    mask5: 5 or6 common neighbors, each occurring12 times.

The two independent rooted-isomorphism searches also verify no isomorphism
between these graphs. No all-permutation search over25! permutations is
claimed; exhaustive refinement and constraint propagation certify that no
omitted branch can contain an isomorphism.

## 6. A precise Frobenius fusion/sign cancellation

If one stops at L=AGL(1,25), one sees two tournament types and two L-even
graph types, represented by masks0 and3. Enlarging to the full normalizer
has two effects:

1. Frobenius fuses the two tournament types into one.
2. Frobenius makes the mask3 graph sign-odd, leaving only mask0.

For mask3, σ(a+bt)=a−bt fixes the graph and swaps the endpoints of exactly
five of its edges. The involution's other edge orbits contribute positive
sign in pairs, so σ has sign−1. This is also checked directly against the
edge-reversal formula on actual labelled edges.

Thus relative-to-global source fusion excess and even-graph sign escape are
both1 for L, but both0 for N. The global equality is real; the intermediate
two-versus-two equality is not already the global census. This example shows
why a prime-power multiplier necklace model must not be silently substituted
for the fully sign-even target.

Connection contract:

- Source: the four-coordinate H-fixed tournament cube.
- Map: quotient first by scalar multipliers, then by the semilinear action.
- Preserved predicates: actual isomorphism under every acting permutation,
  and the explicit edge-reversal sign character.
- Information lost by the scalar-only quotient: Frobenius fusion and odd-sign
  witnesses outside the scalar stabilizer.
- Needed sidecar: the signed action on directed difference classes modulo8.
- Decisive test: the interval skew selector and the adjacent two-coordinate
  graph, which undergo the two compensating changes above.

This is a concrete descent mechanism, not a universal theorem that source
fusion and sign escape always cancel for arbitrary odd subgroups.

## 7. Reproduction, audit, and next targets

Run from this directory:

    python3 verify_affine25.py --output verification.json > verification.out
    python3 -O verify_affine25.py --output verification.optimized.json > verification.optimized.out
    cmp verification.json verification.optimized.json

Both outputs are byte-identical and contain all rooted automorphisms,
normalizer orbits, degree histograms, graph signs, exact signed-cube witness
tables for p mod8 equal to3 and5, and288 comparisons of the compressed
sign formula with full F25 edge signs. All validations use explicit exceptions,
so optimized Python does not disable them. See VALIDATION.md.

Precise next targets:

1. Move from v₂(q−1)=3 to v₂(q−1)=4. The fixed cube grows from4 to8
   coordinates, and the induced Frobenius residue must be retained; no
   one-type formula is predicted without checking it.
2. Test an odd multiplier subgroup smaller than the maximal odd subgroup.
   This preserves the Cayley description but introduces more difference
   coordinates and potentially new full-isomorphism fusion.
3. Seek object-relative fusion/sign descent conditions weaker than ambient
   pronormality. The present family supplies a positive benchmark but proves
   neither pronormality nor a general necessary-and-sufficient criterion.

None of these results settles LRC(14), a Collatz statement, or the unrestricted
odd-group equinumerosity question.
