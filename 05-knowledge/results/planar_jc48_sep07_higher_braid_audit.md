# Independent audit of the three higher-cusp affine braid certificates

**Status: INDEPENDENT ANALYTIC + EXACT SOURCE + REPLAY AUDIT PASS.**
I read the complete [primary proof](planar_jc48_sep07_higher_braid.md)
and standalone source, checked the
actual geometry, rational transport and arbitrary-group eliminations,
and independently checked all frozen raw/compressed witness hashes.
Both complete normal and optimized replays pass and agree byte for byte
with the frozen output. No numerical scout word is used as authority.

## 1. These are the actual curves, including the co-projection

All three literal pairs have `U=t^4+t^3+t^2`. Their second coordinates are

```text
finite7 / infinity9: 2t^6+3t^5-(3/2)t^3-(3/2)t^2,
finite7 / infinity7: 2t^6,
finite9 / infinity7: 2t^6+(12/7)t^5-(6/7)t^3-(6/7)t^2.
```

The source constructs each ordinary resultant `Res_t(U-u,V-v)` and
compares it to the complete literal monic quartic, then substitutes the
parametrization back into it. Its separate symmetric pair calculation
pays birationality: the first divided equation is
`p(p^2+p+1)-(2p+1)q=0`, and at `p=-1/2` its numerator is `-3/8`, so
no pair is lost by the displayed division. The remaining off-diagonal
pair polynomials are finite, with nonzero discriminants and nonzero
resultants against the pair-discriminant factors. Thus the normalization
is generically injective. Monicity of `U` makes its map to its image
finite, and `[C(t):C(U)]=4`; the monic resultant is therefore the actual
minimal quartic equation, with no extraneous component.

The derivative gcd `t` and the gcd `t^2` of `U,V` retain the unique
finite critical normalization point and exclude a second preimage of
the cusp. The finite characteristic calculation removes the actual even
coefficients and finds the first odd terms `9/2`, `-6`, `-44/7` at the
respective orders seven, seven and nine. The infinity chart and its odd
coefficients give the stated infinity types and line contact six.

The tangent ideal is checked before any inference that the pairs are
nodes: its displayed Groebner basis is supported only on the cusp
diagonal. The triple test also retains its necessary denominator sidecar.
The remainder of `V-B` modulo `U-A` is cubic with **constant nonzero**
leading coefficient. Three distinct common preimages would force that
cubic to divide `U-A`; its remainder coefficients generate the unit
ideal. Consequently there is no triple image and no denominator locus
excluded by solving this remainder problem. The node counts three,
four and three agree with the whole rational-sextic genus ledger.

For finite-seven/infinity-seven, the cusp and one node have the same
`u`-projection. The source explicitly keeps

```text
F(0,v)=v^2(v-2)^2,
U=0, V=2 on t^2+t+1=0,
```

so this node is a different target point from the cusp `(0,0)`. The
factor `u^9` in the vertical discriminant must not be reinterpreted as a
finite-nine cusp: it contains the co-projected node contribution. Neither
the local characteristic calculation nor the generic-fibre group argument
requires distinct projections of all singular points.

## 2. Exact rational root transport

The default verifier uses rational real and imaginary parts. Its numerical
root proposals are confined to the separately named production and scout
modes. For every frozen segment it retains all Taylor coefficients of
the literal degree-six-in-`u`, monic-degree-four-in-`v` equation.
The inequalities
`max(|Re z|,|Im z|)<=|z|<=|Re z|+|Im z|` make the displayed Taylor
bound a valid strict Rouche comparison with the linear root term.
It gives one simple actual root in each disk throughout the full base
segment, not merely at its endpoints. Each radius is at most one
thirty-second of every relevant lower separation, so the four disks are
pairwise disjoint and exhaust the quartic's roots.

The endpoint displacement is less than half the old radius. A separate
endpoint Rouche gate isolates the new root within `min(old,new)/64`.
Thus the same root lies in successive disks, and labels cannot be
reassigned by an uncertified nearest-neighbour choice. Initial isolation,
identical labelled base centres for the two loops of a case, and exact
final unordered closure make the two based configuration loops compatible.
Inside disjoint convex disks, the root strands and rational polygonal
strands are homotopic, with the same endpoint corrections.

The actual base path is checked as well: each point lies on the declared
oriented rational edge with a strictly increasing affine parameter, and
all six edges are completed. The orientation-preserving multiplier
`1+i/4` gives exact rational crossing times. The source checks distinct
crossing times, nonzero endpoint real separations, adjacency of the two
crossing strands and nonzero imaginary separation. Its sign is the
inherited below-stem convention, and only adjacent inverse letters are
cancelled. These operations extract the actual braid class from the
certified configuration loop.

## 3. Necessary word constraints in an arbitrary group

Use chronological Hurwitz action

```text
H_i^+(a,b)=(aba^-1,a),   H_i^-(a,b)=(b,b^-1ab).
```

For a conjugate word `A sigma_i A^-1`, fixedness of a tuple is equivalent
to equality of the two entries affected by `sigma_i` after the prefix
`A`. This is valid in every group because the prefix and suffix actions
are inverse automorphisms. It supplies the necessity direction, in
addition to the source's free-word substitution checks of the converse.

Write an original tuple as `(a,b,c,d)`. For finite seven/infinity nine,
the two prefixes are `[3,2,3,2]` and `[-3,-3]`, with central letters one
and two. I independently obtain the equalities

```text
a=b c d c d^-1 c^-1 b^-1,
b=d^-1 c d.
```

Thus all four entries are words in the original positive meridians
`c,d`. The explicit two-generator tuple in the producer is the same
solution, and substituting it fixes both complete words.

For finite seven/infinity seven, the prefixes are `[2,-1,-2,1]` and
`[-1,-2,-3,-1,-1]`, both with central letter two. Their equalities are

```text
b c b^-1=c^-1 b^-1 a b c,
d=c^-1 b^-1 c b c,
```

or equivalently

```text
a=b c b c b^-1 c^-1 b^-1,
d=c^-1 b^-1 c b c.
```

This time the two original positive meridians are `b,c`.

For finite nine/infinity seven, both certified words have the literal
common conjugator `R=[2,3]`. Applying its Hurwitz automorphism once
changes the tuple to

```text
(a, b c b^-1, b d b^-1, b).
```

Each new entry is a positive meridian, and they remain a free meridian
basis. This is a **basis change**; no assertion that `R` itself is an
actual loop in the `u`-plane is required. The actual loop relations are
still the two full certified words. In the new basis their shorter
prefixes are `[2,2,-1]` and `[-1,-2,-1,-3,-1]`. Independent multiplication
of their affected entries gives exactly the same two solved equations
as in the preceding case. Hence two positive meridians in this new basis
generate all four; the inverse basis change recovers the original four.
This argument does not infer group generation from a finite permutation
census or mistake an abstract basis braid for an actual projection loop.

## 4. Actual affine group and all-degree consumer

The monic quartic supplies a continuous global section above a Cauchy
root bound. Over any compact disk filling a certified base loop, that
section can be chosen with a constant sufficiently large vertical
coordinate. Its loop is null in the whole affine complement. Over the
regular base, fibre meridians and section lifts generate the bundle
group. Restoring the finitely many critical fibres kills those section
lifts, while general position lets loops in the full complement avoid
the critical fibres. Therefore the four fibre meridians surject onto the
actual affine fundamental group, and the certified braids fix those
elements exactly. Merely simultaneous conjugacy would not suffice;
the section is the sidecar that pays this marked fixedness.

The two generators obtained in Section 3 are positive curve meridians,
including after the finite-nine basis change. For a transitive action
on `d>1` labels generated by `r` permutations moving at most `delta`
labels each, replace every nontrivial cycle by a tree on its support.
The union graph has exactly the group orbits as its connected components.
It needs at least `d-1` edges and each generator supplies at most
`delta-1`. Nontrivial transitivity forces `delta>=2`; identity generators
then contribute zero without violating the upper bound. Thus
`r(delta-1)>=d-1` in every degree.

For an irreducible whole nonproperness curve, each of the two positive
meridians fixes its `a` actual retained sheets, so `delta<=d-a`.
The inherited odd-cusp node ledger with at least two nodes gives
`d<=2a`. The two-generator inequality is then impossible. This avoids
any unsupported all-degree assertion about groups generated by single
long cycles. The combinatorial bound is sharp for two cycles meeting
in one label, but those abstract actions are not Keller maps.

I also checked the final correction to the node-ledger definition.
For retained branch sets `A,B`, the ledger variable is the overlap
`omega=|(Omega\\A) intersect (Omega\\B)|` of their **deleted complements**.
The actual node fibre has size `|A intersect B|=2a-d+omega`.
Since the affine-line-normalized nodal curve has Euler characteristic
`1-N`, its complement has Euler characteristic `N`, and its smooth
stratum after removing the cusp and nodes has characteristic `-2N`.
Euler integration therefore gives
`1=dN-2aN+sum(2a-d+omega)+n=n+sum(omega)`.
This explicitly pays the definition in the repaired primary text;
calling `omega` retained overlap would be incorrect. The formula,
source and path certificates did not change in that prose repair.

## 5. Exact universe and replay manifest

The finite universe is exactly the three literal coefficient rows,
their full pair/triple/tangent geometry and all six frozen loops.
I independently checked that every witness row contains exactly four
Gaussian rational centres, and that both raw and compressed byte pins
agree with the primary manifest. No larger curve or permutation census
is inferred. The analytic arguments above supply the actual-affine-group
map and the necessary group constraints, beyond merely a passing replay.

I independently replayed the full default normal and optimized verifier.
Both exit successfully and produce exactly the same 2,920 bytes as the
frozen output. Reproduction from the repository root:

```bash
python3 04-computation/planar_jc48_sep07_higher_braid.py
python3 -O 04-computation/planar_jc48_sep07_higher_braid.py
```

All three outputs pass 286,548 always-active gates over
31,807 segments, in the following exact partition:

| Case | Upper/lower segments | Gates |
|---|---:|---:|
| finite7 / infinity9 | 2930 / 6979 | 89,262 |
| finite7 / infinity7 | 2531 / 3745 | 56,578 |
| finite9 / infinity7 | 4760 / 10862 | 140,708 |

SHA-256 pins independently read from the files:

| Artifact | SHA-256 |
|---|---|
| Shared source | `e6b39ca3efa88fe25884331ad8c062f2ec8014d923a47abf70bcbffa315eb5d0` |
| Shared output | `08a41085782dab951ede0882005906d16602083a304a8bfc28de8bf5d27ac088` |
| 79 gzip | `1393ab1970614903a8f7b8b96c550d2a6c67337b9d7afa3221a588167f4035ae` |
| 79 raw JSON | `6a60f48122f2fcdeee3a5eb5b14991b4548cdddda71d176993e0143d93df5379` |
| 77 gzip | `c4ab6a7b28ed6fa34bd8d5071525a805c963e12f365ab9c6018cbca2cf7a4fd9` |
| 77 raw JSON | `405232fd9727a696455a6133f27b8650aed9fa4e4bd75858880ffde893e74ba7` |
| 97 gzip | `1ee79c177a0d6d4783a88ce0613c3bba4067db9caa44a82783c4b683fd3c7932` |
| 97 raw JSON | `0b2ad7bcbc65fd9d2ffbfb87a13a573de135c258eadb7b2032a2eba7b8113448` |

No source or witness correction was needed. The node-ledger definition
was repaired in the primary prose as recorded above. This audit accepts
the actual two-positive-meridian certificates and their whole-support
consumer for the three named curves. Transport to other curves requires
the separately audited higher-odd classification and connected-family
theorem; no broader Jacobian assertion follows from these finite witnesses.
