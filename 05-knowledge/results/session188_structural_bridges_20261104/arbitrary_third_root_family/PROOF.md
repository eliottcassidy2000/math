# Arbitrary descendant families with two small cyclic-root children

This combines the root-family method in ../covered_root_family/PROOF.md,
the original-module decomposition theorem, and the arbitrary-third-block
identities in ../arbitrary_third_block/PROOF.md. It preserves each dependency's
scope and does not alter their frozen files.

Let C be any nonempty named tournament. Let T=C3[A,B,C], with A and B both
transitive of order m=2 or both transitive of order m=3. Let U result from
simultaneous edge-XOR of any family of modules of the original T that
contains the whole root with odd multiplicity. Reduce repetitions modulo
two. There is no bound on the number, depth, or overlap of proper selected
modules inside C or inside the other children.

**Corollary.** The original T and U always have different full named
two-port responses. A fixed battery of two actual queries suffices when
m=2:

    (0,0), (1_A,1_B).

For m=3 a fixed battery of three actual queries suffices:

    (0,0), (1_B,1_C), (1_C,1_A).

These batteries distinguish the original T from each root-containing XOR
target. They do not assert injectivity between every pair of targets in
an arbitrary original-module orbit.

## Proof

Every proper module of C3[A,B,C] is contained in a single child. Indeed,
if a module meets two children, any omitted vertex of the third sees those
two parts in opposite directions. It must therefore contain the third
child entirely. An omitted vertex of either occupied child then sees the
other occupied child and the fully included third child in opposite
directions. Thus no vertex can be omitted, and the module is the root.

Consequently the selected root operation converses the cyclic quotient
exactly once, and all proper operations remain inside its three original
children. Write the target as U=C3^op[A',B',C'] on the same physical blocks.

The original-module decomposition theorem gives H(U)>=H(T). If strict,
the actual query (0,0) separates, since it has both ports preceding all
old vertices and hence equals H. No assertion that crossing reversals
form a legal sequential operation is required.

If H(U)=H(T), the same theorem gives equality of all ordered-path-cover
profiles on every fixed physical induced subset. In particular

q_h(A')=q_h(A), q_h(B')=q_h(B), q_h(C')=q_h(C)

for every h. This is stronger than the hypothesis needed here; the target
children may have changed internal orientations in many ways.

A block-uniform two-port query depends only on the quotient orientation
and the three complete q profiles, by compression into maximal child runs.
Thus its actual T-versus-U contrast is exactly the corresponding contrast
for the quotient converse with all original interiors unchanged.

For m=2, the arbitrary-third theorem gives contrast

    -8H(C)-14q_2(C)+4q_3(C),

which is 4 modulo 8 if |C|>=2 and is -8 for singleton C. It is nonzero.
For m=3, the sum of the two prescribed root-query contrasts is

    242H(C)-283q_2(C)-552q_3(C),

strictly negative if |C|>=2 and equal to242 for singleton C. Hence at least
one of the two actual queries separates. Together with the scalar-strict
case this proves the corollary.

For a laminar selected family, all selected sets remain modules during
reversal, so scalar equality always holds. The scalar query may then be
omitted. This is a statement about the known laminar class, not permission
to omit the scalar check for arbitrary crossing families.

## Boundaries

The result concerns the whole cyclic core, or ordered earlier/later SCC
padding under the separately proved exact localization rule. It does not
assume that these integer or mod-eight root gaps transfer through an
arbitrary strongly connected outer quotient. It does not cover arbitrary
unequal child profiles, a prime quotient of arbitrary size, or arbitrary
legal sequences using new modules that were not original modules.

The proof method is the same scalar-strict/hereditary-equality split as the
companion equal-profile root-family theorem. The new input is the unequal
arbitrary-third-block separator, not a new claim about modular decomposition.

## Reproduction

The checker enumerates every root-containing original-module XOR family
for TT2,TT2,TT4 (including crossing intervals in the third child), and for
TT3,TT3,C3 (including crossing intervals in the first two children). It also
checks declared deep and prime third-block examples and the singleton
boundary. It constructs all targets on the same names and evaluates each
specified battery by literal Hamilton-path recursion.

    python verify.py > normal.json
    python -O verify.py > optimized.json
    cmp normal.json optimized.json
