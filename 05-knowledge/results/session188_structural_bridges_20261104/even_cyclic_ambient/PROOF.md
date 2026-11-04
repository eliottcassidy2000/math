# Even cyclic cores: arbitrary ambient and highest-node corollaries

These are two separate synthesis corollaries. Their inputs are the
original-module decomposition theorem, the highest-changed-node reduction,
the arbitrary-third-block separator, and the independently derived
even-order modulo-eight transfer theorem. None of those frozen packets
is changed here.

## 1. Dependencies and the actual lift

For M=C3[A,B,C], with A and B transitive of order two and C any nonempty
even-order tournament, the physical root query (a,c)=(1_A,1_B) separates
the quotient converse while retaining the original interiors. Its gap is

Delta F=-8H(C)-14q_2(C)+4q_3(C)=4 modulo8.                (1)

This is ../arbitrary_third_block/PROOF.md, Theorem1. Uniform port queries
depend only on the individual child q profiles, so (1) still holds if the
new child interiors merely have those same complete profiles.

The even-order profile theorem ../mod8_transfer/PROOF.md says that for
same-size even tournaments X,Y and any fixed rooted ambient quotient Q,

H(Q[X])-H(Q[Y]) = H(Q[TT2]) (H(X)-H(Y)) modulo8.          (2)

The multiplier H(Q[TT2]) is odd. Indeed Delta q_2=2 Delta H, all
Delta q_k for k>=3 vanish modulo8, and the nonnegative module-run expansion
has c_1+2c_2=H(Q[TT2]). A residue-four gap therefore remains residue four.

Suppose M is a module in an ambient tournament. Let

E={v outside M : v -> M}.

Give both new ports M's original exterior incidences. The actual lifted
query is

    a=1_(E union A),   c=1_(E union B),   x->y.           (3)

Then M together with the two ports is a module in both actual extensions.
Its order is |M|+2=6+|C|, which is even. Thus (2) applies to its complete
port cores and lifts (1) without dropping any cubic or quartic data.

## 2. Corollary A: root-containing families in every ambient quotient

Let T=Q[M] in any fixed rooted tournament Q, with M=C3[TT2,TT2,C] and
C nonempty of even order. Let U be formed by XOR-reversing a family of
modules of the original T, all contained in M, and containing M an odd
number of times. Repetitions are reduced modulo two. Proper chosen modules
may overlap or be deeply nested inside the children.

**Corollary A.** The full named responses of T and U differ. Two fixed
actual queries suffice to distinguish this original T from every such U:
the global scalar query (0,0), and the lifted root query (3).

Proof: every proper module of a cyclic three-child substitution lies in
one child, by the opposite-incidences argument in
../covered_root_family/PROOF.md. Therefore M's quotient is conversed
exactly once, and no exterior incidence changes. The original-module
decomposition theorem gives H(U)>=H(T).

If strict, F_U(0,0)=H(U)>H(T)=F_T(0,0). Otherwise its hereditary equality
conclusion gives equal complete q profiles on each original physical
child. Apply (1) to the local port cores and (2) to the fixed outer quotient.
The actual query (3) has F_T-F_U=4 modulo8, so it separates. This proves
the corollary even when the outer quotient is strongly connected. QED.

The even-order condition on C is essential to the cited transfer argument.
The corresponding odd-core formula has a Hamilton-cycle correction;
no odd-C arbitrary-ambient conclusion is inferred here. For a laminar
family, scalar equality is automatic and the single lifted root query
suffices.

## 3. Corollary B: exclude a highest changed prime C3 node

Let U be any XOR of modules of the original T, with no laminarity or
containment assumption. Use the persistent original strong-module tree
and its highest changed quotient nodes as defined in
../highest_quotient_reduction/PROOF.md. Here strong module means a module
overlapping no other module, not a strongly connected tournament.

**Corollary B.** Suppose a highest changed node M has original prime child
quotient C3, with two physical children of order two and a third child of
nonempty even order. Then F_T!=F_U.

This is a necessary exclusion for full-response collisions, not a claim
that every highest changed prime node has this form.

Proof: if H(T)!=H(U), the scalar query already separates. Otherwise
hereditary q-profile equality holds everywhere. At an originally prime
quotient the decomposition theorem permits only the original quotient
or its converse. Since M is changed, its quotient is conversed. Name the
two pair children A,B so that A->B->C->A originally. Every two-vertex
tournament is transitive, and each child retains its complete q profile.

All highest changed nodes are disjoint, all old arcs outside their interiors
are unchanged, and M's exterior incidences are therefore unchanged. Form
V by making only M's change inside T. Apply the physical lift (3) and (2)
to obtain F_T-F_V=4 modulo8.

For any other highest changed node N, its entire q profile is unchanged
by hereditary equality. In the actual two-port extension, N remains a
module: old incidences are uniform, and both ports see N as M does, hence
uniformly. Replace every such N's interior by its target interior. The
module-run substitution identity preserves the exact Hamilton count at
each replacement because the q profiles agree. These disjoint replacements
turn V into U, so F_V=F_U at the same physical query. Thus F_T-F_U=4 modulo8.
No other highest changed node can cancel the contrast. QED.

Multiple highest changed nodes, including other changed prime nodes, are
explicitly allowed. Descendant changes inside M need not be single
converses and need not be legal in their original listed order. The proof
uses hereditary profile equality, not such a legality assumption.

## 4. Scope and reproduction

Corollary A is original-versus-target separation for a fixed module and
fixed ambient boundary. It does not assert all-pairs orbit injectivity.
Corollary B gives a target-dependent highest-node obstruction, not one
universal fixed battery for every possible highest node in T.

Neither corollary settles arbitrary prime quotient size, arbitrary unequal
cyclic children, or arbitrary legal sequences using modules not present
in the original tournament. The local root identity, all-order mod-eight
transfer, and hereditary decomposition are the proof dependencies; finite
checks below are only implementation controls.

The verifier checks literal global scalar and lifted-root queries in every
declared ambient context, including strongly connected contexts. It checks
all root-containing original-module families in specified even cyclic cores,
including crossing descendants. It also constructs two disjoint highest
changed prime nodes and verifies that the residue-four query survives the
other node's converse.

    python verify.py > normal.json
    python -O verify.py > optimized.json
    cmp normal.json optimized.json
