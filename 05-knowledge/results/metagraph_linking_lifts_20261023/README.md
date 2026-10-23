# Balanced port switches and their first visible electrical layer

October 23, 2026. **PROVED:** the marked four-terminal switch law and its
dyadic visibility criterion. **FINITE-EXACT:** all feasible integer multiples
of the previously authenticated order-eight port switch, their independent
electrical certificate checks, and the displayed linking charges. **OPEN:**
whether a small uniform family of such marked constraints characterizes
tournament metagraph wiring at every order.

## Inheritance and the precise question

The closest proved mechanism is the all-level filtration in
[the full order-eight electrical packet](../metagraph_order8_electrical_20261023/README.md).
Its literal integer lifts retain information lost by parity and by diagonal
energy alone. The hostile object is the actual eight-entry
[one-round port counterfeit](../metagraph_order8_20261023/README.md): it
preserves every stated source and target endpoint/color capacity, reflection
twins, and node tiling counts. The corrected near miss is treating a Smith
group or a list of form ranks as the marked bilinear forms themselves.
The least-used sidecar here is the *same terminal-labelled harmonic basis*
on both candidate wirings.

Incoming [Session38](../session38_incidence_ports_electrical_20261023/README.md)
independently reconstructs the full order-eight graph. The calculations here
use the earlier stable dense merged-node IDs and pin their edge and node
exports. No ID is silently transferred between the two enumerations.

## The all-size four-terminal law

Let A be the root-reduced Laplacian of any connected weighted graph. Select
four distinct nonroot vertices a,b,c,d and write

    u = e_a - e_b,       v = e_c - e_d.

Subtract conductance h from edges ac and bd, and add h to ad and bc.
Assume the resulting conductances remain nonnegative and the new graph
remains connected. Its reduced Laplacian is exactly

    A' - A = h (u v^T + v u^T).                       (1)

The changes to every diagonal entry and every row sum are zero. Thus this
switch preserves vertex degrees, independently of any tournament
interpretation. For arbitrary integer potentials x,y,

    x^T(A'-A)y = h [(u.x)(v.y) + (v.x)(u.y)].         (2)

These assertions follow by expanding the four edge Laplacians. Negative
edge conductance changes give positive off-diagonal changes, which fixes
the sign in (1).

Suppose h = 2^j s with s odd and j >= 1. Define the common spaces

    V_k = {x mod 2 : A x is divisible by 2^k}, k <= j,

and use the previously proved forms q_k(x,y)=x^T A y / 2^k modulo two,
with admissible integer lifts. Since A'=A modulo 2^j, both matrices have
exactly the same admissible lifts through level j. Consequently:

1. Their marked spaces V_1,...,V_j and forms q_1,...,q_(j-1) agree.
2. At the first possible visible level,

       q'_j - q_j = alpha tensor beta + beta tensor alpha,       (3)

   where alpha(x)=u.x modulo two and beta(x)=v.x modulo two,
   restricted to V_j.
3. Every diagonal value q'_j(x,x)=q_j(x,x) still agrees.
4. The difference in (3) has rank two if alpha and beta are independent,
   and is zero otherwise.

For the last statement, independent alpha,beta map V_j onto F_2^2, where
the form has matrix [[0,1],[1,0]]. Dependent binary functionals are either
zero or equal and give zero. Thus visibility is **equivalent** to the two
terminal-difference functionals remaining independent on V_j. Two marked
potentials with (alpha,beta) values (1,0) and (0,1) give the cheapest
decisive test. This is a complete criterion for this particular switch
at its first possible layer, not a universal wiring classification.

Equal ranks of q_j and q'_j do not imply equality of the forms: adding a
rank-two bilinear form can preserve the rank of the whole form. The native
order-eight doubled switch below realizes exactly this failure.

## The entire feasible native port-switch interval

Use a=1238, b=2242, c=2246, d=2315. In literal black counts the inherited
switch is

    ac:-2t,  ad:+2t,  bc:+2t,  bd:-2t.               (4)

Its eight ordinary-port edits have coefficients +/-t. In each reflected
four-entry block the original counts at ac,ad,bc,bd are 2,2,1,3.
Nonnegativity is therefore equivalent to

    -1 <= t <= 2,       t an integer.

Every source-row and target-column change sums to zero. The converse
reflection constraint is linear and was already authenticated for the
inherited switch. Thus **all four values** -1,0,1,2 retain the original
one-round observer's capacities and reflection constraints. This does
not make the three nonzero choices actual tournament metagraphs: the
stronger incidence observer already uniquely specifies t=0.

The four terminals are blue-isolated, so their unit potentials lie in
V_1. For t=1 the direct cross-energy q_1(e_a,e_c) changes by one, while
every diagonal half-energy remains fixed. In particular all source
bucket totals and all diagonal parity tests can agree on a false wiring.

## Independent full-network results

The same complete 3,528-node graph is used in every row. The literal
integer certificates are independently checked using sparse matrix rows,
highest-pivot binary elimination, and the actual saved harmonic lifts.

| t | Exact 2-primary critical group | v_2(tree count) | First changed marked form |
|---:|---|---:|---|
| -1 | (Z/2)^3360 + (Z/4)^2 + Z/32 | 3369 | q_1 |
| 0 | (Z/2)^3358 + (Z/4)^4 + Z/256 | 3374 | native |
| 1 | (Z/2)^3360 + (Z/4)^2 + Z/64 | 3370 | q_1 |
| 2 | (Z/2)^3358 + (Z/4)^4 + Z/8 | 3369 | q_2 |

The t=2 counterfeit preserves the *entire marked half-form*, not merely
its rank. It also preserves the literal mod-four lift basis. In that
common five-vector basis the native and counterfeit quarter-form rows,
encoded as binary integers, are

    native:    [8,12,14,31,8],
    t=2:      [12,12, 7,27,8],
    difference:[4, 0, 9, 4,0].

The difference has rank two and zero diagonal; both total forms have
rank four. Hence their next radical has dimension one in both cases,
but it is a different marked subspace: in the displayed common basis
the native kernel is generated by binary integer 17 (coordinates 0+4),
and the counterfeit kernel by 3 (coordinates 0+1). At the following level the t=2
mode has nonzero q_3 and dies. The native mode survives until q_8.
Keeping ranks while discarding the marked radical would miss the first
decisive obstruction.

For any highest primitive lift x at level e, put y=Ax/2^e. In the
discriminant group coker(A), [y] has exact order 2^e: A^(-1)y=x/2^e,
and x has an odd coordinate. Its linking self-pairing is

    b_A([y],[y]) = x^T A x / 2^(2e) modulo Z.

The supplied generators give self-pairings 15/32,105/256,49/64,1/8
for t=-1,0,1,2 respectively. These are generator-specific rational
numbers, not invariant distinguished labels. The separate Hall/Bott
packet develops the marking and normalization distinction.

## Reproduction and scope

```text
python -B 05-knowledge/results/metagraph_linking_lifts_20261023/generate.py
python -B 05-knowledge/results/metagraph_linking_lifts_20261023/check.py
python -B -O 05-knowledge/results/metagraph_linking_lifts_20261023/check.py
```

`generate.py` uses the pinned previously audited filtration engine and
stores literal lifts in deterministic `certificates.json.gz`. The consumer imports **no
producer code**. It uses the pinned independently written sparse consumer
API from the previous packet, reconstructs the native sparse Laplacian,
applies the four off-diagonal changes directly, and verifies every image
divisibility, parity span, full form, rank, terminal layer, and charge.
It also checks the exact feasible port interval, degree preservation,
the cross-energy identity, and the doubled switch's identical half-form.
Native t=0 is the positive control. The three nonzero choices are hostile
controls to the one-round observer; t=2 is additionally hostile to a
half-form-only or quarter-rank-only observer.

The certificate does not rerun the inherited tournament enumeration.
It binds the stable graph and port witness by LF-normalized hashes, and
consumes their already independently authenticated finite universe.
Normal and optimized consumer replays each pass **465,559 checks**;
the separate small boundary controls pass **16 checks**. The decompressed
literal-certificate SHA256 is
`be1639d01240dd3a5bcfe84cef0f05db2d262e9c11507339680ea19eec59d4ae`.

## Why even the complete dyadic tower is not a universal wiring certificate

The small `boundary_controls.py` packet supplies two exact hostile controls
outside the native tournament family. They stop an unjustified universal
promotion without weakening the native switch result.

First fix the blue path 1-2-3-4. Put one paired-black line on each of edges
12 and 34, or instead on 14 and 23. Both literal graphs have labelled
degrees (3,4,4,3). Their spanning-tree counts are respectively 9 and 17,
verified by cofactor determinants and independent weighted tree enumeration.
Since both counts are odd, V_1=0 and **every** higher binary form and
binary harmonic character is trivial for both graphs. The entire dyadic
tower is blind to this balanced black switch. In this graph class an odd
prime probe is necessary to distinguish the two critical groups. In fact
one ternary potential suffices: x=(1,0,0) on the reduced three coordinates
has native image (3,-3,0), but switched image (3,-1,0). It is harmonic
modulo three only for the native graph. This is a concrete reason to use
both binary and ternary registers, not a numerological identification.

Second, take two unit trees with edge sets

    {12,23,34,15,16},       {12,24,13,35,16}.

Both have the same labelled degrees (3,2,2,1,1,1), but they are not
isomorphic: distances from their unique degree-three vertex to the leaves
are (1,1,3) and (1,2,2). Both determinants are one. Hence even all primary
critical groups and all marked discriminant pairings are trivial. Full
electrical pairing data do not recover graph structure without another
sidecar; loops are also always invisible to a Laplacian. These controls
are statements about graphs, not assertions that either can be a native
tournament metagraph.

The difference between odd and even primes is itself structural. Over an
odd prime field, the diagonal quadratic values determine a symmetric
bilinear form by polarization and division by two. Over F_2 that division
is impossible, and the switch cross-form (3) has zero diagonal despite
being nonzero. This explains the particular value of retained joint
incidence and cross-energy in the even/odd construction.

```text
python -B 05-knowledge/results/metagraph_linking_lifts_20261023/boundary_controls.py
```

The portable conclusion is to retain **marked cross-energy and marked
radicals** along with local incidence capacities. The groups are useful
rejection invariants, but even all their orders and ranks discard which
currents belong to which original ports. Conversely, energy or linking
data alone do not reconstruct a graph: unit trees already all have trivial
critical groups. A universal small certificate remains open.
