# Actual C3 kernel phases as affine ternary codes

24 October 2026. **PROVED, scoped all-size:** affine-fibre classification, exact budget components, constructive bounded-support routes, exact no-return route criterion, three-mode bridge family, and its selected-graph spectrum. **FINITE-EXACT:** the closed checker reconstructs actual owner records and matrices, then checks the stated finite universes. This is an extension of the frozen native equality carrier, not a claim that every linear relation follows from pairwise owner equality or from one global observer.

## 1. Native records, resources, and complete guard scope

Use the actual source S_r=T_r -> C3, with q=8r and r>=1. The rigid eight-owner seed has canonical code29558785 and prefix paths

    P0=(3,0,6,7,1,5,4,2), P1=(3,1,6,7,0,5,4,2).

The endpoint-fixed partner isomorphism is g=(4,0,2,3,1,5,6,7). For each labelled slot i choose a fixed distinct prefix key b_i in {0,1}^r. The complete prefix-pinned menu consists of exactly

    R_i(k)=P_(b_i), (q+k,q+(k+1 mod3),q+(k+2 mod3)),
    k in F3.

There are n slots, with n<=2^r. The terminal owner at the actual indexed clock q determines the phase literally:

    k_i=owner_(R_i)(q)-q in F3.                         (1)

The subtraction is part of the phase chart. If o_i=owner_(R_i)(q) is the raw owner label, then o=q·1+k, so the canonical-phase equation A k=b is equivalent to

    A o=b+q A1 in F3, q=8r.                            (1a)

Using the same right side b on raw owner residues is generally false. For r=1 and row A=(1,0), canonical right side0 means raw right side2: omitting the shift selects first phase1 rather than0. Pairwise difference rows have A1=0, and r divisible by3 has q=0 in F3; those special cases can conceal the error. All raw-owner arithmetic here uses the declared suffix numbering q,q+1,q+2. Under a different owner embedding the phase chart must be transported explicitly; arbitrary owner labels are not automatically F3 coordinates. This raw-owner translation q A1 is separate from the converse-head phase translation2 A1 in section7.

Aut(S_r)=C3 fixes the rigid prefix and rotates the final cycle. Three phase values at one key consume ONE normalized resource, not three. Distinct keys have disjoint four-word source/K/converse reservoirs. Thus every phase assignment on distinct fixed keys is capacity-valid, provided any immutable background avoids those keys. Repeated keys reject the assembly regardless of phase.

Every nonzero phase change at one slot changes exactly three owner-clock occurrences and keeps that same resource continuously occupied. Literal suffix rotations authenticate source and partner transports. Endpoint-fixed copies of g authenticate all prefix partner gluings, so the target isomorphism class is the same for every supplied record. These are actual source and target matrices, not names assigned by fiat. The inherited all-size argument uses ordered strong components, the rigid seed, and endpoint-fixed gluing; finite reconstruction does not substitute for that argument. Source copies and their provenance are included in inputs/ and Sourcearchives.txt.

Now supply a matrix A in F3^(p by n) and b in F3^p, and impose the full-record predicate

    A k=b, with k recovered by (1).                    (2)

This is a declared, genuine relational guard on actual records. Its rows may have arbitrary arity. It is NOT asserted that arbitrary A,b can be compiled from conjunctions of pairwise owner equalities/pins. That smaller class is precisely the voltage/equality subclass in the inherited note. For example k1+k2+k3=0 has nine states and every two-coordinate projection is all nine pairs; no equality/pin system on those three slots has exactly this relation. Nor is (2) advertised as the output of a single target observer.

Retain every other required guard. If the full retained guard factors through the fixed source/target class-incidence data, its entire value is constant on this family: either it accepts every state of (2), or the whole family is empty. This includes non-additive functions of that fixed data, not just a modular or scalar charge check. Phase-sensitive or owner-sensitive guards must be proved constant, appended as additional authenticated constraints, or left outside this theorem. No order-eight numerical harmonic bank is transported to all orders8r+3; no unrestricted electrical-current claim is made.

Atomic packet edges replace a set of full records between accepted arrivals. They are not asserted to be microscopic edges of an earlier metagraph. The graph and its spectrum below are phase-packet graphs on actual supplied records.

## 2. Exact affine-code classification

First perform capacity, constant-guard, and consistency checks. If rank(A|b)>rank(A), there are no accepted assemblies. Otherwise choose one solution k0 and set

    C=ker A, d=n-rank A, F=k0+C.                       (3)

**Affine classification theorem.** The accepted full assemblies are in bijection with F, with precisely3^d states. For any x,y in F, their exact changed-slot cost is wt(y-x), their exact changed-owner-clock cost is3 wt(y-x), and their resource assignment is unchanged.

Proof: (1) is a bijection between each slot's complete menu and F3. Equation(2) is the actual retained predicate, so linear algebra gives(3). At any changed slot all three terminal owners change and no prefix owner changes. Fixed distinct keys prove capacity and continuous resource occupancy at every arrival. This proves both directions, not merely a quotient obstruction. QED.

In particular any finite affine ternary code has such a native full-record realization after choosing r with2^r>=n, with precisely this declared higher-arity predicate. This is a statement about a supplied relational guard, not a native realization of arbitrary pairwise equality, observer, or wiring gadgets.

### Weight enumerator, not a shortest-path enumerator

For any accepted reference x, define W_C(z)=sum_(v in C) z^wt(v). Then the exact joint endpoint-mismatch enumerator is

    sum_(y in F) z^dist_slots(x,y) t^dist_clocks(x,y)
      = sum_(v in C) z^wt(v) t^(3 wt(v))
      = W_C(z t^3).                                  (4)

The number of ordered pairs at slot distance j is3^d times the coefficient of z^j. Equation(4) is the full-record Hamming distance distribution, NOT the shortest-path-distance distribution in a smaller budget graph. An empty family has no accepted reference. A singleton has W=1.

## 3. Full budget graph, exact components, and constructive routes

For a nonnegative integer slot budget B define

    S_B={v in C:0<wt(v)<=B}, H_B=span_F3(S_B), rho_B=dim H_B.

Join distinct accepted arrivals precisely when their difference is in S_B. This is the Cayley graph of the additive group C, translated by k0. It is undirected because S_B=-S_B.

**Budget theorem.** Its components are exactly the affine cosets x+H_B. Thus

    number of components =3^(d-rho_B),
    size of every component=3^rho_B.                   (5)

Two specified endpoints x,y are connected iff y-x belongs to H_B. For an owner-clock budget T replace B by floor(T/3). A singleton is connected at budget0; an empty family is separately rejected and has zero components.

Proof: every packet difference is in H_B, so a route stays in one coset. Conversely a representation y-x=sum_j c_j v_j with v_j in S_B and c_j in F3 is realized by adding each nonzero c_j v_j. Coefficient2 does not change support, and every partial sum is in C. Adding it to x therefore remains in the actual affine accepted space. This gives the reverse implication and(5). QED.

The global connectivity threshold is the least B for which rho_B=d, equivalently the minimum, over bases of C, of the largest support of a basis vector. The endpoint threshold is the least B with y-x in H_B. These are generally neither minimum code distance nor maximum codeword weight. Dimension alone does not determine either threshold: a one-dimensional coordinate code connects at1, while a length-n repetition code requires n.

### Algorithm without enumerating the whole affine fibre

Enumerate the sum_(j=0 to min(B,n)) binom(n,j)2^j ambient phase vectors of support at most B; omit zero for generators. Keep those satisfying Av=0. Use elimination to retain an independent subset of the ORIGINAL low-weight vectors, and calculate its span H_B. Row-reduced output rows need not remain low-weight, so they must not silently replace the actual move generators.

Solve for k0, determine d and rho_B, and use(5). To route x to y, solve for its difference in that retained basis. If no solution exists, their quotient labels in C/H_B differ. If it exists, apply each nonzero coefficient times its basis vector once. This uses at most rho_B events, each at most B slots, with at most3B rho_B owner-clock work. Every step has an authenticated source-C3 action at its changed slots, and every arrival satisfies the original guard and capacity. For fixed B this bounded-support enumeration and field elimination is polynomial in n and the input matrix size; it does not enumerate3^d states. No efficient all-B weight-enumerator algorithm is claimed.

The linear span is sufficient here precisely because the full state set is an affine F3 space and every partial sum is genuinely accepted. The analogous span/lattice test can fail for integer trades constrained inside a binary or capacity box. No such box is being dropped here: the native resources are fixed and distinct, and all three phases are retained at every slot before the explicit linear guard.

These basis routes may return to an earlier FULL RECORD at a slot. They never release or reacquire a normalized RESOURCE: each slot holds its one fixed resource throughout. These are different conditions.

### Full-budget character spectrum

Let chi range over all characters of C (equivalently the3^d elements of its dual group). The entire full-budget phase-configuration Laplacian spectrum is

    lambda_chi=|S_B|-sum_(v in S_B) chi(v).             (5a)

The affine translation by k0 does not affect it. Because S_B=-S_B and the only character values are1,omega,omega^2, let c_j count the generators with chi(v)=omega^j. Then c1=c2 and

    lambda_chi=3c1=(3/2) #{v in S_B:chi(v)!=1}.

Thus every eigenvalue is a nonnegative integer divisible by3. It is zero exactly when chi is trivial on H_B. There are exactly3^(d-rho_B) such characters, agreeing with(5). This is the Laplacian of the supplied phase-configuration packet graph, never the original tournament metagraph Laplacian.


## 4. Exact no-return routes and a sharp separation

Use the strict per-slot condition: after leaving a full record R_i(k), the route never visits it again. Continuous occupancy of the normalized key is still allowed. This is stronger than merely forbidding repetition of an entire assembly.

For endpoints x,y let delta=y-x and S=supp(delta). An endpoint-fixed slot i outside S cannot change at all on a no-return route. For i in S, it either changes directly from x_i to y_i once, or visits the unique third phase x_i-delta_i before reaching y_i. Thus its sequence of nonzero phase increments is either

    (delta_i), or (-delta_i,-delta_i).                 (6)

This is a necessary and sufficient coordinatewise characterization. In particular every no-return route changes a differing slot at most twice, and its total owner-clock work is3(|S|+h), where h is the number of differing slots using the intermediate phase. A direct endpoint packet always gives a no-return route with one event, budget|S|, and minimal possible total work3|S|. For identical endpoints the empty route is optimal.

For exact bounded-budget questions form a finite directed acyclic graph. Give every i in S a stage0,1,2 corresponding to x_i, x_i-delta_i, y_i. Keep only stage vectors whose assembled phases satisfy A k=b; all other slots stay fixed. Join distinct stage vectors when every stage weakly increases. Edge cost is the number of changed coordinates. The start is all0 and the endpoint all2.

**No-return theorem.** No-return full-record routes are exactly the directed paths in this DAG. The minimum no-return budget is its minimum bottleneck path value. Under a fixed budget B, deleting edges of cost>B and using ordinary DAG shortest paths yields the exact minimum number of events. Additive edge cost instead yields minimum owner-clock work after multiplying by3. A route has at most2|S| events. The graph has at most3^|S| vertices; this is an exact finite characterization, not a polynomial-time claim.

Proof: the phase menu has exactly three records. With no return, an endpoint-fixed coordinate cannot leave, and a differing coordinate can use only the two sequences above. These sequences give the stages and a monotone DAG path. Conversely every retained stage vector is accepted and every monotone edge follows those permissible sequences. Costs are literal support costs. The stage sum increases by at least one per event, bounding length by2|S|. QED.

A necessary lower bound is obtained by shortening C to vectors supported in S, then applying the low-weight-span condition there. Additional stage monotonicity is retained by the exact DAG; it is not erased by a generic spanning argument.

### Five resources: connectivity at3, no-return endpoint budget4

Let

    v=(1,0,0,1,1), w=(0,1,1,0,1), C=span{v,w}.

Equivalently impose k2=k3, k1=k4, k5=k1+k2, using one-based coordinates. These are explicit full-record linear relations. There are nine accepted states, realized in S_3 on five distinct resources. The only nonzero codewords of weight at most3 are ±v,±w. Hence the budget3 graph is connected, while budget2 has nine singleton components.

Consider x=0 and y=v-w=(1,2,2,1,0). The two-packet route0,v,v-w has budget3 and never releases any resource, but its fifth full record follows phases0,1,0. Every allowed nonzero budget3 first move changes that fifth phase. Since its endpoints are both0, no no-return budget3 route exists. The direct packet0->v-w has weight4 and is no-return. The sharp endpoint no-return budget is therefore4, strictly above its ordinary connectivity threshold3. Counting full lifts as separate resources would obscure rather than repair this obstruction.

## 5. Three genuine helper modes on2L+1 resources

For L>=1 use old slots(a_j,b_j), j=1,...,L, and an actual helper phase m. All2L+1 prefix keys are fixed and distinct; choose r so2^r>=2L+1. Retain ALL THREE helper records and impose

    b_j-a_j=m in F3 for every j.                       (7)

This is a homogeneous linear full-record guard on the canonical phases extracted from actual terminal owners, of dimension L+1. In raw owner coordinates the same relation is o_bj-o_aj-o_m=-q modulo3, not generally0. It has exactly3^(L+1) states. It uses a genuine ternary relation, not a claim about pairwise old/helper equality. No phase is discarded from the helper menu.

At fixed m each a_j varies freely, and changing one a_j changes both old slots in its pair, at exact cost2 slots/6 clocks. Every old slot has all three phases in every mode. Changing m by delta!=0 changes the helper and at least one member of each pair. Its minimum cost is L+1 slots/3(L+1) clocks; it is attained independently in each pair by choosing delta a_j=0 or -delta. There are2^L minimum bridges for each of the two helper increments.

Consequently, for L>=2:

- Budgets0 and1:3^(L+1) singleton components.
- 2<=B<L+1: exactly three components of size3^L, one per helper phase.
- B>=L+1: the entire family is connected.

The within-mode budget2 graph is the Cartesian product of L triangles. For L=1 the mode bridge already costs2, so budget2 is connected; the three-component statement intentionally excludes that case.

Its complete kernel weight enumerator is

    W_L(z)=(1+2z^2)^L +2z(2z+z^2)^L.                  (8)

The first summand fixes m=0. For each of the other two m values, the helper contributes z and each pair has two weight1 assignments and one weight2 assignment. This proves(8) without numerical extrapolation.

For L=2, S_3 has order27 and five resources; all27 accepted states survive. The enumerator is

    1+4z^2+8z^3+12z^4+2z^5.

Budget2 gives three nine-state components; budget3 is connected and regular of degree12.

## 6. Exact spectrum of a selected phase graph

For general L define a particular Cayley graph on the parameter space(a_1,...,a_L,m) in F3^(L+1). Retain only:

1. The2L within-mode moves changing one a_j by ±1, at slot cost2.
2. All2^(L+1) minimum-mode bridges: delta m=delta in{1,2}, and independently delta a_j in{0,-delta}, at slot cost L+1.

Its degree is d_L=2L+2^(L+1). For L<=2 this is the full budget-(L+1) graph. For every L>=3 it is NOT the full budget graph: a simultaneous same-mode change in two pairs costs4<=L+1 and was not selected. At that budget the missing edges are exactly the same-mode moves changing2 through floor((L+1)/2) pairs; additional bridge patterns enter only at larger budgets. Keeping this distinction is essential.

Let omega be a primitive cube root and chi_(r,t)(a,m)=omega^(r dot a+t m). Write s=wt(r), q=t+sum_j r_j modulo3. Within-mode adjacency has eigenvalue2L-3s. The bridge contribution is

    C(r,t)=sum_(delta=1,2) omega^(delta t)
                     product_j (1+omega^(-delta r_j))
           =2^(L-s)(-1)^s * (2 if q=0, else -1),       (9)

because1+omega^(-delta r_j)=-omega^(delta r_j) when r_j!=0. Thus the combinatorial Laplacian eigenvalue is

    lambda(r,t)=3s+2^(L+1)-C(r,t).                     (10)

Equivalently its complete multiset, merging coincident values, is

    A_s=2^(L+1)+3s-(-1)^s 2^(L-s+1),
       multiplicity binom(L,s)2^s;
    B_s=2^(L+1)+3s+(-1)^s 2^(L-s),
       multiplicity2 binom(L,s)2^s,
    for s=0,...,L.                                    (11)

There is one q=0 character and two q!=0 characters for each r. This proves the multiplicities; their sum is3^(L+1), and A_0=0 is simple. For L=2 the exact spectrum is

    0^1, 9^8, 12^6, 15^12.                            (12)

This example's actual27-by-27 integer Laplacian also satisfies L(L-9I)(L-12I)(L-15I)=0; its traces of powers0,1,2,3 are27,324,4212,56700. These exact identities independently fix the stated multiplicities.

### Exact gap, with locality and degree costs retained

For every L>=1 the selected graph's Laplacian gap is

    gap_L=3+3*2^(L-1).                                (13)

Its multiplicity is4L for L>=2, and6 for L=1. Indeed the odd-s branch B_s is minimized at s=1 and increases strictly along odd s. The positive even-s branch A_s is minimized at s=2 and A_2=B_1+3. The other positive-s branches exceed2^(L+1)+3 and cannot be smaller. Finally B_0-B_1=3(2^(L-1)-1)>=0, with equality only at L=1.

The normalized Laplacian gap satisfies the sharp bound

    gap_L/d_L >=27/40,
    equality exactly at L=4.                          (14)

After multiplying by the positive denominator, this is2^L-9L+20>=0. It holds by checking L=1,...,4 and observing that its successive difference2^L-9 is positive for L>=4. This is a property of the DESIGNED PHASE GRAPH. Its move budget grows as L+1 and its degree grows exponentially. It is not a free bounded-local-move mixing theorem, not an original-metagraph harmonic spectrum, and not an authenticated all-order electrical observer.

## 7. Coordinates and converse must remain actual

A purely mathematical affine monomial reparameterization k=D P z+c transforms the guard into

    A D P z=b-Ac.                                     (15)

It preserves coordinate Hamming costs. General GL(n,3) transformations need not: (u,v)->(u,u+v) sends a weight1 difference to weight2. Equivalence of abstract vector spaces is therefore not enough for this geometry.

Equation(15) is an exact change of coordinate names, not automatic evidence of a native owner symmetry. Literal translations of a phase have the authenticated terminal source rotations above. Relabelled slots carry all their actual keys and records with them. A suffix reflection j->-j is NOT a same-source automorphism of directed C3, and arbitrary sign choices or reflected edge holonomy must not be invented from the code notation. Any native symmetry claim must additionally supply the actual owner permutation, clock map, source/target matrices, guard transport, and resource map.

Converse is authenticated directly. With N=q+3, reverse the path in the transposed source: old clock q+t becomes2-t. The original phase is still recovered at reflected clock2, so the same A k=b applies. If instead the converse chart uses its first owner, ell_i=owner_i^reverse(0)-q, then ell_i=k_i+2. The actual transformed relation is

    A ell=b+2 A 1.                                    (16)

This is a verified affine translation, not a guessed sign reflection. The partner square is K(S_r^op,R_i^reverse)=K(S_r,R_i)^op on the literal owner labels. Distinct resource reservoirs, packet support, and the factor3 clock cost survive. The reversed suffix sits at the FRONT of the converse path; it is not silently treated as the terminal suffix of the original source.

## 8. Reproduction and loss ledger

Run python -B check.py, python -O -B check.py, and python -B verify_manifest.py. No external packages, network, repository, or parent-package imports are used. The checker uses explicit exceptions, never Python assert. The deterministic receipt states the finite universes and exact outputs. The manifest closes every shipped file and rejects unlisted files except bytecode caches; portability is additionally tested after copying only manifest-listed files and the manifest into a fresh /tmp directory.

Preserved by the native-record-to-phase map: the entire supplied record via its fixed prefix sidecar and phase, fixed normalized resource, full linear guard, accepted arrivals, slot cost and actual owner-clock cost. The affine-code map forgets literal labels only if those sidecars are discarded. Keeping merely dimension loses the weight geometry; keeping merely unrestricted span connectivity loses per-slot no-return; keeping only local phase projections loses higher-arity consistency. The all-size results concern the declared phase-packet family and guards that truly factor through its fixed incidence. Pairwise equality universality, a single global observer, arbitrary native affine reflections, microscopic metagraph routing, and transported all-size harmonic/current banks remain outside scope.
