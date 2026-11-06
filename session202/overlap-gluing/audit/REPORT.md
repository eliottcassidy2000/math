# Independent audit: overlap gluing and the first two forward sectors

**PASS, without mathematical correction to the frozen proof.** The overlap-compatible path injection, zero/one-overlap product lower bound, sharp and hostile examples, and both forward-sector rigidity characterizations are valid. These are lower bounds and structure theorems, not an exact additive compiler or eventual-coverage result.

## Frozen inputs

Audited `PROOF.md` SHA-256:
`528959a2f978a62511b46e1a276f8375fab099b02b6151c4a518878bb1af657c`.
Audited `verify.py` SHA-256:
`81d215e7b016e5f6fb42d09bf15b26cfcd39674903e40622dc5b18da62e44725`.
Both remained unchanged during the audit. `INPUTS.json` also records the frozen author receipts.

## The overlap injection

The head-comparison merge is legal: after emitting one head, that vertex points both to its next private-path vertex and to the unchanged other head. Thus whichever head is next chosen is an out-neighbor. Empty input paths and exhaustion of one input cause no exception.

For a shared set S, compatibility fixes the same shared-vertex order in both inputs. Each gap contains two disjoint directed subpaths. Its merged first vertex is one of the original first vertices and therefore receives the required preceding shared arc; its merged last vertex supplies the required following shared arc. When both gap paths are empty, the shared adjacency was present in each input. Reinserting the shared vertices yields a genuine Hamilton path.

Restriction to X and Y recovers both original paths, including their private vertices' gap positions. This proves injectivity across different compatible pairs and across different orders of S. It does not require the shared-vertex order itself to be a directed path. The cases S empty, X=Y, X or Y empty, and one private part empty are valid, including the empty tournament convention.

The zero/one-overlap product conclusion follows because there is one order of S. It prohibits generic additive growth for independently unbounded constituent counts; it does not prohibit constructions that remove or alter constituent subgraphs or exploit larger overlap order sectors.

## Sharpness and the order-profile hostile

The five-vertex supplied graph is a valid strong tournament, has9 Hamilton paths, and induces the two specified triangles with3 paths each. Thus the product lower bound is sharp even for a strong output.

The six-vertex mask91 agrees exactly with rows[54,8,2,5,14,30] under the stated lexicographic pair/bit convention. Its arc0->1 has exactly7 forward and8 reverse Hamilton paths. Therefore majority order need not follow the edge. The proof correctly concludes that a missing sector profile prevents the proposed simplification, not that an additive two-overlap construction has been achieved.

## Forward-sector rigidity

An induced path extends while preserving its order by deterministic legal-gap insertion of each omitted vertex. Restriction recovers the original input, so sector lower bounds survive extension.

If a directed triangle omits at least one of u,v, its three paths merge with the fixed arc path u,v using overlap of size at most one. All outputs retain u before v. If it contains both endpoints, two of its paths retain that order. Since a nontransitive tournament contains a directed triangle, f=1 is equivalent to transitivity and its reverse sector is zero.

If f=2, every triangle must contain both endpoints. There is exactly one nontrivial SCC and it contains u,v; any other nontrivial SCC would supply a forbidden triangle. Every vertex of a strong tournament belongs to a directed triangle: otherwise its out-neighbors cannot reach its in-neighbors or the vertex itself. It follows that the exceptional SCC has the form u->v->W->u. The remaining W is transitive, since a triangle inside it would omit both endpoints. Other SCCs are singletons and provide only a forced source prefix and sink suffix.

In this SCC, a reverse-order Hamilton path must begin at v and end at u, so it is uniquely v,W,u. In a forward-order path, u must be immediately followed by v. Every vertex of the transitive W independently lies before u or after v, with the transitive order forced within each part. Therefore the sectors are(2^m,1), m=|W|>=1. Equality f=2 forces m=1, and the converse with transitive padding is immediate. Both directions and all boundary conditions are covered. The proof makes no claim for larger fixed forward sectors.

## Independent exact controls

The new `check.py` imports no author code. It constructs all1,100 labelled tournaments through order five **including the empty tournament**, recursively enumerates induced Hamilton paths, and checks all254,254 ordered covering pairs. Unlike the author's finite census, it includes empty private parts, containment and equality of the two subsets.

It independently verifies877,822 compatible-path merges, every restriction recovery and injectivity test, and the small-overlap product identity. On the exact author's restricted universe it reproduces187,620 cases and515,524 pairs.

For every oriented edge in that census, it verifies both structural equivalences using SCCs, not merely the total count3 for the f=2 case. There are1,364 forward-one and414 forward-two arc cases. It separately checks the sharp strong graph, the mask91 hostile, and the complete sectors(2^m,1) of the rigid strong family for1<=m<=10.

Independent normal/optimized outputs agree. Both author modes replay exactly and match their frozen receipts. All checks use explicit exception gates, so optimization does not disable validation. These finite controls support the all-order proof and are not an extrapolation argument.

Reproduce from this audit directory:

    python3 check.py > check-normal.json
    python3 -O check.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json NORMAL.json

No unconditional threshold for odd Hamilton-path counts follows from this packet.
