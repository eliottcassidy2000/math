# Sevenfold scaling: a primitive replacement works, and a retained module cannot

This adds a constructive (p,7p) family and a contrasting obstruction to the [target-selection session](README.md). It does not give a universal selector. The distinction is between **replacing** a physical block with a different one and keeping the old block and adding uniform exterior contacts.

## A primitive two-reading replacement

Independent subset path counts give the following literal tiles (outgoing masks):

| Tile | masks | (H,P₂,P₃) |
|---|---|---|
| L₉ | (130,28,281,97,361,71,391,318,43) | (2197,34928,291327) |
| U₇ | (94,60,120,48,96,65,10) | (91,1184,8121) |
| U₈ | (2,8,19,53,195,87,143,47) | (169,3304,30249) |

Set B=L₉→TT₄ and B′=U₇→U₈ (ordinal compositions). Ordered path-slot counts multiply under ordinal composition because every slot concatenates in block order. Thus

    (H,P₂,P₃)(B)  = (2197, 558848, 23597487),
    (H,P₂,P₃)(B′) = (15379, 3911936, 245652129).

The first two coordinates scale exactly by 7. The lower pair is primitive: `H=13³` and `P₂=2⁸·37·59`. The third coordinate does **not** scale by 7, so there is no claimed three-reading or unrestricted iteration law.

## Literal strong witnesses on both prime axes

Let C₃[X,B,1] be a cyclic composition of X, B and one singleton, in the indicated orientation. For each block A write h_A=H(A), R_A=P₂(A)−2h_A (two nonempty slots). The singleton can occur once; X and B are each visited at most twice. There are three single-visit skeletons, one skeleton repeating X, one repeating B, and one repeating both. Therefore

    H(C₃[X,B,1]) =
       3h_X h_B + R_X h_B + h_X R_B + R_X R_B
       = (P₂(X)−H(X))(P₂(B)−H(B)) + 2H(X)H(B).

The composition is strong for any nonempty blocks, irrespective of whether the blocks are internally strong. Replacing B by B′ multiplies the read by 7 for every X. This proves paired physical counts for a variable reader, not just one matched scalar.

For X=TT_s, the lower count is `556651·2^s−552257`. At s=8 it is the prime 141950399; the two strong witnesses have orders 22 and 24. At s=16 the lower prime is 36480127679, with paired strong orders 30 and 32. The prescribed lower target on this reader is a two-coordinate line `p=556651 P₂(X)−552257 H(X)`, with coprime coefficients. When |X|=n≥2, the constraints include

    552257 H(X) ≡ −p (mod 4·556651),
    H(X) ≤ p/(556651n+4394).

These conditions leave a specific finite set of candidate profiles, but do not make any of them physical automatically.

A different exact family takes X a singleton and appends d ordinal singletons to **both** B and B′ before the cyclic wrap. This preserves their H coordinates and doubles their P₂ coordinates for each added singleton. Thus the strong counts are `p_d=558848·2^d+2197` and `7p_d`. For d=7,13,15 the lower primes are respectively 71534741, 4578085013, 18312333461, checked by trial division. Neither family asserts prime infinitude.

## No exact sevenfold host can retain the old block as a homogeneous module

**Theorem.** If M is a nonempty homogeneous module of a tournament T, then (H(T)
e7H(T[M])). Here homogeneous means every exterior vertex either beats all of M or loses to all of M.

First assume M is proper inside a strong tournament S. The quotient obtained by contracting M is strong: project paths and remove repetitions. Each vertex in a nontrivial strong tournament belongs to a directed triangle. Indeed its in- and outneighbor sets are nonempty; if every inneighbor beats every outneighbor, the outneighbors cannot reach the vertex or any inneighbor. So there is a triangle through the contracted vertex. The corresponding induced tournament contains C₃[M,1,1], with count (h+P₂(M)) by the preceding read. Induced extension never decreases H: insert a new vertex into each path at a legal boundary or transition, and different old paths cannot yield the same new one. Thus

    H(S) ≥ h + P₂(M) ≥ (|M|+2)h.

This first bound is sharp for each fixed M in an unrestricted strong host: use C₃[M,1,1] itself. If |M|≥6 it already exceeds 7h. For |M|=5 it also strictly exceeds 7h: P₂(M) is a multiple of four, by pairing complementary subsets in its cover sum and using Rédei oddness, whereas 6h≡2 mod 4. For |M|≤4, h can only be 1,3,5. If h=5 then |M|=4 and P₂=32. One way to see this uses the odd-pack formula: the two directed triangles giving h=1+2·2 contribute 8 each to P₂, in addition to 16 from the all-singleton pack. There cannot be more than two directed triangles on four vertices: the number of transitive triples is the sum of binomial(d⁺,2), at least two for outdegrees summing to six. Thus the strong wrapping floor at h=5 is 37>35. The cases h=1,3 would target the forbidden 7,21.

For arbitrary T, if M meets more than one strong component then it is a union of consecutive whole strong components. A partially met component would have a nonempty cut all in one direction, because a vertex of M outside the component fixes every omitted vertex's uniform relation to M. An omitted component between two met components would see M in both directions. Hence if M spans components (or equals a whole component), the ratio H(T)/h is an ordinal product of the outside counts. This cannot be 7, since 7 is prime and not itself an allowed factor. Otherwise M is proper within a single strong component and the strong bound above applies; multiplying by the other component counts does not weaken that lower bound.

So a universal sevenfold lift must replace the old block, as in the primitive pair here, or expose nonuniform contacts into it. Merely adding vertices that read the old block uniformly cannot finish that route. This says nothing against a nonhomogeneous induced copy.

Replay: `python replay_seven_balance.py`. It verifies the base tile profiles, exact prime examples in both wrappers, and direct literal checks of strong composed tournaments at small parameters. The all-order obstruction and transfer formulas are proofs.
