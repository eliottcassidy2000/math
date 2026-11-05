# Both paired-fan children from one persistent three-port exposure

This is a separate physical-interface upgrade to the audited paired fan. The original two-port packet is unchanged. The new statement is an all-order family of one-generation paired moves, not a returning ternary grammar.

## 1. Actual core, chain, marked root and common ports

Let the core be

`C=[26,16,3,6,12]`, with core vertex a=2.

For m>=1 add a transitive chain Z and a vertex b, with all edges `C -> Z -> b`, all edges `C -> b` except `b -> a`, and the forward order within Z. Label C by0,...,4, label b=5, and label Z by6,...,m+5. Call the graph P_m. Its count is

`H_m=18*2^m+13`.

Keep the SAME ordered actual ports `[0,1,4]` for every m and throughout every move below. The separate distinguished root b is always the old vertex5; it is not silently replaced by an equal-response vertex.

For a left birth L(mask,drop), the newborn loses to the selected current ports and beats every other old vertex. The candidate port list consists of the three old ports followed by the newborn. Dropping index3 forgets only the newborn's exposure; the graph vertex remains present.

The two available words are

- plus2: `L(0,3), L(1,3)`;
- plus4: `L(6,3)`.

Every individual birth retains all three original ports. Thus both choices start from the same marked/exposed graph, rather than requiring different initial pairs.

## 2. Exact graph identities

For plus2, first add a universal source u, then a vertex v which loses only to old vertex0. The second vertex beats u, every other old core vertex, every chain vertex and b. This is the same actual graph as the old paired-fan plus2 gate; only the exposed vertices have changed. Its new core is

`A=[90,16,3,6,12,31,62]`, still marked at core vertex2.

For plus4, add a vertex w which loses only to old vertices1 and4 and beats every other old vertex, including Z and b. Its new core is

`B=[26,48,3,6,44,13]`, again marked at core vertex2.

The output graphs are exactly the same fan construction based on A and B, respectively. Here is the explicit relabeling from the canonical fan labels to the actual birth labels:

- A fan: core0,...,4 stay fixed; its core vertices5,6 map to u=m+6 and v=m+7; its outer root7 maps to old b=5; its chain vertices8,...,m+7 map to6,...,m+5.
- B fan: core0,...,4 stay fixed; core vertex5 maps to w=m+6; outer root6 maps to old b=5; chain vertices7,...,m+6 map to6,...,m+5.

Every edge follows directly from the described births. These identities hold for every m, with the same words and the same original ports.

## 3. All-order counts and physical marked-root profiles

Write h=H(X), and let R_a(X) count ordered spanning two-path covers (U,V), with U possibly empty and V nonempty starting at a. Direct finite certificates for the three cores give

`(H(C),R_2(C))=(13,18)`,

`(H(A),R_2(A))=(41,54)`,

`(H(B),R_2(B))=(43,54)`.

For any nonempty core X and positive chain length m, put S=X->Z. It has H(S)=h, no Hamilton path ending at old a, and `R_a(S)=2^m R_a(X)`: each chain vertex may be assigned independently to either of the two path slots, in its forced increasing order. Appending b with sole outneighbor a gives H(S)-e_a(S)+R_a(S). Paths not using b->a are old Hamilton paths followed by b; those using it correspond exactly to the ordered covers counted by R_a(S). Hence

`H(F_m(X,a))=h+2^m R_a(X)`.

At the physical outer root b, the deletion count and ending count both equal h. Its root-starting two-cover count is h+R_a(S): the path starting at b is either the singleton b, or starts b,a and continues along the a-starting cover slot. Thus its rooted profile (H,e_b,R_b,H(T-b)) is `(H,h,H,h)`.

Applying this exact formula to C,A,B proves

`H(plus2(P_m))=54*2^m+41=3H_m+2`,

`H(plus4(P_m))=54*2^m+43=3H_m+4`.

The same old vertex b=5 has profiles `(H_m,13,H_m,13)`, `(3H_m+2,41,3H_m+2,41)` and `(3H_m+4,43,3H_m+4,43)` respectively. The ports stay `[0,1,4]` in every case. All three cores are strong, and each positive-chain fan is strong: from a one reaches the core and Z, then b, and b returns to a. Equivalently, direct quotient paths verify this property.

## 4. Why this is not an iterated digit grammar

The output core is A or B, not C. Retaining the same exposed vertices does not prove that the numerical rule applies to the outputs again. It does not: at m=1, repeating the plus2 word gives counts `49,149,515`, while repeating plus4 gives `49,151,409`. The second steps are not `3*149+2=449` or `3*151+4=457`.

Thus the theorem supplies both children for every member of one unbounded parent family from a single common persistent three-port exposure. It does not supply a closed state class, arbitrary-target inversion, or a universal ternary induction. The distinction between physical port return and return to the required numerical family is explicit.

## 5. Reproduction

Run `python3 verify.py` and `python3 -O verify.py` in this directory. The verifier checks the core profiles by independent subset DP and literal permutation/cut enumeration, verifies both relabeling identities at positive chain lengths1,...,6, directly recounts every graph and every old-root profile, checks the unchanged port list after each birth, and checks the explicit failures of naive repetition.

The all-order identities follow from the literal construction and the cover formula, not from those finite controls. The standalone graph helper is copied unchanged from the audited paired-fan packet; no previous frozen file is edited.
