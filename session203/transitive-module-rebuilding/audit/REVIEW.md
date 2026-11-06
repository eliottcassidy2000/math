# Independent audit of the three transitive-module repair families

6 November 2026. **PASS in the stated scope.** These are genuine unbounded families of physical H-to-2H+1 rewrites. They do not establish an eventual spectrum threshold, a complete residue class, or a closed induction for arbitrary represented integers.

## Inputs and independent method

The reviewed producer PROOF.md, after correcting its two reproduction commands to python3, has SHA-256 c668cd4f50188050cdf53affd139c4001294c9176532ac273c85c9c6bf6e8596. Its frozen bytes are REVIEWED_PROOF.md. No mathematical correction was required.

The independent check.py imports no producer code. It validates the literal tournaments, core masks, and exact directed old-arc reversals. After inflating a module, it explicitly identifies all old core and module vertices inside the child and verifies that exactly the named old arc changed. For B and C it verifies that vertex 0 is a universal sink at every module size and that its deletion leaves the parent order unchanged.

The coefficient checker uses a different calculation from the producer's product of gap polynomials. It classifies core orderings by the number of allowed and mandatory module gaps. If a is the number allowed and t the number mandatory, the number of ways to occupy exactly r gaps is binomial(a-t,r-t). Summing these counts gives every skeleton coefficient. Cyclic orders fix core vertex 0 and include the wraparound gap, removing rotation exactly once.

Separately, recursive literal enumeration of paths, without subset DP or memoization, verifies the path and cycle counts of 48 inflated graphs at k=0,...,5. This visits 1,126,340 path prefixes. Cycles are counted only from fixed vertex 0, so their normalization is modulo rotation. The normal and optimized Python receipts are byte-identical, and all checks use explicit exceptions.

## Why the finite coefficients prove every module size

Every path has a unique core order and a unique set of occupied module gaps. Assigning the k distinct module vertices to r nonempty ordered gaps gives r!S(k,r) possibilities. Every gap's internal order is forced by transitivity. There is no extra factorial, and the assignment is bijective. The same argument works for cycles anchored at a core vertex, which exists in every graph here.

The finite enumeration therefore gives exact coefficients in a proved all-k formula. It is not polynomial interpolation from the literal sample sizes. The constant skeleton coefficient contributes only at k=0. This distinction is essential for families A and C.

## Verified coefficient lists

All omitted trailing coefficients are zero.

- A parent paths: (25,59,41,8); child paths: (51,119,82,16).
- A parent cycles: (2,3,1); child cycles: (3,5,2).
- B parent paths: (25,45,20); child paths before and after sink deletion: (51,91,40).
- B parent cycles: (1,1); intermediate child cycles: (0); replacement cycles: (3,3).
- C parent paths: (25,53,27); child paths before and after sink deletion: (51,107,54).
- C parent cycles: (1,1); intermediate child cycles: (0); replacement cycles: (3,3).

In every path pair, child minus twice parent is the skeleton vector (1,1). Its evaluation is exactly 1 both at k=0 and k>0. Hence H(child)=2H(parent)+1 for every nonnegative module size.

With delta equal to 1 at k=0 and 0 otherwise, the independent conversion to exponential coefficients gives:

    A parent = 8*3^k + 17*2^k + 1 - delta
    A child  = 16*3^k + 34*2^k + 3 - 2*delta
    B parent = 20*2^k + 5
    B child  = 40*2^k + 11
    C parent = 27*2^k - 1 - delta
    C child  = 54*2^k - 1 - 2*delta.

The cycle coefficients give A's counts 2^k+1 and 2^(k+1)+1, and B/C's parent/replacement counts 1 and 3. Thus all stated parents and final children are strong. The zero-cycle intermediate children in B/C are weak, as required by their universal sinks.

## Physical meaning and limits

A changes one specified old arc and adds one vertex, retaining every old vertex. B and C allow deletion explicitly: reverse an arc incident to old vertex 0, add the new vertex, then remove the resulting universal sink 0. On the final graph this is replacement of one old vertex, leaving all arcs among retained old vertices unchanged. No minimum old-arc edit count is inferred when deletion is allowed.

Deleting a universal sink preserves H bijectively, since every spanning directed path ends at that sink and every path in its deletion extends by that sink. This justifies the same-order count transfer without a formal cancellation or an unrealizable signed sector.

Family B's child is six greater than the next parent, so B is not closed under iteration. C and A likewise provide only the stated sparse parameter families, not arbitrary parent counts. The proof correctly distinguishes physical constructor novelty from novelty of the scalar counts. No Busch bound or tournament bank is needed for these all-k identities.

## Reproduction

    python3 check.py > NORMAL.json
    python3 -O check.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

NORMAL.json retains every core, mask, path/cycle coefficient, allowed/mandatory histogram, exact power coefficients, and all literal controls. SHA256SUMS freezes these audit artifacts and the reviewed proof.
