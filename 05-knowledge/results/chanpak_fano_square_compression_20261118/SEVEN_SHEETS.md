# Seven genuine sheets, even after deleting the variable block

A background research lane found this replacement rail while investigating the prior pointed-set obstruction. The root independently replayed both core profiles with two enumerators, the block-word restrictions, the first heights, and the Paley colored packings. This construction bypasses the old deletion obstruction by **changing the core**, not by cloning or editing the old rail in place. It does not cover the prime-source rail.

## Literal core pair and all-order wrapper

Let X be the transitive order 0,…,5 with arcs (0,5) and (1,3) reversed. Let P be Paley7 with x→y iff y−x is 1,2 or 4 mod 7. Write h for the Hamilton path count, and q for ordered decompositions into two **nonempty** directed paths. Then

    (h(X),q(X))=(27,234),
    (h(P),q(P))=(189,1638)=7(27,234).

The replay gives an independent literal permutation census and subset DP; in the latter q(C)=sum over nonempty proper S of h(C[S])h(C[complement S]). Both cores are strong (each has a directed Hamilton cycle).

Wrap a core C as the cyclic triangle of modules C→A→b→C, with transitive A of order a≥1 and a singleton b. The possible compressed path words are exactly

    CAb, AbC, bCA, CAbC, AbCA, CAbCA.

Indeed only A precedes b and only C precedes A in a quotient crossing; b occurs once, so tracing the forward cyclic word gives the list. The singleton does not inhibit arbitrary assignments to the ordered slots of either of the two variable modules. Counting the six words gives

    H(wrap(C,a)) = 3h+q+(h+q)(2^a−2)
                  = (h+q)2^a+h−q.

Hence the lower h_a=261·2^a−207 and the upper is exactly 7h_a. This yields 315→2205, 837→5859 and 1881→13167, and all widths beyond them. Since h_a is divisible by 9, this new ray does not contain primes.

## Free sevenfold action and actual deletion, not a scalar coincidence

Translations of the seven Paley vertices preserve the core. On either an ordered Hamilton path or an ordered two-cover the action is free, since the first vertex of the first slot moves. Choose the representatives whose first vertex is 0. The ordered representative sets have exactly 27 and 234 elements. Match the X paths and two-covers to the representatives by lexicographic rank independently, then use the translation t to label the seven sheets. In the wrapper, keep both the compressed path word and the precise allocation of A vertices. This is an explicit bijection from lower paths times Z/7 to upper paths, with Z/7 acting by actual translations of the upper core. The matching depends on the reference orders, not the symmetry action itself.

For physical deletion of any subset of A, including all of it, delete any A occurrence whose allocation becomes empty. All successive pairs of remaining module symbols must still be allowed. Validity is solely a test of this projected word. The number and order of C slots are unaltered since any two C occurrences have b between them. In particular the very same core path or two-cover is read after valid deletion; invalid restriction results agree on both sides. This proves equivariance with the actual partial deletion maps (or with an adjoined failure symbol). If A is empty, the sole valid word is bC, and the same one-slot matching applies. No fixed first-path convention is needed.

## Two honest meanings of 1–2–4 on this rail

For the translation label t, write it as alpha^t in F8 where alpha^3+alpha+1=0. Squaring rotates t→2t. The nonzero trace-zero field elements are alpha, alpha², alpha⁴ (since alpha+alpha²+alpha⁴=0); these are precisely the three Paley step-exponent labels. Choosing a line in the trace-zero plane splits the three nonzero vectors there as one plus two, and the affine trace-one coset contributes four. This is a genuine 1+2+4 chart on the sheets; the one-dimensional choice breaks the Frobenius symmetry that cyclically permutes the three trace-zero labels. This polynomial-basis model and the normal-basis model in the main note give two coordinate readings of the same F8, not the same choice of root.

There is also an intrinsic account of six on the cycle-packing side. Under translation of the binary-colored odd-cycle packings of Paley7, the fixed configurations are the empty packing and six colored constant-step Hamilton cycles (three steps 1,2,4 in two colors). Indeed a nonempty invariant packing covers a translation-invariant nonempty vertex set, hence all seven vertices. The only possibility with all its cycles of odd length at least3 is one 7-cycle, which must have a successor map commuting with translation. Every other packing has a free orbit of seven, and there are 26 such orbits. Keeping only the empty fixed configuration gives

    1 + 26·7 = 183 = 7·27−6.

Restoring the six gives 189=7·27. This is a second concrete origin of the earlier deficit of six. The fixed packings here are not fixed Hamilton paths: the [Grinberg–Stanley count](https://arxiv.org/abs/2307.05569) matches cardinality, not the translation action. In particular there is no equivariant bijection between these entire sets.

Replay: `python seven.py | diff -u seven.out -`.
