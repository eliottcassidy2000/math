# Independent audit: marked whole-module context positivity

Status: **PASS after two explicit boundary-scope clarifications.** The final source is bound below. No formulas or verifier data needed correction.

## Frozen source binding

All source files are in `session197-odd-cycle-path-polynomials/`:

- `MARKED_CONTEXTS.md`: SHA256 `885af59f322f7965bcc98b6359bff5a69bb8bec33ca22cd003eb8c0a994eefea`.
- `verify_marked_contexts.py`: `54cb1815560deb053680eec5064cbc96a9706611f07b863922db54b9fe5e9d8b`.
- `marked_normal.stdout` and `marked_optimized.stdout`: `c31cfc72248e397b406749862820b5de4712890cded8cce1131117ba7eb34fb8`.
- Named dependencies: `UNIVERSAL_CONTEXTS.md` (`8615845e7bf06a6cdb1925fdbd8198a77b1a8864aa9c22704b0459678ff2ed51`), `PROOF.md` (`b88c86ac0afa64d3a1bf399c56fa7f8aad046fb93424c70a9ace8d164708183b`), `verify_extension.py` (`3b537a9cbfea66aabe8b1d4cf59c84de45cc41d2ea6f80d053df8c083b1cad76`), and `verify.py` (`10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`).

The marked theorem now states k>=1, with k=0 belonging to the inherited unmarked theorem. The new-apex sentence equating both endpoint counts to H(T) now explicitly assumes nonempty T. If T is empty, W is the directed edge b->a, so that sentence would be false; the displayed polynomial formulas themselves remain correct, giving shifted start profiles x and x+1. This boundary was checked independently.

## Mathematical review

1. Distinct named anchors in distinct labelled bins, with exactly one endpoint requirement per bin, are sufficient for the reserved-bin decomposition. A reserved path cannot contain another named anchor, because that anchor must occur in its own different bin. Choosing those actual disjoint paths leaves an ordinary bin profile on the remainder. The shift x=d-k is therefore necessary and correct, including x=0 and k=|T|.
2. Cutting whole-module runs produces exactly one distinguished local fragment for each reserved bin. Its prescribed endpoint remains the global endpoint. Other fragments are unmarked; all exterior recombination counts depend only on their number, the role pattern and the fixed context. Hence the falling-factorial linear transport acts on the full polynomial, rather than merely on one endpoint scalar.
3. In the monomial test module, end-anchors precede the free chain and start-anchors follow it. No free-chain vertex can immediately precede its end-anchor or follow its start-anchor; other anchors are unavailable. Every reserved local path is consequently its singleton anchor, giving x^m for all m>=0. This establishes the full operator on a spanning family of actual marked inputs, including its constant column.
4. In the skeleton proof, free-chain runs in distinguished open paths have no parity restriction and contribute nabla powers. Crossing odd cycles retain their separate parity projectors and factors two. Anchors and all vertices of distinguished paths contribute no ordinary-bin factor. All unused exterior vertices and packed cycles contribute the displayed X power after the shifts. Fixed labels and canonical rotations eliminate factorial ambiguity. Cutting and restoring the runs are inverse operations.
5. The identity nabla=(delta+Delta)/2 has integral nonnegative basis entries because every displayed delta/Delta coefficient, including the constant boundary cases, is even. Combining this with the inherited positive crossing-cycle kernels proves positivity on every formal basis vector, not only physical inputs.
6. The W cases are exhaustive and disjoint: neither apex in a reserved path, one apex in one reserved path, both in the same reserved path, or one in each of a start and an end path. Their multiplicities are respectively D, kX, kS and seI. Evaluation at zero yields k f(1)+se f(0). The new-apex formulas follow by choosing a possibly empty or nonempty T-path and agree with direct counts.
7. The stated limitations are appropriate: no two endpoint requirements in one bin, shared anchors, coincident bins, arbitrary root transport or arbitrary positive-matrix realizability are inferred.

## Independent computation

The source verifier was rerun normally and under optimization. Both outputs reproduce the frozen receipt byte for byte.

`check_marked_contexts.py` is a separate standalone checker. It obtains all induced path counts by literal vertex permutations, not endpoint subset DP. It then assigns each nonanchor vertex directly to a labelled bin and multiplies the literal path counts in those bins. No source counting, substitution or polynomial routines are imported. Exact rational Newton interpolation recovers the input monomial coefficients, and enough independently evaluated integer points check each complete output polynomial.

The checker reproduces all 7,714 marked context polynomials at 30,130 evaluation points, all 650 explicit W-role polynomials at 2,246 points, and 212 new-apex endpoint-profile values. The latter include the empty-input boundary separately from the nonempty endpoint corollary. Normal and optimized outputs agree.

Independent files: `check_marked_contexts.py`, `marked-independent.json`, and `marked-independent-optimized.json`, adjacent to this report. The source theorem and finite controls pass at the final proof hash above.

Independent checker SHA256: `d9446de257ba232509fc7aa0032d3c2da2dc19ddfea95fb5edcd5bda00863620`.
Independent receipt SHA256: `7a54a7ca7383dcbeeba6c62d6119b966a5a72b3f68d48d41a0efeac6cb64e0d8`.
