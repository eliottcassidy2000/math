# Classify one-sided multiplier contexts

Status: **PROVED and independently audited in the stated one-sided context class.**

This packet constructs, for every m>=1 and every actual rooted input U with exactly one Hamilton path starting at its root,

    n -> (2^m+1)n+2^m,
    n=H(P(U)),

where P(U) is the explicit edge-root-amalgamation context with count H(U)+L(U)-1. The input condition and the parent context are essential: this is not an operation proved for every represented odd n.

The all-order classification is exact within the stated one-sided context class. A core giving a(H(U)+L(U))-start(U) for every rooted U exists precisely when a=1 or a is one more than a positive power of two. The positive cores are the physical marked-edge tournaments E(a).

The factors3 and9 therefore occur, but27 cannot occur in this class, regardless of core size. A visible finite sequence of ternary-looking contexts must not be promoted to an infinite ternary induction. Mixed-mask contexts and target-dependent rebuilding remain outside this theorem.

Read PROOF.md for the construction, necessity argument, and exact scope. It uses the previously proved root-amalgamation formula and universal zero-slack classification. No eventual odd-spectrum threshold follows.

Run:

    python3 verify.py > NORMAL.json
    python3 -O verify.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

The standalone verifier reconstructs six cores and checks 138 graph pairs on 23 rooted inputs by literal subset dynamic programming. Those controls validate the implementation; they are not the all-order proof.

Independent audit: audit/REVIEW.md. The audit verifier checks 450 graph products and 1,099 inherited zero-slack controls, using no producer imports. Reviewed original and packaged proofs are retained; PROMOTION.json records the exact status-only promotion.

Frozen dependencies are in dependencies/. Their published sources, pinned at verified commit 1250ca1a59264cf3c8a3da8c98cd02d0cab80bbe, are:

- session195/zero-slack/PROOF.md, Git blob 8c977f75bdde0a97a75a04b25217103da048a561.
- session197/path-polynomials/ROOTED_ORDERED_AMALGAMATION.md, Git blob a89666a6f1ba8fd4a12999362186181f25af80af.
