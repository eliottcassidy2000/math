# Independent phase-two audit: quantitative weighted stability

**Verdict: ACCEPT, all-order proof, with independent exact rational controls.**
This note is separate from the first-checkpoint global determinant packet.

Let A be a real weighted tournament matrix of order n>=3: diagonal zero,
0<=Aij<=1, and Aij+Aji=1. Write H for the weighted Hamilton-path sum,
u=det(I+A), v=(-1)^n det A, and P=u+|v|. For R obtained by deleting vertices
x,i, write u_R,v_R for the corresponding determinant coordinates. Then

    H-P >= max_x sum_{i!=x} Axi Aix (u_R-|v_R|)
        >= 2^(3-n) max_x sum_{i!=x} Axi Aix
        >= (2^(4-n)/n) sum_{i<j} Aij Aji.                 (1)

## Proof audit

The already-audited signed insertion recurrences express both margins
H-u+v and H-u-v as sums of nonnegative terms. For each fixed x, the old
one-vertex path (i) contributes Axi Aix (u_R-v_R) to one margin and
Axi Aix (u_R+v_R) to the other. Discarding other nonnegative terms, then
using `min(sum a_i,sum b_i)>=sum min(a_i,b_i)`, gives the first line of
(1). This argument uses the minimum of the two signed margins; it does
not select a separate global determinant sign for each deletion term.

For completeness, the needed uniform algebraic lower bound has a short
independent principal-minor proof. For an r>=1 weighted tournament put
K=A-A^T and M=I+K, and let s=1^T M^(-1)1. All principal minors of a real
skew matrix are nonnegative: odd-order ones vanish and even-order ones
are Pfaffian squares. Consequently

    det M >= 1.

Adjoin one index with skew matrix

    L = [ K    1 ]
        [-1^T  0 ].

The Schur complement gives `det(I+L)=det M(1+s)`. Expanding in principal
minors, terms not containing the new index sum to det M. All remaining
terms are nonnegative, and each two-index minor formed by the new index
and one old index equals one. Thus

    s det M = det(I+L)-det M >= r.

The determinant-coordinate identities now give

    u+v = det M / 2^(r-1) >= 2^(1-r),
    u-v = s det M / 2^(r-1) >= r 2^(1-r),
    u-|v| = min(u+v,u-v) >= 2^(1-r).

Here r=n-2 is nonempty precisely because n>=3. The final averaged bound
uses `sum_x sum_{i!=x} Axi Aix = 2 sum_{i<j} Aij Aji` and max>=average.

## Scope and hostile boundaries

At n=3 the root-variance constant is sharp: take A01=1/3, A02=A12=1.
Then H-P=2/9 and the largest root variance is 2/9. No sharpness claim for
the displayed uniform constant at other orders is made.

Order two cannot be included: A01=1/3 gives H-P=0 while the edge variance
is 2/9. The empty graph is not used as a positive-margin remainder.
These are deliberate controls, not omitted exceptional cases.

## Independent replay

`audit_weighted_stability.py` uses no producer imports. It calls the
independent Fraction Gaussian determinant in `audit_feedback2.py`, with
`tournament_kernel.py` as that file's existing local dependency, and a
separate weighted Hamilton subset DP. It checks all ternary edge matrices
through order four, 24 fixed-seed rational matrices at each order five
through eight, and the sharp order-three control: 857 matrices, 6,234
pair-deletion margins, and 3,624 root bounds. The augmented-skew determinant
identity and both cofactor lower bounds are checked directly as well.

    python3 -B audit_weighted_stability.py > weighted_stability.normal.json
    python3 -B -O audit_weighted_stability.py > weighted_stability.optimized.json
    cmp weighted_stability.normal.json weighted_stability.optimized.json

Both receipts have SHA256
`d179728dde319d9ed7b8bd200e89dc4449e772ae90ed82d17b6055fc17166ca3`.
The finite checks support the separately reviewed proof; they are not the
basis for its universal quantifier.
