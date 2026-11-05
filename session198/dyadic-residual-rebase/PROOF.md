# An exact dyadic balanced fiber and a same-order residual replacement

**PROVED.** This is a restricted positive
replacement, not a cofinite or universal odd-coverage theorem.

Use the physical source domain, response `(D,E,U,d)`, and chronological moves
O/Z/J from the published balanced-source normal form, commit
`cf2263bb4449776d1ce13852528d5d55141c9252`.

## 1. Dyadic weak fibers have balanced cores

Every weak source is a strong start-one core, or singleton, followed by a
nonempty transitive sequence of universal sinks. For a strong source,
`U<=2(D-E)<2D`; this follows from the positive grammar: when `E>0`, a last J
has `f=0`, and a last Z has its new `a>=f`, so `f<=a`. The singleton has `U=D`.
Thus every allowed initial core has `1<=U/D<2`.

It follows that a source satisfying `U=2^k D`, `k>=1`, is **exactly**
`O^k(T)` for a balanced source T or a singleton. It cannot be strong, and if
its tail has length r, its core has ratio `2^(k-r)` in `[1,2)`. Consequently
r=k and the core is balanced. The core and the sink tail are canonical from
the strongly connected component decomposition.

This proves an iff characterization of the entire physical fiber, not just
an identity on a selected word family.

## 2. The corresponding balanced fiber

Let W be a balanced source with response `(h,e,h,e)`, and let `k>=1`. Then

    h=(2^(k+1)+1)e

holds iff, canonically,

    W = J(O^(k+1)(T))

for a balanced source T, or singleton, of Hamilton-path count e.

Indeed the canonical inverse of W has profile `(D,E,U,d)` with
`D=e` and `U=(h-e)/2=2^k e`. Section1 gives the result. Conversely the actual
construction has the stated response. The singleton core is allowed even
though W itself is always nontrivial.

## 3. Same-order replacement of the full indicated residual fiber

For every such T and every k>=1, compare the two chronological constructions

    A: T, O^(k+1), J, O, Z^k, J;
    B: T, O^(2k+3), J.

Both are actual strong C1 graphs of order `|T|+2k+4`, and both have

    H=(2^(2k+3)+1) H(T).

Graph A has no balanced root, while B is balanced at its last vertex. Hence
they are not isomorphic. This is a replacement preserving order, H and C1,
not a claimed all-edge isomorphism.

To verify the scalar identity, put `b=2^(k+1)` and `h_0=H(T)`. The intermediate
W in A is balanced with `(h,e)=((b+1)h_0,h_0)`. The exact residual readout from
the published graph dichotomy is

    R_k(W)=(2b-1)h-(b-2)e=(2b^2+1)h_0.

Construction B has the same readout by the elementary O/J law. Its vertex
count agrees because both words add `2k+4` vertices. The graph dichotomy
excludes every balanced root of A, since k>=1 and W is nontrivial.

By Section2, this replacement applies to every physical balanced W with the
specified ratio, not just to graphs initially supplied in the displayed form.
For T a singleton and k=1, the two seven-vertex C1 graphs both have H=33;
for k=2, the two nine-vertex graphs have H=129.

The source endpoint is not freely chosen: membership in the indicated
balanced ratio fiber is an explicit necessary hypothesis. The theorem does
not handle an arbitrary residual W and has no consequence for all large odd
integers by itself.

## 4. Verification and scope

The standalone verifier checks the exact dyadic profile fibers against the
complete physical O/Z/J bank through U=2000; it also checks canonical graph
recovery on every labelled rooted source through order five. Both graph
constructions are counted directly for fifteen pairs, through order fifteen.
Finally their symbolic balanced-input rows are compared for k=1 through40.
The all-order proof is the argument above, not extrapolation of those tests.

Run `python3 verify.py` and `python3 -O verify.py`; the supplied receipts agree.
The dependency statements are frozen in `dependencies/` with pinned provenance
in `DEPENDENCIES.md`. Universal odd coverage and the balanced cofinite lemma
remain open.
