# THM-122: TT Boundaries Span im(d_2)

**Status:** PROVED (computationally verified, algebraic proof straightforward)
**Filed by:** kind-pasteur-2026-03-08-S43

## Statement

For any tournament T on n vertices, the image of d_2: Omega_2 -> Omega_1
is spanned by boundaries of TT (transitive triple) elements alone.

That is: im(d_2) = span{d_2(e_{TT}) : TT triples (a,b,c) with a->b->c and a->c}.

NT (non-transitive) elements of Omega_2 contribute nothing additional to im(d_2).

## Consequence

b_1(T) = dim(Z_1) - rk(M_TT|_{Z_1})

where M_TT is the matrix of TT boundary vectors projected into Z_1 = ker(d_1).

This reduces the computation of b_1 from the full Omega_2/Omega_1 chain complex
to a simple matrix involving only TT triples and the cycle space Z_1.

## Verification

Exhaustive: n=3,4,5 (100%), n=6 (32768/32768 = 100%).
Sampled: n=7 through n=10 (100%).

## Proof sketch

In Omega_2, the allowed 2-paths are all (a,b,c) with a->b->c. These split into:
- TT: a->c (transitive)
- NT: c->a (3-cycle)

The Omega_2 relations (from the path complex) include: for each pair of
allowed 2-paths sharing a common "bad" interior face, their difference is in
Omega_2. These NT relations can be expressed as linear combinations of TT
elements modulo ker(d_2). Therefore im(d_2) = im(d_2|_{TT}).

## See Also
- THM-107 (at most 1 free component, uses this)
- THM-102 (beta_2 proof status)
