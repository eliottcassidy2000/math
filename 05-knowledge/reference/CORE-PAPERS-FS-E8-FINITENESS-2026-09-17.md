# Primary-source audit: Fargues-Scholze finite fibers for E8

**Audit date:** 2026-09-17. **Status:** CITED / CONDITIONAL literature audit.
This records the versions actually retrieved, not an assertion that every
possible manuscript has been located. The target is the map on isomorphism
classes of smooth irreducible representations over `Qbar_l`, with `l != p`,
to conjugacy classes of semisimple Weil L-parameters. A theorem about a
different map, or one assuming finite fibers, cannot settle that target.

## Cotner-Feng: an exact conditional theorem in the requested range

Sean Cotner and Tony Feng, *Local Langlands functoriality for Yu's
supercuspidals II: Fargues-Scholze's parametrization*,
[arXiv:2609.16387v1](https://arxiv.org/abs/2609.16387v1), submitted
2026-09-14; [primary PDF](https://arxiv.org/pdf/2609.16387v1).

**CONDITIONAL, Theorem 10.2.1, printed pp. 48-49:** Let `F` be a
nonarchimedean local field and `G/F` connected reductive, splitting over a
tame extension, with `p` not dividing its absolute Weyl-group order.
Assuming Hypothesis 4.1.4 for cyclic base change at every auxiliary prime
different from `p`, the Fargues-Scholze map preserves depth and has finite
fibers. Its proof treats coefficients `Qbar_l` and `Fbar_l`.

Hypothesis 4.1.4 requires compatibility of parameters with the sigma-dual
homomorphism for irreducible constituents of Tate cohomology of
sigma-invariant representations, including the coefficient Frobenius twist.
Here only its cyclic-base-change instances are assumed.

The introduction and Remark 4.1.2 explicitly retain this hypothesis:
available modular functoriality and Wei's parity result cover auxiliary
primes greater than five; the all-prime extension is described as
forthcoming. Thus this retrieved theorem is not an unconditional solution
for E8. The condition `p>7` concerns residue characteristic, and does not
remove the auxiliary primes `2,3,5` from its assumption.

## Feng: retrieved revision still has a small-prime restriction

Tony Feng, *Modular functoriality in the Local Langlands Correspondence*,
[arXiv:2312.12542v3](https://arxiv.org/abs/2312.12542v3), revised 2024-08-26;
[author-hosted PDF](https://math.berkeley.edu/~fengt/FS-functoriality.pdf).

**CITED:** Theorem 1.3.1 assumes the auxiliary prime exceeds the root-system
bounds for both dual groups; Figure 1 gives the E8 bound as 31. This is
a condition on coefficient characteristic and automorphism order, not on
the residue characteristic `p`. Checking the
[author's current papers page](https://math.berkeley.edu/~fengt/papers.html),
the linked PDF, and arXiv metadata did not locate the forthcoming
all-prime revision on the audit date. This is a bounded retrieval result,
not a proof that no other version exists.

## DHKM: the unconditional bounded-depth input

Jean-Francois Dat, David Helm, Robert Kurinczuk, Gilbert Moss,
*Finiteness for Hecke algebras of p-adic groups*, JAMS 37 (2024), 929-949;
[arXiv:2203.04929v2](https://arxiv.org/abs/2203.04929v2), revised 2022-04-22;
[author-hosted PDF](https://webusers.imj-prg.fr/~jean-francois.dat/recherche/publis/finiteness.pdf).

**CITED, Corollary 3.5(2), printed p. 11:** For each depth bound `r>0`,
there is a wild-inertia cutoff `e` such that the excursion action on the
bounded-depth Bernstein center factors through

```text
Exc(W_F^0/P_F^e, Ghat)_red,
```

and the center over `Z_l[sqrt(q)]` is finite as a module over this algebra
after the corresponding coefficient extension. Corollary 3.5(1) gives
admissibility over the excursion algebra for every finitely generated
smooth `Z_l[sqrt(q)]G(F)`-module. Here `l != p`.

Theorem 1.1 gives finite generation of each compact-open Hecke algebra over
its center for any noetherian `Z_l`-algebra. Together these support
finiteness of FS fibers inside each bounded-depth category. Passing from
all bounded-depth pieces to a whole fiber still needs a depth bound for
that fiber; a union of finite sets need not be finite.

**Different finite map:** Corollary 2.4 concerns pushforward between
parameter GIT quotients for a Weil-stable closed reductive subgroup of
`Ghat`, at a fixed wild cutoff. It must not be identified with the FS map
from all irreducible representations.

## DHKM banal case: avoid a circular use of quasi-finiteness

Dat-Helm-Kurinczuk-Moss, *Local Langlands in families: The banal case*,
[arXiv:2406.09283v2](https://arxiv.org/abs/2406.09283v2), revised 2024-09-23;
[primary PDF](https://arxiv.org/pdf/2406.09283v2).

**CITED boundary:** Theorem 8.2 constructs a quasi-finite morphism between
the spectral and Bernstein coordinate rings from a supplied semisimple
correspondence. Definition 6.19, condition `(C1)`, already requires that
correspondence to have finite fibers. Consequently this theorem cannot
establish finite fibers for the E8 Fargues-Scholze map unless that input
condition has first been justified independently. Its classical-group
applications do not automatically supply the exceptional-group input.

## The E8 prime check

Jessica Fintzen, *Types for tame p-adic groups*, Annals of Mathematics 193
(2021), 303-346; [arXiv:1810.04198v2](https://arxiv.org/abs/1810.04198v2),
revised 2020-11-03; [primary PDF, Table 1, printed p. 8](https://arxiv.org/pdf/1810.04198v2).

**CITED / exact arithmetic:** Table 1 gives

```text
|W(E8)| = 2^14 * 3^5 * 5^2 * 7 = 696729600.
```

Independently, multiplying the standard invariant degrees
`2,8,12,14,18,20,24,30` gives the same integer. Hence for prime `p`,
`p>7` is exactly `p` not dividing this order. Split E8 satisfies the tame
splitting requirement. The table separately lists E8's bad primes as
`2,3,5`, so “good prime” and “prime not dividing the Weyl-group order”
must also be kept distinct.

## Research consequence and scope

The retrieved sources give an unconditional finite-fiber statement after
imposing a depth bound, and a conditional theorem giving the desired
unrestricted conclusion in the E8 range. They do not, as cited here,
supply the remaining unconditional bound on representation depths inside
an arbitrary fixed FS fiber. That is the precise unresolved obligation
for this research attempt; it is weaker than proving exact depth
preservation. No theorem is promoted merely because its missing
hypothesis is expected in forthcoming work.
