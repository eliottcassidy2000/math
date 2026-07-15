---
id: THM-852
title: Klein-four structure and the n=8 refutation of the black self-line SC law
status: CORRECTED - Klein structure and normalized odd-coset lemma proved; all-n count refuted at n=8; involutive direct-realizer inference impossible by the all-n path-square law
source: kind-pasteur-2026-07-15-S128 (cont.13, overnight; owner: prove the 2selfK = SC bijection, check n=8)
depends_on: [THM-466, THM-791, THM-793, THM-849]
related: [THM-781, THM-796, THM-854, HYP-6860, HYP-6885, HYP-6890, HYP-6905]
---

# THM-852 — structure of the black self-line law

> **CORRECTION (opus-S312, 2026-07-15, THM-854 / MISTAKE-150):** the involution
> corollary in (ii) ("when Aut is trivial the unique witness is an INVOLUTION")
> is refuted by direct witness computation: at n=5 the eight witnesses have
> cycle types (1,1,3) and (5,) — orders 3 and 5; at n=6 orders 4 and 5; at
> n=7 orders 6, 10, and 12; at n=8 orders 4, 8, 10, and 12.  No direct
> realizer through n=8 is an involution.  The conditional group lemma is valid
> with the normalizer hypothesis stated below, but its square premise never
> holds for a black kappa-witness:
> σ² realizes the pair-flip p₀ Δ σp₀ ≠ ∅ (THM-854(1)). Witness orders obey the
> parity law: odd-order witnesses can occur only if n ≡ 1, 2 (mod 4)
> (THM-854(3)). Items (i),
> (iii), (iv) — the Klein-four action, orbit granularity, and the n=8
> refutation — are unaffected (my independent n=8 run and codex-S14 concur:
> 404 ≠ 176).

**The dictionary.** For `t` in the quasi-fixed locus `X`,
`T(kappa t)=P union Abar` and `T(t)^op=P^op union Abar`: the complement
tiling is "converse modulo the base path." In general
`cls(g kappa t)=cls(kappa t)^op`; on `X`, quasi-fixedness gives
`cls(kappa t)=cls(t)`, so only there may this be rewritten as
`cls(g kappa t)=cls(t)^op=cls(P^op union Abar)`. Hence quasi-fixed (self-line)
means `cls(P union A)=cls(P union Abar)`, whereas SC means
`cls(P union A)=cls(P^op union Abar)`. The proposed law `2selfK=SC` equated
the sizes of two twisted diagonals of the same Klein-four action off the
grid-symmetric locus; THM-849 refutes that equality at `n=8`.

## The conditional involution lemma

Let `A` be an odd-order subgroup of a finite group and let
`sigma in N(A) minus A` satisfy `sigma^2 in A`.  If `b` is the odd order of
`sigma^2`, then `sigma` has order `2b` and

```text
sigma^b = sigma (sigma^2)^((b-1)/2) in sigma A
```

is an involution.  This is the valid group-extension statement behind the
earlier proposed proof frame.  The normalizer and square hypotheses are
structural data; oddness of `A` does not manufacture them.  Tournament
automorphism groups have odd order, but a torsor of isomorphisms between two
different labelled tournaments is not automatically the nontrivial coset of
an index-two group extension.

For converse, relabelling equivariance supplies the appropriate twisted
extension and this lemma is useful for anti-automorphisms.  It does **not**
apply automatically to the fixed-path all-tile flip.  That operation is
indexed by the observer path `P`, and THM-854 proves the actual law

```text
sigma kappa_P = kappa_(sigma P) sigma.                      (1)
```

Consequently, for a direct realizer `sigma T=kappa_P T`,

```text
sigma^2 T
 =kappa_(sigma P) kappa_P T
 =Flip_(E(P) symmetric_difference E(sigma P)) T.           (2)
```

Thus `sigma^2 in Aut(T)` holds exactly when the two undirected path-edge sets
coincide.  It is an additional path-stabilizer condition, not a consequence
of quasi-fixedness or trivial automorphism group.  The exact THM-854 census
finds **zero** direct black realizers satisfying it at each of `n=5,6,7,8`.
It also finds zero involutive direct realizers.  In particular, the earlier
claim that every odd-size black self-line has a unique involutive witness is
false already at `n=5`.

The Burnside-affine equalities `W=|X|` at `n=5,7` prove triviality of the
automorphism groups only in those two finite censuses.  They do not prove an
all-odd-size theorem.  At `n=8`, THM-849 independently finds all 404 black
endpoints asymmetric while their unique direct realizers still fail (2)'s
square condition.  The missing information is the marked pair `(P,sigma P)`.

## Evidence log

- [x] The Klein four acts on the quasi-fixed set; every black orbit has size
  four for all `n` (THM-849 gives the fixed-point proof).
- [x] Burnside-affine weighted counts are `8,20,88` at `n=5,6,7`.
- [x] The all-size law is refuted at `n=8`: black quasi-fixed endpoints are
  `404`, not `SC(8)=176`; total quasi-fixed endpoints are `412`, not `184`.
- [x] The proposed direct-realizer involution application is impossible for
  every black kappa-witness by THM-854's path-square law; its cycle and
  path-holonomy censuses are exact at `n=5,6,7,8`.
- [ ] Decode `selfK(n)=4,6,44,202` using marked-path holonomy rather than a
  class-only bijection.  THM-849's endpoint-deletion audit rules out the
  simplest unmarked lift of the numerical size-seven residual.
