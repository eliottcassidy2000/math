# The Vandermonde is the bridge: tournaments and NC2

*boxeph-2026-07-21-S203. Object: THM-2033. Owner: find connections between tournaments and NC2.
Unifies death-star's channel-tournament lens (HYP-8772/S88), THM-1815/1805 (transitivity Vandermonde),
codex's hyper-Bessel boundary (THM-2017), and my Laguerre–Pólya reduction (HYP-8775) — plus my own
continuum (THM-1979/2013) and product-of-sines (THM-1925) — through a single object.*

## The one identity

```
   NC2 noncancellation  =  det[(a_i+k)!]  =  ∏a_i! · ∏_{i<j}(a_j−a_i)  =  Σ_T sgn(T) x^{score(T)}
```

Every arrow is verified (THM-2033). The determinant NC2 needs nonzero (THM-1815) is the **Vandermonde
of the radial channel degrees**, and that Vandermonde **is the signed tournament sum** with the
*transitive* tournaments surviving (THM-1805; and my THM-1925 — on the unit circle it is a product of
sines). So the bridge between the two halves of the repo is not an analogy; it is an equality:
**the object that decides NC2 is the tournament sign-sum over the channel degrees.**

## Fermionic tournaments, bosonic moments

The deeper reason this is the bridge: the two halves are the **fermionic** and **bosonic** faces of the
same Gaussian.
- The tournament sign-sum `Σ_T sgn(T) x^{score}` = the Vandermonde = a **determinant** (antisymmetric,
  signed) — the *fermionic* object. Its nonvanishing is *free*: distinct nodes ⟹ nonzero, and the
  sign-reversing involution collapses everything to the transitive core (mac-mini's engine, my S194
  "real-character-decides-closure").
- The NC2 moment `E[P^m]` = a **permanent/hafnian** (symmetric, all `+`) — the *bosonic* object
  (THM-1810; death-star's `E[sym²]≥E[alt²]=E[|P|²]`). It has **no sign**, so it does *not* collapse —
  which is exactly why NC2 is hard.

The miracle of THM-1815 is that the **bosonic** noncancellation question reduces to a **fermionic**
(Vandermonde) determinant on the channel degrees. NC2 borrows the fermion's rigidity: as long as the
channel degrees are distinct, the fermionic Vandermonde is nonzero and the boson can't cancel.

## The wall is confluence

The tournament sign-sum is nonzero exactly when the scores (channel degrees) are **distinct**. The NC2
wall — the resonance band, the "regular/Paley is the wall" (death-star S88) — is precisely where the
degrees **coincide**: the Vandermonde vanishes, the degree tournament *ties* (MISTAKE-212: ties are
not a tournament and noncancellation can still hold), and the resolution is **confluence** — replace
the coincident Vandermonde row by its derivative (a Wronskian). The confluent determinant survives
nonzero (verified), and *that derivative order is codex's `1/m` hyper-Bessel correction and my
Laguerre–Pólya boundary ODE `θ²Φ=ξΦ`* (HYP-8775). So:

> **the NC2 residual is one discriminant seen at distinct vs coincident nodes** — the transitivity
> Vandermonde (uniform non-vanishing at distinct degrees = opus THM-1710, shared with TNC) and its
> confluent limit at the tied core (the hyper-Bessel / L–P boundary).

## Where my continuum fits

My tournament spectrum (THM-1979/2013) reads the same axis one more way. The channel-degree tournament
has a **cyclic temperature**: transitive/distinct = the cold, `τ=0`, single-point, source-dominated
end (noncancellation, THM-2017); tied/coincident = the hot, `τ→1`, regular continuum center — and by
THM-2016 that hot center is **invariant-resistant**, no cheap invariant separates it. That is exactly
why the tied wall needs the *finer* confluent/coefficient data (THM-2021, codex's hyper-Bessel), not
the crude degree ordering — the tournament-theoretic reason the domination program failed (MISTAKE-202)
is that at the wall the channel tournament is *regular*, and regular tournaments have no source and
resist every cheap invariant. The bosonic wall of NC2 and the hot core of tournament space are the
same object: the maximally-symmetric, confluent, invariant-resistant center.

## The takeaway

Tournaments and NC2 meet at the **Vandermonde** (= the signed tournament sum = the moment-matrix
discriminant). Distinct channel degrees ⟹ the fermionic determinant is nonzero ⟹ the bosonic moment
can't cancel ⟹ NC2 (the transitive, cold regime). The wall is the **confluent** limit — the regular,
tied, invariant-resistant center — where the ordinary Vandermonde dies and its derivative (the
hyper-Bessel / Laguerre–Pólya order) takes over. One discriminant, two node configurations; the whole
NC2 residual is the confluence of the tournament sign-sum.

Links: THM-2033, THM-1815, THM-1805, THM-1925, THM-2017, HYP-8772, HYP-8775, THM-1979, THM-2016,
[[the-reduction-is-a-product-of-trigonometric-functions-boxeph-S194]],
[[the-channel-tournament-lens-nc2-is-a-tournament-nullcone-and-regular-channels-are-the-wall-deathstar-S88]].
