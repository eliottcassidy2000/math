---
source: oracle-2026-06-01-S549
status: investigation + verified (doubled primes are the recursion's bridges: H1<->H2, parity, and scale)
tags: [doubled-primes, goldbach, lemoine, bridge, hyperoperation, recursion, lrc, apex, rank-one]
---

# Doubled primes are the recursion's Goldbach bridges (three bridges, one rung)

**Prompt (user):** investigate doubled primes as Goldbach bridges in relation to
this recursion (the hyperoperation tower, S548).

A doubled prime `2q` is a **bridge** in three concrete senses, all inside the
recursion, and they are the same rung seen three ways. Verified
(`doubled_primes_goldbach_bridges_s549.py`).

## (B1) The H1 <-> H2 rung: where additive-prime meets multiplicative-prime

The recursion's gear is `multiplication = repeated addition` (S548). A doubled prime
is the **unique even number where the additive and multiplicative *prime*
decompositions coincide**:

```
2q  =  q + q      (H1, additive: the GOLDBACH DIAGONAL, two equal primes)
    =  2 * q      (H2, multiplicative: a prime DOUBLING)
```

Verified: the doubled primes are *exactly* the even `n` with `n/2` prime, and these
are *exactly* the evens whose Goldbach diagonal `(n/2)+(n/2)` is a prime pair. They
are the rank-one (single-tower, S542) evens with a prime diagonal. (A non-doubled
even `n=2m` has `m` composite: no prime diagonal, and tower rank `>= 2` once `m` has
two prime factors.) So **the doubled prime is the explicit H1<->H2 rung** -- the
fixed point of "repeated addition" where the additive sum `q+q` and the
multiplicative double `2q` are the *same prime-structured* number.

## (B2) The parity bridge: Lemoine's depth-adding bridge node

Read the additive build as a small tree (S548's `H1`):

```
even n = p + q          depth-1 tree: two PRIME leaves p, q          (Goldbach)
odd  n = p + 2q         depth-2 tree: a prime leaf p + a BRIDGE NODE 2q  (Lemoine)
                                                       /    \
                                                      q      q
```

Every odd number is `prime + doubled-prime`, and the **doubled prime is the bridge
node** -- the extra level the tree needs to cross from even to odd. Even parity needs
no bridge (two leaves); odd parity needs exactly one doubled-prime bridge (a node
with two equal prime leaves). Verified Lemoine reps for `9,15,27,45`. So in the
recursion that builds the integers, **the doubled primes are the connectors that flip
parity** -- the even "carrier" every odd is built from.

## (B3) The scale bridge: the LRC apex / halving

In LRC, `n = 2q` is rank-one (`omega(n/2)=1`, S542), and its apex / source-sink
runner has speed `q = n/2` (S530/S547). The doubled prime **bridges `n` down to its
half `q`**: the apex is the halving, and the loneliest config's maximal-odd symmetry
is the `(q,q)` cycle (S547). Verified rank-one + apex `=q` for `n=6,10,14,22,26`. So
the same doubled prime that is the `H1<->H2` arithmetic rung is the `n <-> n/2` scale
rung of the LRC tower.

## The synthesis: doubled primes are the load-bearing connectors

```
                 doubled prime 2q
   H1 additive  ---- q + q (Goldbach diagonal) ----.
                                                    |  the rung
   H2 multipl.  ---- 2 * q (prime doubling) --------'
   parity       ---- the Lemoine bridge node (even -> odd)
   scale (LRC)  ---- the apex / halving n -> n/2 = q  (rank-one)
```

A doubled prime is **rank-one because it is a *single* bridge** -- one prime `q`,
doubled. The whole tower hangs on these rungs:
- pure doublings `2^k` (`n/2 = 2^{k-1}`) are the **2-adic** rank-one bridges (deeper
  towers, but `n/2` not prime -> no Goldbach prime-diagonal);
- doubled *odd* primes `2q` are the **q-adic depth-1** bridges -- the *cleanest*,
  where the additive Goldbach diagonal is a genuine prime pair.

So among all rank-one `n = 2p^k` (S542), the **doubled primes are precisely the
`k=1` (depth-1) prime bridges**, and they are the ones that carry a Goldbach
diagonal. `n = 14 = 2*7` is the canonical LRC bridge: the simplest doubled odd prime.

## Why this matters for LRC

The recursion (S548) said LRC is a cross-level claim (additive runners vs
multiplicative cascade vs exponential count). This investigation pins **where the
levels are bolted together: the doubled-prime rungs.** That is exactly why the
hardest *clean* LRC targets are the doubled primes `n=2q` (rank-one, single bridge):
the conjecture's difficulty concentrates at the rungs that join the additive and
multiplicative floors of the tower, and the doubled primes are those rungs with the
fewest moving parts -- one prime, doubled.

## Open (→ HYP)
- Does the difficulty of LRC(n) scale with how "composite" the `H1<->H2` bridge of
  `n` is -- i.e. with the rank `omega(n/2)` and whether `n` carries a Goldbach prime
  diagonal? Prediction: doubled primes (single bridge) are the cleanest/first
  provable beyond the small range; multi-bridge `n` (`omega(n/2) >= 2`) couple
  several rungs (the S532 channel coupling).
- The Lemoine depth-2 tree for odd `n`: does an odd-`n` LRC (or odd-`n` tournament
  invariant) decompose along its doubled-prime bridge node the way even-`n` splits
  along its Goldbach edge?

## Anchor
`04-computation/doubled_primes_goldbach_bridges_s549.py` (+ `.out`): doubled primes =
evens with `n/2` prime = Goldbach-prime-diagonal = rank-one depth-1 bridges (B1);
Lemoine bridge node (B2); LRC apex/halving (B3). Builds on S546 (doubling=pairing),
S547 (apex/(q,q)), S548 (tower), S542 (rank).
