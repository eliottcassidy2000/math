---
id: HYP-8800
title: LRC14 needs a Frobenius-safe base certificate, not merely a Frobenius analogy
status: >
  PARTIAL / sharply localized. THM-2041 proves that primitive-period and
  Ramanujan projectors are preserved wholesale by every coprime Frobenius
  prime. Three direct NC2 transplants are ruled out: the finite phase-function
  algebra has zero divisors, the LRC danger-count alphabet is bounded, and a
  signed Ramanujan trace need not imply a safe phase. The live constructive
  target is an integral exact-period certificate with a nonzero exposed base
  packet and a positivity/pointwise lift. This is not a proof of LRC(14).
source: codex-2026-07-21-NC2-transfer
related:
  - THM-2022
  - THM-2041
  - THM-1820
  - THM-1830
  - HYP-2758
  - HYP-2978
  - HYP-2979
  - HYP-3024
  - HYP-3036
  - OPEN-Q-108
reflection: 07-reflections/the-nc2-proof-transfers-to-lrc-only-after-a-frobenius-safe-certificate-is-built-codex-20260721.md
---

# HYP-8800 -- the Frobenius-safe LRC certificate target

THM-2022 succeeds through four logically separate services:

```text
algebraic descent
  -> exposed lowest balanced face
  -> p-adic separation of every other layer
  -> Frobenius preservation of the complete tied face sum.
```

LRC already has analogues of the first two services: finite checking and
rational maximizers replace algebraic descent, while the tropical/maximin and
medium-modulus reductions expose active residue packets. THM-2041 now supplies
an exact version of the fourth service for primitive-period/Ramanujan packets.
The absent service is the third together with a pointwise output: LRC has no
Wick factorial that makes off-packet terms pay a strict valuation, and a
nonzero signed trace is not yet a lonely time.

## Conditional finish schema

For a compact LRC residual family, seek integral cyclic certificates
`C_m(S)=ct_q(Theta_q Lambda_S^m)` and infinitely many good coefficient primes
such that:

1. a hypothetical cover forces the relevant `C_m(S)` to vanish;
2. a support/period optimization exposes a base packet with `C_m0(S)!=0`;
3. `Theta_q` is stable under exponent multiplication by the good prime;
4. every off-packet contribution has strictly larger valuation or is removed
   by an exact projector; and
5. the surviving residue reconstructs an actual safe-residue count, endpoint
   gap, or positive lonely-set quantity.

THM-2041 proves item 3 for exact-period projectors and propagates item 2 once
available. Existing LRC packet work supplies much of the bookkeeping for item
4. Items 2 and 5 are the live mathematical content; neither may be replaced
by a signed scalar analogy.

## Three direct transplants that fail

1. **Pointwise danger products.** On a finite phase group, function algebra is
   a product of fields and has zero divisors. The product of complementary
   danger indicators is literally zero for a cover. Frobenius preserves that
   zero identity; there is no DvdK theorem forcing a nonzero face.
2. **Danger-count factorial moments.** The LRC count `X(t)` is at most 13, so
   `binom(X,k)=0` for `k>13`. A Kummer amplification by a prime larger than 13
   exits the alphabet instead of exposing a new moment. This is the opposite
   of NC2's unbounded moment tower.
3. **Raw Ramanujan traces.** A nonzero exact-period trace can be a signed or
   complex cancellation residue with no safe phase. It must retain the
   safe-phase inequality and open/boundary/endpoint-owner data, exactly as the
   HYP-2978/2979 and HYP-3036 quotient audits require.

## Strongest current target

Do not amplify intersection order. Amplify **period packets** after an exact
projector and an honest safe-count/endpoint certificate have isolated them.
For a primitive safe count `D_q`, `0<=D_q<=phi(q)`; therefore a congruence
`D_q != 0 mod p` with `p>phi(q)` proves `D_q>0`. This is the cleanest possible
Archimedean lift. The hard task is to express an exposed compact-residual packet
as such a count, rather than as a merely signed Ramanujan statistic.

The recommended first target is one labelled mod-14 or medium-modulus packet
where the existing exact-period projector, endpoint owner, and safe-residue
deck are all present. Prove its base count nonzero, then use THM-2041 to close
its prime-power descendants at once. This is narrower and more testable than
the earlier unrestricted TTNC proposal.
