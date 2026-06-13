---
id: HYP-2049
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S549
related:
  - HYP-2048
  - HYP-2044
  - HYP-2035
---

# HYP-2049: doubled primes are the recursion's bridges (H1<->H2, parity, scale); LRC difficulty scales with bridge compositeness

**VERIFIED (`doubled_primes_goldbach_bridges_s549.py`):** a doubled prime 2q is a
bridge in three concrete senses inside the hyperoperation recursion (S548):
- (B1) H1<->H2 rung: 2q = q+q (additive Goldbach DIAGONAL, two equal primes) = 2*q
  (multiplicative prime doubling). Doubled primes = exactly the evens with n/2 PRIME
  = exactly the evens with a Goldbach prime-diagonal = the rank-one depth-1 cases.
- (B2) PARITY bridge (Lemoine): even=p+q (depth-1 tree, two prime leaves); odd=p+2q
  (depth-2 tree: prime leaf p + BRIDGE NODE 2q with leaves q,q). The doubled prime is
  the extra level that crosses even->odd; every odd needs one.
- (B3) SCALE bridge (LRC): n=2q rank-one (omega(n/2)=1, S542), apex = q = n/2 (S547);
  bridges n down to its half q. Verified n=6,10,14,22,26.

Among all rank-one n=2p^k, the doubled primes are the k=1 depth-1 PRIME bridges (the
only ones carrying a Goldbach prime-diagonal); pure doublings 2^k are 2-adic rank-one
but n/2 not prime. n=14=2*7 is the canonical LRC bridge.

**CLAIM:** the doubled primes are where the additive (H1) and multiplicative (H2)
floors of the tower are bolted together -- LRC's cross-level difficulty (S548)
concentrates at these rungs, and doubled primes are the rungs with the fewest moving
parts (one prime, doubled). So:
- (A) LRC difficulty scales with bridge compositeness: doubled primes (single bridge)
  are the cleanest/first provable beyond the small range; multi-bridge n
  (omega(n/2)>=2) couple several rungs (the S532 channel coupling).
- (B) odd-n structures (LRC or tournament invariants) decompose along the Lemoine
  doubled-prime bridge node the way even-n splits along its Goldbach edge.

**Files:** `04-computation/doubled_primes_goldbach_bridges_s549.py` (+.out). Reflection:
`07-reflections/doubled-primes-are-the-recursions-goldbach-bridges-s549.md`.
