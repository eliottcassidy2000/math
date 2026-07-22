---
id: HYP-8800
title: LRC14 face-carry-Frobenius transfer audit and cross-problem seed map
status: >
  SYNTHESIS with new exact components; not a proof of LRC(14). THM-2041
  proves whole-packet Frobenius preservation for finite abelian character
  packets. THM-2022 now also proves rational Gamma-radial nullcones, and
  THM-346 gains prime-step tiling-walk congruences. For LRC, the missing
  theorem remains production of a nonzero safe/dual seed and its implication
  to a pointwise lonely time. The characteristic-3 period-14 propagation lemma
  remains open. THM-2043 now settles the raw characteristic-7 question:
  parity-Hasse jets are complete local coordinates, but even the full packet
  plus q-threshold is globally magnitude-blind. THM-2051 has since proved the
  Fejer--BV dissociated exit, and THM-2052 reduces every hypothetical
  counterexample to a rank-eleven bounded relation code in a rational plane.
  The live transfer target is owner-labelled rank-or-Euler termination, not a
  second familywise dissociation theorem.
source: codex-2026-07-21-NC2-transfer
related:
  - THM-2022
  - THM-2041
  - THM-2043
  - THM-2000
  - THM-1605
  - THM-2033
  - THM-1820-lrc-is-a-moment-nullcone-problem-relation-lattice-pairing
  - THM-1830
  - THM-2042
  - THM-404
  - THM-671
  - THM-2051
  - THM-2052
  - THM-873
  - THM-884
  - THM-346
  - LEM-032
  - LEM-033
  - HYP-2758
  - HYP-2974
  - HYP-2978
  - HYP-2979
  - HYP-3024
  - HYP-3036
  - HYP-3402
  - HYP-8801
  - HYP-8841
  - HYP-8846
  - OPEN-Q-108
artifacts:
  - 04-computation/lrc_frobenius_exact_period_projector_codex_20260721.py
  - 05-knowledge/results/lrc_frobenius_exact_period_projector_codex_20260721.out
  - 04-computation/frobenius_exact_period_projectors_codex_20260721.py
  - 05-knowledge/results/frobenius_exact_period_projectors_codex_20260721.out
  - 04-computation/gamma_radial_frobenius_face_codex_20260721.py
  - 05-knowledge/results/gamma_radial_frobenius_face_codex_20260721.out
  - 04-computation/lrc14_parity_hasse_jet_packet_audit_codex_20260721.py
  - 05-knowledge/results/lrc14_parity_hasse_jet_packet_audit_codex_20260721.out
  - 07-reflections/lrc14-from-local-hasse-charts-to-resolved-phase-sheaves-codex-20260721.md
reflection: 07-reflections/the-nc2-proof-transfers-to-lrc-only-after-a-frobenius-safe-certificate-is-built-codex-20260721.md
---

# HYP-8800 -- LRC14 face-carry-Frobenius transfer audit

## 1. The mechanism is four gates, not one slogan

THM-2022 succeeds because four logically separate facts align:

```text
SEED       a characteristic-zero theorem gives a whole-face Q != 0;
SELECTOR   a prime-local valuation kills every channel outside its dilation;
PRESERVER Frobenius turns the complete surviving layer into Q^p;
EXIT       a nonzero moment contradicts nullity.                         (1)
```

The reusable lesson is not merely "use Frobenius" or "use carries." A target
problem must name all four gates. This audit treats a transfer as rigorous
only when its seed and exit are as explicit as its selector and preserver.

For a compact LRC residual family, the corresponding conditional finish
schema is concrete: build integral cyclic certificates
`C_m(S)=ct_q(Theta_q Lambda_S^m)` such that a hypothetical cover forces their
vanishing; expose a base packet with `C_m0(S)!=0`; make `Theta_q` stable under
a good prime; eliminate every off-packet contribution by a proved projector
or valuation; and reconstruct a safe-residue count, endpoint gap, or named
state-lift exit from the residue. THM-2041 supplies the stability step only.

For a scalar-grade channel family

```text
M_m(c)=sum_(|r|=m, gamma dot r=0) multinomial(m;r) w(A(r))c^r,
```

the THM-2022 proof transfers if an exposed face supplies `Q!=0` and the weight
has a prime-block orientation:

```text
w(n)/w(pA0) is p-integral on the chosen side;
w(pA')/w(pA0) is divisible by p at every strict off-face A'.             (2)
```

Kummer kills nondilated allocations, (2) kills dilated off-face channels,
and Lucas/Frobenius preserve the full face. Section 8 of THM-2022 applies
this exactly to `w(A)=(alpha)_A`, proving the same nullcone and Mathieu
property for every positive rational Gamma radial shape `alpha`.

The scalar-grade hypothesis is load-bearing. THM-2022 Section 9 gives a
three-factor Wick-vector example where an off-face channel has lower
`p`-valuation than the proposed scalar-face channel. Higher-dimensional
Gaussian work needs an orthant-exposed vector face or a rank-one grade lock.

## 2. Exact LRC transfer ledger

| Gate in (1) | Best current LRC object | Exact status |
|---|---|---|
| Seed | `B5(S,q)>0` in THM-671, or a negative Toeplitz/Fejer form in HYP-2974 | Each implication to a safe phase is proved once the certificate is present; familywise supply is open |
| Selector | exact-period Ramanujan idempotents; LEM-033 conductor selects one p-adic cross-difference grade | Proved |
| Preserver | THM-2041 finite-abelian packet theorem | Proved at characteristics coprime to the ambient group order |
| Exit | positive safe count, negative Toeplitz form, endpoint-labelled current, or named K33/H7 handoff | Individual exits exist; no universal packet-to-exit theorem |

Several apparent analogues do not instantiate this table:

- **Finite phase-function products:** the function algebra on a finite cyclic
  phase set is a product of fields with many zero divisors. For a cover, the
  relevant product of safe indicators is already pointwise zero; Frobenius
  faithfully preserves that zero. There is no finite-phase DvdK seed theorem.
- **Danger-count factorial order:** the danger count lies in `{0,...,13}`.
  Its falling-factorial/Bonferroni tower terminates, so multiplying the
  intersection order by a prime greater than `13` produces zero rather than a
  new Kummer-separated level.
- LRC speeds are already integers, so THM-2022's algebraic descent gives no
  compactness or bounded-height reduction for the infinite speed family.
- THM-561/HYP-2718's factorial moments are a bounded falling-factorial basis
  transform. They do not create the `m -> pm` multinomial replication that
  Kummer uses.
- HYP-2167's `r+27k` carry is CRT address data, not a valuation penalty on all
  non-dilated channels.
- A primitive phase mask on `a in U_q` (HYP-3036) is not the primitive
  **Fourier-frequency** projector of HYP-2979. Unit dilation preserves both,
  but they cannot be identified.
- At fixed `q`, multiplying all speeds by a unit merely permutes phases:

```text
A_(pS,q)(a)=A_(S,q)(pa).                                  (3)
```

  For a zero-one safe indicator, Frobenius creates no new value. It is a
  gauge, not the new moment level that appears in THM-2022.

## 3. Whole exact-period and conductor packets really are stable

THM-2041 proves the finite-abelian theorem. If `G` is finite abelian,
`char(k)=p` does not divide `|G|`, and `e` is any descended Fourier-packet
idempotent in `k[G]`, then

```text
eX^(p^r)=(eX)^(p^r),
eX != 0  iff  (eX)^(p^r) != 0.                            (4)
```

For `G=C_N`, the exact-order-`d` idempotent is

```text
e_(N,d)=(1/N) sum_(x mod N)c_d(x)u^x.
```

For phase functions, the corresponding Ramanujan energy obeys

```text
E_d(Fr_p f)=E_d(f)^p.                                     (5)
```

The same algebraic theorem applies to parity/conductor character packets on
`G=U_q` from LEM-032/033 at Frobenius primes `p` not dividing `phi(q)`.
LEM-033 is the strongest selector discovered in the unrelated-work pass: for
a possibly different valuation prime `ell`, a conductor of `ell`-depth `beta`
sees exactly difference-valuation grade `alpha-beta` (with the proved
trivial-conductor boundary pair). It selects a complete compatible grade
packet, not a favorite atom.

There is an encoding guardrail. Group-ring Frobenius on `U_q` sends `[u]` to
`[u^p]`, while the physical LRC multiplier is usually translation `u -> pu`.
Both preserve conductor support when `p` is a unit, but they are different
operations. Any propagation proof must state which one it uses.

## 4. Why this still does not prove loneliness

The stored HYP-2979 audit supplies an exact collision. The AP row
`{1,...,13}` and the positive row `{1,...,11,13,84}` both have

```text
E_14(danger count)=216,                                   (6)
```

but their six primitive phase states are respectively `(danger,boundary)=
(0,2)` and `(1,2)`, and their weak-safe energies are `6` and `0`. A nonzero
danger packet therefore does not imply even a weak safe phase at that shell,
much less strict open mass.

The literal safe-indicator seed

```text
Q_(S,N,d)=E_d(A_(S,N))>0                                  (7)
```

would imply a weak witness, but without a restricted familywise theorem it is
just Fourier language for the witness one is trying to prove.

The non-tautological seed candidates are:

1. **THM-671 supply:** prove that every remaining speed set has some
   `q in (Vmax,2Vmax]` with `B5(S,q)>0`. THM-671 already proves the exit.
   Unit dilation preserves a found certificate by permutation, but does not
   supply the modulus or its positive sign.
2. **HYP-2974 rigidity:** prove familywise that every non-AP/GW residual has a
   finite negative Toeplitz/Fejer form, with K33/petal handoffs explicit.
3. **Owner-current seed:** combine HYP-3402 endpoint owners with a conductor
   packet and prove that its first nonzero current cannot be boundary-only.

These are stronger LRC levers than raw Ramanujan energy because each includes
an exit to pointwise safe mass.

There is also a clean modular-to-Archimedean lift once an actual count has
been built: if `D_q` counts primitive safe residues, then
`0<=D_q<=phi(q)`, so `D_q!=0 mod p` for a prime `p>phi(q)` proves `D_q>0`
over the integers. The difficult part is producing this unsigned count with
its safe inequality and endpoint labels intact; a signed trace is not enough.

## 5. New target A: repair the period-14 orbit with characteristic 3

THM-404's famous fragmentation is specific to doubling. On

```text
U_14={1,3,5,9,11,13},
```

multiplication by `3` is the single cycle

```text
1 -> 3 -> 9 -> 13 -> 11 -> 5 -> 1.                       (8)
```

Since `3` does not divide `14`, THM-2041 preserves the complete primitive
period-14 layer in characteristic `3`. Thus the old statement "even n kills
Frobenius propagation" is too broad: even `n=14` kills the **doubling**
propagation, while cubing has maximal orbit connectivity.

The seed must be nonzero as a group-algebra component **after reduction
modulo 3**. The scalar energies displayed in (6), and even the AP weak-safe
energy `6`, vanish modulo `3`; they cannot be the seed. The AP weak-safe
component itself survives modulo `3` despite its isotropic energy, which is
another reason to retain the projected vector rather than only its scalar
pairing.

The exact new theorem target is:

> **Ternary period-14 propagation target.** Construct a polynomial, Fejer,
> `B5`, or owner-current certificate whose endpoint-labelled validity is
> covariant under the cubing step around (8). Then one group-algebra seed
> packet that remains nonzero modulo `3` propagates around the entire AP
> witness orbit without isolating characters.

This target is not proved. In particular, common speed scaling and phase
permutation alone are LRC gauges. The missing content is covariance of the
certificate or the polynomial-sieve implication, not orbit arithmetic.

## 6. Resolved local target B: bad-prime Hasse jets need a resolved-phase sidecar

The good-prime theorem deliberately excludes the most structured period-14
prime. In characteristic `7`,

```text
F_7[C_14]
 = F_7[X]/((X-1)^7(X+1)^7)
 = F_7[X]/(X-1)^7  x  F_7[X]/(X+1)^7.                    (9)
```

So semisimple exact-period shells are replaced by two parity components with
seven-step nilpotent filtrations. Every polynomial has the exact local data

```text
D^(j)F(+1), D^(j)F(-1),        0<=j<=6,
D^(j)F(a)=sum_x binomial(x,j)f_x a^(x-j).                 (10)
```

These Hasse jets are a concrete `parity x seven-depth` side channel. THM-2043
proves that the fourteen coordinates in (10) are complete for a period-14
function after reduction modulo seven. It also proves the sharp limitation.
Let

```text
A   = {1,...,13},
T_n = {1,...,11,13,96+3444n}.                                      (11)
```

Every `T_n` has exactly the AP's owner-indexed residues modulo fourteen,
all raw danger/weak/boundary functions, all Hasse jets, the complete blocked
denominator mask through thirteen, and `q_threshold=14`. Nevertheless the AP
is tight while every `T_n` has the exact strict certificate

```text
C_(41,17)(T_n)
 = min_(v in T_n) (14 min(17v mod 41,41-17v mod 41)-41)
 = 1.                                                               (11a)
```

Even a fixed finite `7`-adic height truncation is insufficient: for every
`k`, one can choose a replacement congruent to `12 mod 14`, divisible by
`12`, congruent to `14 mod 41`, and with lift height `0 mod 7^k`. The same
margin-one certificate (11a) survives. Thus deeper jets cannot reconstruct
height erased before the period-14 pushforward.

The surviving theorem target is therefore not "some Hasse jet is nonzero."
It is a gluing statement

```text
endpoint-owner Hasse module
 + exact height or adaptive resolved (q,a,margin) phase
 + sign/duality implication
 -> safe phase or named boundary/state-lift exit.                    (11b)
```

For bounded heights at most `181`, the owner-labelled residue pair modulo
`13` and `14` recovers the exact speed by CRT. For the unbounded problem the
modulus must adapt, or exact height must remain. THM-2043's eleven-row
HYP-2979 audit also shows that all Hasse depths mix AP/GW, direct-q, petal,
and K33 labels in one seven-row fiber. Endpoint ownership without actual
height does not fix this.

This also corrects older repo language. Because
`Q(zeta_14)=Q(zeta_7)` has conductor `7`, the rational prime `2` is
unramified in that cyclotomic field. The doubling obstruction is
nonunit/non-etale group-algebra behavior, not field ramification.

## 7. Cross-problem transfers found by the unrelated-work pass

### 7.1 Exact transfers

1. **Gamma radial moment problems.** THM-2022 Section 8 proves the complete
   one-complex-variable nullcone and Mathieu property for every rational
   Gamma shape using the block divisibility of `(alpha)_A`. This is the
   strongest direct extension of the NC2 mechanism.
2. **Finite-abelian spectral packets.** THM-2041 preserves any descended
   exact-order, parity, or conductor idempotent at good primes. This is useful
   for cyclic signals and for the LEM-032/033 frame spectrum once the physical
   operation is identified.
3. **Tiling-cube transport.** THM-346 now records, for a mask operator
   `A_M=sum T_u` on `F_2^m`,

   ```text
   A_M^p=A_M mod p  (p odd),       A_M^2=|M|I mod 2.       (12)
   ```

   Thus aggregate `p`-step transport between arbitrary, even non-equitable,
   quotient buckets is congruent to one-step transport.
4. **Positive spanning-tree initial forms.** In the THM-856
   max-plus/deletion-contraction lane, a full minimum-tree initial form is
   automatically nonzero when weights are positive. The THM-2022 lesson is to
   retain that entire initial form instead of naming one Kruskal tree. Signed
   LRC currents lose this automatic seed and need a separate theorem.
5. **Support-measure transform discipline.** THM-2000 separates an atomic
   support from any one scalar transform. Applied here, the period-14 phase
   function is a pushforward of the owner-labelled speed support. Its kernel
   contains every lift direction `v -> v+14`; Fourier and Hasse transforms
   only change coordinates after that quotient. This gives the exact reason
   for THM-2043's no-go and points to a modulus-indexed sheaf of phase views.
6. **Orbit-product monodromy.** THM-1605 turns a local identity on one branch
   into a global contradiction by transporting it over a connected monodromy
   orbit and multiplying all conjugates. An LRC analogue would transport a
   positive resolved-phase margin over an owner packet. It requires a
   connected cover and a multiplicative invariant; neither is yet supplied.
7. **Confluent boundary jets.** THM-2033 replaces a vanishing Vandermonde at
   coincident grades by a nonzero confluent/Hasse determinant. This is an
   exact boundary-to-next-jet device, but it becomes an LRC exit only if the
   determinant retains endpoint ownership and controls a positive phase
   margin rather than merely a signed packet.

### 7.2 Transfers that stop at a guardrail

- **Binary codes/support designs:** the natural characteristic `2` is a bad
  prime for semisimple packet transfer. Odd-prime incidence packets or
  2-primary Hasse jets may help, but scalar enumerators still do not produce
  support realizability.
- **Pollock/Waring:** Freshman's dream amplifies a representation seed only
  while multiplying the number of summands by `p`. Pollock needs fixed arity,
  so the carry-pair ledger remains useful but THM-2022 does not close it.
- **Polynomial irreducibility:** THM-1735's resultant can choose good finite
  places for polynomial seeds. This resembles the prime-selection step, but
  an irreducibility or prime-value exit still needs support/fixed-divisor data.
- **Higher Gaussian dimension:** the multi-factor valuation example in
  THM-2022 Section 9 refutes a scalar-face lift. A vector-valued tropical face
  is a necessary new ingredient.
- **LRC singular series:** its exact relations survive at every prime rather
  than factor into local densities (THM-501/503). There is no Euler-product
  seed for Frobenius to amplify.
- **Rank-two Poisson/Jacobian scaffold:** the incoming THM-2042/HYP-8801
  centralizer problem does not admit p-power amplification of a canonical
  pair. In characteristic `p`, p-th powers are Poisson central, so
  `{S^p,T}=0`, not `1`; Frobenius erases rather than repairs the residual
  symplectic density. The transferable lesson is only to retain the complete
  associated-graded/Moyal correction layer during quantization. Completion
  and de-stabilization still require a characteristic-zero centralizer or
  obstruction theorem.

## 8. Ranked LRC work program

The leverage order after this audit is:

```text
1. active-owner Euler survival or rank-12 gain       (THM-2052/HYP-8841)
2. pointed transport on each rational plane          (HYP-8846)
3. THM-671 resolved-modulus/B5 supply                 (already has the exit)
4. characteristic-3 endpoint-labelled propagation    (repairs the orbit)
5. adaptive {14,27,41,...} resolved-phase sheaf       (global sidecar)
6. LEM-033 conductor-grade + owner-current glue       (selector, seed open)
7. raw Ramanujan energy / characteristic-7 jets       (diagnostic/local). (13)
```

THM-2051 has closed the former familywise Fejer/Toeplitz dissociation item:
absence of a genuine bounded higher relation already gives a strict interval.
THM-2052 starts the surviving relation code at rank eleven; the shortest gap
is now one active independent row or an owner-labelled Euler endpoint. An
unpointed safe direction in the resulting plane has the wrong quantifier.
Characteristic three remains a propagation target; characteristic seven is a
proved local chart and no-go, not a reason to add more raw jet depth.

HYP-8815 remains diagnostic only. MISTAKE-221 corrects its bounded scan to
`L_Q(S)<=M(S)`, proves the pair-sum exact engine, and withdraws its iff claims.
Any useful correlation row must now be independent of the existing relation
code and incident to the active peel/endpoint owners; a scalar score or another
copy of `1+2=3` gives no descent.

## 9. Tournament Analysis and assumption challenge

The exact companion script uses proof carriers as vertices. Its observable is

```text
(strict-LRC implication, endpoint ownership, nonzero-layer retention,
 whole-layer retention, Frobenius stability, linearity, compression).
```

It finds the transitive path

```text
endpoint-labelled primitive packet
 > primitive safe-residue count
 > primitive-period idempotent
 > Ramanujan shell energy
 > raw Ramanujan trace
 > raw denominator blockedness,
```

with score histogram `{0:1,1:1,2:1,3:1,4:1,5:1}`, no directed 3-cycles,
six singleton SCCs, and one Hamiltonian path.

Alternate vertex sets considered were runners, phases, arcs, denominators,
characters, conductor grades, Hasse jets, endpoint owners, tropical faces,
and proof obligations. The chosen quotient preserves the four gates in (1)
and destroys raw circle geometry. The challenged assumptions are:

```text
one surviving atom is needed;              false -- whole packets suffice;
any nonzero harmonic packet is a witness;  false -- (6) is a counterexample;
even n destroys every Frobenius orbit;      false -- (8) is transitive;
a good-prime projector sees the apex;       false -- the apex lives in (9);
a complete apex jet plus q-threshold exits; false -- THM-2043 family (11);
a scalar tropical face controls all Wick factors; false -- THM-2022 (20).
```

The audit's bottom line is therefore precise: THM-2022 supplies LRC with a
robust whole-packet preserver and a better search grammar. LRC still needs a
safe/dual seed theorem. The next proof attempt should be judged by whether it
fills that gate, not by whether it produces another stable arithmetic shadow.
