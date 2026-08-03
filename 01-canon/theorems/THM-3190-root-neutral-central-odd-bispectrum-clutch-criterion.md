---
id: THM-3190
title: "Root-neutral central-odd bispectrum and clutch criterion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  THM-2312's nonzero cubic word-current vector has exact bidegree
  (root-translation character zero, common-scalar charge one).  It is
  therefore invariant under cyclic root reindexing but changes sign under a
  common scalar -1.  On a common equivariant horn carrier where THM-2889's
  quaternionic center acts by that scalar, the full 132-coordinate
  bispectrum vector would detect the central clutch after the known torus
  action is removed.  That common carrier, physical e9 edge, and torus
  normalization are not constructed here.
source: root/2026-08-02
audit: >
  An independent hostile audit rederived the root exponent and 66/66 carry
  split, the common U(1) charge lambda^2 conjugate(lambda)=lambda, the
  THM-2312 nonzero 132-vector handoff, and the conditional -U_tau clutch
  identity.  It verified that common-scalar oddness is not called a literal
  THM-2889 central action before the equivariant-carrier hypotheses.  The
  integer/exponent companion pins both dependencies and checks 1,716 root
  translations, 495 scalar phases, and 3,432 commuting bidegrees.  Normal,
  optimized, and stored outputs agree and the LF hashes match.
depends_on:
  - THM-2312-sparse-root-bispectrum-positive-word-current
  - THM-3187-central-sign-parity-quotient-and-odd-observable-necessity
related:
  - THM-2889-dicyclic-reverse-action-joint-carrier-and-skew-lift-separation
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2420-affine-shell-cross-reference-composition-and-complete-zero-reference-hostile
script: 04-computation/lrc14_root_neutral_central_odd_bispectrum_thm3190.py
output: 05-knowledge/results/lrc14_root_neutral_central_odd_bispectrum_thm3190.out
script_sha256: 4309bb71cbe5fd888f13a178f097826e9c0a8bf99edb566347d58b3c93840fa6
output_sha256: 4ee43881437db8ce70e0571cd0c8cd960d94b21a5dd8d19197a0b45210e9f0a2
hash_basis: LF-normalized bytes
---

# THM-3190 -- root-neutral central-odd bispectrum and clutch criterion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3187 proves that an odd observable is necessary to see THM-2889's
quaternionic central sign.  A seemingly unrelated earlier construction
already supplies the correct parity: THM-2312's sparse-root bispectrum is
cubic and nonzero on every selected positive blocker word.  Its character
indices make it neutral under cyclic root translation, while its unequal
holomorphic/antiholomorphic degree leaves one unit of common scalar phase.

This is the exact combination needed by a central-clutch detector:

```text
root address reindexing:       invisible;
common scalar phase:           retained;
central scalar -1:             sign reversal.                (1)
```

The result is a proved bidegree and a precise conditional bridge.  The
load-bearing common-carrier map to the THM-2889 horn is still open.

## 1. Two commuting actions on the Fourier amplitudes

Let `p=13`, let `zeta` be a primitive thirteenth root, and use THM-2312's
Fourier amplitudes

```text
M_k=sum_(r in F_13) v_r zeta^(-kr).                         (2)
```

For `k,l,k+l!=0 (mod 13)`, put `c=(k+l) mod 13` and

```text
B_(k,l)=M_k M_l conjugate(M_c).                             (3)
```

There are `12*11=132` ordered pairs.  A cyclic root translation by `h`
acts as

```text
M_j -> zeta^(-jh)M_j.                                      (4)
```

Therefore

```text
B_(k,l) -> zeta^((-k-l+c)h)B_(k,l)=B_(k,l),                (5)
```

because `k+l-c` is either `0` or `13`.  Exactly `66` pairs have each carry.

Independently, a common scalar phase `lambda in U(1)` acts by

```text
M_j -> lambda M_j,
B_(k,l) -> lambda^2 conjugate(lambda) B_(k,l)
          =lambda B_(k,l).                                 (6)
```

Thus the bispectrum has the commuting bidegree

```text
(C_13 root character, common U(1) scalar charge)=(0,1).     (7)
```

At the order-two common scalar phase `lambda=-1`, every coordinate changes
sign.  This becomes the THM-2889 central action only under the equivariant
identification hypothesized in Section 3.

### Proof

Equations `(5)` and `(6)` are the direct character and scalar-degree counts.
The two actions commute because root translation is diagonal in the Fourier
index and the common scalar is central.  QED.

The phrase *root-neutral* is load-bearing.  The cubic is not invariant under
an arbitrary common `U(1)` phase; it has charge one under that action.  Calling
it simply gauge-invariant would erase the very coordinate which detects the
central sign.

## 2. A proved nonzero root-neutral, common-scalar-odd word current

Fix any positive THM-2312/THM-2305 canonical word `Q`.  Integrate `(3)` on
that word:

```text
C_(k,l)(Q)=integral_Q B_(k,l)(y)dy.                         (8)
```

THM-2312 proves the quantitative whole-face inequality

```text
sum_(k,l allowed) C_(k,l)(Q)>0,                             (9)
```

and hence the complete vector

```text
C(Q)=(C_(k,l)(Q))_(k,l allowed) in C^132                  (10)
```

is nonzero.  Equations `(5)--(7)` hold pointwise and survive integration.
Consequently `C(Q)` is an already-proved nonzero word current which is
root-translation neutral and odd under a common scalar `-1`.  It has the
formal parity required of a central-odd current; it is not literally a
THM-2889 central representation until the Section 3 hypotheses are supplied.

This does not contradict THM-2312's phase boundary.  On a one-sheet fibre
the bispectrum is blind to the root address precisely because of `(5)`, but
it remains linear in a common scalar phase by `(6)`.  Root address and global
central phase are different coordinates.

## 3. Exact common-carrier clutch criterion

Let `V` be THM-2889's four-channel abstract horn representation.  Write

```text
z=(0,0,-1),        tau=(13,13,+1),
kappa=z tau=(13,13,-1).                                   (11)
```

Suppose a physical horn packet supplies an assignment

```text
Psi: V_endpoint -> (M_j(y))_(j in F_13)                    (12)
```

on one fixed positive word, satisfying all of the following.

1. `Psi` is defined on the same direct and via packet and is equivariant for
   the relevant horn actions.
2. The central element `z` acts on every amplitude by the common scalar
   `-1`.
3. The torus element `tau` induces a known invertible linear action `U_tau`
   on the full `132`-coordinate cubic vector.
4. The physical right edge `e9=(-9,+9,QB)` exists on that packet, so that
   THM-2889's already-composed group identity gives

   ```text
   via=kappa*(direct).                                      (13)
   ```

Then `(6)` implies the exact vector identity

```text
C_via=-U_tau C_direct.                                     (14)
```

Since `C_direct!=0` by `(9)`, the torus-normalized comparison is nonzero and
detects the central sign:

```text
U_tau^(-1)C_via=-C_direct!=C_direct.                        (15)
```

The last inequality uses characteristic zero and `C_direct!=0`.

Keeping the full vector is useful: no uniform choice of the nonzero pair
`(k,l)` is required, and torus mixing among the 132 coordinates is allowed.
Only invertibility and a common selector are required.

### What is and is not proved in the criterion

The implication `(12)--(15)` is proved.  Its hypotheses are not silently
asserted.  THM-2312 and THM-2889 currently live on differently typed
carriers.  Canon does not provide:

- an equivariant identification `(12)`;
- the physical `e9` on the same marked word packet;
- the induced torus action `U_tau` on that packet; or
- a proof that the THM-2312 positive word is preserved by the horn path.

The bispectrum is therefore an odd direct observable candidate, not the
non-co-shifted charged reference of THM-3187/THM-2420.  Root neutrality alone
does not create the missing common carrier.

## 4. Connection contract and cheapest decisive test

```text
source:
  THM-2312's nonzero 132-coordinate cubic current on one exact positive
  owner/clock/word;

target:
  THM-2889's torus-normalized central clutch;

map needed:
  a common-packet equivariant amplitude assignment Psi as in (12);

preserved if the map exists:
  complete word, root-translation quotient, the full cubic vector, and the
  odd common-scalar phase;

destroyed already:
  absolute root address, a uniform selected character pair, and the linear
  ordinary Fourier amplitude;

cheapest test:
  construct e9 on one THM-2312-positive marked word and check whether its
  132-vector transforms by -U_tau rather than +U_tau.                 (16)
```

Failure of `(12)` would be informative: it would identify carrier alignment,
not parity or positivity, as the first failed implication.

## 5. Exact evidence and scope

Run

```text
python 04-computation/lrc14_root_neutral_central_odd_bispectrum_thm3190.py
python -O 04-computation/lrc14_root_neutral_central_odd_bispectrum_thm3190.py
```

and compare LF-normalized bytes with the declared output.  The companion
pins the exact THM-2312 and THM-3187 scripts/transcripts, enumerates all `132`
allowed pairs and the `66/66` carry split, and verifies all root-translation,
global-phase, and commuting-action exponents with integer arithmetic only.

The theorem supplies a nonzero observable of the correct two grades and an
exact conditional clutch detector.  It does not supply the common carrier,
torus normalization, physical right edge, semantic endpoint transport, row
exclusion, or LRC(14).  It makes no NC2/GMC2 implication.

**QED.**
