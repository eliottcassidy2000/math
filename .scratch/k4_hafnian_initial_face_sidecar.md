# K4 hafnian initial-face sidecar and the unbounded cancellation boundary

**SCRATCH PROOF-COMPLETE THEOREM / NO CANON ID / AWAITING REVIEW.**

## Inheritance and question

The closest proved mechanism is THM-3049, which tropicalizes the three
perfect-matching monomials of `K4` and finds the binary/ternary integral
clutch.  Its sharp hostile is

```text
(A_1,A_2,A_3)=(1,1,-2):        v(A_i)=(0,0,0),
A_1+A_2+A_3=0.
```

THM-2290 says why this is the right level: a fixed endpoint-colour query is
the hafnian of its selected pair kernel, but selector or phase erasure is not
a congruence.  THM-2598 supplies the exact `S4/V4=S3` matching torsor and also
warns that the three-matchings quotient does not choose a quartic origin.

The question is what least extra datum, beyond matching valuations, certifies
that the contracted three-term hafnian does not vanish.

## The sharp DVR statement

Let `K` be a discretely valued field, with valuation `v:K^*->Z`, uniformizer
`pi`, valuation ring `O`, and residue field `k` (finite or infinite).  Let

```text
H=A_1+A_2+A_3,                 A_i in K^*.
```

Put

```text
lambda_i=v(A_i),
m=min_i lambda_i,
F={i:lambda_i=m},
alpha_i=bar(pi^(-lambda_i) A_i) in k^*,
sigma=sum_(i in F) alpha_i in k.                         (1)
```

Then:

1. `sigma!=0` if and only if `v(H)=m`.  In particular it certifies `H!=0`.
2. Fix the full labelled data `(lambda_i,alpha_i)`.  Every lift of these data
   has nonzero sum if and only if `sigma!=0`.
3. If `sigma=0`, there are two lifts with exactly the same
   `(lambda_i,alpha_i)`, one with sum zero and one with nonzero sum of larger
   valuation.  Thus valuations plus leading residue units do **not** decide
   actual nonvanishing on the zero-initial-form locus.

The proof is elementary and exact.  Divide by `pi^m`.  Reduction modulo the
maximal ideal gives precisely `sigma`, proving (1).  For (3), choose
`j in F`, lift all terms except `A_j` arbitrarily, and define

```text
A_j=-sum_(i!=j) A_i.
```

Because `sigma=0`, this constructed `A_j` has the prescribed valuation and
leading residue.  It gives an exact zero.  Replacing it by
`A_j+pi^(m+1)` preserves the same first-order data and makes the sum
`pi^(m+1)`.

The zero/nonzero status of `sigma` is independent of the chosen uniformizer.
Replacing `pi` by `u*pi` multiplies every `alpha_i` on the minimal face by
the same nonzero residue `bar(u)^(-m)`.  Thus the natural decision bit is
intrinsic even though a displayed scalar representative uses `pi`.

Consequently, **given the full tropical valuation vector**, the natural
minimal residue sidecar is the one scalar `sigma`; its zero/nonzero bit is
the exact decision quotient for nonvanishing uniformly over all lifts.  It is
not an iff test for one already-fixed amplitude.

## The exact no-go beyond first order

For actual nonvanishing, the canonical object is the filtered contraction
jet.  Starting with `pi^(-m)H in O`, take its successive classes in

```text
O/(pi),  (pi)/(pi^2),  (pi^2)/(pi^3), ... .            (2)
```

The first class is `sigma`.  The sequence reaches a nonzero class exactly
when `H!=0`; it vanishes forever exactly when `H=0`.  No bounded initial
segment decides this uniformly.  For every `N>=1` (in characteristic not
two),

```text
(1,1,-2)                 has sum 0,
(1+pi^N,1,-2)            has sum pi^N,                 (3)
```

while the two triples agree through every coefficient of depth `<N`.
The same construction in the proof above works in arbitrary residue
characteristic.  Finiteness of `k` does not shorten this tower.

Thus the answer to the motivating question has two layers:

```text
robust/noncancelling leading term:  exactly sigma!=0;
mere actual nonvanishing at sigma=0: unbounded higher contraction jet needed.
```

## Why the 24-valued clutch is still too coarse

THM-3049's clutch keeps only

```text
(lambda mod 2, sum lambda mod 3).
```

It does not determine the minimal face `F`.  Over `Q_5`, the two matching
triples

```text
(1,5^2,-2*5^4): lambda=(0,2,4),   sum!=0,
(5^2,5^2,-2*5^2): lambda=(2,2,2), sum=0                 (4)
```

have the same clutch `(000,0)` and the same labelled leading residues
`(1,1,3)`.  The first has a unique minimal term and is robustly nonzero; the
second cancels.  Both are literal `K4` matching triples: put the three values
on edges `01,02,03` and put `1` on `12,13,23`.

Therefore the complete robust sidecar over the 24-valued clutch is

```text
minimal matching face F + its contracted initial-form sum sigma.          (5)
```

If the full matching valuation vector is already retained, only `sigma` is
new.

## Finite residue fields and the roles of two and three

For `k=F_q`, the number of ordered nonzero residue words on a face whose sum
vanishes is

```text
|F|=1: 0,
|F|=2: q-1,
|F|=3: (q-1)(q-2).                                      (6)
```

Hence characteristic two makes the equal two-face maximally cancellation-
sensitive but forbids a three-nonzero-term zero sum over `F_2`; characteristic
three admits the all-ones three-face hostile.  These are additive residue
phenomena.  They are distinct from THM-3049's multiplicative square/cube
divisibility clutches, even though the same primes `2,3` recur.

There is a useful sharp specialization boundary.  THM-3046's quartic
root-difference products satisfy the oriented Pluecker identity

```text
P_1-P_2+P_3=0.
```

For the signed matching triple `(P_1,-P_2,P_3)`, the minimum is therefore
attained at least twice and its initial-face sum is forced to vanish.  This is
an exact hostile on the zero-`sigma` wall, not a nonvanishing mechanism.  An
arbitrary THM-2290 selected `K4` kernel does not inherit this Pluecker
identity.

## K4, selectors, and the V4 torsor

For a fixed THM-2290 endpoint-colour query with selected scalar kernel `B`,

```text
haf(B)=B_01 B_23+B_02 B_13+B_03 B_12.                   (7)
```

Apply (1) to these **already aggregated and selected** kernel entries.  If a
selected pair entry is itself a sum of parallel edges, its valuation and
initial unit must be computed after that inner sum.  Applying (1) before
endpoint-colour selection or parallel-edge aggregation repeats THM-2290's
selector/phase-erasure error.

The `S4` action permutes the three matching terms through `S4/V4=S3`.
Accordingly `(F,sigma)` is equivariant, and the predicate `sigma!=0` is
permutation invariant.  The `V4` kernel acts trivially: this sidecar lives on
the matching quotient and does not choose one of THM-2598's four quartic
sections.  A common vertex-unit gauge multiplies all three matching monomials,
and hence `sigma`, by the same residue unit; nonvanishing is gauge invariant.

This proves no quartic/Keller realization and supplies no endpoint-colour
selector.  It is only the exact additive sidecar once the physical selected
matching triple is already present.

## Exact finite controls

Run

```text
python .scratch/k4_hafnian_initial_face_sidecar.py
python -O .scratch/k4_hafnian_initial_face_sidecar.py
```

The companion exhausts the nonzero triples in
`Z/2^6`, `Z/3^4`, `Z/5^3`, and `Z/7^2`; groups them by full valuation and
leading-residue data; and verifies:

- `sigma!=0` groups contain only valuation-exact nonzero sums;
- every `sigma=0` group below the truncation ceiling contains both a zero and
  a higher-valuation nonzero lift;
- the ceiling groups contain only truncated zeros;
- the finite-field counts (6), the deep hostile (3), and the clutch hostile
  (4);
- all `24` sheet relabellings, the `S3` matching image and `V4` kernel; and
- the common vertex-unit gauge law.

The computation is a control for the general proof, not its source.

Fresh normal and optimized executions both LF-byte-match the stored `23`-line
transcript.  The current evidence hashes are

```text
script SHA-256
  0f8cf744883e537e2a288f4309832531111b32a53a74e4d2f632a6b48f3c7b7c
output SHA-256
  c5659f502635a63e17b1e3b55f3d301ea9f888a720cd97fd5d9e7e84239b3164
semantic SHA-256
  4d3b5150fbef434ec5746dfdc12e90482275b0b5f09ffc07af81109e34f86d3e
```
