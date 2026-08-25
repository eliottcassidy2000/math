---
id: THM-4025
title: "LRC(14) owner-residue odd-dilation semigroup"
status: >
  PROVED ARITHMETIC + SCOPED LRC REDUCTION + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. THM-4003's owner-relaxed residue margin is monotone under common
  odd dilation in the divisibility order. Under the ordinal chart
  phi(n)=2n-1, odd multiplication becomes n star m=2nm-n-m+1; arithmetic
  survivors form a star-upward ideal and closures a divisor-closed set.
  THM-4003 interprets the gate as an LRC failure obstruction only below
  U=91^6. This structures a certificate search but does not prove LRC(14).
source: root + dilation_semigroup_audit / sequence continuation, 2026-08-24
audit: >
  PASS. A Fraction engine and an independent integer-cross-product engine
  agree on 15,864 physical old-strip cells through t=501, the unique
  ordinary-order revival in that universe, 18,030 survivor-divisibility
  implications, and the exact hostile margins. Normal and optimized outputs
  reproduce byte-for-byte; hashes are raw LF bytes.
depends_on:
  - THM-4003-lrc14-scale-two-component-erosion-boundary-strip
related:
  - THM-3756-odd-square-ordinal-berggren-affine-descent
  - THM-4004-lrc14-three-detuned-divisor-comb-profile
script: 04-computation/lrc14_owner_residue_odd_dilation_semigroup_thm4025.py
output: 05-knowledge/results/lrc14_owner_residue_odd_dilation_semigroup_thm4025.out
independent_audit_script: 04-computation/lrc14_owner_residue_odd_dilation_independent_audit_thm4025.py
independent_audit_output: 05-knowledge/results/lrc14_owner_residue_odd_dilation_independent_audit_thm4025.out
script_sha256: 1aa8e55b816cf1638b310d484f62e603e22f5a8cfbab8013d5a9d9933938ad60
output_sha256: 1f3d5e18026a979044e5fbba3b0487405c5d35cefb18dadcd5cd8df687d12c0c
independent_audit_script_sha256: 186bd0714f42cde506b3b7414cdd6e3c9d64e198a2a384efa304e85e60db0807
independent_audit_output_sha256: b2782cd57ff24e462f399d1bf04c1e5c020bfc77d244abd20d568ec85020b858
hash_basis: raw LF bytes
---

# THM-4025 -- the owner-residue odd-dilation semigroup

**PROVED ARITHMETIC + SCOPED LRC REDUCTION + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.** The theorem identifies the operation that a scalar sequence of
THM-4003 cells must retain. The order is odd divisibility, transported to a
natural-number semigroup. It is not ordinary index order, and an exact
minimal hostile shows that ordinary monotonicity is false. LRC(14) remains
open.

## 1. Inheritance pass and typed arithmetic gate

For odd `t`, `U>1`, and `1<=u<=U`, retain THM-4003's first two directed
owner gaps

```text
e_1(t,u)=[3t-4u]^+_(42 gcd(t,u))/(42u),
e_2(t,u)=[9t+16u]^+_(126 gcd(t,u))/(126u).             (1)
```

The other two gaps equal these by reflection. Put

```text
epsilon_i(t,U)=min_(1<=u<=U)e_i(t,u),
B(t,U)=max(0,2/63-epsilon_1(t,U)-epsilon_2(t,U)),
D(t,U)=t(2U-1)/(84U(U-1)).                            (2)
```

Call `D>B` an **arithmetic closure** and `D<=B` an **owner-relaxed arithmetic
survivor**. This terminology does not call a survivor an LRC counterexample.
In the physical old strip with `U<Q=91^6`, THM-4003 proves that every LRC
failure must satisfy `D<=B`; hence `D>B` closes that cell. Beyond `Q`, `(1)`
and `(2)` remain exact arithmetic objects, but this theorem makes no new LRC
interpretation of them.

The closest proved mechanism is THM-4003's owner minimization. Its canonical
near miss is to read the computed cells in ascending odd order. The least-used
sidecar is common dilation of the owner itself.

## 2. Each scaled owner has exactly the same directed gaps

Let `k` be positive and odd. If `d=gcd(t,u)`, then

```text
gcd(kt,ku)=kd.                                         (3)
```

For every integer `x` and positive `M`, least positive residues obey

```text
[kx]^+_(kM)=k[x]^+_M.                                  (4)
```

The nonzero condition is automatic here: the numerators in `(1)` are odd
while their moduli are even. Every numerator, modulus and denominator in
`e_i(kt,ku)` is `k` times its counterpart in `e_i(t,u)`. Therefore

```text
e_i(kt,ku)=e_i(t,u),                 i=1,2,3,4.         (5)
```

The target owner set `{1,...,kU}` contains every scaled owner `ku` from the
source. Minimizing over its additional owners gives

```text
epsilon_i(kt,kU)<=epsilon_i(t,U),
B(kt,kU)>=B(t,U).                                      (6)
```

This is not coordinatewise strict: a new owner may tie the old record.

## 3. The demand decreases by an exact amount

Direct subtraction in `(2)` gives

```text
D(t,U)-D(kt,kU)
 =t(k-1)/(84(U-1)(kU-1))=:Delta_k(t,U).                (7)
```

Thus `Delta_k>0` for `k>1`, and `(6)--(7)` yield the quantitative margin law

```text
[B-D](kt,kU)>=[B-D](t,U)+Delta_k(t,U).                 (8)
```

Consequently

```text
D(t,U)<=B(t,U)  ==>  D(kt,kU)<B(kt,kU)       (k>1),
D(kt,kU)>B(kt,kU)  ==>  D(t,U)>B(t,U).                 (9)
```

Arithmetic survivors propagate upward under odd dilation; closures propagate
downward to every odd divisor point. Closure is not upward under dilation.

## 4. Treating the n-th odd square as n

Define the ordinal chart and transported product

```text
phi(n)=2n-1,
n star m=(phi(n)phi(m)+1)/2=2nm-n-m+1.                (10)
```

Then

```text
phi(n star m)=phi(n)phi(m).                            (11)
```

The operation `star` is commutative and associative with identity one. The
same chart applies if the label is the square `phi(n)^2`, because multiplying
two odd squares multiplies their odd roots. The square therefore need not be
stored, but the transported product must be.

Fix an arithmetic ray `(t_0,U_0)` and let

```text
c(n)=1[D(phi(n)t_0,phi(n)U_0)>B(phi(n)t_0,phi(n)U_0)]. (12)
```

Here `c=1` means arithmetic closure. Equation `(9)` becomes

```text
c(n star m)<=min(c(n),c(m)).                           (13)
```

Thus the survivor indices form a `star`-upward ideal and the closure indices
form a divisor-closed lower set. From scale `phi(n)` to scale
`phi(n star m)`, the guaranteed margin increase is

```text
phi(n)t_0(phi(m)-1)
----------------------------------------------- .     (14)
84(phi(n)U_0-1)(phi(n)phi(m)U_0-1)
```

This is the operation-aware sequence law. A binary list without `star` loses
the theorem.

## 5. Physical and algorithmic consequence

Odd dilation preserves `U<=t` and the old-strip inequality `4U>3t`; it also
preserves physicality upward from `U>=11`. Whenever the target satisfies
`kU<Q`, `(9)` has the following LRC-certificate readout:

- a known owner-relaxed survivor at one scale lets a search skip every odd
  multiple on that ray, because the same gate cannot close there;
- a closure found at a multiple forces arithmetic closure at each odd divisor
  point, and those physical divisor points below `Q` are already closed by
  THM-4003.

This is a pruning theorem, not a closure of the surviving ray. Scales that
are incomparable by divisibility must still be tested separately.

## 6. Minimal hostile to ordinary sequence order

The primitive arithmetic ray `(t_0,U_0)=(5,4)` becomes physical at odd scale
three. At two adjacent odd scales one finds

```text
(t,U)=(45,36), k=9:
epsilon_1=1/966  (owner 23),
epsilon_2=1/2772 (owner 22),
D=71/2352, B=215/7084, D-B=-97/595056;               (15)

(t,U)=(55,44), k=11:
epsilon_1=1/1722, epsilon_2=17/5166 (both owner 41),
D=145/4816, B=8/287, D-B=63/28208.                   (16)
```

So `(45,36)` is an arithmetic survivor but the next ordinary odd scale
`(55,44)` closes. There is no contradiction to `(9)`: `9` does not divide
`11`. In ordinal indices this is the revival from `n=5` to `n=6`, which are
also incomparable under `star`-divisibility.

An exact scan of all `15,864` physical old-strip cells through `t=501`
groups them into `12,757` primitive rays. Two independent engines find just
this one survivor-to-closure revival, and checking every smaller target makes
`(15)--(16)` the minimal physical ordinary-order hostile. The ray word for
odd scales `3,5,...,31` is

```text
CCCSCCSSSSSSSSS.                                       (17)
```

This finite word is telemetry, not an eventual-tail claim.

## 7. Sharp boundaries to stronger claims

- **Closure is not upward.** `(11,11)` closes, while its triple `(33,33)` is
  an arithmetic survivor.
- **Individual minima need not decrease strictly.** One has
  `epsilon_2(21,16)=epsilon_2(63,48)=1/504`.
- **There is no coordinatewise `1/k` gain.** For example
  `epsilon_1(33,33)=1/588>epsilon_1(11,11)/3=1/1008`.
- **Even dilation is outside the theorem.** It destroys the odd-numerator
  noncollision used by THM-4003.
- **The all-`U` arithmetic law is not an all-`U` LRC certificate.** The cited
  failure implication `D<=B` is canonically typed only for `U<Q`.

## 8. Operation and loss ledger

| source | target/map | preserved | destroyed | needed sidecar |
|---|---|---|---|---|
| owner `u` | scaled owner `ku` | every directed gap `(5)` | new unscaled target owners | full target owner minimum |
| cell `(t,U)` | `(kt,kU)` | old strip and quantitative margin order | closure can die | factor `k` and margin |
| odd multiplier `2n-1` | ordinal `n` | multiplication via `star` | ordinary order has no monotonic meaning | operation `(10)` |
| body | scalar bit `c(n)` | this one residue certificate | owner labels, phase, arrivals, actual body | minimizing owners and `B-D` |

The theorem therefore supports a divisor-DAG or semigroup-indexed search,
not a fitted unary recurrence.

## 9. Scope and replay

The result concerns THM-4003's scale-two owner-relaxed arithmetic gate in the
conditional `t>=U` lane. It supplies no body, phase or LRC witness, and an
arithmetic survivor is not a counterexample. The opposite `t<U` lane, actual
eleven-owner compatibility, all cells beyond the pair-height boundary, and
LRC(14) remain open.

Reproduce from the repository root:

```text
python3 -B 04-computation/lrc14_owner_residue_odd_dilation_semigroup_thm4025.py
python3 -B -O 04-computation/lrc14_owner_residue_odd_dilation_semigroup_thm4025.py
python3 -B 04-computation/lrc14_owner_residue_odd_dilation_independent_audit_thm4025.py
python3 -B -O 04-computation/lrc14_owner_residue_odd_dilation_independent_audit_thm4025.py
```

Both normal/optimized pairs reproduce their frozen outputs. **QED.**
