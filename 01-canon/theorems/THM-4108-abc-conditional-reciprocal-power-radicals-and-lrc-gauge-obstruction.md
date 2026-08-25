---
id: THM-4108
title: "ABC-conditional reciprocal-power radicals and the LRC gauge obstruction"
status: >
  PROVED unconditionally for the reciprocal-power separation, primitive
  normalization, and LRC gauge obstruction; CONDITIONAL ON ABC for radical
  saturation; VERIFIED-EXACT + INDEPENDENTLY AUDITED. Assuming the standard ABC schema, the sum and
  nonzero difference of coprime reciprocal powers have logarithmic radical
  exponent tending to one. The local ABC packet of an LRC straddle is
  lift-gauge dependent, excludes AP13's zero boundary, and admits unbounded
  pair-sum-to-determinant ratio. No ABC or LRC(14) proof is claimed.
source: codex-lrc14-abc-exponent-reciprocity-20260825
depends_on:
  - THM-3833-abc-conditional-cube-radical-and-hyperbolic-power-finiteness
related:
  - THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation
  - THM-668-pair-sum-ruler-witness-structure
  - THM-4105-primitive-reciprocal-phase-descent-and-quantitative-arrival
  - THM-4107-gcd-normalized-exponent-tournament-holonomy-and-lrc-blindness
script: 04-computation/reciprocal_power_abc_lrc_thm4107_4108.py
output: 05-knowledge/results/reciprocal_power_abc_lrc_thm4107_4108.out
script_sha256: 977490728ed77150daf580dbc294f7355bc502cf85c750748e8d4ed532d5691b
output_sha256: 6e91f10ea6db2b117937d439bd7d7a9b300ce3fa31a10b92e3811dee79ecbe32
hash_basis: raw LF bytes
audit: >
  PASS after the positive-difference, local-straddle, lift-sign, and all-real-
  delta quantifier repairs. The independent audit verified the separation
  proof, ABC rearrangements, full limits, gauge hostiles, hashes, and literal
  normal/-O/frozen transcript identity.
---

# THM-4108 -- ABC-conditional reciprocal-power radicals and the LRC gauge obstruction

**PROVED unconditionally in its elementary parts; CONDITIONAL ON ABC for the
radical conclusions; VERIFIED-EXACT + INDEPENDENTLY AUDITED.** ABC gives a strong asymptotic radical
statement for coprime reciprocal powers. It does not transfer to LRC(14): the
natural straddle packet is not invariant under the lift gauge and includes an
unbounded exact hostile family.

The ABC schema is exactly the one recorded in THM-3833 and the primary-source
ledger: for every `epsilon>0`, there is `K_epsilon>0` such that every
pairwise-coprime positive triple `A+B=C` satisfies

```text
C<=K_epsilon rad(ABC)^(1+epsilon).                       (1)
```

ABC remains **OPEN**. No disputed IUT-to-ABC implication is used.

## 1. Reciprocal-power separation

Let `a,b>=2` be distinct coprime integers, and put

```text
P=a^b,
Q=b^a,
H=max(P,Q),
D=|P-Q|,
S=P+Q.                                                    (2)
```

Unconditionally,

```text
boxed: D>=H/9,                                            (3)
```

with equality if and only if `{a,b}={2,3}`.

Assume `a<b`. If `a=2`, coprimality makes `b` odd. At `b=3`, `(H,D)=(9,1)`.
For `b>=5`, `H=2^b`, and `b^2/2^b` decreases from `25/32<8/9`; hence
`D/H>1/9`.

If `a>=3`, monotonicity of `log(x)/x` gives `H=a^b`. For fixed `a`, the ratio

```text
b^a/a^b                                                     (4)
```

decreases for `b>=a+1`. At `a=3` its maximum is `64/81<8/9`; for `a>=4`,
it is at most `(1+1/a)^a/a<e/a<3/4`. This proves `(3)` and its equality case.

Also,

```text
gcd(S,D)=1,        a,b of mixed parity;
gcd(S,D)=2,        a,b both odd.                         (5)
```

Indeed, `gcd(P,Q)=1`, so `gcd(P+Q,P-Q)` divides two, and parity decides.

## 2. Conditional radical saturation

The two additive packets

```text
P+Q=S,
min(P,Q)+D=H                                               (6)
```

are primitive, and `rad(PQ)=rad(ab)`. Assuming `(1)`,

```text
rad(S)>=K_epsilon^(-1/(1+epsilon))
        S^(1/(1+epsilon))/rad(ab),

rad(D)>=K_epsilon^(-1/(1+epsilon))
        H^(1/(1+epsilon))/rad(ab).                       (7)
```

If `n=max(a,b)`, the powered term with exponent `n` gives

```text
H>=2^n,
rad(ab)<=ab<=n^2<=(log_2 H)^2.                           (8)
```

Together with `(3)`, this proves that along every sequence of distinct
coprime pairs for which `H -> infinity`, conditionally on ABC,

```text
boxed:
log rad(D)/log D -> 1,
log rad(S)/log S -> 1.                                   (9)
```

The upper bounds are `rad(D)<=D` and `rad(S)<=S`. The comparisons

```text
H/9<=D<=H,                    H<=S<=2H                  (9a)
```

show that `log D/log H` and `log S/log H` tend to one. For the lower bounds,
fix `epsilon`; equations `(7)--(9a)` leave only a fixed
`H^(epsilon/(1+epsilon))` power and a squared logarithm. Letting the
independently fixed `epsilon` tend to zero gives `(9)`.

Equivalently, for every fixed `delta<1`, only finitely many coprime pairs can
satisfy either

```text
rad(D)<=D^delta                 or
rad(S)<=S^delta.                                        (10)
```

Choose `epsilon` with `1/(1+epsilon)>delta`. Under either inequality in
`(10)`, `(9a)` bounds the right side by `C_delta H^delta`: for negative
`delta`, use the lower endpoint of the relevant interval, so the inequality
direction is correct. Equations `(7)--(8)` would then bound a positive power
of `H` by a squared logarithm. Finally, bounded `H` contains only finitely
many base pairs because `max(a,b)<=log_2 H`. This closes every real
`delta<1`, including negative `delta`.

This is a valuation-support theorem. The constants are ineffective here, and
the conclusion does not supply residues, signs, owners, or a common time.

## 3. Noncoprime boundary and valuation defects

For arbitrary positive integer bases, let

```text
G=gcd(a^b,b^a),             A=a^b/G,             B=b^a/G. (11)
```

At every prime `p`, put

```text
Delta_p=b v_p(a)-a v_p(b).                              (12)
```

Then

```text
A=product_p p^max(Delta_p,0),
B=product_p p^max(-Delta_p,0),
rad(AB)=product_(Delta_p!=0) p.                         (13)
```

ABC may be applied to the normalized sum packet. It may be applied to the
normalized difference packet only when `A!=B`; otherwise the difference is
zero and not a positive ABC triple. It gives no saturation unless `rad(AB)`
is separately sub-power in the normalized height. The tie `2^4=4^2`, for
which `A=B=1`, is the minimal hostile: all apparent raw height disappears.

Thus radical support is the zero/nonzero shadow of the full signed valuation
defect; it loses both depth and slot.

## 4. A local pair-straddle packet and its exact obstruction

Let `u,w` be distinct positive speeds and `alpha,beta` integer lifts. At an
oriented local pair-sum straddle, write

```text
q=u+w,
p=alpha+beta,
D=u beta-w alpha,                 0<D<=q/2,
t=p/q,
margin=D/q.                                               (14)
```

THM-4105 proves the exact reciprocal exponent factorization of this signed
determinant. To expose the ABC packet, put

```text
X=u beta,              Y=w alpha.                        (15)
```

If `X,Y>0`, then with `G=gcd(X,Y)`,

```text
min(X,Y)/G + D/G = max(X,Y)/G                            (16)
```

is a primitive ABC triple. This is a lawful conditional consumer. It cannot
imply the LRC pair requirement `q<=14D`, for four exact reasons.

### Zero boundary

The canonical AP13 straddle

```text
(u,w,alpha,beta)=(1,13,0,1)                              (17)
```

is the degenerate packet `0+1=1`, outside positive ABC.

### Lift-gauge obstruction

For every integer `N`, the same physical circle point may be represented by
`t+N`. Its lifts change by

```text
(alpha,beta) -> (alpha+Nu,beta+Nw).                      (18)
```

Consequently,

```text
(X,Y) -> (X+Nuw,Y+Nuw),                                  (19)
```

while `q`, `D`, the margin, and every physical LRC phase remain fixed. ABC
height and radical are therefore not intrinsic to the straddle unless one
retains an arbitrary lift gauge. For AP13 the positive packets, which require
`N>=1`, form the family

```text
13N+1=(13N)+1.                                           (20)
```

At `N=0` this is the recorded zero boundary; negative `N` gives signed terms,
not a positive ABC packet.

### Unbounded local hostile

For every `n>=1`, take

```text
(u,w,alpha,beta)=(n+1,n,1,1).                            (21)
```

Then

```text
q=2n+1,
p=2,
D=1,
t=2/(2n+1),
margin=1/(2n+1),                                         (22)
```

and `(16)` is just the consecutive-integer ABC packet

```text
(n)+(1)=(n+1).                                           (23)
```

Thus `q/D=2n+1` is unbounded and already exceeds `14` for `n>=7`. This is a
local two-speed straddle, not a claim that the pair is a global active
maximizer in a thirteen-speed row. ABC is fully compatible with its pair
margin lying below the LRC(14) target.

### Missing global predicate

Even a favorable radical-height bound for `(16)` says nothing about whether
the other eleven speeds clear `t=p/q`, whether this local straddle is the
global maximizer, or whether the owner packet is physically synchronized.

## 5. Transfer ledger and generated tasks

```text
source                 preserved                    lost
reciprocal powers      primitive height/radical     phase, owner, arrival
valuation defect       prime support and depth      support projection loses depth
LRC straddle           signed local determinant     lift-invariant ABC height absent
gcd-normalized edge    common-dilation class        LRC margin (THM-4107 hostile). (24)
```

The useful next proposals are therefore conditional and owner-aware:

1. Apply ABC only to actual positive three-term relations already found by
   the short-Graver compiler; retain coefficient slots, valuation depths,
   physical phase, and arrival.
2. In the live finite atlas, split low-radical fibres from rich-radical
   fibres. Use low radical only as a conditional height gate; test whether
   rich radical supplies a modulus with a missing signed blocker card.
3. Factor `a^b-b^a` for coprime pair edges, but require one prime to select a
   common owner-safe residue for the whole row. Aggregate radical abundance
   alone is not that selector.
4. Use common-dilation rays and the lift family `(20)` as cheap hostile tests
   for every proposed addition/multiplication invariant.

No statement here proves ABC, validates IUT, supplies an effective cutoff, or
proves LRC(14).

## 6. Exact audit

Reproduce with

```text
python -B 04-computation/reciprocal_power_abc_lrc_thm4107_4108.py
python -B -O 04-computation/reciprocal_power_abc_lrc_thm4107_4108.py
```

The normal and optimized streams equal the frozen transcript. Exact integer
arithmetic checks `705` coprime sum/difference channels, `3,540` valuation
packets, `28,794` normalized straddle packets, the AP13 lift gauge, and the
unbounded consecutive-integer hostile. The ABC implications themselves are
proved symbolically from the explicitly assumed schema, not by computation.
