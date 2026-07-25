---
id: THM-2174
title: "Endpoint-phase scale obstruction to target-preserving radix pumping"
status: >
  PROVED + VERIFIED-EXACT. For a fixed finite core, the THM-2162 endpoint
  phase word is periodic in the added speed W, but its signed cocycle is
  C_r/W on residue class r and therefore retains unbounded scale whenever
  C_r is nonzero. An explicit defect-six-shaped thirteen-speed base-two pump
  preserves a seven-speed core, strict order, primitivity, owner and quotient-
  tie states, three independent height-one crossing relations, and every
  far endpoint phase modulo the full core denominator 1680, while changing
  the exact 1/14-safe-set measure. Thus endpoint phase labels without current
  magnitude are not a target-preserving sidecar for THM-2171. This is a sharp
  carrier obstruction, not a counterexample to LRC(14) or a proof that no
  other finite target invariant exists.
source: codex-2026-07-24-LRC-post-2168-closure
depends_on:
  - THM-2162
  - THM-2171
related:
  - THM-2166
  - THM-2168
  - THM-2169
script: 04-computation/lrc14_endpoint_phase_scale_obstruction_thm2174.py
output: 05-knowledge/results/lrc14_endpoint_phase_scale_obstruction_thm2174.out
script_sha256: 561050e451f5d60526381e3c750da9fd45904cac9ce4af92a9263682bba552ce
output_sha256: ba4c1eba837c247b64ce6d6c407dbf1f49b36bb2edd4e46612fc7a319c3c2e49
hash_basis: working-tree bytes (LF)
---

# THM-2174 -- endpoint-phase scale obstruction

At radius `1/14`, write

```text
G_E={t in R/Z:||et||>=1/14 for every e in E}.         (1)
```

THM-2171 makes bounded relation systems pumpable without losing positivity,
order, distinctness, or primitivity. This theorem proves that even the
natural endpoint-phase sidecar from THM-2162 is not enough to preserve the
LRC target: its phase labels are finite, but its signed magnitude remembers
the unbounded speed scale.

## 1. Endpoint phases versus endpoint current

Let the positive-length components of `G_E` be

```text
I_s=[l_s,r_s],                    s=1,...,K,           (2)
```

with harmless choices at their measure-zero endpoints. Choose an integer
`L_E` divisible by every denominator of every `l_s,r_s`. Let `H` be the
continuous periodic primitive from THM-2162:

```text
H'=1_(-1/14,1/14)-1/7,             H(0)=0.            (3)
```

For a residue `r mod L_E`, define the endpoint numerator

```text
C_E(r)=sum_(s=1)^K [H(r r_s)-H(r l_s)].              (4)
```

THM-2162's exact endpoint identity immediately gives:

> **Phase-scale lemma.** If `W` is a positive integer and
> `W=r mod L_E`, then
>
> ```text
> epsilon_W(G_E)
> :=|G_E intersect {t:||Wt||<1/14}|-|G_E|/7
> =C_E(r)/W.                                         (5)
> ```

Indeed, `Wl_s` and `rl_s` agree modulo one, as do `Wr_s` and `rr_s`,
while the Jacobian factor in the endpoint identity remains `1/W`.

This separates two coordinates that a finite phase word conflates:

```text
endpoint phase label:       r=W mod L_E;
signed endpoint magnitude:  C_E(r)/W.                (6)
```

If `C_E(r)!=0`, then the speeds `r+nL_E` have one common endpoint phase
word but infinitely many distinct endpoint currents. More sharply, inside
that residue class,

```text
C_E(r)/W=C_E(r)/W'  iff  W=W'.                       (7)
```

Thus preserving the exact current on a nonzero class retains the speed
scale itself.

## 2. A seven-speed core with nonzero phase currents

Take

```text
E=(1,2,3,4,5,6,8).                                   (8)
```

Every boundary of `G_E` is a danger-comb endpoint

```text
(14m+/-1)/(14e),                  e in E,             (9)
```

so one may take

```text
L_E=lcm(14e:e in E)=1680.                            (10)
```

Exact interval arithmetic gives

```text
|G_E|=27/70,             K=16,                       (11)
C_E(1)=-27/490,          C_E(2)=-27/245.             (12)
```

Both residue classes therefore have the infinite scale sensitivity in
(7).

## 3. The relation- and phase-preserving pump

For `a=1,2,3`, put

```text
f_(a,-)=3360a+1,             f_(a,+)=3360a+2,
f'_(a,-)=1680a+1,            f'_(a,+)=1680a+2.       (13)
```

The two ordered thirteen-speed rows are

```text
V =E union (3361,3362,6721,6722,10081,10082),
V'=E union (1681,1682,3361,3362,5041,5042).          (14)
```

Both are primitive because they contain `1`. Apply the base-two THM-2171
pump which deletes radix level `j=4` and splices in level `k=5`. The core
coordinates are below `2^4` and remain fixed, while

```text
3360a+r
 =r+2^5(105a)
 |-> r+2^4(105a)
 =1680a+r,                       r=1,2.              (15)
```

Thus the pump sends `V` exactly to `V'`. Its two boundary quotients are

```text
Z_4=(0,0,0,0,0,0,0,210,210,420,420,630,630),
Z_5=(0,0,0,0,0,0,0,105,105,210,210,315,315).        (16)
```

They have the same owner suffix `{8,...,13}` and the same one-based
quotient-tie cut mask `{7,9,11}`.

There are three independent height-one crossing relations. With `e_i`
denoting the coordinate vectors of the ordered row, set

```text
rho_a=-e_1-e_(6+2a)+e_(7+2a),        a=1,2,3.        (17)
```

For each pair,

```text
-f_(a,-)+f_(a,+)=1,
-f'_(a,-)+f'_(a,+)=1,                              (18)
```

and coordinate `1` pays `-1`. Hence every `rho_a` annihilates both rows.
Their private far-coordinate pairs make them independent over every field.
Equation (16) also shows that their carries at levels four and five are
both zero, so (17) is preserved by the exact THM-2171 splice.

These relations have a stronger version of the coefficient shape forced
by THM-2166:

```text
cut carry=1,
far height=1 with every nonzero far coefficient a 7-unit,
core support=1 and core height=1.                    (19)
```

The rows in (14) are not zero-safe, so (19) is a carrier stress test, not
an invocation or converse of THM-2166.

Finally,

```text
f_(a,+/-)-f'_(a,+/-)=1680a.                          (20)
```

Thus the pump preserves every coordinate modulo `L_E=1680`. In particular,
it preserves the complete endpoint phase word of every far comb against
every component of the fixed core `G_E`; this is exactly THM-2171's optional
finite-residue sidecar with `N=1680`.

## 4. The exact current and target still move

The far speeds in each minus column have residue one modulo `1680`, and
those in each plus column have residue two. Equations (5) and (12) give

```text
epsilon_W(G_E)=-27/(490W),       W=1 mod 1680,
epsilon_W(G_E)=-27/(245W),       W=2 mod 1680.        (21)
```

Therefore every paired old/new far comb has identical endpoint phases but
different exact current.

The failure persists for the full thirteen-speed target, not only for a
one-comb marginal. Two independent exact interval calculations give

```text
|G_V|
 =1561405750435498559/10390707539702618590,

|G_(V')|
 =317645844187362436/2113439446871390435,             (22)
```

and

```text
|G_(V')|-|G_V|
 =126112336463776271349/4405994577616688706478598
 >0.                                                  (23)
```

One evaluator merges the complete family of rational danger intervals.
The other forms the complete rational boundary arrangement and sums exactly
the safe open cells. They agree on both fractions in (22). Normal and
optimized Python produce the stored transcript.

## 5. Consequence and exact residual

The example preserves more than the algebraic state of THM-2171:

```text
fixed core;
positive strict order and primitivity;
owner and quotient-tie state;
three independent height-one crossing relations and their carries;
the full coordinatewise residue word modulo the core endpoint modulus;
every resulting endpoint phase label.               (24)
```

Nevertheless the exact safe measure changes. Therefore:

1. endpoint phase labels, even at the full core denominator, are not a
   target-preserving pump sidecar;
2. the missing datum already appears in the one-comb THM-2162 current as
   the scale factor `1/W`;
3. retaining the exact nonzero current on a fixed phase class is not finite
   state compression--by (7), it retains `W` itself.

This sharpens the post-THM-2171 residual. Algebraic feasibility has a finite
terminal, and any chosen finite residue bank can be preserved, but a
counterexample descent still needs a scale-sensitive phase/current theorem
that preserves the **zero-safe predicate**, not merely its endpoint labels.
The theorem does not say that no other finite target invariant can exist,
does not close THM-2168's `(0,0)`, `(0,1)`, `4+3`, or `5+3` profiles, and
does not prove or refute LRC(14). QED.
