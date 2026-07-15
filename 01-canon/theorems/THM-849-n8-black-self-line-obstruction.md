---
id: THM-849
title: The n=8 obstruction to the black self-line equals SC law
status: PROVED FINITE-EXACT n=3..8 + PROVED ALL-n KLEIN AND REALIZER IDENTITIES
source: codex-2026-07-15-S15
depends_on: [THM-549, THM-644, THM-781, THM-793, THM-796, THM-843]
related: [THM-409, THM-587, THM-830, HYP-6885, HYP-6890]
verification:
  - 04-computation/n8_black_self_line_obstruction_codex_S15.py
  - 05-knowledge/results/n8_black_self_line_obstruction_codex_S15.out
---

# THM-849 - the n=8 black self-line obstruction

The proposed identity

```text
2 selfK(n) = SC(n)                                            (1)
```

is false.  It holds at `n=5,6,7`, but at `n=8` its two sides are

```text
2 selfK(8) = 404,                  SC(8) = 176.                (2)
```

The failure is not blue contamination or an automorphism-weighting artifact.
There are 412 size-eight fixed-path tilings whose all-tile flip is isomorphic
to the original tournament.  Eight are blue and 404 are black.  All 404
black tilings carry asymmetric tournaments, so their realizer-weighted and
unweighted counts agree.

This refutes HYP-6890's all-size sharpening.  It does not decide the narrower
odd-size conjecture in THM-793/HYP-6885; `n=9` remains its next test.

## 1. The three different fixed-point questions

Let

```text
X_n = F_2^{S_n},                  |S_n|=binom(n-1,2),
kappa(t)=t+1,                    all-tile or path-complement flip,
sigma(t)_(x,y)=t_(n-y+1,n-x+1), staircase reflection.         (3)
```

Write `T_t` for the fixed-path tournament represented by `t`, and let

```text
q_n:X_n -> I_n,                  q_n(t)=[T_t]                  (4)
```

map to the **ordinary**, not converse-merged, tournament class.  Then

```text
kappa sigma = sigma kappa,
q_n(sigma t)=q_n(t)^op.                                      (5)
```

Define the quasi-fixed endpoint sets

```text
Q_n     ={t:q_n(kappa t)=q_n(t)},
Q_n^B   =Q_n intersect Fix(sigma),
Q_n^K   =Q_n minus Fix(sigma).                               (6)
```

The complement line `{t,kappa t}` is an ordinary self-line precisely when
`t in Q_n`.  Its colour is blue in `Q_n^B` and black in `Q_n^K`.  Therefore

```text
|Q_n^B|=2 selfB(n),                 |Q_n^K|=2 selfK(n).       (7)
```

These are not the same fixed-point problem as

```text
SC(n)=#{C in I_n:C=C^op}.                                  (8)
```

Equation (6) counts observer-path orbits over **all** tournament classes.
Equation (8) counts class orbits fixed by converse.  A black self-line carrier
need not be self-converse; this already occurs extensively at `n=7`.

## 2. What the Klein four-group really proves

For every `n>=3`, `V=<kappa,sigma>` is a Klein four action on `X_n`.
The involution `kappa` has no fixed tiling.  Nor does `kappa sigma`: the apex
cell `(n,1)` is fixed by staircase reflection, while `kappa` toggles its bit,
so fixedness would require `t_(n,1)=1-t_(n,1)`.

The set `Q_n` is `V`-invariant.  Invariance under `sigma` follows from (5):
conversing two isomorphic tournaments preserves their isomorphism.  Hence

```text
every V-orbit in Q_n^B has size 2: {t,kappa t},
every V-orbit in Q_n^K has size 4:
       {t,kappa t,sigma t,kappa sigma t}.                    (9)
```

Consequently the genuine all-size laws are

```text
|Q_n^B| is divisible by 2,
|Q_n^K| is divisible by 4,
selfK(n) is even.                                           (10)
```

The Klein action gives no equality with `SC(n)`.  In particular, at `n=8`,
the black quasi-fixed locus consists of 101 free Klein orbits and hence 202
black self-lines.  The earlier values happened to give 2, 3, and 22 free
black Klein orbits at `n=5,6,7`.

## 3. The exact automorphism-weighted identity

For `t in X_n`, put

```text
R_kappa(t)={g in S_n:g T_t=T_(kappa t)}.                    (11)
```

This set is empty off `Q_n`; on `Q_n` it is a torsor for `Aut(T_t)`.  Thus,
for either colour `c in {B,K}`,

```text
P_kappa^c(n)
 =sum_(t in Q_n^c)|Aut(T_t)|
 =sum_(t in X_n, colour c) |R_kappa(t)|
 =sum_(g in S_n) #{t in X_n of colour c:gT_t=T_(kappa t)}.  (12)
```

There is an equivalent Hamiltonian-path form.  For a tournament class `C`,
let `h_kappa^c(C)` be the number of directed Hamiltonian paths `P` in a
representative `T` for which the path-complement tournament `C_P(T)` is
isomorphic to `T` and the corresponding tiling has colour `c`.  Automorphisms
act freely on ordered Hamiltonian paths, because one fixing a path fixes every
vertex.  THM-781 therefore gives

```text
|Q_n^c|       =sum_(C in I_n) h_kappa^c(C)/|Aut(C)|,
P_kappa^c(n) =sum_(C in I_n) h_kappa^c(C).                 (13)
```

This is the correct weighted pair-counting frame.  It also exposes the missing
step in (1): no Klein or torsor argument converts the path-incidence sum (13)
into the number of converse-fixed classes.

At `n=8`, the endpoint automorphism census is

```text
black: 404 x |Aut|=1,
blue:    6 x |Aut|=1, 2 x |Aut|=3.                         (14)
```

Hence `P_kappa^K(8)=404`, `P_kappa^B(8)=12`, and the total pair sum is 416.
Automorphism weighting changes only the blue count and cannot repair (1).

## 4. Independent twisted Burnside count

For completeness, the verifier computes both tournament class counts without
using the fixed-path isomorphism engine.  Let a permutation `pi in S_n` act on
the unordered vertex pairs.  For `e={i<j}`, attach the orientation-transport
sign

```text
s_pi(e)=1[pi(i)>pi(j)].                                    (15)
```

On every induced pair-cycle `C`, ordinary relabelling is consistent exactly
when `sum_(e in C)s_pi(e)=0 mod 2`.  Converse-twisted relabelling adds one at
each step, so its consistency condition is

```text
sum_(e in C)(s_pi(e)+1)=0 mod 2.                           (16)
```

If all pair-cycles are consistent, their bits are independently chosen and
the fixed count is `2^c(pi)`, where `c(pi)` is the number of induced
pair-cycles; otherwise it is zero.  Calling these fixed counts `F_+(pi)` and
`F_-(pi)`, ordinary and twisted Burnside give

```text
T(n)  =(1/n!) sum_(pi in S_n) F_+(pi),
SC(n) =(1/n!) sum_(pi in S_n) F_-(pi).                     (17)
```

At `n=8`, this independent computation gives

```text
T(8)=6,880,                         SC(8)=176.              (18)
```

## 5. Exact census and the first obstruction

The exhaustive fixed-path census is

| `n` | tilings | `T(n)` | `SC(n)` | `|Q_n|` | blue endpoints | black endpoints | `selfB` | `selfK` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| 4 | 8 | 4 | 2 | 2 | 2 | 0 | 1 | 0 |
| 5 | 64 | 12 | 8 | 8 | 0 | 8 | 0 | 4 |
| 6 | 1,024 | 56 | 12 | 16 | 4 | 12 | 2 | 6 |
| 7 | 32,768 | 456 | 88 | 88 | 0 | 88 | 0 | 44 |
| 8 | 2,097,152 | 6,880 | 176 | 412 | 8 | 404 | 4 | 202 |

The exact finite residual has two suggestive lower-level forms:

```text
|Q_8^K|-SC(8)=404-176=228=T(7)/2,
selfK(8)-SC(8)/2=202-88=114,                              (19)
```

and `114` is THM-793's number of converse-merged black self-loops at `n=7`.
Equation (19) is a verified boundary coincidence, not an asserted recurrence.
It points to lower-level merged-loop transport as a concrete correction term
to test, rather than to an all-size cancellation.  The most direct such test
is negative.  Projecting both endpoint-deletion faces on every one of the 101
black Klein orbits gives

```text
55 orbits: every Klein sheet has exactly one ordinary n=7 self-line face,
46 orbits: neither endpoint face is an ordinary or merged-only self-line,
 0 orbits: a merged-only endpoint face.                            (20)
```

Thus the numerical `114` is not a direct lift of the 114 lower merged loops.
Any genuine recursion must retain additional gluing or marked-deletion data.

The computation enumerates all `2^21` tilings.  Sorted scores are used only as
a necessary rejection filter.  Surviving pairs undergo common directed
equitable refinement followed by complete individualization/backtracking, so
no fingerprint is accepted as an isomorphism certificate.  The isomorphism
engine agrees with direct permutation enumeration on every candidate through
`n=6`.  It also reproduces the established `n=3..7` self-line census and the
four known size-eight blue self-lines.

## Tournament Analysis and preservation boundary

The challenged vertex choice is important here.  Runners, arcs, gaps, tilings,
Hamiltonian paths, realizers, complement lines, and isomorphism classes are
different possible tournament vertices; collapsing them caused the false
proof intuition.  For the local audit, take the four operations

```text
{1,kappa,sigma,kappa sigma}                                (21)
```

as vertices on each free black quasi-fixed orbit.  The pair observable is the
difference in `H`; the switch is its sign, with endpoint-word lexicographic
order as the tie gauge.  Isomorphism, converse invariance, and (6) force every
`H` comparison to tie, and the same argument applies to `C3`.  The verifier
additionally finds the score multiset tied on all 101 size-eight orbits,
consistent with THM-790's flip-score constraint.  The lex gauge therefore
gives score histogram
`{0:1,1:1,2:1,3:1}`, no directed cycles, four singleton SCCs, and one
Hamiltonian path.

This operation tournament preserves the four Klein sheets and their free
action.  It destroys which observer-path orbit lands in which ordinary class,
which is exactly the incidence summed in (13).  The false assumption was that
a free Klein orbit count could be identified with the converse fixed-class
count merely because the first three nontrivial sizes agreed.  The `n=8`
obstruction separates those two objects exactly.  No statement here concerns
the LRC metric stalk or proves LRC(14).  ∎
