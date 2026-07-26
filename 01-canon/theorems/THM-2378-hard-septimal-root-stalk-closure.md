---
id: THM-2378
title: "Hard septimal root-stalk closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the
  hard W=k=1 septimal lane of THM-2367/2372, the pointwise
  thirteen-root stalk is stronger than the signed event-capacity and
  nested-toothpick conclusions. After dividing the three blockers by
  thirteen, the low blocker comb would have to lie in the union of the
  two high blocker combs: on a low-danger/high-safe base phase, the low
  blocker contributes one on all thirteen inverse roots, so exactness
  away from the top-unit absorber would force the guard-danger root
  word, of size three or four, inside the top-unit danger word, of size
  one or two. Independently, if the low quotient has septimal valuation
  below M and both high quotients have valuation above M, every generic
  7^(M+1)-root fibre has exactly N/7 low-danger roots while both high
  masks are constant. The high-safe base set has mass at least 5/7,
  giving an explicit uncovered floor 5/49. Hence the entire hard W=k=1
  lane is empty uniformly. This removes that septimal sublane from both
  strict and repeated-first role branches, but excludes no one of the
  165 thirteen-adic rows by itself, lands no target, and does not prove
  LRC(14).
source: codex-2026-07-25-hard-septimal-root-stalk-closure
depends_on:
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-2372-hard-septimal-signed-stalk-and-toothpick-divisibility
related:
  - THM-1156-tooth-seam-chi7-bipartition
  - THM-2135-root-profile-invoices-and-first-deep-scalar-tail-closure
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2368-owner-pivot-root-fibre-radon-invertibility
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
script: 04-computation/lrc14_hard_septimal_root_stalk_closure_thm2378.py
output: 05-knowledge/results/lrc14_hard_septimal_root_stalk_closure_thm2378.out
script_sha256: 10e7476cc8384d0573d6bd4032d96dd93b8faf22ec7e55b501b71cf953cc4f42
output_sha256: ab4cea7518582b01f5cccfacae49f338edfa1ab1f49312fb58aac01de9a63821
hash_basis: working-tree bytes (LF)
---

# THM-2378 -- the hard septimal root stalk cannot exist

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2372 extracts a mean-zero signed defect, an additive cubic, and a
ninety-eight-fold nested toothpick from the hard lane. The full pointwise
root stalk is more rigid than each of those quotients:

```text
exact lower partition away from three absorbers
  + one low blocker constant on thirteen roots
  + guard root support 3 or 4
  + top-unit root support 1 or 2
  -> low blocker supported by the two high blockers;

strict septimal depth gap
  -> that two-high-blocker support is quantitatively impossible.       (1)
```

The first line uses the `13`-root phase which support-only or integrated
Fourier summaries discard. The second uses the `7^(M+1)` inverse fibre.
The two primes have distinct jobs and no Chinese-remainder averaging is
performed between them.

## 1. Hard-lane notation and exact signed current

Retain the hard `W=k=1` alternative of THM-2367/2372. Thus a canonical
first-depth-one scalar cover has

```text
M=max(nu_7(H),nu_7(q_1),...,nu_7(q_5))>0,

nu_7(H)<M,

nu_7(q_*)=M,

nu_7(q_i)<M                         for q_i!=q_*,

nu_7(c_*)<M,

nu_7(c_3)>M,                        c_* in {c_1,c_2},

nu_7(a),nu_7(b)>M,                  {a,b}={c_1,c_2,c_3}\{c_*}. (2)
```

The six unit labels `H,q_1,...,q_5` are units modulo thirteen, while all
three blockers are divisible by thirteen. Put

```text
D_v={x:||vx||<1/14},

E_H={x:||Hx||<1/7},

mathcal A=D_(q_*) union D_a union D_b,

L
 =1_(E_H)
  +sum_(q_i!=q_*)1_(D_(q_i))
  +1_(D_(c_*)),

F=L-1.                                                       (3)
```

The total mask mass is

```text
2/7+4/7+1/7=1,
```

so

```text
integral_T F=0.
```

THM-2367 proves the load-bearing pointwise identity

```text
F=0                         almost everywhere on mathcal A^c. (4)
```

For completeness, write the lower roles as `(v_i,s_i)`, where

```text
(v_0,s_0)=(H,2)
```

and the other five widths are one. The distributional derivative

```text
mu_L=F'
 =sum_i sum_(m=0)^(v_i-1)
   (
    delta_((m-s_i/14)/v_i)
    -delta_((m+s_i/14)/v_i)
   )                                                        (5)
```

has the exact Fourier transform

```text
mu_L^hat(n)
 =2i sum_(i:v_i|n)
      v_i sin(pi s_i n/(7v_i))                    (n in Z), (6)
```

with value zero at `n=0`. Equivalently, for `n!=0`,

```text
F^hat(n)
 =sum_(i:v_i|n)
    sin(pi s_i(n/v_i)/7)/(pi(n/v_i)).                       (7)
```

Equation (4) says more than an event count. If `X` is the set of distinct
lower endpoint nodes outside `closure(mathcal A)` and

```text
c_x
 =# left walls at x - # right walls at x,
```

then

```text
c_x=0                            for every x in X.             (8)
```

The equivalent signed stalk moments are

```text
M_n=sum_(x in X)c_x exp(-2 pi i n x)=0             (n in Z). (9)
```

If `J=|X|`, the finite block `M_0,...,M_(J-1)` already implies (8) by
the Vandermonde determinant on the distinct roots `exp(-2 pi i x)`.
Thus (9) retains event sign and phase, unlike the unsigned capacity
inequality in THM-2367.

The closure below uses the stronger value identity (4), not merely (8).
This distinction is why it sees an obstruction that the event-current
and cubic quotients do not.

## 2. The thirteen-root stalk forces low-blocker containment

Write

```text
c_*=13C,                 a=13A,                 b=13B.        (10)
```

Fix a base phase `y` outside the finite pullback of every mask endpoint
under all thirteen maps

```text
y -> (y+h)/13,                         h in F_13,
```

and outside the endpoints of the three quotient combs in (10). Write its
thirteen inverse roots as

```text
x_h=(y+h)/13.                                             (11)
```

On this generic set, (4) holds pointwise and every relevant mask is
locally constant. The blocker masks are constant on each root fibre:

```text
1_(D_(c_*))(x_h)=1_(D_C)(y),

1_(D_a)(x_h)=1_(D_A)(y),

1_(D_b)(x_h)=1_(D_B)(y).                              (12)
```

For every thirteen-unit ordinary speed `q`, the translated thirteen-grid
meets a danger arc of length `1/7` in one or two points. For the guard,
it meets the arc of length `2/7` in three or four points. Hence

```text
#{h:x_h in D_q} in {1,2},

#{h:x_h in E_H} in {3,4}.                             (13)
```

Now suppose

```text
y in D_C minus (closure(D_A) union closure(D_B)).     (14)
```

Then the low blocker contributes one to `L(x_h)` for every `h`, while
the two high blockers are absent on the entire root fibre. If
`x_h notin D_(q_*)`, then `x_h notin mathcal A`; equations (3)--(4)
give

```text
L(x_h)=1.
```

The low-blocker contribution is already one, so every one of the five
lower unit masks, including the guard, must vanish at that root.
Consequently

```text
{h:x_h in E_H} subset {h:x_h in D_(q_*)}.             (15)
```

Equation (13) makes (15) impossible:

```text
3<=#{h:x_h in E_H}
  <=#{h:x_h in D_(q_*)}<=2.                           (16)
```

All excluded base phases form a finite null set. If the strict set in
(14) were nonempty, it would contain an interval and hence a generic
phase, so (16) proves, up to immaterial boundary conventions,

```text
D_C subset closure(D_A) union closure(D_B).           (17)
```

This also explains the exact coordinate lost by the event-capacity proof:
the low blocker is constant on the `13`-root fibre, whereas the top unit
has at most two occupied roots and cannot hide the guard's three-root
support.

There is a secondary capacity consequence. On a high-blocker-safe base
phase, (17) makes `C` safe. If both `q_*` and one ordinary lower `u` are
base-danger, their root words are singletons. Exactness on the twelve
`q_*`-safe roots would require at least twelve lower incidences, but the
guard, the distinguished `u`, and the other three ordinary lower masks
have at most

```text
4+1+3*2=11.
```

Thus

```text
D_(q_*) intersection D_u
 subset closure(D_A) union closure(D_B)               (17a)
```

for each of the four ordinary lower labels. The same count with a
base-danger guard gives `3+4*2=11`, hence the analogous guard
containment. There are four ordinary analogues, or five pair containments
total including the guard.

## 3. A septimal anti-shield lemma with a `5/49` floor

The containment (17) is incompatible with the hard valuation gap.

> **Septimal anti-shield lemma.** Let `C,A,B` be positive integers and
> let `M>=1` satisfy
>
> ```text
> nu_7(C)<M<min(nu_7(A),nu_7(B)).                     (18)
> ```
>
> Then
>
> ```text
> measure(
>   D_C minus (closure(D_A) union closure(D_B))
> )>=5/49.                                           (19)
> ```

**Proof.** Put

```text
N=7^(M+1).
```

Then `N|A,B`. On the inverse fibre

```text
y_j=(z+j)/N,                          j in Z/NZ,      (20)
```

the two high masks are constant:

```text
1_(D_A)(y_j)=1_(D_(A/N))(z),

1_(D_B)(y_j)=1_(D_(B/N))(z).                         (21)
```

Write `r=nu_7(C)<M`. As `j` varies, `C(z+j)/N` traverses a
`7^(M+1-r)`-grid, each point with multiplicity `7^r`. The grid order is
divisible by seven. Away from its finitely many aligned endpoints, an
open arc of length `1/7` contains exactly one seventh of the grid.
Therefore

```text
#{j:y_j in D_C}=N/7                                  (22)
```

for almost every `z`.

Let

```text
G=T minus (closure(D_(A/N)) union closure(D_(B/N))).
```

Each danger comb has mass `1/7`, so

```text
measure(G)>=1-1/7-1/7=5/7.                           (23)
```

For `z in G`, all `N` roots in (20) are high-safe, while exactly `N/7`
are low-danger. Disintegration over the map `y -> Ny` gives

```text
measure(
  D_C minus (closure(D_A) union closure(D_B))
 )
 =(1/N) integral_G (N/7) dz
 =measure(G)/7
 >=5/49.                                             (24)
```

This proves the lemma.

In the hard lane, division by thirteen does not change septimal
valuation. Equations (2) and (10) satisfy (18), while (17) says the left
side of (19) is zero. This contradiction proves:

> **Hard-lane closure.** No genuine scalar cover lies in the hard
> `W=k=1` septimal lane of THM-2367/2372.

The quantitative floor `5/49` is robust: endpoint conventions and null
wall coincidences cannot repair the contradiction.

## 4. A useful two-comb sidecar

The same inverse-fibre calculation quantifies the pair containments
proved in (17a). More generally, if

```text
nu_7(u)<M=nu_7(q)<min(nu_7(A),nu_7(B)),
```

then every generic `N=7^(M+1)` fibre has exactly

```text
N/49
```

roots in `D_q intersection D_u`: `D_q` selects one of the seven first
residue classes, and inside it `D_u` selects one seventh of the remaining
grid. Consequently

```text
measure(
 D_q intersection D_u
 minus (closure(D_A) union closure(D_B))
)>=5/343.                                            (25)
```

Thus even one of the pair inclusions

```text
D_(q_*) intersection D_u
 subset closure(D_A) union closure(D_B)
```

would already contradict the hard depths. This `5/343` statement is also
a standalone anti-shield lemma whenever such a pair containment is
available; it does not require the low-blocker argument. In the present
lane, (17a) supplies it. The low-blocker containment (17) is simpler and
gives the stronger floor (19).

The compact-subgroup cross lemma of THM-2135 sends such a pair
containment only to

```text
gcd(q_*,u)|A             or             gcd(q_*,u)|B.
```

That divisor colouring loses the decisive fact that **both** targets are
one full septimal layer above `q_*`. The exact hostile

```text
D_1 intersection D_7 subset D_8 union D_105
```

from THM-2135 lies precisely outside the new hypothesis:

```text
(nu_7(7),nu_7(1),nu_7(8),nu_7(105))=(1,0,0,1),
```

so its two targets are not both strictly above depth one.

## 5. Toothpick divisibility: true boundary and why it is no longer needed

The signed event route in THM-2372 remains correct but is now a weaker
necessary consequence of an empty branch. Its fourteen-fold factor has
the following sharp elementary boundary.

> **Universal wall-absorption lemma.** Every boundary atom of `D_v` lies
> in the closed comb `closure(D_w)` if and only if
>
> ```text
> v|w
> ```
>
> and
>
> ```text
> w/v congruent 0,+1,or -1                 (mod 14). (26)
> ```

Indeed, as the tooth index varies, the target phases form a
`v/gcd(v,w)`-grid. No nontrivial grid lies in a closed arc of length
`1/7`, so `v|w`. The common endpoint phase is then `+/- w/(14v)`,
which is in the closed danger arc exactly for the three residues in
(26).

The tempting unconditional inference `14v|w` is false: the ratios
`14k+1` and `14k-1` are exact boundary hostiles. But if

```text
nu_7(w)>nu_7(v),
```

then `7|(w/v)`, and (26) leaves only residue zero. Thus

```text
14v|w.                                               (27)
```

This is the local toothpick mechanism behind THM-2372. It cannot be
applied to a distributed union of absorbers without the event-current
argument. The root-stalk closure instead proves that the distributed
two-high-absorber configuration cannot exist at all.

## 6. Exact hostile controls

The event-capacity-passing but noncover shield from THM-2367 is caught by
one explicit base phase. In its notation,

```text
q_*=7,

c_*=13,

a=98*13^2*60,

b=13a.
```

After division by thirteen,

```text
(C,A,B)=(1,76440,993720).
```

At

```text
y=1/99
```

the five relevant distances are

```text
||7y||=7/99<1/14,

||2y||=2/99<1/14,

||Cy||=1/99<1/14,

||Ay||=4/33>1/14,

||By||=14/33>1/14.                                  (28)
```

On the thirteen inverse roots, the top-unit danger word has size one,
the guard-danger word has size three, and both high blockers are absent.
More explicitly, their root supports are

```text
q_*:       {0},

guard:     {0,1,12}.
```

Thus roots `1` and `12` are literal pointwise witnesses to (16). The old
event screen passes (`28<=63`), so (28) proves that the new obstruction
is strictly stronger rather than a repackaging of capacity.

THM-2372's forced nested carrier is not a hostile to this closure. It is
a conditional one-factor carrier obtained after forgetting the full
root partition. Its exact circulancy explains why it cannot land a
target, while the present theorem shows that it cannot be completed to
a hard-lane scalar cover.

## 7. Scope and remaining frontier

This theorem removes the hard `W=k=1` septimal sublane uniformly:

- on a strict row it removes both choices of the unique low blocker,
  including the source-only `c_1` role which THM-2367 could not send to
  an excluded target;
- on a repeated-first row it removes both `c_1/c_2` choices inside the
  same hard sublane, where previously only `c_3` had the strict
  shallow/deep carrier;
- it does **not** remove the `k=2` alternatives

  ```text
  (t,b) in {(1,0),(1,1),(2,0),(5,2)},
  ```

  including the saturated `W=7` case `(5,2)`, nor configurations outside
  `M>0`, `nu_7(H)<M`, and the only-`c_3` dominant reduction.

Thus the strict and repeated-first role branches are narrowed but not
globally closed by this theorem. The `165` ledger counts thirteen-adic
rows, not septimal realizations. Therefore no one row is removed merely
by emptying this septimal sublane. The theorem lands no THM-2365 target
and does not solve cancellation in the remaining roles.

The connection to THM-2368 is exact: its missing sidecar was the rooted
phase/event word before chamber integration. Here the `13`-root phase
and blocker-constant bit are retained, and the contradiction occurs
before integration; the all-modes-but-zero-drift hostile therefore does
not apply. In THM-2370 language, no deletion-layer target need be
inferred, because the alleged final hard packet already violates its
pointwise root partition. A tournament on handoff labels would lose
precisely this fibre coordinate.

No scalar row is excluded, the ledger remains `165`, and LRC(14)
remains open.

## 8. Exact companion

The dependency-free companion:

- checks `11,520` generic thirteen-root count cases for both ordinary and
  guard widths;
- checks the support inequality `3>2`;
- checks the exact `N/7` one-comb count on `60` septimal inverse fibres
  and the `N/49` pair count on `360` fibres through depth four;
- records the measure floors `5/49` and `5/343`;
- replays the exact `y=1/99` shield witness, including root support sizes
  one and three;
- audits the universal-wall residues and the THM-2135 hostile boundary.

Run

```bash
python3 04-computation/lrc14_hard_septimal_root_stalk_closure_thm2378.py
python3 -O 04-computation/lrc14_hard_septimal_root_stalk_closure_thm2378.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_hard_septimal_root_stalk_closure_thm2378.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

## 9. Independent hostile audit

An independent read-only audit reconstructed the decisive argument from
the pointwise identity rather than trusting the companion. It checked:

- the finite exceptional set needed to transport the almost-everywhere
  identity to every generic thirteen-root stalk;
- the exact support squeeze `3 or 4 <= 1 or 2`;
- the `N/7` inverse-fibre count, constancy of both higher masks, and the
  resulting `5/49` uncovered floor;
- the explicit `y=1/99` hostile witness; and
- the strict scope `M>0`, `nu_7(H)<M`, `nu_7(c_3)>M`, and
  `c_* in {c_1,c_2}`.

The audit requested four textual repairs: make the endpoint pullback
explicit, distinguish four ordinary pair containments from five total,
state the `5/343` result as a standalone sidecar, and repeat the full
hard-lane scope at the conclusion. All four are incorporated above.
Normal and optimized companions byte-match the stored transcript. QED.
