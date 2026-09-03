# LRC(14): third-tooth competition is a prefix envelope, and nested teeth cast a boundary shadow

**Status:** PROVED elementary prefix-envelope, exact successor-residue,
candidate-sieve, nested-multiple, wall-shadow, and sharp wall-capacity lemmas
+ FINITE-EXACT controls.  The companion prefix-envelope / nested-wall-shadow
theorem is a canon-ready proof candidate; this file records the research interpretation.  The full `420|h`
minority branch and LRC(14) remain open.

## Inheritance pass and concept board

The closest proved mechanism is the shortest addressed cover chain in
[THM-4335](../01-canon/theorems/THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal.md).
The immediately preceding
[address-capacity audit](lrc14-transition-reuse-address-capacity-20260902.md)
proves that a fixed addressed edge occurs on at most one anchor component,
but that one speed pair can occur on unboundedly many components after the
addresses are forgotten.  Its sharp internal hostile is the pair
`{525,1365}`: there are `210` located candidate edges but only two actual
uses in either frozen `h=420` row.

The corrected near miss was to attribute the whole `210 -> 2` collapse to
third-tooth competition.  Exact classification below shows that `208` of the
`210` candidates are not renewal obligations at all: their components are
uncovered or one-tooth-spanned.  The least-used relevant sidecar is therefore
not another pair invariant, but the component's prefix-reach orbit.

The concept board was

```text
located handoff | renewal obligation | prefix-record envelope
centered successor residue | visibility cap | multiplicative nesting
anchor-safe u-wall | residual boundary cover
```

The method cards used were “turn certificate failure into an address,”
“controlled forgetting requires a sidecar,” and “re-evaluate a certificate
after every fibre-changing operation.”  Here the speed-pair quotient forgets
both the component status and the frontier at which a tooth becomes visible.

## 1. Farthest reach is an exact prefix-record process

For an odd speed `w` and tooth address `n`, write

```text
T(w,n)=(a(w,n),b(w,n)),
a(w,n)=(14n-1)/(14w),       b(w,n)=(14n+1)/(14w).     (1)
```

Fix an anchor component

```text
I_k=[L_k,R_k]
   =[(14k+1)/(28h),(14k+13)/(28h)].                  (2)
```

Among the finitely many tail teeth meeting `I_k`, order teeth by the exact
deterministic greedy key

```text
(b(T),-a(T),-w(T),-n(T)).                            (3)
```

Thus farther right wins; at an equal right wall the wider, hence smaller-
speed, tooth wins.  This only fixes genuine endpoint ties and never turns an
open endpoint into an active point.

For a frontier `x`, let

```text
Psi_k(x)=the maximum-key tooth among those with a(T)<x. (4)
```

If `b(Psi_k(x))<=x`, the cover has stopped at a gap.  Otherwise put

```text
x_0=L_k,       T_j=Psi_k(x_j),       x_(j+1)=b(T_j). (5)
```

> **Prefix-envelope lemma.**  Recurrence `(5)` is exactly the
> farthest-reach greedy chain, with the same open-endpoint and tie semantics.
> It reaches beyond `R_k` if and only if the tail teeth cover `I_k`.

Indeed, if some tooth is active at `x`, every already expired tooth has right
wall at most `x`, while every active tooth has right wall greater than `x`.
Hence taking the maximum among all started teeth is the same as taking it
among active teeth.  If no active tooth exists, the prefix maximum has already
expired and the inequality in the stopping rule detects exactly that gap.

There is a useful static form.  Define the **visibility cap**

```text
c(T)=min({b(T)} union {a(Z): Z outranks T and a(Z)<b(T)}). (6)
```

Then `T=Psi_k(x)` precisely when

```text
a(T)<x<b(T)                 and                 x<=c(T). (7)
```

The nonstrict inequality at `c(T)` is forced: a competing open tooth whose
left wall equals `x` is not yet active.  Thus a third tooth truncates a
candidate's visible-frontier interval from the right.

For a located candidate edge `A->B`, actual greedy selection is now an exact
three-stage condition:

1. the component is a genuine covered, nonspanning renewal obligation;
2. `b(A)<=c(B)`, equivalently `B=Psi_k(b(A))`;
3. the prefix orbit `(5)` actually reaches `A`.

These conditions are necessary and sufficient.  Condition 2 handles the
outgoing third-tooth contest; Condition 3 handles an earlier tooth which
jumps over the entire visibility interval of `A`.

## 2. The successor contest is one centered-residue comparison

Let the outgoing tooth be `A=T(u,n)` and put `x=b(u,n)`.  For each competing
speed `z`, take the unique centered residue

```text
e_z = z(14n+1) mod 14u,             -7u<e_z<=7u.     (8)
```

Then the `z`-tooth is strictly active at `x` if and only if

```text
|e_z|<u.                                             (9)
```

In that case its address and its rightward reach are

```text
ell_z=[z(14n+1)-e_z]/(14u),
b(z,ell_z)-x=(u-e_z)/(14uz).                         (10)
```

Now let the candidate successor be `B=T(v,m)`, and retain the usual proper-
crossing determinant

```text
D=um-vn,                    |u-v|<14D<u+v.           (11)
```

For this candidate

```text
e_v=v-14D,
R_v=(u-e_v)/v=(u-v+14D)/v.                           (12)
```

> **Exact third-tooth domination criterion.**  The candidate `B` wins at
> `b(A)` if and only if there is no distinct speed `z` satisfying `(9)` and
>
> ```text
> (u-e_z)/z > R_v,                                   (13a)
> ```
>
> or satisfying equality in `(13a)` with `z<v`.      `(13b)`

Formula `(13b)` is exactly the equal-right-wall rule: the smaller speed has
the wider tooth and wins.  Thus a candidate edge can be checked without
floating point, interval sampling, or a stochastic approximation.  The
companion probe verifies this criterion against literal rational teeth.

## 3. An exact aggregate sieve

Let `E_h(W)` be the located candidate events of the address-capacity audit.
For an event `e`, its handoff window is contained in one anchor component.
Define

```text
O = {e in E_h(W): that component is a renewal obligation},
S = {e in E_h(W): its proposed successor wins (13)},
A = {actual selected greedy edges}.                  (14)
```

Then

```text
A subset S intersect O subset O subset E_h(W).       (15)
```

The remaining reachability condition in Section 1 makes membership in `A`
an equivalence, not merely a bound.

There is also a useful exact local count.  Group candidates by their addressed
outgoing tooth and component, and let `G` be the set of nonempty groups.  At
most one proposed successor can win in a group.  Hence

```text
|E|-|S| = (|E|-|G|) + (|G|-|S|).                     (16)
```

The first term counts excess candidates sharing one frontier.  The second
counts outgoing frontiers whose true winner is not itself a located candidate
in that group.  This is a genuine aggregate third-tooth lemma, though its
right side can still be large.

## 4. What actually happens in the two `h=420` controls

The exact sieve is:

| extra tail | all candidates `|E|` | renewal candidates `|O|` | successor survivors in `O` | actual edges `|A|` | groups `|G|` | all successor survivors `|S|` |
|---:|---:|---:|---:|---:|---:|---:|
| `1287` | 21,176 | 50 | 10 | 6 | 15,409 | 10,239 |
| `9009` | 27,404 | 86 | 11 | 6 | 20,340 | 12,458 |

Both rows have the same component classification:

```text
uncovered 726,              one-tooth-spanned 110,
renewal 4.                                             (17)
```

For `P=1287`, the local loss in `(16)` is

```text
21,176-10,239 = 5,767 + 5,170.                       (18)
```

For `P=9009`, it is

```text
27,404-12,458 = 7,064 + 7,882.                       (19)
```

The pair `{525,1365}` is even more diagnostic.  In each orientation its 105
candidates split as

```text
90 on uncovered components,
14 on one-tooth-spanned components,
 1 on a renewal component.                           (20)
```

That single renewal candidate survives the successor test and is the actual
edge.  Thus, across both orientations,

```text
210 located candidates -> 2 renewal obligations -> 2 actual uses. (21)
```

Local third teeth kill only 39 of the 210 candidates when `P=1287` and 33
when `P=9009`; none kills either actual obligation.  The earlier description
“competition is doing most of the work” is therefore too coarse.  The main
loss in this hostile is the **obligation filter**, not successor domination.
For the complete banks, however, the nontrivial `50 -> 10 -> 6` and
`86 -> 11 -> 6` steps show that both successor and upstream reachability are
real after the obligation filter.

## 5. Nested multiples are invisible third teeth

There is a sharp obstruction to any unconditional competition theorem.

> **Nested-multiple preservation lemma.**  Suppose a closed interval `I` has
> a two-tooth farthest-reach cover by teeth of speeds `u` and `v`.  Adjoin any
> collection of distinct speeds
>
> ```text
> c_i u,                 c_i>=2,             14 does not divide c_i. (22)
> ```
>
> Then the original two-tooth greedy chain is unchanged, whether the
> `u`-tooth is first or second.

At either boundary of an `u`-tooth,

```text
t=(14n+-1)/(14u),
||(c_i u)t||=||c_i/14||>=1/14.                       (23)
```

Thus no `c_i u` danger tooth crosses an `u`-tooth boundary.  Any such tooth
active inside the `u`-tooth ends no later than its parent; an equal right wall
is won by the smaller speed `u`.  If the `u`-tooth is first, this proves it
still wins at `L`, while every multiple is inactive at its outgoing wall.  If
the `u`-tooth is second, every multiple active before its left wall expires
before that wall and cannot beat the first tooth; once inside the `u`-tooth it
cannot beat its parent.  This proves both orientations and any number of
adjoined multiples.

The condition `14 does not divide c` is the exact boundary-crossing condition
in the general interval statement.  Outside the odd-tail branch it cannot be
dropped: at

```text
h=26, k=35, base speeds (3,29),
(3,2)->(29,20)
```

is the two-tooth chain, while adjoining `42=14*3` changes it to

```text
(3,2)->(42,29).                                      (24)
```

Inside the minority branch every multiplier preserving oddness is odd, so it
automatically satisfies `(22)`.

## 6. A full twelve-tail, denominator-complete, half-turn-hostile family

For every `M>=1`, put

```text
h=420M,
u=140M+1,                 v=140M+3,
C=(9,11,13,39,41,43,45,47,49,51),
W_M={u,v} union {cu:c in C}.                          (25)
```

Then `{2h} union W_M` is a primitive thirteen-speed minority row:

- all twelve tails are distinct and odd;
- `gcd(u,v)=1`, hence the complete row is primitive;
- `420|h`, while `9u,11u,13u` pay the additional `9,11,13` divisibility
  necessities in THM-4330;
- the five tails `cu`, `c=39,41,43,45,47`, satisfy

  ```text
  12h<cu<16h,                                         (26)
  ```

  so both inherited half-turn clocks fail.

Take `N=10M` in the primitive coprime family of the address-capacity audit.
For every determinant `d=N+1,...,2N`, each orientation of `{u,v}` is a
literal two-tooth greedy cover.  The nested-multiple lemma preserves all of
them after the other ten tails in `(25)` are adjoined.  Therefore the same
coprime pair occurs at least

```text
2N=20M=h/21                                           (27)
```

times in actual greedy renewal chains of the full twelve-tail row.

This family is not a candidate counterexample.  It has the uniform exact
wall witness

```text
t_M=1/(14u),                 clearance=1/14.          (28)
```

Indeed `u` is tight, every `cu` is safe by `(23)`,
`vt_M=1/14+1/(7u)`, and

```text
||(2h)t_M||=60M/(140M+1)>=1/14.                      (29)
```

The exact runs at `M=1,2,5` retain respectively at least `24,50,122` base
pair edges, while the full rows have `26,50,124` uses of the same pair.
Thus even all twelve odd tails, primitive normalization, `420|h`, the
remaining small-denominator necessities, and failure of both half-turn clocks
do **not** make third-tooth competition bound speed-pair reuse by a constant.

The mechanism also gives a positive replacement.

> **Wall-shadow reduction.**  At every `u`-wall `(14n+-1)/(14u)`, all speeds
> `cu` with `14` not dividing `c` are automatically safe.  Consequently, if
> a minority row is split as
>
> ```text
> W={u} union {c_i u} union R,                         (30)
> ```
>
> then strict failure forces the residual nonmultiple bank `R` to cover every
> `u`-wall which is safe for the anchor `2h`.

This converts a large nested class from useless “extra competitors” into a
small exact boundary-cover obligation.  It is the natural sidecar missing
from a competition-only count.

## 7. Exact wall capacity and the sharp two-residual hostile

The boundary obligation can be counted exactly.  For a positive speed `s`,
write

```text
d_s=gcd(u,s),              q_s=u/d_s,              S_s=s/d_s,
C_q(a)=#{y in Z:-q<y<q and y=a mod 14}.             (31)
```

Then the number of signed points of `X_u` which are strictly dangerous for
`s` is

```text
N_u(s)=d_s[C_(q_s)(S_s)+C_(q_s)(-S_s)].             (32)
```

For a fixed wall sign, danger is equivalent to the existence of an integer
`y` with

```text
y=14S_s n+S_s sigma-14q_s m,              |y|<q_s. (33)
```

The congruence class of `y` determines `n mod q_s` uniquely; each such class
has `d_s` lifts modulo `u`.  This proves `(32)`, including all strict endpoint
semantics.  Equivalently, each sign gives a translated uniform `q_s`-point
grid, so

```text
N_u(s)<=B_u(s):=2d_s ceil(q_s/7).                    (34)
```

Strict failure of a row

```text
{2h,u} union {cu:14 does not divide c} union R
```

therefore requires the aggregate capacity inequality

```text
2u<=B_u(2h)+sum_(r in R)B_u(r).                      (35)
```

For odd `u`, every proper quotient `q_s>1` is odd and at least three, and

```text
N_u(s)<=2u/3,                                        (36)
```

with equality only at `q_s=3`.  For odd `s`, equality further requires
`S_s=+-1 mod 14`.  This already proves a clean one-residual result:

> If `h,u,r` are positive integers, `u,r` are odd, `u` does not divide `r`,
> and `7u` does not divide `h`, then for every finite multiplier set with
> `14` dividing none of its elements, the row
> `{2h,u,r} union {cu}` is safe at some point of `X_u`.

If `u` does not divide `h`, anchor danger and residual danger each have size
at most `2u/3`; if `u|h` but `7u` does not divide `h`, the anchor is safe on
all of `X_u`.  The excluded resonance is exact: when `7u|h`, the anchor is
dangerous on every `u`-wall.

Two residuals can cover, but only in one rigid equality pattern.  Let
`u,r_1,r_2` be positive odd integers, suppose `u` divides neither residual,
and retain `7u` not dividing `h`.  The anchor and the two residual danger
masks cover all `2u` walls if and only if, after exchanging the residuals,

```text
u=3d,                 h=dH,                 r_i=dR_i,
H=+-7 mod 21,
R_1=+-1 mod 42,       R_2=+-13 mod 42,              (37)
```

for an odd positive `d`.  Equality in the three copies of `(36)` forces all
three quotient denominators to be three and their masks to be disjoint.  On
the quotient wall set

```text
{+-1,+-13,+-15} mod 42,                              (38)
```

the two odd residual classes in `(37)` kill the first two signed pairs, while
the anchor class kills the last.  This proves necessity and sufficiency.  The
least positive residual representatives are `(1,13)`.  If both residuals are
required to exceed the anchor `2h`, the least representative hostile is

```text
(h,u,r_1,r_2)=(7,3,29,41).                           (39)
```

It is not an unsafe LRC row; it only shows that a universal two-residual wall
closure would be false.

The gcd concentration is useful on the target branch.  In `(37)`, the anchor,
`u`, both residuals, and every nested multiple all share `d`.  If the complete
row is primitive, then `d=1`; if also `420|h`, this contradicts
`H=+-7 mod 21`.  Thus, away from the exact resonance `7u|h`, a primitive
`420|h` row with at most two residual nonmultiples is safe on `X_u`.

The exact program checks `(32)--(36)` on 14,336 `(u,s)` pairs through odd
`u=63`, including the sharp capacity example `(u,s)=(31,433)`, and audits
the quotient-three classification on 20,384 scaled rows.  Exactly 128 of
those rows cover, all and only the class `(37)`.

## 8. Connection contract and next sharp problem

```text
source:       addressed tail teeth on one anchor component
target:       a prefix-record envelope and its frontier orbit
map:          tooth -> (left wall,right wall); iterate the maximum begun reach
preserved:    exact greedy choice, strict activity, ties, component address,
              successor and reachability
destroyed:    the speed-pair quotient loses component status and frontier;
              the envelope forgets alternative nonwinning covers
sidecar:      renewal/spanning/missing status, visibility cap, centered residue
decisive test: E -> O -> S intersect O -> A
```

The unexpected useful analogy is interval scheduling, not a tournament: right
walls form a transitive prefix order, and manufacturing pairwise orientations
would add no target information.

The next sharp route is now a narrower dichotomy.

1. If many tails are **nonnested**, use `(13)` and the exact visibility caps to
   bound successor and upstream reach slots on genuine renewal components.
2. If many tails lie in one multiplicative nesting class, discard them on the
   corresponding `u`-wall set.  Banks of at most two residual nonmultiples are
   now closed in the primitive `420|h` branch whenever `7u` does not divide
   `h`.  The remaining sharp cases are the total anchor resonance `7u|h` and
   banks of at least three genuinely nonnested residuals.

The full-cover hypothesis must enter one of these two sides.  The family
`(25)` proves that third-tooth competition alone, even with every obvious
minority normalization imposed, cannot close the branch.  The wall theorem
does prove a residual bound through size two away from resonance, but no
counterexample-conditional competition theorem for three or more residuals,
minority closure, or LRC(14) proof is claimed here.

## Reproduction

From the repository root:

```bash
python3 -B 04-computation/lrc14_third_tooth_competition_probe_third_tooth_20260902.py \
  | diff -u 05-knowledge/results/lrc14_third_tooth_competition_probe_third_tooth_20260902.out -
```

The normal, optimized, and hash-seeded replays are required to match the
frozen output byte-for-byte.

Raw-LF SHA-256:

```text
28d170e65faee92103baf279d91b4b7705e10be9223a46b8a47eea7a48486540  script
80ea67fb849fee0c36ccef5b75283d3c4caf90218ab920c66f11113513b43b39  output
```
