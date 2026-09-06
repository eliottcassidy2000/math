# THM-4448 result note: general shore attachment and decoder-pair cones

**Status: PROVED ELEMENTARY + PROVED RELATIVE TO CITED `LRCUpTo13` +
FINITE-EXACT + INDEPENDENTLY AUDITED.**  The attachment lemma and hostile
families are proved arguments.  The tail cells, pair maxima, and bounded body
census are finite-exact.  This closes a cofinal scale cone, not the residual
opposite-scale wedge; `LRC(14)` remains **OPEN**.

## General attachment statement

Write a ten-body as

```text
C=B union hU,   |U|=r,   |B|=10-r,   1<=r<=9,   h in Z_(>=1),
```

where `B` is nonempty and the displayed body speeds are distinct.  Let `T`
be any three distinct positive tails; tails divisible by three are allowed.
Put

```text
D_U={z: ||uz||<1/14 for some u in U}
```

and assume `G_U` is nonempty.  If `delta(U)` is the maximum length of an open
component of `D_U`, apply cited lower-dimensional LRC to the at most `13-r`
moving speeds in `3B union T`.  For a resulting witness `x_0`, set `y_0=3x_0`
and

```text
rho_*(x_0)=min(
  min_(b in B) (||3b x_0||-1/14)/b,
  min_(t in T) 3(||t x_0||-1/14)/t).
```

The closed quotient arc of radius `rho_*` about `y_0` is body-safe and carries
the unique local inverse branch through `x_0`, which remains safe for all
tails.  Moreover

```text
rho_*(x_0) >= rho_r
rho_r=min(r/[14(14-r)max(B)], 3r/[14(14-r)max(T)]).
```

Every open component of `D_U` of length `ell<=1` pulls back under `y -> hy`
to components of length `ell/h`.  Consequently

```text
2h rho_*(x_0) >= delta(U)
```

completes `3(B union hU) union T`; the uniform condition
`2h rho_r>=delta(U)` suffices.  Equality is included: a closed arc cannot be
contained in an equal-length open arc because positive endpoint gaps would be
required.  The proof uses no inverse of a tail modulo three.

For `r=1`, after naming the single physical shore speed `v=hu`,
this recovers the familiar `v>=13 max(B)`,
`3v>=13 max(T)` scale, included as a prior-mechanism sanity check rather than
a novelty claim.

## Exact two-speed constants

For coprime `p<q`, let `delta(p,q)` be the longest component of the strict
danger union for `{p,q}`.  At `r=2`, with `M=max(B)` and `E=max(T)`, the
uniform radius is

```text
rho=min(1/(84M),1/(28E)).
```

The two exact implementations give:

```text
unfiltered p+q<=356:
  19,314 pairs, max delta=15/98 uniquely at (1,14),
  cone 7h>=45M and 7h>=15E;

THM-3818 decoder filter:
  5,855 pairs, max delta=29/196 uniquely at (1,28),
  cone 14h>=87M and 14h>=29E.
```

The decoder filter requires every prime divisor of `p+q` to be `2 mod 3`
with exponent at most two.  Thus `29/196` is not an all-pair bound: `(1,14)`
is the sharp unfiltered countercontrol.  Strict-open merging is also
load-bearing: `(1,13)` has danger gap `1/7`; joining its teeth through their
safe equality singleton produces the false value `15/91`.

THM-3818's separate number `46,837` counts positive-scale divisor-seam
triples `(p,q,s)` with `s>=2`, not ratios.  Together with `5,855` scale-one
triples it gives `52,692` residual seams.

## Critical-tail addresses

For the two ternary-unit tails, the complete open failure cells and owner
vectors are

```text
T=(1,5,11):
  (25/154,31/154)   owners (0,1,2)
  (123/154,129/154) owners (2,1,0)

T=(2,11,20):
  (5/56,3/28)       owners (0,1,2)
  (123/280,129/280) owners (1,2,0)
  (151/280,157/280) owners (1,0,2)
  (25/28,51/56)     owners (2,1,0).
```

Their masses are `6/77` and `11/140`. If `C` is nonempty and `J=[L,R]` is
any closed component of its nonempty safe set `G_C`, including a singleton,
then `J` is wholly
spoiled precisely when there are integers `n_t` such that

```text
n_t/t-3/(14t) < L <= R < n_t/t+3/(14t)
{-n_t t^(-1) mod 3 : t in T}=Z/3Z.
```

Thus `G_(3C union T)` is empty iff every component of `G_C` satisfies this
strict addressed containment.  An equality endpoint immediately frees a
sheet.

The primary census finds no containment for either tail among all `43,758`
ten-subsets of `[18]`.  The cumulative body counts at heights
`13,14,16,18` are `286,1001,8008,43758`; scalar-low counts are
`12,12,15,44` for `(1,5,11)` and `13,13,18,47` for `(2,11,20)`.  This is
**FINITE-EXACT** only.  The height-13 column recovers THM-4442 as a positive
control, and all 572 tail/body cases there also agree with literal physical
safe-set construction.

Full-overlap controls show why scalar mass is insufficient:

```text
T=(1,5,11),  C=(1,2,3,4,7,8,9,13,14,19), overlap=6/77;
T=(2,11,20), C=(1,2,3,4,5,6,7,8,12,14), overlap=11/140.
```

Each body contains the entire failure set but escapes on other components.

## Cofinal prescribed-component obstruction

Let `A={1,...,8}` and `C_N=A union {N,4N}`.  The progressions

```text
T=(1,5,11):  N=53+2310k, y_*=2/11;
T=(2,11,20): N=121+210k, y_*=1/10,
```

for `k in Z_(>=0)` are primitive and have distinguished cross height exactly `N`.
The component through `y_*` is contained in one speed-`N` safe tooth of
length `6/(7N)`.  Since

```text
6/(7*53)<3/154,   6/(7*121)<1/140,
```

the marked component lies strictly inside the displayed tail-failure cell.
Nevertheless the cofinal decoder cone completes every row in both families.
Therefore arbitrarily large cross height cannot certify an arbitrary or
preselected component; a valid proof must select a component existentially
or retain its inverse-branch/address word.

## Reproduction and hashes

From the repository root:

```powershell
python -B 04-computation/lrc14_general_shore_attachment_thm4448.py
python -B -O 04-computation/lrc14_general_shore_attachment_thm4448.py
python -B 04-computation/lrc14_general_shore_attachment_thm4448_independent.py
python -B -O 04-computation/lrc14_general_shore_attachment_thm4448_independent.py
```

Normal and optimized executions are line-identical to their respective
frozen outputs.  Hashes are SHA-256 of raw LF repository bytes:

```text
04-computation/lrc14_general_shore_attachment_thm4448.py
  f2070a2b77d17b2635c73efe53cbddf450e94fc71161a71c2fbc12a55d41488b
05-knowledge/results/lrc14_general_shore_attachment_thm4448.out
  9dfabe053f72d88a3115d2dff171ab085527c9ca21d001cad19c9325737b9667
04-computation/lrc14_general_shore_attachment_thm4448_independent.py
  989ed36ca551e0e31d29071579e13e4a5ba696d2ea3776f8a61d4683f374fb26
05-knowledge/results/lrc14_general_shore_attachment_thm4448_independent.out
  b02e4827a56d453df05836bc8d6ed13eacfb70253d0ffa742d5d96bc0aa5f916
05-knowledge/results/lrc14_general_shore_attachment_thm4448_independent_audit.md
  817d4835492b3b1b8cff7c8580ad846ce208cbfcc499b467a94559378464f8ca
```
