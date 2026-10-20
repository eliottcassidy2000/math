# Discarded reset debt composes, with its multiplier gauge and source domain

2026-10-20. **PROVED scoped elementary identities and implications;
FINITE-EXACT controls.** These are research propositions, not new canon IDs.
Universal Collatz entry and first-exit positive-charge payment remain
**OPEN**. No literature-priority claim is made.

## Inheritance and the coordinate being composed

The inherited mechanism is the exact positive-charge decomposition in
[adaptive budgets, Connections §3](../adaptive_budgets_smooth_interfaces_20261020/CONNECTIONS.md#3-a-positive-charge-collatz-reroute-and-what-it-has-not-proved),
continued in [first-exit charge, §§1–2](../first_exit_complement_games_20261020/README.md).
The canonical hostile is that note's authenticated later descent from55
to53, which is not positive-charge paid. The corrected near miss is
discarding a nonnegative reset dividend from an already descending word.
The least-used useful sidecar here is the exact source pullback, inherited
from [THM-4512, coefficient-descent classes](../../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md)
and [the boundary compiler §1](../collatz_boundary_compiler_20261004.md#1-exact-and-coarse-words-including-the-signed-boundary).

This is the assigned session's Collatz niche. Its local concept board is
**ordered word, forward multiplier, discarded reset debt, original source,
first-exit gap**. The composition law links the first three; pullback and
payment keep the latter two. It is a different debt from the signed-row
conformal residual in [the Oct15 reset-debt packet](../collatz_reset_debt_20261015f.md).

**Concurrent inheritance, integrated before close-out:** incoming commit
`70cb917742` supplies [the forward debt cone and terminal-payment theorem](../charge_palettes_boundary_20261020g/charge.md).
Its forward composition law overlaps R1 below and proves the stronger
family `D<=c(a-1)`. Section3a records its terminal-payment consequence with
an explicit zero-budget boundary. The reciprocal gauge, exact pullbacks,
mixed-cone separation, and source-gap composition here are complementary
clarifications. The overlapping forward result is not claimed as new here.

Use the odd map `U(n)=oddpart(3n+1)` on positive odd integers. A word
`w=(a_1,...,a_m)` has positive integer letters specifying exact valuations.
Write `K_i=a_1+...+a_i`, `K=K_m`, `A=3^m`, `Q=2^K`, and

    B_w = sum_(i=1)^m 3^(m-i) 2^K_(i-1).

In the coordinate `z=n-1`, define

    a_w = A/Q,
    P_w = sum_(a_i=1) 3^(m-i)/2^(K-K_i),
    D_w = sum_(a_i>=3) (1-2^(2-a_i)) 3^(m-i)/2^(K-K_i).

Thus `P_w,D_w>=0`, and the exact endpoint and its positive-charge bound are

    F_w(n)-1 = a_w(n-1)+P_w-D_w,
    E_w(z)   = a_w z+P_w.                               (1)

This follows by composing the individual shifted steps
`z'=(3/2^a)z+2^(2-a)-1`: letter1 emits+1, letter2 emits0, and letters
at least3 emit a negative reset dividend. Empty words have `(a,P,D)=(1,0,0)`.

## 1. Exact composition and the two gauges

**Proposition R1 (PROVED).** Concatenate u followed chronologically by v,
denoted uv. Then

    a_uv = a_v a_u,
    P_uv = a_v P_u+P_v,
    D_uv = a_v D_u+D_v.                                 (2)

Every reset already in u passes through the full slope of v. Substitution
in (1), or splitting the defining suffix sum at the boundary, proves (2).
The identity word supplies a unit and associativity follows from affine
composition. In general the debt product is not commutative.

**Forward gauge.** If the user's multiplier is `lambda=a_w` and the debt
is `delta=D_w`, the predicate

    0<=delta<=lambda-1                                  (3)

is closed under chronological composition. Indeed, with
`s_w=a_w-1-D_w`,

    s_uv = a_v s_u+s_v.                                 (4)

Positive multipliers preserve nonnegative slack. For a sequence of blocks
the total debt is `sum_j D_j product_(k>j) a_k`, bounded by
`product_j a_j-1` when each block satisfies (3).

This gauge is **expansion-only**: (3) forces `a_w>=1`. For a nonempty
Collatz word, `3^m!=2^K`, so it forces `a_w>1`. In particular
`F_w(n)=a_w n+B_w/Q>n` for every positive legal source. It cannot itself
certify a descending block.

**Reciprocal, source-normalized gauge.** If instead

    lambda=1/a_w,       delta=D_w/a_w,

then the identically written condition `delta<=lambda-1` means

    D_w<=1-a_w.                                        (5)

Its exact law is different:

    lambda_uv=lambda_u lambda_v,
    delta_uv=delta_u+lambda_u delta_v.                  (6)

It is again composition-stable. Its slack obeys
`(lambda_uv-1-delta_uv)=(lambda_u-1-delta_u)
 +lambda_u(lambda_v-1-delta_v)`.
In forward coordinates the same fact is

    1-a_uv-D_uv = a_v(1-a_u-D_u)+(1-a_v-D_v).             (7)

This second cone is contraction-only for nonempty words. Changing a to1/a
without also changing D to D/a does not give (6).

**Equality boundary.** Both cones include the empty word on their boundary.
Neither has a nonempty Collatz word on its boundary. In fact

    Q D_w = sum_(a_i>=3) (2^a_i-4) 3^(m-i) 2^K_(i-1)

is divisible by4, whereas `A-Q` and `Q-A` are odd. Thus either weak
inequality is automatically strict for a nonempty word in its cone.
These integrality facts are additional to the general affine proof.

## 2. Source domains compose by pullback, not by forgetting guards

**Proposition R2 (PROVED).** The exact source domain is

    C_w={n>0 odd: n=(Q-B_w)A^(-1) mod2Q}.               (8)

On this domain every intermediate value is positive odd and the word
is literal. Endpoint oddness and backward induction enforce each exact
division; this is the inherited exact-cylinder dictionary. In particular
the modulus is2Q, notQ.

The exact composition domain is

    C_uv = C_u intersection F_u^(-1)(C_v).              (9)

More explicitly, take the least positive representative r_u and put
`s_u=F_u(r_u)`. Every source of u is `n=r_u+2Q_u t`, `t>=0`, with endpoint
`s_u+2A_u t`. If r_v represents C_v, its next guard is

    A_u t = (r_v-s_u)/2 mod Q_v.                       (10)

Both representatives are odd and A_u is invertible modulo Q_v. Hence
there is one residue for t modulo Q_v, giving the single exact cylinder
of modulus `2Q_u Q_v`. All finite positive words have nonempty positive
cylinders, but arbitrary subsets chosen by a controller can have empty
pullbacks. No all-source coverage follows from algebraic composition.

Extra source conditions must also be pulled back. For example an actual
first exit requires, for the *original* n,

    F_(a_1,...,a_j)(n)>=n for every j<m,
    F_w(n)<n.                                         (11)

These are explicit affine inequalities on the same residue class. First
exit is not preserved merely by concatenating two words that are first
exits from their separate starting values.

For source-dependent multipliers and debts, formulas (2)–(7) are pointwise:
evaluate the second block at `F_u(n)`. The second bound must hold there,
not only on an independently selected set of source values.

## 3. Hostiles: order, mixed gauges, and genuine source payment

The first short expansion hostile is

    79 --1,1,1,3--> 101,
    a=81/64, D=1/2>17/64=a-1.

All coefficient prefixes have slope greater than1, so even that stronger
prefix condition does not pay discarded-reset debt. A word containing a
reset cannot have forward slope greater than1 at length at most3, since
its total valuation is at least m+2. This is therefore a minimum-length
hostile among expanding words containing a reset.

Order remains substantive even at fixed length and clock:

| Word | Least exact source | Endpoint | a | D | D<=a-1 |
|---|---:|---:|---:|---:|---|
|11113|223|425|243/128|1/2|yes|
|31111|253|485|243/128|81/32|no|

Moving the reset before the growth letters amplifies its discarded debt.
The different exact cylinders must move with the words; this comparison
does not reorder either source's physical trajectory.

The union of the two stable cones is **not** composition-stable. A letter1
belongs to the forward cone; a letter3 belongs to the reciprocal cone.
Their chronological concatenation13 has

    a=9/16, D=1/2>1-a=7/16,
    19 --1,3--> 11.                                    (12)

It lies in neither cone, although this native first exit is paid:
`D=1/2<19-11=8`. Thus the reciprocal condition is sufficient for a debt
bound, not necessary for positive-charge payment. A strategy alternating
the two cone types needs a separate budget at their joins.

The actual payment test, for a legal block with source n and endpoint h, is

    E_w(n-1)<n-1  iff  D_w<n-h.                        (13)

It is stronger than actual descent. If the reciprocal cone holds and an
actual strict descent is already supplied, then `D_w<=1-a_w<1<n-h`
(oddness gives `n-h>=2`), so payment follows. This does not prove that
every source enters that cone, encounters a coefficient crossing, or has
an actual first exit. The source threshold remains
`n>1+P_w/(1-a_w)` for positive-charge payment when `a_w<1`.

There is a separate useful composition rule **at authenticated descending
sources**, without either uniform cone. Suppose u takes n to h<n, v takes
h to k<h, and both local bounds are paid. Put `g_u=n-h`, `g_v=h-k`.
Then `D_u<g_u`, `D_v<g_v`, and `0<a_v<1` (the unshifted carry is positive).
Therefore

    D_uv=a_v D_u+D_v < g_u+g_v=n-k.                    (14)

So the composite positive-charge bound is paid. More precisely its
debt/gap ratio is at most the larger local ratio, since it equals
`[a_v g_u(D_u/g_u)+g_v(D_v/g_v)]/(g_u+g_v)`.
This cannot prove the unrestricted first-exit hypothesis by subdividing a
first exit into earlier descending pieces: those pieces do not exist
before its first descent. The retained source/history is the obstruction.

## 3a. The forward cone pays a terminal descent

The [concurrent packet](../charge_palettes_boundary_20261020g/charge.md)
proves that, for each fixed `c>=0`, the forward cone
`a>=1, 0<=D<=c(a-1)` is stable by (2). This follows from exactly the
telescoping calculation behind (4). Its useful consumer is a terminal
descent, rather than descent of the expanding prefix itself.

**Proposition R3 (PROVED; concurrent theorem with boundary made explicit).**
Let an authenticated prefix u satisfy `a_u>=1` and
`D_u<=c(a_u-1)` for `0<=c<=1`. If one further actual step of valuation
b gives endpoint h<n from positive odd n>1, its whole positive-charge
bound pays. The same holds for a nonempty terminal packet v whenever
`mu=a_v<1`, `c(1-mu)+D_v<=2`, and its actual endpoint h is below n.

**Proof for one final step.** The whole unshifted word has positive carry,
so `mu*a_u<h/n<1`, with `mu=3/2^b`. Since `a_u>=1`, this forces b>=2.
For c>0, composition gives

    D_total<=c*mu*(a_u-1)+D_v<c(1-mu)+D_v.             (15)

If b=2 then `D_v=0`, so this bound is at most1/4. If b>=3 then
`D_v=1-4/2^b` and

    D_total<c(1-mu)+1-4/2^b<=2-7/2^b<2.               (16)

If c=0, instead `D_u=0` and `D_total=D_v`: this is0 at b=2 and less
than1 at b>=3. Thus in every case `D_total<2<=n-h`, which is (13).
No first-exit premise is needed. No existence of the final valuation or
of a descending endpoint is inferred.

**General terminal packet and the c=0 boundary.** For c>0 the same positive
unshifted carry gives (15), hence the displayed budget proves
`D_total<2<=n-h`. At c=0, the inequality in (15) need not be strict;
we have `D_total=D_v<=2`. Every positive literal reset debt is a
nonintegral dyadic rational, so equality at2 is impossible. To see this,
let i be the first reset index of v. In the numerator expression for
`Q_v D_v` above, the term for index j has two-adic valuation exactly
`K_(j-1)+2`, because `2^a_j-4=4(2^(a_j-2)-1)` with the bracket odd.
These valuations strictly increase with reset index j. The earliest term
is the unique minimum, so

    v_2(Q_v D_v)=K_(i-1)+2<K_v.                        (17)

Consequently the reduced denominator of D_v is
`2^(K_v-K_(i-1)-2)>=2`. If there is no reset, D_v=0. In either case
`D_v<=2` implies `D_v<2`. This proves the terminal theorem also at c=0.

The initial incoming version70cb917742 used the unqualified strict
inequality (15), which fails at c=0, for example for the empty prefix
followed by letter2 at source9. Both sides there are0; the payment
conclusion survives. Before our integration, incoming8dcc679ae3 independently
made this same dyadic-debt repair and logged its mechanism. The current
source and the proof above agree. The repair does not strengthen the
theorem's source-coverage claim.

The concurrent packet also gives an infinite nonvacuity family:
`b=(1,1,1,1,3)` has `a_b=243/128` and debt ratio
`D_b/(a_b-1)=64/115`. Each repeated prefix `b^k` stays in that cone.
Choosing a final sufficiently large valuation gives an exact cylinder
with infinitely many sources whose first descent occurs at its end.
This permits arbitrarily many resets. It does not cover every source;
the same packet retains63 as an unresolved prefix for this cone, even
though its first exit is paid by the direct test (13).

## 4. Reproduction, connection contract, and stopping point

[The exact checker](reset_composition_check.py) compares recurrence data
with independent suffix sums; checks every cut of words of length0..5 over
letters1..4; replays their first three positive cylinder representatives;
checks both cone laws, empty-word equality, and the stated hostiles; and
checks (14) on two consecutive first descents from each odd source3..1001
whose first endpoint is greater than1. Both descents must terminate within
10000 odd steps; hitting the cap is an explicit error, not a false return.
That finite source universe supplies controls only, not a termination proof.
The program uses explicit exceptions and produces the same receipt with-O.

The [saved receipt](reset_composition_receipt.json) records1365 words,
7737 cuts,4095 native cylinder replays, and492 paid descending pairs;
eight first endpoints equal1 and are excluded from the second descent.
Normal and optimized stdout agree byte-for-byte. SHA256:
`73b2ea2dc51c508b6b86cba5df83978ed8282902e4b40e2a4ffbc2d5d229db24`
for the checker and
`b3f0bf64c1b9344f1ff7fcb0b102a6788a8f12e7dd864ff46ec0253c09f8dd80`
for the LF-normalized saved receipt.

    python 05-knowledge/results/trijection_reset_20261020/reset_composition_check.py
    python -O 05-knowledge/results/trijection_reset_20261020/reset_composition_check.py

Connection contract: the source is an authenticated ordered valuation word;
the target is `(a,P,D)` with its exact cylinder; the map takes suffix-weighted
emissions. It preserves the affine endpoint, its positive upper bound, and
the two stated cone predicates under lawful composition. Forgetting the
word loses internal history and first-exit status; forgetting the cylinder
loses the actual integer orbit. Retain the word, original source, and
prefix inequalities when those are the consumer. The shortest decisive
tests are1113 for automatic forward payment and13 for mixed-cone closure.

The structural question about composition is closed in these gauges.
The remaining research question is a source-covering selection of lawful
blocks whose reset debt is bounded by the original descent margin. No
larger word census is needed to establish the composition identities.
