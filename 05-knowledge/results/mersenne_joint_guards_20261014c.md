# Mersenne joint guards: three short repairs and authenticated child terminals

**Status.** PROVED: the run obstruction, complete short-head classification in
the stated template, native exponent phases, all-height paid joins, and the
finite-anchor rank obstruction. FINITE-EXACT: the declared phase intersections
and supplied ROOT certificates. CONDITIONAL: completion of arbitrary members
of the infinite paid families requires their child certificates. Universal
Mersenne and Collatz coverage remain OPEN.

The three additional heads give exact exponent phases
`79 mod128`, `199 mod256`, and `323 mod512`. Each reduces the unchanged
Mersenne source by two exponent bits. Their three least children are supplied
with strict ROOT words, grounding nine selected Mersennes. A compatible
subphase of the third rule also admits an inherited four-bit continuation.
This is not a proof that its arbitrary emitted child reaches ROOT.

## 1. Inheritance and the missing joint guard

The closest mechanisms are
[THM-4555, uniform switches and collisions at minus one](../../01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md),
[THM-4601, corrected two-anchor reduction](../../01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md),
the [variable-depth reset-two rules](collatz_reset2_rules_20261007.md), and
the [general head phase compiler](collatz_general_head_phases_20261007b.md).
The [two-twos grounding packet](collatz_two_twos_grounding_20261007c.md)
keeps finite terminal proofs separate from all-height paid rules. The
[fibre integration](collatz_fibre_integration_20261009.md) retains the odd-count
and normalized-cost defects that deletion rules preserve.

The immediate hostile is the preceding
[paid block bank](collatz_successor_holonomy_20261014b.md): its words
`1^k 2^(2k)` pay every native source, but every nontrivial Mersenne misses
their joint guard. The new calculation makes that incompatibility structural
and then repairs a different, actual post-run branch.

Board: **unchanged source / actual run length / marked affine anchor /
exponent phase / child type / supplied first-hit terminal**. A head that works
at a formal checkpoint must still meet the source's exponent phase. An
available smaller child does not yet discharge its ROOT obligation.

Write `M_E=2^E-1`, and use the positive odd map
`U(n)=(3n+1)/2^{v_2(3n+1)}`, stopped at the first 1. Superscripts on word
letters mean repetition.

## 2. The first two runs cannot pay an odd Mersenne

For `E>=2`, the initial actual word is `1^(E-1)`, ending at

\[
x=2\cdot3^{E-1}-1,\qquad 3x+1=2(3^E-1).
\]

The extra factor two gives the next actual valuation

\[
a=1+v_2(3^E-1)=
\begin{cases}2,&E\text{ odd},\\3+v_2(E),&E\text{ even}.\end{cases}
\tag{1}
\]

For odd `E>=3`, put `k=E-1` and `m=v_2(k)`. Since
`v_2(x-1)=m+3`, the maximal following two-run has length

\[
J=1+\lfloor m/2\rfloor\le k/2. \tag{2}
\]

The intermediate bound `J<=2^(m-1)` is equality for `m=1,2`, while
`m+2<=2^m` suffices thereafter. Combining it with `k>=2^m` gives (2).
After `j<=J` twos, the slope of the entire
prefix is `3^(k+j)/2^(k+2j)`, decreasing with `j`. At `j=k/2` it is
`(sqrt(27)/4)^k>1`. Every nonempty positive valuation word has positive
carry. Thus **every prefix through the actual maximal two-run stays strictly
above the original Mersenne**. Merely prescribing more contracting twos
does not fix the source: the actual native guard has already ended that run.

In particular the previous bank's demand for `2k` twos is impossible here.
This is a run-specific obstruction, not a claim that later descents fail.

## 3. A complete short repair in the one-two shell

Now retain the actual shell `E=3 mod4`, whose maximal two-run has length one.
After the source's initial ones and its single two, and the even child's
initial ones and reset four, the checkpoints are

\[
X=\frac{3^E-1}{2},\qquad
Y=\frac{3^{E-1}-1}{8},\qquad X=12Y+1. \tag{3}
\]

The same `Y` is reached by `M_(E-2)` after `1^(E-3),(2,2)`.
These are actual positive odd checkpoints, because `v_2(E-1)=1`.

For a source head `u` of length `l`, cost `A`, and carrier
`F_u(x)=(Px+B)/Q`, seek a partner `v` satisfying

\[
|v|=l+1,\quad Q_v=Q/16,\quad
F_v(Y)=4F_u(12Y+1)+1. \tag{4}
\]

Its required carry is

\[
B_v=(P+B)/4+Q/16. \tag{5}
\]

We classify only **one- or two-letter heads beginning at least 3**, with
this fixed sibling gap and anchor. Other lengths, gaps and anchors are not
excluded.

For a one-letter head `(a)`, the partner must have two positive letters,
cost `a-4`, and carry `1+2^(a-4)`. A two-letter carry is `3+2^d`, so
`2^d=2^(a-4)-2`. Valuation and positivity force the unique solution
`u=(6),v=(1,1)`, the inherited sporadic rule.

For a two-letter head `(a,b)`, set `x=a-2>=1`, `y=a+b-4>=3`.
Equation (5) requires a three-letter carry

\[
3+2^x+2^y=9+2^d(3+2^e). \tag{6}
\]

* If `x=1`, valuation of the difference from 9 forces `d=2`; then
  `2^e=2^(y-2)-4`, forcing `y=5,e=2`.
* If `x=2`, it forces `d=1` and `2^e=2^(y-1)-4`, hence `y=4,e=2`.
* If `x,y>=3`, it forces `d=1,e=1`, and
  `2^(x-1)+2^(y-1)=8`, hence `x=y=3`.

The remaining partner letter is fixed by the total cost and is positive.
This proves completeness of the following table within the stated template:

| Source head `u` | Partner `v` | Exact exponent phase | Status relative to the declared baseline |
|---|---|---|---|
| `(6)` | `(1,1)` | `E=31 mod64` | inherited sporadic control |
| `(3,6)` | `(2,2,1)` | `E=323 mod512` | additional |
| `(4,4)` | `(1,2,1)` | `E=199 mod256` | additional |
| `(5,2)` | `(1,1,1)` | `E=79 mod128` | additional |

## 4. The phase is exact, and the source remains unchanged

The head's exact native odd residue is

\[
x_u=(Q-B)P^{-1}\pmod{2Q}.
\]

Substituting the actual source checkpoint (3), its native guard is exactly

\[
3^E\equiv2x_u+1\pmod{2^{A+2}}. \tag{7}
\]

Each target lies in the subgroup generated by 3; its order modulo
`2^(A+2)` is `2^A`. This gives the unique phases in the table, with the
displayed minimal periods. The compiler obtains the phase by a two-candidate
binary lift. The script independently tests every odd exponent in each
complete period. There is no replacement of a supplied exponent by the
nearest phase member.

For first-hit safety, the sufficient cutoff

\[
E-1\ge2A+2 \tag{8}
\]

ensures every source prefix through the head stays above `M_E`: a common
lower slope bound is

\[
\frac{3^{E}}{2^{E+1+A}}
\ge\frac{27}{16}\left(\frac98\right)^A>1.
\]

All four least representatives already satisfy (8), so **no positive member
of any displayed phase is discarded**. The cutoff is sufficient, not a
necessary payment criterion.

Put `Z=F_u(X)>1`, and take its actual next valuation `c=v_2(3Z+1)`.
Equation (4) makes the partner endpoint `4Z+1`, with next valuation `c+2`.
Consequently the actual source and child words are

\[
\begin{array}{c|l}
M_E&1^{E-1},2,u,c\\
M_{E-1}&1^{E-2},4,v,c+2\\
M_{E-2}&1^{E-3},2,2,v,c+2.
\end{array} \tag{9}
\]

All reach the same odd endpoint. Final oddness of the partner carrier forces
every earlier valuation to be exact: a dyadic noninteger cannot regain
integrality under an odd numerator multiplier, and an even integer would
make the next step nonintegral. The positive partner checkpoint and endpoint
also rule out earlier ROOT. The last common endpoint itself may be ROOT,
in which case the suffix is empty.

Thus both children are strictly smaller, and a supplied strict ROOT suffix
can be transported. The D2 path is exactly the composition of the D1 path
and the known even reset, with its literal common-future seam retained.
The two alternatives are factorizations, not independent guard coverage.

## 5. Coverage and child types

Declare the comparison bank to be the ordinary even-exponent reset and the
old sporadic reduced word `(2,6,c)`. ROOT at the finitely exceptional smallest
exponents can be checked separately. This is **not** a comparison against
every rule accumulated in the repository.

The new reduced source heads begin `(2,3,6)`, `(2,4,4)`, and `(2,5,2)`.
They miss both declared baseline branches. Their exponent phases are disjoint,
as also follows from their different first post-run valuations. Within the
shell `E=3 mod4`, the old phase has density `1/16`, the new phases have
densities `1/128,1/64,1/32`, and hence

\[
\text{additional density}=\frac7{128},\qquad
\text{union density}=\frac{15}{128}. \tag{10}
\]

These are densities of **exponents within a shell**, not densities of
integer Mersenne values or ROOT certificates. The residual shell fraction is
`113/128`; for example `E=327` matches none of these four heads. That is a
rule-bank miss, not a statement of nonconvergence.

Child types are useful for the next unchanged-source obligation. The D2
child of `E=323 mod512` has exponent `F=321 mod512`, so `v_2(F-1)=6` and
its maximal two-run has length four. The other three rows have
`v_2(F-1)=2`, hence two twos. Thus they must not all be handed to the same
post-run interface.

There is a genuine simultaneous continuation on the four-two child.
Intersecting `F=321 mod512` with the **frozen 225-head bank** from the
[general phase compiler](collatz_general_head_phases_20261007b.md) at `J=4`
gives exactly 24 disjoint phases, of conditional density `429/32768`.
The simplest is its old head `(1,10)` at `F=52545 mod131072`. Therefore

\[
E\equiv52547\pmod{131072}
\quad\Longrightarrow\quad
M_E\rightsquigarrow M_{E-2}\rightsquigarrow M_{E-6}. \tag{11}
\]

Both guards refer to the same integer source and its actual child.
The inherited growth cut is satisfied at every member. The script checks
all 225 phase intersections, pairwise disjointness, and the exact finite sum.
It does not expand the 52,547-bit source to establish (11). The terminal
`M_(E-6)` still requires ROOT grounding. The particular `F=321` hits none
of these 24 phases, which is why its independent supplied certificate matters.

## 6. Finite grounding and an invariant obstruction to infinite closure

The frozen [JSON](mersenne_joint_guards_20261014c.json) contains run-length
encoded, strict actual ROOT words at four selected children. Discovery was
explicitly bounded by 100,000 odd steps per input; production APIs only replay
the supplied data. Exact ranks, costs, peak bit lengths and decoded-word hashes
are checked.

| Supplied child exponent | Odd ROOT rank | Valuation cost | Peak bits |
|---|---:|---:|---:|
| 29 | 162 | 286 | 46 |
| 77 | 333 | 605 | 123 |
| 197 | 979 | 1749 | 312 |
| 321 | 1774 | 3133 | 509 |

Equation (9) grounds the corresponding triples
`29,30,31`; `77,78,79`; `197,198,199`; and `321,322,323`.
The two larger exponents in each triple have the same odd ROOT rank as its
supplied child and costs larger by one and two. Both factorizations export
the identical literal source ROOT word. These are twelve finite certificates,
nine associated with the additional phases; no arbitrary member is thereby
grounded.

The limitation has an exact invariant. A strict common-future receipt with
equal source and child word lengths preserves the total odd ROOT rank whenever
either endpoint is supplied with ROOT. Every receipt in (9) has zero odd-count
defect. But `M_E`, for `E>=2`, has `E-1` initial odd steps followed by at least
one further step, so its ROOT rank is at least `E`.

**Finite-anchor closure bound.** If all receipts in a grounding derivation
have zero odd-count defect, and all supplied terminals have odd rank at most
`R`, then any grounded Mersenne exponent `E>=2` satisfies `E<=R`.
This holds through arbitrary compositions and alternatives. For the four
supplied terminals here, `R=1774`. No amount of repeated zero-defect deletion
alone can turn this finite terminal set into ROOT proofs for unbounded
exponents. This is the finite-anchor form of the inherited fibre invariant,
not a claim against all possible common-future certificates. A route with
nonzero odd-count defect retains a different bill.

## 7. APIs and exact checks

`compile_rule` reconstructs the partner and exponent phase; `audit_rule`
rechecks canonical fields before equality. `guard` tests the supplied
exponent. `materialize_join` has an explicit size cap and produces the actual
two words. `discharge` verifies the child's full strict ROOT word and its
literal seam before exporting a source certificate. `root_word` reads only
the four supplied source-bound rows, with no discovery fallback.

`head_endpoint_residue` retains the full denominator `2Q` before modular
division. The symbolic controls include exponents with phase parameter
`10^30+7` and moduli divisible by 2 and 3. A bounded terminal-valuation query
returns `None` if its requested precision does not decide the valuation;
that result is not interpreted as a missing valuation or failed convergence.

Declared controls: all one-letter heads `3<=a<=64`; all two-letter heads
`3<=a<=64,1<=b<=64`; complete native exponent periods for the four rules;
actual initial runs for `2<=E<=256`; four literal members per phase; the four
terminal words and all twelve transported certificates; and all 225 inherited
head intersections. Hostiles include incorrect phase, wrong child terminal,
corrupted partner, boolean/float aliases, ROOT padding and materialization-cap
failure. Completeness of the short template is proved above, not inferred
from the finite head search.

```text
python -B 04-computation/experiments/mersenne_joint_guards_20261014c.py
python -B -O 04-computation/experiments/mersenne_joint_guards_20261014c.py
```

Both runs reproduce **6,161 checks** and the saved JSON, LF SHA256
`c3cb61b5a2d91b2262b1f69e214f187e91f0945c5006f90c0a4b7de743713ce8`.
