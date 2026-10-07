# Variable-depth rules inside the first-reset-two branch

Status: **PROVED** guarded common-future families and exponent phases;
**FINITE-EXACT** bounded collision search and three first-hit ROOT certificates.
Universal coverage by these rules is **OPEN** and is not asserted.

The named residual input \(2^{1459}-1\) is now fully grounded. A new uniform
two-bit deletion reduces it to \(2^{1457}-1\); an explicitly bounded search
for that one selected child supplies a frozen first-hit certificate. A deeper
six-bit deletion also transfers that certificate to \(2^{1451}-1\). Each
point has 7,347 odd steps to ROOT. The reusable rules apply to infinite
families, but their arbitrary smaller children remain proof obligations.

## 1. Inheritance and the information retained

The closest proved mechanism is
[THM-4555, uniform switches are collisions at minus one](../../01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md).
The preceding implementation is
[terminal lifts, sections 3–4](collatz_terminal_lifts_20261007.md): it removes
the old 128-letter cap for the ordinary reset and for the sporadic
\((2,6,c)\) collision. That note leaves \(2^{1459}-1\) as an explicit
remaining branch. We apply its existing collision theorem to new exact
cores; there is no priority claim for the general theorem or its inverse
grammar.

The live objects are the immutable positive source, an actual valuation word,
its affine carry, the smaller child, the exponent phase, and the authenticated
ROOT suffix. Equality of slopes or word lengths alone is insufficient. In
particular, the variable final valuation and the difference in core lengths
must accompany the collision. A finite search miss is not a nonconvergence
claim, and a successful smaller-child rule is not itself a ROOT proof.

Write \(U(n)=(3n+1)/2^{v_2(3n+1)}\), for positive odd \(n>1\), and stop at
the first 1. For a word \(w\), let

\[
 F_w(x)=\frac{P_wx+B_w}{Q_w},\qquad
 P_w=3^{|w|},\quad Q_w=2^{\sum w}.
\]

Composition is chronological. The exact native source cylinder is
\(x=(Q_w-B_w)P_w^{-1}\pmod{2Q_w}\). Its final oddness forces all its
specified valuations, by reversing the powers of two. Positivity of the
affine partial maps supplies positive intermediate states. First-hit ROOT
typing is an additional condition, retained below.

## 2. A common compiler for both depths

Let reduced heads \(u,v\) start with valuations at least two, with

\[
 |v|-|u|=D\ge1,\qquad Q_u=4Q_v,
 \qquad F_v(-1)=4F_u(-1)+1.                 \tag{1}
\]

Take a positive odd source \(n\), set \(r=v_2(n+1)-1\), and require
\(r\ge D\). Its actual initial word is \(1^r\), with checkpoint

\[
 x=3^r(n+1)/2^r-1.
\]

Require the native cylinder for \(u\), and put \(X=F_u(x)>1\). Then
\(c=v_2(3X+1)\ge1\). The smaller source

\[
 h=(n+1)/2^D-1<n                                      \tag{2}
\]

has the common future specified by

\[
 n:\quad 1^r,u,c,
 \qquad h:\quad 1^{r-D},v,c+2.                       \tag{3}
\]

Indeed, after the child's initial ones its checkpoint is
\((x+1)/3^D-1\). Since \(P_v=3^DP_u\), (1) gives

\[
 F_v((x+1)/3^D-1)=4F_u(x)+1=4X+1.
\]

This positive odd endpoint forces the child's intermediate integrality and
exact valuations. Finally \(3(4X+1)+1=4(3X+1)\), proving (3). The two heads
end above 1, so neither actual prefix can have passed through ROOT: continuing
an exact positive ROOT word would stay at 1 with valuation two. The final
common endpoint may equal 1 and then the suffix is empty. The guard \(X>1\)
avoids padding a source that had already reached ROOT in its head.

This is a source-authentic conditional implication
\(\operatorname{Root}(h)\Rightarrow\operatorname{Root}(n)\), with literal
certificate substitution. The converse also holds for these checked common
futures and can transport a supplied source certificate to the child. Neither
direction assumes an unsupplied trajectory terminates.

## 3. The short rule at exponent 1459

Use

\[
\begin{aligned}
 u&=(2,3,2,2,2,2,1,5,1,1,2),\\
 v_1&=(4,2,1,2,2,2,1,1,3,1,1,1),\\
 v_2&=(2,2,2,1,2,2,2,1,1,3,1,1,1).
\end{aligned}
\]

The first two carriers are

\[
 (P_u,Q_u,B_u)=(177147,8388608,12565973),
\quad (P_{v_1},Q_{v_1},B_{v_1})=(531441,2097152,15017419).
\]

Replacing the initial 4 of \(v_1\) by \((2,2)\) preserves its value at
\(-1\). Both partners therefore satisfy (1), with deletions \(D=1,2\).
For every \(c\ge1\), the completed cores have the common value

\[
 F_{u,c}(-1)=F_{v_i,c+2}(-1)
       =\frac{22777543}{2^{22+c}}.
\]

The exact native checkpoint guard is

\[
 x\equiv4662161\pmod{2^{24}}.                       \tag{4}
\]

For a Mersenne source \(n=2^K-1\), its initial run has length \(K-1\),
and \(x=2\cdot3^{K-1}-1\). Condition (4) is equivalent to

\[
 K\equiv1459\pmod{2^{21}}.                         \tag{5}
\]

This is an exact phase, not a sampled pattern. After division by two the
modulus is \(2^{23}\), where the order of 3 is \(2^{21}\). The elementary
proof is the even-exponent formula
\(v_2(3^d-1)=2+v_2(d)\), obtained by repeated squaring, together with the
odd-exponent obstruction modulo eight. Direct modular exponentiation verifies
the phase representative 1459. It is the least positive representative;
all later checkpoints increase and the seed head ends above ROOT.

Thus (5) gives \(2^K-1\rightsquigarrow2^{K-2}-1\) with the D2 receipt.
The source and child receipts each have \(K+11\) odd letters. The first
reduced letters are \((2,3)\): this whole phase misses both the ordinary
reset rule and the old sporadic \((2,6,c)\) rule.

## 4. A deeper rule at exponent 1457

The exact heads are stored as ordered tuples in
[the frozen data](collatz_reset2_rules_20261007.json). They have lengths
61 and 67, costs 128 and 126, and satisfy (1) with \(D=6\).
The source head starts \((2,2,2,1)\). Its native residue is

\[
 x\equiv420960162085323720610343395041256925057
                 \pmod{2^{129}}.
\]

Exactly as above, its Mersenne guard is

\[
 K\equiv1457\pmod{2^{126}},                        \tag{6}
\]

and it gives \(2^K-1\rightsquigarrow2^{K-6}-1\). Dropping the final source
letter 1 and partner letter 3 from the discovered cores retains the entire
variable-final-letter family, not just the sampled final valuation. Both
receipts have \(K+61\) odd letters. This phase is disjoint from the short
rule's phase and from the two preceding uncapped rules.

The declared finite discovery compared \(K=1457\), all deletions
\(1\le D\le64\), and reduced source depths at most 256, at synchronized
total odd time. It found exactly six deletion choices with a join: depths
62 for \(D=3,4,5,6\), and depth 81 for \(D=1,2\). Every found pair is
also independently checked as an exact rational collision at \(-1\).
This finite minimum is restricted to that comparison universe. It is not
a claim that every useful join must have this form.

## 5. Relative coverage and the remaining guard obstruction

The baseline is frozen, not chosen after seeing these rules: eight adaptive
macros with per-macro cap 128, then the two inverse rows and final frontier
from [uncovered join routes](collatz_uncovered_join_routes_20261007.md),
augmented by the ordinary reset and sporadic \((2,6,c)\) rules from terminal
lifts. Both named sources fail this baseline. Infinite disjoint subfamilies
are

\[
\begin{array}{ll}
 K=1459+729\cdot2^{21}s,&s\ge0,\\
 K=1457+729\cdot2^{126}s,&s\ge0.
\end{array}                                                       \tag{7}
\]

Here is the all-height argument. For either ray \(K\ge1457\); at each
checkpoint after \(d\le896\) ones, the current state has
\(v_2(x+1)=K-d\), \(v_2(x+5)=2\), and \(v_2(11x+19)=3\). The two repeat
patterns and proper sibling rules fail. Every actual prefix of length at
most 128 grows, so a core-bank descent cannot match. A virtual rule would
need more than 128 letters. The to-ROOT branch fails. Each fallback therefore
consumes precisely 128 ones, and the final checkpoint has no proper sibling.
This is the same complete trace argument used in the inherited cap-escape
proof, valid for every member of (7).

The extra factor 729 preserves \(K\pmod{1458}\), and
\(2^{1458}=1\pmod{2187}\). The two source residues are respectively 1 and
1093 modulo 2187; the latter is 40 modulo 81. Neither meets the old inverse
guards, which require 10 modulo 81 or 111 modulo 2187. The reduced prefixes
already exclude both old uncapped rules. Thus every source in (7) acquires
a new paid dependency relative to this declared baseline. We do not claim
the full dyadic phases (5) and (6) are disjoint from every earlier bank.

A cheap residual diagnostic remains useful. For odd \(K>1\), put
\(m=v_2(K-1)\). At the reduced checkpoint,
\(v_2(x-1)=m+3\). Each valuation-two step multiplies \(x-1\) by \(3/4\).
It follows that the initial reduced two-run has exactly
\(\lfloor m/2\rfloor+1\) letters. The next letter is 1 when \(m\) is
even and at least 3 when \(m\) is odd. Thus these prefixes can be arbitrarily
long. This is a depth sidecar, not a proof that no other rule works there.
The specific exponent 1451 misses all three newly frozen guards even though
its point certificate is now grounded. Guard failure must remain separate
from convergence status.

## 6. Grounding the named points, with explicit provenance

The only newly discovered ROOT word is for \(2^{1457}-1\). A declared cap
of 100,000 odd steps returns first ROOT after 7,347 steps, valuation sum
13,102, and peak bit length 2,309. The JSON stores its initial 1,456 ones
as one run and its remaining valuations as nonzero hexadecimal digits.
This is an exact word encoding, not a probabilistic checksum or a ROOT claim
inferred from a residue. The production decoder independently replays the
entire word and checks first-hit ROOT, rank, cost and peak size.

Two substitutions then use this one authenticated suffix:

| Source | Proof construction | Odd rank | Valuation cost | Peak bits |
|---|---|---:|---:|---:|
| \(2^{1459}-1\) | D2 receipt plus supplied 1457 suffix | 7347 | 13104 | 2312 |
| \(2^{1457}-1\) | frozen bounded-discovery word | 7347 | 13102 | 2309 |
| \(2^{1451}-1\) | reverse transport through the D6 receipt | 7347 | 13096 | 2301 |

For reverse transport, remove the checked source prefix from its supplied
ROOT word and prepend the checked child prefix. Both reach the identical
retained endpoint. The compiler verifies exact prefix equality and replays
the output. There is no orbit-discovery call for 1459 or 1451. The unchanged
odd rank follows because both sides of each full receipt have equal length;
the cost change is exactly the deleted run length.

These point certificates do not prove the arbitrary child on (5), (6), or
(7) reaches ROOT. Conversely, no assertion is made that the fixed rule bank
must recognize every point for which another certificate is available.

## 7. Reproduction, APIs and controls

Run from the repository root:

```text
python -B 04-computation/experiments/collatz_reset2_rules_20261007.py
python -B -O 04-computation/experiments/collatz_reset2_rules_20261007.py
```

[Script](../../04-computation/experiments/collatz_reset2_rules_20261007.py),
[frozen data](collatz_reset2_rules_20261007.json), and
[saved output](collatz_reset2_rules_20261007.out) are the complete package.
There are no default file writes. Production entry points are `audit_rule`,
`mersenne_guard`, `receipt`, `frozen_seed_word`, `reverse_discharge`, and
`grounded_named_points`. `discovery_control` is separately labelled and
called only by the finite reproducibility experiment.

Besides the declared collision search, controls cover 90 fixed-run/lift
sources, 96 symbolic exponent lifts, variable terminal valuations 1 through
12, the two baseline seeds, the exact three point certificates, and ten
all-two-prefix parameters. The deep exponent lifts are checked by modular
exponentiation; no huge \(2^K\) is expanded for them. Independent literal
replay uses the preceding implementation's repeated division; the discovery
control uses a lowest-set-bit valuation calculation.

Hostiles reject invalid integer types, forged rule lengths or valuation
fields, corrupted supplied words, insufficient initial runs, and ROOT
padding. During implementation audit, caching a dataclass before validating
its exact fields let a float or boolean alias reuse an integer cache entry.
The cache was removed; every rule now receives exact type validation on
every call. Warm integer-then-float/boolean regressions retain that boundary.

The next obligation is a sound adaptive rule choice or another independently
grounded smaller child at an arbitrary supplied residual source. These
three finite proofs and two infinite guarded reductions advance that
obligation; they do not replace it with an assertion of universal coverage.
