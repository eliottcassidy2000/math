# Both reset two boundary types admit smaller dependency rewrites

2026-10-04. **PROVED** the eight guarded all-height common-future families,
their disjoint densities, and the necessary least-counterexample restriction.
**FINITE-EXACT** the least costs within the stated word search.
Universal Collatz coverage remains **OPEN**. A smaller dependency needs its
own route certificate; none of the density statements counts counterexamples.

Artifacts: [exact search](../../04-computation/experiments/debt_word_join_search_20261004.py)
and [saved output](debt_word_join_search_20261004.out).

## 1. Inheritance and why the search changed

The [reset-two debt note](reset_two_debt_families_20261004.md) proves a
common-future rewrite at odd debt level one. Its short boundary-to-sibling
ansatz fails at level two. That is a failure of a specified word shape.
It does not exclude longer common futures. The present search retains the
full affine carry and permits a longer word, producing a level-two identity
at total source-word cost twelve. Applying the same decoder to the other
terminal boundary type gives four additional disjoint families.

The closest proved mechanism is the exact carry decoder, already used in
[ordered word observations](geometry_collatz_drift_carriers_20261004.md).
The canonical hostile is reset-two source7, which is outside these selected
families. The corrected near miss is treating the failed short ansatz as an
obstruction to every common future. The least-used coordinate is the target
word's length: it must exceed the boundary word's length by the debt level.

Anchor: smaller proof dependencies inside the reset-two kernel. Niche:
decode a target word from its required affine carry. Wildcard: let an ansatz
failure specify which word dimension to enlarge. The live board is
**original source / child / debt level / ordered carry / dyadic guard /
certificate rank**. No tournament is imposed on the partial proof relation.

## 2. A finite algebraic search with an exact decoder

Write U(n)=oddpart(3n+1). For a positive valuation word u of length l, cost S
and carry B, its formal map is

    F_u(n)=(3^l*n+B)/2^S.

At debt level e>=1, the boundary source is Y=2^s*3^e*M+1 with M positive odd
and s in{1,2}, the two terminal dyadic exponents of the retained debt state.
A second positive word v gives the identical affine endpoint F_u(Y)=F_v(M)
if and only if

    length(v)=l+e, cost(v)=S-s, carry(v)=(3^l+B)/2^s.     (1)

Necessity follows by comparing slopes and using unique factorization, then
comparing constants. Sufficiency is substitution. At s=1 every actual Y has
first valuation one. At s=2 it is5 modulo8, so its first valuation is at
least three. Words u outside these respective guards are impossible here.

There is at most one positive word of specified length, cost and carry.
For length d>1 and carry C, its first exponent is

    a=v2(C-3^(d-1)),
    remaining carry=(C-3^(d-1))/2^a.                    (2)

Indeed C=3^(d-1)+2^a*C_tail and C_tail is odd. Reject nonpositive differences
or a<1, subtract a from the remaining cost, and repeat. At length one the
carry must be one, and the remaining cost is the last positive exponent.
These steps prove both uniqueness and the rejection tests. No trajectory
simulation or approximate logarithm is used to select the identities.

The computation enumerates **all 65535 positive compositions of total cost
one through sixteen**, with no length cap. For each pair s in{1,2} and
e in{1,2,3,4} it applies
(1)--(2). The selected first successful costs and their generalized last
exponents are:

| s | e | u acting on 2^s*3^e*M+1 | v acting on M | Least source cost at a=1 |
|---:|---:|---|---|---:|
|1|1|(1,6,a)|(1,2,1,a+2)|8|
|1|2|(1,6,4,a)|(2,1,1,1,3,a+2)|12|
|1|3|(1,2,9,a)|(1,1,1,1,1,2,a+4)|13|
|1|4|(1,14,a)|(2,2,1,3,3,1,a+2)|16|
|2|1|(6,a)|(1,1,a+2)|7|
|2|2|(10,a)|(3,2,1,a+2)|11|
|2|3|(10,a)|(1,1,1,3,a+2)|11|
|2|4|(3,1,11,a)|(1,1,2,1,1,1,2,a+4)|16|

Here a is any positive integer. The least costs are **FINITE-EXACT** within
the declared word model; they do not minimize every possible certificate
representation. Increasing both final exponents by the same amount leaves
the carries unchanged and preserves (1), proving every a in the table.

## 3. Transport back to a smaller original-source dependency

For any r>=1, use the following two words, where 2^e means e copies of the
letter two in this display:

    source n:       (1^r, 2^e, u),
    child (n-1)/2:  (1^(r-1), 2e+s, v).                (3)

Require n to lie in the exact source-word cylinder. For a word w of cost A
and carry B_w, this is the one odd residue

    n=(2^A-B_w)*(3^length(w))^(-1) mod2^(A+1).          (4)

The exact-cylinder lemma in [the return-tail note, section2](modular_return_debt_20261004.md)
proves that (4) is equivalent to the displayed positive valuations. One
direct proof recovers the first valuation by reduction modulo2^(a_1+1),
then divides it out and repeats; the final numerator has valuation exactly A.

To see (3), write n=2^(r+1)*t-1 with t odd. The child's first r-1 steps give
3^(r-1)*2*t-1. A child reset of size 2e+s leaves

    M=(3^r*t-1)/2^(2e+s-1).

The source's first reset-two endpoint is 3*2^(2e+s-2)*M+1. Each further
valuation-two step sends 3^j*2^h*M+1 to 3^(j+1)*2^(h-2)*M+1 while h>=3.
After e such source resets, its boundary is exactly 2^s*3^e*M+1. Apply (1).
Alternatively the composed source and child maps have equal length,
Q_source=2*Q_child, and B_source=2*B_child-3^length; substitution proves
their equality directly. Their exact final odd endpoint transfers (4)
to the child word, proving all required intermediate guards.

The child is positive and strictly below n. Neither displayed word can visit
one before its final edge; the final endpoint need not be one. The child's
last exponent is at least three, which excludes an
earlier visit to the root with its valuation-two self-edge. The source's
preterminal values at s=1 are respectively

    (27M+7)/64, (243M+85)/1024,
    (729M+25)/2048, (729M+7)/16384.

Setting these equal to one requires M=19/9,313/81,2023/729,5459/243,
respectively, none an integer. The four corresponding s=2 equations require
M=5/3,85/9,85/27,8173/2187, again nonintegral. A still earlier root visit
would contradict the subsequent large exponent in u. The sibling relation
has the child preterminal value equal to4 times the source value plus1, or
16 times the source value plus5 in the two rows with last offset four;
its direction alone does not exclude a source preterminal root.

Both prefixes have the same number of odd steps. The source prefix has one
additional halving, so a **supplied child home certificate** transfers with
the same odd rank and ordinary rank increased by one. Matching (4) alone
provides a smaller dependency, not the missing certificate for that child.
Other available actual endpoints and sibling certificates should remain
eligible; the [adaptive-selector hostile](adaptive_boundary_selector_20261004.md)
shows why choosing only one numerically small alternative can lose a proof.

## 4. Arbitrarily long growing prefixes and exact kernel densities

Fix a=1. Every nonempty source-side prefix has coefficient greater than one
for all r at or above the following sharp thresholds:

| s | e | First such r | A guarded source at that r | Smaller child | Common endpoint |
|---:|---:|---:|---:|---:|---:|
|1|1|8|236031|118015|478505|
|1|2|13|73580543|36790271|159293099|
|1|3|15|13581352959|6790676479|24807944825|
|1|4|24|307417686474751|153708843237375|674602512456395|
|2|1|9|445439|222719|903035|
|2|2|16|2629959679|1314979839|4270161683|
|2|3|17|14322237439|7161118719|26161257809|
|2|4|21|10898995412991|5449497706495|21259500274817|

At each displayed r, compare the finitely many integers3^j and2^A for its
prefixes. Increasing r multiplies every fixed-tail coefficient by3/2,
while the new initial coefficients are themselves powers of3/2. At r-1 a
prefix coefficient is below one. This proves the thresholds for all r,
not merely for the displayed examples. Positive carries then show that
every observed actual iterate exceeds its original source. The smaller
dependency can therefore appear before any actual or coefficient descent.

Union over all a>=1 in a row. For each fixed r it is exactly the cylinder
for the prefix ending just before a, of modulus2^(r+c_(s,e)), where

    s=1: (c_1,c_2,c_3,c_4)=(10,16,19,24),
    s=2: (c_1,c_2,c_3,c_4)=(9,15,17,24).

Different r have distinct initial runs of ones; different e have distinct
following runs of twos. After those twos, s=1 has next valuation one while
s=2 has next valuation at least three. Thus the eight rows and all r are
disjoint. Their absolute natural densities are2^-c_(s,e). To justify the infinite union, the
tail r>R is contained in n=-1 mod2^(R+2), whose upper density tends to zero.
This avoids assuming countable additivity for natural density.

The inherited least-counterexample restriction is precise. Immediate
valuation at least two descends. An initial run of ones followed by a reset
at least three joins the smaller child; n=3 is handled directly by3->5->1.
Every positive n has a finite initial run. Hence a least positive non-home
integer, if one exists, must have first reset two. That necessary set has
absolute density1/8, or1/4 among odd integers, as proved in the reset-two note.
Every family (3) lies inside it and supplies a smaller dependency, so a least
counterexample lies outside their disjoint union. The remaining necessary
set has exact density

    removed absolute density = 25041/8388608,
    remaining absolute density = 1/8-25041/8388608=1023535/8388608.

The latter is the density among all integers, or
**1023535/4194304 among odd integers**.
This is a necessary-domain calculation. It is not a density of unresolved
inputs under all known rules, nor a density of nonconvergent inputs.
No claim of disjointness from all prior descent or sibling rules is made.

## 5. Supplied-source selector and certificate transport

The executable `select_join(n)` recovers the boundary coordinates directly:

    r=v2(n+1)-1, t=(n+1)/2^(r+1),
    b=1+v2(3^r*t-1), e=floor((b-1)/2), s=b-2e.

It handles only r>=1 and b>=3, since these are the first-reset-two sources.
The inherited normalization proves the recovered e,s give the actual boundary.
It looks up the selected row, checks its exact original-source prefix cylinder,
and obtains the final exponent from the actual odd frontier. Failure returns
no rule, with no claim that the supplied source fails to converge. Success
returns the smaller child, endpoint, and both guarded words; no root search
is performed.

`transport_certificate(n, child_certificate)` separately validates the supplied
AST and its source identity, checks that it contains the required child word,
retains the rooted common-endpoint suffix, and prepends the original-source
word. It checks the recovered source and both full route ranks. The caller
must supply the child certificate; matching a residue is not permission to
invent it or change its source.

For all10000 odd inputs through20000, an independent literal prefix parser
finds exactly the same60 applicable sources, with no duplicate selected rows.
This is an applicability census, not a home census. Eight additional controls
use the growing examples in section4: a bounded search of at most5000 odd
steps explicitly supplies each child premise; transport then agrees with a
separately encoded full source route. Four malformed inputs and an unrelated
child certificate are rejected. The production transport itself performs no
such search. No claim is made that all60 inputs were newly covered compared
with the repository's other selectors.

## 6. Reproduction and the next boundary

    python -X utf8 -B 04-computation/experiments/debt_word_join_search_20261004.py
    python -O -X utf8 -B 04-computation/experiments/debt_word_join_search_20261004.py

The full source-word search also round-trips all65535 carries. An independent
separator-subset generator and direct carry sums check all words through
cost ten. There are768 literal first-hit boundary/source controls, spanning
all eight rows, a in{1,2,3,7}, r in{1,2,8,13,21}, and four source heights.
The eight growing examples, eight preterminal-root obstructions, and
preceding-threshold hostiles are additional.
Checks remain active under-O.

The source-to-target map is (1), followed by the guarded lift (3). It preserves
the common future and smaller-source certificate transport. Forgetting the
exact cylinder loses integer legality; forgetting the rooted child suffix
loses completion; forgetting the ordered carry loses the target word.

**OPEN:** obtain a total rule for the remaining debt states, or a proved
well-founded way to move them into accepted states. The existence of some
accepted cylinders at both boundary types and four debt levels neither handles
every source there nor proves that every level admits such a word. The decisive next
experiment is adaptive word-length extension at a *supplied unresolved
state*, keeping its actual guards and existing root dependencies, rather
than using the bounded counts as evidence for global coverage.
