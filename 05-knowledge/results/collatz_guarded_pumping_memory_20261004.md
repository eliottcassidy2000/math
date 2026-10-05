# Repeating a certificate consumes prime-adic information

2026-10-04. **PROVED:** exact ternary fuel for repeating an inverse word;
the dual binary fuel formula; every regular language consisting entirely
of first-hit root certificates has bounded odd-step depth and covers only
polylogarithmically many inputs below a height bound; the full first-hit
language is nonregular, unconditionally; the two-anchor refill boundary.
**FINITE-EXACT:** the declared controls below. **OPEN:** universal positive
Collatz coverage. These statements concern a specified word encoding, not
every possible automaton, digit representation, or proof system.

## 1. Why this is a coverage question

The incoming [branch-toll and rank note](collatz_branch_toll_rank_20261004.md)
replaces numerical descent by a proper rank and proposes finite controllers
with unbounded anchor precision. The earlier
[whole-cell obstruction](collatz_partition_cover_20261004.md) rules out
uniformly bounded common-future depths. A tempting repair is to add loops
to a finite certificate grammar and let those loops generate arbitrary
depth. The result here specifies why an unguarded finite grammar cannot
do that, and which arithmetic information makes a guarded loop sound.

Closest proved mechanism: affine word transport with an exact integrality
guard. Hostile: a loop that is legal once but becomes fractional when
repeated. Corrected near miss: a finite control graph with an unbounded
integer register is not a finite-state recognizer. Least-used coordinate:
the ternary distance from the inverse word's rational fixed point.

The older [loops-and-escapes note](collatz_procgen_20260922_loops_and_escapes.md)
uses a different, extended reverse-move system. Its loops are not imported
as loops of the deterministic Collatz map used here. The
[guard synthesis](collatz_guards_20260921_synthesis.md) already emphasizes
retaining cofactors and size; no priority claim is made for the elementary
affine identities or the general automata method.

## 2. Exact ternary repetition fuel

Use the shortcut map on positive integers:

\[
 T(n)=\begin{cases}n/2&n\text{ even},\\(3n+1)/2&n\text{ odd}.
 \end{cases}
\]

Its inverse letters are

\[
 D(x)=2x,\qquad R(x)=(2x-1)/3.
\]

The latter is legal exactly when its value is a positive integer; that
value is automatically odd. Words are written in execution order. For a
word \(w\) of length \(L\), containing \(r>0\) letters R, write

\[
 F_w(x)=\frac{qx-B}{d},\qquad q=2^L,\quad d=3^r,
 \qquad \Delta_w(x)=(q-d)x-B.
\]

The carry \(B\) is a positive integer. Iteration gives the exact identity

\[
 (q-d)F_w^k(x)-B=(q/d)^k\Delta_w(x).                 \tag{1}
\]

Because \(3\nmid q(q-d)\), for an integer starting value \(x\),

\[
 F_w^k(x)\in\mathbb Z
 \iff 3^{rk}\mid\Delta_w(x).                      \tag{2}
\]

An integral endpoint also makes every intermediate inverse letter
integral. Each intermediate value lies in \(\mathbb Z[1/3]\) by inverse
construction, and in \(\mathbb Z[1/2]\) by reconstructing it from the
integral endpoint using the forward affine letters. The intersection is
\(\mathbb Z\). If the endpoint is positive, the same forward reconstruction
makes every intermediate value positive. The prescribed parity is then
automatic at each inverse edge.

Consequently, if \(\Delta_w(x)\ne0\), the **exact integrality budget** is

\[
 k\le\left\lfloor\frac{v_3(\Delta_w(x))}{r}\right\rfloor. \tag{3}
\]

Positivity and a first-hit stopping condition are separate tests. In
particular, integrality alone does not certify an indefinitely repeating
positive path. If \(\Delta_w(x)=0\), then \(x=B/(q-d)\) is a fixed point
of the composite word. This is the exceptional cycle case, not infinite
fuel for a new route to a root.

This is also a compression rule. A verifier may store the word's affine
data, its repeat count, the divisibility (2), the positive endpoint, and
the first-hit condition, instead of expanding every repeated letter.
The general first-hit test must still be supplied; it cannot be erased by
the affine compression.

## 3. Finite-state certificate banks have bounded odd depth

Let \(\mathcal P\) be the language of words in D,R which, when executed
from 1, have positive integer values and never revisit 1 after time zero.
Include the empty word. Reading such a word backwards gives an actual
first-hit T-route from its endpoint to 1. This language is defined without
assuming that every positive integer has such a route.

**Theorem.** If a regular language \(\mathcal L\) is a subset of
\(\mathcal P\), then there is a finite \(K\) such that every word in
\(\mathcal L\) has at most \(K\) letters R. An N-state finite automaton
recognizing \(\mathcal L\) permits \(K\le N-1\).

**Proof.** Discard states that lie on no accepting path. Suppose a cycle
in the remaining automaton contains R. Choose an accepted path
\(u w z\) that traverses that cycle, so \(u w^k z\) is accepted for all
\(k\ge0\). The integer \(x\) reached after u is fixed. All iterates
\(F_w^k(x)\) must be integral. Equation (3) forces
\(\Delta_w(x)=0\), so the block w returns to x.

But an actual inverse path from 1 with no subsequent visit to 1 cannot
repeat any vertex: forward iteration from a repeated vertex would follow
a cycle while its already supplied path also leads to 1. Determinism
forces that cycle to contain 1, contrary to first-hit status. Thus there
is no R-containing cycle on an accepting path. Every R edge crosses
between different strongly connected components; their condensation is
acyclic and has at most N vertices. There are at most N-1 such crossings.
The argument also allows nondeterminism and epsilon edges. QED.

**Quantitative coverage consequence.** If a valid first-hit inverse word
has d letters D and r letters R and endpoint \(n\le X\), then

\[
 d\le\lfloor\log_2 X\rfloor+r.                   \tag{4}
\]

Indeed D doubles, while every legal R on a nonroot positive input x>=2
has \(R(x)\ge x/2\). Thus \(n\ge2^{d-r}\).
For r<=K, put \(H=\lfloor\log_2 X\rfloor+K\). Counting the r+1
gaps of D letters gives the explicit upper bound

\[
 \#\{n\le X:\text{some word in }\mathcal L\text{ certifies }n\}
 \le\sum_{r=0}^K {H+r+1\choose r+1}
 =O_{\mathcal L}((\log X)^{K+1}).                 \tag{5}
\]

Every sound regular first-hit certificate bank therefore has density-zero
input coverage. A finite union of such banks has the same limitation.
This is stronger than saying that a particular bounded search failed.
It does **not** apply to a finite program with arbitrary integer registers,
unbounded valuation guards, a stack, or a different compressed language.

## 4. Nonregularity does not depend on the Collatz conjecture

There are actual first-hit words with arbitrarily many R letters. For
each \(r\ge2\), take

\[
 w_r=D^{3^{r-1}}R^r.
\]

The elementary lifting identity
\(v_3(2^{3^{r-1}}+1)=r\) gives exactly r legal repetitions of R after
the initial doubling block. Its endpoint is

\[
 n_r=\frac{2^r(2^{3^{r-1}}+1)}{3^r}-1>1.
\]

The R portion decreases to this value, so every intermediate value is
greater than 1. The initial D portion also has this property. Thus
\(w_r\in\mathcal P\) for all r>=2. The theorem proves
\(\mathcal P\) is nonregular without any universal-convergence assumption.

Positive control: \(D^{2k+1}RD^b\), k>=1,b>=0, is a sound infinite regular
language. It has only one R. Hostile boundary: DR is the reverse root
cycle \(1\to2\to1\), so it can be pumped integrally but is excluded
from \(\mathcal P\). Omitting first-hit status would invalidate the stated
bounded-R theorem.

## 5. The forward mirror and a sharply located refill event

For a forward parity word of length L with r odd letters, write

\[
 A_w(n)=\frac{3^r n+C}{2^L},\qquad
 \Gamma_w(n)=(2^L-3^r)n-C.
\]

The same calculation, exchanging the denominator prime, gives

\[
 (2^L-3^r)A_w^k(n)-C
  =(3^r/2^L)^k\Gamma_w(n),\qquad
 A_w^k(n)\in\mathbb Z\iff2^{kL}\mid\Gamma_w(n).   \tag{6}
\]

As above, integer endpoints restore integer prefixes and their required
parities. This is the binary mirror of (2), and specializes to the
incoming rational-anchor fuel criterion. The two directions consume
different prime-adic resources even though their formal affine maps
belong to a solvable group.

There is a useful restriction on changing anchors. For distinct rational
anchors \(\alpha,\beta\) with odd denominators, let
\(\kappa=v_2(\alpha-\beta)\). If
\(K=v_2(n-\alpha)\ne\kappa\), then

\[
 v_2(n-\beta)=\min(K,\kappa).                     \tag{7}
\]

This is the exact ultrametric identity. Arbitrarily large new precision
at beta can appear only at the **equality boundary** K=kappa (or at beta
itself). Thus a finite anchor collection does not permit unexplained
precision replenishment: large refills must pass one of finitely many
anchor-distance levels, and the retained cofactor determines the new
depth there.

The boundary is real. With alpha=-1, beta=1 and n=1+2^j, j>=2, old
precision is exactly 1 while new precision is j. Therefore a bound on
pairwise anchor distances does not bound all future refills. The incoming
proper rank can pay for some such switches: for even j>=6,
\((2^{j+1}+1)/3\to1+2^j\) is an actual odd step with lower root-centered
rank. This is a specialization of that note's proved forward comparison,
not a new global decrease claim. At j=4 the source11 is its known exception.

## 6. What to retain in the next coverage construction

The connection contract is:

| Source -> target | Preserved | Lost without an extra coordinate | Decisive test |
|---|---|---|---|
| Inverse block -> affine triple (q,d,B) | Exact endpoint and composition | Positive intermediate path and first hit | Divisibility, endpoint positivity, root check |
| Repeated block -> ternary fuel | Exact number of integral repeats | Which extension pays original source rank | Keep original rank and final child |
| Forward block -> binary fuel | Exact repeated parity word | Behavior after the guard fails | Explicit exits and anchor-switch cofactor |
| Finite proof grammar -> regular language | Word acceptance | Unbounded ternary compatibility | Pump an R-containing automaton cycle |
| Finite anchor atlas -> collision levels | Ultrametric transition law | Unbounded new depth at equality | Keep the equality cofactor; test1+2^j |

Anchor: universal coverage by complementary guarded rules. Niche: exact
compression and limits of finite certificate languages. Wildcard: binary
and ternary fuel as the information lost by an unguarded group quotient.
The live board is original rank / repeat word / binary fuel / ternary fuel /
anchor-switch cofactor / grounded root certificate.

The productive repair is a finite **description** with unbounded checked
arithmetic memory. Every repeat and switch should state how much precision
it consumes, where new precision comes from, and why its final dependency
is below the immutable original rank. The current results make these
obligations explicit. They do not prove that every supplied input has a
paying extension.

## Reproduction and exact scope

Run the [script](../../04-computation/experiments/collatz_guarded_pumping_memory_20261004.py)
with `python -X utf8 -B`, and again with `-O`; compare the saved
[stdout](collatz_guarded_pumping_memory_20261004.out).
Both runs agree after LF normalization; SHA256:
`4bd4b410bbf5cb4692996488d4dd7f28178374e4da1a12683c8afb95c9fa0783`.

Universe: all inverse and forward binary words of lengths1..7, starts1..64,
repeat counts0..5; every inverse binary word through length16; every
R-containing subword of the accepted first-hit paths in that universe;
the explicit unbounded-depth witnesses r2..9; 1023 one-R regular controls;
all ordered distinct pairs of five rational anchors on1001 positive odd
inputs; and the explicit counting bound at X=10^6 for K0..3. Independent
literal rational replay checks affine formulas, and forward replay checks
every finite root certificate. Root-cycle and failed-pump controls remain
visible. The proofs carry the unbounded quantifiers.
