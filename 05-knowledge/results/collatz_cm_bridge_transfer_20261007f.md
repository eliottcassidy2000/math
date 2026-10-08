# A protected-source auxiliary bridge and its exact payment budget

**PROVED:** the elementary sibling-phase compiler, finite exact payment search,
and the infinite family below. **FINITE-EXACT:** the declared controls and two
supplied ROOT discharges. **CITED / USER-TRUSTED INPUT:** the linked modularity
paper's headline and architecture; its theorem is not a premise of any Collatz proof.
**OPEN:** universal paid coverage and ROOT grounding of arbitrary emitted children.

Artifacts: [program](../../04-computation/experiments/collatz_cm_bridge_transfer_20261007f.py)
and [output](collatz_cm_bridge_transfer_20261007f.out). No mathematical priority
claim is made for sibling identities or the inverse-word algebra.

## 1. Source architecture and inherited interfaces

The user supplied the release
[Every Elliptic Curve over a CM Number Field Is Modular, v1.0.0](https://github.com/CaptainSude/elliptic-curves-cm-modularity/releases/tag/v1.0.0).
The complete [release LaTeX](https://raw.githubusercontent.com/CaptainSude/elliptic-curves-cm-modularity/v1.0.0/elliptic-curves-cm-modularity.tex)
was read, without a new audit of its headline. Its sections 2–6 prepare one
auxiliary specialization satisfying simultaneous local conditions, use
integral evaluation lattices for three congruences through coefficient primes
`p,3,p`, and retain a protected Steinberg place through characteristic-5
lifting, patching, and field preparation. The useful design lesson is to
carry the protected condition through every interface, rather than recover
it from a reduced numerical invariant. This is an architectural transfer,
not an identification of Collatz maps with Galois representations.

The closest proved arithmetic mechanisms are the actual sibling relation
`S(n)=4n+1`, the native progression composition in
[carry interfaces, §3](collatz_carry_interfaces_20261004.md), and the
authenticated common-future discharge in
[marked completion, §2](collatz_marked_completion_20261007c.md).
[Grounded child closure](collatz_grounded_child_closure_20261007.md) instead
varies an already grounded source parameter. Here the supplied source stays
fixed: only its auxiliary sibling depth is selected. The canonical hostile
is 7, proved below. The corrected near miss is to mistake solvability of a
congruence for an affordable child or a completed ROOT certificate. The
least-used sidecar is the original-source inequality after the auxiliary lift.

The concept board is **immutable source / auxiliary sibling / native inverse
word / exact phase / original budget / supplied grounding port**.

| Source object | Arithmetic map | Preserved predicate | Lost data if omitted | Required sidecar |
|---|---|---|---|---|
| simultaneous auxiliary specialization | select sibling depth `k` for the fixed `n` | inverse-word integrality | payment to the original `n` | strict inequality (3) |
| integral bridge between marked objects | literal words meeting at `U(n)` | actual common future | native valuations and early ROOT | ordered words and exact source |
| protected local condition | retain `n` through all interfaces | source ownership | a nearby legal source is not this source | recomputed receipt fields |
| final transfer from grounded auxiliary object | supplied child ROOT word | first-hit ROOT | existence of that child proof | authenticated child suffix |

## 2. One auxiliary phase, but a finite payment window

Write the odd map as `U(n)=(3n+1)/2^v2(3n+1)`. Let `n>1` be a supplied
positive odd integer and `w` a nonempty positive valuation word, with

\[
 F_w(x)=(Px+B)/Q,\qquad P=3^L,\quad Q=2^A.
\]

For every integer `k>=0`, set

\[
 S^k(n)=\frac{4^k(3n+1)-1}{3},\qquad
 h_k=\frac{Q4^k(3n+1)-Q-3B}{3P}.                 \tag{1}
\]

**Theorem 1 (exact auxiliary guard).** There is exactly one class
`k=k0 mod3^L` for which `h_k` is integral. At every nonnegative depth in
that class, `h_k` is positive odd and `w` is the actual word from `h_k`
to `S^k(n)`. Neither that word's start nor its intermediate states equal 1.

Indeed, the guard is

\[
 4^k\equiv (Q+3B)\,[Q(3n+1)]^{-1}\pmod{3^{L+1}}. \tag{2}
\]

Its right side is 1 modulo 3. The principal units modulo `3^(L+1)` have
`3^L` elements, and 4 has exactly that order: binomial lifting gives
`v3(4^u-1)=1+v3(u)` for nonzero `u`. Hence the unique phase can be lifted
one ternary digit at a time. Full inverse-word integrality forces each
intermediate inverse to be an odd integer, by successively reducing the
carry identity modulo 3. A nonpositive odd integer stays negative under
every exact positive-valuation step, so it cannot reach the positive
endpoint. A positive path passing through 1 stays there under its only
actual continuation, of valuation 2, and cannot reach `S^k(n)>1`.

**Theorem 2 (complete finite budget test).** Among the legal depths, strict
payment `h_k<n` is equivalent to

\[
 Q(3n+1)4^k < 3Pn+Q+3B.                         \tag{3}
\]

Define

\[
 C=\left\lfloor\frac{3Pn+Q+3B-1}{Q(3n+1)}\right\rfloor.
\]

The admissible positive depths are exactly the members of the class (2)
between 1 and `floor(log4 C)`; if `C<1` the set is empty. The implementation
uses integer bit lengths, not floating logarithms. This is a total finite
decision for the supplied pair `(n,w)`. It neither chooses a replacement
source nor promises success as words vary.

If `a=v2(3n+1)`, the legal bridge has the explicit common future

\[
 U_{(a)}(n)=U_{w\,(a+2k)}(h_k).                 \tag{4}
\]

There is no early ROOT in either prefix: the child first reaches
`S^k(n)>1`, then follows one actual edge. The common endpoint may be ROOT.
Given a supplied first-hit child ROOT word, the compiler verifies its
source, every valuation, and its prefix `w,(a+2k)`, then replaces that
prefix by `(a)`. Its production API performs no ROOT discovery.

## 3. A sharp infinite family inaccessible to direct predecessors

For `L>=1`, take `w_L=1^(L-1),2`, so

\[
 P=3^L,\quad Q=2^{L+1},\quad B=3^L-2^L.
\]

At the single auxiliary lift `k=1`, (1) becomes

\[
 h_L+1=\frac{2^L(8n+3)}{3^L}.                  \tag{5}
\]

**Theorem 3.** For a positive odd source `n`, this word is native exactly
when `3^L | (8n+3)`. On its native domain it is paid exactly when `L>=6`.
Consequently this family emits a smaller child exactly when
`v3(8n+3)>=6`. Among its legal paid depths, the unique smallest child is
obtained at `L=v3(8n+3)`.

For the proof put `t=(8n+3)/3^L`, a positive odd integer. Then

\[
 n-h_L=\frac{(3^L-2^{L+3})t+5}{8}.             \tag{6}
\]

For `1<=L<=5`, the coefficient of `t` is at most `-13`, so this difference
is negative. For `L>=6`, it is positive, starting with `729-512=217`.
Increasing a legal depth multiplies `h_L+1` by `2/3`, proving the minimum
claim. The original source is unchanged during this optimization.

The first paid domain and child are

\[
 n=273+1458t,\quad h_6=191+1024t,\quad t\ge0.   \tag{7}
\]

Every such source is divisible by 3. Yet every ordinary odd Collatz
endpoint is a ternary unit. Thus these sources have **no nonempty direct
odd predecessor word at any depth**, while the auxiliary sibling bridge
does give a smaller common-future child. This is a distinction between
two operations, not a claim of novelty or disjointness from every previous
paid selector. The domain has relative natural density `1/729` among odds.
It is a paid dependency domain, not an all-height ROOT claim.

On the subdomain with odd `t`,

\[
 n=1731+2916u,\quad h=1215+2048u,
 \qquad U_{(1)}(n)=U_{(1,1,1,1,1,2,3)}(h).     \tag{8}
\]

The first source step grows, so this control does not get its payment from
an immediate forward descent. Two independently supplied finite ROOT words
ground `273` via child `191`, and `1731` via child `1215`. They yield 8 and
53 odd steps respectively. The maximum-depth option at `273` has child
`127`, instead of `191`; no minimum over arbitrary common-future diagrams
is asserted. Smaller numerical size does not certify an available ROOT
suffix. The explicit `depth` argument retains the alternative `191` receipt
when that is the child for which a proof has actually been supplied.

## 4. Strong failure boundary and the current residual source

**All-depth hostile.** Source 7 has no smaller child obtainable by taking
any sibling `S^k(7)`, then any inverse actual valuation word. Its possible
smaller positive odds are `1,3,5`; their full future is contained in the
closed set `{1,3,5}`. It cannot reach `S^k(7)>=7`. Source 7 nevertheless
has the actual ROOT word `(1,1,2,3,4)`. Thus this auxiliary mechanism is
not a necessary condition for ROOT and cannot be universally successful.

A separate exploratory beam on the current residual Mersenne source
`n=2^99708993705-1` used sibling depths 1 through 8, inverse letters 1
through 6, depth at most 100, and width 2,000 sorted by accumulated binary
cost. It found no contracting bridge; depths 2,5,8 were already divisible
by 3. This is a bounded heuristic stopping observation, not exhaustive
word failure, and is not part of the reproducible proof controls below.
The parent session subsequently found a deeper, different collision route.

## 5. Verification and boundaries

The standalone script verifies all 84 words of lengths 1–3 over letters
1–4, each on the 32 odd sources 3–65. It independently scans an entire
ternary phase period and compares the complete positive payment list.
Family controls cover 512 members of (7), plus depths 1–24 with 12 legal
native parameters each. Literal replay is independent of the carrier phase
calculation. Malformed bool/float fields, altered source, wrong child proof,
missing guard, unpaid lift, and ROOT-padding inputs are rejected.

Reproduce in the repository root:

```text
python -B -X utf8 04-computation/experiments/collatz_cm_bridge_transfer_20261007f.py
python -B -O -X utf8 04-computation/experiments/collatz_cm_bridge_transfer_20261007f.py
```

Both modes must equal the saved output. The paper provides the protected
bridge design; equations (1)–(8), the source guard, affordability decision,
and ROOT splice are proved here with elementary arithmetic. They do not
derive any individual Collatz conclusion from modularity or a local
congruence alone.
