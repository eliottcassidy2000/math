# A guarded sibling proof grammar: retain intermediate ladder heights, not just a canonical representative

**Status.** PROVED: the signed sibling normal form and clock identities;
the guarded common-future certificate transport; an infinite family beyond
every-step numerical descent; the exact local residue coverage; and the
forward-closure obstruction to gaining new certificates from inverse
ports alone. FINITE-EXACT: the generated-family census and compiled ordinary
root-4 certificates below. Collatz convergence and universal termination
of the proof search remain OPEN. This is a new synthesis of inherited
inverse-fibre identities, with no novelty claim for those identities.

## Inheritance, portfolio, and concept board

Closest proved mechanism: equations B1--B3 in
[the inverse-fibre braid](arithmetic_braids_20260917_collatz.md), especially
`B(4n+1)=B(n)` and the complete valuation-parametrized inverse fibre.
The earlier [bug note](collatz_bugs_20260925.md) connects the same ladder
to the proposed root geometry. Canonical hostile: the inverse-fibre note's
arbitrarily long exponent-one growth blocks. The
[extended-graph note](collatz_mod6_20260917_extended_collatz_scc.md) is a
second hostile: allowing even sources of `3n+1` changes the problem, and
its inverse-density theorem is not a universal deterministic certificate.

The corrected near miss in this experiment is the tempting choice to
keep only the smallest sibling representative. It loses a usable
certificate already at `241`, and keeping just the raw and smallest
representatives still misses `483`. The least-used sidecar is the
intermediate ladder height carrying an already proved certificate.

Anchor: a legal decoder to ordinary root `4`. Niche: a finite proof grammar
with strictly smaller certificate dependencies. Wildcard: test the same
identities on the minus sheet rather than mistaking local algebra for
one-basin convergence. The live board is odd valuation clocks, binary
suffixes, sibling height, certificate dependencies, and actual orbit paths.
The meaningful change of direction was from a single normal form to all
guarded certificate-bearing ladder positions.

Prior-work search recovered the inverse-fibre and extended-graph results
before deriving the rules below. Searches for the explicit `(8n+1)/3`
join schema and its base family did not identify a separate existing
theorem; the underlying algebra is nevertheless an immediate consequence
of the inherited fibre theorem and is not advertised as novel.

## 1. Signed objects, guards, and normal forms

For `s in {+1,-1}` and positive odd `n`, define

    a_s(n)=v_2(3n+s),
    B_s(n)=(3n+s)/2^a_s(n),
    S_s(n)=4n+s.

The map `B_s` uses the odd-to-odd clock. The ordinary map `C_s` sends an
odd integer to `3n+s` and an even integer to `n/2`; one `B_s` step costs
`a_s(n)+1` ordinary steps. In the single-halving shortcut convention it
costs `a_s(n)` steps. These clocks are kept separate throughout.

The exact signed sibling identity is

    3 S_s(n)+s=4(3n+s),
    B_s(S_s(n))=B_s(n),       a_s(S_s(n))=a_s(n)+2.          (1)

Define `rho_s(n)` by repeatedly applying `(n-s)/4`, with the guard

    plus sheet:  n=5 mod8,
    minus sheet: n=3 mod8.

Each guarded replacement is a positive odd integer strictly smaller than
its source. The process removes exactly `floor((a_s(n)-1)/2)` sibling
levels and stops with valuation one or two. It is the unique smallest
positive odd member of the same `B_s` fibre, by the inherited complete
inverse-fibre formula. Write `j_s(n)` for the number of levels removed.

On the plus sheet, `S_+(n)` appends binary digits `01`. The inverse suffix
deletion is legal only when the remaining prefix is odd, exactly the
guard `n=5 mod8`. On the minus sheet, `(n+1)/4` has a carry; it is not
the same binary deletion with the sign ignored.

Set `Q_s=rho_s o B_s`. Because `B_s o rho_s=B_s`,

    Q_s^t(n)=rho_s(B_s^t(n))             for t>=1,
    B_s^(t+1)(n)=B_s(Q_s^t(n))          for t>=0.            (2)

Thus `Q_s` preserves the root-reachability question and the cycle
obstruction. It does not prove convergence merely by changing the object.
The useful step is the guarded certificate rule obtained from the full
fibre, not a claim that the normal form solves the dynamics.

## 2. The common-future rule and a noncircular compiler

Suppose `n,a` are positive odd integers, `j>=0`, and

    B_s(n)=S_s^j(a).                                      (J)

Then

    B_s^2(n)=B_s(a).                                      (3)

Assume an actual finite certificate from `a` to odd root `1` is supplied.
For `a!=1`, its first edge is necessarily `a -> B_s(a)`. Replace this
first vertex by `n, B_s(n)` and retain its remaining suffix. Equation (3)
proves every new edge is legal. If `a=1`, use the known fixed odd step
`B_s(1)=1`, obtaining `n -> B_s(n) -> 1` directly. This is an explicit
certificate transformation, not an assumption that the new source
converges.

For induction, require the additional rank guard

    a<n.                                                  (R)

Certificate dependencies then strictly decrease in the positive-integer
order. The actual orbit may increase: `n -> B_s(n)` is permitted to go
above `n`, because the proof dependency is `a`, not that orbit value.
No edge `n -> a` is claimed. The certificate is a common-future diagram.

The direct forward rule is included by taking `j=0` and `a=B_s(n)<n`.
The canonical rule chooses `a=Q_s(n)`. The full rule keeps every legal
sibling ancestor of `B_s(n)` satisfying (R).

**Clock accounting.** Under (J), the path from `n` to `B_s(a)` takes two
odd steps, versus one from `a`. Its extra shortcut cost is
`a_s(n)+2j`; its extra ordinary cost is `a_s(n)+2j+1`, because
`a_s(B_s(n))=a_s(a)+2j`. These are costs to the displayed common future;
an earlier visit to a chosen root may shorten the final certificate.

**Ordinary root 4 on the plus sheet.** Expand every checked odd edge into
the ordinary `3n+1` step followed by its exact number of halvings. After
odd root `1`, append the legal edge `1 -> 4`, then cut the path at its
first visit to `4`. The resulting finite path consists solely of ordinary
Collatz edges. Its reversal is the requested guarded route from `4`.
There is no appeal to the conjecture, including in the base case.

## 3. Keeping all intermediate heights gives a real gain

There are three different choices of dependency, and they are not
interchangeable under the strict numerical rank.

**Canonical-only loses a forward certificate.** The direct path

    241 -> 181 -> 17 -> 13 -> 5 -> 1

strictly decreases at every odd step. But `Q_+(241)=11`, since the
sibling chain of `181` is `181 -> 45 -> 11`; the canonical-only grammar
has no strict-descent certificate for `11`. Canonicalization preserves
the eventual orbit question while discarding this particular rank proof.

**Raw-plus-canonical still loses an intermediate certificate.** At `483`,

    B_+(483)=725,
    sibling ancestors of 725: 725,181,45,11.

The raw value `725` is larger than `483`, and the canonical value `11`
is not generated. The intermediate value `181<483` has the certificate
above. Rule (J) therefore produces

    483 -> 725 -> 17 -> 13 -> 5 -> 1.

This source is certified only after the missing ladder position is retained.
The example is a concrete reason to store several certificate alternatives
instead of identifying them with the smallest representative.

### Exact finite closure from the sole seed 1

For each odd `n<=100000`, in increasing order, accept it only if a listed
rule references an already accepted smaller source. This computes the
complete inductive closure in that universe, with no orbit-search oracle.

| Permitted decreasing dependencies | Certified odd sources, including 1 |
|---|---:|
| Direct `B(n)<n` only | 252 |
| Canonical `Q(n)<n` only | 274 |
| Both direct and canonical choices | 541 |
| Every sibling ancestor of `B(n)` | 640 |

The first two sets are not nested. The third includes both. The fourth
adds 99 further sources, with first new source `483`. Every one of its
640 certificates is compiled and checked against ordinary edges ending
at root `4`. The first ungenerated source remains `7`.

These counts measure exactly these restricted grammars. They neither
count all convergent inputs in the range nor suggest that the remaining
inputs diverge.

## 4. A whole infinite certificate family, with exact clocks

Let positive odd `a` have a supplied root certificate. Choose `j>=0` with

    j=1-a mod3,
    m=S_+^j(a),
    W_j(a)=(8m+1)/3.                                      (4)

Since `S_+(a)=a+1 mod3`, the guard makes `m=1 mod3`, so `W_j(a)` is a
positive odd integer. In fact `W_j(a)=3 mod16` and `W_j(a)>a`.
Direct calculation gives

    3W_j(a)+1=2 S_+^(j+1)(a),
    B_+(W_j(a))=S_+^(j+1)(a),
    B_+^2(W_j(a))=B_+(a).                                 (5)

Thus (J) certifies every such `W_j(a)` from the smaller certified `a`.
Its two odd exponents are exactly

    1, a_+(a)+2j+2.

The extra shortcut cost versus `a -> B_+(a)` is `2j+3`; the extra ordinary
cost is `2j+4`. Integrality, parity, and both clocks are explicit.

With `a=1` and `j=3ell`, this becomes

    N_ell=(32*64^ell-5)/9,
    N_ell -> (16*64^ell-1)/3 -> 1,
    ell=0,1,2,... .                                       (6)

The first sources are `3,227,14563,932067,...`. Every first odd step
increases, so none can belong to the direct-every-step-descent grammar.
All belong to the full sibling grammar. This proves an infinite proper
extension, beyond the finite census. It is a certificate family, not a
cover of the positive odd integers.

The signed version is `W_{s,j}(a)=(8S_s^j(a)+s)/3`, guarded by
`j=1-sa mod3`. It satisfies the same two-step join and the same clock
formula with `a_s(a)`. On the minus sheet it transports the supplied
basin certificate; it does not move all basins to root `1`.

## 5. Exact local coverage, and why inverse ports do not enlarge this family

For `n>1` on the plus sheet,

    Q_+(n)<n iff n=1 mod4 or n=3 mod16.                    (7)

**Proof.** If `n=1 mod4`, then `a_+(n)>=2`, so
`Q_+(n)<=B_+(n)<n`. If `n=3 mod4`, its first exponent is one and
`B_+(n)=(3n+1)/2>n`. A sibling reduction exists precisely when this image
is `5 mod8`, equivalent to `n=3 mod16`. Its first reduced ancestor is
`(3n-1)/8<n`, and subsequent reductions are smaller still. Otherwise
`Q_+(n)=B_+(n)>n`. QED.

There exists *some* smaller sibling ancestor iff the smallest is smaller,
so the full ladder has the same local applicability condition (7), though
it certifies more sources. This covers five of the eight odd classes
modulo 16.

A valid descending inverse port is

    n=5 mod6:  P(n)=(2n-1)/3<n,     B_+(P(n))=n.           (8)

Given a certificate for `P(n)`, delete its forced first odd edge. This is
sound and removes one odd step, one shortcut step, or two ordinary steps.
Together (7) and (8) cover 18 of the 24 odd classes modulo 48. Their
uncovered classes are exactly

    7,15,27,31,39,43 mod48.                                (9)

These are local minima of the stated rank rules, not counterexamples to
Collatz. A smaller endpoint being available does not ensure that it
already has a certificate.

### Forward closure is the decisive obstruction

Let `C` be the infinite family generated from seed `1` by the full sibling
rule (J)+(R). Then

    n in C implies B_+(n) in C.                           (10)

**Proof by induction on a finite certificate derivation.** It holds at
the seed. Suppose the last rule certifies `n` from previously certified
`a<n`, with `B(n)=S^j(a)`. If `j=0`, its successor is `a`, already in
`C`. If `j>=1`, the induction hypothesis gives `B(a) in C`, while

    B(a)<=(3a+1)/2<4a+1<=S^j(a)=B(n).

The direct forward rule therefore certifies `S^j(a)` from `B(a)`.
This gives (10). QED.

Consequently **no finite inverse-orbit port can enlarge this already
closed family**. If an already certified `y` has a finite forward path
to `n`, then repeated use of (10) already places `n` in `C`. This remains
true when the port is discovered by a decreasing inverse formula.

The same proof applies to the raw-plus-canonical family, because its
rules include direct forward descent. It explains the exact no-gain
censuses: adding (8) leaves `541` or `640` unchanged, respectively.
The coordinating lane's sound two-step inverse port

    n=13 mod18: y=(8n-5)/9<n,
    y -> (4n-1)/3 -> n,     exponents (1,2),               (11)

also leaves the full-ladder count at `640`. Locally, (11) removes residues
`31,103,139 mod144` from the lifts of (9), improving applicability from
`3/4` to `19/24` of odd residues; (10) shows why this does not automatically
improve the generated root-certified family.

Thus the next productive rule must supply a new forward/common-future
certificate, a new independently certified seed, or a genuinely stronger
rank argument. More inverse ports alone cannot close this particular gap.

## 6. Hostile controls and the remaining global obligation

On the minus sheet, `Q_-(5)=7` and `Q_-(7)=5`. Its descending inverse
port sends `7` to `5`, where it stops. The signed family gives

    13 -> 19 -> 7 -> 5 -> 7 ...,

because it transports the basin of `a=5`. This confirms that the local
algebra respects distinct basins instead of proving they coincide.

Sibling normalization also does not eliminate long rising prefixes. For
every `L>=1`, put `n0=2^(L+2)-1`. For `1<=j<=L`,

    Q_+^j(n0)=B_+^j(n0)=3^j*2^(L+2-j)-1>n0.              (12)

At every displayed step the remaining power of two is at least four,
so the odd value is `3 mod4`, never `5 mod8`; no sibling compression is
available. The forward exponent is one. This proves that no fixed number
of canonical steps forces numerical descent for every source.

The rank in (R) guarantees termination of a *supplied dependency proof*.
It does not guarantee that each integer admits a rule application leading
to certified smaller dependencies. The residual source `7` already shows
the difference. A universal root decoder still needs a total rule-selection
or well-founded proof showing that every source eventually enters the
certified family. Nothing in (1)--(12) supplies that universal quantifier.

## 7. Reproduction

    python 04-computation/experiments/creative_sibling_20260925.py

Output: [creative_sibling_20260925.out](creative_sibling_20260925.out).
All arithmetic is exact, with explicit exceptions rather than removable
assertions. Universes: both signs, every positive odd input through `20000`,
normalization, valuation guards, descent residues and eight intertwining
clocks; complete increasing-source closure through `100000`; all 640
compiled ordinary root-4 certificates; 20 certified bases and eight
parameters each in (4); the minus controls; and (12) for `L<=20`.
The all-parameter claims have proofs above; the computation is a control,
not an extrapolation. No shared index, theorem reservation, or canon file
was modified by this lane.
