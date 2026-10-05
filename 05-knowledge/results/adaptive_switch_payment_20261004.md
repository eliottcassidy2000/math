# Payment through an exhausted-pattern switch

**PROVED:** the source-relative payment theorem, exact maximal-run families,
their compiled actual prefixes, and constructed completed families at every
depth. **FINITE-EXACT:** the stated experiments. **OPEN:** coverage of arbitrary
positive sources by a paid controller. Exhaustion, a larger search budget, or
a new anchor is not itself a payment proof.

Artifacts:
[program](../../04-computation/experiments/adaptive_switch_payment_20261004.py)
and [saved output](adaptive_switch_payment_20261004.out).

## 1. Inheritance and the switched object

The closest mechanisms are the
[paid portrait controllers](collatz_paid_portrait_controllers_20261004.md),
especially their original-source payment and anchor-refill sections, and the
[recursive dependency kernel](collatz_recursive_dependency_kernel_20261004.md).
The latter supplies
\[
 H(n)=\frac{729n+669}{1024},\qquad n\equiv155\pmod{2048},
\]
a strictly smaller common-future dependency, with exact repeat fuel
\[
 \Delta(n)=295n-669,\qquad
 H^m\text{ legal}\ \Longleftrightarrow\
 v_2(\Delta(n))\ge10m+1.
\]
The actual word \(12\) supplies
\[
 G(x)=\frac{9x+5}{8},\qquad x\equiv11\pmod{16}.
\]
It increases every positive legal source and repeats exactly
\(\lfloor(v_2(x+5)-1)/3\rfloor\) times.

The canonical hostile is an expanding stage after an already-paid reduction.
Its missing coordinate is the unchanged original source. The corrected near
miss is to call every numerical increase a graph-rank failure: the inherited
rank can decrease even when the integer grows. The useful sidecar is the
odd cofactor remaining when one anchor's fuel has been exhausted.

The live board is **original source / ordered affine carry / binary fuel /
ternary integrality / graph rank / grounded endpoint**. We study
\[
 N\xRightarrow{H^m}c\xrightarrow{G^q}y.
\]
The first arrow denotes common-future dependencies, while the second is an
actual forward word. No arrow type is erased. A completed endpoint proof can
be transported back to \(N\), but an unproved endpoint remains an obligation.

## 2. A small domination proof pays all heights

Let
\[
 L(x)=H^{-1}(x)=\frac{1024x-669}{729},\qquad
 K(x)=G^2(x)=\frac{81x+85}{64}.
\]
Their difference is
\[
 L(x)-K(x)=\frac{6487x-104781}{46656}.
\]
Since \(104781/6487<17\), both increasing maps preserve \([17,\infty)\)
and \(L(x)>K(x)\) there. Induction gives
\[
 L^m(x)>K^m(x)\qquad(m\ge1,\ x\ge17).
\]
This proof does not commute the two maps; it compares their compositions
pointwise using monotonicity.

**Uniform payment theorem.** Suppose \(m\ge1\) kernel repetitions are legal
from \(N\), and \(q\) actual copies of \(12\) are legal from \(c=H^m(N)\),
where \(1\le q\le2m\). Then
\[
 111\le c<G^q(c)\le G^{2m}(c)<L^m(c)=N.
\]
Indeed the inherited kernel gives \(c\ge111\), while \(G\) increases positive
inputs. Thus the switched stage may grow its own input and still pay the
unchanged original source.

For the inherited proper graph rank
\[
 R(1)=(0,0),\qquad
 R(n)=\left(3^{v_2(n-1)}
 \left(\frac{n-1}{2^{v_2(n-1)}}\right)^2,\ v_2(n-1)\right),
\]
payment also gives \(R(y)<R(N)\). The kernel source has \(N\equiv3\pmod4\);
its energy is \(3(N-1)^2/4\). Every positive odd \(y<N\) has energy at most
\(3(y-1)^2/4\), strictly below that value. This rank is inherited from
[collatz_branch_toll_rank_20261004.md](collatz_branch_toll_rank_20261004.md).

The proof is conditional on the actual guards. It does not say that a
controller can always find a legal next pattern, or that every paid child
has already been certified to reach ROOT.

The concurrent [adaptive credit potential](adaptive_credit_potential_20261004.md)
extends this grouped comparison to arbitrary guarded interleavings:
\(H\) earns two credits, \(G\) spends one, and
\((x+5)(9/8)^{\text{credits}}\) decreases at \(H\) and is invariant at \(G\).
That theorem supplies the more general controller payment rule. The
families here supply concrete all-height cases where both repeat patterns
are actually exhausted, as well as explicit rooted exits and sharp hostiles.

### Exact payment outside the uniform allowance

Write
\[
 A=1024/729,\quad \mu=9/8,\quad \alpha=669/295.
\]
Then
\[
 N=A^m(c-\alpha)+\alpha,\qquad
 y=\mu^q(c+5)-5.
\]
The exact ordered composite has positive integer carrier
\((P,Q,B)\), with
\[
 P=3^{6m+2q},\qquad Q=2^{10m+3q}.
\]
Its necessary and sufficient numerical payment condition is
\[
 (Q-P)N>B.
\]
In particular \(P<Q\) is necessary. Equivalently,
\[
 y-N=(\mu^q-A^m)(c+5)+(A^m-1)(\alpha+5).
\]
At \(q=3m\),
\[
 P/Q=(531441/524288)^m>1,
\]
so every such legal switch grows past \(N\). This refutes a uniform
three-block-per-kernel allowance for numerical payment. Section 4 retains
the separate graph-rank boundary.

## 3. Exact families where both patterns are exhausted

Fix \(m\ge1,\ q\ge2,\ e\in\{1,2,3\}\). Let \(t_0\) be the least positive
odd solution to
\[
 295\,2^{3q+e}t_0\equiv2144\pmod{729^m}.
\]
It is unique modulo \(2\cdot729^m\). For every \(s\ge0\), put
\[
 t=t_0+2\cdot729^m s,\qquad
 c=2^{3q+e}t-5,\qquad
 N=\frac{669+1024^m(295c-669)/729^m}{295},
\]
\[
 y=2^e3^{2q}t-5.
\]
All three are positive odd integers. Ternary integrality follows from the
displayed guard; division by 295 is integral because
\(1024\equiv729\pmod{295}\). Positivity and the inverse-kernel identity
follow as in the inherited repetition theorem.

Since \(q\ge2\), the power of 2 in the first term below is at least seven:
\[
 295c-669=295\,2^{3q+e}t-2144,\qquad2144=2^5\cdot67.
\]
Consequently
\[
 v_2(\Delta(N))=10m+5,\qquad v_2(\Delta(c))=5,
 \qquad v_2(c+5)=3q+e.
\]
The \(H\) run is exactly maximal at \(m\), and the subsequent \(G\) run is
exactly maximal at \(q\). The three values of \(e\) exhaust the possible
remaining binary fuel for a maximal \(G\) run. Conversely every positive
source with exactly this grouped pair of maximal runs, with \(q\ge2\),
has this parameterization.

For fixed \(m,q,e\), these are affine progressions. The source period is
\[
 2^{10m+3q+e+1};
\]
the terminal and final periods are respectively
\[
 2^{3q+e+1}3^{6m},\qquad 2^{e+1}3^{6m+2q}.
\]
Thus \(q=2m\) gives an explicit family of fully exhausted switches paid at
every parameter height. The three residual classes occupy seven eighths of
the full \(H^mG^{2m}\) source cylinder; the remaining eighth permits more
than \(2m\) copies of \(G\), beyond the uniform allowance.

These families lie inside the already-paid \(H\) input cylinder. The advance
is an exact switch/closure mechanism and an actual forward certificate,
not additional source coverage beyond that inherited cylinder.

### The compiled actual word

Use the inherited
\[
 v=(1,2,1,1,1,2),\qquad W=(3,2,1,1,1,2).
\]
Substituting the first edge of the child's actual \(12\) word gives the
source prefix
\[
 w_{m,q}=v\,W^{m-1}(3,2)(1,2)^{q-1}.
\]
It has \(6m+2q\) odd edges, cost \(10m+3q\), and reaches \(y\) exactly.
The source and endpoint identities are checked using both the composed
operators and literal valuation replay. There is no early ROOT: the
switched states increase from \(c\ge111\), and the inherited splice proof
retains the source's earlier nonroot states.

At \(q=2m\) this is a genuine actual descent after \(10m\) odd edges with
halving cost \(16m\). The dependency derivation compresses how this word was
obtained; it does not replace its source guard or first-hit conditions.

## 4. Positive and hostile boundaries

The least residual-one row at \(m=1,q=2\) gives
\[
 N=241819,\qquad c=172155,\qquad y=217885.
\]
Thus \(c<y<N\): insisting that the new pattern decrease its own input
would reject a valid payment against the original source.

At \(m=1,q=3,e=1\), the corresponding row begins
\[
 110747\xRightarrow H78843\xrightarrow{G^3}112261.
\]
The integer grows beyond its original value, but the proper graph rank
decreases because the final binary precision increases. Numerical and
graph-rank payment must not be conflated.

For a genuine failure of both payments, use \(e=2\):
\[
 1159323\xRightarrow H825339\xrightarrow{G^3}1175143,
\]
\[
 R(1159323)=(1008020624763,1),\qquad
 R(1175143)=(1035719040123,1).
\]
More generally, every \(q=3m\) row with \(e=2\) or \(e=3\) has both source
and endpoint congruent to 3 modulo 4. Its numerical growth therefore
strictly increases the graph rank at every height. These are failures of
this proposed payment rule, not failures of Collatz convergence; the earlier
kernel dependency \(c<N\) remains available.

## 5. A constructed ROOT exit at every depth

For \(q=2m,e=1\), use the least odd address \(t_0\) from section 3. Choose
\(b=2a\ge4\) with
\[
 4^a\equiv2\cdot3^{4m+1}t_0-14
       \pmod{3^{10m+1}}.
\]
The target is 1 modulo 3. The elementary principal-unit order of 4 gives
one exponent class modulo \(3^{10m}\), so there are infinitely many
positive choices. The script lifts the exponent one ternary digit at a time.

Set
\[
 t=\frac{2^b+14}{2\cdot3^{4m+1}}.
\]
This is an integer, is odd, and lies in the required class
\(t_0\pmod{2\cdot729^m}\). To verify this, the exponent congruence says
\[
 2^b+14=2\cdot3^{4m+1}t_0+3^{10m+1}z.
\]
Both first terms are even and the modulus is odd, so \(z\) is even;
division shows \(t=t_0+729^m(z/2)\). Independently \(b\ge4\) gives
\(v_2(2^b+14)=1\), so \(t\) is odd and \(z/2\) must be even. Positivity
then gives a nonnegative parameter on the canonical positive row.

Its switched endpoint is
\[
 y=2\cdot3^{4m}t-5=\frac{2^b-1}{3},
\]
which has the one-edge first-hit route \((b)\) to ROOT. Hence
\[
 w_{m,2m}(b)
\]
is a completed first-hit source certificate with ranks
\[
 \text{odd rank}=10m+1,\qquad
 \text{ordinary rank}=26m+b+1.
\]
This is a completed constructed family at every \(m\), with no orbit-search
premise. It is still not a certificate for arbitrary free row parameters.

There is a precise comparison with the previously constructed completed
family in section 4 of the
[recursive dependency kernel note](collatz_recursive_dependency_kernel_20261004.md).
Those sources have exact fuel \(10j+1\), whereas every source constructed
here has fuel \(10m+5\). The two completed families are therefore disjoint
for every \(j,m\ge1\). This adds rooted members relative to that specific
completed-family bank; it is not a disjointness claim against all previous
certificates or an enlargement of the underlying paid \(H\) input guard.

At \(m=1\), one solution is \(b=962\). The 961-bit source was expanded
and independently replayed; its ranks are \((11,989)\). For \(m=2\) the
least displayed phase has \(b=668613680\); the script deliberately does
not expand that source. Higher controls use full-denominator modular
division and an independent inverse-letter reader, with exact ROOT-codec
rank checks.

## 6. Interface, accounting, and reproduction

A Switch record contains only the supplied source and the two chosen
counts. Its verifier checks exact source guards, the common-future terminal,
the actual \(12\) guard, the ordered affine summary, and payment against
the unchanged source. It never infers ROOT status from exhaustion.
The certificate adapter consumes a supplied exact endpoint word and
exports the inherited first-hit ROOT grammar with source identity intact.

The separate test helper that discovers endpoint proofs by forward
observation is not part of either verifier. Four such controls consume
317 explicitly counted odd observations. The constructed families in
section 5 require none.

Reproduce:

    python -B 04-computation/experiments/adaptive_switch_payment_20261004.py
    python -B -O 04-computation/experiments/adaptive_switch_payment_20261004.py

Finite universe: \(m=1,\ldots,12\), \(q=2m\) or \(3m\), residuals
\(e=1,2,3\), and row parameters \(0,1,7,31\). All 288 sources are checked by
literal source and switched-prefix replay. There are 144 paid controls
and 144 failures of numerical payment; the 96 designated residual-two/three
hostiles also strictly increase graph rank. There is one expanded constructed
ROOT family control, plus 36 symbolic plans at \(m=1,\ldots,12\), exponent
lifts \(0,1,7\), with 144 independent modular-reader comparisons.
Malformed types, corrupted sources, a wrong endpoint proof, an invalid
exponent phase, and a too-small expansion cap are rejected.

All checks remain active under optimization. The compression retains the
two counts and ordered carriers; actual integer sizes and optional word
export costs are not hidden. No prior controller, global index, or existing
artifact is modified by this package.
