# Information-preserving observers still need their source and action markers

2026-10-14. **PROVED:** the finite Fisher identities, verifier-relative
equivalence, cycle-phase obstruction and cocycle repair below.
**FINITE-EXACT:** the declared exhaustive finite controls and literal ROOT
receipt. **CITED:** the classical Fisher and exceptional symmetric-group
interfaces. **OPEN:** universal Collatz convergence. No statistical model
here supplies a missing ROOT certificate.

## 1. Inheritance and the object being compressed

The closest mechanism is the separation of cylinder mass from affine carry
in [valuation information spectrum](valuation_information_spectrum_20261005.md).
The canonical hostile is two valuation words with the same length and
total cost but different native source cells. The corrected near miss in
this session is assigning a normalized cocycle around a cycle without
checking its holonomy. The retained sidecars are the original source,
the representation on which a group acts, and the phase within a cycle.

Our board is source / score / verifier / action / phase / terminal.
The anchor is an exact criterion for when a compression preserves a given
finite verification task. The niche is the outer automorphism of $S_6$,
which preserves a group without preserving its usual point action. The
boundary-function reading supplies a separate approach-marker hostile.

For a connection, the source is a finite marked record $x$, the target is
a statistic $T(x)$, and the preserved predicate must be specified. A
parameter score, a Boolean verifier and a pair relation are different
predicates. The examples below keep their maps explicit.

## 2. Fisher loss is exactly conditional score variance

Let $\Omega$ be finite, $p_\theta(x)>0$ a continuously differentiable
probability law near a fixed parameter, and $T:\Omega\to\mathcal T$ a
parameter-independent statistic. Put

\[
s_\theta(x)=\partial_\theta\log p_\theta(x),\qquad
I_X(\theta)=\sum_xp_\theta(x)s_\theta(x)^2.
\]

Differentiating each finite fibre sum gives

\[
s_T(t)=\mathbb E_\theta[s_\theta(X)\mid T=t],\qquad
\boxed{I_X-I_T=\mathbb E_\theta\operatorname{Var}(s_\theta\mid T).}       \tag{1}
\]

Thus equality holds at this parameter exactly when the score is constant
on every fibre. This is the elementary finite form of Pollard,
[*A note on insufficiency and the preservation of Fisher information*,
Theorem 7, p. 6](https://arxiv.org/pdf/1107.3797).
The proof here uses positive finite support; no claim about changing
support, singular models, or global sufficiency follows from a pointwise
equality in an arbitrary model.

There is a constructive verifier-specific form. Fix a positive base law
$p$ and a supplied Boolean predicate $V$. Define the exponential tilt

\[
p_\theta(x)=\frac{p(x)e^{\theta V(x)}}{\sum_yp(y)e^{\theta V(y)}}.
\]

Its score at zero is $V-\mathbb E_pV$. Consequently

\[
\boxed{I_X(0)=I_T(0)\iff
 V\text{ is constant on every }T\text{-fibre}
 \iff V=\bar V\circ T\text{ for a Boolean }\bar V.}                    \tag{2}
\]

Equation (2) audits compression against a specified task. It does not
compute an unknown $V$: the verifier or its values are premises. For
every perturbation to preserve Fisher information, $T$ must be injective.
If two points merge, choosing a Boolean predicate that separates them
violates (2); the converse follows from invertibility on the image.
For a relation $R(x,y)$, the corresponding exact requirement is constancy
on the fibres of $T\times T$. A one-object score need not test that.

In particular, $T=V$ is a sufficient one-bit statistic for this tilt.
Storing an **asserted** verdict bit is not an independently checkable proof
certificate. The finite information criterion assumes an authenticated
source model and evaluator; effective verification still needs a checked
map, proof witness or replay. Constancy on fibres establishes that a
well-defined $\bar V$ exists, not a complexity bound for computing it.

### An actual Collatz guard lost by a sufficient statistic

Let $w=(a,b)\in\{1,2,3\}^2$, with truncated geometric law

\[
p_z(w)=\frac{z^{a+b-2}}{(1+z+z^2)^2},\quad z>0.
\]

The cost $A=a+b$ is sufficient for $z$: conditionally on $A$, every order
has the same probability, independently of $z$. At $z=1/2$, both Fisher
informations are $208/49$. Nevertheless the odd Collatz carriers are

\[
U_{(1,3)}(n)=\frac{9n+5}{16},\quad n\equiv19\pmod{32};\qquad
U_{(3,1)}(n)=\frac{9n+11}{16},\quad n\equiv29\pmod{32}.                \tag{3}
\]

Take $V(w)$ to mean that $w$ is the actual strict two-step word at source
19. Within the declared nine-word universe this is true only for $(1,3)$.
For its verdict tilt, the full, cost-compressed and lost informations are

\[
\frac{180}{2401},\qquad \frac{148}{7203},\qquad \frac8{147}.
\]

These exact values distinguish statistical sufficiency for $z$ from
verification sufficiency for the original source. The cheap repair is
to retain the ordered word and its source, or an independently audited
equivalent guard representation. The native odd cell for any word with
carrier $(Pn+B)/Q$ is $n\equiv(Q-B)P^{-1}\pmod{2Q}$; strict replay also
rejects every earlier ROOT hit.

The positive receipt is literal:
$19\xrightarrow{(1,3)}11\xrightarrow{(1,2,3,4)}1$.
`prepend` consumes the supplied first-hit ROOT word at exactly 11 and
exports $(1,3,1,2,3,4)$ at 19. A valid certificate at 17 is rejected at
the wrong port. This is the same source-owned interface as
[marked completion](collatz_marked_completion_20261007c.md), not a
probabilistic proof of a new basin.

## 3. The outer automorphism preserves the group, not its point action

On six labelled points, a **syntheme** is a perfect matching of $K_6$;
there are 15. A **pentad** partitions the 15 pairs into five synthemes;
there are six. Relabelling points acts on these six pentads. After fixing
their order, this gives the classical outer automorphism
$\Phi:S_6\to S_6$.

The primary reference is Howard--Millson--Snowden--Vakil,
[*A description of the outer automorphism of $S_6$, and the invariants of
six points in projective space*, §1.4](https://arxiv.org/html/0710.5916v1).
The user's [automorphism overview](https://en.wikipedia.org/wiki/Automorphisms_of_the_symmetric_and_alternating_groups)
was an entry point; the construction is not claimed as new.

The checker enumerates all 15 matchings, all six pentads and all 720 point
permutations. Their pentad images are all distinct, and composition is
checked along all five adjacent-transposition generators. This certifies
the finite group map without identifying the two underlying actions.
In the chosen orders, $(0\ 1)$ becomes $(0\ 1)(2\ 3)(4\ 5)$.
It fixes four original points but no pentad. Therefore no bijection
between the two six-element carriers intertwines these actions: such a
bijection would preserve the number of fixed points.

This is a precise pair-observer hostile. A transposition identifying one
swapped point-pair is transported to three swapped pentad-pairs. Group
order, multiplication and an invertibly transported probability model
remain intact; the original vertex-pair meaning does not. In particular,
all Fisher information of any finite model on group elements survives
the bijection $\Phi$, even though the predicate “fixes this original
vertex” requires the original action or its transported predicate.
`ActionReceipt` stores and checks both permutations. No tournament or
plane-graph relation is asserted merely because both carriers have six
elements.

## 4. Doubling cycles require phase, including primitive-root cases

The finite arithmetic interface inherited from
[THM-4631, level factorisation](../../01-canon/theorems/THM-4631-level-factorisation-of-the-compressed-twisted-transfer-operator-of-the-px-plus-1-maps-fast-modes-are-the-doubling-cycles-on-the-units-mod-p-i.md)
is multiplication by $2^{-1}$ on $(\mathbb Z/p^i\mathbb Z)^\times$.
It has the same cycles, with reversed orientation, as doubling. Their
length is $L_i=\operatorname{ord}_{p^i}(2)$ and their number is
$\varphi(p^i)/L_i$. This section uses only that finite action; it does
not treat the twisted complex operator as a probability law.

Let $U$ be the units modulo any odd $m>1$, and let $T$ record only the
doubling-cycle index. Choose two successive points in one cycle, put
$f=1,-1$ there and zero elsewhere, and set
$p_\theta(a)=(1+\theta f(a))/|U|$ for $|\theta|<1$.
Every cycle mass is constant, so

\[
I_T(0)=0,\qquad I_U(0)=2/|U|>0.                       \tag{4}
\]

Retaining cycle index **and phase** is bijective and restores equality.
For $m=7$ the cycles are $(1,2,4)$ and $(3,6,5)$, and (4) is $1/3$.
For $m=3$ there is only one cycle, but the lost information is 1.
Thus a primitive-root condition removes cycle labels, not the phase
obligation. No Artin-density statement or integer ROOT conclusion is
needed for this example.

### Local correction: a cycle cannot be normalized without holonomy

For a directed cycle $v_i\to v_{i+1}$ with nonzero edge weights $w_i$,
a nonzero solution to $\chi(v_i)=w_i\chi(v_{i+1})$ exists iff
$\prod_iw_i=1$. Necessity follows by multiplication; sufficiency follows
by assigning one nonzero value and propagating once around the cycle.
On an acyclic forward-orbit component the propagation is consistent.

The earlier wording in THM-4631 statement 1 omitted this condition.
For shortcut Collatz $1\to2\to1$, its character weights have product
$\zeta^{2\alpha+\beta}$. Taking $m=2,\alpha=0,\beta=1$ gives $-1$;
only the zero cocycle satisfies both retained edges. First-hit ROOT
characters instead prescribe $\chi(1)=1$ and impose recursion only at
$n\ne1$, on the ROOT-killed graph. Finite-horizon characters have a
horizon index and are not automatically stationary cocycles.

The authorized correction was made in THM-4631 statement 1 and its
[current result discussion](graph_theory_as_common_carrier_20261013.md).
It changes this existence/terminal interpretation only; the stated
operator factorization was neither changed nor broadly re-audited here.
The checker tests all binary charge assignments on cycles of length
1 through 7 at character moduli 1 through 7, including the killed-root
positive control.

## 5. Certainty and boundary limits do not replace a terminal proof

A constant family $\delta_n$ has zero Fisher information and zero
entropy: the current state is known. This does not determine a future
event without its transition law. The two finite maps in the checker
have the same certain initial state, with one reaching the marked root
and the other fixed away from it. This is a logical hostile, not a
claimed counterexample to Collatz. Parameter-varying singular Dirac
families are outside the positive finite-model hypothesis of §2.

The user-named
[mathematical archive](https://github.com/micah541/kaczynski-mathematical-works/tree/main)
is a modern reconstruction, as its
[README](https://raw.githubusercontent.com/micah541/kaczynski-mathematical-works/main/README.md)
and [errata](https://raw.githubusercontent.com/micah541/kaczynski-mathematical-works/main/ERRATA.md)
explain. Its boundary and harmonic sections were read, including the
corresponding full blueprint pages. No Lean build or axiom audit was run,
so the archive's reported formal status is not a verification by this
packet. The original primary paper, Kaczynski,
[*Boundary functions for bounded harmonic functions*, Trans. AMS 137
(1969), pp. 203--209, especially p. 203](https://www.thetedkarchive.com/uploads/tedkarchive/t-k-ted-kaczynski-the-mathematical-work-of-ted-kac-4.pdf),
distinguishes an asymptotic limit along an arc from an unrestricted
boundary limit.

The useful elementary transfer is the missing approach marker. On the
upper half-plane, the bounded continuous function
$f(x,y)=x^2/(x^2+y^2)$ has along $x=cy$, $y\downarrow0$, limit
$c^2/(1+c^2)$. The same boundary point gives 0, $1/2$ or $4/5$ along
three marked rays. This example is **not harmonic** and invokes no
harmonic representation theorem. The boundary point by itself does not
encode the approach. No analytic extension of Collatz or inference from
these arc limits to ROOT is asserted.

Likewise the companion
[continued-fraction code](continued_fraction_collatz_codes_20261014a.md)
can retain a whole trajectory while its irrationality exponent is always
two. Losslessness of the full carrier does not make every scalar of that
carrier sufficient for the target predicate.
For the [parity code](collatz_parity_cf_boundary_20261014a.md), the same
constant exponent statistic has zero Fisher information for every
nonconstant finite verdict tilt, despite injectivity of the full code.

## 6. Exact controls and the next obligation

Run [the standalone checker](../../04-computation/experiments/information_observers_collatz_20261014a.py)
using `python -B -X utf8` and `python -B -O -X utf8`.
The [saved JSON](information_observers_collatz_20261014a.json) records
8,181 always-active checks: 1,954 combinations of all Boolean predicates
and all partitions on one through five points; the nine-word Collatz
model and strict source-owned splice; all 720 outer actions and their
five generator edges; six unit-group phase examples; 1,778 exact
cycle-holonomy tests; 80 marked boundary approaches; and typed hostiles.
All probabilities and tangents in the checker are exact rational values.
Normal and optimized output agree byte-for-byte; the LF-normalized SHA256
is `641a721e226f4984d3db2c06a37540c6911e2c7cd2377706f81981b2539abed8`.

The next concrete obligation is to specify the actual certificate
verifier $V$, then prove a proposed representation makes $V$ constant on
its fibres, retaining the source and ROOT boundary. Equation (2) gives
a finite diagnostic for that task. It cannot manufacture the verifier's
missing positive values. The cycle and $S_6$ examples show two distinct
repairs: retain phase within an orbit, and retain which action interprets
a group element. Neither repair alone proves global convergence.
