# The recovered 13/4 cycle: endpoint identification and a faithful gain lift

Status: **PROVED** elementary gain-graph propositions below; **INHERITED,
independently audited** valuation-cycle and rotation-system results cited by
path. The gain decoration is a new, explicitly chosen representation. It is
not a newly discovered arithmetic iteration, nor a proof about Collatz or
tournament Hamiltonian counts.

## 1. What the earlier cycle actually was

The exact earlier object is
[the Session 32 valuation-cycle note](../../session32_square_monster_20261023/VALUATION_CYCLE_AND_QUOTIENT.md),
independently checked in
[its audit, section 5](../../session32_square_monster_20261023/monster-audit/AUDIT.md).
For

\[
 N=39,046,339,332,000=2^5 3^7 5^3 7^4\,11\,13^2,
\]

the partial map on the prime support, \(p\mapsto v_p(N)\), has arrows

\[
 13\longrightarrow2\longrightarrow5\longrightarrow3
 \longrightarrow7\longrightarrow4,
 \qquad11\longrightarrow1.
\]

There is no prescribed arrow out of 4. Identifying the terminal 4 with the
input 13 gives **five quotient states**:

\[
 [13=4]\to2\to5\to3\to7\to[13=4].
\]

Thus the recalled display `4 -> 2 -> 5 -> 3 -> 7 -> 13 -> 4` can be used
only after distinguishing the two representatives of the joined state and
the final representative reset. It is not six iterations of the original
valuation map. A representative change contributes no extra base-state
period. Literal integer equality would be false.

The earlier correction also remains essential: the Monster 2A centralizer
is \(2.B\), not \(B\). Using the involution-class index changes the core to
\(N/2\), so \(2\to5\) becomes \(2\to4\). The same quotient then has a
two-cycle \([13=4]\to2\to[13=4]\), with the other arrows feeding it.
The five-cycle is therefore sensitive to the chosen ratio normalization.

## 2. A gain sidecar makes 13/4 an exact invariant

**Proposition (endpoint debt).** Take any positive labels
\(x_0,\ldots,x_r\), initially on an indexed path with distinct positions.
Decorate edge \(i\to i+1\) by

\[
 g_i=x_{i+1}/x_i.
\]

Identify the endpoint positions to make a cycle of length \(r\), retaining
these edge gains. Its multiplicative holonomy is

\[
 H=\prod_{i=0}^{r-1}g_i=x_r/x_0.                 \tag{1}
\]

There exist nonzero vertex values \(z_i\) satisfying
\(z_{i+1}=g_i z_i\) on the quotient cycle iff \(H=1\). Indeed, starting
from \(z_0\ne0\), propagation gives \(z_r=Hz_0\), whereas the identified
endpoint requires \(z_r=z_0\). Conversely this condition closes the
propagated values. The zero assignment always exists but carries no
information and is excluded.

Applied to the recovered chain,

\[
 (g_0,g_1,g_2,g_3,g_4)
   =\left(\frac2{13},\frac52,\frac35,\frac73,\frac47\right),
 \qquad H=\frac4{13}.                           \tag{2}
\]

The representative reset from the terminal label 4 back to the input label
13 has reciprocal gain \(13/4\). This is an exact role for the requested
ratio: it pays the transport defect created by identifying two unequal
numeric endpoints. Appending this reset as a **new** edge creates a
six-edge lifted loop with gain product one; it does not change the original
valuation-map period or supply the missing original arrow out of 4.

This construction retains an output value that the bare set quotient loses.
With the marked starting label \(x_0\), the indexed gains recover every
original label by \(x_i=x_0\prod_{j<i}g_j\). Without that anchor they
retain only ratios and lose a common multiplicative scale. Without the
indexed seam, cycle rotation also loses the original starting position.
For the normalized \(N/2\) chain \(13\to2\to4\), the period is two but
the same endpoint debt \(13/4\) remains. Thus period and endpoint debt are
distinct invariants: neither determines the other.

## 3. Observer changes preserve the debt; a logarithm cut cannot erase it

For an oriented gain graph, change the nonzero local frame at vertex \(i\)
by \(h_i\), so the represented coordinate is \(y_i=h_i z_i\). The new
edge gain is

\[
 g_i'=h_{i+1}g_i/h_i.
\]

On a closed cycle the \(h_i\) telescope, so \(\prod g_i'=H\). This is
the precise observer-frame invariant. Taking \(h_i=x_0/x_i\) on the
first \(r\) path positions makes all but the closing edge have gain one,
and the closing edge have gain \(H\). This gauge is rational when all
labels are rational. It shows that the entire quotient defect can be
stored at one marked seam without removing it.

Over the complex numbers, choose any \(\lambda\) satisfying
\(\lambda^r=H\). Recursively set
\(h_0=1\), \(h_{i+1}=\lambda h_i/g_i\). Equation (1) gives
\(h_r=h_0\), so this is a legitimate cyclic frame change and every
transformed gain is \(\lambda\).

For (2), choose

\[
 \lambda=(4/13)^{1/5}\exp(2\pi i/5).
\]

The transformed transport is a uniformly rotating, contracting spiral.
After five steps it has turned once but has radius multiplied by \(4/13\):
it has returned to the base state, **not** to the transported amplitude.
This is a genuine monodromy model for the newly defined gain lift, not a
claim that prime valuations themselves evolve by a complex spiral.

The interpolation

\[
 z(t)=\exp\!\left[\frac t5\bigl(\log(4/13)+2\pi i\bigr)\right],
 \qquad0\le t\le5,
\]

crosses the negative-real ray at \(t=5/2\). The principal-log cut changes
the displayed imaginary part by \(2\pi i\); a continuous argument keeps
the full turn. Neither convention changes the real displacement
\(\log(4/13)\). In general, all logarithms of \(H>0\) are
\(\log H+2\pi i k\), so branch choices cannot turn \(H\ne1\) into zero
logarithmic displacement. This open spiral has an argument change of
\(2\pi\); do not call it a winding number of a closed curve without
adding a closing path.

No covering map from the original valuation chain to the quotient is
claimed: its endpoint fibre has size two while the other fibres have size
one. The gain lift is extra structure that remedies this loss, not an
automatic covering supplied by the quotient.

**Exact integer-lift stopping corollary.** Start this newly defined gain
transport at a positive integer amplitude \(s\), and require every
intermediate amplitude to remain an integer. It can complete exactly
\(v_{13}(s)\) turns of the five-edge quotient cycle. The successful-edge
clock is therefore exactly \(5v_{13}(s)\), and the next first edge fails
integrality. The normalized two-edge cycle with gains \(2/13,2\) has the
same number of complete turns but clock \(2v_{13}(s)\).

Proof: during turn \(j\ge0\), the five successive amplitudes are
\(s4^j x_i/13^{j+1}\), with \(x_i\in\{2,5,3,7,4\}\). Each \(x_i\)
is coprime to 13, so this turn is integral iff \(13^{j+1}\mid s\).
The two-edge version uses \(x_i\in\{2,4\}\) and has the same divisibility
condition. Thus a finite source-specific valuation pays for precisely the
stated number of turns. No positive integer supports arbitrarily many
integral turns in this gain system. The quotient still cycles forever if
the amplitude is discarded; unrestricted rational amplitudes also continue
forever. This is a stopping theorem for the explicitly defined gain lift,
not for the original valuation map or Collatz. It separates a base-state
cycle, an exact source budget, and the number of actual transport edges.

## 4. What mod 8 and mod 9 retain

The full original prime-support map does not descend modulo 9:
\(11\equiv2\pmod9\), but \(v_{11}(N)=1\not\equiv5=v_2(N)\pmod9\).
Any modulus identifying 13 and 4 divides 9; modulo 3 the inputs 2 and 5
agree, but their images 5 and 3 disagree. Thus no nontrivial such ring
quotient carries the full map. Even the selected five arrows cannot come
from an integer polynomial modulo 9, by the same modulo-3 compatibility
failure at inputs 2 and 5. These are inherited corrections.

For the new gain lift, blindly reducing all fractions modulo 8 or modulo 9
is also invalid: denominators include 2 and 3. There is, however, a lawful
valuation-plus-unit representation. For a prime \(p\), write each nonzero
rational \(g=p^a u\), with \(a\in\mathbb Z\) and numerator/denominator
of \(u\) coprime to \(p\). Then

\[
 \rho_{p,k}(g)=(a,u\bmod p^k)
 \in\mathbb Z\times(\mathbb Z/p^k\mathbb Z)^\times
\]

is a group homomorphism, with addition in the first coordinate and
multiplication in the second. This follows directly from the unique
valuation decomposition, including negative valuations.

For the complete debt:

\[
 \rho_{2,3}(4/13)=(2,5),\qquad
 \rho_{3,2}(4/13)=(0,1).                         \tag{3}
\]

So the mod-8 register with its integer valuation detects this cycle defect;
the mod-9 register does not. The latter equality is possible because
\(4/13\equiv1\pmod9\) after the complete product is formed. It does not
repair the failed valuation-map quotient, nor permit divisions by 3 in
individual ordinary residue calculations.

Neither register, nor their pair, is faithful on all rational gains:
\(73\ne1\) has \(\rho_{2,3}(73)=\rho_{3,2}(73)=(0,1)\).
Every positive integer \(1+72j\) has this same pair. Exact prime-valuation
vectors (all primes retained, together with the sign if negative gains are
allowed), or the exact rational gain, preserve the
remaining scale. Congruence agreement is a guard, not transport closure.

## 5. The recovered geometry has a different, explicit sidecar

[The inherited 5-6-7 rotation-carrier theorem](../../geometry_567_flag_carriers_20261004.md)
uses directed edge incidences (darts), edge reversal \(\alpha\), and a
local cyclic order \(\sigma\). The face permutation is
\(\psi=\sigma\alpha\). For a cubic map whose faces have size \(q\),
with \(d\) darts,

\[
 \chi=d\left(\frac12+\frac13+\frac1q-1\right).
\]

The excess is \(1/30,0,-1/42\) at \(q=5,6,7\), respectively. This is
the precise spherical/Euclidean/hyperbolic boundary. A scalar edge gain is
not a rotation system, and neither a five-cycle nor its labels determine
this geometry.

The inherited hostile is especially relevant to observer frames. One
Heawood graph supports four Paley-equivariant rotation systems. Equal
point/line rotation signs give seven hexagonal faces and genus one;
opposite signs give three 14-gonal faces and genus three. Global reversal
preserves genus; reversing only one of the two local conventions changes
it. The graph, its 8-dimensional cycle space, and its common Paley symmetry
do not select which surface is intended.

Thus the actual shared mechanism is **composition of retained local
transitions**, with different targets: rational gains preserve amplitude
transport, while dart permutations preserve face incidence. Neither can
be substituted for the other. The next useful application to tournament
observer records is to attach the actual transition map to a shared port
and test its composed return, rather than counting repeated printed labels.

A second mod-9 comparison has an equally precise universe boundary.
[THM-085, the all-permutation forward polynomial](../../../../01-canon/theorems/THM-085-f-omega-mod9-universal.md)
concerns \(F_T(\omega)=\sum_{\pi\in S_n}\omega^{\operatorname{asc}_T(\pi)}\),
with \(\omega^2+\omega+1=0\). Its inherited divisibility by 9 for
\(n\ge6\) does not pass to the Hamiltonian-path subset. Every directed
Hamiltonian path has \(n-1\) forward edges, so that restricted sum is
\(h(T)\omega^{n-1}\). A transitive six-vertex tournament has one such
path; the restricted sum is \(\omega^5=-1-\omega\), which is not divisible
by 9 in \(\mathbb Z[\omega]\). All-permutation phase cancellation and
source-specific Hamiltonian admissibility are distinct predicates; retaining
one does not retain the other.

## 6. Controls and boundary

The companion `check_gain_lift.py` uses exact fractions and explicit errors,
not optimization-sensitive assertions. It checks the recovered valuation
graph and both normalizations, their quotient periods and modular
obstructions; telescoping and rational seam gauges for every label word of
length 2 through 6 over `{1,2,3,4,5}`; and the two register homomorphisms on
all ordered pairs of fractions with numerator/denominator in `1..8`.
It includes the positive zero-debt control, the joint-register hostile 73,
and an exhaustive 720-permutation transitive-order-six subset control.
The integer-clock corollary is replayed for all sources `1..512` and for
`13^k b`, `0<=k<=12`, `b in {1,2,3,5,7,11,17,19}`.
Reproduction:

```text
python -X utf8 -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/history/check_gain_lift.py
python -O -X utf8 -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/history/check_gain_lift.py
```

This does not establish a transfer from the Monster-order integer to a
tournament obstruction, a Collatz source, or a curvature regime. It gives
a specific reversible sidecar for the proposed quotient and an exact reason
why a branch reset or a finite congruence register can conceal its debt.
