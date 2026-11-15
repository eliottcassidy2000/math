# Exact windows: the finite-state boundary and two-radix criterion

[Session 25](../session25_spectrum_wiring_20261021/session25-spectrum/PROOF.md)
adds an apex to a transitive spine whose contact word repeats $1^L0^L$
exactly $t$ times. For $Q=2^L$, its actual Hamiltonian count is
$$
h(L,t)=1+\frac{Q^2(4^{Lt}-1)}{(Q+1)^2}+\frac{2t}{Q+1},
\qquad L,t\ge1.
$$
The inherited Session 25 construction supplies arbitrarily large actual
counts in each parity-admissible residue class. Calibrating its modulus
changes $L$; this is not finite-state decoding of the absolute target.

## Exact count and an unconditional nonregularity theorem

The map $(L,t)\mapsto h(L,t)$ is injective. Its counting function is
$$
\#\{h(L,t)\le X\}=Y\log Y+(2\gamma-1)Y+O(\sqrt Y),\qquad Y=\log_4 X.
$$
In particular the canonical digit language of this family is not regular
in any integer base.

*Injectivity.* Put $D=Lt$. The bounds
$(4/9)4^D\le h(L,t)\le4^D+2D/3$ separate successive shells $D$.
For fixed $D$, set $F(Q)=Q^2/(Q+1)^2$. For $q\ge2$,
$$
F(2q)-F(q)=\frac{q^2(4q+3)}{(2q+1)^2(q+1)^2}>\frac1{4q}.
$$
Indeed $2q+1\le(5/2)q$, $q+1\le(3/2)q$, and
$4/[(5/2)^2(3/2)^2]>1/4$. For two choices $L_1<L_2$ of divisors of
$D$, the bulk increases by more than
$(4^D-1)/(4\cdot 2^{D-1})=2^{D-1}-2^{-D-1}$ because
$2^{L_1}\le2^{D-1}$. This exceeds $2D/3$ for $D\ge2$
(the left increases by more than one per step). The correction term can
decrease by no more than $2D/3$. Thus $h(L,D/L)$ increases with $L$ for
fixed $D$. For $D=1$ there is one pair.

*Counting.* Every shell $D\le\lfloor Y\rfloor-1$ is fully included
and none with $D\ge\lfloor Y\rfloor+2$ is included. Each shell
contains exactly $\tau(D)\le2\sqrt D$ values. Thus the count differs
from $\sum_{d\le Y}\tau(d)$ by $O(\sqrt Y)$. The hyperbola split gives
$\sum_{d\le Y}\tau(d)=2\sum_{j\le\sqrt Y}\lfloor Y/j\rfloor
-\lfloor\sqrt Y\rfloor^2=Y\log Y+(2\gamma-1)Y+O(\sqrt Y)$.

*Nonregularity.* A nonempty regular language has either exponential
cumulative growth in word length or cumulative growth $\Theta(n^k)$
for an integer $k\ge0$. Trim a deterministic automaton to states on
successful paths. If a recurrent strongly connected component has a
vertex with two internal outgoing transitions, there are two return
loops with different first symbols. Repeat the loops to a common
length; arbitrary concatenations exhibit exponential growth. Otherwise
the recurrent components are simple cycles. A successful path has one
of finitely many skeletons with a repetition count for each traversed
cycle. If the largest number of cycles on one successful path is $k$,
there are at most $O(n^k)$ words of length at most $n$; fixing a
skeleton with $k$ cycles gives $\Omega(n^k)$ distinct such words by
varying the repetition counts independently. In base $b$, length at
most $n$ means value below $b^n$, so this particular language grows as
$\Theta(n\log n)$, neither exponential nor an integer power.

This concerns exactly the calibrated family, not the whole one-apex
spectrum. It rules out any finite automaton recognizing its **exact
outputs**, even in base two.

## A sufficient two-radix criterion

Let $A$ be a set of odd positive integers with arbitrarily large elements
in each parity-admissible progression. This includes the full strong
Hamiltonian-count spectrum by Session 25. Then the following are equivalent:

1. $A$ contains all sufficiently large odds.
2. $A$ is eventually periodic.
3. The same exact set $A$ has regular canonical digit languages in two
   multiplicatively independent integer bases.

The first property gives the second and third by recognizing parity and
finitely many exceptions. The third implies the second by
[Cobham's theorem](https://research.ibm.com/publications/on-the-base-dependence-of-sets-of-numbers-recognizable-by-finite-automata).
The second implies the first because any late missing odd makes a
whole odd progression missing after the threshold (double the period
if needed), contradicting the progression property.

Two different sufficient subclasses, one binary and one ternary,
do not suffice. One physical subclass with the progression property
and two exact recognizers would suffice. In particular it is not
necessary to recognize all tournament constructions. However, the
existing Session 25 compiler itself cannot be used as either exact
regular language by the preceding theorem.

An explicit finite gate follows from a *certified* $d$-state DFA in
canonical base $b$ for the finite missing-odd language: an accepted
word of length at least $d$ repeats a state; pumping the loop preserves
its leading nonzero digit and gives arbitrarily long accepted words.
Thus every missing odd is below $b^d$. Qualitative cofiniteness alone
does not give that effective gate.

One radix plus progression coverage and ordinal-product closure is
insufficient: odd positives with $7^j$ omitted for every $j\ge1$ are
multiplicatively closed, meet all admissible progressions arbitrarily
high (powers of 7 are logarithmically sparse), and form a base-7
regular set with infinitely many holes.

## What the two real holes teach

The [all-order 7/21 proof](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md)
reduces candidate witnesses to strong cores by SCC multiplication;
the Moon–Busch floor caps the orders; the finite strong spectra
finish the proof. Its positive dual would put **all candidate holes**
in one finite band. The cap on witness order for one fixed count is
not a cap on every missing count. Nor can any fixed modulus distinguish
7 or 21 from all realizable counts.

The repository's [Collatz/Mahler analogy](../procgen_continuous_20260925_natural_boundary_mahler_harmonic.md)
makes a useful distinction: two independent exact recognizers force
global rigidity, whereas mixed transfer equations or finite residue clocks
do not supply those recognizers. A binary/ternary tournament recursion
needs bounded absolute carry, physical witness ownership, and exact
compilation in both directions. That is a precise obligation, not a
claim that the corresponding Collatz problems transfer.
