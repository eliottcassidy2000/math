# Coverage through earlier reroutes, guarded group actions, and retained arithmetic memory

2026-10-04. **PROVED:** the scoped constructions and obstructions linked
below. **FINITE-EXACT:** their explicitly declared controls and certificate
demonstrations. **OPEN:** guaranteed completion for every positive integer.
A smaller dependency is an induction step; it is not a supplied proof that
every free-parameter child reaches 1.

The main coverage gain is to choose where a construction starts, rather
than require every construction to act at the same reset. This produces
new all-height cancellations inside the current critical set. Burnside's
theorem and the three quadratic maps lead to two exact supporting ideas:
restore the integer guards behind a finite group action, and retain the
unbounded precision consumed when a word is repeated.

## 1. Direct coverage progress

The [early-reroute theorem](collatz_early_reroute_20261004.md) starts a
ternary sibling construction at the original source, or at an actual
initial valuation-one checkpoint, then applies the half-child reset
switch. It compares only the final child with the immutable original
source. Its exact least inverse clock is

\[
 D(k)=\min\{d\ge1:3^d>2^{d+2k-1}\}.
\]

The earlier sibling clock paid a factor \(4^k\); the final switch reduces
the factor to \(4^k/2\). A proof gives all necessary position and ternary
guards, a necessary-and-sufficient applicability test within this entire
grammar, and a finite source-height bound on the search. The selector is
not a guessed depth cutoff.

Two concrete families, for every integer \(t\ge0\), are

| Original source n | Smaller child h | Actual valuation words to a common future |
|---|---|---|
| \(89211+93312t\) | \(62655+65536t\) | source \((1)\), child \((1^5,2,3)\) |
| \(2424699+2519424t\) | \(2018303+2097152t\) | source \((1,2)\), child \((1^9,2,6)\) |

Every source is27 modulo96, a local minimum of the incoming
[proper graph rank](collatz_branch_toll_rank_20261004.md). Both source
and child have \(v_2(n-1)=1\), so the strict numerical inequality also
strictly decreases that rank. These are actual cancellations within its
residual problem.

The overlap analysis is explicit. Two parameter refinements of the first
row each have at least \(77/78-\varepsilon\) of their members outside
the named earlier binary, ternary, small-predecessor and composed banks,
where
\[
 \varepsilon=\frac{3^{-55}}{1-3^{-30}}.
\]
Their combined added density among odd integers is at least
\(2(77/78-\varepsilon)/(64\cdot3^7)\). The second row has at least
\(25/26\) of its members outside the entire strengthened original-source
ternary bank. These are family-relative and bank-relative coverage
statements; they do not measure completed root proofs or compare against
every earlier repository operation.

Modulo19 remains useful without demanding that every new rule preserve
the residue identically. For the two rows the exact maps are respectively
\[
 h\equiv8n+2\pmod{19},\qquad h\equiv13n+17\pmod{19}.
\]
They are permutations, and all19 source phases occur. Keeping the rule
with the residue retains the phase information while allowing more
construction choices.

## 2. Reroute when a whole infinite family resists

The previous after-reset construction missed every odd-exponent Mersenne
input \(2^p-1\). The new rules handle the following whole exponent
progressions:

\[
\begin{array}{c|c}
 p&\text{scope}\\\hline
5+6t&\text{inherited small inverse rule}\\
75379+354294t&\text{new sparse family inside }1\pmod6\\
36969036897+62762119218t&\text{new sparse family inside }3\pmod6
\end{array}
\]

These are not all exponents in the latter two classes. The large final
family is checked symbolically through its exact ternary exponent guard;
no billions-bit integer is expanded.

The opposite prime tower gives a sharp new hostile. The entire early
grammar misses every \(3^a\), \(a\ge1\), and more generally it misses
\[
 n=3^a u,\quad 3\nmid u,\qquad 9^a\ge32u^3.
\]
This is a source-size versus ternary-valuation obstruction at all indices
and depths, not a finite failed search.

Dividing this hostile among other constructions immediately reduces it:

| Exponent of \(3^a\) | Available reduction |
|---|---|
| a even | Immediate actual odd descent |
| \(a\equiv1\pmod4\) | Inherited reset-at-least-three half-child switch |
| \(a\equiv7\pmod8\) | Actual three-step descent |
| \(a\equiv3\pmod8\) | Residual of these three tests; other rules may apply |

An [arbitrary-checkpoint compiler](collatz_checkpoint_reroute_20261004.md)
then breaks into that last domain. Unlike the initial-run construction,
it retains the ordered carry of an arbitrary actual source prefix. Its
inverse depth may be smaller than that prefix's length; imposing the
full prefix-depth ternary guard would incorrectly discard valid moves.

The first successful example initially gave a modulus of 131072. Cutting
its redundant common-future suffix and allowing the last valuation to
vary reduces the modulus by a factor of 64. For every integer \(t\ge0\),

\[
 \boxed{n=155+2048t,\qquad h=111+1458t<n.}
\]

The first six actual valuations of n are \((1,2,1,1,1,2)\), leading to
\[
 x=\frac{729n+925}{256}=4h+1.
\]
Since \(3(4h+1)+1=4(3h+1)\), the odd Collatz values \(U(x)\) and
\(U(h)\) are equal. This uses the inherited sibling identity at a new
checkpoint. It needs no conjectured convergence and gives a strict
numerical induction step. It also lowers the graph rank from section 1:
the source is 3 modulo 4, so its energy is \(3(n-1)^2/4\), whereas
every smaller positive odd child has energy at most \(3(h-1)^2/4\).
The even-t half has no ordinary descent through the seven-step join;
the odd-t half already descends at that seventh step. The unified guard
dispatches both, and the even-t half is the nontrivial common-future gain.

In particular, \(\operatorname{ord}_{2048}(3)=512\) and
\(3^{483}\equiv155\pmod{2048}\), so the entire exponent progression
\[
 \boxed{a\equiv483\pmod{512}}
\]
has this reduction for \(n=3^a\). It occupies one of the 64 exponent
classes modulo 512 inside \(a\equiv3\pmod8\). This is coverage of that
exponent family, not a claim to have rooted every resulting child.

The incoming [paid portrait controllers](collatz_paid_portrait_controllers_20261004.md)
then supplied a complementary unbounded family. Their power-of-three
rules follow \((1,2)^q\), \(q\ge2\), and a sufficiently large final reset.
Every covered exponent is 11 modulo 32; every exponent in our new family
is 3 modulo 32. The two domains are therefore disjoint. If \(\delta\)
denotes that note's proved relative density within \(a\equiv3\pmod8\),
their combined reduction coverage is exactly
\[
 \delta+\frac1{64}\ \approx\ 0.08258927719932.
\]
Thus this explicitly compared union handles about 8.2589 percent of the
residual **exponent class**. It is not an integer-density or a root-proof
percentage. The incoming unbounded repetition parameter also explains
why its controller lies outside the bounded-depth obstruction below.

The complete per-checkpoint search has a finite height bound. The number
of checkpoints needed for an arbitrary source still has no proved bound.
The small hostile inputs remain useful: the specified checkpoint grammar
finds no reroute before ordinary descent for 7, 27 or 703. Ordinary descent
does close these particular inputs, so a failed grammar test is not an
unresolved Collatz trajectory.

## 3. Burnside's idea becomes an explicit two-layer calculation

Burnside's finite-group theorem concerns solvability for order
\(p^a q^b\). For the actual affine group generated by the inverse
shortcut operations modulo \(m>1\), \(\gcd(m,6)=1\), we can prove more
directly:

\[
 G_m=(\mathbb Z/m\mathbb Z,+)\rtimes\langle2,3\rangle_m,
 \qquad G_m'=\mathbb Z/m\mathbb Z,
 \qquad G_m''=1.
\]

The layers are multiplier and translation. This group has derived
length two even when its order has more than two prime divisors; at19
its order is342. The [guarded affine-lift theorem](collatz_affine_guarded_lifts_20261004.md)
then restores what this solvable group forgets.

For every supplied positive odd hub \(u\) coprime to3, every
\(g\in G_m\) has a legal inverse Collatz representative from u. It
constructs infinitely many odd ancestors above u, with the entire modular
affine action g, exact integer guards, and a first-hit route to u. A
supplied root certificate for u splices to each ancestor. A finite Cayley
walk chooses the operation pattern; an unbounded ternary discrete
logarithm supplies the initial doubling exponent. Large exponents and
endpoints are stored symbolically with independently checked residue
readers.

This proves existence of certified representatives for all modular
targets. It does not identify an arbitrary prescribed integer with one
of those representatives, nor make them smaller than the supplied hub.
Both the construction and its quantifier boundary are essential.

The audit also corrected a printed CRT step in the literature. At95,
\(3=2^{31}\) and \(\operatorname{ord}_{95}(2)=36\), so the global
multiplier group has36 elements, although the local multiplier groups
have orders4 and18. Their product72 incorrectly treats the two shared
generator exponents as independent. The exact retained correlation is
\(\chi_5(a)=\chi_{19}(a)\). The component note gives the primary source,
independent enumeration, and correction scope; the separate
arithmetic-progression sufficiency theorem is not implicated by this
specific failure.

## 4. The quadratic maps connect exactly at the level of words

The [quadratic component and rank atlas](quadratic_escape_rank_atlas_20261004.md)
recovers the three specified integer graphs. Their finite preperiodic
sets coexist with infinitely many escaping components. For \(x^2-c\),
\(c=0,1,2\), each escaping component has a unique least positive label
b with b+c nonsquare. Repeated positive integer square-root stripping
finds it. This is a terminating component classifier, and its terminal
labels need not represent a chosen periodic basin. That distinction
also matters when reducing Collatz to rank-local minima.

There is an exact further connection. For a valuation word w, retain
its matrix, slope and carry:
\[
 M_w=\begin{pmatrix}3^r&B_w\\0&2^A\end{pmatrix},
 \qquad\lambda=3^r/2^A,\qquad J=\lambda+\lambda^{-1}.
\]
Word doubling is matrix squaring, and therefore
\[
 \boxed{\lambda\mapsto\lambda^2,\qquad J\mapsto J^2-2.}
\]
This realizes the first and third quadratic maps on repeated word
presentations. It gives no analogous operation for \(x^2-1\), and no
conjugacy of Collatz starting integers to a quadratic orbit.

The affine anchor stays fixed while the raw gap and carry grow together.
For word12, \((3^r,B,2^A)=(9,5,8)\); doubling gives \((81,85,64)\).
The new factor17 cancels, leaving the same anchor -5. New factors in a
repeated presentation need not mean a new fixed point. Conversely the
same slope counts do not determine the anchor: word21 has anchor -7.
The actual source guard must still be checked;27 admits12 but not1212.

## 5. Why the recursive object needs arithmetic memory

The [guarded-pumping theorem](collatz_guarded_pumping_memory_20261004.md)
makes the cost of repetition exact. For an inverse shortcut word with
length L, \(r\ge1\) odd inverse steps and carry B, put
\[
 q=2^L,\quad d=3^r,\quad\Delta(x)=(q-d)x-B.
\]
If \(\Delta(x)\ne0\), its exact integer-repeat budget is
\[
 k\le\left\lfloor v_3(\Delta(x))/r\right\rfloor.
\]
The forward mirror consumes binary precision. Switching between rational
anchors can replenish precision only at an explicit equality boundary;
the next cofactor must be retained there.

Consequently every regular language made entirely of valid first-hit
root words has bounded odd-step depth, and covers at most a polynomial
in \(\log X\) inputs below X. The full first-hit language is
unconditionally nonregular, witnessed by explicit rooted families.
This does not rule out a finite program with unbounded integer registers.
It explains why a loop in a finite proof diagram needs a checked guard,
rather than an assumption that a successful word can repeat freely.

Independently, no finite collection of eventually nondecreasing height
charts can strictly rank every chosen forward macro of uniformly bounded
length. Arbitrarily long increasing one-runs force a repeated chart.
The valuation ranks and common-future moves under study escape those
hypotheses; they retain precisely the additional information the
obstruction excludes.

## 6. Current proof obligation and evidence

A [negative-cycle shadow theorem](collatz_negative_cycle_shadow_20261004.md)
gives the next exact boundary. Fix bounds R and S on the source and child
word lengths. If
\[
 n\equiv-5\pmod{2^{\lfloor3R/2\rfloor+1}},\qquad 3^S\mid n,
\]
there is no smaller-child common-future join within those bounds whose
affine child map \(h=\lambda n+b\) has \(b>-1\). Individual halving
exponents are unrestricted. Every row of the checkpoint grammar satisfies
this intercept condition.

The source's checked prefix imitates the signed cycle
\(-5\to-7\to-5\), while the ternary divisibility restricts the child
intercept. Together these force any putative paying join to have
\(\lambda\ge1\), contradicting its source payment. Unlike the earlier
\(-1\) binary shadow, the \(-5\) shadow is attained by positive powers
of 3 to every finite binary precision. Thus every such bounded-depth bank
misses an infinite exponent progression inside \(a\equiv3\pmod8\).
For example, source depth at most 6 misses \(3^a\) with
\(a\equiv11\pmod{256}\) once \(a\ge S\).

This explains why the new exponent family does not close the whole
problem. It points toward growing depth, a child construction outside the
stated intercept condition, or another well-founded proof move. It does
not exclude those alternatives or an adaptive controller with unbounded
arithmetic state.

For this particular shadow the retained memory is especially concrete.
An actual valuation block \((1,2)\) satisfies
\[
 U^2(n)+5=\frac98(n+5),
\]
and the number of consecutive legal blocks at a positive odd source is
exactly \(\lfloor(v_2(n+5)-1)/3\rfloor\). Each block consumes three
binary digits. The extra one in the guard enforces an odd endpoint:
mere integer divisibility would incorrectly accept the block at n=3,
where its formal endpoint is 4. A productive next controller must account
for the cost of changing patterns at this finite fuel boundary.

The live board is **original rank / checkpoint carry / ternary guard /
repeat fuel / modular action / grounded component**. Anchor: cover every
positive source by a paying dependency. Niche: legal compressed group
lifts. Wildcard: quadratic dynamics of repeated word presentations.
The research moves are to move the operation before changing its bound,
retain the quotient's missing coordinate, and classify the exact
exceptional domain before interpreting density.

The missing theorem is still global: every unresolved source must admit
a verified lower-ranked dependency, or its continued extension process
must have a proved well-founded rank. Current constructions discharge
new infinite domains and expose narrower residuals. Neither finite
solvability, modular reachability, repeated-word identities nor a larger
finite census supplies that last implication.

All six component packages have hand proofs, declared hostile controls,
independent audits, and matching normal/optimized Python outputs. The
experiments use exact integer or rational arithmetic. The incoming
critical-rank and negative-cycle gap audit was integrated from the live
shared branch; its newly recorded repeated-word cancellation is used
only in the scope stated above. The later incoming paid-controller package
was also integrated, with its complementary exponent domains combined by
the explicit modulo-32 disjointness proof.
