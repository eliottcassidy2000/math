# Exact basin and clock observers for the rank-one pair chain

**PROVED:** the elementary evaluation, terminal-label, finite-observer, and Haar-null ROOT-basin statements below. **FINITE-EXACT:** the supplied cycles and receipts in the declared control universe. **NOT AUDITED HERE:** the probability proof of THM-4610. Its file labels the theorem proved, while [HYP-9244](../hypotheses/HYP-9244-polya-trichotomy-for-haar-coalescence-the-debt-lattice-rank.md) explicitly records that its audit was pending. No conclusion here depends on promoting that status. HYP-9244's general classification and [HYP-9245](../hypotheses/HYP-9245-rank-one-integer-laws-for-the-base-p-collatz-maps-basin-boundaries-and-repunit-orphans.md) remain open in their stated scopes.

The missing information is not one more drift statistic. An affine pair state, its evaluation at a particular source, the eventual terminal cycle, and the clock phase are different objects. The exact examples below separate them. This is a small integration and certificate diagnostic, not a new Haar-coalescence theorem or new ordinary ROOT coverage.

## 1. Inheritance and the map

The closest mechanism is [THM-4610 — rank-one Haar coalescence for the base-p Collatz maps](../../01-canon/theorems/THM-4610-rank-one-haar-coalescence-for-the-base-p-collatz-maps-on-z-p.md), with the source discussion [zp_rank_one_coalescence_20261008](zp_rank_one_coalescence_20261008.md). For a prime $p$, put $P=p+1$ and

\[
C_p(x)=
\begin{cases}
x/p,&x\equiv0\pmod p,\\
(Px+p-i)/p,&x\equiv i\pmod p,\quad1\le i<p.
\end{cases}
\]

For $p=2$, this is the shortcut Terras map. On positive integers it preserves positivity. The distinguished cycle is $1\to2\to\cdots\to p\to1$.

The canonical hostile is a genuine positive cycle avoiding 1. The corrected near miss is the former sheet-blind reading of multiplicative debt: [MISTAKE-589 and MISTAKE-590](../../01-canon/MISTAKES.md) retain additive branch data and accessibility obstructions. The least-used sidecars here are exact reference evaluation and the terminal cycle's phase. A branch/norm marker transports a finite guard; it does not authenticate a terminal.

The live objects are: formal affine frame; ordinary source pair; finite local observation; witnessed terminal cycle; and the odd-step versus shortcut clock. Their interfaces are made explicit below.

## 2. Formal absorption and numeric equality are different events

Write the inherited affine relation as

\[
u=P^k v+\frac{A}{P^{\max(0,-k)}},\qquad k,A\in\mathbb Z.
\]

The exact evaluation numerator is

\[
\boxed{J(k,A;v)=
(P^{\max(k,0)}-P^{\max(-k,0)})v+A.}                 \tag{1}
\]

When the relation evaluates to an ordinary integral $u$, it satisfies

\[
J=P^{\max(0,-k)}(u-v).
\]

Thus $J=0$ is precisely numeric equality. Formal identity is the stronger condition $k=A=0$. If $k\ne0$, equation (1) has at most one rational solution $v$; that singleton is Haar-null, but can be exactly the ordinary input being studied. If $k=0,A\ne0$, it has none. Consequently the distinction disappears almost everywhere for a fixed frame, and cannot be discarded pointwise.

The code derives a step independently from the rational branch affine maps. With branch slopes $m_i\in\{1,P\}$, constants $b_i=0$ or $p-i$, and source/reference digits $j,i$,

\[
M'=m_jM/m_i,\qquad e'=(m_je+b_j-M'b_i)/p.
\]

It retains the integer frame lattice and checks against direct ordinary edges. No random-digit theorem is used.

Three exact one-step controls are:

| base | source pair | common next state | frame after the step | witnessed basin |
|---|---|---:|---|---|
| 2 | $(10,3)$ | 5 | $k=-1,A=10$ | distinguished cycle |
| 3 | $(30,7)$ | 10 | $k=-1,A=30$ | cycle through 7, avoiding 1 |
| 11 | $(7711,642)$ | 701 | $k=-1,A=7711$ | cycle through 642, avoiding 1 |

The initial frame in each case is $k=0,A=u-v$. After equality the two actual branch digits agree forever, so $k=-1$ stays fixed. Since the common state $v$ is positive, $A=pv\ne0$. These pairs really have merged, although their inherited formal frames never enter $(0,0)$. This is a pointwise specialization, not a counterexample to a Haar-almost-everywhere equivalence.

The supplied non-ROOT cycle for $p=3$ is

\[
7,10,14,19,26,35,47,63,21.
\]

The script supplies and authenticates all 57 edges of the $p=11$ cycle whose minimum is 642. It does not discover cycles by an unbounded search or claim these are the only positive cycles. The nontrivial $p=3,p=11$ cycles were already reported in THM-4610's computation; our new use is the exact pair-state and terminal observer test.

## 3. A complete terminal label on the certified periodic domain

Suppose a finite checked path enters a supplied primitive cycle

\[
\mathcal C=(c_0,\ldots,c_{\ell-1}),\qquad
C_p(c_i)=c_{i+1\bmod\ell}.
\]

Canonicalize the order by starting at the least vertex. If a path of length $s$ ends at $c_j$, define its label

\[
\boxed{(\mathcal C,\ j-s\pmod\ell).}                \tag{2}
\]

This is unchanged by extending the checked path along the cycle. Two sources with such receipts merge at equal shortcut times **if and only if** their labels agree. If the labels agree, after both paths have entered the cycle they occupy the same vertex at every sufficiently large equal time. Conversely an equal-time meeting forces all subsequent states equal, hence the same cycle and phase. This theorem applies to the supplied eventually periodic domain; it supplies no terminal finder for an unresolved source.

Reaching distinguished 1 asks only whether $1\in\mathcal C$. It forgets the phase in (2). Therefore:

* equal-time coalescence does not imply reaching 1 without a grounded terminal;
* two grounded sources need not coalesce at equal times;
* a checked merge plus a checked ROOT suffix does ground the original source.

On the distinguished cycle, sources 1 and 2 have different phases for every prime base, despite both visiting 1. More concretely, the inherited [grounded child-port factory](collatz_grounded_child_ports_20261007f.md) supplies

\[
151\xrightarrow{(1,1,10)}1,\qquad
75\xrightarrow{(1,2,8)}1.
\]

Both accelerated odd ROOT ranks are 3. Their shortcut ROOT times are respectively 12 and 11. Since the shortcut ROOT cycle has length 2, their phases differ; the shortcut trajectories never meet at equal times. The script expands each supplied valuation $a$ to one odd shortcut step followed by $a-1$ even steps, and checks strict first-hit ROOT. This is an actual map between clocks, not an identification of the clocks.

In general, two sources already certified to reach 1 at times $s,t$ under $C_p$ merge synchronously exactly when $s\equiv t\pmod p$. Allowing unequal times gives a different, coarser grand-orbit relation. In particular, a synchronous comparison must not silently replace the odd-step clock by the shortcut clock.

## 4. Every finite local observer can miss the exact equality event

Fix any prime $p$, positive integer $c$, and depth $K\ge1$. Use the fixed frame

\[
k=-1,\quad A=pc,\qquad u=(v+pc)/(p+1).
\]

At $v=c$, one has $u=v=c$. For every positive integer $s$, take

\[
v'=c+(p+1)p^K s,\qquad u'=c+p^K s.                 \tag{3}
\]

Both are positive integers, satisfy the same frame, and are unequal. The two reference values $v,v'$ agree modulo $p^K$, hence have the same first $K$ branch digits. Indeed, within one branch the difference is multiplied by a $p$-adic unit and divided by $p$, so agreement loses exactly one available digit per step. The same is true for the two source values $u,u'$.

Thus no classifier using only this frame and finitely many source/reference residue digits can decide numeric equality for all its ordinary specializations. This is a narrowly specified finite-observer obstruction. It is not an impossibility theorem for algorithms with access to the exact integers; equation (1) decides the event immediately.

The [marked normalized Gaussian norm](collatz_square_norm_transfer_20261007f.md) gives a stronger explicit control. For $K\ge2,d\ge1$, put

\[
u'=5+3^d2^K,\qquad v'=5+3^{d+1}2^K.
\]

They satisfy $u'=(v'+10)/3$. Both agree with 5 modulo $2^K3^d$, have the same markers modulo 4 and 3, and have

\[
q(u')\equiv q(v')\equiv q(5)\pmod{2^K3^d},
\qquad q(n)=(n^2+1)/2.
\]

For example, the divisibility follows from $q(x)-q(5)=(x-5)(x+5)/2$, with the extra even factor in $x+5$ restoring the division by 2. The imported marked-norm guard decoder independently checks the collisions. Retaining both finite norm observers still loses the exact equality event. Retaining the **full** positive norm does not: $n=\sqrt{2q-1}$ recovers the original source uniquely. This distinguishes finite precision loss from the false claim that the full positive norm itself is noninjective.

## 5. The Haar ROOT event has measure zero

Each residue branch of $C_p$ is a bijection onto $\mathbb Z_p$. Its inverses at $y\in\mathbb Z_p$ are

\[
x_0=py,\qquad x_i=\frac{py-p+i}{p+1},\quad1\le i<p.
\]

The $p$ points have distinct residues $i\pmod p$. It follows inductively that

\[
|C_p^{-N}(1)|=p^N.
\]

The set of $p$-adic points that ever reach distinguished 1 is the countable union of these finite sets, and is Haar-null. The same holds for the basin of any supplied finite cycle. This elementary fact needs no coalescence theorem. In particular, **Haar-almost every $p$-adic point never reaches 1**, even if every positive integer did. Positive integers themselves form a Haar-null subset.

There is no contradiction between a Haar coalescence result and this fact: it compares two moving trajectories, not either trajectory with a fixed authenticated terminal. Nor may ordinary natural-density basin experiments be read as Haar basin masses. The numerical basin percentages in the inherited note remain numerical statements about their finite ordinary universes.

The $p=3,p=11$ controls therefore block only the proposed implication from shared rank/drift/coalescence inputs **alone** to distinguished ROOT convergence, when those inputs are used identically across the base-$p$ family. They do not exclude a stronger argument adding a source-specific ROOT witness, a well-founded paid dependency, or information special to $p=2$. The general classification and tail laws in HYP-9244/9245 are not needed here.

## 6. Receipt and reproduction boundaries

The source object is an exact ordinary path or a typed affine frame with its marked reference. The target observer is either equation (1) or the terminal label (2). The preserved predicates are actual numeric equality and, on supplied periodic receipts, exact synchronized coalescence. A finite local quotient loses the exact reference; forgetting the terminal loses ROOT status; forgetting clock phase loses equal-time information. The repair stores those coordinates explicitly. No norm/source substitution, fabricated suffix, or orbit discovery is permitted by the APIs.

The cycle verifier rejects duplicate states and checks the closing edge. `root_path` rejects wrong valuations and ROOT padding. `terminal_label` checks an actual path into the supplied canonical cycle. Its path need not be a first-entry path; its label is independent of that choice. The `rational_C` inverse-layer controls explicitly operate in rational points of $\mathbb Z_p$, not only positive ordinary integers.

Reproduction:

```powershell
python -B -X utf8 04-computation/experiments/collatz_basin_observer_20261009.py
python -B -O -X utf8 04-computation/experiments/collatz_basin_observer_20261009.py
```

The declared universe is: $p\in\{2,3,11\}$, $k=-3,\ldots,3$, $A=-12,\ldots,12$, reference $v=1,\ldots,30$, filtered only by positive integral source evaluation; 120 post-merge steps for each of the three named pairs; collision depths $1,\ldots,12$ at $c\in\{1,5,7\}$; mixed norm depths $K=2,\ldots,12,d=1,\ldots,8$; inverse layers through depth 6, 5, 3 for bases 2, 3, 11; the supplied cycles and factory words; and malformed-receipt controls. Optimized execution uses explicit checks, not removable assertions. The saved output records the exact count. These finite controls support the displayed elementary proofs; they are not evidence of universal ordinary convergence.
