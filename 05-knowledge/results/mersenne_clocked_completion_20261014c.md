# Mersenne completion: joint guards, fixed-head rigidity, and the clock a repeating rule owes

**PROVED scoped:** the fixed-head coincidence classifier, necessary clock bill and extension of inverse-guard failures to an exponent phase. **FINITE-EXACT:** the declared controls, complete inverse search through length 64, and companion supplied ROOT certificates. **CONDITIONAL:** every ungrounded family transports ROOT only from its actual child. The giant Mersenne obligation and universal Collatz remain **OPEN**.

Artifacts: [checker](../../04-computation/experiments/mersenne_clocked_completion_20261014c.py), [saved output](mersenne_clocked_completion_20261014c.json), [new joint guards](mersenne_joint_guards_20261014c.md), [child dependencies](mersenne_child_dependencies_20261014c.md), and [Wilson lift trees](wilson_defect_trees_20261014c.md).

## 1. Inheritance and the live concepts

The closest proved mechanism is [THM-4555, uniform switches are collisions at minus one](../../01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md), followed by the source-specific two-anchor rules in [THM-4601](../../01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md). The decisive hostile is the exact common future of 7 and 3: it has nonzero clock defect and cannot be promoted from one source to a fixed-head all-height family. The corrected near miss is the [MISTAKE-588 repair retained by the fibre ledger](collatz_fibre_integration_20261009.md): normalized absorption preserves both residual clocks, with explicit pre-ROOT scope. The least-used sidecar here is the *clock cost of returning to the same exponent phase*.

Our concept board is: unchanged source / joint native guards / child type / two clocks / grounded terminal / free lift coordinate. The anchor is the residual first-reset-two branch; the niche is a complete algebraic classifier rather than a longer undirected search; the wildcard is the Wilson product lift, where the number of genuinely free coordinates can be counted exactly.

The [previous successor packet](carmichael_successor_blocks_20261014b.md) gives a useful operational analogy. The source 561 has no one-prime Carmichael extension, yet the pair 401 and 3041 succeeds together. What transfers is checking the *whole coupled guard* while keeping the source fixed. It is not a map from those primes to Collatz valuations, nor a guarantee that an arbitrary failed step has a successful paired repair.

## 2. A complete fixed-head classifier

Let $U$ be the odd accelerated map, let $M_E=2^E-1$, and let

\[
X_E=U^{E-1}(M_E)=2\cdot3^{E-1}-1.
\]

Fix an exponent drop $D\ge1$ and two finite positive valuation heads $u,v$. We use reduced heads (first letter at least 2) or the empty head, so their attachment after the maximal initial ones run is unambiguous. Write

\[
F_u(x)=\frac{p_u x+b_u}{q_u},\qquad
p_u=3^{|u|},\quad q_u=2^{A_u},\quad A_u=\sum u.
\]

The empty head has $(p,q,b)=(1,1,0)$; appending $a$ changes the carrier to $(3p,2^a q,3b+q)$. Define the two exact integers

\[
C=q_vp_u3^D-q_up_v,\qquad
R=q_u(b_v-p_v)-q_v(b_u-p_u).
\]

For every $E\ge D+2$, rational endpoint equality is equivalent to

\[
\boxed{2\cdot3^{E-D-1}C=R.}\tag{1}
\]

This follows by cross-multiplying the affine formulas. It gives a complete three-way decision:

1. If $C=R=0$, the equality is an identity in $3^E$.
2. If $C=0\ne R$, it never holds.
3. If $C\ne0$, it holds at at most one exponent: $R/(2C)$ must be exactly $3^k$ for an integer $k\ge1$, and then $E=D+1+k$.

This is a source-line specialization of the inherited collision mechanism, not a new independent route around THM-4555. Its extra executable output is the complete single-exponent candidate, including cases with no integral or positive candidate.

### What an identity preserves

In case 1, unique factorization applied to $C=0$ forces

\[
A_u=A_v,\qquad |v|=|u|+D.
\]

For the full words $1^{E-1}u$ and $1^{E-D-1}v$, the normalized two-clock defect is

\[
\delta=(A_u-A_v,\ D+|u|-|v|)=(0,0).\tag{2}
\]

The coordinates are exactly those of the [source-authenticated fibre ledger](collatz_fibre_integration_20261009.md): residual Terras time and total number of odd steps. Thus a finite bank with fixed drops and fixed post-ascent heads has only finitely many possible *nonzero-defect* Mersenne coincidences. Its repeatable identities stay inside the old zero-defect fibres. A universal grounding construction using this representation must let some heads/drops vary without a finite bound, leave the Mersenne child type, or provide independent terminal proofs. This statement does not exclude any of those alternatives.

### Algebra is followed by the native guard

Equation (1) is not a valuation certificate. Every step must have the stated **exact** power of two, stay positive and odd, and precede the first occurrence of 1; the join itself must exceed 1. The checker exposes classification and authentication as separate operations.

The identity $u=(3),v=(2,1),D=1$ illustrates the distinction sharply. It is a true affine identity but has no native Mersenne source: the first reset of $M_E$ is

\[
1+v_2(3^E-1)=
\begin{cases}2,&E\text{ odd},\\3+v_2(E)\ge4,&E\text{ even}.
\end{cases}
\]

For a positive isolated control take $u=(2,3),v=(),D=1$. Here $C=-5$, $R=-30$, and (1) permits only $E=3$. Its authenticated join is the inherited $7\to11\to17\to13\to5\leftarrow3$, with defect $(5,3)$. For an algebraic-but-nonnative control, $u=(2,5),v=(5,1),D=1$ also singles out $E=3$ but uses wrong valuations. No new source is substituted to make that guard pass.

## 3. A repeating phase cycle must carry growing clocks

Suppose an arithmetic phase has exponents $E_j=E_0+jP$, where $P\ge1$ and $E_0\ge2$. Assume:

- $M_{E_0}$ has an independently supplied strict first-hit ROOT certificate;
- for every $j\ge1$, an authenticated pre-ROOT common-future receipt joins $M_{E_j}$ to $M_{E_{j-1}}$;
- every receipt has the same normalized defect $(a,b)$.

The receipts may have variable words; no fixed-head hypothesis is imposed here. By actual prefix alignment and suffix transport, every member is rooted. Write its residual Terras time as $s_j=\sigma_T(X_{E_j})=\sigma_T(M_{E_j})-(E_j-1)$ and odd-step count as $o_j$. The alternative normalization $\sigma_T(M_{E_j})-E_j$ is one smaller and has the same defect. The inherited defect-addition rule gives

\[
s_j=s_0+ja,\qquad o_j=o_0+jb.
\]

The initial Mersenne ascent alone gives $o_j\ge E_j$. Every Terras step decreases a state by at most a factor of 2, so

\[
s_j\ge\log_2 X_{E_j}>(E_j-1)\log_2 3.
\]

Divide both inequalities by $j$ and let $j$ grow. Consequently

\[
\boxed{b\ge P,\qquad 2^a>3^P.}\tag{3}
\]

For the second condition, the limit first gives $a\ge P\log_2 3$; equality is impossible for positive integral $P$ and integral $a$. Both coordinates matter. At $P=1$, $(a,b)=(1,1)$ fails the time condition and $(2,0)$ fails the odd condition. The smallest possible integer pair is $(2,1)$. Passing this necessary bill supplies neither a receipt nor its source guards.

If $b<P$, a finite contradictory lift can be printed without searching an orbit:

\[
j=\left\lfloor\frac{o_0-E_0}{P-b}\right\rfloor+1,
\quad o_0+jb<E_0+jP.\tag{4}
\]

For example, the rooted input $M_3=7$ has odd rank 5. A claimed period-one zero-odd-cost completion from that base fails by $j=3$: it would assign rank 5 to $M_6$, which already owes at least six odd steps. This rejects that proposed completion, not Collatz convergence for $M_6$.

This phase-cycle formulation repackages the prior fibre invariant as a design constraint for finite-state completion. It is not an obstruction to one paid deletion, or to grounding arbitrarily large sources by words whose defects grow. It also requires the return to the *same phase*; applying a rule once and landing outside its guard does not establish the hypotheses.

## 4. What the new blocks actually add

The [joint-guard packet](mersenne_joint_guards_20261014c.md) proves three further short paired templates on the $E\equiv3\pmod4$ branch:

| Native exponent guard | Paid Mersenne child | Child's initial two-run |
|---|---|---:|
| $E\equiv323\pmod{512}$ | $M_{E-2}$ | 4 |
| $E\equiv199\pmod{256}$ | $M_{E-2}$ | 2 |
| $E\equiv79\pmod{128}$ | $M_{E-2}$ | 2 |

These guards are disjoint from the inherited one-letter $E\equiv31\pmod{64}$ rule. They add relative mass $7/128$ *within that first shell*; this is not additional mass against the repository's entire adaptive rule bank. The least instances have supplied child ROOT certificates at exponents 321, 197 and 77, from which the corresponding sources and intermediate children are strictly grounded. The arbitrary family members retain a child ROOT obligation.

The first row also composes with an inherited long-two rule on the simultaneous guard $E\equiv52547\pmod{131072}$, yielding $M_E\leadsto M_{E-2}\leadsto M_{E-6}$. This is a useful paired-block repair with exact overlap information. It is still zero-defect transport, and its general child remains ungrounded. The least grounded child $M_{321}$ is outside the tested 225-head long-two bank; its separate ROOT word cannot be replaced by an assertion of that bank's coverage.

The [child-dependency packet](mersenne_child_dependencies_20261014c.md) supplies a second precise boundary: a fixed contracting inverse word has at most one possible Mersenne-to-smaller-Mersenne exponent pair. None of the complete 953 short contracting valuation words or 9,840 nonempty words in the declared three-generator bank has such a pair. Their arithmetic descent on ordinary children remains valid. The obstruction is closure in the Mersenne child type, not a failure of the inverse maps.

The packet also transports existing grounded evidence through the exceptional exponent 18269 and the even reset, grounding $M_{18267}$ by a checked suffix splice. That is a discharged terminal, not merely another labelled deletion edge.

**The actual giant obligation was tested as well.** Among the 1,949 inherited fan labels, the three new phases match 4, 7 and 15 sources respectively. All 26 children $E-2$ already lie in that same fan. Its minimum $E=99708991741$, the parameter-23 source $99708993705$, and the top $99708993713$ miss all three guards. These are exact small-integer label/phase calculations using the corrected inherited fan, not a replay of its enormous orbits. Thus these short rules add useful source families and compression, but they do not ground or escape the outstanding giant component.

### A source-specific obstruction to short ordinary-child exits

The attempted alternative was a smaller **ordinary odd** predecessor, without requiring a Mersenne child. Fix the actual minimum $E_*=99708991741$. For each $1\le L\le64$, the checker decides all positive valuation words of length $L$ and inverse slope $2^A/3^L<1$. A native word with carrier $(P,Q,B)=(3^L,2^A,B)$ must pass

\[
2^{E_*+A}\equiv Q+B\pmod{3^L}.\tag{5}
\]

**FINITE-EXACT result:** none passes, for any of these 64 lengths. The search uses only residues, never constructs $2^{E_*}$, and is exhaustive by the following recursion. Set $c_L=\lfloor\log_2(3^L)\rfloor$, exactly computed as an integer bit length. Start with $M_{E_*}\bmod3^L$. At a state with residue $n\bmod3^r$ and used cost $c$, try every $1\le a\le c_L-c-(r-1)$. If $2^a n-1$ is divisible by 3, pass to its quotient modulo $3^{r-1}$; otherwise reject that branch. Reverse the resulting inverse-step list to obtain the chronological forward word. This is a complete enumeration of native contracting words, with impossible branches pruned. A node-budget interruption raises an error and is not reported as an empty complete search.

The completed run visits 28,087,511 prefix nodes. Its formal universe contains

\[
\sum_{L=1}^{64}\binom{c_L}{L}
=7\,818\,118\,225\,914\,629\,120\,771\,112\,104
\]

words; that count describes the universe covered by pruning, **not** individually evaluated words. Independent direct carrier enumeration through length 12 checks 68,722 words and agrees. Small-source controls include the positive native inverse $31\rightsquigarrow27$ with forward word $(1,2)$, so the algorithm is not an always-empty test.

**PROVED extension to arbitrary valuation cost.** A positive word of length $L$ satisfies

\[
\frac BQ\le(3/2)^L-1,
\]

by the recurrence $B'/Q'=(3B/Q+1)/2^a$ and $a\ge1$. If $Q>P$, then $Q/(Q-P)\le P+1$, whence

\[
\frac{B}{Q-P}\le\big((3/2)^L-1\big)(3^L+1)
<(9/2)^L<2^{3L}.
\]

Therefore, whenever $E\ge3L+1$, the inverse endpoint $(QM_E-B)/P$ is larger than $M_E$. The equality $Q=P$ is impossible for $L\ge1$. Thus at this source, allowing arbitrary large valuations cannot rescue a smaller predecessor omitted by the contracting-word search.

Finally, every guard (5) for $L\le64$ depends only on $E\bmod2\cdot3^{63}$, by Euler's theorem for powers of 3. The least positive member of the phase below is already much larger than the size cutoff. Combining the exact base calculation with these two elementary lemmas gives

\[
\boxed{E\equiv99708991741\pmod{2\cdot3^{63}}
\quad\Longrightarrow\quad
\text{no }1\le h<M_E\text{ odd has }U^j(h)=M_E
\text{ for }1\le j\le64.}\tag{6}
\]

Here $2\cdot3^{63}=2289122546861674989771899392854$. This is a guarded **obstruction to short inverse exits**, not a ROOT obstruction. Longer inverse words, two-sided common-future receipts and forward descent remain available. Inverse descent alone is not universal even among already rooted inputs: none of $1,3,5$ reaches 7. The decisive positive target is still a grounded forward or common-future exit for the unchanged giant source.

## 5. Wilson triples and the right abstract transfer

The [Wilson companion](wilson_defect_trees_20261014c.md) proves the user's factorial observation for every $n\ge2$: the remainder is 2 exactly at 3 and 4; other primes give $-1$ and other composites give 0. The extra congruence for a Wilson prime is modulo $p^2$, with a fixed ordinary factorial.

There is an actual reversible ternary tree nearby: at $p=3$, labelled residue sections with product $-1\pmod{3^k}$ have three admissible next-digit lifts. In general there are $p^{p-2}$ children, with an explicit inverse decoder. The *unchanged factorial section* is a marked path within this much larger tree. It is not interchangeable with a newly corrected section.

This supplies a concrete common lens across the three lanes. A coupled extension has an aggregate defect, free coordinates, and a pivot guard. The Carmichael block chooses prime factors under simultaneous divisibilities; the Wilson tree chooses lift digits under one linear equation; the Collatz receipt must satisfy the native word guards and reach the actual child. Only the Wilson construction here proves a complete free-address bijection. None of these maps preserves the other two problems' arithmetic predicates automatically.

There is also an explicit bijection of **formal encodings** with the old Collatz ternary alphabet. Name its inverse generators

\[
G_1(n)=(2n-1)/3,\quad G_5(n)=(8n-5)/9,
\quad G_{17}(n)=(2048n-2363)/2187.
\]

Send the letters $(G_1,G_5,G_{17})$ to the free Wilson digits $(0,1,2)$ at $p=3$, apply the section-lift encoder, and invert using its free-digit decoder. A formal word of length $k$ becomes a section at depth $k+1$; deleting its last letter is exactly reduction of that section modulo $3^k$. The first three images are $(1,8),(4,2),(7,5)$. This is a prefix-preserving isomorphism between the **formal three-letter word tree** and the **labelled Wilson-section tree**, with a literal inverse. It preserves the chosen address, not Collatz integrality, numeric size, or ROOT status.

The actual giant source is a decisive hostile to dropping that qualification. Its exponent is 1 modulo 6 and 703 modulo 1458, so all three inverse generators fail their source guards; the corresponding three Wilson lifts still exist. Forgetting the source turns a pruned arithmetic tree into a full formal tree and fabricates coverage. The checker verifies all 1,093 addresses through word length six, retaining prefix reduction and inverse recovery.

The next decisive Collatz target is therefore **a source-native, variable-depth exit carrying nonzero clock defect into grounded evidence**, or an independently grounded ordinary child outside the Mersenne subcategory. A repeating exponent-phase construction must pay (3); a fixed-head proposal must pass the complete classifier (1) and the native guards. These are precise filters for allocating search effort, rather than a claimed proof of universal coverage.

## 6. Reproduction and audit scope

Run `python -B -X utf8 04-computation/experiments/mersenne_clocked_completion_20261014c.py`, then repeat with `-O`. The declared classifier control universe contains 36 empty/reduced heads of length at most two, first valuation 2 through 6 and second valuation 1 through 6, with drops 1 through 3: 3,888 comparisons. Independent successive rational evaluations check every allowed exponent through 24. The all-exponent classifier finds four identities, 13 singleton candidates and 3,871 impossible pairs. Exactly two singleton candidates are native; 11 fail the separate guards. The even-reset control uses all even exponents 4 through 98. Fourteen malformed/source/ROOT-padding API cases and the absent reset-three phase are rejected.

The later source-specific controls include the 1,949 inherited fan labels, all contracting inverse lengths through 64 for the unchanged giant source, an independent full-carrier path through length 12, all source exponents 2 through 32 at inverse lengths through six, explicit budget exhaustion, and all formal ternary addresses through length six. The 28,087,511 pruned-search nodes are reported separately from the explicit validation-check counter.

Normal and optimized output agree at **187,695 exact checks**, plus those explicitly counted search nodes. LF-normalized JSON SHA256: `1cdc908345d49a1b62b4dd653ba4ad63c7a1d0cd7b843da7b9f1098db18b4565`. None of the finite counts is used to promote an infinite ROOT statement. The clock/classifier core and final inverse/tree additions independently passed proof/API review; a peer reproduced the complete optimized output, and a separate modular kernel independently reproduced the inverse-search verdict. The three companion packets also passed independent proof/API and optimized-output audits.
