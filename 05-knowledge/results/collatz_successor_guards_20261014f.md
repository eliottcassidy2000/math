# Predecessor divisibility, actual Collatz prefixes, and the missing new-prime guard

**Status:** PROVED elementary transfers and obstructions below; CITED Dirichlet prime supply; FINITE-EXACT controls in the companion JSON. Full Giuga/Korselt intersection existence and universal Collatz convergence remain OPEN. No ROOT conclusion is inferred from prime density, a modular return, or a paid prefix alone.

## 1. Inherited interfaces and the proposed transfer

The current arithmetic interface is the complete factor certificate, not an unlabelled congruence. The audited [grand-circuit typing note](collatz_mod6_20260921_grand_circuit_typing.md), §1, identifies a composite Giuga power-sum counterexample with a squarefree composite satisfying, for **every** prime factor p,

\[
 p\mid N/p-1,\qquad p-1\mid N-1.
\tag{1}
\]

The second condition is Korselt's condition. The [successor-block compiler](carmichael_successor_blocks_20261014b.md) retains the source N, its prime support, its conductor \(\lambda=\operatorname{lcm}_{p\mid N}(p-1)\), and the slack \(e=(N-1)/\lambda\). Its simultaneous two-prime repair of561 is an arithmetic completion, not a Collatz orbit operation. The [Collatz holonomy packet](collatz_successor_holonomy_20261014b.md) already proves that even full affine quotient identity and the initial residue can agree for growing and declining literal words.

The live objects here are complete versus marked prime support, the actual integer source, the 2-adic native word, and the independently supplied terminal word. The useful transfer retains the conductor and slack; the obstruction uses the least-used coordinate, the new prime's **own** two tests in(1).

## 2. A complete Korselt certificate gives a source-owned descent test

For every positive odd N>1,

\[
 U(N)<N\quad\Longleftrightarrow\quad N\equiv1\pmod4,
 \qquad U(N)=\frac{3N+1}{2^{v_2(3N+1)}}.
\tag{2}
\]

Indeed the other class has valuation1 and grows, while valuation at least2 gives \(U(N)\le(3N+1)/4<N\).

For a completely certified odd Carmichael source, any factor \(p\equiv1\pmod4\) forces \(N\equiv1\pmod4\). If such a factor is \(1\pmod8\), the first valuation is exactly2. The precise first-rise complement consists of an odd number of factors, all \(3\pmod4\). For example \(8911=7\cdot19\cdot67\) passes all Korselt tests but has first valuation1. This is a missing exit for this rule, not a nonconvergence certificate.

The conductor/slack coordinate gives the exact number

\[
 R=v_2(N-1)=v_2(\lambda)+v_2(e).
\tag{3}
\]

In fact the following identity needs no primality assumption. If R>=3, put \(r=\lfloor(R-1)/2\rfloor\). Then the actual word \(2^r\), meaning r letters of value2, has endpoint

\[
 y=1+\frac{3^r(N-1)}{4^r},\qquad 1<y<N.
\tag{4}
\]

At prefix j<r the current value minus1 has valuation R-2j>=3; therefore \(3n+1=4+3(n-1)\) has valuation exactly2. All these prefixes decline. At the endpoint, R odd gives \(y\equiv3\pmod4\) and next valuation1; R even gives \(y\equiv5\pmod8\) and next valuation at least3. R=2 is handled by its actual first valuation at least3. R=1 is the explicit open branch of this dispatcher. Formula(4) stops before ROOT and never inserts ROOT padding.

The code authenticates all source factors before using(3). Its `descent` returns a literal receipt or `None` for R=1. `discharge` requires a supplied child ROOT word and checks the exact integer seam. The fixed561 control uses

\[
 (2,4,1,1,1,3,4,1,3,1,2,3,4)
\]

with odd rank13 and cost30; its first receipt is561→421. Among the inherited sixteen Carmichael sources at most100000, fifteen receive a paid prefix and8911 does not. Notably15841→8911 by(2,2): landing at a smaller Carmichael source need not preserve this rule's ready-to-descend branch.

**Separate finite discharge of the missing named source.** A bounded discovery with an explicit100000-step cap found the first descent

\[
8911\xrightarrow{(1,1,1,3,1,1,1,5)}3569.
\]

The checker retains that word and verifies every earlier odd checkpoint is at least8911. The supplied strict ROOT suffix from3569 is

\[
(2,5,1,2,1,2,1,1,1,1,1,3,1,1,1,4,2,2,4,3,1,1,5,4).
\]

Production only replays these data. The complete8911 word has odd rank32 and cost64. Concatenating the authenticated15841→8911 receipt gives an explicit15841 ROOT word of rank34 and cost68. Thus all sixteen named sources have a paid exit if this singleton is added, while561,8911 and15841 have the separately recorded full ROOT receipts in this packet. The generic R=1 rule is still unresolved: one supplied singleton does not authenticate its other inputs.

**An infinite first-rise subfamily, with the same word.** The exact carrier and native cell of the eight-letter first-descent word are

\[
(P,Q,B)=(6561,16384,9425),\qquad n\equiv8911\pmod{32768}.
\]

For every integer t>=0 it therefore gives the actual relation

\[
8911+32768t\longmapsto3569+13122t,
\qquad n-U^8(n)=5342+19646t>0.
\tag{10}
\]

The seven proper prefix pairs(P,Q) are
\((3,2),(9,4),(27,8),(81,64),(243,128),(729,256),(2187,512)\).
Each has P>Q and positive carry, so its value is greater than the unchanged source for every positive member. Thus(10) is the **first** strict descent at exactly eight odd steps, without earlier ROOT. Every source is3mod4, so this supplies a genuine infinite guard inside the first-rise branch; it needs no Carmichael hypothesis. It does not prove that infinitely many sources in this ray are Carmichael, nor that every child3569+13122t is rooted. Only t=0 is discharged by the supplied suffix above. The exported `ray_descent` retains that distinction.

The ray has **no positive Mersenne member**. For E>=15, \(2^E-1\equiv32767\pmod{32768}\), not8911. For1<=E<15 the native congruence would require \(2^E=8912\), which is not a power of two. Thus this exact new guard does not discharge the earlier giant Mersenne source obligations.

## 3. Every finite native word survives two marked Giuga/Korselt tests

Let q be prime with

\[
 q\equiv257\pmod{330},\qquad N=33q.
\tag{5}
\]

Then N is squarefree and the two tests in(1) both pass at each marked prime3 and11. Explicitly they require q2mod3, q4mod11 and q7mod10, exactly(5). At the remaining prime q they would require respectively

\[
 q\mid32,\qquad q-1\mid32.
\tag{6}
\]

Since q>=257, both fail. No member of this family is a full Giuga or full Korselt completion.

Now prescribe **any nonempty finite word** w of positive valuations. Write its affine carrier as

\[
 F_w(n)=\frac{Pn+B}{Q},\qquad P=3^{|w|},\quad Q=2^{\sum w}.
\]

Its exact native odd-source congruence is

\[
 n\equiv r_w=(Q-B)P^{-1}\pmod{2Q}.
\tag{7}
\]

Final oddness forces all the intermediate valuations: reverse the affine steps, or induct on the last denominator. The only additional issue for a strict first-hit word is reaching1 before its end. Taking N>=Q excludes that: every nonempty prefix with cost Q_j<=Q has positive carry and value at least \(3^jN/Q_j\ge3^j>1\).

On N=33q, (7) asks \(q\equiv r_w33^{-1}\pmod{2Q}\). This and(5) agree modulo2, so CRT gives one class

\[
 q\equiv a_w\pmod{165\cdot2Q}.
\tag{8}
\]

It is reduced: the class is odd and congruent257 modulo165, a unit. The compiler computes it by writing q=257+330j and solving

\[
 j\equiv\frac{r_w33^{-1}-257}{2}\,165^{-1}\pmod Q.
\tag{9}
\]

Here the dyadic residue representative in the numerator is odd, so the displayed division is integral; the inverse of165 is taken modulo Q.

**CITED prime supply.** Dirichlet's theorem gives infinitely many primes in the reduced progression(8), and thus infinitely many beyond the strict-word cutoff N>=Q. See [Dirichlet's original1837 paper, translated by Ralf Stephan](https://arxiv.org/abs/0808.1408). Combining this input with the proved CRT calculation gives infinitely many positive squarefree N passing all marked tests at3 and11 and realizing the prescribed word w, while always failing both tests at q.

In particular words \(1^k\) give arbitrarily long initial rising runs within this marked family. This rules out a uniformly bounded first-descent horizon from those two marked guards. It does **not** construct one positive integer with an infinite rising word: the source changes when the word changes.

## 4. What analytic prime input can and cannot complete

The prime-supply step above is unconditional. Conditional least-prime bounds can make a search effective with a sharper size estimate; they do not remove(6). For a precise primary reference on the relevant analytic task see [Bach and Sorenson, Explicit Bounds for Primes in Residue Classes,1996](https://doi.org/10.1090/S0025-5718-96-00763-6), also [the author's institutional record](https://digitalcommons.butler.edu/facsch_papers/1121/). No particular numerical bound from that paper is used here.

This is the exact boundary of the proposed RH/GRH transfer: analytic hypotheses about prime supply may discharge a prime-in-progression obligation, but the current construction still must satisfy the new prime's algebraic guard, bind the resulting source, and supply a ROOT suffix. The obstruction does not exclude a different use of analytic information. It excludes promoting these retained marked guards to a full certificate.

The current [S32 collector result](prime_power_bias_collatz_collector_20261014.md), Proposition2, is similarly an ensemble law over initial residues. It is not a stopping-time or orbit-convergence theorem. The known ROOT cycle visits1 and2 with equal time frequencies even though the limiting input-ensemble proportions modulo3 are1/3 and2/3. We inherit only this correctly scoped distinction.

[THM-4644](../../01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md) counts invariant connection sets under a specified finite multiplier action. Its Wieferich orbit splitting changes the number of available invariant choices. It supplies no map to an actual Collatz source, word, or terminal proof. Neither that finite action nor an AP prime count replaces the source and new-factor coordinates in this packet.

## 5. Exact controls, API, and remaining target

The [script](../../04-computation/experiments/collatz_successor_guards_20261014f.py) exports `KorseltSource`, `audit_source`, `predecessor_state`, `descent`, `discharge`, `ray_descent`, `compile_partial`, `audit_plan`, and `realize_partial`. Production paths take supplied exact primes and words; prime search occurs only in the explicitly capped test driver. Type checks precede equality, so bool/float aliases, altered source factors, and malformed phases are rejected.

The finite universe is the sixteen inherited Carmichael sources; every odd n from3 through8191 for(4); all84 words on valuations1..4 of lengths1..3; and three additional rising words of lengths4,8,12. Formula(10) has literal controls t=0..127 and exact checks of all seven prefix slope inequalities. Each of the87 word controls supplies an independently trial-divided prime, verifies every marked guard and every missing new-prime guard, replays by repeated division, and reverses its word with exact division by3. For example:

| Word | Supplied q | Actual source33q | Endpoint |
|---|---:|---:|---:|
| (1) |587|19371|29057|
| (2) |257|8481|6361|
| (1,2,3) |16427|542091|228695|
| twelve ones |2292767|75661311|9816778151|

The three finite ROOT controls are supplied separately and stop at their first1. The JSON records each source and exact word, rank, cost, and word hash. Wrong child suffixes, wrong integer seams, and padding after ROOT are rejected. No finite control supplies a universal terminal theorem.

The next source-preserving obligation is concrete: cover the remaining Korselt complement R=1 outside(10), or ground a new child in(10) with an independently authenticated dependency. Enlarging a marked prime bank must recheck the newly added prime, not merely provide it in a reduced progression. A source-changing prime-existence theorem and a source-preserving ROOT receipt remain different interfaces.

Reproduction, from the repository root:

```text
python -B 04-computation/experiments/collatz_successor_guards_20261014f.py
python -B -O 04-computation/experiments/collatz_successor_guards_20261014f.py
```

The companion [JSON](collatz_successor_guards_20261014f.json) is the deterministic LF-normalized output. Normal and optimized runs agree on **4,851 checks**, LF SHA256 `2d04ba1afb459ac43fff935eba72d829c5eb73f76709e380d66da242b21fe38b`.
