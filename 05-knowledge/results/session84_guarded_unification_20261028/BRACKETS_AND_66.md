# Two exact routes to66: odd-square brackets and the Berggren root

Status: **PROVED inherited bracket results and elementary synthesis; FINITE-EXACT checks.** The bracket classification and its multiplier window are inherited, not newly discovered here. A parameter plateau is not a dynamical attractor.

## 1. Recover the original meaning of {2,3,11}

For j>=1 define B_j=((2j-1)^2,(2j+1)^2]. A prime p has a unique such bracket. The original question asks whether2p stays in p's own bracket, not whether2p is prime.

The qualifying primes are exactly **{2,3,11}**. Indeed the simultaneous inequalities require

\[
2(2j-1)^2<(2j+1)^2.
\]

For j>=3 the reverse weak inequality already holds, since the difference is4j²-12j+1>0. Thus only B_1=(1,9] and B_2=(9,25] matter. The primes retained there are2,3 and11, respectively.

Sources: [original bracket theorem](../procgen_numerology_20260923_odd_square_brackets.md), [pairing/transition continuation](../procgen_brackets_20260924_pairings_transitions.md), and [later arithmetic-clock recovery](../arithmetic_clock_tournaments_20261021/session24-arithmetic/README.md). No prime-density heuristic is used in this theorem.

## 2. The exact stability interval

Let E_lambda be the primes whose scaled values lambda*p remain in the same bracket as p, with lambda>0. Then

\[
E_\lambda=\{2,3,11\}\quad\Longleftrightarrow\quad25/13<\lambda\le25/11.
\]

**Necessity.** Retaining11 gives9/11<lambda<=25/11. Excluding5 then forces lambda>9/5, because its alternative exclusion lambda<=1/5 is impossible. Excluding13 now forces lambda>25/13.

**Sufficiency.** On the displayed interval,2,3,11 remain. All other primes in B_1 or B_2 leave above the upper endpoint. In B_3, the least prime is29 and49/29<25/13. For j>=4, the ratio of bracket endpoints is at most81/49<25/13 and decreases with j. Hence every prime in those brackets also leaves above the upper endpoint.

The lower endpoint is strict because13*25/13=25 is still in B_2; the upper endpoint is included because11*25/11=25 remains. The first boundary adds13, while just above the second boundary11 is lost. These are actual guard changes.

For integer multipliers k>=2 the complete permitted pairs (k,p) are

\[
(2,2),(3,2),(4,2),(2,3),(3,3),(2,11).
\]

No larger bracket can work once doubling fails. The common doubling set therefore has product66 throughout the entire stated real-parameter cell, not just at the isolated value2. This explains robustness of the support selection under that particular rule. It says nothing about an iteration converging to66.

## 3. A second exact construction gives the same product

The primitive Pythagorean root is(3,4,5). Its even-prime semiprime image is obtained from

\[
(q-2,\ 2\sqrt{\sigma(2q)},\ q+4)=3(3,4,5).
\]

The first and last coordinates both give q=11; the middle gives sigma(22)=36. Thus22→14 is a square-sum aliquot edge. The prime3 has sigma(3)=4 and is coprime to22, so multiplication by3 lifts this to

\[
66\to78,\qquad\sigma(66)=144.
\]

The full guarded Berggren map and its completeness proof are developed in the companion Berggren note. Here the exact agreement is

\[
\prod_{p\in E_2}p=66=3\cdot(2\cdot11).
\]

This is a genuine equality of two explicitly defined constructions. It is not a proof that bracket dynamics and the Berggren tree are conjugate, nor a rule identifying their larger families. The branch choice, even-prime constraint and coprime square-sigma multiplier are all part of the construction.

## 4. Where F=U+S sits, and where it does not

Use the original **count** convention on D(N)={d|N:1<d<N}:

\[
F=\tau(N)-2,\quad U=\omega(N)-[N\text{ prime}],\quad
S=2^{\omega(N)}-1-[N\text{ squarefree}].
\]

Here S is not sigma, which we writeΣ when needed. The underlying prime set is contained in the squarefree set; F=U+S is not a disjoint partition.

For66, (F,U,S)=(6,3,6), so its defect F-U-S is-3. Doubling one of its prime exponents gives132,198 or726, each of shape p²qr and hence balanced, with(F,U,S)=(10,3,7). Their sigma values are336,468,1596, all nonsquares. The precise exponent-promotion obstruction is proved in the FUS note: a square-sigma source cannot keep that property when a prime of exponent1 is promoted to exponent2.

The existence of other balanced square-sigma sources, such as364=2²*7*13 with sigma784, prevents a false global incompatibility claim. Balance is a constraint on the divisor-count shape; square sigma is a constraint on evaluated prime weights. Their intersection and the operations transporting them must be studied separately.

## 5. Checks and scope

Run `python3 check_brackets_66.py` and the same command under `-O`. The standard-library checker directly tests all primes through100000 for the exact bracket set, checks rational multiplier boundaries, enumerates the integer-multiplier pairs, and verifies the66/promotions arithmetic by divisor enumeration. All-height conclusions use the finite bracket-reduction proof above, not a prime search cutoff.

This packet treats recurring numbers as prompts to look for a specified map. It records both the map and what it does not preserve. That is the useful common pattern with the source-retaining aliquot carrier and the guard-aware tournament/branch observations.
