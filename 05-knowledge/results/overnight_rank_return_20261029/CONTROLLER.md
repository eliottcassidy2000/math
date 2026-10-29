# Integer termination and five failure cylinders for the cheaper ternary rule

**2026-10-29. PROVED C1 and C3 by elementary inequalities; PROVED C2 using
the explicit finite arithmetic certificate in controller.json; FINITE-EXACT
finite-height return cover. Same-agent audit, not independent review.
Universal Collatz descent remains OPEN.**

The inherited [source-refinement budget](../affine_tree_recurrences_20261028/REFINEMENT.md)
controls tails of the `(d,r)=(3,3)` grammar under odd Haar measure, but cannot
exclude a null exceptional integer. An integer-valued decreasing potential
now excludes every nonterminating positive-integer **single block**. It also
reduces all non-descending blocks to five explicit source cylinders.

## C1 Every positive integer block terminates

For positive odd n write a=v2(3n+1), S(n)=(3n+1)/2^a. Start a queue q=1,
read the actual chronological valuations, and update

\[
q'=q-1+3\mathbf1_{3\mid a}.
\]

Stop at the first q=0. The same result holds with any initial positive integer
q. Put f(n)=n^3+5 1_{n=3} and W(n,q)=4^q f(n). Then every executed step obeys

\[
\boxed{W(S(n),q')\le {125\over128}W(n,q).} \tag{1}
\]

**Proof.** No positive odd integer maps to 3: 3n+1=3*2^a is impossible modulo
3. For n>=7 the ratio is

\[
4^{-1+3\mathbf1_{3\mid a}}2^{-3a}(3+1/n)^3
\le (3+1/n)^3/32\le1331/1372<125/128.
\]

The coefficient bound uses a>=1 off the branch event and a>=3 on it. For
n=1,3,5 the ratios are respectively 1/4,125/128,1/500. This proves (1),
with equality only at n=3. The queue cannot cross below zero before hitting
zero. Thus W is a strictly decreasing positive integer while the parser
runs, and the parser must stop. For its length J,

\[
(128/125)^J\le4^{q_0}f(n_0).
\]

For q0=1 this is J=O(log n0). Every encountered integer satisfies
n_j^3<=4 f(n0); the bit lengths remain O(log n0). These are single-block
complexity bounds, not a bound on the number of blocks to reach 1.

**Hostiles and sharp scope.** At n=3 the block is 3->5, so termination does
not imply descent. On positive rationals with odd numerator and denominator,
n=1/5 is fixed with valuation 3; its queue increases by two forever. The
integer lower bound in the proof is essential, and a Haar-null set argument
alone could not supply it. This is a grammar-dependent potential, not a
decreasing function of the Collatz value alone.

## C2 Exactly five cylinders fail to descend

Let R(n) be the first-block endpoint. For every odd n>1,

\[
\boxed{R(n)\ge n\iff n\equiv3\pmod4\ \text{or}\ n\equiv125\pmod{128}
\ \text{or}\ n\bmod4096\in\{1357,2269,2877\}.} \tag{2}
\]

All failures in (2) strictly increase. The extra fixed source is R(1)=1.

| Exact word | Source residue | Endpoint |
|---|---|---|
| (1) | 3 mod 4 | (3n+1)/2 |
| (3,1,1,1) | 125 mod 128 | (81n+179)/64 |
| (3,3,1,1,1,1,1) | 1357 mod 4096 | (2187n+16177)/2048 |
| (3,1,3,1,1,1,1) | 2269 mod 4096 | (2187n+12289)/2048 |
| (3,1,1,3,1,1,1) | 2877 mod 4096 | (2187n+9697)/2048 |

**Proof.** Let m count branch letters. Queue balance gives J=3m+1.
With m=0 the word has one nonbranch letter; a=1 increases, while a>=2
decreases at every n>1. With m=1 the only tree shape has minimal word
(3,1,1,1). With m=2 there are exactly the three minimal words in the table.

Replacing any valuation by a larger valuation with the same branch status
decreases the resulting affine map on positive inputs. For m=1, a nonminimal
word has slope at most 81/128 and intercept at most 179/64. Hence it
decreases for n>358/47. The least possible source with m=1 is at least 13:
below 13 none has its initial valuation divisible by 3. For m=2, a
nonminimal word has slope at most 2187/4096 and intercept at most
16177/2048. It decreases for n>32354/1909<19. Direct checks of the odd
sources below 19 show none has two branch letters. Thus the only possible
failures for m<=2 are precisely the five minimal words, each of slope >1.

It remains to exclude every m>=3, rather than infer that from a sweep.
When all states of a block are >=33, the proof of C1 improves its factor to

\[
\rho={(3+1/33)^3\over32}={31250\over35937},\qquad4\rho^{10}<1.
\]

Here J>=10 and the terminal queue is zero, so
R(n)^3<=4 rho^J n^3<n^3. If a state is <=31, its entire forward orbit
lies in the following explicit closed set U:

```text
1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 35 41 47 53 61 71
91 103 107 121 137 155 161 167 175 233 251 263 283 319 325 377 395
425 433 445 479 577 593 719 911 1079 1367 1619 2051 2429 3077
```

Membership of each S(u) in U is checked by 53 exact integer substitutions.
Therefore any source n>3077 whose block visits this set descends. The last
boundary consists of the 1539 odd sources 1<=n<=3077. The explicit rows
`finite_boundary` in [controller.json](controller.json) give each source,
endpoint, step count and branch count. Recomputing them using C1's terminating
algorithm shows that every row with m>=3 descends. This is a bounded
arithmetic proof certificate for the remaining cases; it invokes no
unproved Collatz convergence claim.

Finally, for a word w define J=length(w), D=sum(w), and recursively
C_empty=0, C_(w,a)=3C_w+2^D. Its source cylinder is

\[
n\equiv(2^D-C_w)3^{-J}\pmod{2^{D+1}}, \tag{3}
\]

and its endpoint is (3^J n+C_w)/2^D. Reverse substitution through the word
shows that (3) enforces every specified valuation, not just the total D.
This yields the table and proves (2).

The cylinders are disjoint. Their odd Haar mass is
1/2+1/64+3/2048=1059/2048. This mass belongs to the original one-block
source law; it is not a fresh failure probability after each return.

## C3 A controller that preserves its original source

After a failed block reset the queue, but retain the original source n0
and concatenate the actual words. Test descent below n0 at each block
boundary. A terminal word and its earlier block cuts determine an exact
interval intersected with one cylinder (3). At a cut write
g=2^D-3^J. The endpoint is below n0 exactly when g>0 and n0>C/g.

Consequently an accepted word has lower guard
L=max(3,floor(C_final/g_final)+1). Every earlier cut with g>0 gives the
upper guard n0<=floor(C/g); earlier cuts with g<=0 fail automatically.
Take the minimum of these upper bounds, or infinity if there are none.
Together with (3), this is an exact certificate of the first successful
block boundary, preserving source, chronological cuts, carry and both clocks.
Different accepted words have disjoint guards because the native dynamics
and parser are deterministic.

[controller.py](controller.py) constructs such cells for every odd source
3<=n0<=10001. A separate arithmetic pass covers that entire finite universe
exactly once using only the cells' residue and interval guards. This finite
cover does not certify reference descent for all integers. C1 establishes
that R is total; convergence of its iterations to 1 is equivalent to
ordinary positive-odd Collatz convergence and remains OPEN.

## Reproduction and remaining target

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/controller.py --write
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/controller.py
```

The default extra hostile sweep compares the five-cylinder predicate with
native execution for every odd source 3<=n<=100001. It is additional
validation, not the justification for the all-size quantifier in C2.

The next precise problem is to control returns through the five expanding
cylinders while keeping the marked source. C2 supplies a finite **failure
menu**, not a closed five-state dynamical system: future residue precision
can increase. The old 81/27 tournament control is the analogous warning
that a short list of observed types need not preserve all transporters.
