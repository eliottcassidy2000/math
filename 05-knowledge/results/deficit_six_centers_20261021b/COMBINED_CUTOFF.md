# Combined deficit-six cutoff: order 35, ten residual orders

October 21, 2026. **PROVED from the explicitly inherited analytic inputs and
the completed order-ten certificate:** every strongly connected,
self-converse tournament with black path-orbit deficit six has order at most
**35**. The only still-possible orders above the established small centers are

\[
\boxed{12,14,15,18,20,21,24,25,28,35}.                       \tag{1}
\]

Every coefficient at these ten orders remains **OPEN**. Thus

\[
R_6(z)=z^5+z^6+z^7+O(z^{12}),\qquad \deg R_6\le35.          \tag{2}
\]

The associated trivalent paired-black counting series has minimal eventual
recurrence order **18**, with a uniform guaranteed onset **$n\ge48$**.
This is neither an earliest-onset claim nor a completed center census.

## 1. Source lineage and the operation being composed

The source is a literal strong self-converse tournament $T$. Let $H$ count
directed Hamilton words, $A=\operatorname{Aut}(T)$, $h=H/|A|$, and let $b$ count
the word orbits whose positional reversal is an anti-automorphism. Write
$\Delta=h-b$. We assume $\Delta=6$; the analytic obstructions also apply when
$\Delta\le6$.

Three separately scoped inputs are combined:

1. **PROVED, audited general closing-cycle theorem:**
   [cycle-cut obstruction, §5](../paired_black_leaf_reduction_20261020/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md)
   and its [explicit prefix source](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md).
   Although the first note's title concerns deficit two, its §5 explicitly
   states the general odd-order theorem. No deficit-two-only restriction is
   imported here.
2. **PROVED, independently reviewed incoming odd-shore and score-run bounds:**
   [current proof, §§1–3 and 8](../odd_shore_deficit_centers_20261021/session22-centers/PROOF.md),
   with the [third review](../odd_shore_deficit_centers_20261021/session22-centers-review/THIRD_REVIEW.md),
   received at commit `6b12a62f6b`. Its standalone cutoff is 47 because it
   retains the possible $q=1$ lane. Its finite controls are not center censuses.
3. **FINITE-EXACT:** this folder's frozen [original checker](check.py) and
   [original receipt](receipt.json) exclude all 128 necessary normalized
   order-ten frames. They also consume the inherited complete small tables
   and independently replay the three surviving sources. The original
   quadratic cutoff 248 remains a valid weaker theorem; its bytes are not
   rewritten by this follow-up.

The new composition uses input 1 to choose a Hamilton cycle to which input 2
applies with $q\ge2$, then uses input 3 to remove order ten. This preserves the
source and an actual Hamilton cycle. It does not replace the exact cut index
by a conveniently chosen period. The source's Hamilton ports and targets
remain in the frozen receipt; the combined order list alone forgets them.

## 2. Why a suitable cycle always has $q\ge2$

For a marked Hamilton cycle let $R$ be its successor, and put

\[
m=|A\cap\langle R\rangle|,\qquad q=n/m.
\]

Every tournament automorphism group has odd order, so $m$ is odd. The $n$
cycle cuts meet exactly $q$ path orbits: the unique vertex map taking one
ordered cut to another is a power of $R$. The number of blue cut orbits is
zero or $\gcd(2,q)\le2$. Hence

\[
q\le\Delta+2\le8.                                         \tag{3}
\]

At even order, oddness of $m$ already implies $q\ge2$. At odd order $n\ge7$,
take any Hamilton cycle. If $R$ is not an automorphism, its exact index is
already greater than one. Otherwise the cycle presents $T$ as an odd
circulant. The inherited closing-prefix theorem supplies another Hamilton
cycle with a nonblue cut. That new cycle's successor cannot be an
automorphism: an automorphic successor would give circulant coordinates in
which every positional reflection is anti, making every cut blue.

The closing coordinate is indispensable. In the inherited explicit families,
F1, F2 and F3 have return increments $3$, $k$ and $2$, respectively, each
belonging to the stated connection set. Unit-step reindexing preserves this
closure, and the four order-seven witnesses also close. Thus the theorem
supplies an actual cycle with a nonautomorphic successor, rather than merely
a nonblue bypass word. This reasoning is independent of the assumption
$\Delta=2$ used elsewhere in the older note.

This is an **existential choice of cycle**, not a statement about every
Hamilton cycle. The checker retains Paley7 as a hostile to that stronger
claim. Its standard cycle has $q=1$, while the closing word
$(0,1,2,4,5,6,3)$ has $q=7$ and is nonblue. Both use the same literal source.

Therefore every still-unknown center, whose order is at least ten, admits a
marked cycle with

\[
2\le q\le8,\qquad n=qm,\qquad m\text{ odd}.                \tag{4}
\]

## 3. Odd-shore balance and score runs reduce the parameter set

Here is the incoming mechanism, restated with the chosen-cycle quantifier.
The $q$ orbits of $R^q$ have odd size $m$; each induces a regular tournament.
The marked cycle supplies a directed perfect matching from each shore to the
next. If a proper subset of one shore is removed and $m\ge3$, the intact
shores are strong and their directed chain closes through any surviving
vertex of the altered shore. Hence the complement is strong. This is stronger
than an unmarked generic connectivity estimate.

Suppose $m\ge9$. Take a directed eight-cycle $S$ inside one shore and another
shore $B$. Every vertex of the first shore beats the same number $k$ of
vertices of $B$, so $\#(S\to B)=8k$. If every vertex in $B$ had exactly four
incoming arcs from $S$, the same total would be $4m$, forcing $2k=m$, impossible
for odd $m$. Some vertex therefore supplies at least five exits or at least
five entries. Entries are converted to exits by the actual anti-isomorphism
of the whole source. The strong complement supplies one common Hamilton
suffix. The half-supported eight-cycle cuts occupy at least eight black
orbits, including distinct converse mates, contrary to $\Delta\le6$.
Consequently

\[
m\in\{1,3,5,7\}.                                         \tag{5}
\]

The incoming score-run refinement applies when $q\ge5,m\ge5$. A full-score
class is a union of shores. If it consists of one shore, four-cycle cuts give
four black orbits with a four-vertex constant-score run absent from every
marked-cycle cut. If it consists of two, three or four shores, choose the
attachment terminal in another shore of that same score. The resulting run
has length five, while the cyclic shore-score word has runs of length at most
four. Converse mates have the corresponding complementary-score run, so they
too are disjoint from all cut orbits. Adding these four to the at least
$q-\gcd(2,q)$ black cut orbits gives $\Delta\ge8$.

Thus each full-score class has at least five shores. As $q\le8$, the whole
tournament must be regular. This is impossible for even $q=6,8$, because $m$
is odd and a regular tournament has odd order. For $q=7,m=7$, regularity also
gives a contradiction: regular connectivity is at least $\lceil49/3\rceil=17$,
so a removable sixteen-cycle yields at least eight entries or exits and
$\Delta\ge14$.

Combining (4), (5) and these exclusions gives exactly the following necessary
parameter rows at orders at least ten:

| $q$ | Allowed $m$ producing $n\ge10$ | Orders |
|---:|---|---|
| 2 | 5, 7 | 10, 14 |
| 3 | 5, 7 | 15, 21 |
| 4 | 3, 5, 7 | 12, 20, 28 |
| 5 | 3, 5, 7 | 15, 25, 35 |
| 6 | 3 | 18 |
| 7 | 3, 5 | 21, 35 |
| 8 | 3 | 24 |

In particular $n\le35$. These are possible parameter rows, not constructed
tournaments or isomorphism-class representatives. Regularity at the rows
$(q,m)=(5,5),(5,7),(7,5)$ is an additional necessary sidecar, not an existence
assertion.

## 4. The exact order-ten input and the ten open orders

Order ten forces $q=2,m=5$. In marked cycle coordinates there are four
circulant connection sets on each five-shore. The cross rule
$a_i\to b_j\iff f(j-i)=1$ has $f(0)=1,f(-1)=0$, leaving three free bits.
Thus every candidate lies in the 128-frame universe already completely
checked by the frozen producer. No extra reflection or self-converse filter
was imposed while generating that universe.

The producer enumerated all score-compatible maps for the full automorphism
and anti-isomorphism groups, used exact subset DP for $H$, and counted every
blue word once through its unique anti-involutive positional reversal. The
64 self-converse frames have defects
$1528,1628,2908,2926,3064,3114$, with multiplicities $8,8,12,12,12,12$.
None has defect six. The source and literal word/target sidecars remain in
the original receipt.

The new checker pins that producer and receipt, reconstructs every normalized
source from its seven bits, and replays its stored maps and native ports.
It **consumes** the previously proved exhaustive counts; it does not label
this lightweight replay as a new census. Removing ten from the preceding
table gives exactly (1). The inherited complete small tables supply one
center at each of orders five, six and seven, and none at eight or nine.
Equation (4) excludes eleven. This proves (2), while leaving every coefficient
listed in (1) undetermined.

## 5. Minimal recurrence order eighteen, safe onset forty-eight

Put $p=1-z-z^3$, $u=p(z^2)$ and $A=1/p$. The inherited SCC/center law and the
proved three-presentation classification give

\[
G(z)=\frac{z^5}{2}\bigl(A(z)^2-A(z^2)\bigr)
      +A(z^2)R_6(z)+z^{10}A(z^2)^2(1+z+z^3).
\]

Over $D=p^2u^2$, its numerator is

\[
N=\tfrac12z^5(u^2-p^2u)+p^2uR_6+z^{10}p^2(1+z+z^3).       \tag{6}
\]

Since $\deg R_6\le35$, $\deg N\le47$. The denominator has degree eighteen,
constant coefficient one, and the homogeneous coefficient recurrence is
valid at every $n\ge48$.

There is no denominator cancellation for any actual $R_6$. The factors $p,u$
are squarefree and coprime. At a root of $p$, $u=2z^4\ne0$, and (6) is
$z^5u^2/2\ne0$. At a root of $u$, it is
$z^{10}p^2(1+z+z^3)\ne0$: simultaneous vanishing of the last factor would again
force $2z^4=0$. Thus the minimal eventual recurrence order is exactly eighteen.
This denominator fact was inherited; the new combined center bound improves
the guaranteed onset.

The exact formal control $R_6=z^5+z^6+z^7+z^{35}$ has numerator degree 47 and
leading coefficient $-1$. The homogeneous relation fails at 47 and holds
from 48. The top coefficient is **synthetic** and is not a claimed order-35
center. Therefore the control proves the sharp uniform algebraic bill for
the stated degree bound, not the earliest recurrence onset of the actual
counting sequence.

## 6. Reproduction and status

```text
python -B 05-knowledge/results/deficit_six_centers_20261021b/combined_check.py --verify
python -B -O 05-knowledge/results/deficit_six_centers_20261021b/combined_check.py --verify
```

The checker uses always-active exact checks and records LF SHA-256 pins for
the incoming proof/review, both closing-cycle sources, the inherited SCC
formula inputs, and the frozen local producer/receipt. Its 16,889 gates cover
those pins, every stored order-ten source/map/port, the exact $(q,m)$ support,
the Paley7 cycle-choice hostile, and formal recurrence controls. This is a
lightweight **certificate consumer and algebra check**, not a census of the
ten remaining orders.

Normal and optimized runs agree with [combined_receipt.json](combined_receipt.json),
whose LF SHA-256 is
`9a1bdb64fdc717cb3af1aa9c2c65eef5b59978786f161969fbe1147411d8e048`.
The original producer and receipt are unchanged; the README routes the
strengthened result. Independent root review of the complete proof and
certificate-consumer API passed. The root's independent optimized
`-B -O --verify` replay also passed all 16,889 checks and matched the saved
receipt hash above. No remaining-order census is asserted.
