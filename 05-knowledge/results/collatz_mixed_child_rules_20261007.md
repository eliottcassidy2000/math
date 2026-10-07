# Paid exits after the smaller child leaves the Mersenne family

Status: **PROVED** guarded inverse-prefix composition, sharp six-bit deletion
budget, and displayed infinite mixed phases. **FINITE-EXACT** native-cell
controls and selected point certificate transport. Universal closure under
arbitrary children remains **OPEN**.

The new step is a genuine smaller obligation for the mixed source. Returning
from that source to a larger old Mersenne state is only an intermediate
prefix. Composing that prefix with the previously proved six-bit deletion
can finish at an obligation below the mixed source. This works after at most
35 inverse-five blocks, with an exact all-positive failure at 36. The
concrete stopped inverse-chain child arising from exponent 1459 is also paid
and fully grounded using the already stored exponent-1451 certificate.

## 1. Inheritance, baseline and retained state

The closest proved mechanisms are the variable-depth collision in
[reset-two rules](collatz_reset2_rules_20261007.md), the inverse negative-five
branch in [collision DP](collatz_collision_dp_20261007.md), and immutable-source
payment in [adaptive switch payment](adaptive_switch_payment_20261004.md).
These are applications of the guarded common-future interface of
[THM-4555, uniform switches are collisions at minus one](../../01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md).
The inherited inverse-one rule also appears in
[complement routing, section 5](collatz_complement_routing_20261004.md).
No novelty claim is made for either elementary inverse branch.

The active objects are the original source, its smaller dependency, the
actual connecting word, ternary divisibility, the affine intercept, and the
first-hit ROOT suffix. The canonical hostile is to compare a final child
only with the larger intermediate state. That comparison does not prove
the child is smaller than the current mixed source. The decisive control
below is the otherwise identical 35-block versus 36-block composition.

Let \(U(n)=\operatorname{oddpart}(3n+1)\), stopping at the first 1. Define
the guarded inverse operations

\[
 I_5(x)=\frac{8x-5}{9},\quad x\equiv4\pmod9,
 \qquad J(x)=\frac{2x-1}{3},\quad x\equiv2\pmod3.
\]

They are positive odd integers smaller than their positive odd inputs, and
their actual return words are respectively \((1,2)\) and \((1)\).
For repeated applications the exact guards and formulas are

\[
 I_5^r(x)=8^r\frac{x+5}{9^r}-5,
 \quad 9^r\mid x+5;
 \qquad
 J^e(x)=2^e\frac{x+1}{3^e}-1,
 \quad 3^e\mid x+1.                               \tag{1}
\]

Here \(r,e\ge0\), and positivity is retained. The guards ensure all
intermediate states are positive odd integers. These operations transport
a supplied proof, or produce a strictly smaller dependency; they do not
create an unsupplied ROOT certificate.

## 2. Compose first, then compare with the mixed source

Suppose there is a checked common-future receipt from \(m\) to

\[
 q=\frac{m+1}{2^D}-1<m,\qquad D\ge1.              \tag{2}
\]

Both actual words and their identical endpoint are part of this premise.
Write \(X=m+1\). For a legal inverse-five depth \(r\ge1\), let

\[
 h_r=I_5^r(m)=\frac{8^r(X+4)}{9^r}-5.
\]

The source \(h_r\) follows \((1,2)^r\) to \(m\). Prepend that word to
the source side of the receipt (2), retaining the same child-side word.
This gives an exact common future. It is a **paid** receipt relative to
\(h_r\), rather than merely relative to \(m\), if and only if

\[
 (2^D8^r-9^r)X>2^{D+2}(9^r-8^r).                \tag{3}
\]

This follows by subtracting \(q\) from \(h_r\), with denominator
\(2^D9^r\). All source identities and actual-word guards are already
required, so (3) is a payment test, not an independent legality test.

If \(2^D8^r\le9^r\), payment is impossible at every positive \(X\).
Otherwise the least integer height that passes is exactly

\[
 X_{\min}=1+\left\lfloor
 \frac{2^{D+2}(9^r-8^r)}{2^D8^r-9^r}
 \right\rfloor.                                  \tag{4}
\]

The strict inequality and intercept matter at its boundary. The function
`payment_threshold` implements (4), with `None` for the impossible case.
`inverse_prefix_exit` independently checks the connecting source and actual
word, then compares the child with this new source. It returns an unpaid
result rather than recycling the old comparison.

### The sharp six-bit budget

For the preceding deep collision, \(D=6\). Exact integer arithmetic gives

\[
 64\,8^{35}>9^{35},\qquad64\,8^{36}<9^{36}.
\]

The threshold (4) increases with \(r\) while the coefficient is positive:
putting \(t=(9/8)^r\), it is the strict threshold
\(256(t-1)/(64-t)\), whose derivative is positive on \(1<t<64\).
At \(r=35\) its least integer is 6780. Consequently, for every
\(X\ge8192\),

\[
 h_r>q\quad\Longleftrightarrow\quad1\le r\le35.  \tag{5}
\]

For every \(r\ge36\) the coefficient in (3) is negative and its right side
is positive, so failure holds for **all** positive \(X\), not only at sampled
sources. This is a sharp limit of the given deletion, not an obstruction to
using a different deletion, a different common future, or a ROOT proof.

## 3. A family that returns to the Mersenne type with a smaller obligation

Impose the stronger binary phase

\[
 K\equiv1459\pmod{2^{126}},\qquad K\ge1459.       \tag{6}
\]

Then \(m=2^{K-2}-1\) is in the preceding deep rule's phase
\(K-2\equiv1457\pmod{2^{126}}\). That rule supplies (2) with
\(q=2^{K-8}-1\). For any \(1\le r\le35\) with
\(9^r\mid2^{K-2}+4\), (3) is satisfied. Thus

\[
 h_r=\frac{8^r(2^{K-2}+4)}{9^r}-5
       \rightsquigarrow 2^{K-8}-1<h_r.             \tag{7}
\]

The source word is \((1,2)^r\) followed by the exact deep-collision source
word at \(m\); the child word is that collision's original child word.
All proper prefixes avoid ROOT, since the connecting prefix ends at
\(m>1\), and both collision words already have first-hit typing.

These are infinite families even when the inverse-five depth is maximal.
For each \(r\ge1\) and \(u\in\{1,2\}\), choose

\[
 K\equiv1459\pmod{2^{126}},\qquad
 K\equiv4+u3^{2r-1}\pmod{3^{2r}}.                 \tag{8}
\]

CRT supplies infinitely many positive exponents. Since \(K\) is odd,
elementary three-adic lifting yields

\[
 v_3(2^{K-2}+4)=1+v_3(K-4)=2r.
\]

Hence (8) has exactly \(r\) legal inverse-five turns. The paid part of this
family is precisely \(r\le35\) under the fixed six-bit deletion. No pointwise
ROOT claim for every \(q\) in (7) is included.

The full inverse child has \(v_2(h_r+1)=2\) and
\(v_2(h_r+5)=3r+2\): it leaves the Mersenne initial-one-run state and has
exactly \(r\) initial \((1,2)\) blocks. Those blocks alone only return to
the larger \(m\). Equation (7), including its smaller final obligation,
is the additional inference.

### Independent finite controls beyond astronomical exponents

The full deep native cell at \(m_0=2^{1457}-1\) contains

\[
 m\equiv m_0\pmod{2^{1585}}.
\]

The 1585 bits retain 1456 initial ones and the deep head of cost 128,
including its last oddness bit. Combining this cell with
\(m+5\equiv u9^r\pmod{3^{2r+1}}\) supplies finite-bit exact controls for
every tested \(r\). The experiment replays both sides at \(r=1,\ldots,36\),
both units, and two parameter lifts: 140 paid cases and four genuine
unpaid cases at 36. This checks the whole native-word interface, not a
rounded comparison of enormous powers.

`MixedPlan` separately stores (6)–(8) without expanding \(2^K\). For a
requested modulus \(M\), the source reader reduces \(2^{K-2}\) modulo
\(9^rM\) **before dividing by** \(9^r\), retaining the required precision.
The child is read as \(2^{K-8}-1\pmod M\). Materialization has an explicit
bit cap. Modular reads do not themselves supply missing ROOT proofs.

## 4. A stopped inverse chain also gets an exit

At \(K=1459\), put \(X=2^{1457}\) and \(m=X-1\). Maximal alternating
applications of the two elementary inverse rules give

\[
\begin{aligned}
 m&\xrightarrow{I_5}h=(8X-13)/9,\\
 h&\xrightarrow{J^5}g=(256X-2315)/2187,\\
 g&\xrightarrow{I_5}z=(2048X-29455)/19683.           \tag{9}
\end{aligned}
\]

The exact guards are
\(v_3(m+5)=2\), \(v_3(h+1)=5\), and \(v_3(g+5)=3\).
The final state satisfies \(z\equiv1\pmod9\), so neither \(J\) nor
\(I_5\) applies there. This statement concerns those two inverse guards,
not all possible rules. The states modulo 2187 are respectively
1093, 242, 1885, and 1918.

The whole chain and its exact stop persist for

\[
 K=1459+2^{126}3^{10}s,\qquad s\ge0,              \tag{10}
\]

with \(X=2^{K-2}\). The order of 2 modulo \(3^{11}\) is
\(2\cdot3^{10}\), so this phase preserves all numerator residues needed
through the nine ternary divisions and two final guard digits. Reducing
successively after each exact division proves the guards from the finite
base residue; the script checks this same precision schedule independently.

The actual connecting word from \(z\) to \(m\) is

\[
 W=(1,2),1^5,(1,2),\qquad
 F_W(z)=(19683z+27407)/2048=m.
\]

Again, following this word alone does not pay against \(z\). But composing
with the deep deletion gives the smaller obligation \(q=X/64-1\), because

\[
 z-q=\frac{111389X-625408}{64\cdot19683}>0
 \quad\text{for every }X\ge6.                    \tag{11}
\]

All sources in (10) satisfy this height condition. Thus a frontier where
the two simple inverse operations stop is not a frontier of the composed
receipt language. The proof changes the final obligation; it does not
declare the stopped state solved just because its forward path is known.

## 5. Any mixed inverse word has the same exact payment interface

The concurrent [child normal-form package](collatz_child_normal_forms_20261007.md)
supplies three guarded inverse generators, in application order:

\[
 G_1(m)=(2m-1)/3,\quad
 G_5(m)=(8m-5)/9,\quad
 G_{17}(m)=(2048m-2363)/2187.
\]

Their actual return words are \((1)\), \((1,2)\), and
\((1,1,1,2,1,1,4)\). A generated word \(w\), with binary and ternary
costs \(A,B\), has the exact carrier

\[
 h=\frac{q m-C}{p}=\frac{qX-D_0}{p},\qquad
 q=2^A,\ p=3^B,\ D_0=q+C,
 \quad m\equiv Cq^{-1}\pmod p.                   \tag{12}
\]

The normal-form theorem guarantees that this native congruence is sufficient
for every intermediate positive odd state. The actual return word is the
concatenation of the generator words in **reverse application order**.
Prepending that return word to the deletion receipt (2) gives

\[
 h>q_{\rm child}=X/2^D-1
 \quad\Longleftrightarrow\quad
 (2^Dq-p)X>2^D(D_0-p).                            \tag{13}
\]

This is the general counterpart of (3). All source and native guards remain
premises. The exact least integer height is

\[
 1+\left\lfloor\frac{2^D(D_0-p)}{2^Dq-p}\right\rfloor
 \quad\text{when }2^Dq>p,                         \tag{14}
\]

and no positive height can pay when \(2^Dq\le p\). To see the latter claim
and its scope, the carriers obey the stronger envelope

\[
 p-q\le C\le17(p-q),\qquad 0\le D_0-p\le16p.    \tag{15}
\]

For one generator its carry is \(g(p_g-q_g)\), \(g\in\{1,5,17\}\).
On appending that generator,

\[
 p'-q'=q_g(p-q)+(p_g-q_g)p,\qquad
 C'=q_gC+g(p_g-q_g)p.
\]

These preserve (15). Equality \(D_0=p\) holds exactly for words containing
only \(G_1\), including the empty word. All nonempty words compatible with Mersenne
inputs begin with \(G_5\) or \(G_{17}\), so their right side in (13) is
strictly positive. The unfulfilled budget case therefore cannot be repaired
merely by increasing the source height.

The necessary coefficient budget \(p/q<2^D\) is additive in logarithms:
\(B\log 3-A\log 2<D\log 2\). Production uses exact integer comparison,
not rounded logarithms. That budget forgets order; the carry and native
source phase restore information it loses. For example, the application
words \((17,5^{34})\) and \((5^{34},17)\) both have
\((A,B)=(113,75)\), but their exact six-bit height cuts are respectively
2726 and 3244. Equality of the cost vector is not equality of the paid rule.

### Every funded mixed word gives an infinite deep-phase family

Fix \(D=6\). The least ratio \(p_g/q_g\) of the generators is
\(2187/2048\), attained by \(G_{17}\). Exact integer comparisons give

\[
 2187^{63}<64\,2048^{63},\qquad
 2187^{64}>64\,2048^{64},\qquad
 9\,2187^{62}\ge64\cdot8\,2048^{62}.
\]

Thus every word passing the six-bit budget has length at most 63, and the
only one of length 63 is \(G_{17}^{63}\). Its height cut is 45594.
This is a bound for a fixed deletion budget, not a bound on the depth of a
controller that later obtains another deletion.

Every nonempty word beginning with \(G_5\) or \(G_{17}\) has exactly one
native Mersenne exponent phase
\(E\equiv E_w\pmod{2\cdot3^{B-1}}\), with \(E_w\) odd, by the
normal-form theorem. Intersect this with \(E\equiv1457\pmod{2^{126}}\).
The common parity is the sole shared factor, so CRT gives an infinite phase
of period \(2^{126}3^{B-1}\). Set \(K=E+2\), \(m=2^E-1\), and retain
the actual source \(h\) in (12). Whenever the budget passes, (13) supplies

\[
 h\rightsquigarrow2^{E-6}-1<h.                  \tag{16}
\]

The exact cut (14) is retained by the compiler. In fact every such cut is
automatically met on the declared deep phase \(E\ge1457\): from (15),
the positive integer gap \(64q-p\ge1\), and \(B\le7\cdot63=441\),

\[
 X_{\min}\le1024\,3^{441}+1<2^{710}<2^{1457}.
\]

This proves an infinite paid family for **every** mixed word that begins
with a compatible generator and passes the fixed budget. It does not say
every supplied integer lies in one of their guards, or that the child in
(16) is already grounded.

`word_payment_threshold`, `mixed_word_exit`, `word_exit_phase`, and
`mixed_word_residues` implement these statements. A symbolic exponent
phase is not an expanded actual word: its collision side contains
\(E-1\) initial ones, which can be enormous. The new modular reader
retains the full factor \(p\) before division, and no huge Mersenne integer
is created in its phase controls. The earlier explicit materialization cap
remains in force.

Sixteen selected mixed words test pure-branch boundaries,
different orders, the stopped-chain word, and the 63/64-letter boundary.
Two independently generated native lifts per word give 26 paid and six
unpaid receipts. These are exact finite controls of the all-length proof,
not an exponential enumeration of all qualifying words.

## 6. Grounded points and precise remaining obligation

For \(K=1459\), the child \(q=2^{1451}-1\) already has a frozen first-hit
certificate in the preceding package. Supplying that certificate discharges
both (7) at \(r=1\) and (11):

| Source | Odd ROOT rank | Total valuation cost |
|---|---:|---:|
| \(h=(8\cdot2^{1457}-13)/9\) | 7349 | 13105 |
| \(z=(2048\cdot2^{1457}-29455)/19683\) | 7356 | 13113 |

Every exponent is replayed, the first ROOT is enforced, and the compiled
words are also compared with the independently authenticated connecting
prefix followed by the old 1457 word. No new ROOT orbit is discovered by
this package. This identity check does not confuse the larger intermediate
with the smaller obligation used by the new receipt.

For arbitrary parameters, the remaining obligation is still
\(\operatorname{Root}(2^{K-8}-1)\), or another valid smaller-child rule at
that source. Nor does the fixed six-bit exit handle maximal inverse-five
depths 36 and above: the exact failure (3) requires a stronger deletion or another pattern.
This is a proved partial closure mechanism with a quantified boundary,
not universal closure under the children of the preceding rules.

## 7. Reproduction and adversarial controls

```text
python -B 04-computation/experiments/collatz_mixed_child_rules_20261007.py
python -B -O 04-computation/experiments/collatz_mixed_child_rules_20261007.py
```

[Script](../../04-computation/experiments/collatz_mixed_child_rules_20261007.py)
and [saved output](collatz_mixed_child_rules_20261007.out) use exact integer
arithmetic and make no default file writes. The declared controls include
the 144 pure-five native-cell cases, 32 mixed-word native lifts,
exact phase/residue checks through depth 40 with
two ternary units and three lifts, 32 symbolic stopped-chain phases, and the
two supplied-certificate discharges. There is no further orbit-discovery
control or hidden terminal oracle.

Invalid integer types, native-guard failures, wrong connecting sources,
insufficient materialization budgets, and truncated child certificates are
rejected. The 36-block cases deliberately retain the exact common future
while refusing to label it a smaller-child receipt. This separates an
algebraically valid composition from its required original-source payment.
