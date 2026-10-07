# Explicit two-anchor exponent phases and a smaller-child exit

**Status: PROVED** for the parameterized guards, exact common-future receipts,
strict smaller-child payment and declared entry-bank comparison below.
**FINITE-EXACT** for the listed implementation controls. **CONDITIONAL** for
ROOT completion: the smaller child still needs a supplied first-hit certificate.
This is not a universal Collatz proof or an audited promotion of every claim in
the incoming two-anchor papers.

The inherited mechanism is the affine clearing and head identity in
[THM-4601, two-anchor reduction](../../01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md),
parts (iii)–(iv). We independently derive exactly that template here. Its two
Mersenne examples, \(K=1889,5249\), are inherited. The new transfer is a uniform
phase for both shell parities at every \(J\ge3\), its exact comparison with the full entry bank of
[three inverse generators](collatz_child_normal_forms_20261007.md), and the
source-preserving composition of the resulting exits. The prior
[mixed-child payment compiler](collatz_mixed_child_rules_20261007.md) explains
why a requested deletion depth must be authenticated, not treated as an
available budget. The canonical hostile is the distinction between a paid
dependency and a ROOT proof. The least-used sidecar here is the exponent's
ternary residue, which distinguishes old entry coverage inside a binary phase.

## 1. The consumed affine identity, checked independently

Write \(F_w(x)=(3^{|w|}x+B_w)/2^{\sum w}\), with words in chronological
order. Let \(n=2^K-1\). Its first \(K-1\) odd valuations are \(1\), reaching

\[
x=2\,3^{K-1}-1.
\]

The deletion children \(h_D=2^{K-D}-1\), \(D=3,4\), similarly reach

\[
y_D=(x+1)/3^D-1.
\]

The relevant exact carriers are

\[
F_{222}(x)=\frac{27x+37}{64},\qquad
F_{411}(y)=\frac{27y+89}{64},\qquad
F_{2211}(y)=\frac{81y+143}{64}.
\]

Consequently both child clearings have value \((x+63)/64\), and

\[
F_{222}(x)-1=27\bigl(F_{411}(y_3)-1\bigr)
 =27\bigl(F_{2211}(y_4)-1\bigr).
\]

The map \(F_2(z)-1=3(z-1)/4\) preserves this relation. After \(J\ge3\)
source twos and the remaining \(J-3\) child twos, let the values be (X,Y).
Then \(X-1=27(Y-1)\). Finally,

\[
F_{10}(X)=\frac{3X+1}{1024},\qquad
F_{3113}(Y)=\frac{81Y+179}{256}=4F_{10}(X)+1.
\]

Put \(Z=F_{10}(X)\) and \(c=v_2(3Z+1)\). The identity

\[
3(4Z+1)+1=4(3Z+1)
\]

gives the common-future words

\[
\begin{aligned}
w&=1^{K-1}2^J\,(10,c),\\
v_3&=1^{K-4}(4,1,1)2^{J-3}(3,1,1,3,c+2),\\
v_4&=1^{K-5}(2,2,1,1)2^{J-3}(3,1,1,3,c+2).
\end{aligned}
\]

Both sides have \(K+J+1\) odd edges. The source's total valuation cost exceeds
the child's by \(D\). These formulas are affine identities; the following
guard is what makes them actual positive first-hit receipts.

## 2. A complete explicit phase for every \(J\ge3\)

Let \(e=K-1>0\). For even \(e\), the elementary lifting identity is

\[
v_2(3^e-1)=2+v_2(e).
\]

If \(v_2(e)=2J-1\), then \(v_2(x-1)=2J+2\). There are exactly \(J\)
successive twos after the ones-run, followed by a valuation at least \(3\).
The endpoint is

\[
X=1+\frac{3^J(3^e-1)}{2^{2J-1}}.
\]

The next valuation is exactly \(10\) if and only if

\[
3^{K-1}\equiv
R_J:=1+255\,2^{2J+1}3^{-J-1}
\pmod {2^{2J+10}}. \tag{1}
\]

Indeed, the numerator of \(3X+1\) is

\[
3^{J+1}(3^e-1)+2^{2J+1},
\]

and it must be \(2^{2J+9}\) modulo \(2^{2J+10}\). This proves both directions
of equation (1), including the extra bit that excludes a valuation greater than \(10\).

For \(b\ge3\), the powers of \(3\) modulo \(2^b\) are exactly the residues

\[
u\equiv1\text{ or }3\pmod8,
\]

and the order is \(2^{b-2}\). To see this, the displayed valuation identity
gives that order; every power belongs to the indicated subgroup, which has
the same size. Here \(R_J\equiv1\pmod8\), so there is a unique phase

\[
\boxed{K\equiv k_J\pmod {P_J},\qquad P_J=2^{2J+8}.} \tag{2}
\]

Moreover \(v_2(R_J-1)=2J+1\), so equation (1) itself forces

\[
v_2(K-1)=2J-1.
\]

There is no separate shell assumption to guess. The implementation finds the
unique phase by testing the two exponent lifts at each successive binary bit.
It uses only modular powers, never a real logarithm or the enormous source.

| \(J\) | Least \(k_J\) | Period \(P_J\) |
|---:|---:|---:|
| 3 | 1889 | 16384 |
| 4 | 5249 | 65536 |
| 5 | 50689 | 262144 |
| 6 | 67585 | 1048576 |
| 7 | 90113 | 4194304 |
| 8 | 11304961 | 16777216 |

The period is minimal, not merely sufficient. Each phase occupies exactly

\[
\frac{1/2^{2J+8}}{1/2^{2J}}=\frac1{256}
\]

of its shell \(v_2(K-1)=2J-1\). No claim is made that this one head covers
the other \(255/256\).

### The even-shell companion closes the missing parity

The second inherited head pair is

\[
u=(1,10),\qquad v=(1,1,2,2,3).
\]

Its exact carriers are

\[
F_u(X)=\frac{9X+5}{2048},\qquad
F_v(Y)=\frac{243Y+283}{512}=4F_u(X)+1
\]

when \(X-1=27(Y-1)\). Use these heads in place of \((10)\) and
\((3,1,1,3)\) in Section 1, leaving both clearing words unchanged.

For \(v_2(K-1)=2J-2\), the two-run has length \(J\), then its next valuation
is \(1\). The following valuation is \(10\) exactly when

\[
9X+5\equiv2^{11}\pmod {2^{12}}.
\]

Substituting the same formula for \(X\) gives the unique phase

\[
\boxed{\begin{aligned}
3^{K-1}&\equiv1+1017\,2^{2J}3^{-J-2}\pmod {2^{2J+11}},\\
K&\equiv k^{(1)}_J\pmod {2^{2J+9}}.
\end{aligned}} \tag{2a}
\]

The numerator is \(3^{J+2}(3^{K-1}-1)+7\,2^{2J}\); thus the coefficient
\(1017=1024-7\) is exact. Lifting again forces
\(v_2(K-1)=2J-2\), so this is an iff native head phase without a separate
shell guess. The least residues for \(J=3,4,5\) are \(6129,52545,517889\),
with periods \(32768,131072,524288\). It occupies exactly \(1/1024\) of its
even shell. Together (2) and (2a) provide an explicit paid head phase in
**every shell \(v_2(K-1)=m\ge4\)**. Their relative sizes are stated separately.

## 3. Positivity, exact child words and first ROOT

Every positive exponent in phase (2) has \(e\ge2^{2J-1}\ge32\). The smallest
source-prefix slope after the ones-run, twos and head is

\[
\frac{3^{e+J+1}}{2^{e+2J+10}}>1. \tag{3}
\]

For the least possible \(e=2^{2J-1}\), verify inequality (3) at \(J=3\) by the exact
integer inequality \(3^{36}>2^{48}\). Increasing \(e\) increases the ratio.
Passing from \(J\) to \(J+1\) at these least values multiplies it by

\[
\frac34\left(\frac32\right)^{3\cdot2^{2J-1}}>1.
\]

Every nonempty source prefix through \(Z\) consequently has slope greater
than \(1\) and positive carry, so every one is greater than the immutable
source \(n\). In particular, \(Z>n>1\). The final join may descend; we do not
claim that it also grows.

The affine child prefix ends at the positive odd integer \(4Z+1>1\).
For a word with positive valuation letters, a nonintegral dyadic intermediate
can never become integral later: \(v_2(3z+1)=v_2(z)<0\) when \(v_2(z)<0\).
An even integer intermediate makes the next value nonintegral. Thus its odd
integer final value forces every intermediate to be an integer of the correct
odd parity and all stated valuations to be exact. Positivity also follows
from each affine map's positive slope and carry applied to \(h_D>0\).

Such a valid child prefix cannot have visited \(1\), since an actual odd step
at \(1\) stays at \(1\), whereas its endpoint \(4Z+1\) is larger. The final

\[
U(Z)=U(4Z+1)
\]

therefore produces strict first-hit receipts, including the case when this
common endpoint is \(1\). Both \(h_3\) and \(h_4\) are positive and strictly
smaller than \(n\). The preferred stronger exit is \(h_4\).

For the even-shell companion, \(e\ge2^{2J-2}\ge16\) and the last head slope
is \(3^{e+J+2}/2^{e+2J+11}>1\). Its base case is the exact inequality
\(3^{21}>2^{33}\), and the same induction applies. All positivity and
first-hit arguments therefore hold unchanged. Each side now has
\(K+J+2\) odd edges, and the cost difference is still \(D\).

Completion is a separate typed operation. Given an independently supplied
first-hit ROOT word for \(h_D\), its initial segment must be \(v_D\). Replace
that segment by \(w\) and keep its remaining suffix. `discharge` authenticates
the child word, the common future, and the resulting first-hit source word.
It does not obtain a missing child proof by orbit search.

## 4. Exact coverage beyond a declared entry bank

The baseline is all nonempty words, of arbitrary length, over

\[
G_1(n)=\frac{2n-1}{3},\quad
G_5(n)=\frac{8n-5}{9},\quad
G_{17}(n)=\frac{2048n-2363}{2187},
\]

with their native positive odd integer guards. This is **not** a comparison
against every historical selector or against the full Collatz graph. A longer
word cannot start unless its first generator is legal.

For \(n=2^K-1\):

* \(G_1\) never starts: its source must be \(2\pmod3\), whereas a Mersenne
  source is \(0\) or \(1\pmod3\).
* \(G_5\) starts exactly when \(K\equiv5\pmod6\).
* \(G_{17}\) starts exactly when \(K\equiv733\pmod{1458}\).

For the last assertion, the native congruence is

\[
2^{K+11}\equiv37\pmod {2187}.
\]

Direct modular evaluation gives \(2^{744}\equiv37\), and \(2\) has order

\[
2\cdot3^6=1458
\]

modulo \(3^7\), by \(v_3(4^j-1)=1+v_3(j)\). These entry classes are disjoint.

Within any phase (2) or (2a), its \(729\) parameter residues give all odd residues of

\[
K\pmod {1458},
\]

because \(P_J\) is a power of \(2\). Exactly \(243\) admit \(G_5\), exactly
one admits \(G_{17}\), and the remaining \(485\) admit neither. Therefore

\[
\boxed{\frac{485}{729}}
\]

of **each compiled head phase** is outside the entire old entry bank, yet has
the paid deletion receipt of Sections 1–3. The same old entry split holds
inside every shell \(v_2(K-1)=m\ge4\). However, the guarantee for the new
single head is only

\[
\frac{485}{729\cdot256}=\frac{485}{186624}
\]

of each odd shell \(m=2J-1\ge5\). All these are ordinary relative natural
densities in explicit arithmetic progressions of **exponents** \(K\), not
densities of integer sources \(2^K-1\).

The even-shell guarantee is \(485/(729\cdot1024)=485/746496\) of each
\(m=2J-2\ge4\). The factor \(485/729\) is the within-head-phase exclusion
fraction in both cases, not the fraction of an entire shell newly covered.

The inherited points \(1889,5249\) are already \(G_5\)-covered. A genuinely
new point for the declared baseline is \(K=18273\). The whole progression

\[
\boxed{K=18273+11943936s,\qquad s\ge0} \tag{4}
\]

has \(J=3\), satisfies the head guard, and avoids all old entry words.
This is \(K=1889+16384(1+729s)\); hence both the binary phase and every old
entry predicate are fixed.

The even-shell least phase already has an uncovered point \(K=6129\).
Its independently preserved excluded progression is

\[
\boxed{K=6129+23887872s,\qquad s\ge0.} \tag{4a}
\]

## 5. An actual change of residual type, and a composable further exit

Every phase (2) or (2a) has \(K\equiv1\pmod {16}\). Its four-bit child has exponent

\[
E=K-4,\qquad v_2(E-1)=v_2(K-5)=2.
\]

Therefore the authenticated exit sends an arbitrarily long two-run \(J\) to
a Mersenne child with exactly **two** twos after its initial ones-run. This
normalizes the run type while decreasing the source; it does not solve the
remaining two-twos type.

On either excluded progression (4) or (4a), \(K\equiv3\pmod6\), so \(E\equiv5\pmod6\).
The four-bit child now admits \(G_5\), producing the still smaller integer

\[
g=G_5(2^{K-4}-1)=\frac{2^{K-1}-13}{9}.
\]

Its actual forward word ((1,2)) reaches the four-bit child. Prepending this
word to \(v_4\) yields a source-preserving common-future receipt from the
original \(n\) to \(g\). This composition does not recycle a route back to
the original \(n\): it returns to the already strictly smaller \(h_4\).
ROOT completion still requires a supplied proof for \(g\).

The restricted \(+1\) barrier in
[HYP-9240](../hypotheses/HYP-9240-plus-one-barrier-for-all-deletion-depths.md)
concerns deletion-clock collisions during the two-run. It does not exclude
all other arithmetic reductions. For example, \(G_5(1)=1/3\) at the rational
anchor, and \(F_{12}(1/3)=1\); its integer native guard pays one third of
every long-run shell before reading the run. Conversely,

\[
G_{17}(1)=-35/243
\]

is negative at that rational limit even though every positive odd native
integer source has a positive child. Positivity of an excluded rational
anchor is not positivity on the actual guarded integer domain. Neither
observation refutes the restricted deletion-clock conjecture.

The remaining obligation is explicit: to ground these new families, prove
ROOT for their actual \(h_4\) or \(g\) children, or supply their individual
certificates. A statement about an entire child family is an infinite schema,
not a finite list of assumed seeds.

## 6. Exact implementation and reproduction

The owned script is
[collatz_completion_anchor_20261007b.py](../../04-computation/experiments/collatz_completion_anchor_20261007b.py).
It exports `phase`, `contains`, modular `source_residue`, bounded
`materialize`, `g5_extend`, and supplied-certificate `discharge`. Canonical
phase fields, exact integer types, valuation words, source ownership, strict
payment and first-hit boundaries are authenticated. Symbolic readers keep
their source parameter; they never substitute another member of the phase.

The finite controls comprise:

* both clearing carriers and the head identity, independently checked by
  exact rational substitution;
* both branches at \(J=3,\ldots,32\), with five selected parameters each, including a large
  symbolic parameter; all \(729\) old-entry parameter residues per phase;
* independent old native-guard counts in shells \(m=4,\ldots,20\);
* full literal receipts for \(K=1889,5249,18273,6129\), for both deletion depths;
* composed \(G_5\) receipts at \(K=18273,6129\), and \(128\) modular controls
  along each of (4) and (4a);
* a one-less-bit hostile, forged phase and receipt fields, malformed numeric
  types, an explicit materialization cap and ROOT-padding rejection.

No selected child ROOT orbit is discovered in this experiment. The small
transport control uses the supplied empty proof of ROOT and the explicit
word \(3\to5\to1\); it checks the API boundary, not new completion coverage.

```powershell
python -B 04-computation/experiments/collatz_completion_anchor_20261007b.py
python -B -O 04-computation/experiments/collatz_completion_anchor_20261007b.py
```

Both runs must equal the saved
[output](collatz_completion_anchor_20261007b.out). No assertions, randomness,
floating logarithms, hidden convergence oracle or import-time census is used.
Normal and optimized runs agree on 1,008,509 explicit checks. Their normalized
LF output SHA256 is
e38bd4926e3a667f83a240ecc7f936aad3f8633e2b7a40bddc8d209570d14235.
