# Guarded inverse debt: a descent certificate and an exact clock boundary

**PROVED:** the positive-word debt criterion, cycle budget, and guarded
source family in §§2–4. **INHERITED:** the reachability distinctions in
§1. **FINITE-EXACT:** the accompanying small controls. Classical positive
Collatz, and universal reachability in the literal-extra graph below,
remain **OPEN**. The solved larger union graph is not identified with
either of them.

## 1. Four objects that must remain separate

Let the ordinary native map be
\(C(n)=n/2\) for even \(n\), and \(C(n)=3n+1\) for odd \(n\).

| Object | Allowed native moves | Inherited status |
|---|---|---|
| Ordinary forward Collatz | \(n\to C(n)\) | All-positive convergence open |
| \(A\), all genuine reverses | Ordinary edges in both directions; predecessors \(2n\), and \((n-1)/3\) only for \(n\equiv4\pmod6\) | Connectivity to 1 is exactly ordinary convergence |
| \(L\), literal extra moves | Ordinary forward edges; even \(n\to2n\); odd \(n>1\), \(n\equiv1\pmod3\), \(n\to(n-1)/3\) | Universal reachability open |
| \(B\), the union | All \(A\) edges plus the extra odd inverse just specified | Explicit descent to 1 is proved |

These are recovered from
[endpoint-context Collatz models](../../endpoint_context_recursion_20261113/COLLATZ_MODELS.md),
the independently integrated
[literal-graph note](../../fill_tree_matrix_20261113/collatz/NOTE.md), and
[spanning-tree-fill's union-graph proof](../../spanning_tree_fill_20261113/collatz/PROOF.md).
The first source proves even the equality of the minimum possible peak
between ordinary forward reachability and model \(A\). The union \(B\)
has a strict macro-parent with peak at most \(12n+4\); the literal graph
\(L\) does not inherit that proof. For example \(40\to13\) and
\(5\to10\) used by union macros are not legal \(L\) edges.

The extra odd inverse at \(7\to2\) has an **even** quotient. It is not a
genuine Collatz inverse, and it is not the minus-map step \(3n-1\).
Separately the tournament clock in
[CLOCK_BOUNDARY.md](CLOCK_BOUNDARY.md) implements some steps of
\((3H-1)/2\) on positive odd \(H\), equivalently an ordinary plus Collatz
odd step followed by one halving on the **negative** integer \(-H\).
This is a separate arithmetic convention,
not an edge of positive ordinary Collatz merely because it has a similar
formula.

The inheritance board is: genuine reversibility / literal guards /
affine translation carry / positive-source domain / clone resources.
The older
[commutator clock](../../carmichael_collatz_carry_20261019/clock/COMMUTATOR_CLOCK_MIXING.md)
shows why finite modular affine words need not have legal integer lifts.
[THM-4503, the two-chain parity-poset bridge](../../../../01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md)
retains source residues when word order changes.
[THM-4581, affine coalescence](../../../../01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md)
already supplies a Haar/density-one result, explicitly not a theorem about
every fixed positive source. None of these is re-promoted to pointwise
termination here.

## 2. Odd returns in the literal graph

Between two consecutive odd vertices, an \(L\)-path first uses one of
the two odd actions

\[
U(x)=3x+1,\qquad V(x)=(x-1)/3,                           \tag{2.1}
\]

where \(V\) requires \(x>1\), \(x\equiv1\pmod6\). Both outputs are
positive even integers. Subsequent even doublings and halvings remain
on one dyadic chain. If the path reaches an odd vertex again, its endpoint
is the odd part of that output; the extra up-and-down even excursions
can be deleted. Thus the reduced odd-return word consists of actions
\(U\) or \(V\), each followed by \(k\ge1\) halvings. In an actual
integer return, \(k\) is exactly the output's 2-adic valuation.

The next lemma applies to a larger real domain, so its proof does not
hide a parity assumption. Write \(E(x)=x/2\).

**Lemma 2.1 (no positive balanced return).** Consider a finite cyclic
word of \(U,V\) with one or more halvings between consecutive actions.
If its number \(b\) of \(V\) actions is at least its number \(a\) of
\(U\) actions, it has no periodic orbit with every intermediate value
strictly positive.

**Proof.** Use the class of affine contractions

\[
\mathscr D=\{x\mapsto\alpha x-\beta:
                   0<\alpha\le1,\ \beta\ge0\}.
\]

It is closed under composition and contains every halving. If \(D\) is
in this class, a cyclically adjacent \(U\)-then-\(V\) pair, including
the intervening \(D\), reduces exactly to

\[
V\circ D\circ U(x)
 =\alpha x-\frac{1-\alpha+\beta}{3}\in\mathscr D.         \tag{2.2}
\]

As long as both letters remain on the circle, an adjacent \(UV\) pair
exists. Cancel pairs using (2.2), retaining all the intervening affine
maps. Since \(b\ge a\), eventually no \(U\) remains. The remaining
\(V\) maps are also in \(\mathscr D\). The total resulting map has slope
strictly smaller than one, because the original word contains a positive
number of halvings. Its intercept is nonpositive, so it has no strictly
positive fixed point. A putative positive periodic orbit can be rotated
to any stage of the cyclic cancellation, and each replacement retains its
positive endpoints. This contradicts the final contraction. ∎

The intermediate-positivity hypothesis matters. For instance
\(U/2\circ V/2\) is the affine map \(x\mapsto(x+1)/4\), with positive
fixed point \(1/3\). Its intervening value after \(V/2\) is \(-1/9\).
An endpoint-only positive fixed point therefore does not constitute a
positive return.

## 3. A stronger finite-path descent criterion

**Theorem 3.1 (inverse debt paid implies strict descent).** Take any
nonempty word in \(U,V\), each action followed by one or more halvings,
with \(b\ge a\). If its
successive real values are positive, its final value is strictly smaller
than its initial value. In particular this holds for every legally
guarded positive integer \(L\)-path of that form. Its endpoint may be
even; completing the final odd return only decreases it further.

**Proof.** The whole affine map is \(F(x)=\lambda x+c\), where

\[
\lambda=3^{a-b}2^{-K}<1
\]

and \(K\ge1\) is the total number of halvings. Suppose
\(F(x_0)\ge x_0>0\). Its fixed point
\(x_*=c/(1-\lambda)\) satisfies \(x_*\ge x_0\). Every prefix of the word
is an affine map with positive slope, so all prefix values at \(x_*\)
are at least the corresponding positive values at \(x_0\). This supplies
the all-positive periodic orbit forbidden by Lemma 2.1. ∎

Thus the integer difference

\[
\#(\text{forward }3n+1)-\#(\text{extra odd inverse})
\]

is a useful sufficient descent certificate when it has returned to zero
or become negative. It is not a necessary condition for descent:
\(3\to10\to5\to16\to8\to4\to2\to1\) descends using two \(U\)'s and
no \(V\). Nor is there a proof that the budget is eventually paid on
every source or on a specified greedy path.

**Cycle corollary.** Every \(L\)-cycle containing an odd vertex uses
strictly more \(U\) actions than \(V\) actions, after deleting even-chain
backtracking. Purely even cycles such as \(2\to4\to2\) are outside this
statement; they exist because even doubling was explicitly allowed.
The ordinary root cycle \(1\to4\to2\to1\) has one \(U\), no \(V\), and
is consistent with the corollary. This does not classify the remaining
positive cycles or establish termination.

## 4. An explicit guarded source family with paid debt

For \(r\ge1,q\ge1\), put

\[
R_r=\frac{4^r-1}{3},\qquad n=2\,4^r q+R_r.
\]

There is a literal \(L\)-path

\[
\boxed{n\xrightarrow{U}4^r(6q+1)
 \xrightarrow{2r\text{ halvings}}6q+1
 \xrightarrow{V}2q\xrightarrow{E}q.}                   \tag{4.1}
\]

Every displayed guard holds: \(n\) is positive odd, \(6q+1>1\) is odd
and is one modulo three, and the intermediate halving inputs are even.
The path has \(2r+3\) native edges and peak \(3n+1\). It has one \(U\)
and one \(V\), so (4.1) is a concrete positive instance of Theorem 3.1.
Its endpoint is specified in advance, not merely some unspecified smaller
integer. For \(r=1\), it reads

\[
8q+1\to24q+4\to12q+2\to6q+1\to2q\to q.
\]

At \(q=0\), the forward path from \(R_r\) reaches 1; it must stop there,
because the extra move \(1\to0\) is not permitted. This is the sharp
positive-domain boundary of the formula.

The family concerns \(L\), not \(A\): its \(V\) step is precisely the
nonclassical even-quotient inverse. It gives a controlled source lift
and descent whenever that arithmetic shape is present, without claiming
that all positive sources have such a shape or reach one.

## 5. The tournament-side resource is now explicit

[CLOCK_BOUNDARY.md](CLOCK_BOUNDARY.md) proves the complete zero-site law

\[
H(\Gamma[M])=60P_M(3)-17H(M).
\]

For the inherited \(\Pi\) chart, consuming the final dyadic clone misses
the desired minus-Collatz transfer by exactly \(18H(M)\). A fresh native
root-child birth restores the clock, with exact count law
\(H\mapsto2H+35H(M)\) at the boundary; the full reusable circuit is
\(H\mapsto3H-2H(M)\). It therefore pays a physical resource cost rather
than providing a free continuation of the minus map.

The common lesson is precise. An affine word needs its admissible positive
source in the arithmetic graph; a scalar tournament identity needs its
marked nonempty module in the construction graph. In each case the lost
sidecar can be restored and the cost computed. Neither restoration is a
proof of arbitrary-source or arbitrary-target coverage.

## 6. Controls

Run from the worktree root:

```bash
python3 -B 05-knowledge/results/target_collatz_transitions_20261113/collatz/check.py
python3 -B 05-knowledge/results/target_collatz_transitions_20261113/collatz/clock_boundary.py
```

The first checker uses exact rational arithmetic on all \(U,V\) words of
length at most five, with each intervening halving count in \(1,2,3\),
and tests the positive-domain debt statement. It checks literal guards,
source formulas, and the named forbidden transitions independently.
These finite controls supplement the all-word proof; they do not infer
that an infinite orbit terminates. The second checker certifies the
fixed-skeleton clock boundary and its arbitrary-module law separately.
