# A 425 source rebase and its exact infinite-family boundary

**PROVED construction and formulas; FINITE-EXACT complete ready-value
inverse at425; Chan–Pak remains OPEN.**

## 1. An actual source return after the uniform reset

Use the append ordering of V(k,1) in [COMPOSITION](COMPOSITION.md):
prefix vertices0..k−1, B=k, A=k+1, ports a=k+2,b=k+3,
root r=k+4, first sink x=k+5, second sink y=k+6, apex z=k+7.
Delete b and x, preserving every remaining arc, and call the resulting
graph W_k. Append a vertex v beaten only by z, and call the result R_k.
This retains k+6 of the k+8 vertices of V(k,1), including its old root
and one old port. It sacrifices the other old port.

For every k>=1,

    H(W_k)=11·2^k+9,
    H(R_k)=24·3^k+53·2^k−3,
    (s_v,q_v,F_v−H(R_k),delta_v)=(11·2^k+8,11·2^k+9,1,0).

In particular the uniform141 reset now has the continuation

    141 → 59 → 59 → 59 → 425 → 53 → 425 → 799 → 1277.

The second425 is genuinely1-ready. The last two steps are the usual
two successive births, first beaten only by v and then beaten only by
the previous newborn; their final count is3·425+2. The last target is
not asserted to be ready again. The intermediate sink steps are weak,
as already declared in the uniform construction.

The 425 rebase has ten old vertices, retains eight, and has nine final
vertices. Its retained old labels are(0,1,2,3,4,6,8,9); its newborn
outmask in this retained order is127. All retained arcs are checked.
This is a source-compatible continuation, not merely a second scalar
realization of425.

## 2. Why it is ready, and why it only preserves the count at k2

In W_k the old apex z has exactly one Hamilton path starting there:

    z,B,A,a,r,prefix0,...,prefix(k−1),y.

Indeed z can only go to B. Vertex y only beats z, already used, and
therefore must be last. Each of B,A,a is then forced to choose the next
displayed vertex; the remaining transitive prefix has only its increasing
Hamilton path. The displayed path closes by y→z, proving W_k strong.
Inserting v between z and B gives a spanning directed cycle in R_k.
The indegree-one birth identity says the new root's gap is the old
root's start count, and gives q_v=H(W_k), s_v=H(W_k)−1.

For independent all-k count proofs, treat the first prefix vertex at k1
as a placeholder and inflate it to TT_k. The positive ordered-run
coefficients (A1,A2,A3,...) are

    W: (31,11,0,...),       R: (175,125,24,0,...).

The exact placeholder recurrence and independent literal DP certify
these complete finite lists. Substitution into
H=sum_j A_j·j!·S(k,j) gives the displayed formulas, including the
positive3^k channel in R. This channel is not inferred from a fit.

The count difference from V(k,1) is

    f(k)=24·3^k−49·2^k−20.

Here f(1)=−46, f(2)=0, and
f(k+1)−2f(k)=24·3^k+20>0. Thus the rebase preserves the count
exactly at k2; it lowers it at k1 and raises it at every k>=3.
The operation is uniform and always creates a1-ready source, but it
does not give a uniform same-count return. Its target addresses have
a new positive third-run coordinate.

## 3. Complete ready-image information at425

The inherited complete inverse in [OBSTRUCTION, O2–O3]
(../../ready_return_pentagon_20261103/OBSTRUCTION.md) reduces every
strong1- or2-ready root to a classified one-/two-start parent, or one
of five exceptional2-ready tails. A newborn of target count425 has
parent count at most(425+b)/2, so Busch's published bound f(11)=225
forces parent order<=10. The reproduced complete parent grammar has
counts, by orders3..10,

    one-start: 1,3,8,21,55,144,377,987;
    two-start: 0,1,5,18,60,191,589,1775.

The inverse finds425 in the1-ready image and excludes it from the
2-ready image. Enumerating all one-start parents, rather than stopping
at the first value witness, finds exactly ONE rooted strong1-ready425
class. It is the nine-vertex returned source above. None of the ten
nine-vertex deletion cards of V(2,1) is isomorphic to that graph or its
converse. Thus the exhibited eight-retained return is retention-optimal
among replacements of this V(2,1) into any strong1-/2-ready425 source,
including readiness after reversing all arcs, with every retained arc
unchanged. This optimum permits any number of auxiliary insertions and
deletions: the final strong ready graph is necessarily the specified
nine-vertex graph. Weak padded targets are outside this optimum's scope.

Consequently the usual2-ready+4 recipe cannot begin
from any strongH425 source, at any order. This is an interface
obstruction, not a claim that1279 is unattainable or that other+4
operations are impossible. The explicit1-ready witness above is
verified independently of that inverse.

## 4. The actual source return improves every future cost

The incoming [exact source-order theorem]
(../../h13_paired_descent_20261111/dynamics/PROOF.md) uses
(H,E,U,d), where E counts paths ending at the root, U is the ordered
two-path-cover sum with the first path starting there, and d is the
root-deletion count. The actual old root and new root above have states

    old V425: (425,4,2542,41),
    new R425: (425,15,1491,53).

Their positive-frame difference in q=(H,U,U−d,H+U−E) is

    q_old−q_new=(0,1051,1063,1062).

Consequently every common future O/Z/J program has old count at least
new count. Pure O preserves425; the first Z or J after any O-prefix
creates a strict decrease, and subsequent operations preserve it. The
theorem also permits arbitrary fixed right-root amalgams. This is a
cost comparison for matched future programs, not a target-surjectivity
claim: a cheaper continuation may have a different desired target.

There is a stronger fixed-context statement. In the seven-coordinate
order(H,E,S,U,V,d,K) of incoming [ADDITIVITY section5]
(../../h13_paired_descent_20261111/modular_cycles/ADDITIVITY.md), the
old and new states are

    (425,4,41,2542,449,41,2848),
    (425,15,52,1491,426,53,1544).

Their H,U,V,K differences are nonnegative and their E,S,d differences
are nonpositive, exactly that theorem's sufficient cone. Hence inserting
the rebase into any common fixed rooted-amalgamation context on either
or both sides also weakly reduces the count, with the stated sidecar
inequalities preserved. A subsequent common O/Z/J program preserves
the resulting q comparison. No necessity of the seven-coordinate cone
or arbitrary intervening left contexts after each root operation is
claimed.

Reproduce all construction checks with `python -B ready_rebase.py`;
add `--inverse` to rerun the complete ready-value inversion through425.
The first mode checks k1..6, exact placeholder coefficients, actual
persistent arcs, strongness, every claimed profile, and two-birth
continuations. The second mode writes a separate complete-inverse receipt.
