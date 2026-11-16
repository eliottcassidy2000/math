# A quantitative phase cutoff and two continuations at height three

This is a follow-up to the [real-boundary packet](../real_boundary_signatures/README.md).
It does not prove Chan–Pak or exclude an unrestricted height-three continuation.

## A finite phase-crossing theorem

Let a target ray be A·2^s+a and a physical reader ray F·2^b+D,
with A,F positive integers and s,b nonnegative. Write the positive terms
as u·2^α and v·2^β with u,v odd. If d=D−a≠0, put w=v₂(|d|).
When α=β, the equation (u−v)·2^α=d has at most one solution and α≤w.
When α<β, α=w and d/2^w<u−v. When α>β, β=w and d/2^w>u−v.
The strict inequalities make the two unequal-valuation cases mutually
exclusive; in the available one, the remaining exponent is forced.
Thus distinct-intercept rays meet in at most two points. Also,
`A·2^s ≤ max(u,v+1)|D−a|`: either α≤w, or β=w and
u·2^α=v·2^w+d. Equal intercepts instead require equal odd mantissas
and give a fixed shift of exponent.

### Application to the old middle target 756·2^s−31

By the [finite-core reduction](../real_boundary_signatures/COMPACTIFICATION.md),
normalized hosts at g≤100 and |H−2E+g|≤Cg are a finite set of cores
and padded reader rays, with |D+g|≤Cg. By the
[aligned certificate](../real_boundary_signatures/ALIGNED.txt) and
its all-order stripping argument, D=−31 and odd mantissa 189 are
impossible together in this deletion band. The crossing theorem
therefore proves eventual exclusion for every fixed C.

We can also give an effective cutoff without enumerating physical cores.
Let K*=floor((1+C/2)·100), J*=3+floor(log₂ K*), and B=2^J*·10!.
After reduction the deletion has at most J* singletons. Each nonsingleton
component of size n has path count at least ceil(5^((n−1)/3)) by Busch.
The product of those lower bounds is at most 100, so n≤9. For a component
of size n its ordered two-slot count (empty slots permitted) is at most
(n+1)n!=(n+1)!: choose the first slot's length in a permutation.
In an ordinal sum the two-slot counts multiply, and a singleton contributes 2.
The greatest product of (n+1)! over admissible nonsingleton size lists is
10!, by the finite integer recursion in check.py.
The core's F=E+O−g is bounded by this full two-slot count B: remove the
root from a weighted root-start path, retaining the selected path as the
first slot and a Hamilton path on its complement as the second.
This map is injective. The core's unpadded H is at most 100+2B by the
weighted even-path cut formula, and a root with no contacts has H=g.

Consequently any admitted height satisfies one of
```
756·2^s ≤ max(189,B+1)·(100(C+1)+31)    [padded]
756·2^s ≤ 131+2B                       [unpadded or sink].
```
**In particular, at C=1 all heights s≥31 are excluded.**
The practical implication is that in this normalized bounded-defect
route, either the deletion budget or the relative defect cap has to
be exceeded at all sufficiently high prescribed heights. The special
intercept exclusion is inherited from the exact aligned scan, not from
a statistical fit.

## A different signature: ordinary sinks

Appending an ordinary sink preserves H and g. For any weighted
root-start path in the original tournament, the new sink can be last in
its selected path or last in its complementary path. The two choices
are bijective, have opposite selected parity, and have the same weight
as the original selection. This includes empty complement. Thus the
new pair is (E+O,E+O); m≥1 sinks give
E'=O'=2^(m−1)(E+O). Two pre-padding profiles can be equalized by
ordinary sinks exactly when g and oddpart(E+O) agree. This is a
different reader from the near-sink append used above.

The background search supplied and check.py directly verifies these
two exact pairs at the unresolved height:
- 6015 and 6017 give (E,O,g)=(12217,12217,729).
- 6017 and 6019 give (E,O,g)=(16102,16102,1185).

Each witness has 10 vertices before padding. The two middle witnesses
cannot be identified by sink padding, since their deletion counts differ.
The first pair can be closed by finding 6019 at (g,oddpart(E+O))=(729,12217),
the second by finding 6015 at (1185,8051). I searched all root contact
masks over 115,350 distinct labeled order-nine deletions drawn from a
seeded random mixture of mutations around the four witnesses and
uniform orientations. This gives 79 keys containing at least two target
counts, none containing three. It is not an exhaustive census of order
nine or of feasible keys. Reproduce with
`g++ -O3 -std=c++17 search.cpp -o /tmp/phase_scan && /tmp/phase_scan 150000`.

## A parity invariant recovered from the archive

For every root, O−E equals the number of Hamilton paths ending at that
root, not merely a nonnegative signed count in a finite census.
In G=T−r, separate root-start paths by first successor:
O−E=H(G)−Σ_(r→u)(O_u(G)−E_u(G)).
Induction identifies the subtracted terms with the paths of G ending
at outneighbors. The remaining paths extend by appending r. The
singleton start supplies the base. So 0≤O−E≤g, and for a nonempty
deletion equality on the right means r loses to every vertex of the
terminal deletion SCC (each vertex there can end a spanning path of
that SCC, by opening a Hamilton cycle when nonsingleton).
Ordinary sinks force the left equality; near-sinks force the right.
Credit: [archived head parity](https://github.com/eliottcassidy2000/math/blob/7a324a5a6e44f5d5e55c4ffd8876e419735be992/04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/head_parity_check.py);
the [Klein continuation note](https://github.com/eliottcassidy2000/math/blob/1b66550eae98f33b5f35f6ef66c2aa639e216c99/05-knowledge/results/path_splice_rows_continuation_20261116/CONTINUATION.md)
supplies the ordinary-sink change of perspective.
