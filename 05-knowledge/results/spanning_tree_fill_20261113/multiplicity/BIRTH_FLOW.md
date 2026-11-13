# Exact one-node count flow needs the gap and cut response

Status: **PROVED** weighted flow identity and literal preparation; **FINITE-EXACT**
complete small-source replay. The endpoint/two-path quadratic kernel is
shared with and derived independently in [matrix/PROOF.md](../matrix/PROOF.md);
we do not claim that kernel as a separate new result.

## 1. Exact labelled flow from a complete source row

For a tournament T on n vertices, let a newborn w beat exactly the subset
encoded by x∈{0,1}ⁿ. Let s_i and e_i count full paths starting and ending at i,
and let W_ij count ordered covers by two nonempty paths, the first ending at i
and the second starting at j. Partitioning paths according to w's position
is the exact quadratic kernel

\[
Q_T(x)=H(T)+\sum_i(s_i-e_i)x_i+
\sum_{i,j}W_{ij}(1-x_i)x_j.
\tag{1}
\]

Define `B_T(y)=Σ_x y^{Q_T(x)}`. Then, with ordinary source classes and actual
labelled output multiplicity L,

\[
\boxed{L_{n+1}(k)=[y^k]\sum_{[T]:|T|=n}\frac{n!}{|\operatorname{Aut}(T)|}B_T(y).}
\tag{2}
\]

Every labelled child on[n+1] has a unique parent obtained by deleting label
n+1, and a unique cut x. This proves(2), including the normalization: there
is no extra factor n+1. Equivalently, this is an exact flow between old and
new H labels, with the old classes carrying their orbit-stabilizer weights.
The row sums and first moment are

\[
\sum_k L_{n+1}(k)=2^{\binom{n+1}{2}},\qquad
\sum_k kL_{n+1}(k)=(n+1)!2^{\binom n2}.
\]

For the first moment, each ordered spanning path fixes n arcs and leaves
`binom(n+1,2)-n=binom(n,2)` free. These mass identities do not prove that
individual target coefficients are nonzero.

Ordinary child counts cannot be obtained by treating all cuts or all labelled
children as one class each. Parent-cut orbits correspond to isomorphism
classes of a child with one distinguished vertex. A child can have several
vertex orbits and hence several parent-cut presentations. Its ordinary
contribution becomes1 only after this marked-child multiplicity is removed.

## 2. Minimal small-source obstruction to a coarse transition state

At order5, compare these actual sources:

- `T_end=C3→TT2`, with rows `(26,28,25,16,0)`;
- `T_mid=1→C3→1`, with rows `(30,20,24,18,0)`.

Both have H=3, automorphism order3, and the same ordered nontrivial core word
(C3). They even share the entire ordered-slot polynomial
`P_T(z)=z^5+2z^3`: ordinal multiplication gives
`P_C3(z)P_TT2(z)=(z^3+2z)z^2`. Here P_T(m) counts ordered m path bins,
with empty bins allowed. This unmarked polynomial still loses the selected-cut response. Their32 newborn cuts have the following complete output multiplicities:

| child H | end-padded source | centered source |
|---:|---:|---:|
|3|4|4|
|5|6|6|
|9|4|6|
|13|3|6|
|15|4|2|
|17|3|0|
|23|3|0|
|25|3|0|
|27|1|2|
|29|0|6|
|33|1|0|

Thus neither `(n,H)` nor `(n,H,core word,automorphism order)` determines the
next-source distribution. The missing coordinate here is the singleton gap
vector. Even target support differs: one allows17 and the other29.
This is the first such H-only per-source discrepancy through order5; lower
orders have only identical or converse-paired sources in each H-fibre.
Converse preserves the all-cut distribution by complementing the cut and
reversing every arc.

This is a failure of a **common per-source rule** on the coarse state. It
does not deny that an aggregate conditional transition table can be defined
from a complete labelled source ensemble. Computing that table still needs
the source response distribution omitted by the coarse state.

## 3. A literal same-count preparation changes the next exact target

Keep the vertex names of T_end. Move singleton3 ahead of its C3 by reversing
exactly its three arcs with vertices0,1,2. The result is

```
T_prepared = (18,20,17,23,0).
```

It is isomorphic to T_mid and still has H=3. All five physical vertices are
retained; this is an ordinal-component relocation, so exact H preservation
follows from the product law. It is a count-neutral simultaneous preparation,
not a claim that each of the three individual arc flips is neutral.

Now apply the **same named newborn cut x=9**, with w beating vertices0 and3:

```
T_end + w:      (26,60,57,16,32,9), H=17;
T_prepared + w: (18,52,49,23,32,9), H=29.
```

The physical newborn incidences agree in both cases. The difference in future
count is12. Reversing the preparation selects17 instead of29. This gives an
actual use for the gap sidecar: moving a neutral singleton component can
change the target achievable by the next fixed cut without changing the
current H or the nontrivial core inventory.

## 4. Reproduction and scope

```
python3 05-knowledge/results/spanning_tree_fill_20261113/multiplicity/birth_flow.py
```

The script uses all18 inherited ordinary classes at orders3,4,5 and all464
newborn cuts. It independently compares(1) with literal subset path counting,
then recovers every labelled count at order6 from(2), agreeing exactly with
the inherited full row. `birth_flow_receipt.json` retains the two complete
histograms and the converse source. These are small controls of the new data
loss and exact flow; no large census or all-target coverage claim is made.

The padding forest in [PROOF.md](PROOF.md) already retains the gap vector.
To upgrade it from count-preserving padding to arbitrary target-producing
births, retain the cut response (1) as well. A fair enumeration of this richer
object space is complete for its actual objects; proving that a prescribed
odd target appears remains the target-support obligation.
