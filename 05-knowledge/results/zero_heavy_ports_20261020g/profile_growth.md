# Exact quadratic port response, and a six-vertex obstruction to scalar growth

**2026-10-20. PROVED:** the port-tensor formula and its mean/submodularity
consequences. **FINITE-EXACT:** the minimal reinforced-profile collision
through order six. **OPEN:** a closed smaller state for unrestricted
unlabelled tournament growth.

## Inheritance and the operation being preserved

The source is a tournament T; the operation adds one vertex x with arbitrary
incidences. The target observable is the entire **unordered directed**
spanning path-cover profile, including H=c_1. This is not the undirected
colour-response profile from the complementary game.

The closest structural mechanism is the exact clone dichotomy in
[weighted-core inheritance, Propositions 1–3](../weighted_core_transport_20261020/inheritance.md).
The closest count mechanism is [THM-3121, path-cover walk-content substitution](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).
Uniform substitution consumes different data from arbitrary attachment.
[THM-3235, observer twins, section 6](../../../01-canon/theorems/THM-3235-transitive-tournament-blowup-saddle-and-hypergeometric-decimation.md)
already proves that equal complete profiles can split under uniform lifts.
The present hostile strengthens this distinction for **pair-free cores and
one-vertex growth**, while adding an exact operation-specific repair.

The corrected near miss is deleting x from a Hamiltonian path and assuming
that the remaining word is a path. The middle deletion can leave two paths
whose joining arc points backwards. These are the orphan paths already
retained in [THM-012b, insertion decomposition](../../../01-canon/theorems/THM-012b-insertion-decomposition.md).
The least-used relevant data are the actual two component endpoints.

## 1. The quadratic formula

Write z_u=1 for u->x, and z_u=0 for x->u. Let c_d count unordered directed
spanning covers of T by d nonempty paths. Set out-of-range c_d=0; T here
is nonempty. Define:

- e_d(u): covers by d paths having a component ending at u;
- s_d(u): covers by d paths having a component starting at u;
- B_d(u,v): covers by d+1 paths, together with an ordered pair of **distinct**
  components, the first ending at u and the second starting at v.

Then, for every d>=1,

\[
\boxed{
c_d(T+x)=c_{d-1}+d c_d+
 \sum_u(e_d(u)-s_d(u))z_u+
 \sum_{u\ne v}B_d(u,v)z_u(1-z_v).
}                                                       \tag{1}
\]

**Proof.** Delete x from a cover of T+x. If x was a singleton, the
remaining cover has d-1 components, giving c_{d-1}. If x was an endpoint
of a nonsingleton component, the remainder has d components: it can be
appended to an endpoint u with weight z_u, or prepended to a start v with
weight 1-z_v. Summing the latter constant term gives d c_d. If x was
internal, deletion splits its component into two distinct nonempty paths
in a (d+1)-cover, with weight z_u(1-z_v). These three cases are disjoint;
the reverse insertions are unique. No assertion about the arc u->v is
made or needed. QED.

Thus the one-step response has degree at most two, even though the parent
can have arbitrarily many vertices and arbitrarily complicated paths.
The small fact is that x has at most two neighbours in a path cover.
It limits the **number of ports consumed**, not the information needed to
compute their coefficients or to update them indefinitely.

For H specifically this is

\[
 H(T+x)=H(T)+\sum_u L_u z_u-
                  \sum_{u<v}C_{uv}z_uz_v,              \tag{2}
\]

where L_u=e_1(u)-s_1(u)+sum_v B_1(u,v) and
C_uv=B_1(u,v)+B_1(v,u)>=0. All profile coordinates are therefore
submodular set functions of the set of old vertices beating x: their
mixed second Boolean differences are -C_uv (with the appropriate d).
This is a weighted cut response with possible linear terms, not an
intrinsic orientation of a new tournament on the ports.

For independent fair attachment bits,

\[
\boxed{\mathbb E c_d(T+x)
 =c_{d-1}+d c_d+\frac{d(d+1)}4c_{d+1}.}                 \tag{3}
\]

Indeed sum_u e_d(u)=sum_u s_d(u)=d c_d, and
sum_uv B_d(u,v)=d(d+1)c_{d+1}. In particular
E H(T+x)=H(T)+c_2(T)/2. The full parent profile controls the
**mean over labelled attachments** exactly; it does not control their
distribution or their unlabelled child classes.

## 2. Smallest reinforced-profile hostile in the audited universe

Outgoing-row bitsets on vertices 0,...,5 give

```
A = (2,12,17,37,11,23)
B = (2,36,9,19,7,29).
```

Both are pair-free, with

```
(c_1,...,c_6) = (29,131,162,77,15,1)
sorted scores = (1,2,2,3,3,4)
|Aut| = 1
minimum feedback arcs = 3.
```

They are neither isomorphic nor converses. Hence even their complete
automorphism cycle indices agree. Their unit-weight versions are literal
zero-heavy cores at feedback budget 3 (deficiency D=7).

Over all 64 attachments of one vertex, however,

\[
\max H(A+x)=159,\qquad \max H(B+x)=147.                 \tag{4}
\]

Both maxima occur at z=(1,1,1,0,0,0). Their linear coefficient vectors
in (2) are respectively

```
(84,48,46,26,38,20), (86,36,44,36,38,22).
```

Their complete coefficient arrays and response histograms are in
[the exact output](profile_growth_check.out). Their means are both 189/2,
as required by (3). The joint response of full path profile **and actual
canonical-chain lengths** also differs. Forgetting the ports is the first
failed implication; adding feedback, scores, symmetry, or total order does
not repair it. No such nonconverse pair exists at a smaller order among
pair-free cores with this exact reinforced signature.

## 3. Clone information is retained separately

For a canonical chain of T, its binary attachment word splits it into
maximal constant runs. Construct the quotient on these runs and x,
retaining run lengths and every incidence; then take its maximal
transitive homogeneous chains and sum their lengths. This exactly gives
the canonical weighted core of T+x.

This is inherited, rather than a new growth theorem here. The stronger
clone dichotomy says that if x joins any old chain, the update is simply
insertion at a cut of that chain, with the old core unchanged. Otherwise
there is no recompression after the run split. Vertices in two distinct
old chains cannot merge: deleting x from such a new transitive homogeneous
module would contradict the old maximal chains.

The present checker independently confirms this rule on 4,054 attachments
of every class of orders 1–6. Equation (1) and this structural update have
different jobs. Their combination preserves the immediate numeric and
clone outputs, but still needs child port tensors, automorphism orbits,
and parent ownership to iterate as an unlabelled recurrence. No closure of
those data has been proved.

## 4. Proof repair in the older insertion note

The original mechanism paragraph of THM-1900 wrongly said that every
Hamiltonian path restricts to a Hamiltonian path after vertex deletion.
Already a directed triangle a->x->b->a has the path (a,x,b), whose deletion
leaves the invalid word (a,b). The neutrality theorem survives, with the
following proof retaining all paths.

Along each old HP, write 1 for vertices beaten by x and 0 for vertices
beating x. The number of allowed insertion positions is

    [first bit=1]+#(0->1 transitions)+[last bit=0]
      = 1+#(1->0 transitions).

Different old HPs give disjoint non-orphan extensions. If the beaten set P
is a suffix of the SCC order, every old HP has a single insertion cut;
the whole new tournament is an ordered join with x at that cut, so there
are also no orphans. Thus H is unchanged.

If P is not such a suffix, either a strong component is split by P, or
two unsplit SCCs appear in the wrong order. In the first case, a Hamilton
cycle in that SCC has a P->complement transition; remove a different
cycle edge and concatenate HPs of the other SCCs. In the second case,
concatenate arbitrary component HPs. Either construction gives an old HP
with a 1->0 transition, hence at least one additional **non-orphan**
extension. All orphan contributions are nonnegative. Therefore H strictly
increases. This uses the classical strong-tournament Hamilton-cycle
theorem (Camion), also explicitly inherited by
[THM-4526](../../../01-canon/theorems/THM-4526-shaved-tournaments-path-plus-chord-and-the-largest-unavoidable-graphs.md).

## Reproduction and limits

```
python 05-knowledge/results/zero_heavy_ports_20261020g/profile_growth_check.py
python -O 05-knowledge/results/zero_heavy_ports_20261020g/profile_growth_check.py
```

The universe is every isomorphism class through order six, generated by
one-vertex augmentation; the row counts are 1,1,2,4,12,56. The path-profile
DP is checked independently by cuts of every vertex permutation and by
explicit unordered covers. All 4,054 extensions check the full formula
and clone update. The two hostile cores additionally receive a full
permutation isomorphism/converse, automorphism and feedback audit, and all
128 child H values are counted by literal seven-vertex permutations.
Assertions are not used. The bounded computations support the stated
proofs and locate a counterexample; they do not settle A000568 growth or
all-order black-core connectivity.
