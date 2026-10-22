# Independent audit: forest asymptotics and merged-metagraph masses

Date: 2026-10-22. Verdict: **PASS**.

This audit accepts the basic ensemble identity already independently audited
in `session31-mirror/HAMILTONIAN_FOREST_AUDIT.md` and does not duplicate its
tridiagonal proof. It audits the later asymptotic and merged-mass addenda of
`session31-root/HAMILTONIAN_FOREST_SIZE_BIAS.md`. No novelty claim is made.

## 1. Third-order expansion, by an independent route

Write `R_n(t)=E[kappa F]/E[F]` and `v=t^2/(4(1+nt/2)^2)`.
The audited real polynomial identity gives

`R_n(t)=sum_{l=0}^{floor((n-1)/2)} binom(n-1-l,l) v^l`.

For fixed `t>0`, `v=n^-2(1+2/(nt))^-2`. The first four terms are

`1 + (n-2)v + ((n-3)(n-4)/2)v^2 + ((n-4)(n-5)(n-6)/6)v^3`.

For every allowed `l`, `binom(n-1-l,l)<=n^l/l!`. Thus the omitted tail is
bounded by `(nv)^4 exp(nv)/24=O(n^-4)`. Expanding just these four rational
terms proves

`R_n(t)=1+1/n-(3/2+4/t)/n^2+(4/t+12/t^2-10/3)/n^3+O_t(n^-4)`.

This obtains the source's coefficient without using cosine power sums.
The original logarithmic derivation also passes: the half-spectrum sums
are `(n-2)/4` and `(3n-8)/16`, respectively, for `n>=3`.

As an extra exact regression, retaining `l=4` gives the next coefficient

`[n^-4] R_n(t) = 85/24 + 26/t - 4/t^2 - 32/t^3`.

The independent checker verifies the expansion symbolically and, at
`n=8,16,32,64` and `t=1/2,1,2,7/3`, checks the entire finite remainder with
rational arithmetic. It uses the explicit upper bound

`sum_{l>=4}(nv)^l/l! <= (nv)^4 / (24(1-nv/5))` when `nv<1`.

The restriction to fixed positive `t` is material. At `t=0`, `R_n=1`
exactly. For `t=tau/n` with fixed `tau>0`, the same matching expression gives
`n(R_n-1) -> tau^2/(tau+2)^2`, not one. Thus the source is correct not to
claim uniformity down to zero.

## 2. Orbit accounting and the mass identities

The group `Aut(T)` acts freely on ordered directed Hamiltonian paths: an
automorphism fixing a spanning ordered list fixes every vertex. Hence
`h(T)=H(T)/|Aut(T)|` is an integer. Relabeling any ordered path to the fixed
path `1->...->n` identifies its automorphism orbit with exactly one fixed-path
presentation of the unlabelled class. There are `2^d` presentations in total,
where `d=(n-1)(n-2)/2`.

Reversal of an ordered path is a bijection between paths of `T` and `T*`;
also `Aut(T)=Aut(T*)` as permutation groups. Therefore paired converse
classes have equal `h`. With weight `h` on a self-converse node and `2h` on
a distinct pair, the average of the two forest polynomials is precisely
what turns a merged sum back into the unmerged class sum.

More generally, every isomorphism-invariant statistic `Phi` satisfies

`sum_classes H(T) Phi(T)/|Aut(T)| = (1/n!) sum_labelled H(T) Phi(T)`.

Since `kappa=2^(n-1)H/n!` and the labelled universe has `2^(n(n-1)/2)`
elements, this equals `2^d E[kappa Phi]`. Substituting the audited forest
identity proves the merged polynomial mass identity. Taking
`Phi=||(A-A^T)1||^2`, whose unbiased expectation is `n(n-1)` and whose
covariance with kappa is `-2(n-2)`, proves the energy mass

`2^d(n^2-3n+4)`.

No additional incidence or independence assumption is hidden in these
steps. The graph inventory must include all isomorphism classes on the
fixed order `n`, exactly once before converse merging, with the stated
weights. Arbitrarily filtered, duplicated or differently weighted
inventories do not satisfy the same normalization automatically.

## 3. Converse averaging is not optional

A directed triangle plus a universal source has

`Fout(t)=1+6t+12t^2+9t^3`.

Its converse, a directed triangle plus a universal sink, has

`Fout(t)=1+6t+12t^2+7t^3`.

Both have `H=3`, automorphism order three and `h=1`; the pair contributes
weight two times `1+6t+12t^2+8t^3`. Replacing the average by either
representative breaks the cubic mass. A self-converse class has no such
ambiguity because isomorphism invariance already equates the two values.

## 4. Independent finite-exact test of the addenda

`audit_bridges.py` begins with fixed-path orientations, builds each full
permutation orbit, and checks that the orbits cover every labelled
tournament. It obtains automorphism orders from orbit sizes, independently
counts Hamiltonian paths by permutations, and counts rooted forests by
acyclic parent functions. No source-worker code is imported; neither a
determinant nor subset-DP path counts are used.

For `n=2,3,4,5,6`, the respective numbers are:

- Unlabelled classes: `1,2,4,12,56`.
- Converse-merged nodes: `1,2,3,10,34`.
- Fixed-path presentations / total node mass: `1,2,8,64,1024`.
- Signed-energy masses: `2,8,64,896,22528`.

Every class satisfies `presentation_count=H/|Aut|`. Every merged
polynomial mass and energy identity holds exactly. The complete class and
merged inventories are saved in `audit_bridges.json`, including rational
averages without rounding.

At `n=6`, for example, the merged forest-mass polynomial has ascending
coefficients `[1024,15360,93184,285696,442560,277056]`.

## 5. Strong information-loss witness independently certified

The root researcher supplied the two order-five masks 40 and 41 in lexicographic
unordered-pair encoding, with a zero bit orienting the smaller vertex toward
the larger and a one bit reversing it. The audit regenerates the tournaments
from those integers and independently verifies

`Fout=Fin=1+10t+39t^2+70t^3+49t^4`

for both, but obtains respectively

- Mask 40: `H=15`, `|Aut|=3`, `h=5`.
- Mask 41: `H=13`, `|Aut|=1`, `h=13`.

Both classes are self-converse, and both have signed score energy eight.
Thus even the ordered pair of complete in- and out-forest polynomials,
together with score energy and self-converseness, does not determine the
Hamiltonian count or the node mass. This is stronger than a collision of one
special value. It cannot be repaired merely by adding the converse forest
polynomial. A separate Hamiltonian/automorphism or incidence sidecar remains
necessary.

The audit's mask convention inside the enumeration is the opposite of the
root collision convention; it explicitly complements all bits before testing
the supplied examples. This convention change is recorded in the JSON.

## 6. Reproduction and limits

Run from the workspace root:

`python3 session31-bridge-audit/audit_bridges.py`

The original category checker can be replayed without touching its output:

`python3 session31-categories/check_categorical_parameters.py > session31-bridge-audit/category_original_replay.json`

These identities constrain a metagraph inventory but do not certify its
realizability, recover its colored incidence, or prove a Hamiltonian-spectrum
gap theorem. The directed branching forests here are not the signed linear
forest coefficients of the earlier loop expansion.
