# Coarse-mean certificates with unresolved Collatz terminals

2026-10-15. **PROVED:** a source-preserving signed-terminal representation
for every coarse current mean; an upper certificate improving under every
terminal extension; and a sufficient mean/variance energy criterion that
does not require grounding every source. **FINITE-EXACT:** successful
certificates on six windows through X=2^20, including 2,962 sources still
unresolved by the successful controller at X=2^20. **OPEN:** uniform
acceptance over all large dyadic windows.

This is the next step after
[the sharp refinement envelope](collatz_clock_mesh_20261015a.md). It uses
the exact native odd and even edges, unlike the doubling-only hostile in
[the duality/glue note](clock_duality_glue_20261015a.md). The general
principle of retaining signed boundary forcing was already present in
[the Poisson patch calculus](collatz_poisson_patch_20261005.md).
The new specialization is the actual coarse-mean weight row, whose norm
can certify energy even with nonzero unresolved terminal mass.

Other current routes remain distinct: [Mersenne guards](mersenne_clocked_completion_20261014c.md)
prove the giant minimum has no smaller odd predecessor through64 steps;
[signed fibres and first-hit floors](prime_power_orbit_boundary_20261014d.md)
retain the source and sign. Those ROOT-entry obligations are not discharged
by an aggregate energy certificate.

## 1. Keep the mean as a weighted source observable

Let T(n)=n/2 for even n and (3n+1)/2 for odd n, and let t be the
ROOT-supported first-hit parity label, zero on sources not reaching1.
Then t(1)=1 and t(n)=-t(T(n)) for every positive integer n.
The identity at1 is consistent because the terminal shortcut cycle
1->2->1 has even length. No convergence assumption is made.

For a physical interval I=[a,b) of length L=b-a, define the finite row

    A_I(n) = -L w(n)                           (1<=n<a),
             -(b-n)w(n)+W(n)                   (a<=n<b),
             0                                (otherwise),

where w(n)=3[n=2 mod3]-1 and W(n)=0,-1,+1 by residue. Exactly,

    2 M_I := 2 sum_(n in I) J(n) = sum_n A_I(n)t(n).           (1)

This follows by summing 2J(n)=W(n)t(n)-H(n) over I. The weights
on n<a are indispensable: they carry the incoming current. Replacing
them with zero resets the source and changes the quantity being proved.

Choose any finite, individually authenticated prefix for each source n
in the support, stopping at endpoint e(n) after d(n) shortcut steps.
If a prefix reaches1, stop it there. Push the row by

    B_I(e)=sum_(n:e(n)=e) A_I(n)(-1)^d(n).

The endpoints may coincide, the prefix depths may differ, and endpoints
may be much larger than b. Only exact native edges are permitted.
Writing K_I=B_I(1) and R_I=sum_(e!=1)|B_I(e)| gives

    2M_I=K_I+sum_(e!=1) B_I(e)t(e),
    |2M_I-K_I|<=R_I,
    |M_I| <= (|K_I|+R_I)/2.                                 (2)

This representation exists for every finite source row, including the
zero-depth controller. Existence alone proves no small bound.

## 2. Changing terminal controllers cannot worsen this bound

**Theorem T1 (PROVED).** Extend any unresolved terminal e by one exact
step, transporting its whole coefficient to T(e) with a minus sign and
combining coincident terminals. Keep1 absorbed into K_I. Then

    Bnd_I:=|K_I|+R_I

cannot increase. Therefore any sequence of finite terminal extensions,
chosen independently for each row or shared across rows, retains (2)
and weakly improves its upper bound.

**Proof.** Moving to a previously unused terminal preserves absolute
mass. At a collision, |x+y|<=|x|+|y|. At ROOT, replace K by K+x
and remove residual x; again |K+x|<=|K|+|x|. These are all cases. ∎

The theorem is deterministic. It neither assumes random terminal signs
nor estimates residual mass by the number of unresolved sources. Large
opposing weights can cancel, and many unresolved sources can share one
terminal. Conversely a lone unconstrained terminal can simply move with
unchanged norm, so there is no uniform strict contraction for arbitrary
rows. The clock lengths L and addresses a,b stay fixed while d(n) changes.

Combining T1 with M1--M2 of the mesh note gives, for any physical partition,

    S_X <= Z_term := sum_I Bnd_I^2/(4L_I)
                        + sum_I V_*(a_I,L_I),               (3)
    E_X <= (2 sqrt(Z_term)+sqrt(X))^2.                       (4)

Thus no unknown terminal sign needs to be evaluated. A finite certificate
consists of exact native prefixes or an equivalent checked transport,
the physical weight rows, terminal coefficient sums, and the rational
test `A=2X^2-4Z_term-X>=0`, `A^2>=16X Z_term`.

### Physical splitting is a different refinement

For exact means, M2 proves monotone upper bounds under every interval
split. Independent rowwise absolute-value bounds in (3) can discard
cross-row correlations, so they are not automatically monotone under
physical splitting. Retain the parent upper bound and take the smaller
of valid certificates, or retain the joint terminal rows and maximize
over the same terminal cube. For the latter, the pointwise M2 inequality
survives the common supremum; forgetting shared terminals would break
that argument. Terminal extension at a fixed partition is monotone
without this additional optimization.

## 3. Exact successful certificates before full grounding

Each controller below follows every contributing source n<2X for at most
D shortcut steps, stopping early at1. The native prefix is reconstructed
using arbitrary-precision integers. No parity label for an unresolved
endpoint is used in computing the certificate.

| k, X=2^k | mesh b | depth D | unresolved sources | nonzero unresolved terminals | Z_term/X^2 |
|---|---:|---:|---:|---:|---:|
| 10 | 32 | 96 | 49 | 25 | 0.296254426 |
| 12 | 64 | 128 | 25 | 17 | 0.132839782 |
| 14 | 128 | 160 | 77 | 20 | 0.326241773 |
| 16 | 256 | 160 | 1076 | 115 | 0.311158444 |
| 18 | 512 | 192 | 1789 | 167 | 0.323872687 |
| 20 | 1024 | 224 | 2962 | 229 | 0.399489714 |

Every row passes the exact square test and therefore certifies
E_X<=2X^2. These are the first successful depths in the declared tested
lists, not proved minimal depths. For instance k20 at depth192 fails
this sufficient test, with Z_term/X^2=1.137066579; extending the same
source rows to224 makes the certificate pass. No source is replaced.

The output records both the known signed contribution K_I and residual
l1 weight R_I for every physical block. An independent audit grounds
sources by first descent, reconstructs actual J and block sums, checks
their inclusion in every terminal enclosure, compares with the separate
native Haar census, and checks controller monotonicity block by block.
An additional formal audit assigns {-1,0,1} values to terminal variables
and compares the weight row with direct summation of the resulting
source current. These assignments are algebraic controls, not claimed
Collatz labelings.

```sh
python3 04-computation/experiments/collatz_terminal_mesh_20261015a.py 05-knowledge/results/collatz_terminal_mesh_20261015a.json
python3 04-computation/experiments/collatz_terminal_mesh_20261015a.py 05-knowledge/results/collatz_terminal_mesh_large_20261015a.json --large
python3 04-computation/experiments/collatz_terminal_mesh_audit_20261015a.py
```

Data: [small-window controllers](collatz_terminal_mesh_20261015a.json),
[larger controllers](collatz_terminal_mesh_large_20261015a.json),
[independent audit](collatz_terminal_mesh_audit_20261015a.json).

## 4. A sharper open target, and what would count as progress

For all dyadic X>=2^34, construct finite native terminal controllers on
the square-root mesh such that

    sum_I (|K_I|+R_I)^2/(4L_I) <= (4/9) X^2.                (5)

Together with the proved fine-scale allowance this suffices for the
amplitude energy premise. It is weaker than requiring every contributing
source to reach1, as the finite witnesses above demonstrate. It is
stronger than merely providing a complete grammar for transporting rows.

Three concrete helper questions now have an exact predicate to test:

1. Can a variable-depth controller bound the weighted residual R_I in
   aggregate, rather than the unweighted number of escaping sources?
2. Can the joint terminal-row Gram matrix exploit correlations enough to
   give a stable cone under the native odd/even transfer? Marginal l1
   cancellation is only a sufficient relaxation.
3. Can a source-addressed block/glue grammar certify those cone inequalities
   across scales, retaining its seam cost and physical length? The golden
   conditional-comparison hostile and the Q projective-scale example are
   mandatory controls for any proposed closure.

There is no evidence here of a uniform strict per-step contraction, and
no proved upper bound on the necessary controller depth at all scales.
The six successful windows are finite evidence for this route, not a
proof of (5) or of universal ROOT entry.

## 5. Joint terminal Gram control improves the mean bound

The previous rowwise relaxation forgets that the same unresolved terminal
has the same label in different blocks. Retain this information. Include
ROOT as terminal column0, with z_0=1, and let v_I be the complete signed
terminal row B_I, including its known ROOT coefficient. Define

    G=sum_I v_I v_I^T/(4L_I),
    C(P)=z^T G z,   z_0=1, |z_e|<=1.                       (6)

This is the exact source coarse energy, for the actual terminal labels.
The finite PSD matrix is an authenticated object, not a postulated
two-coordinate covariance closure.

**Theorem T2 (PROVED).**

    C(P) <= ||G||_(entrywise 1)
         <= sum_I (|K_I|+R_I)^2/(4L_I).                    (7)

Under every common native terminal extension, including ROOT absorption,
the entrywise norm in (7) cannot increase. Keep the physical partition
fixed for this statement.

**Proof.** The first inequality follows by taking absolute values in the
quadratic form, and the second by the triangle inequality before summing
the outer products. A common terminal extension is a signed identification
matrix A with exactly one entry +1 or -1 in each old-terminal row; the
ROOT row maps to ROOT with sign+1. The new Gram matrix is A^T G A.
Every old matrix entry contributes, with a sign, to precisely one new
entry. Therefore `sum|A^T G A|<=sum|G|`. ∎

This explains why retaining the joint matrix is useful: opposite
contributions can cancel across blocks before absolute values are taken.
It is the concrete source-preserving counterpart of asking the Q--Pell
analogy to retain correlations instead of only one marginal statistic.
No Q sampling hypothesis is used in T2.

### A small signed graph sometimes solves the entire terminal cube

Put an undirected signed edge e--f whenever G_ef!=0, with sign sgn(G_ef).
Do not orient this symmetric observable into a cosmetic tournament.
If there are signs sigma_e with sigma_0=1 and
`sigma_e sigma_f=sgn(G_ef)` on every edge, the graph is balanced and
the upper bound (7) is attained exactly at z=sigma. Thus a sign witness
solves the full relaxed cube maximum, however many terminals there are.
Conversely attaining the entrywise bound requires this condition on all
nonzero off-diagonal entries. Zero entries impose no condition.

For an unbalanced triangle with edges of weights |G_ef|, every corner
assignment leaves at least one edge unsatisfied. Relative to the
entrywise bound it loses at least four times the smallest edge weight.
Summing this deduction over edge-disjoint negative triangles is valid.
Since G is PSD, its maximum on the terminal cube is attained at a corner
(successively maximize each convex coordinate function). This yields
a rigorous optional improvement, not an assumed random-sign model.

The fixed-source gains are already strict:

| k | depth | unresolved sources | rowwise S upper/X^2 | joint S upper/X^2 |
|---|---:|---:|---:|---:|
| 10 | 89 | 122 | 0.485919505 | 0.419628024 |
| 12 | 110 | 266 | 0.509313870 | 0.458170826 |
| 14 | 134 | 359 | 0.540767211 | 0.480941281 |

In all three cases the rowwise sufficient test fails while the joint
test certifies E_X<=2X^2. At k12, depth112, the signed graph on44
columns is balanced, so its recorded sign witness solves the relaxed
cube maximum over43 free labels without enumerating2^43 assignments.
At depth128 the graph is unbalanced; all2^17 assignments were enumerated
both by Gray-code updates and by independent direct quadratic evaluation,
as a control on the triangle bound. The exact maximum numerator is
348638364, below the entrywise bound349335076 and the triangle
bound349137956, all with denominator256 for the coarse energy.

The script checks each actual common extension against separately
recomputed source prefixes at depth D+1 and verifies matrix-norm
monotonicity. It saves the full integer Gram matrix, denominator,
terminal addresses, sign witnesses, and negative-triangle certificates.
The triangle improvement need not itself be monotone under a greedy
choice of triangles; retain the best previous valid bound when using it.

```sh
python3 04-computation/experiments/collatz_terminal_gram_20261015a.py
```

[Joint certificate data](collatz_terminal_gram_20261015a.json).
The strengthened open target is now
`||G_X||_(entrywise 1)<=(4/9)X^2` after suitable native common
terminal extensions on every large dyadic window. This is a sufficient
condition, weaker than the rowwise condition (5); a direct bound on
the joint cube maximum could be weaker still. None is proved uniformly.
