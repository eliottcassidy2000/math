# Shared alignment: integration of the concurrent Session42 packet

23 October 2026. Incoming commits `cf7c4fc327` and `cce2774f3c` were fetched
after this packet's main derivations. **PROVED scoped** statements below
are separated from finite checks and open transfers. No priority claim is
made for results independently obtained in both sessions.

## 1. Overlap and what it changes

The incoming [cyclic attachment paper](../session42_shared_lattices_square_colors_20261023/lattices/README.md)
independently proves the same two-congruence HNF description, one-sided
kernel formulas, composite alignment cokernel and actual cyclic-nine
subcarrier. Its Goursat formulation explains why each individual attachment
is functional after the correct kernel quotient; it does not make the
original port relation functional. Our any-degree CRT elimination adds
exact pairwise fill-in and hidden-witness counts. The full-group common
kernel factorization and nested-cell marking theorem are separate extensions.

The incoming [founder-color report](../session42_shared_lattices_square_colors_20261023/square-colors/RESULTS.md)
independently establishes the same sharp all-founder tail at 71. It adds
forest multiplicities and weighted joint partition functions. Our report
adds the first shared-marginal failure at 22 and the pinned three-site
proper-color repair at 35. These are distinct coloring predicates.

The incoming [correspondence paper](../session42_shared_lattices_square_colors_20261023/correspondence/CORRESPONDENCE_AND_OBSTRUCTIONS.md)
independently proves the only-empty/complete observer-independence boundary
and connectedness of the full incidence graph. Its incidence weights are
the same two-sided order counts used in our packet, after the declared
class merging. Our quantitative adjacent-swap jump and order-seven census
refine that shared picture. Distinct observer-distribution rows do not
imply injectivity on arbitrary mixtures of source classes.

The incoming [square-word reconstruction](../session42_shared_lattices_square_colors_20261023/SQUARE_WORD_RECONSTRUCTION.md)
retains the exact ordered edge sums and an affine starting coordinate. It
proves at most one Hamiltonian realization at odd order and at most two at
even order, classifies the alternating-block twins, and certifies their
first square-sum occurrence at 36. Its indexed common-port equations give
a precise route to the coupled integer decoder, rather than a coloring
analogy.

## 2. Odd chart groups cannot carry a nontrivial mod-eight character

Every finite tournament automorphism group has odd order: an involution
would exchange a pair and reverse its directed edge, and Cauchy's theorem
excludes even group order. Its subgroups and quotients also have odd order.
Consequently every group homomorphism from such a chart group to the
additive group Z/8 is trivial: its image order divides both an odd integer
and 8.

In the incoming character criterion, an additive chart function
`Lx mod 8` is well-defined on cyclic charts of orders `m_i` only when
`8 divides L_i*m_i` for each i. With odd `m_i` this forces `8 divides L_i`.
It is therefore identically zero before any further compatibility test.
Likewise, the pure alignment cokernel formed from odd cyclic groups has
odd order and no nontrivial character to Z/8.

This does **not** make the native mod-eight guard constant. That guard may
be a nonhomomorphic function of actual owners, a path-position statistic,
or a relation involving further data. It must retain those coordinates or
remain an exact boundary filter. A converse branch can introduce an even
action; the odd-group statement does not include it by omission.

## 3. A complete joint Gaussian decoder

The incoming [Gaussian paper, section 4b](../session42_shared_lattices_square_colors_20261023/gaussian/RESULTS.md)
shows why independently choosing representatives of D4 point-orbits loses
the alignment needed for addition. Here is an exact repair.

Let D4 act on ordered pairs of **real** vectors by the same signed
coordinate permutation. For `z=(x,y)` and `w=(u,v)`, the following six
polynomials separate all simultaneous ordered-pair orbits:

```
N_z=x^2+y^2,    Q_z=x^2*y^2,
N_w=u^2+v^2,    Q_w=u^2*v^2,
D=x*u+y*v,     K=x^3*u+y^3*v.
```

All six are invariant under the simultaneous action. To prove separation,
first use `(N_z,Q_z)` to put z uniquely in the chamber `x>=y>=0`.

- If `x>y>0`, the equations `xu+yv=D` and `x^3u+y^3v=K` recover `(u,v)`:
  their determinant is `-xy(x^2-y^2)`, which is nonzero.
- If `x>0,y=0`, D fixes u and N_w fixes `|v|`; the remaining reflection
  of the second coordinate identifies its two signs.
- If `x=y>0`, D fixes `u+v`, and N_w fixes `uv`. The unordered roots u,v
  are recovered; their exchange is precisely the remaining chamber symmetry.
- If z is zero, `(N_w,Q_w)` classifies w under the remaining full D4 action.

Thus every case retains exactly the residual stabilizer, proving the claim.
This is a separating set over the reals, not a claim that these polynomials
generate the invariant ring or constitute a minimal separating set.

The mixed quartic is necessary for this particular five-invariant truncation.
Take `z=(2,1)` and compare `w=(-2,4)` with `w'=(2,-4)`. Both inputs have

```
(N_z,Q_z,N_w,Q_w,D) = (5,4,20,64,0),
```

but K is respectively -12 and 12. Their sums are `(0,5)` and `(4,-3)`:
both have norm 25, yet their quartics are 0 and 144. These are exactly the
two norm-25 orbit types appearing in the incoming single-vector obstruction.
Even individual quartics **plus the additive norm cross term** can lose the
orbit of the sum. The mixed invariant restores the ordered pair, and hence
every equivariant operation on it, including addition.

## 4. The common map is a relative-frame quotient

More generally, fix individual orbit types of z and w under a finite group
G, with stabilizers H_z and H_w. Every ordered pair can be brought to
`(z,g w)` by a simultaneous action. Two such pairs agree up to that action
exactly when their relative g values belong to the same double coset

```
H_z \ G / H_w.
```

This gives a bijection of simultaneous pair-orbits with those double cosets.
For D4 the counts are 8 for two generic types, 3 for axis/axis, 2 for
axis/diagonal, 3 for diagonal/diagonal and 1 when one vector is zero.
The literal eight signed permutations give the finite controls.

This is an exact structural connection to shared attachment charts: both
retain relative alignment under one simultaneous group action. In the
cyclic network specialization, the retained object is the shared incidence
cokernel. The six Gaussian polynomials do not automatically become complete
tournament observers; constructing such a separator would require the
actual cell action, port predicate and realizable frame constraints.

The [checker](check_joint_alignment.py) independently enumerates every ordered
pair in the integer square `[-4,4]^2`, comparing the six-polynomial partition
against literal D4 orbits. The all-real theorem rests on the case proof,
not that finite census. A separate agent reviewed the proof, its boundary
stabilizers and the real/ordered-pair scope before publication.

## 5. Final incoming synthesis through 7e36e38161

The final synchronization added the concurrent
[cyclic-attachment packet](../cyclic_attachment_lattices_20261023/README.md).
Its transporter/HNF/cokernel results overlap the two independently developed
descriptions already credited above. It also authenticates full cyclic
automorphism groups for every odd consecutive circulant cell. Its
[nonabelian extension](../cyclic_attachment_lattices_20261023/NONABELIAN.md)
records a **CITED** polynomial feasibility method for explicitly supplied
group-coset constraints and proves a complete prefix-transversal bank. We
therefore do not classify explicit full-group feasibility as unresolved.
The remaining boundary is implicit native-group authentication, compact
port observations and the effects of additional non-coset restrictions.

Its [charge theorem](../cyclic_attachment_lattices_20261023/CHARGE.md) gives
restricted-Fourier constancy tests for integer-valued local tables, retaining
mixed carries. Such a table need not be a homomorphism to Z/8. Consequently
section 2's trivial-character result does not invalidate that stronger test.
Its [carry construction](../cyclic_attachment_lattices_20261023/CARRIES.md)
retains an actual kernel rotation in the nonsplit `(d,d^2)` cyclic extension
after a quotient return. This is a concrete reason to retain both our quotient
frame and its kernel lift. A bound on Fourier term count alone is not a bit
complexity theorem for exact cyclotomic zero testing.

The incoming [coupling operator](../session42_shared_lattices_square_colors_20261023/coupling/COUPLING_OPERATOR.md)
is exactly our incidence-flow transition P. It supplies a stronger all-size
result: graph-shaped Walsh modes have eigenvalues `mu(H)^2`, and the exact
spectral gap is **8/9 for every n>=3**, before and after complement/converse
merging. The [sign-moment theorem](../session42_shared_lattices_square_colors_20261023/sign-moments/THEOREM.md)
identifies the extremal nonconstant shapes as P3 and C4, with isolates.
The P3 shape already appears in our degree/triangle observable, as the
following explicit map shows.

Put `s_uv=(-1)^g_uv`, `S_v=sum_(u!=v) s_uv=n-1-2d_G(v)`, and

```
u_P3(G) = (sum_v S_v^2 - n(n-1)) / (n(n-1)(n-2)).
```

Expanding the squares shows that this is the average Walsh character over
all two-edge stars: each unordered pair of neighbors occurs twice. Our
earlier expected cyclic-triangle count is

```
m(G) = sum_v d_G(v)(n-1-d_G(v))/6
     = binom(n,3)/4 * (1-u_P3(G)).
```

The incoming spectral theorem therefore gives the exact all-size identity

```
P (m-binom(n,3)/4) = (m-binom(n,3)/4)/9,       n>=3.
```

Thus the degree/triangle observable is an extremal slow mode of this exact
coupling. This is a proved algebraic corollary of the cited incoming theorem,
not a fitted spectral coincidence. It is not a native Hamiltonian-path
metagraph spectral theorem. Nor does a constant gap give a uniform total
variation mixing time from arbitrarily rare graph classes.

The incoming [palindromic-clock theorem](../session42_shared_lattices_square_colors_20261023/PALINDROMIC_SQUARE_CLOCKS.md)
adds a quotient Hamiltonian-path test with a marked central endpoint and
first occurrence at 48, distinct from the first twin at 36. Forgetting that
marked endpoint falsely accepts its `(40,4)` control. This reinforces our
indexed-port requirement; a piecewise block involution is not automatically
an equality in the signed-affine compiler without retaining its branch.

Finally, the incoming [Fibonacci addition theorem](../session42_shared_lattices_square_colors_20261023/fibonacci-addition/RESULTS.md)
and [square-count distribution](../session42_shared_lattices_square_colors_20261023/fibonacci-addition/SQUARE_COUNT_DISTRIBUTION.md)
supply a genuine Fourier/addition bridge: canonical colored pair counts
depend on the one common phase `{N*phi^(-2)}`, with a logarithmic error bound.
Their nine limiting coordinates are functions of that same phase, not nine
independent draws. Every legal branch still contains the same-colored square
edge `2+7=9`, so these are ancestry/addition colors, not proper graph colors.
Branch selectors must also remain common when several target sums are queried.
These canonical-distribution statements give neither uniform per-source
positivity nor a map to prime-restricted or native tournament problems.

Three agents made this final bounded comparison without rerunning or
rebranding the incoming research. Our frozen computational inputs are local
to this packet; these additional inherited theorems are linked with their
scope, not silently absorbed into a new empirical claim.
