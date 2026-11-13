# Theta arms, Cartan forms, and gasket scales

Status: PROVED IN NOTE for the elementary identities and classification below; classical names and the cluster consequence are CITED. Finite verification is `check.py` / `check.json`. This is an explicit map between objects, not an identification of their path-count predicates.

## T1. A rank-one change turns a theta network into a Cartan star

Let Theta(p,q,r) consist of three internally disjoint paths between hubs u and v, of positive integer edge lengths p,q,r. Parallel edges are allowed if more than one length is 1; the observed square-sum cases are simple. Give every edge unit conductance. Ground v, let L be the grounded Laplacian, and put

`C = L - e_u e_u^T`, with `N = p+q+r-2` rows.

Then C has diagonal 2, and off-diagonal -1 along the star with arm lengths p-1,q-1,r-1. An empty arm is permitted. Eliminating the three positive-definite path blocks A_(p-1), A_(q-1), A_(r-1) gives the one-by-one Schur complement

`chi(p,q,r) = 1/p + 1/q + 1/r - 1`.

Thus

```
det C = pq + pr + qr - pqr;
inertia C = (N,0,0) if chi>0,
            (N-1,0,1) if chi=0,
            (N-1,1,0) if chi<0,
```

where inertia lists positive, negative, zero eigenvalue counts.

Proof. The path block A_(s-1) has determinant s and its inverse entry at the root end is (s-1)/s; both follow by the determinant recurrence d_j=2d_(j-1)-d_(j-2), d_0=1,d_1=2, or by linear interpolation on the arm. Subtracting these three entries from the root diagonal 2 gives chi. Block congruence proves inertia and the determinant product formula. Equivalently, among assignments with u=t,v=0 the theta energy is minimized at linear interpolation along every arm, with minimum `t^2(1/p+1/q+1/r)`. The rank-one subtraction removes t^2. This proves the same threshold electrically.

When chi=0 a positive integer kernel vector has root value M=lcm(p,q,r) and value M(s-j)/s at distance j on the length-s arm; divide by the gcd of all coordinates to make it primitive. This states the equality boundary explicitly.

## T2. The complete triple classification

For 2<=p<=q<=r, chi>0 exactly for

`(2,2,r) (r>=2), (2,3,3), (2,3,4), (2,3,5)`.

It is zero exactly for `(2,3,6), (2,4,4), (3,3,3)`. Every other triple has chi<0 and exactly one negative direction. If p=1, C is the path Cartan matrix A_(q+r-1), so chi>0.

Proof. If p>=3, the sum of reciprocals is at most 1, with equality only at (3,3,3). If p=2 and q>=4, equality is possible only at (2,4,4). If p=2,q=3, the sign is the sign of 6-r. The q=2 case is always positive. These alternatives exhaust sorted triples.

| Arms p,q,r | N | det C | Positive/negative/zero | Classical star |
|---|---:|---:|---|---|
| 1,4,4 | 7 | 8 | 7/0/0 | A7 |
| 2,3,3 | 6 | 3 | 6/0/0 | E6 |
| 2,3,4 | 7 | 2 | 7/0/0 | E7 |
| 2,3,5 | 8 | 1 | 8/0/0 | E8 |
| 2,3,6 | 9 | 0 | 8/0/1 | affine E8 |
| 2,3,7 | 10 | -1 | 9/1/0 | E10 |

The E-series identification is classical; see Miranda–Morrison, [Embeddings of Integral Quadratic Forms](https://www.math.colostate.edu/~miranda/preprints/eiqf.pdf), printed p.16. They use the opposite sign for the definite root lattices; the table here uses diagonal +2. By Fomin–Zelevinsky's [finite-type classification](https://arxiv.org/abs/math/0208229), a cluster algebra from an acyclic orientation of this star is finite type precisely in the positive cases. This consequence concerns the star quiver, not an arbitrary orientation of the original theta graph and not a tournament.

The old [THM-4838 reciprocal-triangle theorem](../../../../01-canon/theorems/THM-4838-von-staudt-clausen-as-reciprocal-sum-thresholds-orbifold-euler-characteristics-triangle-groups-2-3-p-and-the-twelve-behind-the-owners-sets.md) already identifies the signs for {2,3,5}, {2,3,7} and distinguishes their Bernoulli and divisor mechanisms. The present map makes the same chi the actual Schur complement of a network occurring in the square-sum computation.

## T3. The square-sum connection and its exact loss

The [fixed-endpoint reduction](../squares/PROOF.md) leaves, through Q25, the theta kernels `(1,4,4)` and `(2,3,6)`. Thus their Cartan shadows are A7 and affine E8. The map is: retain three arms, ground one hub, and subtract one root diagonal unit. It preserves the reciprocal threshold and determines determinant/inertia; it forgets the outside forced forest and the designated Hamiltonian endpoints.

The retained outside forest is essential. Q25 and Q26 have the same `(1,4,4)` residual kernel but respectively one and two Hamiltonian completions in the exhibited endpoint problems. Moreover Theta(2,3,r), for r=4,6,8, residual interior demand 1 at every internal vertex, and hub demands (2,2), has the same arm parity and two residual degree completions in all three cases, while C is respectively positive definite, semidefinite, and indefinite. These are residual edge choices after forced edges have been accounted for, not ordinary Hamiltonian paths of the bare theta graph. Cartan type is therefore not a Hamiltonian feasibility test. The correct combined state is `(arm parity/lengths, all residual demands, outside connectivity partition)`. The endpoint-relocation theorem in the square note supplies a genuine operation on the discarded cycle-defect sector.

## T4. The electrical {2,3,5} and {2,5,7} families

There is a second exact source of these small triples. Define the first refinement of a q-corner gasket (q>=2) as q copies of K_q: cell i has one outer vertex a_i and the midpoint vertices y_ij (j!=i), with y_ij=y_ji shared by cells i and j. All edges have conductance 1. If the outer potentials are a_i and S=sum a_i, the unique harmonic interior is

`y_ij = (a_i+a_j+S)/(q+2)`.

The outer current at i is `q/(q+2) * (q*a_i-S)`. Thus eliminating all midpoints multiplies the K_q boundary Laplacian by exactly q/(q+2). Proof: substitute into the midpoint equation, whose degree is 2(q-1), then sum the q-1 outgoing currents at each outer vertex. Positivity of the grounded interior Laplacian gives uniqueness.

Under iteration, conductance scales by q/(q+2) per refinement and resistance by (q+2)/q. With q equal-mass cells per level, the product of cell-count and resistance factors is q+2. The geometric, cell-count, and this mass-resistance scale triple is therefore `(2,q,q+2)`: `(2,3,5)` for the triangle, `(2,5,7)` for five corners, and `(2,7,9)` for seven corners. These are derived network scales; no spectral-limit theorem is needed or claimed. The prime triple `{2,3,7}` is not the triangle's scale triple.

This scale construction and T1 share exact Schur elimination, but their scalar outputs are different functions. Likewise E8^3 is a rank-24 even unimodular lattice with norm-2 roots, so it does not become the rootless Leech lattice merely by having rank 24. The [paired-port audit](../ports/PROOF.md) retains the actual code and distance needed to test that analogy.

## Reproduction and stopping boundary

`python -B 05-knowledge/results/triples_ports_sieve_20261112/triples/check.py`

The checker constructs matrices rather than inserting their determinants, compares exact elimination on all sorted triples 1<=p<=q<=r<=12, checks every zero-case null vector, and independently eliminates the first-level gasket for q=2..9. It does not enumerate all square-sum graphs or prove Chan–Pak. The new actionable questions are which graph births preserve these marked states and when cycle-defect promotion can be scheduled repeatedly with a preserved tournament target.
