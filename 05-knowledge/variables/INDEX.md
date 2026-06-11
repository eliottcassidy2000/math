# Variable Registry — Index

Alphabetical index of all mathematical quantities in this project.
Each variable has its own file with full cross-references.

**How to use:** Find a variable below, click to its detail file for equations, values, and links.
**How to add:** Create `variables/VARNAME.md` from the template at bottom, add a row here.

---

## Quick Lookup

| Variable | File | One-line description | Key equations |
|----------|------|---------------------|---------------|
| `A[i][j]` | [tournament-matrix.md](tournament-matrix.md) | Adjacency matrix of tournament T | A[i][j] + A[j][i] = 1 |
| `alpha_k` / `a_k` | [alpha-k.md](alpha-k.md) | Count of disjoint k-tuples of odd cycles | H = 1 + 2*sum(alpha_k) |
| `B[i][j]` | [signed-matrix.md](signed-matrix.md) | Signed adjacency: B = 2A - J | B[i][j] in {-1, +1} |
| `bc33` | [bc33.md](bc33.md) | Disjoint 3-cycle pairs | bc33 = alpha_2 at n=7 |
| `bc35_w` | [bc35-w.md](bc35-w.md) | Weighted (3-cycle, 5-cycle) disjoint pairs | Appears in c_2 at n=9 |
| `B_M` / `W_d` | [merged-bucket-size.md](merged-bucket-size.md) | Merged tiling bucket size and Hamming-layer transport | `sum B_M=2^m`; `row(W_d)=B_M*C(m,d)` |
| `bucket_balance` | [tiling-bucket-balance.md](tiling-bucket-balance.md) | Half-line balance for finite bucket maps and tiling quotients | Lean: `selfHalf+crossHalf=fiber*moves`; tiling: `2*self+cross=|bucket|*|M|` |
| `c_k` | [fourier-coefficients.md](fourier-coefficients.md) | Coefficients in W-polynomial expansion | W(T,r) = sum c_k * r^k |
| `D_k` | [D-k.md](D-k.md) | Degree-2k Fourier component of W | W = sum D_k; D_0 = c_0 |
| `d(T)` | [d-tournament-determinant.md](d-tournament-determinant.md) | Normalized tournament determinant det(I+S)/2^(n-1) | floor 1 ⟺ local order; ceiling ⟺ DRT switching class; E[det]=A000085 |
| `F_k` | [forward-counts.md](forward-counts.md) | Count of HPs with exactly k forward edges | F(T,x) = sum F_k * x^k |
| `fwd(P)` | [forward-edges.md](forward-edges.md) | Number of forward edges in HP P | fwd(P) + bwd(P) = n-1 |
| `B_N(x)` / `c_L` | [good-cut-bucket-polynomial.md](good-cut-bucket-polynomial.md) | Good-cut bucket generating polynomial and connected-run cover counts | `B_N=B_{N-1}+Σ c_L x^L B_{N-L-1}` |
| `g(τ)` | [good-cut-count.md](good-cut-count.md) | Number of cuts crossed by at least one upward tile | `g≠1`; `g(reflect τ)=g(τ)` |
| `H(T)` | [hamiltonian-paths.md](hamiltonian-paths.md) | Total Hamiltonian path count | H = tr(M); H = I(Omega,2); H odd (Redei) |
| `I(G,x)` | [independence-poly.md](independence-poly.md) | Independence polynomial of graph G | I(Omega(T), 2) = H(T) (OCF) |
| `M[a,b]` | [transfer-matrix.md](transfer-matrix.md) | Transfer matrix entry | M symmetric; tr(M) = H; M[a,b] = sum_S (-1)^|S| E_a B_b |
| `mu(C)` | [mu-C.md](mu-C.md) | Weight of odd cycle C in OCF | mu(C) = I(Omega(T-v) restricted, 2) |
| `Omega(T)` | [omega-graph.md](omega-graph.md) | Conflict graph of odd cycles | Vertices = odd directed cycles; edges = shared vertex |
| `delta_proj` | [projection-defect.md](projection-defect.md) | Loss/residue under support, deletion, or cycle-space projections | `loss_v + keep_v = |V(Omega)|` |
| `r_core(s)` | [single-core-cycle-count.md](single-core-cycle-count.md) | Odd-cycle count for a single-core insertion into a transitive tournament | r_core(s)=Σ_{1...0 pairs} 2^{max(gap-1,0)} |
| `rank_res` | [residue-rank.md](residue-rank.md) | Independence rank of a projected support residue | `rank_res(v)=max{k: alpha_k(R_v)>0}` |
| `s_e` | [edge-variable.md](edge-variable.md) | Centered edge variable | s_e = A_e - 1/2; s_e in {-1/2, +1/2} |
| `S(T)` | [signed-hp-permanent.md](signed-hp-permanent.md) | Signed HP permanent | S(T) = sum_P prod B[P_i][P_{i+1}] |
| `t_k` | [cycle-counts.md](cycle-counts.md) | Number of directed k-cycles in T | t_3 = C(n,3)/4 iff regular |
| `W(T,r)` | [W-polynomial.md](W-polynomial.md) | W-polynomial / weighted HP count | W(T,r) = sum_P prod(r + s_e) |

---

## Template for new variable files

```markdown
# Variable: NAME

**Symbol:** `name`
**Type:** integer / polynomial / matrix / graph
**Defined in:** `01-canon/definitions.md` line N; THM-XXX

## Definition
[Precise mathematical definition]

## Values at small n
| n | Range | Notes | Source script |
|---|-------|-------|--------------|
| 3 | ... | ... | `04-computation/...` |

## Equations involving this variable
- EQ-1: [equation] — source: THM-XXX
- EQ-2: [equation] — source: HYP-YYY

## Relationships to other variables
- Related to [VAR](link): [how]

## Tags
#tag1 #tag2 #tag3
```
