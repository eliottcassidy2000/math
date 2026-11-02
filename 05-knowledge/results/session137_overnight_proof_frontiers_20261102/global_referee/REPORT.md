# Independent adversarial audit: global tournament determinant bound

**Verdict: ACCEPT.** The candidate in `../discriminant/GLOBAL_PROOF.md` proves the claimed global bound, including the weighted extension. Its ordinary equality cases are correct. The stronger proper-cycle gap and its packing equality criterion are also correct. This verdict is based on an independent reconstruction, not another reviewer's approval or finite extrapolation.

## 1. Independent reconstruction and delicate cases

Let `u=det(I+A)`, `v=(-1)^n det A`, and `delta=H-u`. For `n>=1`, set `M=I+A-A^T` and `D=det(M)/2^(n-1)`. Because `M+J=2(I+A)` and `M-J=-2A^T`, the matrix determinant lemma gives

```
u+v=D>0,                 u-v=sD>0,
s=1^T M^(-1)1=||M^(-1)1||^2>0.
```

The same assertions apply to real weighted tournament matrices. At the empty matrix, `H=u=v=1`: the two positivity quantities are 2 and 0, not both strictly positive. Crucially, the induction needs only `delta(empty)=0`, not the false empty inequality `H>=u+|v|`.

Choose a root `x`. Let `C` run over simple directed cycles containing `x`, with cyclic rotations identified. Put `R_C=T-V(C)`, and let `w_C` be the product of its arc weights. For ordinary tournaments every such cycle has length at least three; weighted matrices additionally allow length two. Direct permutation expansion and OCF give

```
u(T)=u(T-x)+sum_C (-1)^(|C|-1) w_C u(R_C),
v(T)=-sum_C w_C v(R_C),
H(T)=H(T-x)+2 sum_(C odd) w_C H(R_C).
```

These formulations independently verify the candidate's endpoint convention: writing the cycle as `x,p1,...,pell,x` means `x->p1` and `pell->x`. Hence the old path begins in the out-neighborhood and ends outside it in the ordinary case. Fixing `x` removes cyclic-rotation multiplicity. In `v`, a full cycle cover with `r` cycles has sign `(-1)^r`, so the new cycle always contributes a minus sign, regardless of length.

For either `epsilon=+1` or `-1`, subtraction yields

```
delta(T)+epsilon v(T)
 =delta(T-x)
  +sum_(C odd) w_C [2 delta(R_C)+u(R_C)-epsilon v(R_C)]
  +sum_(C even) w_C [u(R_C)-epsilon v(R_C)].
```

Induction on the weaker assertion `delta>=0`, together with algebraic positivity, makes every term nonnegative. This proves both signed margins and therefore `H>=u+|v|=max(1,s)D>=D`. The argument works for every chosen root and does not use strong components or any computational cutoff.

For the weighted extension, OCF is obtained by independently orienting each unordered edge with its given probability. Each odd-cycle family uses each unordered edge at most once, as does each Hamiltonian word. Expectations therefore multiply correctly. Determinants themselves are **not** claimed to commute with this randomization: weighted two-cycles are treated directly in determinant expansion. A weighted two-cycle `x,i,x` is even and contributes the beneficial even-cycle term above, with weight `A_xi A_ix`.

## 2. Ordinary equality and the stronger cycle gap

For integer tournament matrices, `u>0` is an integer, so `u>=1`; also `u±v>=1` on every nonempty induced subtournament. Set `C_all(T)` to the number of all simple directed cycles, counted modulo rotation, and `C_proper(T)` to the number not spanning all vertices.

First, the unsigned subtraction identity gives

```
delta(T)=delta(T-x)
 +sum_(C odd through x) [2 delta(R_C)+u(R_C)]
 +sum_(C even through x) u(R_C).
```

Inducting with `u(R)>=1`, including `u(empty)=1`, proves `delta>=C_all`. Now apply the signed identity: the root-avoiding contribution is at least `C_all(T-x)`, every proper cycle through the root contributes at least one, and spanning-cycle contributions are nonnegative. Thus

```
H(T)>=u(T)+|v(T)|+C_proper(T)              (nonempty T).
```

If two vertex-disjoint triangles exist, put the root in one of them. Its odd-cycle summand has a remainder containing the other triangle, so `2 delta(R)>=2 C_all(R)>=2`. Consequently the gap in this strengthened inequality is at least two.

Conversely, if there are no two vertex-disjoint triangles, there are no two vertex-disjoint directed cycles: each cyclic induced subtournament contains a triangle. Write `c_o,c_e,c_n` for odd cycles, even cycles, and spanning cycles. Only single nontrivial cycles can occur in the determinant and OCF collections, so

```
H=1+2c_o,       u=1+c_o-c_e,       v=-c_n.
```

Therefore `H-u-|v|=c_o+c_e-c_n=C_proper`. This proves the exact equality criterion: equality in the strengthened bound iff there are no two vertex-disjoint triangles.

For the original bound, every nontransitive tournament of order at least four has a triangle through a root with nonempty remainder, which makes both signed margins strictly positive. The remaining equality cases are transitive tournaments and C3. Only transitive tournaments satisfy `H=D`.

### Optional weighted equality classification

The candidate deliberately refrains from a classification, but its identity proves one. If `n>=3` and any edge probability is fractional, root at one endpoint. The corresponding positive-weight two-cycle has a nonempty remainder, so both signed margins are strictly positive. Thus for `n>=3`, equality `H=u+|v|` forces a binary tournament and has exactly the ordinary equality cases. For `n=1,2`, equality holds for every weighted tournament matrix. At order two, writing the edge probability as `p`, one has `H=1`, `u=1-p(1-p)`, `v=-p(1-p)`. Equality `H=D` occurs exactly for binary transitive matrices at every positive order; for order two, `D=1-2p(1-p)`.

## 3. Counting-complexity companion

The four nonnegative gaps admit canonical polynomial-length #P certificates. The full recurrence proof, exact fixed-width binary encoding, complexity bounds, and standalone verifier/counter are separated into `COUNTING_GAPS.md` and `counting_gaps.py`. This is a special constructive #P result, not a claim that every nonnegative GapP function is in #P.

## 4. Independent exact computational checks

`check.py` counts Hamiltonian words literally by permutations, enumerates cycles with minimum vertex first, and computes determinants independently by integer Bareiss elimination. Rational weights are scaled by a common denominator; every comparison is exact. It checks all vertex recurrences, the separate odd-cycle collection expansion, empty cases, all three inequalities, both ordinary equality criteria, the stronger packing equality, and the fractional-edge strictness statement.

Completed census/sample:

* All 1,100 labelled tournaments of orders 0 through 5 (including empty).
* 100 random order-six, 20 random order-seven, and 3 random order-eight tournaments.
* Explicit two-C3 examples of orders 6, 7, and 8.
* 180 rational weighted matrices, 30 at each order 1 through 6, entries from `{0,1/4,1/2,3/4,1}`.

All checks pass. Every audit gate uses an explicit `require`/raise and remains active under Python optimization. Normal and `-O` runs produce byte-identical `check.out` and `check.optimized.out`. The script prints its own SHA-256. These checks are diagnostics, not the universal proof.

## 5. Exact primary sources and literature boundary

1. Darij Grinberg and Richard P. Stanley, *The Redei–Berge symmetric function of a directed graph*, arXiv:2307.05569 (10 July 2023), Theorem 1.39 (p.14), Lemma 6.5 (p.52), and equation (42) (p.55): <https://arxiv.org/pdf/2307.05569>. These give the odd-cycle-family formula after specialization. Their specialization naturally counts complement Hamiltonian paths; reversal gives the same count for a tournament, and loops do not affect vertex-distinct paths.
2. John Irving and Mohamed Omar, *Revisiting the Redei-Berge Symmetric Functions via Matrix Algebra*, Electronic Journal of Combinatorics **32**(4) (2025), #P4.43, DOI 10.37236/13841, Corollaries 19–20 (pp.14–15): <https://www.combinatorics.org/ojs/index.php/eljc/article/download/v32i4p43/pdf/>. Corollary 20 directly states the specialized counting formula and attributes it to Grinberg–Stanley Theorem 1.39 and Lemma 6.5. The older arXiv version is titled *A Survey of Redei-Berge Symmetric Functions via Matrix Algebra*: <https://arxiv.org/html/2412.10572v2>.
3. Alexander Divoux and Shouda Wang, *The Maximum Number of Hamiltonian Cycles in Tournaments*, arXiv:2610.23755v1 (16 October 2026), especially Section 7 and Section 10: <https://arxiv.org/html/2610.23755v1>. This is nearby current literature involving determinants and Hamiltonian-path counts, but its normalized Seidel determinant ratio and asymptotic extremal bounds are different from the pointwise lower bound audited here.

The targeted searches and inspected sources did not identify an already-stated equivalent `H>=det(I+A)+|det A|` inequality. This is a bounded search observation, **not** a claim of novelty or priority. The inherited OCF must remain attributed.

## 6. Separate inherited error in THM-1950

The canonical file's prose claiming that disc is super-multiplicative over strong components is false. Its correct two-block law has factor `(1+s1 s2)/2`, which can be below one. This does not affect the new proof or the stated max-kernel reduction.

A concrete strong order-seven witness has adjacency rows

```
0011101
1001100
0101000
0000100
0010000
1111100
0111110
```

It has `(H,u,v,D,s)=(31,5,1,6,2/3)`. Thus the ordered strong-component sum `C=>C` has `D=26<36=D(C)^2`. Direct 14-by-14 Bareiss elimination confirms `det(I+K)=212992=26*2^13`. The incorrect super-multiplicativity prose should be corrected separately and must not be copied into the closure statement.
