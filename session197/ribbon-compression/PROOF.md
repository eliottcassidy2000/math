# Logarithmic witnesses for bounded-backward-span path spectra

Session 197. Independent structural companion to `../ribbon/PROOF.md`. No global spectrum assertion is used or proved.

A span-`w` tournament is supplied with an ordering in which every backward edge has endpoint distance at most `w`. Assume `w>=1`. Let `H(T)` count its directed Hamiltonian paths.

## 1. Path merging is supermultiplicative

**Lemma.** For every vertex partition `V(T)=V_1 union ... union V_k`,

`H(T) >= product_i H(T[V_i])`.

**Proof.** First take two disjoint directed paths `P,Q`. Merge them by repeatedly comparing their two current heads and outputting the head that points to the other. When one list is exhausted, append the other. Two consecutive outputs from one list follow a path edge. At a switch between lists, the previous head pointed to the current head by the comparison at which the previous head was chosen. The output is therefore a directed path preserving the orders in both input lists.

This operation is deterministic. Applying it repeatedly to one chosen Hamiltonian path from each part yields a Hamiltonian path of their union, preserving the path in every part. Its restriction to each part recovers the input, so different tuples yield different outputs. QED.

Consequently, `k` vertex-disjoint directed triangles imply `H(T)>=3^k`; unused vertices can be made singleton parts. No odd-cycle expansion identity is needed.

## 2. Reduction by strong components

The strong components of a tournament have a unique total ordering in which all intercomponent edges point forward. Every Hamiltonian path visits them in that order, each as one contiguous block. Thus `H(T)` is the product of the counts of the strong components.

For `H(T)>1`, delete every one-vertex strong component. The remaining tournament `T*` is nonempty, this deletion preserves `H(T)`, and every remaining component has at least three vertices. For `H(T)=1`, use a single-vertex witness instead; no empty-graph path-count convention is needed. Deleting vertices and retaining the supplied relative ordering preserves the span bound, since distances between remaining endpoints can only decrease.

Every vertex of a nontrivial strong tournament belongs to a directed triangle. To see this for a vertex `v`, let `O` and `I` be its out- and inneighbors. Both sets are nonempty. If no edge points from `O` to `I`, all edges point from `I` to `O`, and no vertex of `O` can reach `v`, contradicting strong connectivity. Any edge `o -> i` produces `v -> o -> i -> v`.

Hence every vertex of `T*` lies in a directed triangle.

## 3. Locality forces a packing

If a directed triangle has ordered labels `a<b<c`, there are only two cyclic orientations:

- `a -> b -> c -> a`, where the backward edge gives `c-a<=w`;
- `b -> a -> c -> b`, where the two backward edges give `b-a<=w` and `c-b<=w`.

In either case its index diameter is at most `2w`.

Let `N=|V(T*)|` and `L=4w+1`. Take the centers

`1, 1+L, 1+2L, ... <= N`.

Choose a directed triangle containing each center. All vertices of a chosen triangle lie within distance `2w` of its center. The intervals of radius `2w` about successive centers are disjoint, because their centers differ by `4w+1`. Thus these triangles are vertex-disjoint. Their number is `ceil(N/L)`.

Combining with path merging gives

`H(T)=H(T*) >= 3^ceil(N/(4w+1))`.

## 4. Logarithmic witness theorem

**Theorem.** For every integer `h>1` realizable by a span-`w` tournament, there is a realization on at most

`(4w+1) floor(log_3 h)`

vertices. The realization is obtained simply by deleting all singleton strong components from any existing witness. Its Hamiltonian-path count is unchanged and its backward span does not increase.

The integer floor can be calculated without floating point as the unique `k` satisfying `3^k<=h<3^(k+1)`.

This is a uniform all-order theorem. It does not assert that a given integer is realizable. At widths 2 and 3 the respective bounds are `9 floor(log_3 h)` and `13 floor(log_3 h)`.

**Algorithmic corollary.** For every fixed `w`, membership in the span-`w` spectrum is in NP when the target integer is given in binary. A word describing a witness has `O(log h)` vertices and hence `O(log h)` bits for fixed width. The exact frontier transfer in the companion packet verifies its Hamiltonian-path count with a constant-size state vector and polynomially many bit operations, saturating at `h+1` if desired. The value `h=1` has the single-vertex witness.

A direct bounded word search is consequently finite without any appeal to stabilization: it suffices to inspect orders up to the displayed bound. Its crude running-time bound is polynomial in the *numerical value* of `h`, with an exponent depending on `w`. This is not a polynomial-time claim in the binary input size.

## 5. A quantitative width-1 contrast

The companion packet proves that the complete width-1 spectrum is the multiplicative monoid of the sequence

`h_0=h_1=h_2=1`, `h_n=h_(n-1)+h_(n-2)+h_(n-3)` for `n>=3`.

For `n>=3`, `h_n>=3^((n-2)/3)`. Indeed the sequence is nondecreasing, so `h_n>=3h_(n-3)`, and the three base cases `n=3,4,5` satisfy the inequality.

For `X>=3`, any product `product_(n>=3) h_n^(e_n) <= X` therefore satisfies

`sum_(n>=3) (n-2)e_n <= 3 log_3 X`.

Set `A(X)` equal to the number of distinct width-1 spectral values at most `X`. Each value has at least one such exponent vector, and distinct values cannot come from the same vector. With `N=3 log_3 X`, for any `t>0`, the number of admissible vectors is at most

`exp(tN) product_(k>=1) (1-exp(-tk))^(-1)`.

Taking logs of the product and expanding `-log(1-z)` gives

`sum_(m>=1) 1/[m(exp(tm)-1)] <= (1/t) sum_(m>=1) 1/m^2 = pi^2/(6t)`.

Choosing `t=pi/sqrt(6N)` proves

`A(X) <= exp(pi sqrt(2 log_3 X))`.

Thus the width-1 spectrum has zero natural density and is smaller than every positive power of `X` asymptotically. The all-order width-1 obstruction is much stronger than a finite list of exceptional small values. It supplies no corresponding upper bound at widths 2 or 3.

## 6. Verification

`verify_compression.py` independently checks, on all width-1 order-7, width-2 order-7, and width-3 order-6 inputs, that deleting singleton SCCs preserves the path count, the reduced graph remains within the claimed span, every remaining vertex is in a local triangle, the center triangles are disjoint, and the numerical lower bound holds. It also checks the head-merge operation and its injectivity for every split of sampled tournaments of order at most 8. Proofs above establish arbitrary-order validity; these are implementation and edge-case checks.
