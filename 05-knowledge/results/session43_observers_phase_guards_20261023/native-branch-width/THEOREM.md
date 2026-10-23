# Exponential native coset width with a constant-state transfer compiler

Status: proved on the explicitly constructed regular hosts and their contiguous cyclic-cut path families. This is an extension of Section5 of the Session43 native observer packet, sharpened by the Session43 root synthesis. It does not count all Hamilton paths in these hosts or prove a complexity lower bound for arbitrary representations.

## 1. Actual regular tournament carrier

Take any odd t>=3. Let the outer tournament be the regular circulant tournament on cells 0,...,t-1, oriented i->j when j-i mod t lies in {1,...,(t-1)/2}. In particular, 0,1,...,t-1 is an actual outer Hamilton path.

Each cell is a directed C3 with arcs u->u+1. Whenever the outer arc is i->j, put

    (i,u) -> (j,v) iff v-u mod3 belongs to {1,2}.

Use the opposite orientation on all remaining cross pairs. Every vertex has one internal out-neighbor, two out-neighbors in each of its (t-1)/2 outer-outgoing cells, and one in each of its (t-1)/2 outer-incoming cells. Its score is therefore (3t-1)/2. The host is a genuine regular tournament of order3t.

For phases a=(a_0,...,a_(t-1)) in F3^t, concatenate the actual internal path in cell i starting at a_i, in the natural cell order. Its boundary from i to i+1 is valid exactly when a_(i+1)-a_i belongs to {0,1}. All internal arcs are already valid. Thus the valid phase set is

    S_t = { a : Delta(a) belongs to {0,1}^(t-1) },
    Delta(a)_i = a_(i+1)-a_i.

Delta is onto and its kernel is the diagonal common-phase group of order3. The family therefore has exactly 3*2^(t-1) distinct literal Hamilton paths.

## 2. Exact lower bound on ANY affine-coset union

**Theorem.** The minimum number of affine subgroup cosets of F3^t whose union is exactly S_t is 2^(t-1).

Proof. Every coset C=a+H in an exact union must itself lie inside S_t. Its image Delta(C) is a coset of the subgroup Delta(H), contained in the binary box {0,1}^(t-1). If Delta(H) is nonzero, some coordinate projection of it is the whole F3. The corresponding coordinate of Delta(C) then contains all three values, contradicting the binary box. Hence Delta(H)=0, so H is contained in ker Delta and |C|<=3.

Since |S_t|=3*2^(t-1), at least 2^(t-1) cosets are necessary, even if overlaps are permitted. The 2^(t-1) complete diagonal fibres Delta^-1(d), for binary d, achieve this bound. QED.

This lower bound allows arbitrary affine subgroup cosets, not merely the particular branches initially supplied to the compiler. It therefore defeats a universal polynomial bound on explicit disjoint-coset branch expansion, even on these regular hosts. It does not defeat compact constraints, transfer matrices, circuits, or other representations: S_t itself has a linear-sized list of local constraints.

## 3. Three states preserve shared owners without branch expansion

Let M be the3-by3 integer matrix M_ab=1 exactly when b-a belongs to {0,1} mod3. Equivalently M=I+Q, where Q_ab=1 when b=a+1.

For r=t-1, (M^r)_alpha,beta counts paths in this family with actual start phase alpha and actual last-cell phase beta. The transfer state is the shared actual phase in the current cell; it is not independently resampled at the two ends of a boundary. A pin at an intermediate actual owner can be inserted as a diagonal0/1 matrix at its stage. Several pins, exact integer clocks, and additional finite guard states can be retained by a product-state transfer, at the corresponding cost in state size.

Write delta=beta-alpha mod3. Fourier diagonalization gives the exact count

    A_r(delta) = [2^r + 2 cos(pi*r/3 - 2*pi*delta/3)]/3.

Indeed the eigenvalues are2,1+omega,1+omega^2, with omega=exp(2*pi*i/3), and 1+omega=exp(pi*i/3). Equivalently, entirely over integers,

    A_r(delta) = [2^r + c_(r-2delta mod6)]/3,
    (c_0,...,c_5) = (2,1,-1,-2,-1,1).

The correction has period6. This is an exact third-root-of-unity effect, not a numerical analogy.

Normalization is essential: A_r(delta) fixes the initial phase as well as the endpoint difference. If ONLY the endpoint difference is fixed, leaving the common phase free, the answer is3 A_r(delta). Summing over delta gives2^r for a fixed initial owner and3*2^r with that owner unpinned.

For regular odd-t hosts r is even, so only the even-r subsequence is used; the full period6 identity also applies to arbitrary chain lengths before closing the regular host construction.

## 4. Prime-cell extension

Let p be any odd prime and let E be a nonempty proper subset of F_p. Put D=E+1. Replace each cell by the regular circulant p-vertex tournament with steps1,...,(p-1)/2, which contains the standard cyclic-cut Hamilton paths. For every outer arc i->j, orient cross pairs by v-u in D.

Writing d=|D|, every vertex has score

    (p-1)/2 + [(t-1)/2]d + [(t-1)/2](p-d) = (pt-1)/2.

Thus the host is regular, and the contiguous cyclic-cut family is exactly

    S_(p,t,E) = Delta^-1(E^(t-1)),
    |S_(p,t,E)| = p |E|^(t-1).

The same proof shows the minimum exact affine-coset cover has size |E|^(t-1): every nonzero subgroup image has some coordinate equal to all F_p, impossible inside proper E. The p-state transfer matrix is M_ab=1_(b-a in E). Its eigenvalues are the finite character sums sum_(e in E) omega^(ke), and the exact endpoint count is the inverse discrete Fourier sum. Integer matrix powering avoids numerical evaluation entirely.

Primality matters to this particular lower-bound proof: for composite moduli a nonzero coordinate subgroup need not be the whole cyclic group and may fit inside E. No composite extension is asserted without examining E's subgroup-coset structure.

## 5. Converse, kernel lifts, and guards

Converse sends every actual accepted path to its reversed path in the converse host. The common phase is a retained kernel lift; fixing a first-cell owner kills it, whereas recording only endpoint difference does not. Any further guard, including a labelled mod8 electrical predicate, must be evaluated on its actual lifted state or shown to factor through a smaller state. The transfer construction does not assume that arbitrary such guards are additive, phase-invariant, or finite-state with bounded memory.

The exponential coset-width theorem is therefore a useful two-sided result: a particular explicit branch representation is necessarily large, while the native shared-state path representation and exact count remain small.

## Checks and limits

The checker constructs actual regular hosts of orders9,15,21,27 for p=3, plus two order15 p=5 examples. It enumerates all phase tuples, verifies every literal path boundary, degree balance, transfer counts, and pinned/unpinned normalizations. Independently it constructs ALL subgroups and affine cosets in F3^3 and F5^3, checks the largest contained coset has size p, and verifies the diagonal fibres attain the exact cover bound. The integer period6 formula is checked for r=0,...,30. These finite checks supplement, not replace, the all-size proofs.
