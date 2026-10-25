# Rooted-forest EGF equals the small mirror base

## Status and scope

An exact structural identification, with branch guards and explicit errors. The underlying Lambert-W and rooted-forest facts are classical; this document does not claim a new classical theorem. The mirror series, its coefficients, convergence, asymptotics, and tail estimates are RECOVERED from Session31 tree-moments/TREE_MOMENT_BRIDGE.md; that source has stronger quantitative bounds and Padé enclosures. The complete mirror-to-rooted-forest bridge is also RECOVERED: Session31/SYNTHESIS.md §1 states it explicitly, and HAMILTONIAN_FOREST_SIZE_BIAS.md connects these coefficients to the expected directed-forest determinant in a tournament ensemble. This packet is a compact independent proof/control, NOT a new Session56 discovery. The triangle rooted/unrooted distinction is a compatibility guard. See SOURCE_PINS.json. No conclusion about CM, MRB defects, or transcendence is strengthened here.

## The identity

Put a_0=1 and a_m=(m+1)^(m-1)/m! for m≥1. Then

F(z)=sum_{m≥0} a_m z^m = exp(T(z)) = T(z)/z,
T(z)=z exp(T(z))=−W_0(−z),

where T(0)=0 and the quotient at zero is its removable value 1. Thus the user's coefficients 1, 3/2, 8/3, 125/24, 54/5, 16807/720 are precisely a_1 through a_6, not a_0 through a_5.

Proof: Lagrange inversion applied to T=z exp(T) gives [z^k]T=k^(k-1)/k!, k≥1. Dividing by z gives [z^m]F=(m+1)^m/(m+1)!=(m+1)^(m-1)/m!. Alternatively, applying inversion directly to exp(T) gives the same result. The defining equation yields exp(T)=T/z and consequently F=exp(zF).

## Mirror branch theorem

For positive real n, the equation r^(1/r)=n^(1/n), r>0, is equivalent to h(r)=h(n), where h(x)=log(x)/x. Its derivative is (1−log(x))/x². Hence h increases on (0,e), decreases on (e,∞), is zero exactly at 1, has maximum 1/e at e, and tends to 0 from above at infinity.

For n>e set z=log(n)/n, so 0<z<1/e. There are exactly two positive solutions: n and one r in (1,e). Since T(z) belongs to (0,1), log(F)=T and log(F)/F=T/exp(T)=z. Therefore the distinct mirror is exactly

r_n=F(log(n)/n)=sum_{m≥0} [(m+1)^(m-1)/m!] (log(n)/n)^m.

For 1<n<e the same F(log(n)/n) equals n, the diagonal solution, and the distinct mirror is −W_{−1}(−z)/z>e. At n=e the branches coalesce and there is no distinct solution. At 0<n≤1 there is only the diagonal positive solution. For n>e the W_{−1} formula returns n. In particular, r_4=2 exactly; n=2 has distinct mirror 4, whereas its principal series returns 2. For integer n≥3 the principal series always gives the distinct mirror.

## Combinatorial and determinant identification

Let f_m be the number of forests on the labeled vertex set [m], with one distinguished root in every component; isolated vertices count as rooted one-vertex components. Adjoining a new vertex 0 and joining it to every distinguished root is a bijection from these forests to labeled trees on {0,...,m}. Cayley's formula yields f_m=(m+1)^(m-1) for m≥1; the empty forest gives f_0=1. Thus F is exactly their exponential generating function, not just a function that also uses W.

Equivalently, the complete-graph Laplacian is L(K_m)=mI−J. The eigenvalues of I+L are 1 on the constant vector and m+1 on its orthogonal complement, so det(I+L(K_m))=(m+1)^(m-1)=f_m for m≥1. The matrix-forest theorem identifies this determinant with rooted forests for general graphs as well.

This does NOT identify unrooted forest counts with these coefficients. For K_3, the unrooted counts by edge cardinality are 1,3,3, totaling 7. Their root-multiplicity weights yield 1,6,9, totaling 16=det(I+L). Therefore a score-fibre model counted by unrooted forests cannot be replaced by this determinant without accounting for the product of component sizes. This is a guard, not a new claim establishing such a score-fibre model.

## Convergence, asymptotics, and explicit errors

For every m≥0,

a_{m+1}/a_m=(1+1/(m+1))^m < e.

The ratio tends to e, giving radius 1/e. Stirling's asymptotic and (1+1/m)^(m-1)→e give

a_m ~ e^(m+1)/(sqrt(2π) m^(3/2)).

In particular the series converges absolutely even on |z|=1/e. At z=1/e, continuity from below gives T=1 and F=e. For 0≤z<1/e and M≥0 the nonnegative remainder satisfies

0≤F(z)−sum_{m=0}^M a_m z^m ≤ a_{M+1}z^(M+1)/(1−ez).

Indeed every subsequent term ratio is at most ez. An exact-arithmetic version replaces e by any explicit rational upper bound E≥e with Ez<1.

Using the standard Stirling lower bound m!≥sqrt(2πm)(m/e)^m and (1+1/m)^(m-1)≤e gives a_m e^(−m)≤e/(sqrt(2π)m^(3/2)). Hence uniformly for |z|≤1/e and M≥1,

|F(z)−sum_{m=0}^M a_m z^m|≤2e/sqrt(2πM).

This follows by comparing sum_{m>M}m^(−3/2) with its integral from M to infinity. It remains valid at the branch point, where the geometric bound is unusable. On 0≤z<1/e one can take the smaller of these two bounds.

## Exact verification and limitations

verify.py checks the initial series through degree 24 using the independent formal identity F=exp(zF); exhaustively enumerates all forests of K_m for m≤6; checks the exact coefficient ratios; and supplies exact rational enclosures using the positive atanh series for log(2). Normal and python -O outputs agree. All acceptance checks use explicit exceptions, not assert.

For log(2), summing j=0,...,K−1 of 2(1/3)^(2j+1)/(2j+1) gives a lower bound, with tail at most 2(1/3)^(2K+1)/((2K+1)(1−1/9)). Also e<11/4 follows from summing through 4! and bounding the remaining factorial tail by 1/100. These give rational z bounds at n=4 and 16 and monotone polynomial/tail enclosures. The n=4 enclosure contains the independently exact value 2. Float decimal conversions in the certificate are illustrations, not certified endpoints; separate outward_lower_30dp and outward_upper_30dp fields are rigorously outward-rounded by integer arithmetic. Finite tests are controls, not proofs of the universal identities above.

## Primary sources

- NIST DLMF §4.13, https://dlmf.nist.gov/4.13, for W branches and the tree function. The series-to-mirror argument is explicitly derived here.
- Pavel Chebotarev and Elena Shamis, Matrix-Forest Theorems, arXiv:math/0602575v2, https://arxiv.org/abs/math/0602575v2, for det(I+L) and rooted forests. The complete-graph specialization and the added-root bijection are explicitly derived here.
