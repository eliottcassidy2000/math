# Polynomial tournament-count tails: sharp indefinite spectral realizations

**Status:** elementary proof plus exact rational tests. This applies classical indefinite-inner-product and moment-matrix mechanisms to polynomial counting sequences. It is not a claim of novelty for the general Jordan or moment theory.

## 1. Construction

Let P(n) be a real polynomial of degree d>=1 with leading coefficient a>0. Define c_j=Delta^j P(0) for0<=j<=d, and c_j=0 for j>d. Here Delta P(n)=P(n+1)-P(n), so c_d=d! a>0. Newton's formula gives P(n)=sum_(j=0)^d c_j binom(n,j) for every nonnegative integer n.

On the real vector space with basis e_0,...,e_d, put

N e_i=e_(i+1) for i<d, N e_d=0,
A=I+N,
G_ij=c_(i+j),
[x,y]_G=x^T G y.

Then A^T G=G A. Indeed (G N)_ij=c_(i+j+1)=(N^T G)_ij, including the boundary because all c_j beyond d vanish. G is nonsingular: its entries with i+j>d vanish and its entire reverse diagonal is c_d. Thus

det G=(-1)^(d(d+1)/2) c_d^(d+1).

The G-adjoint of A is therefore well-defined, and A is G-selfadjoint. Expanding (I+N)^n yields

[e_0,A^n e_0]_G=sum_(j=0)^d c_j binom(n,j)=P(n).

All spaces and matrices here are finite-dimensional, so there are no unbounded-operator domain issues.

## 2. Exact signature

Continuously multiply c_0,...,c_(d-1) by t in[0,1], leaving c_d fixed. The determinant formula stays nonzero throughout. Inertia cannot change along a continuous path of real symmetric nonsingular matrices. At t=0, G=c_d J, where J reverses the order of the d+1 basis vectors. Each paired pair of coordinates contributes one positive and one negative eigenvalue; if d+1 is odd the middle coordinate contributes one positive. Therefore

inertia(G)=(ceil((d+1)/2), floor((d+1)/2),0).

This proves existence with exactly floor((d+1)/2) negative directions.

## 3. The negative index and dimension are minimal

Consider ANY scalar realization P(n)=[v,B^n v] in a vector space with a nondegenerate real symmetric or complex Hermitian form, where B is symmetric for that form on a common invariant domain containing all iterates of v. A bounded indefinite-selfadjoint operator is an important special case. Assume its ambient negative index is defined as the supremum of dimensions of negative-definite subspaces.

The Gram matrix of v,Bv,...,B^d v is

H_ij=[B^i v,B^j v]=[v,B^(i+j)v]=P(i+j), 0<=i,j<=d.

For our construction, let V_ji=binom(i,j) for j<=i and0 otherwise. Its columns give A^i e_0 in the original basis. V is triangular with diagonal1, and

H=V^T G V.

Thus H has the same inertia and nonzero determinant as G. This identity concerns only the prescribed moments, so it holds as an identity for the Gram matrix of any realization. The vectors v,...,B^d v must be independent, their span must be nondegenerate, and it contains a negative-definite subspace of dimension floor((d+1)/2). Consequently every realization requires at least d+1 total dimensions and at least floor((d+1)/2) negative directions. Our construction attains both lower bounds simultaneously.

For a degree0 positive constant, the same formulas reduce to a one-dimensional positive realization. For negative leading coefficient the signature is reversed; that case is outside the positive counting-tail hypothesis here.

## 4. What 'Lorentzian' does and does not mean here

For degree1 or2, the minimum negative index is1. Hence a one-negative-direction realization exists, of signature(1,1) or(2,1). For degree at least3, a one-negative-direction ambient space cannot realize these same scalar moments. This is stronger than merely saying the moments are not positive-definite: it quantifies their minimum signature cost.

In particular, the degree1 sequence associated with the inherited h=63 tail is P(n)=2(n+1). It has

c_0=2, c_1=2,
G=[[2,2],[2,0]],
A=[[1,0],[1,1]],
v=(1,0)^T.

One checks A^T G=G A, [v,A^n v]_G=2(n+1), det G=-4 and signature(1,1). By contrast, a positive Hilbert scalar moment realization is impossible because P(0)P(2)-P(1)^2=2*6-4²=-4.

This is a real bilinear-form statement. It does NOT identify G with an even unimodular integral lattice, II_(25,1), a Leech cusp, or a scattering operator. Even integral P(n) only gives an integral G with determinant constrained as above, not unimodularity. The explicit h=63 G has determinant-4, so it is not an integral unimodular rank2 lattice in this basis or any integral unimodular change of basis. Real signature alone does not supply arithmetic lattice structure.

If a tournament-count sequence is only EVENTUALLY polynomial, shift to a starting order where the polynomial identity holds and apply this theorem to that tail. Shifting does not change its degree or positive leading coefficient. This does not automatically realize a finite exceptional prefix with the same negative index.

### Application to the earlier repunit counting theorem

Session27's proved ordered-strong-component formula makes the ordinary isomorphism-class count N_h(n) and reciprocal-automorphism-weighted count W_h(n) eventually polynomial of degree w(h), where w(h) is the maximal number of nontrivial actually attained strong factors whose Hamiltonian counts multiply to h. This application inherits that theorem and its attainment/gap inputs; it does not reprove the strong spectrum.

For h=31,63,127 the earlier result gives w(h)=1. Thus each of these polynomial tails has minimal negative index1, even though63 is composite. For h=255=3*5*17 it gives w(h)=3, so that tail requires at least2 negative directions. The comparison h=257 has degree1 and needs only1. This is a precise bridge between the strong-component factorization degree and the minimum signature of a scalar spectral observer.

The statement concerns N_h and W_h at fixed h as the order grows, not the labelled count n!W_h(n), not a count of strong classes alone, and not automatically a merged-metagraph node count. Minimal signature is an invariant of this observer realization, not a cause of rarity and not an identification with the actual tournament score lattice or II_(25,1). Provenance: Session27 repunit report, sections 'A proved measure-aware explanation', 'Application to31,63,127' and 'Hostile controls and limits'.

## 5. Repeated poles survive through an indefinite Jordan block

For |z|<1,

[v,(I-zA)^(-1)v]_G = sum_(n>=0)P(n)z^n
= sum_(j=0)^d c_j z^j/(1-z)^(j+1).

Because c_d>0, the pole at z=1 has order exactly d+1. A is a single nontrivial Jordan block when d>=1, although G-selfadjoint. There is no contradiction with the positive Hilbert spectral theorem: indefinite self-adjointness does not impose its diagonalizability or positive spectral-measure conclusion. This is the precise survivor of the earlier positive-moment no-go.

The associated functional on polynomials also has the concrete distributional representation

L(f)=sum_(j=0)^d c_j f^(j)(1)/j!,

since L(t^n)=P(n). It is a finite-order distribution supported at1, not a positive measure. In distribution notation it is sum_j (-1)^j c_j delta_1^(j)/j!. A general signed measure representation is a separate question and is not supplied merely by this indefinite operator realization.

## 6. Classical provenance and verification

The construction is a finite Hankel/Krylov realization with the standard reversing-matrix signature of an indefinite-selfadjoint Jordan block. For modern primary literature on the classical flipped-orthogonal canonical-form framework, see Akgul, Minenkova and Olshevsky, *Existence of flipped orthogonal conjugate symmetric Jordan canonical bases for real H-selfadjoint matrices*, https://arxiv.org/abs/2203.09877 . This note proves the needed special case directly; it does not attribute a new general canonical-form theorem to the present session.

The accompanying standard-library Python check uses exact Fraction arithmetic on32 polynomials, four choices at every degree1 through8. It checks Newton coefficients, A^T G=G A, the determinant formula, inertia by rational congruence elimination, H=V^T G V, moments through n=2d+6, and two larger Hankel null directions. Failure checks remain active under python -O; ordinary and optimized JSON outputs agree byte for byte. These finite checks supplement the all-degree proof, not replace it.
