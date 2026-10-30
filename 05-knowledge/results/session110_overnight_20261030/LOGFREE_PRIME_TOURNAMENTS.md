# Logarithm-free prime tournaments in the ramified Pisano bank

2026-10-30. **PROVED from a cited discrepancy input, with independent proof review and exact finite cyclotomic checks.** This is an existence theorem, not a claim of an implemented optimizer, a sharp constant, or global literature priority.

## The result

Let p be an odd prime and a an integer with v_p(a²+4)=1. Put

    O=Z_p[mu],  mu²=a mu+1,  delta=2mu-a,  g=mu⁴,
    R_k=O/(p^k),  N=p^(2k).

Necessarily p>=5. There exists a recursively balanced family of g-invariant Cayley tournaments T_k on R_k with skew adjacency matrices K_k such that

    sqrt(pN) <= ||K_k|| <= 200 sqrt(pN).

The lower bound holds for **every** tournament in this specified bank. The constructed balanced tournaments are prime for k>=2. The estimate is therefore optimal in order of magnitude in p and N, without claiming that200 is remotely sharp.

For fixed p, ||K_k||/N tends to zero: these are quasirandom prime tournaments with the prescribed growing-order arithmetic symmetry. Fully coherent lexicographic towers have linear norm and do not satisfy this conclusion. Balanced lifting is a different, noncoherent operation.

## 1. What is inherited, and what is added

Session109 supplies the ramified orbit lengths and the p-child attachment structure. The present session proves the ideal-module prime criterion and constructs balanced lifts. The earlier random balanced argument gives O(sqrt(pN log N)); the argument here removes that logarithm by combining a new-conductor period bound with block-preserving discrepancy rounding.

The concurrent [pin and reset packet](../pinned_discrepancy_reset_tails_20261029/README.md), at commit11f85d93cc29c48409d218be0c292b3d9c57a4ff, removes a logarithm for fixed bounded orbit size. Our orbit size grows as p^k, so that theorem cannot be applied unchanged. The reduction in the number of distinct Fourier rows and the new period bound are essential additional steps.

### Cited discrepancy input

We use the arbitrary-real-vector partial-colouring lemma in [Lovett–Meka, Theorem4](https://arxiv.org/html/1203.5747#S2), not a set-incidence-only theorem. It controls row changes from an arbitrary cube point when the sum of exponential threshold costs is at most one sixteenth of the number of variables, while placing at least half the coordinates near cube faces. Compactness gives the exact-face existence version.

The following bounded-entry consequence is also derived in the concurrent packet. For an M-by-t real matrix with entries bounded by L, M>=t, a sign vector can be chosen with discrepancy at most

    18 L sqrt(t log(16M/t)).                              (1)

For completeness, partial colouring at n remaining variables costs at most
4L sqrt(n log(16M/n)). Summing over n<=t/2^s and using
sqrt(1+s log2/log(16M/t))<=1+s/8 gives a multiplier
10+(11/2)sqrt2<18. This calculation is not a new discrepancy theorem.

## 2. Orbit and Fourier coordinates

Write D=a²+4. Then delta²=D, O is a ramified quadratic DVR, and

    g-1=a mu² delta,
    v_delta(g^j-1)=1+2v_p(j).

The last formula follows by the first nonzero term in the pth-power binomial expansion, since p>=5. An element z with valuation d<2k has orbit length

    p^floor((2k-d)/2).

There are B_k=p^k-1 antipodal orbit pairs. Reduction R_k to R_(k-1) has p child pairs above every nonzero parent pair, and p-1 additional kernel pairs.

Use the primitive character

    psi(x+y delta)=exp(2 pi i y/p^k).

The pairing (s,x)->psi(sx) is nondegenerate: the delta coefficient of
(s0+s1 delta)(x0+x1 delta) is s0 x1+s1 x0. Thus every additive character
is chi_s(x)=psi(sx), and its conductor is2k-v_delta(s). Fourier eigenvalues
are constant on g-orbits of characters and opposite at antipodal characters.

The two new conductor levels2k-1 and2k have exactly

    M=(p-1)p^(k-1)

antipodal character-orbit representatives in total. There is no need to
treat all N Fourier characters as unrelated rows.

## 3. New-conductor orbit periods

**Lemma.** For a character parameter s of delta-valuation e=0 or1 and any free child orbit O=Hz, H=<g>, its period obeys

    |sum_(x in O) chi_s(x)| <= sqrt(p) N^(1/4).

Consequently the real coefficient of an oriented antipodal packet is at
most B=2 sqrt(p) N^(1/4).

### The analytic calculation

The p-adic logarithm has

    log g=delta ell,  ell in Z_p^×.

Indeed g-1 has delta-valuation1 and all later logarithm terms have larger
valuation. Norm(g)=1 makes log g anti-invariant under conjugation, so it
lies in delta Z_p. The logarithm and exponential converge here because
p>=5.

Let d=v_delta(z), h=floor((d+e)/2), L=k-h, and write
sz=p^h(x+y delta) with residual delta-valuation0 or1. The coefficient-of-delta phase of sz g^t, after removing p^h, is

    f(t)=y cosh(delta ell t)+x sinh(delta ell t)/delta.

If x is a unit, its reduction is affine with unit slope. Hensel lifting
therefore makes f a permutation modulo every p^L, and the complete
nontrivial character sum is zero.

Otherwise p divides x and y is a unit. The restricted integral power series

    F(t)=(f(t)-y)/p

has reduction

    (x ell/p)t + y(D/p)ell² t²/2  modulo p.

Its quadratic coefficient is nonzero. Integrality and disappearance of
higher terms follow from factorial valuations; p>=5 is used here.

For q>=1, expand the squared magnitude of the complete sum of F modulo
p^q by shifts d0. For any nonzero d0 modulo p^q, divide
F(t+d0)-F(t) by p^v_p(d0). The resulting series has affine reduction with
unit slope, hence is a permutation modulo p^(q-v_p(d0)). Its complete
sum vanishes. Only the zero shift remains, so the squared magnitude is
exactly p^q. For f modulo p^L this gives magnitude p^((L+1)/2), also
when L=1.

Summing over all t modulo p^k repeats the reduced phase p^h times.
The stabilizer of z has order p^ceil(d/2), so the orbit period is at most

    p^(h+(L+1)/2-ceil(d/2)) <= p^((k+1)/2).

For free child variables d<=2k-3, hence L>=1. Any constant-phase endpoint
with L=0 belongs to the separately treated kernel; it is not hidden in
this estimate.

### Exact finite checks

[ramified_periods.js](portable_checks/ramified_periods.js) checks the norm
identity algebraically in the prime-power cyclotomic ring. It forms
integer autocorrelation coefficients and tests divisibility by
1+X^(p^(k-1))+...+X^((p-1)p^(k-1)); no floating tolerance is used.

The checks include p=5 at k=1,2,3, and p=13,17 at k=1,2. At the two new
character types the maximal squared orbit periods are respectively
[1,25], [5,125], [25,625] for the golden cases; [13,2197] at p=13,k=2;
and [17,4913] at p=17,k=2. These are finite controls, not the proof of
the all-depth lemma.

## 4. Rounding while preserving every child sum

Suppose a real matrix has n columns in disjoint blocks and entries
bounded by B. Let x in [-1,1]^n have integral block sums of the correct
parity. Existing integral pins may be retained.

Set y=(x+1)/2, so each block sum is an integer. First take y dyadic.
At precision2^q, every block has an even number of odd numerators.
Pair them within their block. Round the two entries of each pair in
opposite directions; the block sum remains unchanged.

At that precision the paired matrix columns are A_i-A_j, bounded by2B.
Apply (1) to these t<=n/2 columns. The y-discrepancy costs at most

    36 B 2^(-q) sqrt(t log(16M/t)).

Convert back to x and sum over q>=1. When M<=n and t<=n/2<=M,
monotonicity gives total error at most

    72 B sqrt((n/2)log32).                               (2)

An arbitrary feasible y is approached by same-block-sum dyadic points
in the integral hypersimplices. The finite set of sign outcomes has a
constant subsequence, proving the exact existence claim by passage to
the limit. Pins never enter a pair. This argument establishes existence,
not an implemented real-arithmetic optimizer.

## 5. Induction without accumulating an error at old frequencies

Fix a parent orientation. In each of its p-child blocks, set

    x_i=parent_sign/p.

The block sum is the required parent sign. The union of its child
orbits is the full reduction preimage of the parent orbit. Every new
character is nontrivial on the kernel, so its coefficients sum to zero
over that block. The actual fractional starting contribution is
therefore zero at every new Fourier row.

There are

    n=p(p^(k-1)-1),   M=(p-1)p^(k-1),

with0.8n<=M<=n for p>=5,k>=2. Combining the period lemma with (2),

    rounding error
    <=72*(2 sqrt(p) N^(1/4))*sqrt((n/2)log32)
    <190 sqrt(pN).

Choose any legal kernel tournament. Its contribution is bounded by
p²-1<=sqrt(pN) for k>=2. Thus all new eigenvalues have magnitude less
than191 sqrt(pN).

At old characters, each child orbit contributes p times its parent
period, and the child signs sum to the parent sign. Therefore the old
eigenvalues are **exactly p times** the parent eigenvalues. Since sqrt N
also grows by p, the normalized old bound is preserved, not increased.

For the base R_1, the unipotent clock forces a product of two
p-vertex circulant tournaments. Applying (1) with L=2 and
M=t=(p-1)/2 gives each factor norm below43 sqrt p. The product norm
is below43 p^(3/2)=43 sqrt(pN). Hence200 closes the induction at every
level.

Balanced means a signed sum constraint, not that all children equal
their parent. Fully coherent lifting instead multiplies old eigenvalues
by p² and retains a macroscopic spectral mode.

## 6. Why every balanced lift at k>=2 is prime

By the independently proved ideal-module criterion, a composite bank
would have an ideal cut delta^r, 1<=r<2k.

For r<=2k-2, take a unit parent difference. Its full reduction fibre
lies in one nonzero coset of delta^r and meets every child orbit.
Module constancy would force all p child signs to agree, contrary to
their sum being plus-or-minus1.

For r=2k-1, take a parent difference of valuation2k-3. Its parent orbit
is a singleton. Each of its p child orbits has length p, and its points
are pairwise distinct modulo delta^(2k-1): for1<=j<p the displacement
(g^j-1)x has valuation2k-2. Thus every top ideal coset meets every
child orbit once. Constancy again forces all child signs equal.
These contradictions exclude every proper ideal cut.

The first level is correctly excluded: its bank consists of products.

## 7. The unavoidable lower scale

For each nonzero residue a in F_p, its fibre contains p^(k-1) unit
orbits of length p^k. Its signed charge is p^k t_a, where t_a is an
odd integer. Opposite residues have opposite charges. The residue-zero
charge vanishes by skewness.

For the p-1 nontrivial characters factoring through F_p, Parseval gives

    sum |lambda|²=p^(2k) p sum_(a!=0) t_a² >=p^(2k) p(p-1).

At least one eigenvalue therefore has magnitude at least sqrt(pN).
This obstruction applies to every bank state, including the constructed
prime states. It is an order-of-growth lower bound, not a claim that
equality is attained.

## 8. Boundaries and consequences

- The proof requires simple ramification. For dual numbers
  R=(Z/q)[epsilon]/epsilon² with g=1+epsilon and q=p^k, every unit
  first-coordinate fibre is one orbit and has charge plus-or-minus q.
  Parseval instead forces
  ||K||>=sqrt(q³ phi(q)/(q-1))=Omega_p(N^(3/4)), N=q².
  Coherent products with low-norm circulant factors attain that order.
  The simple-DVR theorem cannot be transferred unchanged to this family.
- Prescribed child pins require a feasible block-sum fractional point
  and its actual Fourier contribution. The proof does not reset that
  conditioned contribution to zero.
- The theorem supplies existence at the correct scale. The portable
  code verifies character-sum inputs, not a certified implementation of
  the rounding optimizer.
- For any comparison class with nuclear norm at most C N log N, an
  edge-edit route from one of these matrices needs at least
  [N(N-1)/(200 sqrt(pN))-C N log N]_+/4 flips.
  This follows from Frobenius norm squared N(N-1), nuclear norm at least
  Frobenius-squared/operator-norm, and nuclear change4 per reversed edge.
  It is conditional on the comparison-class bound, not a statement about
  every possible tournament grammar.
- No choice of tournament signs licenses choosing branches of the
  deterministic Collatz map. Source weights, guards and clocks remain
  separate arithmetic obligations.
