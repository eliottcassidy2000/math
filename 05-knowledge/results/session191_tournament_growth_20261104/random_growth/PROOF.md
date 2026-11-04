# Session 191: exact random-path moments and almost-sure factorial growth

**Status:** PROVED, self-contained below; independent finite checks in
`check_moments.py`. This is a uniform-random model, not a
classification of possible path counts and not a proof of the {7,21} spectrum.

## Inheritance and scope

The exact second moment is inherited, not new: see
`01-canon/theorems/THM-589-metagraph-H-variance-closed-form-nud-succession.md`
and its corrected odd-composition identity. Its invalid unweighted recurrence
was withdrawn on 2026-10-14; we do not use it. The orthogonal even-forest
expansion and the leading variance terms are inherited from
`01-canon/theorems/THM-4761-random-tournaments-the-forest-identity-is-the-walsh-chaos-decomposition-of-kappa.md`
and, earlier,
`01-canon/theorems/THM-848-hamiltonian-path-walsh-needle-stalk.md`.

The additions here are an elementary uniform **telescoping tail bound** and a
self-contained **complete-convergence / almost-sure normalization theorem**.
A targeted search of tournament/path theorem and result files found no prior
almost-sure formulation; this is not a literature-priority claim. The proof
does not use the cited Moran central limit theorem or independence across n.

Live concepts: signed successions; even path forests; orthogonal Walsh layers;
quadratic score energy; almost-sure size growth. The hostile controls are n=1,2
(zero variance), n=3,4 (only the quadratic layer), and the known false
unweighted recurrence for W. The score layer is the less-used sidecar that
upgrades the second-moment statement to a strong law.

## 1. Exact moment identity

Let T be a uniformly random labeled tournament on n>=1 vertices: all unordered
pairs are oriented independently and fairly. Let H(T) count directed
Hamiltonian paths, treating a path as its full vertex permutation. Put
mu_n=n!/2^(n-1), kappa=H/mu_n, and m=floor((n-1)/2). Define

    T_(n,j) = ((n-2j)!/n!)
              sum_(r=1)^min(j,n-2j) 2^r C(j-1,r-1) C(n-2j,r),
              1<=j<=m.

Then

    E H = mu_n,
    E H^2 = n! S_n / 4^(n-1),
    S_n = n! + sum_(j=1)^m (n-2j)!
                       sum_r 2^r C(j-1,r-1) C(n-2j,r),
    Var(kappa) = sum_(j=1)^m T_(n,j).

**Proof.** Each specified path requires n-1 independent orientations. For the
second moment fix one permutation as (1,...,n). For a second permutation sigma,
let F_i and R_i indicate whether i,i+1 occur consecutively in forward and
reverse order. Relative to 2^(-2n+2), their joint path probability has weight

    product_(i=1)^(n-1) (1+F_i-R_i).

Indeed a common forward edge doubles the probability, a reverse common edge
makes it zero, and all other edges impose distinct constraints. Expand the
product. Each selected run of consecutive reference edges must appear as one
contiguous block of sigma, either all forward or all reverse: a mixed direction
would require two incoming or two outgoing path edges at an internal vertex.
For a run of length ell, the two orientations have signed sum 1+(-1)^ell.
Thus odd runs cancel and even runs contribute 2. With 2j selected edges and r
even runs, their half-lengths form a composition of j into r positive parts,
counted by C(j-1,r-1). A fixed list of lengths has C(n-2j,r) placements in the
n-1 reference-edge slots. Contracting the runs leaves n-2j blocks, hence
(n-2j)! permutations. The empty selection contributes n!. Sum and then multiply
by n! choices of the first path. This proves every displayed identity. QED.

This also explains the inherited odd-composition expression: the surviving
contracted blocks have odd vertex sizes, with weight 1 at size 1 and weight 2
at every odd size >=3.

## 2. A uniform, one-sided remainder bound

Write (n)_j=n(n-1)...(n-j+1), and A_(n,j)=2^j/[j!(n)_j]. For 1<=j<=m,

    0 < T_(n,j)
      <= A_(n,j) (n-2j)/(n-j)
      <= A_(n,j)-A_(n,j+1).

Consequently, for 1<=J<=m+1,

    0 <= sum_(j=J)^m T_(n,j) <= 2^J/[J!(n)_J].             (1)

In particular, for n>=3,

    2(n-2)/[n(n-1)] <= Var(kappa) < 2/n.                  (2)

The lower bound is an equality at n=3,4, and strict at n>=5.

**Proof.** Replace 2^r by 2^j and use Vandermonde in its weak-composition form:

    sum_r C(j-1,r-1) C(n-2j,r) = C(n-j-1,j).

This bounds T_(n,j) by

    (2^j/j!) (n-2j)(n-j-1)!/n!
      = A_(n,j) (n-2j)/(n-j).

Now A_(n,j+1)/A_(n,j)=2/[(j+1)(n-j)]. The second inequality is exactly
j>=2/(j+1), valid for every j>=1. Telescoping gives (1), and retaining the
strictly positive terminal A_(n,m+1) gives the strict upper bound in (2).
The j=1 term is 2(n-2)/[n(n-1)]. QED.

In particular, truncation to any fixed number of j terms has a rigorous
O(n^(-J)) remainder. Thus finite rational expansions really do prove

    Var(kappa) = 2/n - 14/(3n^3) - 58/(3n^4) + O(n^(-5)),

including the cancellation of n^(-2), rather than merely suggesting them by
formal fixed-shape power counting. For instance, the error after j=1 is at
most 2/[n(n-1)]. Ordinary Chebyshev gives the explicit size-wise bound
P(|H/mu_n-1|>=epsilon) <= 2/(n epsilon^2).

## 3. Orthogonal quadratic layer

For i<j put X_ij=+1 if i->j and -1 otherwise, and X_ji=-X_ij. These are
independent Rademacher signs for unordered pairs. Expanding each path gives

    kappa = (1/n!) sum_sigma product_a (1+X_(sigma_a,sigma_(a+1))).

As above, all odd-length path components cancel. Write

    kappa = 1 + Q_n + R_n,

where Q_n is its homogeneous degree-2 polynomial in the independent signs,
and R_n is the sum of homogeneous degrees >=4. The even-forest coefficient
for a support with 2j edges and r components is, up to its orientation sign,
2^r (n-2j)!/n!. Distinct supports are distinct orthogonal Walsh characters.
Counting these forests with r components gives

    #(supports) = (n)_(2j+r) C(j-1,r-1)/(2^r r!),

so the squared L2 norm of the degree-2j layer is exactly T_(n,j). Therefore

    E Q_n = E R_n = E(Q_n R_n) = 0,
    E Q_n^2 = 2(n-2)/[n(n-1)] <= 2/n,
    E R_n^2 = sum_(j>=2) T_(n,j) <= 2/[n(n-1)].           (3)

For a concrete structural description, if d_v=outdeg(v)-indeg(v), then

    Q_n = 1 - sum_v d_v^2/[n(n-1)].

This follows by summing the signed two-edge paths: their total is
n(n-1)-sum_v d_v^2. Thus degree-2 score-energy fluctuation dominates the
random path count, while every longer even forest together has summable
mean-square size.

## 4. Elementary quadratic fourth-moment lemma

If Q=sum_(a<b) c_ab X_a X_b for independent fair signs, then

    E Q^4 <= 15 (E Q^2)^2.                              (4)

**Proof.** Replace the signs one at a time by independent standard Gaussians.
Conditional on all other variables, Q=A X+B. Its fourth conditional moment
changes from A^4+6A^2B^2+B^4 to 3A^4+6A^2B^2+B^4, so never decreases.
The second moment is unchanged. At the end, write the Gaussian quadratic form
as G=g^T M g, with symmetric zero-diagonal M and trace zero. Diagonalizing gives
G=sum_a lambda_a(g_a^2-1). Direct fourth-moment expansion yields

    E G^2 = 2 sum lambda_a^2,
    E G^4 = 12(sum lambda_a^2)^2 + 48 sum lambda_a^4
          <= 60(sum lambda_a^2)^2 = 15(E G^2)^2.

This proves (4). QED.

## 5. Complete convergence and the strong law across sizes

**Theorem.** Let (T_n)_(n>=1) be ANY random sequence such that each T_n has the
uniform distribution on labeled tournaments of order n. No assumption about
dependence between different orders is needed. Then, for every epsilon>0,

    sum_(n>=3) P(|H(T_n)/mu_n-1|>epsilon) < infinity,

and consequently

    H(T_n)/(n!/2^(n-1)) -> 1 almost surely.               (5)

This applies, in particular, to the nested induced subtournaments on [n] of
one countably infinite tournament with independent fair orientations.

**Proof.** By (3)-(4), E Q_n^4 <= 60/n^2. Since kappa-1=Q_n+R_n, Markov gives

    P(|kappa-1|>epsilon)
      <= 960/(epsilon^4 n^2) + 8/(epsilon^2 n(n-1)).      (6)

Both series converge. Borel-Cantelli, first for each positive reciprocal
integer epsilon and then their countable intersection, proves (5). QED.

In particular, almost surely,

    log H(T_n) = log(n!/2^(n-1)) + o(1)
              = n log n - (1+log 2)n
                + (1/2)log(2 pi n) + log 2 + o(1),
    H(T_n)^(1/n)/n -> 1/(2e),
    log H(T_n)/(n log n) -> 1.

The theorem is stronger than sparse-subsequence interpolation of an O(1/n)
variance bound. The missing ingredient in that weaker route was the quadratic
sidecar, not independence between sizes. These random statements do not
constrain exceptional deterministic tournaments enough to resolve the global
{7,21} question, and make no claim for the random-regular model.

## 6. Reproduction and finite controls

Run `python3 check_moments.py` from this directory (SymPy is required for the symbolic expansion). The standard-library checks
use (a) the even-run formula, (b) an independently enumerated succession sum
over all permutations through n=9, (c) dynamic programming for all odd
compositions through n=9, and (d) Hamiltonian-path subset DP over every labeled
tournament through n=6. SymPy checks the rational expansion through n^(-6).
The exact tail inequality is tested for every admissible cutoff through n=200.
The script generates `checks.json` as a machine-readable certificate; the mathematical proof is
not limited to the tested sizes. The script retains the composition's number
of parts, so it cannot silently inherit the withdrawn scalar recurrence.
