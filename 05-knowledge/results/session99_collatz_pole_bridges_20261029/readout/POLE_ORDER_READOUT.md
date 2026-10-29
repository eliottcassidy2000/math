# Reading the Collatz step and its clock from tournament singularities

**PROVED:** exact readout conjugacy, amplitude clock identities and hostile
controls. This is a reformulation of arithmetic inside a proved singularity
atlas, not a proof of Collatz convergence or an identification with analytic
continuation/dilation dynamics. No originality claim is made for elementary
factorial valuation identities.

## 1. The singularity data we actually possess

For the Session98 normalized tournament function

    B(z)=sum_n [n! A000568(n)/2^(n choose2)] z^n,

the pole at p_k=4^k has order d_k=3k+1 and leading scaled amplitude

    A_k=lim_(z->4^k)(1-z/4^k)^(3k+1)B(z)
       =(3k)!/(3^k k!) * 2^(3k^2+k).                    (1)

Both the order and amplitude are inherited proved formulas, not values
fitted to computations. Their mechanism is the unique maximal 3-cycle core
at each symmetry deficit. The pole order has no relation to an order of an
automorphism group unless a separate map is supplied.

## 2. A literal arithmetic conjugacy, with its construction exposed

Let oddpart(m)=m/2^v2(m), for positive integers m. On the subset of poles
p_k whose index k is positive and odd, define the singularity readout

    Phi(p_k)=4^oddpart(ord_(p_k) B).                     (2)

Then, with S(k)=(3k+1)/2^v2(3k+1),

    Phi(4^k)=4^S(k).                                    (3)

Thus k ->4^k conjugates odd Syracuse dynamics to (2). This is an exact map:
the old source index, the new index, and the dyadic clock are recoverable.
Every positive odd source reaches 1 iff every pole in this odd-index subset
reaches the fixed pole4 under Phi. The pole at1 (index0) is outside this
positive-odd domain. Ordinary Collatz steps are recovered by one multiplication
step followed by a=v2(d_k) separate divisions; their cost is a+1, not one.

Proof is substitution of d_k=3k+1. Its very shortness is the limitation:
we have DEFINED a readout using pole order and odd-part extraction. We have
not shown that analytic continuation, residue calculus or the core-insertion
operator performs this extraction or contracts Phi. A conjugate reformulation
does not remove the original orbit problem.

The genuine core-growth operator in ../poles/POLE_TRANSPORT.md sends
k to k+1, pole4^k to pole4^(k+1), and order d to d+3. This differs from
Phi: for k=3, core growth gives index4 whereas Phi gives index5.
Their clocks and sources must not be conflated.

## 3. The amplitude ladder independently exposes the exact valuation clock

Remove the explicit dyadic scale from (1):

    D_k=A_k/2^(3k^2+k)=(3k)!/(3^k k!).

This counts permutations of cycle type3^k. Cancelling the factorial ratio gives

    D_(k+1)/D_k=(3k+1)(3k+2),
    A_(k+1)/A_k=2^(6k+4)(3k+1)(3k+2).                  (4)

For odd k, the cofactor3k+2 is odd, hence

    a(k)=v2(3k+1)
        =v2(D_(k+1))-v2(D_k)
        =v2(A_(k+1))-v2(A_k)-(6k+4),                 (5)
    S(k)=d_k/2^a(k).

This is an exact bridge from adjacent leading pole amplitudes to the native
Collatz division clock. Removing the known scale and the adjacent odd
cofactor, rather than merely recognizing the numeral3, is what makes it work.
It does not say adjacency of poles equals adjacency of a Collatz orbit.

The odd-source guard is essential. At even k=2r the increment instead is

    v2(D_(2r+1))-v2(D_(2r))=1+v2(3r+1),              (6)

because3k+1 is odd and3k+2=2(3r+1). Using (5) for even sources would read
the wrong factor's clock. Both branches are checked explicitly.

## 4. Binary carries are the residual arithmetic coordinate

Write s_2(n) for binary digit sum, with s_2(0)=0. Legendre's factorial identity
v2(n!)=n-s_2(n) gives

    v2(D_k)=2k+s_2(k)-s_2(3k).

For H(k)=s_2(k)-s_2(3k), odd k therefore satisfies

    a(k)=2+H(k+1)-H(k).                                 (7)

The multiplication-by3 binary carries determine the deviation from the
baseline clock2. In the amplitude formula the large quadratic scaling is
completely explicit; H is the residual digit information, O(log(k+1)).

Equation (7) is a finite difference along consecutive k, NOT a coboundary
along Syracuse steps k->S(k). Summing it on a Collatz trajectory does not
telescope: the terms H(k+1) are not generally the next source's H(S(k)).
This is the exact missing alignment, not a technical nuisance. For example
S(3)=5, while3+1=4. Nor do the pole radius or amplitude supply a monotone
height: S(3)=5 increases both, while S(5)=1 decreases both.

The repository's THM-4507,
`01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md`,
finite-polynomial-valuation potential obstruction
is relevant as a warning about restricted height corrections, but it is not
silently applied to H: binary digit sums are not that theorem's declared
finite polynomial-valuation feature vector. No impossibility theorem for
all digit-based potentials is claimed here.

## 5. A decisive same-mechanism hostile

The all-ell-cycle subseries derived in ../poles/POLE_TRANSPORT.md has poles
2^((ell-1)k) of order ell*k+1 for each fixed odd ell>=3. The same definition
as (2), using that pole grid, therefore reads the odd map

    k -> oddpart(ell*k+1).

For ell=5 it reads the positive cycle

    13 ->33 ->83 ->13,

with division clocks1,1,5, disjoint from the cycle1->3->1. All the analogous
meromorphic continuation, geometric-slot, and positive-coefficient properties
remain available. They therefore cannot, by themselves, force convergence
to the orbit of1. This tests the proposed inference using the same analytic
construction, not an unrelated numerical counterexample.

An odd affine multiplier also has the same finite Haar valuation-cylinder
masses2^(-sum a_i): its powers are units modulo2. Thus sharing a geometric
valuation law or a ternary-tree residual-slot law does not fix integer orbit
behavior. The multiplier, exact affine carries, source cylinder and order
of operations remain essential.

## 6. Reproducibility and provenance

Run `python3 check_pole_readout.py` and its optimized equivalent. Active checks
compare independent factorial and ladder computations through k=300, actual
large amplitudes through k=50, exact dyadic clocks and digit identities, and
selected finite native Collatz paths. Selected paths are controls, not an
all-source theorem. The5x+1 cycle and wrong-parity use are hostile controls.

The source is Session98 `limits/ALL_POLES_CORE_EXPANSION.md`, pinned at
`d9d7e33349653b6b65dbad477841c9ed518c081a`. The construction here is distinct
from the classical2-adic shift conjugacy and from the base-conversion
simulation of Stérin and Woods, [The Collatz process embeds a base conversion
algorithm](https://arxiv.org/abs/2007.06979). That paper supplies contextual
motivation for preserving binary/ternary conversion data, not a proof of
our readout theorem or of convergence. Every identity used here has a direct
proof above.
