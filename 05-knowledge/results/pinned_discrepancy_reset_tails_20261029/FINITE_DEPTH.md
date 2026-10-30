# Finite reset depth retains a critical tail

2026-10-29. **PROVED elementary / FINITE-EXACT controls / OPEN uniform depth.**
Companion to [the one-block residual theorem](README.md). Same-agent review;
the exact cylinders are inherited mechanics, not a new parity-vector theorem.

## Inherited mechanism and the quantity actually bounded

[Session29, §§3–5](../session29_framed_profiles_20261022/session29-collatz/REPORT.md)
already proves exact valuation cylinders and explains why exponent-density
cannot be spent as integer-source mass. We retain its cylinder and ROOT
guards, but now sum the original weights `w(j)=1/[j(j+1)]` at source indices
j=(n+1)/2. This is a different measure, not another claim of independent
valuations along a fixed orbit.

For b>=1 let sigma_(b,z) be the residual after b first-reset blocks, killed
as soon as j=1. The discount is z to the **actual number of odd steps**.
Let S_(b,z)(J) be its mass at endpoints greater than J. Every positive
source has a finite first-reset block, so this is defined independently
of eventual Collatz termination. The construction iterates the actual map
on the actual residual, never a refreshed endpoint law.

## Exact word contributions give an all-depth obstruction

A block has valuation word `(1^k,a)`, k>=0,a>=2. A string of b blocks has
odd-step length L=sum(k_i+1) and division clock A=sum(k_i+a_i).
On its unique source class modulo D=2^A, its endpoint has the form

    h=(3^L j+C)/D.                                      (17)

For clarity, a single native letter a maps j to
`(3j+2^(a-1)-1)/2^a`. Thus an exact integer compiler starts with
`M=1,C=0,D=1` and updates

    (M,C,D) -> (3M, 3C+(2^(a-1)-1)D, 2^a D).

Its unique source residue is `j=-C M^(-1) mod D`. Integrality at the
final step forces all previous congruences, by backwards reduction;
the congruence for one letter forces its exact valuation. This also
follows by converting Session29's odd-source class modulo2^(A+1) to j.
Every fixed word loses only finitely many sources to early ROOT arrival,
because every intermediate affine map has positive nonzero slope.

Along a fixed progression the weight is asymptotic to j^(-2), so the
source-weighted endpoint tail has exact leading coefficient

    lim_(J->infinity) J * [this word's discounted tail]
        = z^L 3^L / 4^A.                               (18)

The lattice spacing contributes one factor D; the endpoint slope
contributes 3^L/D. Neither can be discarded. In particular cylinder
frequency alone, 1/D, is not the tail coefficient.

**Proposition B.** For each fixed b>=1 and 0<z<=1,

    liminf_(J->infinity) J S_(b,z)(J)
        >= [z/(4-3z)]^b.                               (19)

**Proof.** The actual b-block words partition the surviving sources.
For any finite bank of these words, sum (18). Its contribution is a
lower bound for the full residual. Exhaust all words. The coefficient
sum is a nonnegative tensor product:

    [sum_(k>=0,a>=2) z^(k+1) 3^(k+1)/4^(k+a)]^b
        = [z/(4-3z)]^b.

This is Fatou's lower-bound direction; it does not require a uniform
tail majorant or exchange a limit with an unbounded signed sum. QED.

At z=1 the lower coefficient is1 at **every finite reset depth**. Thus
no envelopes `S_(b,1)(J)<=c_b/J` valid at all J can have c_b<1, let alone
c_b tending to zero. The one-block theorem proves equality of the
coefficient at b=1; this packet deliberately claims only the lower bound
at b>1. A proof of equality at all finite depths would still not license
interchanging reset depth and source cutoff.

This does not obstruct convergence. Total residual mass can tend to zero
while every finite-depth tail has coefficient1, with its onset escaping
to larger J. Conversely a good asymptotic coefficient alone need not
control the finite sources. The appropriate target is a two-parameter
bound with its onset retained, rather than a decreasing coefficient.

## An independent but degrading upper bound

There is also a finite-scale obstruction stronger than coefficient1.
The original source n=27, j=14, has weight1/210 and reaches j=456 at
its tenth reset. Therefore `455 S_(10,1)(455)>=455/210=13/6`.
Any common envelope `S_(b,1)(J)<=C/J` valid at **all** b and J must
have C>=13/6. This does not impose that lower bound on an envelope
whose onset J_b is allowed to exceed455; (19) still imposes C>=1 there.

Put gamma=log_2(3). For every j>=1, with ROOT held at1 for this estimate,

    B(j)<=R(j)=3^k ell <= j^gamma,

because j=2^k ell and ell>=1. Hence every b-block endpoint is at most
j^(gamma^b). Discounting and killing can only reduce its mass. Therefore,
for J>=1 and any 0<z<=1,

    S_(b,z)(J)
      <= 1/(floor(J^(gamma^(-b)))+1)
      <= J^(-gamma^(-b)).                              (20)

The first expression is the exact original prior tail beyond
J^(gamma^(-b)); it makes no assumption about transported equidistribution.
This is an elementary fallback, not an improvement over the stronger
one-block O(1/J) estimate. Its exponent deteriorates with depth and so
does not close uniform residual control.

The corresponding b-block trial flow is nevertheless summable. At block
i, its length is at most `1+log_2(j_i)<=1+gamma^i log_2(j_0)`.
Consequently its undiscounted occupation mass is at most

    b/2 + ((gamma^b-1)/(gamma-1)) * sum_(j>=2) w(j)log_2 j < infinity.

Thus the signed-residual identity applies at every finite b with no
unstated summability assumption. Its cost is not uniform in b.

## Why a depth-uniform tail is already a major arithmetic target

For this strictly positive original prior, the following are equivalent:

1. `lim_(J->infinity) sup_(b>=0) S_(b,1)(J)=0`.
2. Every individual positive orbit is bounded.

If one reset orbit is unbounded, its original atom w(j)>0 appears beyond
every cutoff at some depth, contradicting1. Conversely, retain finitely
many original sources with omitted total mass below epsilon. If their
orbits are bounded, the union of their values has a finite maximum.
Beyond it all mass comes from the omitted sources and is below epsilon,
uniformly in b. Bounded reset endpoints also bound the intervening native
steps: a bounded endpoint set is finite and each deterministic reset
block has finite length and height. This proves the equivalence for the
actual Collatz orbit, not only the sampled map.

Moreover an all-depth bound C/J would force, for each original j,

    sup_b B^b(j) <= 1+C/w(j)=1+Cj(j+1).

Test at J=h-1 whenever a surviving endpoint is h. Thus this apparently
modest residual estimate would supply a quadratic height bound for every
source. Boundedness still permits a hypothetical nonterminal cycle; it
does not alone prove ROOT arrival. These are elementary quantifier
checks explaining the difficulty of a uniform bound, not new proofs of
boundedness or convergence.

An abstract two-step hostile shows why a C/J target can be unnecessarily
strong. Let1 be terminal, send every odd j>1 to1, and send every even
j=2k to `2^(2k)+1`. Every source terminates within two steps. Nevertheless
at J=2^(2k), its one-step tail is at least `w(2k)`, so
`J S_1(J)>=2^(2k)/[2k(2k+1)]` is unbounded. The same example excludes
every uniform power-tail envelope. It is not the Collatz map; it refutes
the abstract inference from termination to a power bound on excursions.
Source-dependent cutoffs and signed defects may avoid this stronger demand.

## A sharper next target and its falsification rule

The viable helper question is now: can one give explicit, independent
onsets J_b and residual masses m_b such that

    S_(b,1)(J) <= min(m_b, C/J) for every J>=J_b,
    m_b -> 0,

with a controlled finite-region contribution below J_b? The constant C
must be at least1 by (19); nothing here proves such an envelope. To avoid
merely renaming convergence, J_b and the finite-region payments must be
derived without using unproved whole-orbit stopping times. A more flexible
route retains signed residual credit and certifies its total positive
mass directly. Either route must carry source weights, residue phases and
the true clock through every refinement.

The tournament parallel is precise only at the level of proof method:
its fractional barrier isolates the non-roundable pin debt; (19) isolates
the non-contracting Collatz tail coefficient. A theorem that permits
choosing tournament signs does not permit choosing Collatz branches.

`check_words.py` tests 819 words from k=0..2,a=2..4 at depths1..3, both
discounts1 and3/4, perturbs a terminal guard bit, and checks literal
source/endpoint equality through four resets on10,000 original sources.
It also independently constructs the actual625-vertex coherent ramified
tournament and verifies both invariant-space identities behind (16).
Normal and optimized receipts must agree; proofs, not the finite census,
establish (19) and (20).
