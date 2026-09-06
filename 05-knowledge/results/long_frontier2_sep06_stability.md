# Regular-root obstructions to a global three-atom stability modulus

**Status: PROVED analytic family and barrier statements; FINITE-EXACT actual
hostile; [independent audit passed](long_frontier2_sep06_stability_audit.md).** The best global distance coefficient
is positive and at most

    kappa4=2(sqrt(2)-1)(sqrt(3)-1)/9=0.067383416186...

The matching global lower bound is **OPEN**. The separate finite-length
inequality `N(R-K3)>=C` for N>=6 is still **OPEN**. This note supplies no
finite-N optimizer and no actual Laurent or LRC consequence.

## 1. Inheritance, domain, and live board

Keep the exact finite real-list domain `p1=p2=1`, `E=(1-p4)/2>0` from
[THM-4454 / sharp global signed-root duplication stability](../../01-canon/theorems/THM-4454-sharp-global-signed-root-duplication-stability.md).
With u=sqrt(2), v=sqrt(3), z=1/v and the three largest positive entries
a>=b>=c>=0, zero padded, put

    J=(5-8p3+3p4)/(3(1-p4)), c_*=(13-8u)/3,
    d2=2-u(a+b), R=(J-c_*)/d2,
    K3=4(u-1)v/[3(1+v)],
    Delta=2-2(a+b+c)/v,
    M=p4-2p3/v+1/3.

All three denominators d2, Delta, M are positive on this actual domain.
The nearest mechanisms are
[THM-4455 / minimizing-sequence rigidity](../../01-canon/theorems/THM-4455-three-atom-minimizing-sequence-rigidity.md),
[THM-4456 / finite-length asymptotics](../../01-canon/theorems/THM-4456-sharp-finite-length-signed-root-stability-asymptotics.md),
and the [sharp local modulus](long_frontier_sep06_local_modulus.md).
The incoming [continuing4 synthesis](continuing4_20260906_synthesis.md)
retains the global modulus as open; its moment and coefficient reductions
reinforce the need to keep the information lost by a relaxation.

The board is finite length, the one-atom energy boundary, ordered root
matching, local three-atom response, regular higher multiplicities, and
signed normalization dust. The canonical hostile is the
[length-four/five competing boundary](long_frontier_sep06_finite_dimension.md).
The corrected near miss is promoting a local sharp coefficient to a global
one. The least-used sidecar here is the integer number of nonzero roots.
Targeted current-repository searches found no duplicate of the four-root
constant below; no external priority claim is made.

## 2. A general normalization lift

Let q be any finite real list with sum q_i^2=1 and sum q_i=S. For an integer
n>=2 define

    alpha_n=[S+sqrt(n(n+S^2-1))]/(n+S^2),
    beta_n=(1-alpha_n*S)/n,
    r^(n)=(alpha_n*q_i)_i concatenated with n copies of beta_n. (1)

Then alpha_n>0, and `sum r_i=sum r_i^2=1` exactly. Indeed alpha_n is the
positive root of

    (n+S^2)alpha^2-2S alpha+1-n=0,

which is equivalent to alpha^2+n beta^2=1 after the first-moment equation.
Also alpha_n->1, beta_n=O(1/n), and n beta_n^2->0. The resulting lists
converge in permutation square norm to q; their third and fourth power
sums and ordered positive coordinates converge to those of q. Whenever
all denominators of a specified quotient are positive at q, that quotient
therefore converges without a singular substitution. The dust has limiting
signed first moment 1-S, which need not be zero.

This proves that an arbitrary unit-square root configuration is a lawful
limiting object, rather than asserting that its own uncorrected first sum
is one. For the positive regular configurations below, S>1 and beta_n<0.
Their lifts have positive energy and all required denominators positive.

## 3. Every regular multiplicity and the four-root barrier

Take an integer m>=4 and q consisting of m positive entries x=1/sqrt(m).
Its limiting power sums are p3=x, p4=x^2; the three leading coordinates
are all x. Direct cancellation gives

    J=(5-3x)/[3(1+x)],
    R=4(u-1)/[3(1+x)],
    (R-K3)/Delta = 2(u-1)/[3(1+v)(1+x)].             (2)

There is no singular cancellation in this range: x<=1/2, so E,d2 and
Delta are strictly positive. Since x decreases strictly as the integer m
increases, the final expression increases strictly with m. Thus m=4 is
the unique minimizing multiplicity among every regular limit with m>=4,
and its value is exactly kappa4. The actual lifts (1) prove the global
upper bound, rather than a surrogate without p1 normalization.

The local distance coefficient A=4.053382428... is therefore much larger
than any possible global coefficient. The diffuse m->infinity limit gives
only the weaker bound `(4(u-1)/3-K3)/2`; the discrete four-root limit is
the stronger obstruction.

The m=3 configuration is deliberately excluded from this continuity
argument: both Delta and M vanish there. Substituting x=z into a canceled
formula from noninteger multiplicities does not give the actual local
coefficient of a normalized three-root lift; that coefficient is A by
the inherited local theorem.

For completeness the same regular configurations give

    M=(z-x)^2,
    (R-K3)/M = 4(u-1)/[3(1+z)(1+x)(z-x)].           (3)

The denominator `(1+x)(z-x)` strictly decreases for 0<=x<=1/2, since its
derivative is z-1-2x<0. Thus these moment quotients strictly decrease as
m increases, and their infimum is the diffuse barrier

    kappaM <= 4(u-1)/(1+v)=0.6064507456... .         (4)

To realize this upper bound by actual lists, for each m choose a sufficiently
large lift parameter n that its quotient differs from the regular value
by less than 1/m, then let m increase. This diagonal choice retains both
actual normalizations. No uniform estimate in both parameters is assumed.

Both global infima are strictly positive, though no explicit lower value
is proved here. For the distance quotient, if an actual sequence had
`(R-K3)/Delta->0`, then Delta<=2 would imply R->K3. The proved local
liminf bound would give `(R-K3)/Delta>=A-o(1)`, a contradiction. For M,
the uniform bound `M<=4/3+2/sqrt(3)` gives the same implication and the
local lower limit 3A yields the contradiction. Consequently

    0 < inf_actual (R-K3)/Delta <= kappa4,
    0 < inf_actual (R-K3)/M <= 4(u-1)/(1+v).         (5)

Neither right endpoint is claimed to be a global sharp coefficient.

## 4. An entirely rational actual hostile

For every integer n>=1 set

    r_n=(n,n,n,n,1,-2,...,-2)/(2n+1),                (6)

with exactly n copies of -2. Its actual length is n+5 and no coordinate
vanishes. Its numerator sum is 2n+1 and square sum is
`4n^2+1+4n=(2n+1)^2`. Hence p1=p2=1 exactly. Its largest three positive
entries are n/(2n+1), also at n=1 where positive ties cause no ambiguity.
Its energy and all matched distances are positive.

The complete power sums are

    p3=(4n^3+1-8n)/(2n+1)^3,
    p4=(4n^4+1+16n)/(2n+1)^4,
    J=(7n^3+32n^2+62n+34)/[3(3n^3+8n^2+6n-2)].    (7)

Thus r_n tends to four positive halves, with vanishing tail square mass
and signed dust sum tending to -1. Formula (2) proves that its distance
quotient tends to kappa4. In particular the concrete length-1005 list
r_1000 satisfies the exact rationally enclosed inequality

    74921/10^6 < (R-K3)/Delta < 74922/10^6 < 1/10.  (8)

This is an actual counterexample to global coefficient 1/10, without
using an asymptotic endpoint in place of an eligible finite list. It is
not a counterexample to the unrelated finite-length N>=6 question.

## 5. A precise obstruction to forgetting integer multiplicity

If formula (2) is extended to a formal noninteger number of equal roots,
take x=4/7 and formal multiplicity m=49/16. It has formal p2=1, p3=x,
p4=x^2 and retains the three selected coordinates x. The relaxed distance
quotient is

    21*kappa4/22 < kappa4.                          (9)

This is **not an actual root list**. Equality p4=p3^2 in Cauchy applied to
the vectors `(r_i)` and `(r_i^2)` at p2=1 would force every nonzero actual
entry to equal p3=4/7. Its multiplicity would have to be 49/16, which is
impossible. Equivalently, after selecting three roots 4/7, the remaining
square mass is only 1/49, so it cannot be supplied by even one additional
root 4/7. A count-free weighted-atom relaxation erases this obstruction.

The connection contract is: actual normalized roots map to a limiting
unit-square configuration via (1); square norms, moments and ordered
roots survive, while the first moment is carried by dust. Passing further
to fractional multiplicities loses the integer root count. Formula (9)
is the cheapest exact hostile to that second reduction. A proof of the
candidate global kappa4 lower bound must retain the ordered tail capacity
or an equivalent multiplicity constraint. A quartic majorant with equality
at 0 and 1/2 is a possible next instrument, not a theorem supplied here.

## 6. Verification and stopping scope

The [source](../../04-computation/long_frontier2_sep06_stability.py) checks
the symbolic cancellation and normalization identities, the formal
multiplicity obstruction, and exactly four literal rational rows
n=1,10,100,1000. Ordinary polynomial multiplication truncated after degree
four independently recovers D/E on every actual row. Rational interval
arithmetic with integer-square-root certificates proves the named hostile.
These finite controls do not infer a global lower bound or an optimizer.

```bash
python3 -B 04-computation/long_frontier2_sep06_stability.py
python3 -B -O 04-computation/long_frontier2_sep06_stability.py
```

Normal and optimized runs pass **40 always-active gates** with byte-identical
[frozen output](long_frontier2_sep06_stability.out). The raw LF hashes are

```text
source SHA256 bf0734ba632810b307b0a3dde8b788fcf6dafb4f226ed0c27308d86b7d086f8e
output SHA256 cd9662c07adeca90b24086fd029186f1db13987cf240e20480d651d5a10c588c
semantic SHA256 9d80ce2955f4555bda812c1b4bfbc0439efb081f989e05265d4c061e7fd62ffd
```

The global kappa4 and moment-barrier lower bounds, exact finite-N
optimizers, and the N>=6 coefficient inequality remain open. This bounded
pass changes the permissible global coefficients and identifies the
integer-multiplicity obstruction; it does not reopen the proved local or
finite-length asymptotic theorems.
