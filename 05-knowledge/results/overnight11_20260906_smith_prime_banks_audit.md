# Independent audit of the full p11 and p13 prime-order Smith laws

**Status: PASS, proof-level audit with independent exact certificates. No
mathematical repair requested.** The complete p11 and p13 partitions, their
all-depth/all-lift scope, the intermediate p13 exceptional ideal, and the
stated relative prime minimality are accepted. No full law at general primes,
other multiplicities, close-pair depths, or larger node sets is inferred.

Producer: `overnight11_20260906_smith_prime_banks.md` and its source, output,
and `_certificates.json`, frozen in `C:/w/overnight11_jets`. This referee
imports no producer code and does not modify any producer artifact. Its JSON
input is the pinned finite polynomial certificate, whose entries are each
independently proved below rather than treated as trusted conclusions.

## 1. Observer, normalization, and inherited endpoints

The domain is the complete coefficient module of degree below 3m, and the
observer contains Hasse derivative orders 0 through m-1 at each of three
nodes, where m=(p+1)/2 and p is 11 or 13. Translation, permutation of node
banks, and a p-adic unit change of variable are unimodular on these complete
modules. Their coefficients remain p-integral and their diagonal factors are
units. Thus normalization to 0,p^e,p^e a preserves the p-local Smith form
for every a in Z_p with a and a-1 units. This does not normalize away a
nonunit scale or change the derivative orders.

The m rows at zero give an identity in columns 0 through m-1 and zeros
afterward. Integral row clearing therefore leaves exactly the residual
matrix in producer equation (4). The complete Hasse convention is essential:
the determinant is the product of pair differences to power m^2, with no
extra derivative factorial content. Consequently the total determinant
valuation is 3m^2 e. At e=0 it is a unit, so every factor is a unit.

The terminal two residual ideals are legitimate imports from the independently
audited ninth theorem,
`05-knowledge/results/overnight9_20260906_jets_deuring.md`.
It gives L=(3m-1)e-sigma for e>=1 and hence D_(3m-1)=3m^2e-L.
The residual ideals are D_(m+r), because the first m factors are zero.
This gives exactly E11=91e+sigma,E12=108e at p11 and E13=127e,E14=147e
at p13. No intermediate ideal is recovered from the largest factor alone.

## 2. Complete weight bands and exact polynomial identities

For a residual r-minor every determinant term has the same p-scale
exponent W e, where W is the column sum minus the derivative-order sum.
All remaining coefficients form an integral polynomial in a. The least
column sum uses the first r residual columns, while the greatest derivative
sum uses the r largest members of the multiset in which each order occurs
twice. Therefore

    w_(2s)=3s^2,       w_(2s+1)=3s^2+3s+1.

There is one minimal row set for even r and two for odd r. All choices
outside that minimum have integral weight loss at least one. A weight-one
loss is exhausted by either a one-step column-sum increase or a one-step
row-order decrease. The referee independently enumerates all row subsets
and all column subsets into their first two weight levels, then takes their
permitted products. It obtains exactly the producer's 66 index sets:
26 at p11 and 40 at p13, including every required next-band minor.

Every retained polynomial identity is independently certified by a finite
degree argument. If k of a minor's rows lie at node a, then each determinant
term assigns those rows k distinct columns. If their derivative orders sum
to J, every monomial exponent lies between

    lo=sum(k smallest retained columns)-J,
    hi=sum(k largest retained columns)-J.

Thus P(a)/a^lo has degree at most hi-lo. The referee computes the literal
integer determinant with full-pivot fraction-free elimination at each of
a=1,...,hi-lo+1, and checks the given polynomial at those distinct points.
Together with the verified monomial factor and degree bound, these equalities
prove the entire polynomial identity. They are not sample-based inference.
This differs from the producer's symbolic Laplace expansion and its five
evaluation controls.

All 66 identities pass, using 1,092 determinant values. Their complete integer
coefficient gcds, all constants K_(p,r), and every displayed factorization
through q_0,...,q_5 pass separately. The exact next-band content tables are

| p,r | content valuation 1 | content valuation 2 |
|---|---:|---:|
| 11,7 | 4 | 2 |
| 11,8 | 4 | 1 |
| 13,8 | 0 | 5 |
| 13,9 | 0 | 6 |
| 13,10 | 4 | 1 |
| 13,11 | 6 | 0 |

The derivative/degree shift is also exact: increasing both c and j by one
multiplies an entry by (c+1)/(j+1), so a determinant is multiplied by the
product of column factors divided by the product of row factors. This is
a rational scalar identity between integral polynomials; it need not be a
p-adic unit. All twenty shifts whose index sets are retained at both m6
and m7 are independently checked. The operation transports the packet
shape, not an unchanged Smith form or unchanged p-content.

## 3. All-lift attainment and the exceptional intermediate ideal

For each ordinary rank, removing the exact common p-content from the
minimal polynomials leaves a packet with at least one unit at every
admissible residue. This is verified on the complete finite fields, hence
holds at every p-adic lift, not just the integer representatives. The p11
s4 packet reduces to a unit times a^2-a+1, which has no root modulo eleven;
the first p13 s4 packet reduces to a nonzero constant. These are consistent
with the more general direct unit-packet checks.

The p13 rank-eleven packet is the sole exceptional rank needed here. The
referee independently multiplies the two factored degree-five polynomials
and verifies the exact identity

    q2+a^2 q1
      =13(a-1)(a^2-a+1)(4a^4-2a^3-a^2-2a+4).

The common admissible zero set modulo thirteen is exactly {2,7,12}, and
the right-hand factor after dividing by thirteen has residues 2,1,12 there.
Thus at every lift of these classes both q1,q2 are divisible by thirteen,
but they cannot both be divisible by its square. At every other admissible
class at least one is a unit. This proves the exact joint valuation tau.
Restoring K_(13,11) adds one, giving E11=91e+1+tau.

As a hostile to retaining only q1, the independent controls find

    a=1536,1099,662:  v13(q1)=3, v13(q2)=1,

above residues 2,7,12 respectively. All 507 lifts of these three residues
modulo 13^3 have joint valuation one. The algebraic divided identity, rather
than this lift bank, proves the all-lift assertion.

For a target correction zero, every minor has weight at least w_r. For
correction one, the minimal band has the required p-content and all other
weights cost at least e>=1. For correction two, the minimal band has
valuation at least two, the complete next band contributes at least e+1,
and all omitted weights contribute at least 2e. These are all at least two
when e>=1. This explicitly handles the shallow e=1 boundary. The same
argument with the divided companion applies at exceptional p13 rank eleven.
At every rank a minimum-weight polynomial attains the bound. Both directions
of the ideal valuation statement are therefore proved.

Subtracting consecutive residual ideal valuations gives the complete lists
in producer (1)--(2). Their shallow ties are genuine and retained. The
referee also checks that all adjacent factor gaps at e=1 and all their
affine e-slopes are nonnegative, which certifies ordering at every e>=1.

## 4. Full Hasse matrices and the actual precision hostile

The independent matrix engine begins with the original 3m by 3m Hasse
matrix modulo p^B, using B one greater than its determinant valuation.
It removes unit pivots over Z/p^B. Whenever no unit entry remains, it
divides the entire residual matrix by p and decreases the available
modulus, recording the new valuation layer. Each unit pivot is a unimodular
operation; each common division increments all remaining factor valuations.
The determinant bound guarantees enough precision to identify every factor.
This is modular layer peeling, not the producer's rational minimum-valuation
pivot algorithm, and it does not consult the predicted partition.

All 296 matrices pass: the full admissible parameter bank modulo p^2 at
e=1, every residue at e=0 and e=3 with separate fourth-digit lifts, all six
node permutations with a signed affine unit change, and the two precision
hostiles. Only after elimination does the engine compare against factors
derived independently from cumulative ideal formulas and the inherited
largest-loss theorem.

In particular, p13 nodes (0,13,26) and (0,13,39) have respectively

    (0^7,1,2,4,5,7,8,11,12,13,14,16,16,18,20),
    (0^7,1,2,4,5,7,8,11,12,13,14,15,17,18,20).

Both are Deuring-ordinary, have identical pairwise depth one, and have
largest exponent twenty. Modulo 13^16 the sums of capped exponents are
141 and 140. Since each factor p^s contributes p^min(s,16) to the kernel
of this square observer, the claimed kernel orders 13^141 and 13^140 are
correct. This is an actual observer distinction, not only a selected-minor
distinction.

## 5. Symmetry, minimality scope, and remaining question

The AP bit tau is intrinsic. For residues {0,1,a}, the three midpoint
possibilities are a=2,a=-1,a=1/2, giving {2,12,7} modulo thirteen.
It is stable under common affine units and all node permutations. The
referee checks the complete six cross-ratio transforms, as well as literal
matrix permutations. The Deuring roots are exactly {2,6,10} at p11 and
none at p13, so no higher unit digit enters either claimed law.

The relative least-prime statement is properly narrow. At p3,m2 the
three-node two-jet theorem THM-4429 is metric-only; at p5,m3 the independently
audited seventh three-jet law is metric-only; at p7,m4 the eighth full law
depends only on its Deuring bit; p11 is covered here. Thus p13 is the least
odd prime exhibiting failure of the Deuring bit to determine the full
partition in this exact equilateral prime-order family. This does not claim
minimal node diameter, minimal multiplicity for all primes, or any full
partition theorem at arbitrary primes. Those scope boundaries are accepted.

The failed implication is exact top factor plus metric implies all
intermediate ideals. The strongest survivor remains the ninth top-loss law;
the added coordinate here is the midpoint packet in the eleventh residual
ideal, with its integer content and divided companion. The next genuinely
open object is which higher-prime minimal packets acquire common roots and
which additional weight bands are then required.

## 6. Frozen reproduction

    python -B 04-computation/overnight11_20260906_smith_prime_banks_audit.py
    python -B -O 04-computation/overnight11_20260906_smith_prime_banks_audit.py

Both runs pass **260,363 always-active gates**, with byte-identical LF output.
No repository or producer writes occur. The independent source SHA256 is
`c7e7909f60b94782731d0881e7838a74276bfdf9ce331746ea622c0ce643cd14`;
its output SHA256 is
`5773f57284b42b8c5eabd982fec10210f24ccdc9bb80af14b6772fd02eee3f3b`.

Audited producer source SHA256:
`ed102132d643cbdea698f3705576dced1ac51f6a26cd9a4659390da5e5a453ff`.
Producer output SHA256:
`1ef4498d01078e2085c65c025a73a00c5a342ade1aed188b74e9a45bfc18d621`.
The standalone coefficient certificate SHA256, checked before reading, is
`05f79d3dba323be53a1308ea174f61cc2cb79dff053c61cd4c93ec80a15751ec`.

**Filing:** root integrated these audited artifacts in the eleventh checkpoint;
reproduction commands above are relative to the repository root. Outside-worktree
locations preserve author provenance, not the present reproduction location.
