# From isothermal stresses to tournament families

**PROVED IN THIS NOTE:** the stress comparator, exact parabola and power-curve
images, and their counting formulas. **CITED:** Brouwer's classification of
locally transitive tournaments. **FINITE-EXACT:** the marked-seed census in
`tournament_receipt.json`. Chan–Pak remains OPEN.

## A literal pairwise observable

At sites z_i let q_i be the complex deviatoric stress from `STRESS.md`. Define

    i -> j iff Im(conj(q_i)*q_j)>0.                       (1)

This is the oriented commutator coefficient of the two traceless symmetric
stress matrices, up to a fixed sign convention. The sites are the vertices;
the plane orientation fixes the gauge. Zero stress difference or parallel
deviators give genuine ties, excluded from the tournament universe below.
No artificial tie-breaking is used. Hydrostatic pressure is invisible.

For q_i=h_i*lambda_i, with lambda_i real,

    Im(conj(q_i)*q_j)
       =lambda_i*lambda_j*Im(conj(h_i)*h_j).              (2)

Thus changing the stress parameters performs a **vertex switch** on a fixed
rank-two determinant tournament: reverse exactly the pairs whose lambda
signs differ. Magnitudes and p0 are lost. This construction preserves the
sign of the actual observable, not physical work, energy or Hamiltonian-path
multiplicity before taking the tournament.

There is a sharp regularity boundary: if the field is required to have
distinct principal stresses **everywhere on one connected domain**, lambda
has constant sign and (2) is independent of all five parameters. The varying
switch families below permit isotropic loci between sample sites, while
requiring every sampled comparison to be nondegenerate. Five free parameters
do not by themselves produce five observable tournament degrees of freedom.

The rank-two local-transitivity theorem is inherited from
[THM-2294, anchored-plucker-tournament-and-kakeya-address-bank](../../../01-canon/theorems/THM-2294-anchored-plucker-tournament-and-kakeya-address-bank.md).
Its short proof applies directly: each in- or out-neighborhood lies in an
open semicircle of vector directions, and angular order makes that
neighborhood transitive.

Choose projective representatives in one open half-plane and sort their
angles. Then (1) has the exact binary-word form

    i -> j (i<j) iff s_i=s_j,                            (3)

where s records which representative was negated. Conversely every such
word is a rank-two realization. Brouwer's
[1980 primary paper](https://ir.cwi.nl/pub/6805/6805D.pdf) proves that these
are all locally transitive tournaments, with isomorphism given by rotations
of the antiperiodic word `s complement(s)`. The resulting count is A000016,
equivalently

    a(n)=sum_(d|n, d odd) phi(d)*2^(n/d)/(2n).

This connects the stress construction to a **specific complement-feedback
shift register**, not to arbitrary computation or FRACTRAN universality.

## Fixed sites versus moving sites

At n fixed sites, each lambda_i is a linear form on the four amplitude
parameters. Their zero hyperplanes divide parameter space into at most

    2*sum_(j=0)^3 binom(n-1,j)

strict sign chambers, interpreting out-of-range binomial coefficients as
zero. Global sign reversal leaves (2) unchanged, so there are at most
`sum_(j=0)^3 binom(n-1,j)` labelled tournament outputs. This is the standard
central-hyperplane recurrence: adding a hyperplane adds the regions of the
induced arrangement in one lower dimension. The bound is attained below.
It is polynomial in n despite the exponentially many possible switches.

Moving the sites is a different operation. In fact **one fixed net and one
fixed stress field realize every locally transitive tournament**, as sites
vary, even with globally distinct principal stresses. Take

    Z=exp(z), h=exp(2z), lambda=1, q=exp(2z), p=0.

For ordered `0<alpha_1<...<alpha_n<pi` and desired switch bits s_i, choose
`z_i=x_i+i(alpha_i+pi*s_i)/2`, with arbitrary real x_i. Then q_i has argument
`alpha_i+pi*s_i`, giving exactly (3). The entire field is compatible and
nonisotropic. The isothermal coordinate is locally regular; its global
injectivity is not needed.

A second realization, useful for the curve tests below, takes on the first quadrant

    Z=z^2/2, h=z^2, lambda=|z|^2-1,
    q=z^2(|z|^2-1), p=-(1/2) Re(z^4).

Choose ordered angles `0<alpha_1<...<alpha_n<pi` and put
`z_i=r_i exp(i alpha_i/2)`, with r_i=2 or1/2 according to the desired sign.
Then arg h_i=alpha_i, and (2) is exactly (3). All sampled stresses are
nonisotropic and all comparisons are tie-free. The field has an isotropic
circle between some sample regions; the prescribed isothermal net extends
through it, while the stress alone does not specify axes on that circle.

This distinguishes three quantities: the five-dimensional field space, the
four-dimensional sign arrangement at fixed sites, and the unbounded family
of tournaments available when the sites themselves are construction data.

## The positive parabola gives an exact infinite family

Fix the same net h=z^2 and sample

    z(t)=t+i t^2, t>0.

Its h arguments `2 arctan(t)` increase inside (0,pi), so the unswitched base
is transitive. The amplitude along this curve is

    lambda(t)=A t^4+(A-2 Im B)t^2+2 Re(B)t+C.             (4)

By Descartes' sign rule it has at most three positive sign changes. Every
word with at most three changes at any prescribed ordered positive sites
is attainable. For zero, one or two changes, put a root in each corresponding
gap and use the resulting polynomial of degree at most two. For three
changes, take the three positive gap roots and a fourth root equal to the
negative of their sum. The monic quartic then has no cubic coefficient and
belongs to (4). Its overall sign can be chosen freely.

Thus the labelled output count is exactly

    1,2,4,8,15,26,42,64,93,130,176,232,...               (5)

for n=1,2,... . These are sign pairs, not isomorphism classes. The first
missing sign pair occurs at five sites: the alternating pattern, which
would require four positive roots.

Compress each maximal constant-sign run in (3) into a transitive block.
With an odd number r of runs the quotient is the regular locally transitive
tournament on r vertices. With even r, the first and last runs form a single
transitive module (last before first); merging them gives the odd quotient
on r-1 vertices. Consequently the **complete unlabelled parabola image** is

    TT_n, or C3[TT_a,TT_b,TT_c], a,b,c>=1.               (6)

Every positive triple occurs: choose a three-run sign word with those block
lengths. The cyclic quotient has no homogeneous pair, so its maximal
transitive blocks and their cyclic order are intrinsic. Two triples give
isomorphic tournaments exactly when they differ by a cyclic rotation.

Let

    L_k(n)=1/k * sum_(d|gcd(n,k)) phi(d)*binom(n/d-1,k/d-1)

for n>=k, and zero otherwise. This counts positive length-k compositions
modulo rotation, by Burnside's lemma. For prime k it simplifies to
`[binom(n-1,k-1)+(k-1)*1_(k|n)]/k`.
The parabola family therefore has exactly **1+L_3(n)** classes at every n.
The first twelve values are

    1,1,2,2,3,5,6,8,11,13,16,20.

This is a genuine efficient enumeration of a defined infinite subclass of
A000568. It does not change the count of all tournaments.

## Cubics and every higher integer power add precisely one core

Now sample `z(t)=t+i t^m`, for any fixed integer m>=3. The phase remains
strictly ordered, and

    lambda(t)=A(t^(2m)+t^2)-2 Im(B)t^m+2 Re(B)t+C.       (7)

It has at most four positive sign changes. Four are actually possible for
every m>=3: use

    f_m(t)=t^(2m)-3t^m+t^2-t/100+1/1000000.

At `t=10^-8,1/200,1/50,1,2` its signs are `+,-,+,-,+`. The same signs follow
for every m>=3 by bounding the small powers by their m=3 values; at2 the
leading term dominates. Thus by placing any desired numbers of sites in
these five sign intervals, every positive five-run block pattern occurs.
Patterns with fewer runs occur by using fewer intervals.

The complete image **as sites and stress parameters vary** is exactly

    TT_n, C3[TT_a,TT_b,TT_c],
    or C5[TT_a1,...,TT_a5], all block lengths positive.   (8)

There are exactly **1+L_3(n)+L_5(n)** isomorphism classes. The first twelve
values are

    1,1,2,2,4,6,9,15,25,39,58,86.

The receipt also records a separate rational cubic witness with coefficient
3 replaced by10. It checks the same five-sign phenomenon; it is not the
uniform-in-m proof above. We do not claim every four-change word is available
at a *fixed* set of cubic sample sites.

The striking threshold is geometric: the quadratic curve permits the C3
core, while every integer power m>=3 permits exactly the C3/C5 palette.
Increasing m alone introduces no C7 core in this experiment. More complicated
sampling curves can introduce more crossings, and arbitrary sites recover
all locally transitive tournaments as shown above.

## What this does and does not do for odd path counts

No locally transitive tournament has H=13. Here is a short all-order reduction.
A nontransitive locally transitive tournament is strong: otherwise a cyclic
triangle in one strong component is contained in the in- or out-neighborhood
of a vertex of a different component, a contradiction. The exact local
spectra through order five are

    n=1,2: {1}; n=3: {1,3}; n=4: {1,5};
    n=5: {1,9,11,15}.

Every strong tournament of order at least six has H>=15, using Moon's strong
subtournament theorem, insertion monotonicity and the complete order-six
minimum in [THM-1370, h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
Therefore13 is absent at all orders. The five-dimensional stress family does
not directly yield every allowable odd tournament count.

It may still supply **marked seeds for operations that leave this family**.
That is a different question. The exact census tests every locally transitive
isomorphism type and every nonempty selected set through order fourteen, with
actual automorphism quotienting. It retains

    delta=H-sum_{a<b in O} W_ab, U=sum_{a in O} U_a.

Unit seeds have delta=1 and yield the inherited tail counts `1+2^t U` for
t>=1. There are328 marked unit-seed types in this finite universe. Of these,
324 lie in (6), and four lie in the C5 palette. The first two, at order nine,
have H=705 and U=6288 or2718, with odd parts393 and1359. Two more at order
eleven have H=54123 and U=358600 or236287. This is a
specific next construction opportunity, not a coverage theorem or a novelty
claim about prior repository discovery. Their recursive behavior is treated
separately rather than inferred from one successful unit defect. The
compiled continuation in `round_extension.json` extends this census;
[SEEDS](SEEDS.md) gives the all-length families and exact obstructions.

## Connection contract and controls

Source: actual elastic stress samples. Target: a labelled tournament.
Map: (1). Preserved: sign of the commutator; parameter changes become actual
vertex switches. Lost: stress magnitude, hydrostatic pressure, physical
distance, isotropic intermediate loci and displacement periods. Sidecars:
the sample positions, h values, lambda values, and the field's period data.
Hostiles: ties, the square sensor failure, the missing alternating word on
the parabola, and the all-order H13 obstruction.

The finite engine compares the binary-word quotient with graph isomorphism
through order six and with the entire labelled tournament universe through
order five. Endpoint-cover counts are independently checked by permutation
and cut enumeration through order seven. Normal Python computes the original
census through14. Optimized replay checks its prefix through8; the separately
compiled engine rechecks all levels through14 before extending to18. Its
optimized replay also reaches14; the expensive15..18 extension ran once.
See `replay_receipt.json`. Universal claims follow from proofs, not cutoffs.

## Exact enumeration across all sizes

For prime k=3 or5, positive k-part compositions modulo cyclic rotation give

    L_k(n)=[binom(n-1,k-1)+(k-1)1_(k divides n)]/k, n>=k,
    L_k(n)=0, n<k.

Burnside proves the formula: the identity fixes every composition, and each
nonidentity rotation fixes exactly the constant composition when k divides n.
No reflection quotient is taken. The previously proved sampling images give

    P_n=1+L_3(n),          Q_n=1+L_3(n)+L_5(n), n>=1.

Consequently

    P(z)=z/(1-z)+(z^3/(1-z)^3+2z^3/(1-z^3))/3,
    Q(z)=P(z)+(z^5/(1-z)^5+4z^5/(1-z^5))/5.

These are all-order laws, with eventual quasipolynomial periods dividing3
and15, respectively; retain the denominators and small-order boundary.
At n18 the counts are47 and523, compared with7286 unrestricted locally
transitive classes. At n64 they are652 and119785. These are family sequences,
not new values or replacements for A000568.

`aspect_sequences.py` independently enumerates cyclic compositions through18
and matches the binary-word class census through14. It also evaluates the
classical general tournament Burnside formula through30 as a scale control;
that formula and those known terms are not claimed as new work.
