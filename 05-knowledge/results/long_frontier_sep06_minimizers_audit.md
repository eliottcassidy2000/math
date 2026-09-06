# Independent audit of the three-atom minimizing-sequence classification

**Status: INDEPENDENT ANALYTIC AND SOURCE AUDIT — PASS.**

Auditor: `three_ray_geometry`, separate from the proof/source author, root.
This audit accepts the complete
[minimizer classification](long_frontier_sep06_minimizers.md), relative to
**THM-4454**,
[sharp global signed-root duplication stability](../../01-canon/theorems/THM-4454-sharp-global-signed-root-duplication-stability.md).
The auditor previously proved the regional packing part of that theorem;
root and another independent reader audited that regional proof. The
present boundary analysis and classification are new root-authored work.

## 1. Statement, normalization and exact scope

The theorem concerns arbitrary sequences of finite real lists with
p1=p2=1 and E=(1-p4)/2>0. Lengths and all sign patterns may vary.
The quotient R retains both E and the distance d2 in its denominator.
Its approach to K3 is equivalent to square-norm convergence, after
permutation and zero padding, to three positive entries 1/sqrt(3).

I checked the exact squared-distance formula and the treatment of fewer
than three positive entries. The chosen target slots are the largest
positive coordinates, with padded zeros available. Thus the distance
decomposes into the three matched squared errors and every remaining
squared coordinate. There is no ambiguity from repeated entries.

The theorem preserves the signed first-moment condition. It concludes
only the signed sum of the remaining entries, not the separate positive
and negative sums. Its actual mixed-dust example demonstrates the
distinction. No claim is transported to actual Laurent return rows.

## 2. Independent checks of the two boundary expansions

For the one-positive-atom boundary, let delta=1-a^2. The bounds

    |sum_tail r_i^3|<=delta^(3/2),
    sum_tail r_i^4<=delta^2

hold independently of length and signs. Expanding a=sqrt(1-delta)
therefore gives exactly the stated expansions of p3,p4,D,E. Since an
eligible list has delta>0, division yields J->1 and

    R->sqrt(2)-2/3>1/2>K3.

For the two-equal-positive-atom boundary, I independently used

    t=a+b=sqrt(2-2delta-w^2),
    a^3+b^3=(t^3+3tw^2)/4,
    a^4+b^4=(t^4+6t^2w^2+w^4)/8.

The displayed cubic, quartic, distance and g expansions follow.
Every omitted tail term is bounded by delta^(3/2) or delta^2;
the analytic expansion errors are uniform for delta,w^2>=0 with
eta=delta+w^2 small. The denominator delta+w^2/2 is at least eta/2.
Consequently the quotient error is O(sqrt(eta)), in particular o(1),
uniformly as the ratio delta:w^2 varies.

The response coefficients are exactly

    Ldust=(28sqrt(2)-32)/3,
    Lsplit=(64-44sqrt(2))/3,
    Ldust-Lsplit=24sqrt(2)-32>0.

The comparison Lsplit>1/2 is equivalent to 125>88sqrt(2),
with positive sides and squared comparison 15625>15488.
The convex-combination form therefore proves the claimed strict gap
from K3. The formal boundary lists are used for expansion only and
are not misdeclared as satisfying p1=1.

## 3. Compactness, equality extraction and the third atom

I checked both closed envelope domains and their exact zero sets against
THM-4454. In the complementary domain the zeros are (z,z) and (1,0).
In the packing domain strict fixed-tail concavity and the two endpoint
factorizations leave precisely (z,z,z) and (h,h,0). The two expansions
above exclude the unwanted boundary limits for a minimizing sequence.

Since E and d2 are bounded, R->K3 implies F_actual->0. Every sequence
has a subsequence with convergent top two positive coordinates and a
further subsequence lying on one side of b=z. The envelope argument
therefore leaves (z,z) as the only possible top-two subsequential limit.
This is finite-dimensional compactness and does not assume norm
compactness of the varying lists.

For the third-atom step, the signed estimate

    sum_tail(r_i^3-Cr_i^4)<=m f_C(c)

is valid for positive, negative and zero tail entries. The derivative
of f_C on [0,sqrt(m)] has a fixed positive lower bound eventually.
The comparison with F3 gives exactly

    0<=epsilon*m*(sqrt(m)-c)<=F_actual-F3->0.

No global nonnegativity of this auxiliary F3 outside the packing
region is assumed; only its continuity and its value zero at the
three-atom point are needed. Thus c->z, and p2=1 recovers all lost
norm: the square mass outside the three selected positive coordinates
tends to zero.

The converse follows from dimension-independent continuity of p3,p4
and the top-two positive sum on the unit square-norm ball. I checked
the k-Lipschitz gradient estimate for p_k, k=3,4, and the sqrt(2)
bound for the top-two sum. The limiting energy is 1/3 and the limiting
d2 is strictly positive, so the quotient converges to K3. The exact
distance decomposition proves the remaining equivalence in the statement.

## 4. Negative-coordinate count and mixed-dust hostile

The total negative magnitude is at least max(a+b+c-1,0), while its
square mass is at most Delta3. Cauchy--Schwarz therefore gives the
displayed bound

    N_minus*Delta3>=
      max(sqrt(3)-1-(sqrt(3)/2)*Delta3,0)^2.

In particular the negative-coordinate count diverges for every minimizing
sequence. This argument permits any additional positive dust.

I checked the hostile construction for all n>=4: the defining quadratic
is strictly increasing for a>=0, is below one at a=1/2 and above one
at a=z. The stated rational endpoint estimate is already strict at n=4.
Its unique a_n lies in (1/2,z). Both dust magnitudes are smaller than
a_n, so the selected top-three coordinates are the actual leading
positive entries. The two moment normalizations hold exactly.

The dust square mass tends to zero and a_n tends to z, so the converse
applies. Meanwhile its positive dust sum is n and its negative magnitude
is n+3a_n-1. Both diverge while their signed difference tends to
1-sqrt(3). This is a valid actual hostile to stronger separate-sum claims.

## 5. Source audit and independent replay

I read the complete
[exact source](../../04-computation/long_frontier_sep06_minimizers.py).
Its quotient identities, parameterized boundary derivatives, radical
sign gates, and exact algebraic mixed-dust construction match the proof.
All gates use explicit exceptions and survive optimized Python. No
finite sample is used to infer the classification or a compactness claim.

The finite universe is exactly the seven actual algebraic lists indexed
by n=4,5,6,8,12,20,50, with no exclusions, together with symbolic
identities and formal boundary controls. The quadratic formula keeps
the positive root, all multiplicities and both dust signs.

Independent normal and optimized replays match the frozen output
byte-for-byte. All **91 always-active gates PASS**. Independently
computed raw file hashes and the output's semantic hash are:

```text
source SHA256 6d3be5780138c24f58312733270981d3375eae93b9281b26710529bf955dc7cd
output SHA256 19d6c54228a78db3b7231e30e9d6897bdc9a71ba904dd5159b9c861ab2040b20
semantic SHA256 306ae986875172f1ab4c69e9779d02cb0dc82e3c3cfdf4fdd144361797e9869a
```

No mathematical correction, scope defect or source defect was found.
