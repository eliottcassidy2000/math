# Exact arithmetic reduction for a finite O/Z tail bank

Status: independently audited elementary all-order reduction, not an occupancy proof or a general finite-bank impossibility proof. No finite physical bank is supplied that covers a cofinite set.

For a fixed reachable edge-grammar state s=(h0,R0,e0,d), append any O/Z word. Encode its chronological letters as binary digits b_1,...,b_n from least significant to most significant, with O=0,Z=1. Thus m=sum_i b_i 2^(i-1). Terminal O letters can be deleted for scalar output, so the nonzero encoded integers enumerate all nonconstant scalar tails, without a guessed length bound.

Define Q(0)=0 and

    Q(2m)=Q(m),
    Q(2m+1)=Q(m)+m+(m mod 2).

Then the scalar output is exactly

    H_s(m)=h0+R0*m-d*Q(m)-e0*(m mod 2),    m>=0.

This is an identity about actual O/Z descendants of the specified source. A hypothetical nonnegative state cannot be substituted for a physical seed without proving its attainability.

## Proof

Let b_i indicate whether step i is Z. At each step R_i=2R_(i-1)-d*b_i, while e_i=d*b_i. Hence

    R_(i-1)=2^(i-1)*R0-d*sum_(j<i) b_j*2^(i-1-j).

The increment in h at step i is b_i*(R_(i-1)-e_(i-1)). Summing gives the formula with

    Q(m)=sum_(1<=j<i<=n) b_i*b_j*2^(i-1-j)
         +sum_(i=2..n) b_i*b_(i-1).

Deleting a zero low bit does not change this pair sum, yielding Q(2m)=Q(m). Inserting a one low bit contributes m from the first sum and the low bit of m from the adjacent-pair sum, yielding the odd recurrence. The initial e0 contributes only at step1.

Every Q(m) is even: adjacent pairs contribute 1+1, and all longer-distance pairs have an even coefficient. Induction also gives 0<=Q(m)<=m-1 for m>=1. Therefore, for m>=1,

    H_s(m)>=h0+(R0-d)*m+d-e0>=h0+(R0-d)*m.

At m=0, H_s(0)=h0, so the weaker bound H_s(m)>=h0+(R0-d)*m holds for every m>=0. Since reachable unit-start seeds have R0>d, this gives a complete finite inverse-search bound for each target, m<=(H-h0)/(R0-d). It does not imply the target has a preimage.

## Precise remaining covering obligation

For a proposed fixed finite bank B of physically verified states, eventual odd coverage by that bank's O/Z tails is equivalent to

    for every sufficiently large odd N, there exist s in B and m>=0
    such that N=h_s+R_s*m-d_s*Q(m)-e_s*(m mod2).

This equivalence concerns only the O/Z tails of that fixed bank. Such coverage would suffice for eventual coverage of the unrestricted tournament spectrum; it is not known to be necessary for that global conjecture. The fixed-bank equivalence is exact and all-order. It is not an inference from dense finite image computations. The identity JZ=OJ allows response-preserving word normalization, but does not by itself reduce every source to an O/Z tail of one fixed finite bank: the oriented rewrite JZ -> OJ gives no reduction of the subword JO, and this identity alone supplies no bound on the possible prefix states.

The root-only specialization is H(m)=1+2m-Q(m). Its checked finite image contains249077 odd values<=999999 and24904080 odd values<=99999999. These are exploratory finite arithmetic computations, not a proof of an asymptotic density, noncofiniteness, or any finite-bank impossibility. In particular, approximate half-density for one bank member says nothing decisive about a finite union of different images.

## A necessary reciprocal-width condition on a proposed bank

For each seed put g_s=R_s-d_s>0. The inverse bound shows that its scalar image contains at most X/g_s+O(1) distinct values at most X. Consequently a finite bank satisfying

    sum_(s in B) 1/g_s < 1/2

cannot cover all sufficiently large odd integers: its union has at most X*sum_s(1/g_s)+O(1) values, whereas a cofinite odd set has X/2+O(1). Repetitions only strengthen the upper bound. This excludes such banks, not every finite bank. In particular it gives no obstruction to a bank containing the rooted edge, for which g=1.
