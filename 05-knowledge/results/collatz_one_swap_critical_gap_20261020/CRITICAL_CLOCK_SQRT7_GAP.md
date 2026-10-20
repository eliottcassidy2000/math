# A square-root-seven content gap at the critical one-swap clock

Date: 20 October 2026. Status: **PROVED elementary refinement**, with an independent audit supplied alongside it. This is separate from the sharp all-length norm theorem. Its analytic consequences alone do not assert the final positive one-swap exclusion.

## 1. Exact setting and conclusion

Inherit the definitions, hypotheses, cut-safe parameterization, and proofs from [the sharp sparse-norm theorem](../collatz_one_swap_sparse_norms_20261020/SHARP_ONE_SWAP_CONTENT_BOUND.md). Thus p>=2, q>p, gcd(p,q)=1, and w is obtained by one cyclic adjacent swap of unequal bits from the zero-intercept lower mechanical word. Write

    D=2^q-3^p,    g=gcd(|D|,B(w)),    q=pk+r.

The predecessor proves g<=3^p-2, and any positive integer cycle must have q=ceil(p log_2 3).

**Theorem (critical-clock gap).** At every coprime critical shape q=ceil(p log_2 3),

    g < 7^(p/2).                                       (1.1)

**Necessary exponential condition.** Any positive integer cycle in this one-swap family therefore satisfies

    0 < q log 2-p log 3 < exp(-p/9).                   (1.2)

More directly, D/3^p<(7/9)^(p/2)<exp(-p/9). Thus no numerical logarithm computation is needed to test the stronger content bound g^2<7^p on a supplied word.

The key new coordinate is the edit's legality. At k=1, the two types having norms close to 3^p are not physical adjacent swaps. Moreover the entire allowed type range is bounded by r. That makes an elementary second-moment estimate sufficient; a variance improvement or linear-form theorem is not needed for (1.1).

## 2. Exact admissible-type ranges at k=1

Use the predecessor's residue labels

    a_j=(-jr mod p),    a=a_(t+1),

where the moved one is t+1 and the unchanged anchor is its predecessor t. Let A=a_t and a'=a_(t+2). The neighboring gap corrections are

    delta_(t+1)=(r+a-A)/p,
    delta_(t+2)=(r+a'-a)/p,

with A=(a+r mod p) and a'=(a-r mod p). For k=1, a move is legal exactly when its destination-side correction equals 1. Consequently

    epsilon=-1 is legal iff a>=p-r,
    epsilon=+1 is legal iff a<r.                       (2.1)

In particular a=0,- and a=p-1,+ are never legal at k=1. These are exactly the two nonfree endpoint types whose norms approach 3^p. The other endpoint types a=0,+ and a=p-1,- have norm modulus M=1.

This argument applies equally to the original cyclic crossing cut, because the residue model roots at the preceding unchanged one. It does not presume every formal signed type is an actual word.

## 3. The critical residue occupies at most three quarters of the period

At p=2 the critical length is q=4, not coprime to p, so there is no instance. For p>=3 the critical length satisfies p<q<2p, hence k=1. One may prove the strict upper bound without decimals: 3^p<2^(2p-1) holds at p=3 and is preserved by induction.

Also

    r/p <= 3/4.                                       (3.1)

For p>=7, the elementary inequality 3^5<2^8 gives

    r=ceil(p log_2 3)-p
      < (3/5)p+1 <= (3/4)p.

For p=3,4,5,6, the critical residues are respectively 2,3,3,4, which satisfy (3.1). The noncoprime p=6 shape may be retained in this arithmetic check without being counted as a word in the theorem.

## 4. The sparse norms are uniformly below 7^(p/2)

Let rho=2^(1/p) and z_j=rho exp(2 pi i j/p). The predecessor supplies the exact universal norm gates

    g divides M_+,    Q_+(Z)=Z^(a+1)-Z^a+1,
    g divides M_-,    Q_-(Z)=Z^(p-a)-Z+1,

where M is the positive 2,3-free part of Res(Z^p-2,Q). These gates retain prime-power multiplicities.

The nonfree endpoint types are forbidden by (2.1), and the free endpoints already have g=1. It remains to consider 1<=a<=p-2. In each polynomial the three exponents are distinct modulo p, so averaging over the p conjugates cancels all cross terms.

For a legal plus type, a+1<=r. Therefore

    (1/p) sum_j |Q_+(z_j)|^2
      =1+4^(a/p)+4^((a+1)/p)
      <1+2*4^(r/p)
      <=1+2*4^(3/4)
      =1+4 sqrt(2)<7.                                  (4.1)

The last strict inequality is equivalent to sqrt(2)<3/2.

For a legal minus type, p-a<=r, and p>=3. Thus

    (1/p) sum_j |Q_-(z_j)|^2
      =1+4^(1/p)+4^((p-a)/p)
      <=1+4^(1/3)+4^(3/4)<6<7,                         (4.2)

using 4^(1/3)<2 and 4^(3/4)=2 sqrt(2)<3.

AM-GM applied to the p nonnegative values |Q(z_j)|^2 gives

    |Res(Z^p-2,Q)| < 7^(p/2).

Removing powers of 2 and 3 cannot increase its positive modulus. As g divides that modulus, (1.1) follows. Nonvanishing of the resultant and its complete prime-power divisibility gate are supplied by the predecessor.

## 5. Exponential clock consequence and its limits

For a hypothetical positive integer cycle, D=g. Hence

    D/3^p < (sqrt(7)/3)^p = exp[-(p/2) log(9/7)].

The elementary integral inequality log(1+x)>x/(1+x), x>0, yields log(9/7)>2/9. Therefore D/3^p<exp(-p/9). Finally

    q log 2-p log 3 = log(1+D/3^p)<D/3^p,

proving (1.2).

An effective two-log lower bound and a rigorously certified continued-fraction computation can consume this exponential condition. Those are separate dependencies, not silently asserted here. The condition itself does not prove that all remaining small p are impossible.

The restriction to critical shapes matters for the square-root-seven constant. The all-length sharp family can have g=3^p-2, much larger than 7^(p/2). The legality lemma does apply to every coprime shape p<q<2p, but its quantitative residue bound r/p<=3/4 was obtained from the critical-clock condition.

## 6. Exact finite controls and the small-weight boundary

`verify_sqrt7_gap.py` uses no external mathematical libraries and no repository mathematical-code imports. It checks:

- every coprime critical shape with 2<=p<=512, using q=bit_length(3^p), an exact integer computation;
- every unequal cyclic adjacent swap occurrence, including the crossing cut;
- the complete type legality equivalences in (2.1);
- 4r<=3p, g^2<7^p, and g<D;
- separately, the full small-weight universe 2<=p<100, to make that finite exclusion available as an explicit boundary certificate.

The full receipt records 312 critical shapes, 93,878 swaps, and 312 original-cut crossings. The p<100 subreceipt is exact and exhaustive within the stated primitive-mechanical one-swap universe, not a census of arbitrary words. Combined with the predecessor's critical-length theorem, its g<D checks exclude every length for each weight 2<=p<100.

Reproduce:

    python3 verify_sqrt7_gap.py > sqrt7-receipt.replay.json
    python3 -O verify_sqrt7_gap.py > sqrt7-receipt.optimized.replay.json
    cmp sqrt7-receipt.json sqrt7-receipt.replay.json
    cmp sqrt7-receipt.json sqrt7-receipt.optimized.replay.json

The independent audit also checks the norm reduction and legality argument. These computations corroborate, but do not replace, the unbounded proof in Sections 2-5.
