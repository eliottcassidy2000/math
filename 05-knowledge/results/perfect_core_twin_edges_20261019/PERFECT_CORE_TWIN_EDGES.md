# Perfect cores turn twin-prime pairs into actual aliquot edges

19 October 2026. Elementary consequences of the repository's guarded prime-extension formula. No historical novelty or prime-infinitude claim is made.

The inherited identity from [the guarded prime-extension note](../aliquot_guarded_cycles_20261014e.md) is, for p prime and p not dividing C,

    s(Cp)=(sigma(C)-C)p+sigma(C).

The prime and coprimality conditions are essential. The formula is not an iteration rule after the marked cofactor ceases to satisfy them.

## 1. Perfectness is exactly the translation case

**Theorem.** A fixed positive integer core C is perfect if and only if its guarded prime-extension rule is the translation

    Cp -> C(p+2).

More strongly, if there are two distinct primes p_1,p_2 not dividing C and a common number h such that

    s(Cp_i)=C(p_i+h),      i=1,2,

then sigma(C)=2C and h=2.

Proof: subtract the two affine identities. Their nonzero difference gives sigma(C)-2C=0. Substitute back to obtain Ch=sigma(C)=2C. The converse is immediate.

This gives an exact meaning to the deficient/perfect/abundant core triple. The cofactor branch has slope sigma(C)/C-1, respectively below, equal to, or above one. Perfect cores are precisely the neutral branches, with their positive carry fixing the translation at +2.

## 2. A literal twin-prime/aliquot bijection

Take the perfect core C=6. For every prime p>=5,

    s(6p)=6(p+2).

Therefore

    (p,p+2) both prime, p>=5
      <->
    an actual aliquot edge 6p -> 6q with p,q>=5 prime.

The correspondence is bijective, with q=p+2. In particular, the twin-prime conjecture is equivalent to infinitely many such edges in the aliquot graph. At a finite cutoff X, the number of these edges whose source is at most X is exactly the number of twin-prime pairs with lower prime p>=5 and p<=X/6.

This is an exact reformulation, not a prime-counting estimate. It does not follow from global cofactor convergence or supply a proof of twin-prime infinitude.

For any other fixed perfect C, the same equivalence holds after excluding the finitely many prime endpoints dividing C. Thus the finite exception does not change the infinitude question. No odd perfect number is assumed to exist; C=6 already supplies the full statement.

## 3. Every fixed-6 prime lane is short, independently of twin-prime infinitude

Although there might be infinitely many such edges across different starts, every run restricted to numbers 6p with prime p>=5 leaves that class after at most two aliquot steps.

Indeed p is 1 or 2 modulo 3. If p≡1, then p+2 is a multiple of 3 greater than 3 and is composite, so the first step leaves the class. If p≡2, the first successor can be prime, but then p+4 is a multiple of 3 greater than 3, so the second step leaves. Both affine steps are used only while their input cofactors are prime and coprime to 6.

Thus infinitely many admissible two-vertex segments would coexist with a uniform bound on the length of each single prime-guarded run. Quantifiers over different sources and quantifiers along one trajectory are distinct. Nothing here describes the subsequent aliquot dynamics after leaving this class.

A numerical agreement outside the domain does not restore the prime guard. For example s(54)=66 happens to equal 6(9+2), but 9 is composite, so this is not one of the prime-to-prime edges in the bijection. Conversely s(150)=222 differs from 6(25+2)=162. The compiler must retain the source predicate even when an unguarded example happens to match its affine formula.

## 4. Relation to this session's descent work

The dyadic-only class C=2^a has deficient cores and contracting guarded affine slopes; its associated [auxiliary cofactor-map dynamics were globally classified](../dyadic_cofactor_global_descent_20261019/COFACTOR_MAP_GLOBAL_CLASSIFICATION.md), with the aliquot interpretation retained only on the odd-prime domain. Adding the perfect core 6 produces a neutral +2 branch whose legal transitions are controlled by twin primality. An abundant fixed core produces an expanding branch. This is a direct arithmetic connection between the three slope regimes and the original deficient/perfect/abundant triple, with the native prime domain retained throughout.

The potential proof target supplied by this bridge is a lower bound for the number of neutral-core prime-to-prime edges. A theorem bounding the length of one such lane does not establish that lower bound. Any proposed transfer back to twin primes has to estimate the number of distinct legal sources, rather than merely classify their finite trajectories.
