# Session110 checkpoint 1: verified prime structure and arithmetic observers

2026-10-30. **The overnight session is active.** This checkpoint preserves proved and independently checked milestones while the larger reproducibility packet and further branches are assembled. It does not claim completion of the overnight request.

## 1. Exact reconstruction from two prime-power aliquot probes

Write s(n)=sigma(n)-n. Fix a prime p and t>=1. Given a=s(n), b=s(p^t n), put

    e=v_p(b-a),
    R=1+p+...+p^(t-1),
    D=(b-p^t a)/R.

Then the source is uniquely

    n=(1+p+...+p^e)D-a.

There is no coprimality guard on n. To prove this, write n=p^e u with p not dividing u, and D=sigma(u). Prime-power multiplicativity gives

    b-p^t a=R D,
    b-a=p^e R [u+p s(u)].

Both R and the bracket are positive units modulo p. Thus the observed difference recovers the occupied exponent e, and the first identity recovers D. Substituting in sigma(n)-a gives n. This extends the coprime attachment identity in Session108; it is not an inverse of s(n) from one value alone. Any two distinct known positions along a prime-multiplication observation ray similarly determine its original source.

**Independent finite checks:** the producer checked12000 probes with a divisor sieve; a separate auditor checked3600 probes with direct divisor-pair enumeration. Both normal and optimized Python executions agree, with explicit checks retained under optimization.

### A genuine composite-multiplier boundary

The two distinct sources387152 and437510 have

    s(387152)=s(437510)=362986,
    s(3871520)=s(4375100)=5275324.

Their factorizations are2^4*24197 and2*5*67*653. Exact divisor sums, verified independently by trial division and divisor enumeration, are750138,800496,9146844,9650424 respectively.

For general m with at least two prime divisors, write n=h u, with h supported on primes dividing m and gcd(u,m)=1. Then

    b-ma=[sigma(mh)-m sigma(h)] sigma(u)=K_m(h) sigma(u).

If m=product p^t and h=product p^e, expanding the prime-power sums gives K_m(h)>=p^(e+t) for every p dividing m. Therefore the observed positive number D=b-ma bounds every occupied exponent. Enumerating those signatures and checking the remaining coprime core gives the complete fibre.

For the multiplier10 witness the only surviving signatures are h=16 and h=10, with K values68 and37. The complete fibre is exactly the two displayed sources. Finiteness alone is elementary from ordinary aliquot fibre bounds; the useful refinement is this explicit small candidate set and exact inverse test, not a new claim of finiteness or polynomial-time factorization.

## 2. Ramified tournament prime factors

Let mu^2-a mu-1=0, let p be an odd prime with v_p(a^2+4)=1, and let O=Z_p[mu] with uniformizer delta=2mu-a. Set R_n=O/(delta^n), g=mu^4. Consider labelled Cayley tournament connection sets S with S disjoint from -S, S union -S=R_n minus zero, and gS=S.

The inherited order-lifting argument gives

    A_(2k)=2^(p^k-1),
    A_(2k+1)=2^(((p+1)/2)p^k-1).

Here A_n counts labelled bank states on p^n vertices, not all tournament isomorphism classes.

Every composite bank state has a canonical partition into ideal cosets. Indeed a vertex-transitive tournament is strongly connected. Its maximal proper modules are disjoint and canonical. Translations make the block through zero an additive subgroup H; g preserves H. The identity

    g=(a^2+1)+a(a^2+2)mu

has mu coefficient congruent to -2a modulo p, a unit. Thus H is closed under mu and is an O-ideal. Conversely, an ideal whose nonzero cosets have constant orientation supplies a module partition.

The canonical outer quotient is prime. If P_n counts prime labelled states and A_0=1, unique quotient/kernel gluing gives

    A_n=sum_(r=1)^n P_r A_(n-r),
    P_n=A_n-sum_(r=1)^(n-1) P_r A_(n-r).

No additive splitting or section multiplicity is needed: the quotient and kernel connection sets directly define the source-labelled lift.

For the golden prime5:

    A_1=4, A_2=16, A_3=16384, A_4=16777216,
    P_1=4, P_2=0, P_3=16320, P_4=16646400.

All16 depth-two states are the same graph type T5[T5], where T5 is the cyclic five-vertex tournament. The full automorphism group has order5^6. Independent audit tested all4800 pair-module closures across the16 states and all120 permutations for each five-vertex factor.

At depth four,130560 composite states are not coherent under reduction by p. Thus failure of that single attachment test does not imply primality. The full ideal filtration is the necessary carrier. These counts are supported by a general proof, signed-cut constraints, and small exact controls; no exhaustive scan of all2^24 depth-four tournaments is claimed.

Classical modular decomposition and wreath-product automorphisms are credited. The scope extension is their application to the explicit ramified bank, not discovery of those general mechanisms.

## 3. Observed arithmetic codes can hide all true dependencies

For any finite binary matrix A on selected prime rows, CRT and Dirichlet produce distinct primes whose p+1 valuation parities project to A while their FULL squareclass columns are linearly independent. Add private witness-prime rows forming an identity matrix and prescribe exact valuations by congruences p=-1+l^e modulo l^(e+1). Each progression is coprime to its modulus, so Dirichlet provides infinitely many choices.

An explicit false-Fano bank is

    (7,2,41,19,109,1019,389),

observed only on prime rows2,3,5. The projected rank is3 and the full rank is7. There are16 projected square-syndrome selections but only the empty true square-sigma selection. The inherited genuine Fano bank is retained as a positive control. This concerns row projection, not arbitrary full prime-matroid realizability.

## Scope and continuation

The arithmetic and ramified branches have independent proof reviews and exact checks. Additional verified signed-residual work retains the original Collatz source law and negative credits; its certificate is finite and does not prove convergence. General-field enumeration, balanced spectral constructions and larger graph-isomorphism checks continue in this same overnight session.

No unrestricted A000568 census, Collatz convergence theorem, Wall–Sun–Sun classification, or global originality claim is made.
