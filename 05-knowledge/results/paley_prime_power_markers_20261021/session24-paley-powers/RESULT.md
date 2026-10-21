# Minimum signed-regular Paley markers over every prime power

Status: PROVED elementary certificate theorem conditional only on a CITED classical difference-preserver theorem; FINITE-EXACT implementation controls below. No literature-novelty claim.

## Statement

For every prime power q congruent to 3 modulo 4 other than 3 and 7, the minimum size of a subset with trivial setwise stabilizer in the full automorphism/anti-automorphism group of the Paley tournament P_q is exactly three. For q=3 or 7 no such subset exists.

More constructively, write q=p^e with e>1 (necessarily odd), and choose t of degree exactly e over F_p. If 3 does not divide e, every such t gives a regular set {0,1,t}. If 3 divides e, set r=p^(e/3). The set is regular if and only if

    t^r != 1/(1-t) and t^r != (t-1)/t.

There always exists a full-degree t satisfying both inequalities. Thus this is a finite, deterministic certificate search: enumerate field elements and test degree and two equations, rather than enumerate subsets or group elements.

## Inheritance and limits

The exception classification was already proved in THM-4743, `01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md`. Its G2 control already supplies a regular subset for P27. The prime-only predecessor `session24-paley-certificates/RESULT.md` proves the fixed certificate {0,1,3} at every prime p>=11 and its affine triple count. This note adds the all-prime-power minimum-size marker argument, not the exception classification itself. Bounded searches of canon and results for Carlitz, full-degree Paley markers, and signed-regular triples found no existing version of this explicit criterion; no exhaustive novelty search is asserted.

Do not transfer the affine regular-orbit counts to extension fields. Frobenius is a genuine extra automorphism, and it fixes every prime-subfield point. In characteristic 3 the expression {0,1,3} is not even a three-set. The valid hostile {0,1,2} is fixed pointwise by Frobenius in every extension considered here. The repaired ordinary/signed distinction in THM-309 also remains necessary: nonsquare slopes reverse orientation.

Concept board: affine kernel; field-generator degree; faithful three-point action; cyclic Galois image; Frobenius–Möbius obstructions. The kernel and faithfulness sidecars are separate obligations. No centroid division by 3 is used. No claim about Hamiltonian-path gaps follows.

## 1. The full group, with primary source checked

A. Bruen and B. Levinger, *A theorem on permutations of a finite field*, Canadian Journal of Mathematics 25 (1973), 1060–1065, Theorem 1, DOI https://doi.org/10.4153/CJM-1973-113-4. Publisher PDF:
https://www.cambridge.org/core/services/aop-cambridge-core/content/view/916D5CE868A79BA0171F8F1CACBA4C66/S0008414X00051154a.pdf/div-class-title-a-theorem-on-permutations-of-a-finite-field-div.pdf

This paper proves the Carlitz–McConnel theorem directly. Theorem 1, page 1060, states that for a multiplicative subgroup K of proper order d dividing q-1, a map f has every nonzero difference ratio (f(x)-f(y))/(x-y) in K exactly when f(x)=b+a*x^(p^j), where a is in K and q-1 divides d(p^j-1). No exception is stated. The original Carlitz article is identified by its reference 2: *A theorem on permutations in a finite field*, Proc. Amer. Math. Soc. 11 (1960), 456–459. We rely on the directly checked 1973 primary proof, not an uninspected rendering of the 1960 paper.

Take d=(q-1)/2 and K the squares. An ordinary Paley automorphism preserves the quadratic character of every difference, so its difference ratios lie in K. The divisibility condition is automatic because p is odd. Conversely square-slope affine-semilinear maps preserve the character. Since -1 is a nonsquare, negation reverses all arcs, and composing it with any anti-automorphism gives an ordinary automorphism. Therefore

    Aut±(P_q) = AΓL(1,q) = {x -> a*x^(p^j)+b: a!=0, b in F_q, 0<=j<e}.

This identification is in the actual field coordinates, not just an abstract group isomorphism. In particular the order is e*q*(q-1).

## 2. Direct affine audit in every odd characteristic

For S={0,1,t}, t not 0 or 1, an affine map is determined by its images of 0 and 1. Checking the five nonidentity permutations of S gives the following necessary and sufficient conditions:

- Swapping 0 and 1 and fixing t: t=1/2.
- Swapping 0 and t and fixing 1: t=2.
- Swapping 1 and t and fixing 0: t=-1.
- Either three-cycle: t^2-t+1=0.

These statements remain valid when characteristic 3 causes conditions to coincide. Every exceptional parameter has degree at most two. Therefore a parameter of odd degree e>=3 has trivial affine setwise stabilizer. This proof avoids the invalid characteristic-three centroid argument.

## 3. Two faithful maps, neither omitted

Let H be the full affine-semilinear setwise stabilizer of S for such a full-degree t.

The projection AΓL(1,q) -> Gal(F_q/F_p) has affine kernel. Its restriction to H is injective by Section 2. Consequently H is cyclic and its order divides e.

Separately the action H -> Sym(S)=S3 is faithful. Indeed a pointwise fixer first fixes 0, forcing b=0, then fixes 1, forcing a=1, and finally fixes t. Since t generates F_q over F_p, the remaining field automorphism is identity. Hence |H| also divides 6. Since e is odd, |H| is 1 or 3, with 3 possible only when 3 divides e. This separately proves that no hidden pointwise semilinear fixer or order-nine lift is present.

## 4. Exact obstruction and counting

Suppose 3 divides e and put r=p^(e/3). A nontrivial H has order three and its Galois image is the unique order-three subgroup. One of its two nonidentity elements projects to x -> x^r (invert the other if necessary). Its faithful action on S is a three-cycle.

For the cycle 0 -> 1 -> t -> 0, the element must be

    f(x)=(t-1)*x^r+1,

and the last image condition is precisely t^r=1/(1-t). For the reverse cycle 0 -> t -> 1 -> 0, it must be

    f(x)=-t*x^r+t,

and the condition is t^r=(t-1)/t. Conversely either equation produces a nonidentity semilinear set stabilizer. This proves the iff criterion, including both possible Frobenius powers via inversion.

Clearing the nonzero denominators gives two nonzero polynomials of degree r+1. Thus at most 2(r+1) field elements satisfy an obstruction. Every proper subfield has degree a proper divisor d of e, hence d<=e/3. The union of proper subfields therefore has cardinality at most

    sum_{d=1}^{e/3} p^d < p*r/(p-1) <= 3*r/2.

Here the sum over all integer d is an upper bound; it does not assert that a subfield exists for every d. Since q=r^3 and r>=3,

    q - 3*r/2 - 2*(r+1) = r^3 - 7*r/2 - 2 > 0.

So some full-degree t avoids both obstructions, even at q=27. When 3 does not divide e, full-degree elements exist (for example, a generator of the multiplicative group), and Section 3 already proves every one works.

## 5. Minimum and exceptions

The prime-field predecessor proves {0,1,3} regular for each prime p>=11 congruent to 3 modulo 4, and proves the failure of all subsets at p=3,7. Together with Sections 1–4 this supplies a triple for every other admissible q.

No subset of size below three is regular: the empty set is fixed by the nontrivial group; a singleton {u} is stabilized by x -> 2u-x; a pair {u,v} is stabilized by x -> u+v-x. These reflections are nonidentity anti-automorphisms. Hence three is the exact minimum. The prime exceptions may also be checked without enumeration: at q=3 complementation reduces to singletons; at q=7 it reduces to sizes at most three, and the five parameters t=2,3,4,5,6 are exactly the three midpoint exceptions 2,4,6 and the two quadratic roots 3,5.

## 6. Exact full-semilinear triple enumeration

This formula counts the full signed group, not just its affine subgroup. Let e>1 be odd and define

    G_e = sum_{d|e} mu(e/d)*p^d.

This is the number of full-degree field elements, by Möbius inversion over subfields. A normalized regular triple must have full-degree parameter: if t lies in a proper subfield, a nonidentity Frobenius fixes 0,1,t pointwise. Conversely Sections 3–4 completely characterize the obstruction for full-degree parameters.

If 3 does not divide e, the number N_e of regular normalized parameters is G_e. If e=3u, put D=1 when p=3, D=2 when p is 1 modulo 3, and D=0 when p is 2 modulo 3. For each divisor d of e let g=gcd(u,d) and define

    M_d = D                 if d/g=1,
    M_d = p^g+1             if d/g=3,
    L_e = sum_{d|e} mu(e/d)*M_d.

These are the only two possibilities for d/g. Then

    N_e = G_e - 2*L_e,
    number of regular labelled three-sets = q*(q-1)*N_e/6,
    number of free signed three-set orbits = N_e/(6*e).

Proof of M_d: On F_(p^d), the automorphism sigma:x->x^(p^u) has order d/g. If it has order one, sigma(t)=1/(1-t) reduces to t^2-t+1=0. Since d is odd, this has exactly D roots. If the order is three, d=3g; sigma is either x->x^(p^g) or its inverse. Put h=p^g and rho(t)=1/(1-t), so rho has order three. In the first case the cleared equation is (1-X)*X^h-1=0. It has degree h+1 and derivative -X^h, nonzero at every root. Frobenius commutes with rho, so iterating three times shows every root satisfies X^(h^3)=X. Thus all h+1 distinct roots lie in F_(h^3). For inverse sigma, applying its inverse converts the equation into X^h=rho^(-1)(X)=(X-1)/X. The polynomial X^(h+1)-X+1 likewise has nonzero derivative at its roots (a derivative zero would force X^h=1, incompatible with the equation) and all its h+1 distinct roots lie in the same field. Therefore M_d=h+1 in both cases.

Möbius inversion now counts full-degree roots of one obstruction as L_e. The two obstruction sets intersect only where rho(t)=rho^(-1)(t), which is t^2-t+1=0; such parameters have degree at most two and cannot have degree e. Hence precisely 2*L_e full-degree parameters are bad. Every unordered triple has six ordered first-two-point normalizations, yielding the labelled count. Dividing by the full signed group order e*q*(q-1) is legitimate for the regular triples and gives the orbit formula.

For e=3 the orbit count simplifies to

    (p^3 - 3*p - 2 + 2*D)/18.

The exact controls give respectively one, eight, and eighteen free signed triple orbits at q=27,243,343. The formula is a proved consequence of the same marker analysis; finite counts are checks, not its justification.

## Exact controls and reproduction

`verify.py` uses only Python's standard library. It builds fields by polynomial quotient, certifies each defining polynomial by exhaustive trial monic divisors up to half degree, and verifies Frobenius and inverses. For every full-degree parameter in q=27,243,343 it enumerates the six normalized image pairs for every Frobenius power, checks the kernel, stabilizer orders, exact obstruction equivalence, and the Möbius count. A second exhaustive path scans every a!=0, b, and j for one regular witness and one hostile in each field. The prime-subfield pointwise-Frobenius hostile is explicitly checked. This is an exhaustive finite group check, not independent computational discovery of the full tournament automorphism group: that identification uses the cited theorem.

    python3 session24-paley-powers/verify.py > session24-paley-powers/normal.json
    python3 -O session24-paley-powers/verify.py > session24-paley-powers/optimized.json
    cmp session24-paley-powers/normal.json session24-paley-powers/optimized.json

All checks use RuntimeError rather than removable assert statements. Output records the basis convention, irreducible polynomial, all bad full-degree parameters, a witness, and hostile stabilizers. `receipt.json` hashes sources and outputs. No repository file or Git state was modified by this worker.
