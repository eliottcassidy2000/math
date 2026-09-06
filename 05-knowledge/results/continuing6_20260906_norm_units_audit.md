# Independent referee: prime-divisor complementary norm units

**PASS: the all-parameter prime-divisor theorem, its strengthened valuation threshold, and the carried-boundary consequence are accepted without a mathematical repair.** The audit also independently verifies the entire producer certificate, actual rational norms on a separate small rectangle, and a genuine modular common-root hostile whose rational rows remain coprime.

The reviewed candidate is [continuing6_20260906_norm_units.md](continuing6_20260906_norm_units.md), with source `continuing6_20260906_norm_units.py`. No producer implementation is imported. This audit was conducted outside the repository; parent owns filing and promotion.

## 1. Exact statement accepted and normalization

For integers H>=1,z>=1, put S=z+2H and

`Phi(t)=sum_(j=0)^H H! S! t^j/[j!(z+2j)!(3H-3j)!]`,

`Psi(t)=Phi_(2H,2z)(t)`.

Both are monic. In the H-dimensional rational quotient by Phi, define N_H(z) as the determinant of multiplication by -Psi. Thus

`N_H(z)=(-1)^H Res(Phi,Psi)`.

This identity does not require Phi to be squarefree or to have real roots. The accepted sufficient condition is: some prime p divides2S-1, p>4H, and e=v_p(2S-1)>floor(6H/p). It implies

`v_p([t^j]Phi)=0` for every j, and

`v_p([t^j]Psi)=e-floor((6H-3j)/p)>0` for j<2H.

Consequently both rows are p-integral, Psi reduces to t^(2H), and

`N_H(z) modp=(-1)^H Phi(0)^(2H) !=0`.

The conclusion is rational coprimality and hence absence of a common complex root. It does not determine the real norm sign. The theorem holds for every positive integer z, including even z; the producer's odd-z main bank is a finite validation choice, not a theorem restriction.

## 2. Independent analytic derivation

A convenient independent coefficient construction starts from c_H=1 and lowers by

`c_(j-1)/c_j=j(z+2j)(z+2j-1)/[(3(H-j)+1)(3(H-j)+2)(3(H-j)+3)]`.

Multiplying these ratios gives the same rational row and also the finite-product expression

`c_j=prod_(k=j+1)^H k * prod_(k=z+2j+1)^S k / prod_(k=1)^(3H-3j) k`.

This preserves all numerator and denominator factors without replacing rational factorial ratios by integers prematurely.

Because p is an odd divisor of2S-1, the standard residue of S is(p+1)/2. The last2(H-j) factors ending at S lie entirely among nonzero residues: their length is at most2H<(p+1)/2. Also3H<p, so all remaining numerator and denominator factorials of Phi are units. The constant term is included.

For j<2H, the doubled upper interval has length4H-2j, between2 and4H<p. It contains2S-1 and no other multiple of p. Its valuation is exactly e, even when S is much larger than p. The factorials(2H)! and j! are units. The remaining factorial index6H-3j is less than3p/2 and hence less than p^2; its valuation is precisely floor((6H-3j)/p). This proves the exact doubled formula. No neglected second power of p can enter that denominator.

The strict threshold makes all these nonleading valuations positive. Reduction of a monic resultant is legitimate over the local ring Z_(p), and

`Res(Phi,t^(2H))=prod_(Phi(rho)=0) rho^(2H)=Phi(0)^(2H)`.

The even exponent removes the constant-term product sign; the factor(-1)^H remains from the definition of N. This proves the actual norm statement, rather than only a coefficient statistic.

Since p>4H, floor(6H/p) is0 or1. Thus the stated simpler branches p>6H with one endpoint factor, or p>4H with at least a square endpoint factor, are valid. The second branch is not silently promoted to all endpoint primes above4H.

## 3. Carried-boundary consequence and its limits

The inherited [continuing5_20260906_complementary_norm_jets.md](continuing5_20260906_complementary_norm_jets.md) supplies the exact leading-jet formula and the equivalence between order r-1 at x=-r and N_(h-r)(2r+1)!=0. Its inverse-carry treatment and normalization remain load-bearing. The present theorem changes only the complementary nonvanishing input.

With H=h-r,z=2r+1, one has S=2h+1 and2S-1=4h+1. Thus the claimed condition

`p|4h+1, p>4(h-r), v_p(4h+1)>floor(6(h-r)/p)`

has the correct scale. H=0 is separately covered by N_0=1. When4h+1 is prime, the simple inequality p>6(h-r) is equivalent, for integral r, to r>=ceil(h/3). The complementary heights then reach floor(2h/3).

The elementary proof of infinitely many primes1 modulo4 is valid: for any finite list, a prime divisor of(2 times their product)^2+1 is outside it, odd, and has an element of order4 in its multiplicative group. Hence it is1 modulo4. This establishes unbounded complementary heights analytically; no finite prime bank is used as the infinitude argument.

The accepted claim is exact boundary multiplicity. It is not a sign theorem, a complete root classification, or noncancellation at arbitrary interior phases. The old H<=6 sign certificates are not weakened or replaced.

## 4. Independent exact paths

The referee builds coefficients by the lowering recurrence, checks them against finite falling products, and computes modular polynomial gcds and resultants by Euclidean division. This differs from the producer's direct factorial coefficients and quotient multiplication matrix. It compares every complete producer certificate row, including its norm residue.

- The independently generated universe contains all492 eligible(H,z,p) triples with H1..24, odd z1..101 and all possible endpoint prime divisors. The endpoint is at most297, so the prime sieve through301 is complete. There are292 cases with H>=7 and four square-prime extension cases.
- All16 periodic controls are replayed, including z larger than p.
- Six additional exact controls include even z, higher valuation, and parameters outside the producer bank.
- For every H1..4,z1..12, a literal rational Sylvester matrix gives the actual norm. These48 norms agree with the independent Euclidean resultant modulo101 and103, and every eligible certificate prime in that rectangle has norm valuation zero. This also validates the resultant sign convention independently of the monomial case.
- The prime-boundary ceiling equivalence is checked exactly for all67 prime endpoints with h1..200. This is a control for the elementary integer inequality, not the all-height proof.

All37331 gates are active under optimization. Normal and -O outputs are byte-identical actual LF.

## 5. Two distinct boundary hostiles

For H=z=1,p5, the literal rows are

`Phi=t+1, Psi=t^2+10t+1`.

The endpoint valuation e1 equals floor(6H/p), so the strict hypothesis fails. Psi mod5=t^2+1, not a monomial. Its rational response at the Phi root is-8, while the characteristic norm is+8. These are different quantities and their signs are retained correctly.

The larger hostile is H26,z3,p109. Here p>4H but e1=floor(6H/p), again exactly the excluded boundary. The original rational rows are p-integral; their reduced gcd is

`t+104=t-5 in F_109[t]`.

Both literal rows vanish at t5 modulo109, and their resultant is zero there. Hence deleting the valuation threshold can genuinely destroy the modular-unit conclusion, not merely its proposed proof.

However, reducing the same rational rows modulo157 gives gcd1 and

`Res(Phi,Psi) mod157=N_26(3) mod157=65`.

Monicity and integrality at157 therefore prove their rational gcd is1. This hostile is not a counterexample to rational noncancellation. As a **FINITE-EXACT extra consequence**, the inherited boundary iff gives order0 at(h,r)=(27,1), even though the only endpoint prime109 does not satisfy the new arithmetic certificate. This single point demonstrates that the sufficient condition is not necessary; no broader all-prime or all-height conclusion is inferred.

## 6. Reproduction and frozen pins

```
python continuing6_20260906_norm_units_audit.py
python -O continuing6_20260906_norm_units_audit.py
```

Filed layout prefers the sibling results certificate; the outside fallback uses `C:/w/continuing6_20260906_norm_units_certificate.json`. The certificate is read only for comparison, never imported as a mathematical oracle. The source contains all independent constructions above.

Reviewed producer pins:

- Source:`b63551e61022fb2058d0a65c2b07115346128490ae5fb43e0e69b1f46200ceb7`.
- Candidate report:`24ca90c72183b60f4c59f164e02bb27c499af4789921d533c8572c84973bc985`.
- Certificate:`058354c40113b27150af49e164d8537b6f75eb6352dcc9fe9020229b47184479`.

Independent referee pins:

- Source:`f483a3cca43854059ab5b403c22372bab392bffe8643171b7ed328a2b37a39ea`.
- Normal and optimized output:`58d2f070f2e6436be1ad0dc3c5c3fbfa7ba5fe3137a1ddfd2f04c04bede17053`.

Promotion basis: the full coefficient-valuation proof, actual resultant transfer, boundary normalization, infinitude argument, scope boundaries, and independent exact controls all pass. No repair is requested. Parent may promote the producer to independently audited proof status while retaining this reviewed candidate hash as lineage.
