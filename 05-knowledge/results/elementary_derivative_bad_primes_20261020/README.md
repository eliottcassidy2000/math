# Elementary bad primes for every nontrivial integer derivative

2026-10-20. **PROVED**, with standard-library exact controls and a completed [independent mathematical audit](audit/REVIEW.md). The [concise replacement proof](THEOREM_W_REPLACEMENT.md) supplements the incoming theorem without silently rewriting its source.

## Result

For every integer \(c\in\mathbb Z\setminus\{0,1\}\), infinitely many odd primes \(p\nmid c\) satisfy
\[
v_2(\operatorname{ord}_p(c))=1,
\qquad\text{equivalently}\qquad
c\bmod p\in-K_{\max}(p),
\]
where \(K_{\max}(p)\) is the unique odd-order subgroup of \(\mathbb F_p^\times\) of order equal to the odd part of \(p-1\).

For \(|c|\ge2\), more precisely: **for every odd prime \(\ell\) not dividing \(c+1\), every prime divisor \(p\) of**
\[
N_\ell(c)=\frac{c^\ell+1}{c+1}
=c^{\ell-1}-c^{\ell-2}+\cdots-c+1
\]
**has \(\operatorname{ord}_p(c)=2\ell\).** In particular \(p\equiv1\pmod{2\ell}\), \(p\ne\ell\), and different values of \(\ell\) give disjoint sets of prime divisors. Taking the odd primes \(\ell>|c|+1\) is an immediately transparent infinite family; the weaker condition \(\ell\nmid c+1\) is more economical for finite controls.

This proves the infinitude assertion in incoming THM-4724, Theorem W(ii), without Zsigmondy, Dirichlet, quadratic reciprocity, Kummer theory, or Chebotarev. It does **not** establish a density or prescribe the exact value of \(v_2(p-1)\).

## Elementary proof

For \(|c|\ge2\), let \(\ell\) be an odd prime with \(\ell\nmid c+1\).

1. \(N_\ell(c)\) is a positive odd integer greater than one. For \(c=a\ge2\), it is \((a^\ell+1)/(a+1)>1\); for \(c=-a\le-2\), it is \((a^\ell-1)/(a-1)=1+a+\cdots+a^{\ell-1}>1\). Modulo 2 the alternating polynomial is 1 when \(c\) is even and is a sum of the odd number \(\ell\) of ones when \(c\) is odd.
2. Evaluation at \(c=-1\) gives
   \[
   \gcd(N_\ell(c),c+1)=\gcd(\ell,c+1)=1.
   \]
   Also \(N_\ell(c)\equiv1\pmod c\), so its prime divisors do not divide \(c\).
3. Choose any prime \(p\mid N_\ell(c)\). It is odd, \(p\nmid c(c+1)\), and \(c^\ell\equiv-1\pmod p\). Thus \(\operatorname{ord}_p(c)\) divides \(2\ell\) but not \(\ell\). Since \(\ell\) is prime, this order is either 2 or \(2\ell\). Order 2 would give \(c\equiv-1\pmod p\), contrary to \(p\nmid c+1\). Hence the order is exactly \(2\ell\).
4. An order \(2\ell\) element of \(\mathbb F_p^\times\) forces \(2\ell\mid p-1\); therefore \(p\ne\ell\). There is also a direct check: Fermat's congruence and \(\ell\nmid c+1\) give \(N_\ell(c)\equiv1\pmod\ell\).
5. There are infinitely many allowable odd primes \(\ell\), since only finitely many divide the fixed nonzero integer \(c+1\). A prime \(p\) cannot arise from two different allowable \(\ell\), because \(\operatorname{ord}_p(c)\) is uniquely determined. This proves infinitude.

For \(c=-1\), every odd prime works, since the order is 2. For \(c=1\), the order is always 1, so no odd prime works. The integer \(c=0\) is excluded: it is not a nonzero dilation, and no prime satisfies \(p\nmid0\).

To check the coset formulation, write \(p-1=2^e m\) with \(m\) odd. If \(x=-k\) with \(k\in K_{\max}(p)\) of odd order \(d\), then \(x\) has order \(2d\). Conversely, if \(x\) has order \(2d\) with \(d\) odd, its order-two power \(x^d\) is the unique involution \(-1\), so \((-x)^d=1\) and \(-x\in K_{\max}(p)\). Thus membership in \(-K_{\max}(p)\) is exactly the condition that the multiplicative order has 2-adic valuation one. Equivalently, \((-x)^m=1\).

## Prime-power corollary

For every odd prime \(p\nmid c\) and every integer \(r\ge1\),
\[
v_2(\operatorname{ord}_{p^r}(c))
=v_2(\operatorname{ord}_p(c)).
\]
Indeed, reduction of unit groups has a kernel of order \(p^{r-1}\). If \(d=\operatorname{ord}_p(c)\) and \(d_r=\operatorname{ord}_{p^r}(c)\), then \(d\mid d_r\mid d\,p^{r-1}\), so \(d_r/d\) is odd. Consequently every bad prime constructed above is the base of an entire bad prime-power tower:
\[
c\bmod p^r\in-K_{\max}(p^r)\quad(r\ge1),
\]
where \(K_{\max}(p^r)\) is the unique maximal odd-order subgroup of \((\mathbb Z/p^r\mathbb Z)^\times\). The order is \(2\ell p^j\) for some \(0\le j\le r-1\), not necessarily \(2\ell\).

One can define \(K_{\max}(p^r)\) directly as the inverse image of \(K_{\max}(p)\) under reduction: its order is the odd integer \(p^{r-1}m\), and every odd-order subgroup of the unit group lies in it.

This transports the complete-graph failure to **ring circulants on \(\mathbb Z/p^r\mathbb Z\)**, not circulants on the additive group of the finite field \(\mathbb F_{p^r}\). Write \(H_r=K_{\max}(p^r)\). For every nonzero residue \(x\), including nonunits, no \(h\in H_r\) satisfies \(hx=-x\): if it did, iterating through the odd order \(d\) of \(h\) would give \(x=h^dx=(-1)^dx=-x\), which is impossible in an odd-order additive group. Thus the \(H_r\)-orbits on all nonzero residues come in disjoint antipodal pairs. Choose one orbit from each pair to obtain an \(H_r\)-invariant tournament connection set \(S\). If \(c=-k\) with \(k\in H_r\), then \(1_S(cx)=1_S(-x)=1-1_S(x)\) for every nonzero residue, so the derivative is the full nonzero connection set and represents the sign-odd complete graph \(K_{p^r}\). No prime-power type-bijection theorem is asserted or needed here.

## The two distinct consequences for circulant rules

For an odd prime \(p\nmid c\), define
\[
D_c(S)=\{x\in\mathbb F_p^\times:
1_S(x)\mathbin{\mathrm{XOR}}1_S(cx)=1\}.
\]
A tournament connection set satisfies \(S\sqcup(-S)=\mathbb F_p^\times\).

**1. Infinitely many actual odd-valued outputs.** If \(S\) is \(K_{\max}(p)\)-invariant and \(c\in-K_{\max}(p)\), then \(1_S(cx)=1-1_S(x)\), so \(D_c(S)=\mathbb F_p^\times\). Such invariant tournament sets exist: \(K_{\max}(p)\) is odd-order, so one can choose one coset from every antipodal coset pair. Every such \(S\) maps to the complete graph \(K_p\), which is sign-odd, not a sign-even graph. A vertex transposition has edge-reversal sign \(-1\). This is **not** a degree-parity assertion: \(K_p\) is Eulerian for odd \(p\), since its vertex degrees are \(p-1\). Therefore every fixed \(c\in\mathbb Z\setminus\{0,1\}\) fails sign-even-valuedness at infinitely many primes.

**2. Failure of a uniform derivative type bijection was already independent.** The published, independently audited
[Boolean-affine obstruction](../affine_circulant_ideal_obstruction_20261020/AFFINE_CIRCULANT_IDEAL_OBSTRUCTION.md)
proves that for every prime \(p\equiv1\pmod4\), \(p>5\), no Boolean-affine multiplier-equivariant map to sign-even circulant graphs is surjective on isomorphism types. Every \(D_c\), when \(p\nmid c\), is such a Boolean-affine equivariant rule. Thus \(D_c\) cannot be a type bijection at any such admissible prime, whether its failure is non-even-valuedness or insufficient type coverage. This consequence needs no bad-prime infinitude theorem.

For completeness, a suitable admissible prime \(p\equiv1\pmod4\), \(p>5\), exists for every fixed nonzero \(c\) by an elementary argument: any prime divisor of \((10c)^2+1\) is different from 2 and 5, does not divide \(c\), and has \(-1\) as a square, hence is 1 modulo 4.

The exception \(c=1\) has \(D_1(S)=\varnothing\) for every input. It never has the complete-graph failure, but it is not a type bijection uniformly in \(p\); at \(p=7\) there are two target types. The Boolean-affine obstruction also excludes it at every \(p\equiv1\pmod4\), \(p>5\).

Neither argument proves that **all** fixed nonlinear local dilation rules are impossible. The incoming sentence “Theorem W shows that no local rule does this” overstates its derivative-family scope. Nonlinear AGL-equivariant algorithms are not excluded, and an [all-prime derivative-path construction](../circulant_necklace_bijection_20261020/CIRCULANT_NECKLACE_BIJECTION.md) is already available in the repository.

## Exact finite controls

The positive examples include the factor-16 residual family, an eighth power, and a sixteenth power:

| \(c\) | \(\ell\) | \(N_\ell(c)\) | chosen prime \(p\) | \(\operatorname{ord}_p(c)\) |
|---:|---:|---:|---:|---:|
| 2 | 5 | 11 | 11 | 10 |
| 3 | 3 | 7 | 7 | 6 |
| 4 | 3 | 13 | 13 | 6 |
| 16 | 3 | 241 | 241 | 6 |
| 256 | 3 | 65281 | 97 | 6 |
| 65536 | 3 | 4294901761 | 193 | 6 |
| -2 | 3 | 7 | 7 | 6 |
| -3 | 3 | 13 | 13 | 6 |
| -4 | 5 | 341 | 11 | 10 |
| -16 | 7 | 17895697 | 29 | 14 |
| -256 | 7 | 282578800148737 | 29 | 14 |
| -65536 | 7 | 79229371458530977775699951617 | 29 | 14 |

Additional complete factorizations are \(65281=97\cdot673\) and \(4294901761=193\cdot22253377\). Every displayed factor is prime, and every factor in these two factorizations gives order 6.

The exact controls are not the infinitude proof. The main standard-library checker validates the table, all factors just listed, all 77 prime divisors in its bounded quotient universe, the distinct-order disjointness mechanism, explicit subgroup coset membership, and prime-power lifts through exponent 4. It also constructs invariant tournament sets and checks the complete derivative on every nonzero residue, including nonunits, for seven ring-circulant controls at moduli \(49,121,169,343,841,1331,2197\). The boundary checks cover all 94 odd primes below 500. Its negative controls show why the filter and prime exponent matter:

- Without \(\ell\nmid c+1\): \(c=2,\ell=3,N=3,p=3\) has order 2, not 6.
- Negative-sign counterpart: \(c=-4,\ell=3,N=21,p=3\) again has order 2.
- With composite odd exponent: \(c=3,\ell=9,N=4921\) has prime divisor 7 and order 6, not 18, despite \(\gcd(N,c+1)=1\).
- \(c=1\) never lies in the negative odd subgroup, while \(c=-1\) always does, throughout the tested odd-prime universe.

Reproduction from this packet, comparing rather than overwriting the pinned receipts:

```sh
tmp=$(mktemp -d)
python3 check_bad_primes.py > "$tmp/controls.json"
python3 -O check_bad_primes.py > "$tmp/controls.optimized.json"
cmp controls.json "$tmp/controls.json"
cmp controls.json "$tmp/controls.optimized.json"
python3 audit/check_cyclotomic_primes.py > "$tmp/independent.json"
cmp audit/receipt.json "$tmp/independent.json"
```

The [producer verifier](check_bad_primes.py) has its normal [exact receipt](controls.json). The [independent checker](audit/check_cyclotomic_primes.py), its original [receipt](audit/receipt.json), and its [written audit](audit/REVIEW.md) were supplied by an independent reviewer and copied unchanged. The independent finite checks cover 89 admissible pairs \(c\in[-20,20]\setminus\{-1,0,1\}\), \(\ell\in\{3,5,7\}\), using trial division for the least prime factor and exact tests for the proper order divisors.

## Integration notes and scope

Incoming source: [THM-4724, §3, Theorem W](../../../01-canon/theorems/THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md).

Recommended local repair: replace the universal infinitude argument of W(ii) with the alternating-quotient proof, change “every integer \(c\ne1\)” to “every nonzero integer \(c\ne1\),” and remove the Chebotarev dependency for this assertion. Any retained finer arithmetic-progression statements need their own separately scoped proofs; this packet neither imports nor retracts them. Correct the §4 all-local-rules sentence to the derivative-family conclusion, and distinguish the already published Boolean-affine theorem.

Inheritance / concept board: maximal odd subgroups and order valuation; alternating quotients and denominator exclusion; affine-rule obstruction versus odd-output witnesses; reduction modulo prime powers. The nearest proved obstruction is the Boolean-affine theorem. The hostile boundary is \(c=1\); the corrected near miss is \(c=0\) in a nonzero-dilation theorem; the least-used sidecar is the odd kernel of unit reduction, which supplies the full prime-power tower.
