# Independent audit of the exact one-digit Deuring precision law

**PASS: full analytic proof audit and separate exact controls. No repair
required.** The theorem and adjacent-jet refinement in
[the producer report](overnight9_20260906_jets_deuring.md)
are accepted on the stated equilateral, three-node, prime-order domain.

The audited producer source SHA256 is
`e7e4846f4bdd9cb5b2174571365192436c1578bd110dd71f6dde41a53faa3ab6`;
its frozen LF output SHA256 is
`13b6475b7110129aee0e7c92bafb6aa995465fc3e153e965e1cfe176eb576f08`.
The referee changed no producer files.

## Statement and scope accepted

Let p be any odd prime, k=(p-1)/2, m=k+1, and let three integer nodes
have all three pairwise p-adic depths equal to e. Normalize them by
translation and a p-adic unit change of variable to `0,p^e,p^e a`,
where a,a-1 are units. Put `sigma=[H_p(a)=0 mod p]`.

For e>=1 the sharp largest Smith exponent is

```text
L=(3m-1)e-sigma.
```

At e=0 all exponents are zero. The determinant and penultimate
determinantal valuation are correctly
`D_(3m)=3m^2 e` and `D_(3m-1)=(3m^2-3m+1)e+sigma`.
At a Deuring-zero residue all three normalized jets of order k-1 are
units. Their cost is exactly `(3m-2)e`, tying the top cost precisely
at e=1 on that branch.

The multiplicity m=(p+1)/2, three-node condition, and equilateral depth
are essential scope restrictions. This is not an arbitrary-node packet
formula, a close-pair theorem, or a full intermediate Smith partition.
The eighth-round complete p7 four-jet partition remains a separate
result. The proof handles every higher unit lift, not merely integer
representatives of the residue classes.

## Exact reciprocal normalization

The inherited actual inverse-denominator equality is
[THM-4443](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md),
already independently audited in the preceding rounds. Complete Hasse
jets make translation and unit coordinate changes unimodular. For
normalized differences u,v, expansion of `((u+T)(v+T))^(-m)` at order k
gives exactly the numerator and unit denominator in producer equation(4).

With the producer's P,R,T polynomials, the three top coefficients are
literally

```text
at0: P(a)/a^p,
at1: R(a)/(1-a)^p,
ata: (-1)^k T(a)/[a(a-1)]^p.
```

All denominators are units. Their common scale is `(2m+k)e=(3m-1)e`.
The binomial reductions to the Deuring polynomial and its reflected
and reciprocal transforms are correct. Thus ordinary residues make
all three numerators units; at a Deuring-zero residue all are divisible
by p. No extra factorial factor, orientation sign, or power of p has
been omitted.

## Independent all-prime ODE and Wronskian derivation

The coefficient ratio is

```text
c_(j+1)/c_j = (k+j+1)(k-j)/[(j+1)(2k-j)].
```

Its denominators are nonzero for 0<=j<k. Multiplying out gives exactly
the coefficient equation for
`a(1-a)P''+(-2k-2a)P'+k(k+1)P=0`, including the constant and leading
coefficient boundaries. With p=2k+1, the operator
`L0=a(1-a)D^2+(1-2a)D+k(k+1)` therefore satisfies
`L0P=pP'`. Reflection gives `L0R=-pR'` with the indicated minus sign.

The mod-p reflection relation makes `(P-R)/p` integral. P and R have
the same exact leading coefficient, so the divided companion G has
degree at most k-1. Subtracting the two characteristic-zero ODEs gives
`L0G=P'+R'`; reducing modulo p gives `L0G=2P'` and `L0P=0`.

For `W=PG'-P'G`, direct differentiation yields

```text
[a(1-a)W]' = P L0G-G L0P =2PP'=(P^2)' mod p.
```

The degree bound is decisive: `deg(a(1-a)W-P^2)<=2k=p-1`.
Every nonconstant monomial of this degree range has nonzero derivative
in characteristic p. Hence this polynomial is constant. Its value at
zero is `-P(0)^2=-1 mod p`, because
`P(0)=binom(p-1,k)=(-1)^k mod p`. The exact congruence is therefore

```text
a(1-a)(PG'-P'G)=P^2-1 mod p.
```

At an admissible Deuring zero this becomes
`a(1-a)P'G=1 mod p`. Both P' and G are units. For every lift of that
residue, `P-R=pG` has valuation exactly one. Since P and R are already
divisible by p, their simultaneous minimum valuation is exactly one;
adding the third numerator cannot decrease it below one or raise it.
This proves the joint cap on all higher lifts. It simultaneously proves
the simplicity needed for the individual-numerator Hensel hostile.

The low-degree condition is used correctly. A zero derivative could
otherwise conceal pth powers, so neither a higher-genus nor an
arbitrary-node extension follows by formal analogy.

## Actual precision, adjacent jets, and the p3 boundary

Every lower-order normalized reciprocal coefficient is integral.
For l<k its cost is at most `(2m+k-1)e`. The top cost is
`(2m+k)e-sigma` with sigma in {0,1}; it dominates or ties all those
costs for e>=1. The attained inverse law therefore gives the exact
largest exponent. This proves the main theorem without needing the
stronger adjacent-jet assertion. The unit determinant separately
discharges e=0, where inserting the positive-depth formula would be
wrong on the zero branch.

The adjacent assertion also passes. With

```text
A(a)=[X^k]((X-1)(X-a))^k,
B(a)=[X^(k-1)]((X-1)(X-a))^k,
```

the claimed equation `(k+1)B=-a[kA+(1-a)A']` is exact over the
integers. Independently, its coefficient of a^j follows from

```text
(k-j+1)binom(k,j-1)^2+j binom(k,j)^2
   =(k+1)binom(k,j-1)binom(k,j).
```

The producer's homogeneity and root-translation derivation is valid
as well: reciprocal reversal sends the coefficient of X^(k+1) to B/a.
At A=0 modulo p, a,1-a,k+1 and A'=P' are all units, so B is a unit.

For the cubic `F_a=X(X-1)(X-a)`, its kth power has degree
`3k<2p-2`. Modulo p, translation of the coefficient of degree p-2
therefore receives contributions only from degrees p-2 and p-1.
The result is exactly `B-A x_i`. At a zero of A it equals the same
nonzero B at every endpoint. The Frobenius reciprocal identity
multiplies this by a unit, proving that all three actual adjacent
coefficients are units. This does not require a new root condition.

At p3,k1 the polynomials are
`P=2+2a`, `R=2a-4`, `G=2`. The Wronskian constant remains correct,
and the adjacent order is zero. There is no exceptional division by
three in the proof. As a sharp individual/global hostile, a=-1 makes
P vanish identically while R=-6 has valuation one; the actual largest
loss is still `5e-1` at e>=1. More generally the unit derivative permits
successive choices of Hensel digits killing P to arbitrary order, while
its reflected companion remains of valuation one.

## Independent exact controls

The [referee source](../../04-computation/overnight9_20260906_jets_independent_audit.py)
uses only the standard library, with no producer import, SymPy, or
integer Smith routine. It reconstructs P from its adjacent coefficient
ratio, reflects by Horner composition, and checks the entire integer ODE
and mod-p Wronskian identities. Direct powers of the cubic separately
recover both coefficients A and B. Actual reciprocal jets are computed
by recursively solving their defining polynomial inverse, rather than
using the producer's negative-binomial formula.

The bank contains thirteen odd primes from3 through103, including
ordinary-only and mixed residue behavior. At primes through31 it checks
all2,255 admissible lifts modulo p^2, including ordinary classes as
well as every Deuring-zero class. All joint caps agree with the theorem.
One hundred four literal reciprocal observers, spanning four depths,
include every observed order, negative a, and higher unit lifts. Twenty
full Hasse matrices are reduced with independently checked p-integral
Gaussian pivots; their largest exponents and full determinant valuations
match. These are controls for the all-prime proof, not its justification.

Both modes pass **7,269 explicit gates**, with byte-identical LF output:

```text
python -B 04-computation/overnight9_20260906_jets_independent_audit.py
python -B -O 04-computation/overnight9_20260906_jets_independent_audit.py

source ce5b80b93a3fa39650e4cb01b6eaf92223476ea4eff012b147df501ef774c29a
output 0c8b33e3a0e421fee84c8c809fa72f594b2654a317015f82cc8cdfbc8322a1c3
semantic bank a7cbcf5152cb41411251b9d85850d90b84da29f1786fbdfd02cc7fbc3e07b169
```

The eighth report's primary-source supersingular interpretation was
independently checked in the preceding audit. This extension uses no
additional literature theorem or external priority claim. Its novel
mathematical content is the divided reflected companion and its exact
joint cancellation cap, not the classical Deuring polynomial itself.
All files remain outside the repository for parent-managed integration;
no Git or shared proof surface was edited by this audit.

**Filing:** root integrated these audited artifacts in the ninth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
