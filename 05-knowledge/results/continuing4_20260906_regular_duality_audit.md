# Independent audit of regular reflection and paired characteristic divisors

**Status: PASS -- independent analytic proof audit and complete FINITE-EXACT identity reconstruction.** The all-height paired-divisor and boundary-reflection laws are accepted. The six specified regular families have first nonzero constant-term moment exactly g or2g, relative to the inherited proved/CITED real-root supplier. All-height residual positivity and general Laurent two-rung noncancellation remain OPEN. No mathematical repair was required; the producer's stdout and displayed reproduction-path nits were corrected before final filing.

## 1. Scope and inherited boundaries

Read the complete [producer proof](continuing4_20260906_regular_duality.md), source and full rational JSON, the incoming [carried singular-block proof](third_20260906_laurent.md), and [THM-4436 / complete factorial-row simple negative roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md). The latter is analytically proved relative to its explicitly CITED finite-preserver theorem. This audit introduces no new literature dependency or external-priority claim.

The new incoming [all-height beta-step note](long_frontier_sep06_allh.md) and its audit were also read. Its two sign-propagation hostiles have incorrect interior factorial amplitudes. They defeat the stated qualitative propagation rule, while retaining the genuine first rows and carry/top anchors, but they do not refute the present exact factorial family or subsume its paired divisors and reflection. No inference between different rows is made from separate real-rootedness.

The three scopes must remain separate: the divisor law is algebraic for every integer H>=1 and formal parameter z; the sign certificate has H=1,...,6 only; the actual first-return interpretation further requires integer z>=1 and gcd(z,3H)=1.

## 2. Actual fibres and normalization

For support(-3H,z,6H+3z), let g=3H+z. At mass m the literal equation is

```text
g*n_beta+3g*n_gamma=3H*m.
```

At m=g or2g, solve this equation for n_beta while n_gamma ranges over every possible nonnegative count, and retain n_alpha=m-n_beta-n_gamma>=0. This gives exactly the producer's progressions

```text
(z+2j,3H-3j,j),          0<=j<=H,
(2z+2j,6H-3j,j),        0<=j<=2H.
```

The independent source performs this literal enumeration before constructing the rows. Dividing each multinomial by binom(g,H), or binom(2g,2H), gives precisely the stated monic rising-product Phi and Psi. The monomial factors are X=alpha^z beta^(3H) and X^2, with tau=alpha^2 gamma/beta^3. There is no negative-index channel in this new regular fibre; this follows from complete nonnegativity enumeration and does not delete a term from the old carried row.

The real-root supplier is typed as A=1,B=3, height H, x=z, r=0 and base remainder0. It gives H distinct strictly negative roots of Phi for every allowed actual parameter. A=2 with unbounded base remainder would be inapplicable and is not used.

If gcd(z,3H)=1, then gcd(g,3H)=1 and the charge equation forces g|m. The j=0 channel makes g admissible. Without the gcd hypothesis, z=3H gives support(-3H,3H,15H), first support return2, and all-unit coefficients give CT(f^2)=2; g=6H is not the first return. The audit retains this hostile at every certified H.

## 3. Paired divisibility: formal lifting is sufficient

Fix m in1,...,2H, set ell=ceil(m/2), and delta=z+m. In Phi's coefficient of t^j the consecutive factors are z+2j+1 through z+2H. They contain exactly one vanishing factor if j<ell, and none if j>=ell. In particular the constant coefficient has a simple delta zero, while the coefficient of t^ell at delta0 is nonzero. Thus

```text
Phi_(z=-m)=t^ell*a(t),             a(0)!=0.
```

The corresponding factors for Psi are2z+2j+1 through2z+4H. They contain exactly one zero when j<m and none when j>=m. The factor's derivative is2, so the zero is still simple.

Hensel lifting over C[[delta]] splits off the degree-ell factor because t^ell and a are coprime. With delta=epsilon^ell and t=epsilon*v, the leading small-root equation is a(0)*v^ell+A=0, where A is the nonzero derivative of Phi's constant coefficient. Its ell roots are nonzero and distinct, so each lifts by the formal implicit recursion. At every small root, Psi's terms of index j<m have epsilon order at least ell+j, and its remaining terms have order at least m>=ell. Every small response eigenvalue therefore has epsilon order at least ell.

Any j-fold product of small eigenvalues has order at least ell*j. Their symmetric functions are coefficients of the multiplication operator on the lifted small factor, hence belong to C[[delta]], and have delta order at least j. Chinese remainders supply a complementary operator of size H-ell with regular formal coefficients; when ell=H this complementary part is empty. Any contribution to the kth full characteristic coefficient uses at least max(0,k-H+ell) small eigenvalues. This proves the required valuation at z=-m.

The distinct rational factors z+m are coprime. Grouping the two values m=2ell-1,2ell proves the entire displayed paired divisor in Q[z]. Cancellation can only increase these lower orders; exact all-height multiplicities are not asserted. This argument does not require real roots at the negative formal parameter or a physical interpretation there.

Giving t weight2 and z weight1, the first row's coefficient degree is2(H-j), and the response coefficient degree is4H-2j. Monic reduction preserves these weighted bounds. The multiplication entry in row i,column j has degree at most4H+2j-2i, so index sums cancel in each principal minor and deg c_(H,k)<=4Hk. The divisor degree is k(k+1). Thus the residual degree is at most4Hk-k(k+1), whose maximum over k is3H^2-H. No unproved residual positivity enters either argument.

## 4. Reflection retains the actual scalar normalizations

At old height h and1<=r<=h, put H=h-r,z=2r+1. After setting j=r+k, the old first coefficient is

```text
(2h+1)! H! /[(3H-3k)!(2r+2k+1)!k!],
```

which is exactly Phi_(H,z)'s coefficient. Below r it vanishes. The same factorial substitution at e=2r+j gives Psi_(H,z)/(4h+2)! and vanishes below2r, including the old inverse carry. This proves both coefficient identities for every h,r, including H=0 where the complementary polynomials are1.

Reflecting u to1/u and swapping the endpoint coefficient names gives

```text
(alpha_new,beta_new,gamma_new)=(gamma_old,beta_old,alpha_old),
tau_new=tau_old,
X_old*tau^r=X_new,
g=3h-r+1=3H+z,
binom(g,2h+1)=binom(g,H),
(2g)_(4h+2)/(4h+2)!=binom(2g,2H),
gcd(g,6h+3)=gcd(z,3H).
```

Thus the actual normalization reflects, not only the complementary coefficient list. The boundary construction directly supplies odd z>=3. Other positive z in the regular theorem are justified by Section2's direct complete-fibre construction, not by pretending every such z is an old negative-integer boundary.

The zero-root block must still be excluded from the specialized-response comparison. Independently, at old h=1,x=-1 the generic polynomial quotient response specializes to **-1/90720**, while the raw specialized Laurent row is t^2/720 and vanishes at t=0. The difference is the old inverse-carry ratio after cancellation. At a nonzero complementary root there is no ambiguity: the generic remainder agrees with the specialized Laurent response, and the factor t^(2r)/(4h+2)! is positive on a nonzero real phase. No singular-block equivalence is claimed.

## 5. Complete independent certificate reconstruction

The independent engine uses literal integer multinomial fibres, rational Horner evaluation in the companion operator, clears a common denominator, and computes the **integer Berkowitz characteristic polynomial**. The producer instead uses symbolic monic polynomial division and Faddeev--LeVerrier over Q[z]. No producer module is imported.

For each H the audit evaluates at every integer z from1 through3H^2-H+1. There are respectively3,11,25,45,71,103 points, totaling258 matrices. It divides each literal characteristic coefficient by the analytically proved paired divisor and compares with the entire supplied rational residual polynomial. Both sides have the established residual degree bound. Agreement at strictly more points than that bound proves each polynomial identity, rather than serving merely as a positive sample bank.

All21 polynomials and all833 strictly positive rational coefficients match. There are1253 characteristic matches,4275 complete fibre terms and111 primitive parameter controls; the latter also independently exclude every smaller positive return mass. The complete coefficient law, local zero prerequisites, reflection identities, both scalar normalizations, literal negative-root controls and the singular-carry hostile are retained.

For H<=6 and actual z>=1 every divisor factor is positive, so every nonleading characteristic coefficient is positive. Therefore det(wI-M)>0 for real w>=0. At a root rho of Phi, Psi(rho) is a real eigenvalue of M by the typed real-root supplier. It must be strictly negative. This is a normalized real response statement, not a sign assigned to a raw complex constant term.

Consequently the first nonzero constant term is g unless tau is a first root, in which case it is2g. Both alternatives occur for every admissible parameter by choosing nonzero coefficients that realize a generic tau or a negative first root. As a separate literal control at H=z=1, all coefficients1 give CT(f^4)=8; choosing alpha=beta=1,gamma=-1 cancels every moment below8 and gives CT(f^8)=-224. Both g and2g are below the endpoint width3g. No proof for H>=7 is inferred.

## 6. Freeze and promotion basis

[Independent source](../../04-computation/continuing4_20260906_regular_duality_audit.py) and [output](continuing4_20260906_regular_duality_audit.out) pass **14,026 always-active exact gates** normally and with optimization. Captured output bytes are identical and contain LF newlines; no post-normalization was applied. The certificate is read only as pinned data, and an explicit path supports filing it separately from the computation source.

```text
python -B 04-computation/continuing4_20260906_regular_duality_audit.py --certificate 04-computation/continuing4_20260906_regular_duality_certificate.json
python -B -O 04-computation/continuing4_20260906_regular_duality_audit.py --certificate 04-computation/continuing4_20260906_regular_duality_certificate.json
```

```text
independent source 12635f429915f42a9ad839a570d4a092c33140b1f3b8db829d9c743c9cf4ea83
independent output 77b774e85301a0012c24bbd45623abe9a9515993ddc37934a71c62b4d7a3cfb3
producer certificate 0d5a65f03fc4f4295f3db38bca8609375cfea8805f21499a28e9a5e0d9a1ccd4
producer source as read d649033b472372558d671a3f3d1b9a2b6c80ca4edcde28fbc9384a02cf87ccce
candidate report as read 3d6b264e11043e9c8766a3d85b0b2e2457fb6c0ec52d4301db13717fb4adb3a9
```

The complete analytic audit plus distinct bounded-degree identity reconstruction supports promotion to **PROVED analytically / FINITE-EXACT / INDEPENDENTLY AUDITED**, with the cited real-root dependency and all stated open boundaries retained. The root integrator owns the filed status change and candidate-hash lineage. This audit made no repository or Git changes.
