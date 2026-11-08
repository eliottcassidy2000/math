# Exact two-J fixed-length law for every odd modulus

Status: PROVED, with fresh independent review recorded in AUDIT.md when published.
This is a reconstruction of Session225 reasoning from retained conversation, not restoration of unavailable original files. No priority claim is made.

## Model and scope

Work modulo a fixed odd positive integer q. Start at t=(D,E,U,d)=(1,0,2,1), with
O(D,E,U,d)=(D,0,2U,d),
Z(D,E,U,d)=(D-E+U,d,2U-d,d),
J(D,E,U,d)=(D-E+U,D-E,D+U,D),
and F=D-E+U.

For L>=2 choose uniformly one of the binomial(L,2)*2^(L-2) words of length L having exactly two Js. For every residue r,
Pr(F=r)=L_q(r)/q+O_q(1/L).
For fixed q this also holds in total variation. The q=1 case is trivial.

Put g=gcd(d,U,q), h=gcd(D,U,q). The exact factor is
L_q(r)=q^(-3) sum_(d,U,D mod q) g*h
  * (3/4*[g|D]+1/4*[g|D+1])
  * (3/4*[h|r]+1/4*[h|r+d]).

This is a fixed-length theorem, NOT a height-truncated multiplicity theorem. It gives no scalar coverage, cutoff, finite-bank completeness, or proof of the Chan-Pak conjecture.

## 1. Normalized graph and finite permutation group

After one O/Z letter, E belongs to {0,d}; coincident values are counted once. Fix d, g=gcd(d,U,q), and D=c mod g. These define a broad invariant class. Its directed O/Z graph is two-in/two-out regular, counting edges with multiplicity. If d!=0, each E=0 target has the two O predecessors (one from each E slice), and each E=d target has the two Z predecessors; solving for U uses invertibility of 2. If d=0, O and Z are permutations of E=0. Therefore all strongly connected components (SCCs) are closed, and uniform measure on each SCC is stationary.

Do not identify a broad gcd class with an SCC. It may split. In each SCC with d!=0 the two E slices have equal size, and O maps either slice bijectively onto its E=0 slice. To see this, O is a permutation on the whole E=0 slice, so its restriction to a finite closed SCC is bijective there. Its other predecessor at each target lies in the same closed SCC. The d=0 case has just one slice.

On E=0 put ell=ord_q(2). The realizable finite permutation group contains
S(D,U)=(D,2U),
A(D,U)=(D+U,U+eta), eta=-d/2.
Here A is Z followed by O^(ell-1). Inverses are realizable by positive powers of these finite permutations. For integers m,
A^m(D,U)=(D+mU+eta*m*(m-1)/2,U+m*eta).

Reduce D-c,U,d by g, with reduced modulus Q=q/g. Then gcd(d,U,Q)=1 in the reduced coordinates. For a selected prime power p^k||Q choose a CRT idempotent n, equal to 1 on p^k and 0 at every other prime-power component. The power A^n is A on the selected component and identity elsewhere. This follows from the power formula, since 2 is invertible. Conjugating this localized A by any power of S remains localized. No isolation of S itself is asserted.

For p!=3 put B=SAS^(-1)=(D+U/2,U+2*eta). Comparing A B with B A gives a pure D translation by +/-3*eta/2. Thus if d is a unit on this component, arbitrary localized D translations are available.

For p=3, there exists m with 2^m=-1 mod 3^k: ord_(3^k)(2)=2*3^(k-1), by the elementary valuation identity v_3(4^(3^j)-1)=j+1. Conjugating localized A by S^m gives (D-U,U-eta) on the selected component. Composing it with A gives a pure D translation by +/-eta. If d is a unit these again generate arbitrary local D translations. Outside the selected component the conjugation word is identity; it does not require S^m to be a global reflection.

## 2. Line-transversal lemma

Fix a,b modulo Q and V=D-aU-b. On the E=0 slice of every SCC, V is uniformly distributed modulo Q.

At a prime where d is a unit, localized D translations prove the statement conditional on all other prime coordinates. At a prime where d is nonunit, U is a unit. Along localized A powers,
V(A^m)=V+mU+eta*m*(m-1-2a)/2.
The difference at m,n equals
(m-n)*(U+eta*(m+n-1-2a)/2).
The second factor is a unit, so this polynomial permutes the selected prime-power ring. The local action has period dividing p^k. Averaging over it makes V uniform on that component, conditional on the other components. Successively averaging these localized operations proves the joint assertion modulo Q.

Consequently each SCC E=0 slice meets D=aU+b in exactly fraction 1/Q of its vertices. Uniform U on the gcd stratum, placed on this line, therefore gives SCC masses proportional to SCC sizes. Its Cesaro projection is uniform over the entire broad class even when that class is disconnected. The same holds for a line on E=d, since O carries it bijectively to a line on E=0 and preserves SCC membership.

The lemma also shows that from every normalized starting state, the projected output F is uniform on its invariant D-coset modulo g. On E=0 use a=-1. For E=d, its bijection under O to E=0 writes the OLD readout as D+U_new/2-d, another allowed linear form. Both slices give the same coset since g divides U,d.

This is output uniformity. It does not assert independent uniform D,U within every SCC.

## 3. Initial and first-J laws

For d=1, arbitrary localized D translations and A moves give transitivity: adjust U componentwise using A, then correct D. Thus the whole d=1 class is one SCC. It is aperiodic because O fixes states with E=U=0. Its stationary law has independent uniform D,U and independent fair E in {0,1}.

Write old D=d, old U=u, old E=e. Immediately after J, put V0=d+u, uniform independently of d,e. After the first following O:
E=0, U_new=2V0, D=U_new/2-e.
After the first following Z:
E=d, U_new=2V0-d, D=U_new.
Both preserve gcd(d,U,q) and leave U_new uniform at fixed d. Apply the line-transversal lemma: the respective D cosets are -e and 0 modulo g. Averaging the two fair choices (e and the first O/Z) gives:
- d,U independently uniform modulo q;
- g=gcd(d,U,q);
- D uniform on coset 0 mod g with weight 3/4 or -1 mod g with weight 1/4;
- E independently fair in {0,d}, with coincident values counted once.

Thus relative to independent uniform (d,U,D), the density is
w=g*(3/4*[g|D]+1/4*[g|D+1]).

## 4. Second J and the finite sum

Put h=gcd(D,U,q) before the second J. After J, the first O gives D coset -E modulo h; the first Z gives coset 0. Arbitrary-state projected output uniformity from section 2 applies regardless of the finer SCC. Average the preceding fair E choice in {0,d}: the output is uniform on coset 0 mod h with weight 3/4, or -d mod h with weight 1/4. Integrating against w proves the displayed finite-sum formula.

## 5. Fixed-length convergence

Let P be the fair O/Z transition matrix on the finite state space and Pi its Cesaro projector. Exactly-two-J words are uniform over three O/Z block lengths a,b,c>=0 with a+b+c=L-2. Their distribution is the coefficient average of the product of three resolvents (I-zP)^(-1), separated by deterministic J maps.

Finite stochastic matrices have semisimple peripheral eigenvalues. Each peripheral resolvent pole is simple. The first block is restricted to the initial aperiodic d=1 SCC, so its only peripheral pole is at 1. Hence the triple product has at most an order-three pole at 1 and at most order two at other peripheral points. At 1 the order-three coefficient is exactly the three-projector cascade. Coefficients of order-two terms are O_q(L), while binomial(L,2) is quadratic; eigenvalues strictly inside the unit circle give exponentially smaller terms with fixed polynomial factors. This proves the asserted O_q(1/L) limit.

## 6. Branch factorization and correlation

Define
K_(ab),q(r)=q^(-3) sum_(d,U,D mod q) g*h*[g|D+a]*[h|r+b*d],
for a,b in {0,1}. CRT factors EACH branch:
K_(ab),q(r)=product_(p^k||q) K_(ab),p^k(r mod p^k).
Therefore
L_q=9/16*K_00,q+3/16*K_01,q+3/16*K_10,q+1/16*K_11,q.

The common branch choices survive CRT. One must not replace this mixture of products by a product of averaged local factors. For example L_15(0)=5179/3375, whereas L_3(0)L_5(0)=5134/3375, a difference of 1/75. Explicit prime-power formulas, if included, are supplied and independently audited in KERNELS.md.

## Inheritance and loss ledger

The retained Session223 all-3-adic argument and THM-4916 odd-prime argument motivate the cascade. This reconstruction supplies its all-odd-modulus proof directly, without depending on unavailable temporary files. The hostile is a split gcd class, rather than a primitive-root hypothesis. Source: full normalized O/Z/J states. Quotient: gcd/coset data and scalar residue. Preserved: this fixed-length projector law. Lost: SCC labels, joint D/U dependence, real height and history. Restoration sidecar: line transversals and the shared branch label. Height-conditioned independence and scalar coverage remain OPEN.
