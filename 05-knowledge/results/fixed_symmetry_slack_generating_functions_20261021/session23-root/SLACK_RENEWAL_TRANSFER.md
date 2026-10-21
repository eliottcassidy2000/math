# A shared renewal law for two different tournament constraints

Status: PROVED from the independently audited fixed-symmetry-slack compiler and the stated renewal proof. The compiler proof is pinned at SHA256 22cb6aca899d9a485b2432d37a90d13040dfc6353ee41529a243d8489aa974a1; independent transfer review passed.

## Exact source and distinction
The fixed-h limit laws in RENEWAL_LIMIT_LAWS.md use a finite alphabet of nonrigid strong SCCs and unbounded neutral SCC gaps. The symmetry-slack compiler in session23-connections/SYMMETRY_SLACK_RATIONALITY.md instead uses canonical maximal W-modules, W=SEQ(K1,C3), and bounded short decorations. These carriers are different: C3[K1,W,K1] is strongly connected of symmetry slack3 for every |W|>=2, so bounded slack permits unbounded strong tournaments. No fixed-h finite-atom assertion is being reused outside its scope.

The audited compiler theorem gives, for fixed d>=1,

F_d(z,y)=sum_c z^(a_c)y^(b_c)(A(z,y)-D_d(z,y))^(ell_c),
A(z,y)=1/(1-z-yz³),

with finitely many actual configurations, nonnegative integer weights and a finite polynomial D_d. Here y marks the actual minimum feedback f, not a display-order feedback proxy. For d0, F0=A directly.

## Ordinary and merged-node theorem
For every nonempty fixed slack d, choose a tournament isomorphism class uniformly at order n. Such classes exist for all sufficiently large n: adding an ordered K1 prefix preserves f and moved support s. Put rho+rho³=1, mu=rho²/(1+3rho²), sigma²=rho²(1+rho²)/(1+3rho²)³. Then

E f=mu*n+O_d(1), Var f=sigma²*n+O_d(1),
(f-mu*n)/sqrt(sigma²*n) -> Normal(0,1).

In particular s=3f-d implies s/n ->3mu in probability, and (s-3mu*n)/sqrt(9sigma²*n) ->Normal(0,1). The limiting moved-vertex fraction is approximately0.5827620121.

Proof: let ell be the largest number of long slots in the finite compiler. For a nonempty layer ell>=1 since the layer is infinite; if every term had ell0 it would be a polynomial. Subtracting the finite head D_d changes neither the pole location nor the principal part of A. Finite nonnegative top-weight summation prevents cancellation of the leading moving pole. Hence exactly the moving-root/characteristic-function proof of Theorem C in RENEWAL_LIMIT_LAWS.md applies, with n^(ell-1) in the coefficient amplitude. The same implicit root derivative constants arise, even though the preserved tournament predicate differs.

The actual automorphism order obeys log|Aut(T)|=f(T)log3+O_d(1). In each canonical decorated quotient the group order is3^(sum t_i)|K_0|: the internal rotations form the kernel and every word-preserving quotient automorphism lifts. The exact compiler gives f=sum t_i+C_0 with 0<=C_0<=d and |K_0|<= (4d+1)!. Therefore the logarithmic automorphism order has limit mu log3 per vertex and Gaussian variance constant sigma²(log3)². This asserts only convergence in probability/distribution for logAut; no O(1) variance correction is inferred merely from adding an arbitrary bounded random variable.

The corresponding converse-merged fixed-slack nodes have the same limits. The compiler's self-converse series has all poles at modulus at least sqrt(rho), while the ordinary series has a positive rho pole. Thus the self-converse contribution is exponentially smaller up to a polynomial factor. All statistics above are converse invariant, so the two-to-one quotient on the non-self-converse part preserves their law asymptotically.

## Self-converse variance doubling
On each parity where self-converse fixed-slack classes occur at arbitrarily large orders, the feedback law instead satisfies

E f=mu*n+O_d(1), Var f=2sigma²*n+O_d(1),
(f-mu*n)/sqrt(2sigma²*n) ->Normal(0,1).

The finite compiler pairs every long word with its reversal, with at most one palindromic long word. A palindrome itself is a half-word followed by an optional center letter and the reversed half. For each fixed choice of finite decoration and optional center, y marks twice the half triangle count plus a bounded integer; z marks twice the half size plus a bounded integer. For a middle palindrome center of size h in {0,1,3} with triangle count e in {0,0,1}, respectively, its long half contributes z^h y^e [A(x,y²)-D_floor((d-h)/2)(x,y²)], x=z²; define D_q=0 for q<0. Removing bounded short words only changes finite heads. Equivalently, after fixing vertex parity and putting x=z², the moving denominator is1-x-y²x³. The characteristic exponent per n is−(1/2)log R(exp(2t)), giving mean coefficient mu and variance coefficient2sigma². Finite supported mixtures have bounded differences of means, as in Theorem D of RENEWAL_LIMIT_LAWS.md. Unsupported parities make no assertion.

## What the connection preserves and loses
Map: each of two distinct canonical tournament decompositions produces finitely many decorated renewal word slots over lengths1 and3. Preserved by the moving pole: triangle/feedback fugacity and the linear-size count. Lost: short exceptional decorations, actual edge directions, path positions, group action beyond its bounded factor. The common constants do not identify the two tournament families or their metagraph wiring. The dependence of mirrored halves, not an independence approximation, explains variance doubling.

This is a transfer of a proved mechanism, not evidence that every bounded tournament statistic has this limit. A different unbounded alphabet changes the root, a periodic alphabet requires residue-class treatment, and negative top-weight cancellations or growing defect invalidate the stated finite positive argument.
