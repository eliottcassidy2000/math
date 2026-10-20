# An infinite 2,3-free near-critical two-edit norm barrier

20 October 2026. Status: **PROVED for the explicit near-critical family, independently audited.** This addendum does not modify the frozen two-edit checkpoint. It strengthens a norm-method boundary, not an integer-cycle classification.

## 1. Exact theorem and scope

Write alpha=log_2 3. For every integer a>=2, put

    p=2*3^a,
    L=ceil(p alpha),
    q=min{m>=L : m=1 mod 6},
    r=q-p,    n=p-r=2p-q.

Start with the zero-intercept lower mechanical word of shape (p,q). Perform the overlapping two-plus chain whose follower has residue p-1 and whose leader has residue p-1-r, moving the leader first. The residue coordinates and norms are those of the [two-edit sparse-form theorem](../collatz_two_swap_norm_barrier_20261020/symbolic/TWO_SWAP_SPARSE_FORMS_AND_BARRIER.md).

Its canonical polynomial and positive norm modulus are

    Q_p(Z)=3-Z^(p-1)+Z^n-Z^(n-1),
    R_p=Res(Z^p-2,Q_p),
    M_p=the positive 2,3-free part of |R_p|.

**Theorem.** These are actual legal two-edit paths at coprime positive-clock shapes. For every a>=2,

    gcd(R_p,6)=1,    M_p=|R_p|,

and along this explicit infinite sequence,

    M_p^(1/p) -> 3.                                    (1.1)

Consequently, for any constants C>0 and 0<c<3, the proposed bound M_p<=C c^p fails for all sufficiently large members of this family. In particular there is no eventual uniform stripped-norm base below 3 covering all actual nearby positive-clock two-edit shapes.

There are three essential qualifications:

1. These q lie within six of p alpha but are not required to equal ceil(p alpha). The theorem does **not** establish the corresponding eventual barrier when restricted to exact critical clocks alone.
2. It bounds the available norm modulus, not the actual carry content g. A large norm can have little or no common factor with the clock.
3. Some or all members could be easily excluded individually because their clock D=2^q-3^p already exceeds M_p. The result obstructs this particular uniform norm shortcut, not every asymptotic two-edit argument or a smaller content bound.

No statement about the distribution of fractional parts of 3^a log_2 3 is used. The bounded adjustment of q is what makes the arithmetic family explicit.

## 2. The shapes and ordered paths are genuinely legal

The definition gives

    0<q-p alpha<6,    q=1 mod 6.

Since p is divisible by 6 and has no prime factors other than 2 and 3,

    r=1 mod 6,    gcd(p,r)=gcd(p,q)=1.

The elementary inequalities 3^2>2^3 and 3^5<2^8 give

    3/2<alpha<8/5.

For p>=18,

    p<q<p alpha+6<(8/5)p+6<2p,
    r>p/2.

Also n=2p-q is positive and congruent to 5 modulo 6. Therefore n>=5, so every nonconstant monomial in Q_p has positive exponent and its constant coefficient is exactly 3.

At k=1 the source gaps are 1 or 2. Consecutive source-one residues differ by -r modulo p. The follower residue A=p-1 has next-one residue

    B=p-1-r.

Because A>=r, their source gap is 1. Because 0<=B<r, the gap after the leader is 2. Thus they occur in a cyclic source block 110. Move the leader one place right, then the follower one place right:

    110 -> 101 -> 011.

Both swaps are actual and preserve labelled cyclic order, even if the chosen block crosses the original cut. The sparse-form formula gives

    1+(Z-1)(Z^(p-1)+Z^(p-1-r))
      =3-Z^(p-1)+Z^(p-r)-Z^(p-r-1)

modulo Z^p-2, exactly Q_p.

The two small indices are deliberately excluded: a=0 gives (p,q)=(2,7), and a=1 gives (6,13), both outside p<q<2p. They are not used as k=1 examples.

## 3. No factor 2 or 3 survives in the norm

### Modulo 2

Reduce Z^p-2 modulo 2 to Z^p. Since Q_p(0)=3 is odd,

    Res(Z^p-2,Q_p)=Res(Z^p,Q_p)=Q_p(0)^p=1 mod 2.

Hence R_p is odd.

### Modulo 3

In characteristic 3, with p=2*3^a,

    Z^p-2=(Z^2+1)^(3^a).

The two roots zeta of Z^2+1 lie in F_9 and satisfy zeta^2=-1. In particular they have order 4. On the other hand,

    (zeta-1)^2=zeta,
    (zeta-1)^4=-1,
    (zeta-1)^8=1.

Thus zeta-1 has order exactly 8. Every power zeta^r has order dividing 4, so

    zeta^r != zeta-1                                   (3.1)

for every integer r, with no distribution or congruence assumption on r needed here.

Reducing the sparse polynomial modulo 3 gives

    Q_p(Z)=Z^(n-1)(-Z^r+Z-1).

Neither root zeta is zero, and (3.1) says neither is a root of the second factor. Therefore Z^p-2 and Q_p have no common root over the algebraic closure of F_3. Their resultant is a nonzero element of F_3, proving 3 does not divide R_p. Repeated roots in Z^p-2 do not affect this nonvanishing argument.

Together with the parity calculation, this proves M_p=|R_p| exactly for every a>=2.

## 4. The norm root limit is 3

The analytic argument is supplied explicitly to show why bounded clock adjustment is harmless. In the field Z^p=2, multiply Q_p by Z:

    H_p(Z)=Z Q_p(Z)=3Z-2+(Z-1)Z^n.

Since the absolute norm of Z is 2,

    |Norm(H_p)|=2|R_p|.                                (4.1)

Let rho_p=2^(1/p) and beta_p=n/p. The bounded adjustment of q implies

    beta_p -> beta=2-alpha,
    rho_p^n -> 2^beta=4/3.

The real number beta is irrational: rationality of log_2 3 would imply equality of a positive power of 2 and a positive power of 3.

At the p conjugates Z_j=rho_p exp(2 pi i j/p), write theta_j=2 pi j/p and phi_j=n theta_j. The logarithmic integrands converge uniformly to the continuous function

    f(theta,phi)=log|3 exp(i theta)-2
                      +(4/3)(exp(i theta)-1)exp(i phi)|.

To justify uniformity, for every rho>=1,

    max_theta |rho exp(i theta)-1|/|3rho exp(i theta)-2|
       =(rho+1)/(3rho+2)<=2/5.

Since 0<n<p, rho_p^n<2. The perturbing term therefore has magnitude at most 4/5 of the linear term, uniformly in p and theta. No zero or logarithmic singularity approaches the integration domain. At the limit the sharper ratio bound is 8/15.

The p points (theta_j,phi_j) become uniformly distributed on the two-dimensional torus along this sequence. For any fixed integers (h,l) not both zero,

    (1/p) sum_j exp(i(h theta_j+l phi_j))
       =1 if p divides h+l n, and 0 otherwise.

When l!=0, divisibility for arbitrarily large p would force an integer sequence h/p+l beta_p to converge to the irrational number l beta, which is impossible. For l=0 and h!=0, divisibility also eventually fails. Trigonometric-polynomial approximation therefore gives convergence of averages for every continuous torus function, including f.

For fixed theta, the constant term 3 exp(i theta)-2 dominates the rotating term. Averaging f over phi gives

    log|3 exp(i theta)-2|.

For example, this follows by expanding the real part of log(1+w) as an absolutely convergent series: every nonconstant Fourier term has average zero. Averaging next over theta gives log 3 by the same argument, because 2/3<1. Thus

    (1/p) log|Norm(H_p)| -> log 3.

Equation (4.1) removes only log 2/p, and Section 3 removes no prime factors at all. This proves (1.1) for the stripped norm itself.

Finally, if M_p<=C c^p held eventually with c<3, then taking p-th roots would give limsup M_p^(1/p)<=c, contradicting (1.1).

## 5. Exact controls and the remaining frontier

`certify_stripped_barrier.py` checks F_9 directly as pairs over F_3, verifies both order statements and every exponent class modulo 8, and constructs all actual paths for a=2 through 8. It checks the coprime shapes, the ordered 110 chain, positive clocks, and the unchanged monomial degrees. It also computes exact symbolic resultants for a=2,3,4. It records 112 explicit checks. Normal and optimized runs agree byte-for-byte.

The first two norms are

| a | p | q | r | n | exact M_p |
|---|---:|---:|---:|---:|---:|
| 2 | 18 | 31 | 13 | 5 | 261006559 |
| 3 | 54 | 91 | 37 | 17 | 57871986057582534182735809 |

The first member has actual content 17; the other six tested members through a=8 have content 1. Neither observation is promoted to an infinite content claim.

The finite rows are controls, not evidence substituted for the infinite finite-field and torus proofs. The [independent audit](audit/INDEPENDENT_NEAR_CRITICAL_AUDIT.md) separately checks the argument and small-index boundaries.

Reproduce from this directory:

    python3 certify_stripped_barrier.py > receipt.replay.json
    python3 -O certify_stripped_barrier.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json

The exact-critical stripped-norm asymptotic question from the predecessor remains open. So do unbounded two-edit integer exclusion and a useful smaller bound on actual carry content. The successful change here is precisely to tune the nearby clock modulo 6, retain legal two-edit paths, and eliminate the problematic prime factors before invoking the norm limit.
