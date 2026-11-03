# A cubic norm recognizes the tribonacci guard chain

**Status.** A finite-root theorem for every positive integral cubic norm shell; a complete description of the nonnegative norm-two shell; an explicit identification with the transitive node's partner in the complement-line metagraph; exact modular periods at three primes and all their powers; and counterexamples to a one-orbit norm criterion in degrees four and five. No claim of literature priority is made.

The theorem classifies **endpoint triples**, not all tournaments with those triples. It neither solves a general Thue–Mahler equation nor identifies an unspecified earlier “305” or “Lucas-47” construction.

## 1. The physical source and the algebraic integer

Use [S164's one-vertex guard operation](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/05-knowledge/results/session164_emergent_20261103/profile_algebras/GUARD_ALGEBRA.md), Section 5. The guard q has unique in-neighbor r. Its state is

- s: Hamilton paths starting at q;
- u: Hamilton paths not starting at q;
- d: Hamilton paths of X−q starting at r.

Add a new vertex a with q→a and a→X−q, then mark a. The already-proved physical transfer is

    C(s,u,d) = (u,s+u+d,s),
    C^(-1)(s,u,d) = (d,s,u−s−d).                         (1)

This step preserves strength when its source is strong. The canonical triangle has state (1,2,0).

Let α be the real root of f(x)=x³−x²−x−1 and put

    z(s,u,d) = (α²−α)s + αu + d.

The coefficient change (s,u,d)↦(d,u−s,s) is unimodular, and

    z(Cv)=αz(v),       Norm(α)=1.

Thus the exact integral invariant is

    N(s,u,d) = Norm(z(s,u,d))
      = d³ + 2d²s + d²u + 2ds² − 2dsu − du²
        + 2s³ − 2su² + u³.                             (2)

In particular N(1,2,0)=2. The invariant is not inferred merely from matching characteristic polynomials: the displayed coefficient map intertwines the **literal physical transfer** with multiplication by α.

## 2. Every positive norm shell has finitely many descent roots

For a positive integer K define

    R_K = {(s,u,d)∈Z_≥0³ : N(s,u,d)=K, u<s+d}.

**Theorem 1.** Every nonnegative integral solution of N(v)=K has a unique representation

    v = C^j r,         r∈R_K, j≥0.                       (3)

Moreover, every r=(s,u,d)∈R_K satisfies

    s+u+d ≤ floor((1792K)^(1/3)).                         (4)

Consequently (3) is an effective, finite catalogue of recurrence orbits for every fixed K. An empty R_K proves the entire nonnegative norm shell empty.

**Proof.** There is exactly one real root α: f has critical points −1/3 and 1, and both critical values are negative. Its other roots are a conjugate pair β, β̄. Exact substitution gives

    11/6 < α < 15/8,
    β+β̄=1−α,       ββ̄=1/α,
    (Im β)² = 1/α−(α−1)²/4
             > 8/15−49/256 > 1/3.                     (5)

The polynomial f is irreducible over Q by the rational-root test. For a nonzero integral v, the quadratic polynomial defining z cannot vanish at β. Also z(v)>0 for nonzero nonnegative v. Hence N(v)=z(v)|z_β(v)|² is a positive integer.

Write L=s+u+d, Q=|z_β(v)|², and assume u≤s+d. Since α²−α>1 and α>1, z(v)≥L. Set

    E = Re z_β(v) − Im z_β(v)/(2 Im β)
      = ((2+2α−α²)/2)s − (α/2)u + d.

Using u≤s+d gives

    E ≥ ((2+α−α²)/2)s + (1−α/2)d
      ≥ (1−α/2)(s+d) > (s+d)/16 ≥ L/32.               (6)

Here the second inequality uses α(2−α)>0. By Cauchy–Schwarz and (5),

    E² ≤ (1+1/(4(Im β)²))Q < (7/4)Q.

Therefore every nonzero integral v in this cone satisfies

    N(v) > L³/1792.                                    (7)

This proves (4), with a harmless extra boundary candidate if 1792K is a cube.

If u≥s+d, the inverse in (1) is nonnegative and its coordinate sum is u≤L, strictly unless s=d=0. Equality can occur only if s=d=0; such a nonzero vector instead has inverse (0,0,u), whose sum equals u, so this edge case needs one more step. That inverse lies in the strict root cone and has the same norm. Thus inverse iteration terminates: ordinarily the sum decreases, and the only equal-sum move immediately reaches a root. The final triple belongs to R_K. Since C is invertible and has nonnegative entries, all forward iterates of every root are valid, and no two roots or their forward chains can meet. This proves existence and uniqueness in (3). ∎

The equal-sum edge case matters for general K, for example C(0,0,1)=(0,1,0). It does not occur on the norm-two shell, because u³=2 has no integral solution.

## 3. The norm-two shell is exactly one chain

**Theorem 2.**

    {v∈Z_≥0³ : N(v)=2} = {C^j(1,0,0) : j≥0}.            (8)

The initial states are

    (1,0,0), (0,1,1), (1,2,0), (2,3,1), (3,6,2),
    (6,11,3), (11,20,6), (20,37,11), …

**Proof.** Theorem 1 reduces all possible roots to L≤15, because 15³<3584<16³. There are 816 nonnegative triples with L≤15; exactly 576 also satisfy u<s+d. Complete exact enumeration of these 576 candidates gives

    R_2 = {(1,0,0)}.

This finite enumeration is part of the proof, not a heuristic control. It is fully reproduced in `check.py`. A compact independent check is the list of minimum norms in the strict root cone at each coordinate sum L=1,…,15:

    1, 7, 4, 19, 11, 32, 17, 53, 26, 73, 44, 98, 53, 136, 76.

At L=1 the two possible root triples are (1,0,0), of norm 2, and (0,0,1), of norm 1. Every L≥2 in this list has minimum greater than 2. ∎

The same finite-root algorithm gives the complete catalogues for K≤10:

| K | R_K |
|---|---|
| 1 | {(0,0,1)} |
| 2 | {(1,0,0)} |
| 3 | empty |
| 4 | {(1,1,1)} |
| 5,6 | empty |
| 7 | {(1,0,1)} |
| 8 | {(0,0,2)} |
| 9,10 | empty |

These are arithmetic catalogues, not claims that every listed root is itself a physical tournament interface.

**Scaled recognition.** For nonzero v∈Z_≥0³, let g=gcd(s,u,d). Then

    N(v)=2g³  ⇔  v=gC^j(1,0,0) for a unique j≥0.         (9)

Indeed N is cubic, C is unimodular and hence preserves gcd, and Theorem 2 applies to v/g. Thus the normalized norm is a complete recognition test for this scaled family of endpoint states. It is not a complete recognition test for the underlying graphs, nor a proof that every scaled state is physically attainable by an authorized construction.

## 4. The literal complement-line metagraph bridge

Let T*_n have vertex set {0,…,n−1}, arcs i→i+1, and j→i whenever j≥i+2. This is the repository's all-backward full tiling tournament. Its top vertex q=n−1 has unique in-neighbor n−2. Adding the next guard by C creates **exactly** T*_(n+1), with no relabeling needed.

The directed edge T*_2 has state (0,1,1), and T*_3 is the canonical triangle. Therefore

    state(T*_n, n−1)=C^(n−1)(1,0,0),       n≥2.          (10)

For n≥3 these graphs are strong, as witnessed by the directed cycle 0→1→…→n−1→0. The nonnegative norm-two shell consists of the formal initial state (1,0,0), the edge state, and this entire strong physical tail. The formal first state is not being declared a valid unique-in-neighbor guarded singleton.

In the complement-line metagraph, the transitive class has the all-forward tiling as its unique presentation; complementing the non-path tiles gives exactly T*_n. Thus (10) attaches the norm theorem to an actual all-order metagraph node family. It is not a construction in the different single-arc-flip graph.

The prior results must be kept separate from the new norm recognition. [THM-4731, Statement 6](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md) already records the earlier Hamilton-count law and proves rigidity for n≥4, self-converseness, and the blue multiplicity law. With h_1=h_2=1, h_3=3 and h_n=h_(n−1)+h_(n−2)+h_(n−3), these are

    H(T*_n)=h_n,
    tc(T*_n)=h_n                    (n≥4),
    bm(T*_n)=h_(floor(n/2)+1)        (n≥4).

At n=3 the triangle has nontrivial automorphisms, so tc=H is not asserted. The new statements here are the explicit integral norm, the finite-root theorem, and the complete norm-two endpoint-shell classification, together with their bridge to that already-known node family.

## 5. Exact modular clocks, including the split prime 47

Let M be the matrix of C. The initial triangle has norm 2, so its algebraic image is a unit modulo every odd modulus. Therefore its state period equals the multiplicative order of α, equivalently the order of M.

The scalar count H=s+u has the same period modulo every odd modulus, since the three successive observations have matrix

    [1 1 0]
    [1 2 1],       determinant = 2.                     (11)
    [3 3 2]

Thus three consecutive counts recover the whole state modulo any odd modulus.

**Theorem 3.** For k≥1, the exact periods of both the canonical C-state sequence and its Hamilton-count sequence are

    13·3^(k−1)    modulo 3^k,
    31·5^(k−1)    modulo 5^k,
    46·47^(k−1)   modulo 47^k.                          (12)

**Proof.** Direct exact matrix powers show that the respective orders modulo 3, 5 and 47 are 13, 31 and 46. For 47 there is also the explicit factorization

    f(x)=(x−5)(x−17)(x−26) mod 47,

with root orders 46,23,46. The first lift matrices (M^t−I)/p modulo p, for the indicated base order t, are

    p=3:  [1 0 0]       p=5:  [3 1 0]       p=47: [ 0 8  8]
          [0 1 0]             [1 4 1]             [16 8  8]
          [0 0 1]             [1 0 3]             [ 0 8 39].

Each is nonzero. For odd p, if X=I+p^a A with A nonzero modulo p, the binomial expansion gives X^p=I+p^(a+1)A modulo p^(a+2). Consequently the order multiplies by p at each lift, proving (12). Equation (11) transfers the state period to H. ∎

There is a useful warning about finite norm quotients. Modulo 3 and 5, the norm-two shell has 13 and 31 points, respectively, so in each case it is one C-orbit. Modulo 47 it has 46²=2116 points and splits into **46 C-orbits of length 46**. Indeed the three root evaluations give an isomorphism to F_47³; their product must equal 2, so any two nonzero coordinates determine the third. Every point is a unit and hence has the same period 46. The integer nonnegative norm-two shell is one chain, but reducing its norm equation modulo a prime need not retain that uniqueness. The role of 47 in this statement is its displayed splitting, not an inferred identification with a different Lucas construction.

## 6. Why degrees four and five need more than the norm value

For k=4,5 put f_k(x)=x^k−(1+x+…+x^(k−1)), and work in the free integer algebra Z[x]/(f_k). Its norm means the determinant of multiplication in the basis 1,x,…,x^(k−1). Multiplication by x is the usual nonnegative k-step companion transfer. These are algebraic recurrence examples; no new physical guard realization for k=4 or 5 is being asserted.

There are strictly positive coefficient vectors of equal norm which are not on the same companion orbit:

| k | First coefficient vector | Second coefficient vector | Common norm |
|---|---|---|---|
| 4 | (1,1,2,2) | (1,2,2,1) | 3 |
| 5 | (1,1,1,1,2) | (1,1,1,2,1) | 2 |

Here is an exact certificate of the separation, not a finite failure to find a connecting power.

For k=4 let p=1+x+x² and ε=1+x. Then

    first=x²p,       second=εp,
    Norm(p)=3,       Norm(ε)=1,
    ε^(-1)=2−x+2x²−x³.

At the positive real root α_4, which lies between 19/10 and 2,

    α_4 < 1+α_4 < α_4².

Thus ε(α_4) is not any integral power of α_4. Because p is invertible over Q (its multiplication determinant is 3), an equality connecting the two vectors by any integral power of x would imply precisely such a power identity, which is impossible.

For k=5 use p=1+x and ε=1+x²+x³. Then

    first=x⁴p,       second=εp,
    Norm(p)=2,       Norm(ε)=1,
    ε^(-1)=1+3x+4x²+x³−2x⁴.

At the positive real root α_5∈(19/10,2),

    α_5³ < 1+α_5²+α_5³ < α_5⁴.

For the upper inequality, x⁴−x³−x²−1 is positive at 19/10 and increasing for x≥19/10. Again ε is not an integral power of x, and cancellation of p proves distinct orbits. This gives, in particular, a strictly positive **quintic norm-two counterexample** to importing Theorem 2 unchanged into degree five. In degree four the norm changes sign under one companion step; the displayed equal positive norms remain separated even when arbitrary integral powers are allowed.

A norm collapse is therefore a strong source of necessary arithmetic structure, but a one-orbit conclusion needs its own argument. The cubic theorem supplies that argument through its positive inverse cone and finite root catalogue.

## Reproduction and proof boundaries

Run `python3 check.py` and `python3 -O check.py`. The outputs are identical; tests use explicit exceptions rather than removable assertions.

The checker independently computes multiplication determinants to certify (2), checks the polynomial C-invariance on a degree-determining grid, exhausts all 576 norm-two descent candidates, gives all root catalogues K≤10 using the proven bound, and independently counts Hamilton paths by subset/endpoint DP for T*_n through n=13. It also certifies the three base periods, first lift witnesses, norm-fibre sizes, and the higher-degree norm/inverse identities. Prime-power controls through k=3 support the all-k proof; they do not replace it.

The finite catalogue in Theorem 2 and its bound are part of the exact proof. The physical DP is an independent control for the general edge-by-edge construction. No endpoint classification here is a classification of all graphs, full path-cover profiles, all metagraph nodes, or all exact odd Hamilton counts.
