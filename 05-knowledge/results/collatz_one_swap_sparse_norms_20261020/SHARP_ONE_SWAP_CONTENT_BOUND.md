# A sharp all-length content bound for one-swap mechanical Collatz words

Date: 20 October 2026. Status: **PROVED elementary reduction, sharp uniform bound, exact endpoint families, and periodicity; FINITE-EXACT symbolic and word audits.** The all-weight positive one-swap exclusion remains **OPEN** at a single critical length per weight. No Collatz convergence, general positive-cycle exclusion, or historical-priority claim is made.

## 1. Inheritance and the question that actually remains

For a chronological binary word w of length q with p ones, write

    B(w) = sum_(i:w_i=1) 2^i 3^(number of ones after i),
    D(q,p) = 2^q - 3^p,
    g(w) = gcd(|D|, B(w)).

The shortcut branch is (3^p x+B)/2^q. Its fixed point is B/D, and its reduced denominator is |D|/g. Integer-cycle recognition is g=|D|; positivity additionally requires D>0. Rotation preserves g. These facts, including the full prime-power criterion and the chronological-edge recognition identity, are inherited from [Local carry and prime-power gates](../collatz_local_prime_power_gates_20261020/LOCAL_CARRY_AND_PRIME_POWER_GATES.md) and the canonical predecessor `05-knowledge/results/collatz_rotation_defect_ideal_20261020/LOCAL_ROTATION_DEFECT_IDEAL.md`.

The balanced-word mechanism is THM-4716, `01-canon/theorems/THM-4716-balanced-collatz-cycles-are-three-smooth-ladders-and-are-free.md`. Its freeness theorem is not extended to one-swap words: the canonical hostile `0000101` has p=2, q=7, and g=7. The newer hostile below has p=3 and g=25, so the full valuation vector remains essential.

A necessary preliminary correction to the original fixed-small-p target is that **every** positive integer cycle, not just a one-swap word, satisfies q<=2p. If its positive odd entries are n_1,...,n_p and the gaps between odd steps are a_1,...,a_p, then

    2^q = product_j (3+1/n_j) <= 4^p.

Consequently the inherited q<=400 exact one-swap census already gives an all-length exclusion for p<=200. Rechecking p=3,4,5 would not be substantive new exclusion progress. This report instead proves a uniform content bound, resolves all-length prime-power content classes at small weights, and narrows the all-weight obstruction to q=ceil(p log_2 3).

Concept board: mechanical odd-step gaps; cut-safe local edits; algebraic norms; full content valuations; arithmetic periodicity; positive critical clocks. The anchor is a uniform bound, the niche is exact periodic content, and the wildcard is the sparse norm representation. The map and its losses are explicit below: one-swap words map to a polynomial family and then a sparse norm; the norm preserves an upper bound in divisibility but loses k and the common modular root. Retaining D and P restores the exact gcd.

## 2. A polynomial family covering every cyclic adjacent swap

Assume p>=2, q>p, and gcd(p,q)=1. Write

    q = p k + r,    k>=1,    1<=r<p,    gcd(p,r)=1.

Let c be the zero-intercept lower mechanical word

    c_i = floor((i+1)p/q)-floor(ip/q),    0<=i<q.

Index its ones cyclically by j=0,...,p-1, with j=0 the one at position q-1. Their positions are

    s_j = j k + ceil(jr/p)-1 modulo q.

Extend j to all integers. The gap from one j-1 to one j is

    k + delta_j,    delta_j=ceil(jr/p)-ceil((j-1)r/p) in {0,1}.

Every unequal cyclic adjacent swap moves exactly one indexed one by epsilon=+1 or -1. Anchor the rotated word at its preceding one, indexed t, so the moved one is t+1. That anchor is unchanged, including when the original swap crosses positions q-1 and 0. The offsets of the other ones before the edit are

    m k + eta_m,    eta_m=ceil((t+m)r/p)-ceil(tr/p),    0<=m<p.

Only the m=1 offset changes, by epsilon. The swap is allowed exactly when the destination is zero:

- epsilon=-1: k+delta_(t+1)>=2;
- epsilon=+1: k+delta_(t+2)>=2.

All types are allowed for k>=2. At k=1 these conditions must be retained; a formal polynomial with a forbidden move is not a binary-word instance.

Put

    h=max(0,-eta_1-epsilon),
    P_(p,r,t,epsilon)(X)
      = sum_(m=0)^(p-1) 3^(p-1-m)
          2^(h+eta_m+epsilon*1_(m=1)) X^m,
    F_(p,r)(X)=2^r X^p-3^p.

Here h is 0 or 1, P has integer coefficients and degree p-1, and

    P(2^k) = 2^h B(w rotated to the unchanged anchor),
    g(w) = gcd(|F(2^k)|, P(2^k)).                         (2.1)

This parameterization covers each swap occurrence, not just distinct output necklaces. It explicitly includes crossing-cut edits and swaps returning to a balanced rotation. Rotation and the factor 2^h preserve content because D is odd.

## 3. The universal sparse norm, independent of r

Let z be an algebraic root of z^p=2. The polynomial Z^p-2 is Eisenstein at 2, hence irreducible. Put x=3 z^(-r). Then F(x)=0. Because gcd(p,r)=1, integers u,v exist with pu+rv=1, giving

    z=2^u (3/x)^v.

Thus Q(x)=Q(z), and F is irreducible over Q. In particular the resultant of F with the nonzero degree-(p-1) polynomial P cannot vanish. There is no hidden zero-resultant or cyclotomic exception under the stated coprimality hypotheses.

Define

    a_j = p ceil(jr/p)-jr = (-jr mod p) in {0,...,p-1},
    A=a_t,    a=a_(t+1).

The a_j run through every residue 0,...,p-1. At x=3 z^(-r),

    P(x) = 2^h 3^(p-1) z^(-A)
            [S(z)+(2^epsilon-1) z^a],
    S(z)=1+z+...+z^(p-1)=1/(z-1).                       (3.1)

The last equality uses z^p=2. Multiplying by z-1 gives, up to powers of 2 and a monomial, the two universal polynomials

    Q_(p,a,+)(Z)=Z^(a+1)-Z^a+1,
    Q_(p,a,-)(Z)=Z^(p-a)-Z+1.                          (3.2)

For the minus case the exact relation is

    2-Z^(a+1)+Z^a = Z^a [Z^(p-a)-Z+1] modulo Z^p-2.

Let

    R_(p,a,epsilon)=Res_Z(Z^p-2,Q_(p,a,epsilon)),
    M_(p,a,epsilon)=the positive 2,3-free part of |R|.  (3.3)

None of these resultants is zero. For interior types Q has degree <p. The endpoints reduce to Z, 3-Z^(p-1), 3-Z, or 1, which likewise cannot vanish at a root of Z^p-2.

**Theorem A (uniform norm gate).** For every allowed one-swap word above,

    g(w) divides M_(p,a,epsilon).                       (3.4)

The bound depends only on p, the residue a of the moved one, and its direction. Changing r merely permutes the p edit labels a.

There are two complete proofs, both preserving prime powers.

**Resultant proof.** Taking norms in (3.1), and using |Norm(z-1)|=1 and |Norm z|=2, gives the exact identities

    |Res_X(F,P)|
      = 2^[r(p-1)-A] 3^[p(p-1)] |R_(p,a,+)|   (epsilon=+1),
    |Res_X(F,P)|
      = 2^[r(p-2)]   3^[p(p-1)] |R_(p,a,-)|   (epsilon=-1).

For the second identity, eta_1=(r+a-A)/p and h=1-eta_1. Both displayed powers of 2 have nonnegative exponents. An integer Bezout identity from the Sylvester determinant puts Res(F,P) in the ideal (F,P). Evaluation at 2^k shows that g divides this resultant. Since gcd(D,6)=1, its 2,3 factors can be removed, proving (3.4).

**Modular proof.** Work modulo g. Since 2 and 3 are units, define

    z=2^u(3/2^k)^v,    pu+rv=1.

The clock relation gives z^p=2 and z^r=3/2^k. Also z-1 is a unit: (z-1)S(z)=1. The exact carry identity (3.1) therefore gives Q(z)=0. The integer resultant identity for Z^p-2 and Q again implies g divides R. This is an argument in Z/gZ, not an argument restricted to fields, so it retains every prime-power depth.

The irreducibility hypothesis must not be silently removed from a generic fixed-weight resultant argument. For example, at p=2,r=0, the unswapped repeated balanced family has F=X^2-9 and P=X+3; its resultant is zero because X+3 is a common factor. That example is outside this theorem and is not claimed to be a one-swap counterexample.

## 4. A sharp bound valid for every p and every length

**Theorem B (sharp uniform content bound).** Under the hypotheses above,

    g(w) <= 3^p-2.                                      (4.1)

For each p>=2, equality is attained by infinitely many allowed one-swap words of coprime shape. Thus no smaller bound depending only on p can replace the right side.

**Proof of the upper bound.** The p roots of Z^p-2 are z_j=rho exp(2 pi i j/p), with rho=2^(1/p). For plus types 1<=a<=p-2, the exponents 0,a,a+1 are distinct modulo p. Orthogonality gives

    (1/p) sum_j |Q_+(z_j)|^2
       =1+rho^(2a)+rho^(2a+2)<9.

The arithmetic-geometric mean inequality therefore yields

    |R_+|=product_j |Q_+(z_j)|<3^p.

For minus types 1<=a<=p-2, the exponents 0,1,p-a are distinct modulo p, and the same argument uses

    (1/p) sum_j |Q_-(z_j)|^2
       =1+rho^2+rho^(2(p-a))<9.

The four endpoint values are explicit:

    |R_(p,0,+)|=2,
    |R_(p,p-1,+)|=3^p-2^(p-1),
    |R_(p,0,-)|=3^p-2,
    |R_(p,p-1,-)|=1.

All are less than 3^p. There are no interior types when p=2. Consequently M<3^p for every type. M is a positive odd integer, so M<=3^p-2; (3.4) proves the bound. The construction achieving equality is given in Section 5. ∎

This genuinely uses the one-edit hypothesis. The primitive two-one word with smaller cyclic distance 17 and length 55 has content 49 by the inherited two-one formula, whereas 3^2-2=7. It is not a one-swap perturbation of the balanced necklace. Repetition is another hostile: (01)^3 has p=3,q=6,g=37>25 and is outside the coprime-shape hypothesis.

**Corollary B1 (single critical length).** If such a word encodes a positive integer 3x+1 cycle, then

    q=ceil(p log_2 3).                                  (4.2)

Indeed 0<D=g<=3^p-2 gives

    3^p < 2^q <= 2*3^p-2 < 2*3^p,

so p log_2 3 < q < p log_2 3+1. The logarithm is irrational, since 2^q=3^p is impossible for positive integers p,q. There is exactly one integer in this interval.

Equivalently, the word is excluded whenever 2^(q-1)>3^p, a statement checkable using integers alone. It is also excluded whenever the unique critical length is not coprime to p. The all-weight problem remains open at the coprime critical lengths; (4.2) does not itself make their content proper.

**Corollary B2 (inherited finite census upgraded correctly).** The existing exact census of every cyclic adjacent swap with q<=400 now implies all-length positive-cycle exclusion for every 2<=p<=252, because 3^252<2^400 while 3^253>2^400. This uses the new uniform theorem plus the old finite census, not a new all-q enumeration. At p=253 the critical length is 401 and gcd(253,401)=1, so the old census alone stops there. No general Collatz conclusion follows.

## 5. Two exact endpoint classes and sharpness

The large endpoint norm bounds are accompanied by exact gcd formulas, not just bounds.

**Theorem C (endpoint content).** At the endpoint a=0,epsilon=-1,

    g(w)=gcd(3^p-2, 2^k 3^(r-1)-1).                   (5.1)

At a=p-1,epsilon=+1,

    g(w)=gcd(3^p-2^(p-1), 2^(k+r)-3^(r+1)).           (5.2)

**Proof.** For any positive modulus m coprime to 6, simultaneous vanishing of F(2^k) and P(2^k) is equivalent to the modular root equations in Section 3.

In the first case Q=3-z, so the root is z=3. Its two defining equations become

    3^p=2,    2^k 3^(r-1)=1 modulo m.

Conversely these two equations imply the clock equation and force the constructed z to equal 3: its p-th and r-th powers agree with those of 3, and pu+rv=1 recovers z itself. Thus the sets of common divisors agree, including all prime powers, proving (5.1).

In the second case Q reduces to 3-z^(p-1). Together with z^p=2 this forces z=2/3. The two defining equations are then

    3^p=2^(p-1),    2^(k+r)=3^(r+1) modulo m.

The same converse proves (5.2). Each modulus on the right is already coprime to 6. ∎

For sharpness in Theorem B, take r=1 and the first endpoint. Let L be the multiplicative order of 2 modulo 3^p-2. Any positive multiple k of L with k>=2 gives

    q=pk+1,    gcd(p,q)=1,
    g=gcd(3^p-2,2^k-1)=3^p-2.

The move is admissible because k>=2. Infinitely many such k exist. These sharpness words generally have enormous positive D and a proper content; they are not integer cycles.

Small sharpness examples, obtained by exact arithmetic:

| p | k | q | g=3^p-2 |
|---|---|---|---|
| 2 | 3 | 7 | 7 |
| 3 | 20 | 61 | 25 |
| 4 | 39 | 157 | 79 |
| 5 | 24 | 121 | 241 |
| 6 | 121 | 727 | 727 |
| 7 | 396 | 2773 | 2185 |
| 8 | 117 | 937 | 6559 |

In particular the content is not uniformly bounded as p varies, and squarefree-only tests do not suffice even at fixed p=3.

## 6. Exact content periodicity, with its missing sidecar retained

**Theorem D (periodic all-length content).** Fix p,r,t,epsilon and put M=M_(p,a,epsilon). Then for every allowed k,

    g_k=gcd(M, 2^(pk+r)-3^p, P(2^k)).                  (6.1)

If M=1, then g_k=1 and take period 1. If M>1, the formal gcd sequence has period dividing ord_M(2). Restrict it to the admissible k values when describing actual words. In particular every k>=2 is admissible.

The equality follows from g|M and (2.1). Both polynomial evaluations depend only on 2^k modulo M, proving the period. This is exact for the complete gcd, not just its prime support. The order is a guaranteed period and need not be the least period. `catalogue.json` records the least period of the finite residue sequence for all 68 parameter classes with p=3,4,5; those minimality claims are exact finite calculations, not additional symbolic theorems.

A resultant upper bound can be strictly larger than every attained content. For example, for p=4,r=1,a=3,epsilon=+1, M=73 but (5.2) gives

    g=gcd(73,2^(k+1)-9)=1 for every k,

as the nine powers of 2 modulo 73 omit 9. Likewise p=5,a=2,epsilon=+1 has M=35, but the complete period tables attain only g=1 or 5. The missing factor 7 is not restored by retaining the prime support of M; the same k and the clock must be retained.

For p=3 the only potentially nonunit types are:

| r | type (a,epsilon) | exact nonunit classes of k |
|---|---|---|
| 1 | (2,+1) | g=23 for k=4 mod 11 |
| 2 | (2,+1) | g=23 for k=0 mod 11 |
| 1 | (0,-1) | g=25 for k=0 mod 20; g=5 for k=4,8,12,16 mod 20 |
| 2 | (0,-1) | g=25 for k=13 mod 20; g=5 for k=1,5,9,17 mod 20 |

Every other residue and every other p=3 type has g=1. Apply the k=1 admissibility test before interpreting a table entry as a word. As a concrete prime-square hostile, p=3,q=41,r=2,k=13, with ones at positions 13,27,39, has g=25. It is the zero-intercept mechanical word with its last one moved left once; the original ones were at 13,27,40.

## 7. Verification receipts and reproducibility

The main verifier imports no repository mathematical code. It independently computes carries from positions and affine composition, generates each actual cyclic swap, anchors at the preceding unchanged one, and checks the polynomial identity, type bound, exact gcd cap, and both endpoint formulas.

`receipt.json` records:

- 432 exact symbolic comparisons: 2<=p<=10, every coprime 1<=r<p, every preceding-one anchor and both directions. The ordinary carry resultant and universal sparse resultant have identical 2,3-free parts.
- 154 universal norm values: every a and both directions for 2<=p<=12.
- 15,096 actual swap occurrences: 2<=p<=12, p<q<=180, gcd(p,q)=1, including 1,120 original-cut crossing swaps, 374 admissible k=1 swaps, and 2,240 balanced-return occurrences.
- 68 complete periodic classes for p=3,4,5, with 2,690 exact two-period replays.
- The sharpness examples above, including lengths far beyond the direct finite universe.

Run from this directory:

    python3 verify_resultant_content.py --catalogue catalogue.replay.json > receipt.replay.json
    python3 -O verify_resultant_content.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json
    cmp catalogue.json catalogue.replay.json

An independent audit in `audit/` checks a different parameter universe and the modular proof. Its receipt and report are delivered separately alongside this note. These finite checks corroborate the symbolic proof; they do not supply its all-p or all-length quantifiers.

## 8. Honest remaining frontier

1. The universal one-swap positive-cycle exclusion is open at q=ceil(p log_2 3) when gcd(p,q)=1. The theorem eliminates every other positive clock length, not this remaining infinite family.
2. The exact norm constants lose the clock and the choice of common modular root. Keeping D and P in (6.1), or an equivalent common-root sidecar, is essential.
3. The sharp bound concerns one move of one odd position by one slot. It is not asserted for two or more swaps, a bounded number of arbitrary gap changes, noncoprime shapes, or arbitrary words.
4. Fixed-p periodicity is a full-content theorem. It neither proves eventual trajectory entry nor imposes a finite upper bound on p.
5. The finite all-length p<=252 consequence inherits the predecessor's exact q<=400 universe and its primitive-mechanical one-swap restriction. It is not a census of all parity necklaces.
