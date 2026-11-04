# A smaller affine profile for the full two-cover response

## 1. Physical setting and normalized profiles

Let Q be an arbitrary tournament on m nonempty child blocks A_i. Attach actual ports x->y with two prescribed incidence bits constant on each whole child. Write F_h for the complete number of ordered covers of this augmented tournament by h nonempty directed paths. In particular F_1 is the ordinary physical Hamilton response and F_2 counts all ordered two-path covers, with every port position retained.

The five-profile theorem supplies the normalized child coordinates modulo eight

    r_i = q(A_i)/H(A_i)
        = (0,1,2-2s_i+4b_i,2c_i+4b_i,0,...),

where S=(s,b,c)=(1,0,0), P=(0,0,0), E=(0,1,0), O=(0,1,1), C=(0,0,1). These stand for singleton, pair, even order at least four, odd order with even Hamilton-cycle parity, and odd order with odd Hamilton-cycle parity. Every H(A_i) is odd, and its residue is retained as an actual multiplicative factor. These normalized vectors are not all literal tournament profiles.

## 2. Every two-cover run coefficient is even

For a positive child run multiplicity vector k, let B(k) be the number of ordered pairs (W_1,W_2) of nonempty directed run words in the augmented quotient. Each port occurs exactly once across the pair, child i occurs k_i times, equal adjacent child labels are forbidden WITHIN each word, and every within-word adjacency is an actual directed arc. Equal child labels can occur on opposite sides of the path boundary. There is no arc requirement across that boundary.

An actual ordered global two-cover determines this pair of run words. For each child, reading its run slots first along W_1 and then along W_2 supplies one ordered child path cover. Conversely these choices reconstruct the actual two-cover uniquely. Thus

    F_2 = product_i H(A_i) · sum_k B(k) product_i r_(i,k_i) (mod8).  (1)

The involution (W_1,W_2) -> (W_2,W_1) preserves all defining conditions. A fixed pair would have W_1=W_2. This is impossible: the port x occurs exactly once in the pair, whereas two identical words would contain either zero or two copies of x. Hence EVERY integer B(k) is even, independently of k, the quotient, and the port incidences. This is a statement about abstract run coefficients, not merely the familiar evenness of the final ordered-cover count.

## 3. Exact affine formula

Let B_0=B(1,...,1), B_i=B(i2), and C_i=B(i3), with unlisted children visited once. All three kinds of coefficients are even. Child coordinates of index at least four vanish modulo eight. A term with at least two repeated children has a factor four from the two child profiles and another factor two from B(k), so vanishes. The 4b_i corrections also vanish against every even coefficient. Equation(1) therefore becomes

    F_2 = product_i H(A_i) · [ B_0
           +2 sum_i (1-s_i) B_i
           +2 sum_i c_i C_i ]                 (mod8).             (2)

This formula is separately valid for EACH full physical two-cover response. There is no dependence on any b_i and no quadratic or mixed type term. For nontrivial children it is affine in just the m odd-cycle bits c_i. Quotient-converse contrasts inherit the same formula whenever corresponding child profiles, and hence the odd scalar product, are retained.

The type identifications P/E and O/C are not literal equalities of entire profiles. Equation(2) says that changing b_i while keeping s_i and c_i fixed is invisible to this normalized two-cover functional modulo eight. In the five-state notation this identifies P with E, and C with O.

## 4. Exact extraction from small actual model graphs

At a fixed quotient and fixed block query, let G_0 be the actual two-cover response when all children are singletons. Let G_i^(2) and G_i^(3) replace only child i by a transitive pair or triple. Their exact child profiles give

    B_0=G_0,
    B_i=(G_i^(2)-G_0)/2,
    C_i=(G_i^(3)-3G_i^(2)+2G_0)/6.                         (3)

The first numerator is divisible by four and the second by twelve, since B_i and C_i are even. Thus 1+2m complete physical two-cover counts on known models with at most m+4 augmented vertices determine (2), for arbitrary original child sizes. These are offline counts on known model graphs, not oracle-query counts on an unknown original tournament.

The scalar product of actual H(A_i) residues must be restored. For canonical C3 children that factor is three, and omitting it would give the wrong individual response residues.

## 5. Scope, combination and boundary

This is an all-order structural compression of the two-cover functional, not a non-annihilation or reconstruction theorem. It does not assert that any quotient-converse contrast is nonzero, locally or after placing the core in an ambient. Exterior vertices may be incorporated as singleton quotient children without changing the formula.

The complete Hamilton functional has the separate, richer formula in coarse_mod8_affine/PROOF.md: it is affine in 2m nontrivial type bits and has singleton-pair quadratic terms. The present reduction uses a different involution and does not infer its coefficient parity from that Hamilton formula.

For same-order actual tournaments the previously established modulo-eight profile identities give Delta q3=Delta q2 in odd order and Delta q3=0 in even order; every higher coordinate vanishes. Consequently the present two-cover compression can be used together with a separately certified nonzero Hamilton contrast to CHECK an all-cover ambient certificate. It does not supply that nonzero contrast or a vanishing two-cover contrast by itself.

The two-path exchange proof does not generalize to arbitrary h by pretending that all h! permutations act freely on abstract run words. Identical port-free path words can occur. A concrete boundary has one transitive child, both port incidence bits zero, and exact three-cover model values whose extracted triple-run coefficient is 9, odd. The verifier retains the actual profiles and an explicit fixed three-cover word pattern. Only the h=2 statement is asserted here.
