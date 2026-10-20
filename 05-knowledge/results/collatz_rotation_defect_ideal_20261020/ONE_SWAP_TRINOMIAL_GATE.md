# A sparse three-term gate for one-swap perturbations of balanced necklaces

Date: 20 October 2026.
Status: PROVED elementary reduction and explicit hostile; the proposed family exclusion remains OPEN. This supplements the local rotation-defect ideal theorem without changing its proof.

## 1. The false shortcut and its strongest useful survivor

The tempting statement “one adjacent swap from a primitive balanced word is still free” fails already at

    c=0001001, q=7, p=2,
    w=0000101, obtained by swapping positions 3 and 4 (zero-based),
    D=119, B(w)=112, gcd(D,B(w))=7.

Thus the reduced denominator is 17, not 119. The word is still excluded as a 3x+1 integer cycle, because 7<119. Freeness is stronger than the actual integer-cycle exclusion goal.

A cheap hostile probe covered all cyclic adjacent swaps of every lower mechanical word with coprime q,p and 3<=q<=80: 52,792 directed swap occurrences, including repetitions of the same output necklace. There were 3,821 occurrences with nonunit content. The maximum observed content was 11,243, at (q,p)=(66,35). This is a finite probe, not a census of distinct necklaces or evidence of a uniform bound.

The useful unresolved target is therefore: show that the defect content is a proper divisor of D for a precisely specified unbalanced one-swap family with D>1. The local arithmetic below produces a much narrower object to study than a general orbit.

## 2. Exact perturbation formula, with the crossing cut retained

Let c be the zero-intercept lower mechanical word of coprime shape (q,p), with q>=4 and 0<p<q, explicitly

    c_k = floor((k+1)p/q) − floor(kp/q),   0<=k<q.

Obtain w by swapping unequal adjacent cyclic bits at positions i and i+1 modulo q. Put epsilon=+1 if the swap changes 10 to 01, and epsilon=−1 for 01 to 10.

For each rotation cut r, let c_r and w_r be the rotations starting at r. There is exactly one exceptional cut,

    r_star=i+1 modulo q,

at which the swapped positions lie at the two ends of the linear word. Retain it as exceptional; do not treat it as an adjacent linear swap.

For every other r put j_r=i−r modulo q, with 0<=j_r<=q−2, and let s_r count the ones strictly after position j_r+1 in c_r. The local swap identity gives

    B(w_r)=B(c_r)+E_r,
    E_r=epsilon 2^j_r 3^s_r.                         (2.1)

The two words have unchanged total counts. All statements apply independently of whether w itself is balanced.

## 3. Intercept neighbors become explicit trinomials

Order balanced rotations by intercept t=0,...,q−1. The cut r_t is t p^(-1) modulo q. The intercept-shift lemma (the elementary ingredient of THM-4716) gives

    B(c_(r_t))−B(c_(r_(t+1)))=M_t=2^k_t 3^h_t,
    k_t = [(q−1−t)p^(-1) mod q]−1,
    h_t = p−1−floor((k_t+1)p/q),                     (3.1)

for 0<=t<=q−2. The residue in brackets lies in 1,...,q−1.

Let t_star=(i+1)p modulo q. For every edge t→t+1 whose endpoints avoid t_star, combine (2.1) and (3.1):

    B(w_(r_t))−B(w_(r_(t+1)))
      = M_t + E_(r_t) − E_(r_(t+1)).                 (3.2)

There are at least q−3 such edges. Each is a signed sum of three 2–3-smooth monomials, with fully prescribed exponents and signs. None is zero, by carry injectivity and primitivity (coprime q,p implies every output word of that shape is primitive).

Put D=2^q−3^p and g=gcd(|D|,B(w)). Then

    g divides every trinomial in (3.2).

Consequently

    gcd(|D|, all retained trinomials)

is an explicit upper bound in divisibility for g. If it is a proper divisor of |D|, w cannot be an integer 3x+1 cycle. If it is 1, w is free. The exceptional vertex and any disconnectedness matter: the retained edges alone are not claimed to recover g exactly. Adding enough other local-defect edges to make a connected graph invokes the full local-ideal theorem and restores equality.

For a hypothetical integer cycle, D must divide all these trinomials. Since 2 and 3 are invertible modulo D, each may be normalized into a three-term equation in powers of 2 and 3 modulo D. These linked equations, rather than a visual circle or an unconstrained trace, are the next arithmetic target.

## 4. Boundary and decisive next test

The all-near-balanced freeness claim is refuted above. The proper-divisor target for D>1 is not established here. The reduction also says nothing about integer sources with aperiodic parity sequences or about entry into the previously certified return bank.

A productive next proof would identify two of the prescribed trinomials whose resultant or Bezout combination is a 2–3-unit, or is uniformly too small to be divisible by D, under a stated set of edit locations and slope conditions. Any claimed uniform argument must retain the exceptional cut and the actual dependence of the three exponents on the same (q,p,i). Treating the terms as independently chosen smooth numbers discards the mechanism.

Reproduce the exact perturbation and trinomial controls with `python3 verify_one_swap_gate.py`; normal and optimized outputs must match `one-swap-receipt.json`. The verifier checks every permitted swap occurrence for q<=40, and separately repeats the q<=80 hostile probe. It does not import mathematical code from the repository or assume convergence of any positive trajectory.
