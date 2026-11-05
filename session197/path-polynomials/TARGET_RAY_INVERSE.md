# A complete inverse test for a prescribed ternary ray

Session197. **PROVED finite inverse algorithm for a fixed physical module bank and a specified pair(A,C).** This strengthens the cancellation-only lattice criterion. It still does not solve the unrestricted Hamilton spectrum.

## 1. Fixed language and target

Use the physical quotient and nonempty bank from `RATIO_LATTICE.md`. Fix nonempty leaf graphs P,R and bank graphs G_1,...,G_k. The variable modules U,V are nonempty ordered joins of bank copies. The root module is Tr_m for every m>=1. Given integers A>0 and C, the target is the exact whole family

    H_m=A3^m+C.

Write h_i=H(G_i), b_i=B_(G_i)(2), and let

    A0=2(b_P-h_P)(b_R-h_R)>0,
    C0=(b_P-3h_P)(b_R-3h_R).

The cut bound shows that C0 is never zero: for a nonempty X, B_X(2)-3H(X) is-1 at K1 and positive at every larger tournament. In particular a zero prescribed offset cannot occur in this language.

If U has u_i copies and V has v_i copies, define

    w_i=u_i+v_i, z_i=u_i-v_i.

Multiplicativity and the audited coefficient formula give

    A=A0 product_i b_i^(w_i),
    C=C0 product_i h_i^(w_i).                              (1)

The binary-mode condition is exactly the sign/lattice condition from `RATIO_LATTICE.md`: it is automatic if both leaf e-values vanish, impossible in its stated zero/same-sign cases, and otherwise Mz=t for the prime-valuation matrix of the bank ratios and the positive target kappa.

## 2. Necessary and sufficient integer constraints

A physical solution exists if and only if all of the following hold:

1. The leaf sign case permits cancellation.
2. X=A/A0 and Y=C/C0 are positive integers.
3. There is w in Z_(>=0)^k with product b_i^(w_i)=X and product h_i^(w_i)=Y.
4. There is z in Z^k satisfying |z_i|<=w_i and z_i=w_i modulo2, and satisfying Mz=t in the opposite-sign leaf case.
5. The recovered counts u_i=(w_i+z_i)/2 and v_i=(w_i-z_i)/2 each have positive total.

These constraints are necessary by (1), cancellation and literal copy counting. Conversely they produce nonnegative integer u,v with both modules nonempty. Build their ordered joins and substitute them into the quotient. The exact physical coefficient theorem gives precisely(A,0,C), proving the prescribed ray at every positive m.

Unlike the cancellation-only existence test, one must **not** add common padding after this step: padding changes A and C. Nonemptiness is enforced by item5, not repaired silently.

## 3. A finite algorithm with an explicit complete bound

Every b_i>=2. Thus every possible solution obeys

    sum_i w_i <= floor(log_2 X).                            (2)

This bounds the total number of bank copies independently of graph order or a guessed search cutoff. If X is not a positive integer, reject before using the bound.

Enumerate the finite factorizations in item3. A direct recursive procedure repeatedly divides the remaining X by b_i and the remaining Y by h_i, stopping when either divisibility fails; b_i>=2 guarantees termination even when h_i=1. For each w, enumerate z_i in{-w_i,-w_i+2,...,w_i}, check the relevant lattice equations and nonempty totals, and return the corresponding graph recipe if one passes. If none passes, the target ray is absent from this language at every order.

Smith normal form can reject an impossible unrestricted lattice equation before this bounded search. Its particular unrestricted solution is not by itself enough: it may violate the required bounds, parity or nonempty allocations. The checker searches every bounded allocation, not merely the first Smith solution.

This is a terminating exact inverse algorithm; no polynomial-time claim is made. A witness graph for any m has explicitly bounded order

    m+2+|P|+|R|+sum_i w_i |G_i|.

The bound is for the prescribed ray's source construction. It does not delete padding inside a live graph or transport a simplified state through a later wrapper.

## 4. Complete offset-65 classification in the worked bank

Take P=K1,R=Tr4 and bank(K1,W), where W=[8,1,3,6] has(h,b)=(5,32). Here A0=30,C0=-13 and the cancellation equation has the unique difference vector

    z_K=4, z_W=-1.

Fix C=-65. Equation(1) forces exactly one W copy in total, w_W=1. The difference vector therefore places it in V. The singleton copy count must satisfy w_K>=4 and w_K=4 modulo2. Write w_K=4+2t, t>=0. Then

    U=Tr_(4+t),
    V=Tr_t -> W  (with no prefix when t=0),
    A=30*32*2^(4+2t)=15360*4^t.

Consequently the complete prescribed-offset family in this bank is

    H_m=15360*4^t*3^m-65,  t>=0,m>=1,                      (3)

with order m+15+2t. In particular15360*3^m-65 is realized, while30720*3^m-65 is impossible: its required singleton total is five, which cannot have difference four. The latter failure is not visible from unrestricted rational-span or lattice membership alone.

The first actual source of (3) has H=46015; at t=1 its eighteen-vertex source has H=184255. Both graphs are reconstructed and independently counted in the verifier. Every member of this particular family remains divisible by5, so these examples are arithmetic inverse certificates, not prime-spectrum progress.

### Nonempty-module boundary control

For P=R=Tr2 and bank(K1), cancellation is automatic, A0=18,C0=1. Target(A,C)=(72,1) uses two singleton copies, one in each variable module, and gives the actual217 source. Target(36,1) would provide only one singleton copy in total, so it is impossible with both U,V nonempty. Allowing an empty input would silently change the theorem's domain and is not permitted.

## 5. Reproduction and dependencies

Run `python3 verify_target_ray_inverse.py` and `python3 -O verify_target_ray_inverse.py`. The checker implements the finite factorization/allocation algorithm, checks all successive dyadic amplitudes in the worked range against the all-order parity classification, replays three physical source graphs, and verifies the nonempty boundary and invalid arithmetic cases. Explicit exceptions remain active under optimization.

The substantive dependencies are `RATIO_LATTICE.md` (SHA256 `a6f9c0b09d715865d01bbd62d95f2261c59e1f4a7eea2ed05cc09ac90aff60a3`) for multiplicativity and exact cancellation, and `UNIT_OFFSET_RIGIDITY.md` (SHA256 `3c116086bfd697eaeaf5fcbf46308f1e31cdf653eb2a651053f01bce46ce374b`) for the physical coefficient formula and nonzero boundary statistic. The new point here is that fixing the amplitude bounds every copy count, turning the inverse problem for the complete ray(A,C) into a finite exact integer feasibility test.
