# Integral torsion and exact admissibility counts in two-fibre Seifert gauge

Status: PROVED as the standalone integer-lattice theorem below; exact finite
controls recorded separately. No classification of analytic torus fillings or
six-manifolds is asserted. The use of the word gauge means exactly the two
specified presentation moves, with the exceptional points labelled.

## Inheritance and the new question

The repository's August 24 Hopf/S6 intake proves that, for coprime orders
`a,b`, the signed clutch number

    p = ab ell0 - b ell1 - a ell2

is complete for the displayed Seifert presentation modulo the two standard
moves `(1,a,0)` and `(1,0,b)`. MISTAKE-480 correctly retracts the inference
that a raw-coordinate kernel automatically represents lost analytic moduli.
The new question is what changes when coprimality is removed. The answer is
an actual finite torsion coordinate in the specified presentation quotient,
not an invented analytic parameter.

The paper's source is https://alpo.ge/s6.pdf, especially its opening Seifert
comparison and Section 7.5. The download on 2026-10-19 is byte-identical to the
repo's 2026-08-24 artifact, SHA-256
`283bba102dd1d5dc346af81b28145bdaaea6654398d5032e76e97bafb9a858f2`.
The paper's global complex-S6 claim remains under audit.

## Theorem A: the missing coordinate is exactly cyclic torsion

Let `a,b >= 2`, `g=gcd(a,b)`, and let

    G = Z(1,a,0) + Z(1,0,b) subset Z^3,
    Q = Z^3/G,
    phi(ell0,ell1,ell2) = ab ell0-b ell1-a ell2.

Then:

1. `Q` is noncanonically `Z direct_sum Z/g`.
2. `phi/g` induces a surjection `Q -> Z`; its kernel is `Z/g`.
3. The torsion kernel has a distinguished generator represented by
   `tau=(0,a/g,-b/g)`. Its order is exactly `g`.
4. Every triple has a unique representative with
   `0 <= ell1 < a` and `0 <= ell2 < b`.
5. A signed integer `p` occurs as `phi(ell)` iff `g` divides `p`; when it
   occurs, exactly `g` gauge classes have that signed value, before imposing
   admissibility.

### Proof

The matrix of the two gauge generators is

    K = [[1,1],[a,0],[0,b]].

Its first Smith divisor is one, and the gcd of its two-by-two minors is
`gcd(a,b,ab)=g`. Hence `coker K = Z direct_sum Z/g`.
The coefficients of `phi` have gcd `g`, so `phi/g` is primitive. It kills
both columns of `K` and gives a surjection on the quotient. Its kernel has
rank zero and is the torsion subgroup just computed.

More explicitly, `g tau=(0,a,-b)=(1,a,0)-(1,0,b)`. If `r tau` is in `G`,
comparison of the first coordinates in
`r tau=x(1,a,0)+y(1,0,b)` gives `y=-x`; comparison of the second gives
`r/g=x`. Thus `g` divides `r`, proving exact order `g`.

Subtracting integral multiples of the two gauge vectors puts `ell1,ell2`
in the displayed residue intervals. Two representatives in those intervals
can differ by a gauge vector only if both multiples are zero. This proves
existence and uniqueness. Finally every fibre of the surjective map
`Q -> Z` is a torsor for its order-`g` kernel. This proves all five claims.

## Theorem B: a primewise count of admissible classes

Call a class admissible if

    gcd(ell1,a)=gcd(ell2,b)=1.

This is well-defined under the gauge moves. Let `C_(a,b)(p)` be the number of
admissible gauge classes with signed clutch value `p`. For each prime `q`
dividing `ab`, write

    A=v_q(a), B=v_q(b), P=v_q(p), with v_q(0)=infinity.

Then `C_(a,b)(p)` is the product of the following local factors:

- If `A != B`, the factor is zero unless `P=min(A,B)`. If equality holds,
  the factor is `1` when `min(A,B)=0`, and otherwise it is
  `q^(min(A,B)-1)(q-1)`.
- If `A=B=e>0`, the factor is zero for `P<e`, is
  `q^(e-1)(q-2)` for `P=e`, and is `q^(e-1)(q-1)` for `P>e`.

All these are statements for an individual signed integer `p`, not just
asymptotic densities.

### Proof

By Theorem A, an admissible class corresponds to a pair of unit residues
`x mod a`, `y mod b` satisfying

    b x+a y = -p mod ab,                                  (1)

and then `ell0=(p+b x+a y)/(ab)` is uniquely determined. Chinese remainders
factor this count primewise.

Fix `q`. If `A>B`, write `a=q^A a0`, `b=q^B b0`, with `a0,b0` units.
After division by `q^B`, the congruence is

    b0 x+q^(A-B)a0 y = -p/q^B mod q^A.                    (2)

There are no solutions unless `q^B` divides `p`. For each unit
`y mod q^B`, equation (2) determines a unique `x mod q^A`. That `x` is a
unit iff `v_q(p)=B`, because its reduction modulo `q` is
`-b0^(-1)p/q^B`. Therefore the count is `phi(q^B)` when `P=B`, and zero
otherwise; interpret `phi(1)=1` when `B=0`. The case `B>A` is symmetric.

If `A=B=e`, division by `q^e` gives

    b0 x+a0 y = -p/q^e mod q^e.                           (3)

Again there is no solution if `P<e`. For any `x mod q^e`, equation (3)
determines `y`. If `P>e`, the forbidden residue `x=0 mod q` is exactly the
same as the residue forcing `y=0 mod q`, leaving `q-1` permitted residues.
If `P=e`, those are two different forbidden residues, leaving `q-2`.
Each permitted residue has `q^(e-1)` lifts. This gives the stated factors.

## Corollaries and sharp boundaries

1. If `gcd(a,b)=1`, then `C_(a,b)(p)=1` exactly when `gcd(p,ab)=1`, and it is
   zero otherwise. This recovers the earlier coprime result.
2. An admissible class with `|p|=1` exists iff `gcd(a,b)=1`. Thus coprime
   exceptional orders are forced by the cyclic homology-sphere clutch
   condition in this grammar. This does not select `(3,4)` uniquely.
3. `C_(a,b)(0)=0` unless `a=b`. For `a=b=m`, it equals `phi(m)`.
   Indeed any unequal prime valuations force a zero factor; if all agree,
   the local factors multiply to `phi(m)`.
4. If `v_2(a)=v_2(b)=e>0`, then admissibility forces `v_2(p)>e`.
   Merely requiring `g|p` is not enough: the local factor `2-2` vanishes.
5. The count is periodic in `p` modulo `ab`, and over one period its sum is
   `phi(a)phi(b)`, because each unit residue pair contributes once.

## The clutch group loses the finite gauge coordinate

For an admissible triple the two-by-two relation matrix inherited from the
paper's presentation is

    R = [[a,-ell1],[b,ell2-b ell0]].

Its determinant is `-p`. Since `gcd(a,ell1)=1`, the gcd of all entries is
one. Therefore its cokernel is `Z/|p|` when `p != 0`, and is `Z` when `p=0`.
This proof requires no coprimality assumption on `a,b`.

For `a=b=5` and `p=10`, the exact count is three, with normalized triples

    (1,1,2), (1,2,1), (2,4,4).

They represent distinct classes under the specified gauge moves, all with
clutch cokernel `Z/10`. Even after swapping the two equally labelled orders,
the third normalized residue pair differs from the first two. No claim that
these are distinct analytic fillings, distinct unmarked Seifert spaces, or
distinct six-manifolds follows without the corresponding equivalence proof.

## Connection contract

Source: labelled integral Seifert twist triples.
Target: their cyclic clutch presentation, and then its abelian group.
Map: quotient first by `G`, then evaluate signed `p`, then forget its sign
when recording the finite cyclic group.
Preserved predicate: the exact order and cyclicity of the displayed clutch
cokernel; the condition `|p|=1`.
Lost information: order-`g` gauge torsion for noncoprime orders, followed by
orientation/sign information. Admissibility restricts the torsion fibre by
Theorem B.
Sidecar: normalized labelled residues `(ell1 mod a,ell2 mod b)`, or a
chosen origin plus the torsion class in the fibre of `phi/g`.
Decisive hostile: `(a,b,p)=(5,5,10)` above. Positive control: `(3,4)` has
`g=1`, so the old corrected completeness theorem is recovered exactly.

## Reproducibility and non-consequences

`verify_geometry_extensions.py` exhausts all 277,729 signed-p residue
classes modulo `ab` for ordered `2<=a,b<=32`, checks all 961 gauge Smith
matrices, and compares direct residue-pair counting with Theorem B.
The proof is unrestricted in `a,b,p`; the computation is a finite hostile
control, not its logical basis.

This advances the arithmetic classification of the filling presentation. It
does not establish the period construction, proper toric quotient, analytic
gauge equivalence, global gluing, homology, or a complex structure on S6.
