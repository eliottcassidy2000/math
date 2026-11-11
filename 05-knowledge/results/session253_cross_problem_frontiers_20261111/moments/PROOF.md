# Carried Laurent collisions, Smith precision, and the Chan–Pak interface

Status: **PROVED IN THIS NOTE** relative to the scoped inherited Laurent row formulas and Session252 interface identification. New exact tests support the formulas, not their all-height quantifiers. No independent canon promotion, external priority claim, or resolution of the general Laurent first-return problem or JC(2).

## 1. Main results and exact inheritance

The main new result is a full local Smith partition at every carried Laurent collision. The mechanism also yields an actual algebra isomorphism to the Chan–Pak five-dimensional interface at a ternary collision, together with a sharp obstruction to transporting all five interface directions across the deformation.

Sources (paths relative to repository root):

* `05-knowledge/results/third_20260906_laurent.md`, §§2–3: normalized carried rows, polynomial remainder, and small-root collision at x=-r. Its characteristic-factor divisibility and endpoint-39 closure are inherited, not new here.
* `01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md`: individual factorial-row roots. Not needed for the new formal algebra; needed for actual positive-parameter noncancellation consumers.
* `01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md`: real-core signed duplication, a positive mechanism whose general carried extension remains open.
* `05-knowledge/results/continuing1_20260906_jets_projective.md`, §§2–3: complete projective jets, reciprocal denominators, and module splitting. The new result is a parameter-DVR module computation, not an unqualified p-adic specialization of this theorem.
* `05-knowledge/results/session252_fivefold_bridges_20261111/algebra/PROOF.md`, §1: the Chan–Pak positive interface algebra is the centralizer of E13, equivalently End of a dual-number module of lengths 2 and 1. Seven-state arithmetic realization and concatenation remain additional constraints.
* `01-canon/theorems/THM-4415-source-normal-row-fourteen-same-row-response-rank-obstruction.md` and `THM-4426-source-normal-row-fourteen-weight-eighteen-memory-repair.md`: JC same-row rank obstruction and restored two-memory split conic.

Inherited claims are used in their stated scopes. We do not re-audit those claims. The concept board is: inverse carries; Eisenstein collisions; full Smith modules; centralizers and liftability; positive return signs; and JC memory conics. The canonical hostile is deleting the lower carry. The corrected near miss is mistaking the common five-dimensional special-fiber algebra for five liftable directions. The least-used sidecar is the integral lattice inside the generic commuting algebra.

## 2. Balanced Smith theorem for a totally ramified collision

Let k be a characteristic-zero field, A=k[[δ]], and

  P(t)=t^r+a_(r-1)t^(r-1)+...+a_0,

where all a_j lie in δA and a_0 lies in δA minus δ²A. Set O=A[t]/P. Then O is a DVR with uniformizer t, residue field k, and δ=t^r times a unit.

For a nonzero response F in O, let ν=ord_t(F)=qr+s with 0<=s<r. Multiplication by F on the free rank-r A-module O has **exact Smith exponents**

  q repeated r-s times, and q+1 repeated s times.                 (1)

In particular, its determinant has δ-order ν, its largest inverse-denominator exponent is ceil(ν/r), and

  length_A ker(M_F mod δ^n) = sum_i min(n,α_i)
                           = min(nr,ν).                        (2)

Here length is appropriate over arbitrary k; when the residue field has finite cardinality Q, cardinality is Q raised to this length. We do not assign a finite cardinality to the characteristic-zero residue field.

**Proof.** Eisenstein irreducibility gives a domain. In O the relation writes δ times a unit as -t^r: the constant coefficient divided by δ is a unit, and the other terms contribute multiples of t. Thus its unique maximal ideal is (t), and the finite complete local ring is a one-dimensional regular local domain, hence a DVR. Alternatively successive t-adic division constructs the valuation directly. Write F=t^ν u, with u a unit. Multiplication by u is an A-module automorphism. The quotient O/F O is therefore O/(t^ν). Since δ^n O=t^(nr)O, its quotient modulo δ^n has k-dimension min(nr,ν). This determines the elementary divisors uniquely and gives (1)–(2). QED.

**Boundaries.** The Eisenstein condition is essential: P=t²-δ², F=t gives Smith (0,2), not (0,1). A determinant does not determine inverse precision: diag(1,δ,δ) and diag(1,1,δ²) both have determinant δ², but largest losses 1 and 2. These are new-claim hostile tests, not inherited audits.

## 3. Applying the theorem to the actual lower carry

Use the inherited normalized rows for h>=1:

  p_(h,x)(t)=sum_(j=0)^h (2h+1)! (x+h)_(h-j)
                    t^j / ((3h-3j)!(1+2j)!),

  q_(h,x)(t)=sum_(e=-1)^(2h) (2x+2h)_(2h-e)
                    t^e / ((6h-3e)!(2+2e)!).                 (3)

Parentheses with subscript mean falling factorials. The inherited exact lower-carry division makes R=q mod p a polynomial in Q[x,t]. Fix 2<=r<=h, put δ=x+r, and take the rank-r small-root factor P of p over Q[[δ]]. Hensel factorization is legitimate because

  p_(h,-r)(t)=t^r a(t), a(0)!=0.

All lower coefficients of P are δ-divisible and its constant term has order one, since the constant term of p does and the complementary factor has unit constant term. Hence P is Eisenstein.

The literal row coefficients give

  q_-1=δ B(δ), B(0)!=0;
  q_e is δ-divisible for 0<=e<2r;
  q_e is integral for e>=2r.                               (4)

In O, δ/t is integral and has t-order r-1. The lower carry δ B/t therefore has order r-1; every regular term in (4) has order at least r. Consequently

  ord_t(R)=r-1.                                              (5)

There is no cancellation at that first order. Applying (1), the **complete small-block partition** is

  (0,1,...,1), with r-1 ones.                                (6)

Its cokernel is (A/δ)^(r-1), its determinant has exact order r-1, and its inverse has exactly a simple δ-pole. Thus arbitrarily large collision multiplicity does not create arbitrarily deep worst precision loss on this block.

This strengthens the inherited small-root valuation/characteristic divisibility to a complete lattice statement. It does not establish positivity at the physical domain x>=1. The collision x=-r is an algebraic boundary outside that physical domain.

**Complement firewall.** When r<h, (6) concerns the Hensel small block only. Extra Smith factors may arise in the complementary block. Section 8 now proves its invertibility for h-r<=3, uniformly in r; arbitrary h-r remains open. For r=h there is no complement, so the entire h-by-h multiplication matrix has partition (0,1^(h-1)), and the multiplicity of x+h in its determinant is exactly h-1, for all h>=2. The tests find invertible complements for 2<=r<=h<=6, but this is FINITE-EXACT only.

**Carry deletion.** If the actual e=-1 term is removed, every surviving term has order at least r, so every small-block Smith exponent is at least one. In particular its special fiber is zero. This operation changes the local module, not just a harmless endpoint coefficient.

## 4. Exact Chan–Pak bridge at a ternary collision

Modulo δ, O becomes k[t]/t^r, and R becomes c t^(r-1), c!=0. For r>=2 multiplication by R is rank one and square zero. Its Jordan type is (2,1^(r-2)), so

  Cent(M_R mod δ) = End_(k[ε]/ε²)
                (k[ε]/ε² ⊕ k^(r-2)),

with dimension r²-2r+2. For r=3 this is exactly the Session252 five-dimensional interface algebra. After choosing a Jordan basis, R becomes E13 and its centralizer is

  [[d,o,j],[0,e,i],[0,0,d]].                                 (7)

The multiplication rule is literally Session252 (A1), including the cross term o I in the (1,3) coordinate. This is an algebra isomorphism based on the same nilpotent module, not a comparison of equal dimensions.

**The liftability obstruction is equally concrete.** The centralizer of M_R over A is exactly multiplication by O. To prove this, let K=k((δ)), L=Frac(O). The value group of K(R) contains r=ord_t(δ) and r-1=ord_t(R), hence contains 1. Its ramification index over K is r, forcing K(R)=L. Therefore the generic commuting algebra is L. An element α of L gives an A-linear endomorphism of O precisely when αO is contained in O; evaluating at 1 shows α is in O. Thus Cent_A(M_R)=O.

Specialization of liftable commuting endomorphisms has dimension r, whereas the special-fiber centralizer has dimension r²-2r+2. At r=3 two of the five directions are new only at the singular fiber. In the basis (t²,t,1), the liftable reduction is

  [[a,b,c],[0,a,b],[0,0,a]],                                 (8)

so it is cut out inside (7) by e=d and i=o.

An explicit model is P=t³-δ and R=t². With

  T=[[0,1,0],[0,0,1],[δ,0,0]],

we have T³=δI, M_R=T², and M_R(0)=E13. The liftable algebra consists of polynomials in T; its special fiber is (8). The coefficient e independent of d, or i independent of o, generally cannot be lifted while continuing to commute with the carried response.

**Loss ledger.** Source: actual carried Laurent small block. Target: Chan–Pak ambient interface algebra. Map: specialize, take centralizer, choose Jordan basis. Preserved: multiplication, nilpotent type, and cross term. Lost: physical positive parameter, coefficient signs, tournament realization, and deformation liftability. Sidecars: actual cut masks/seven-state arithmetic on the Chan–Pak side, and the O lattice with conditions (8) on the Laurent side. No new all-height Chan–Pak covering theorem follows.

## 5. Connection to complete projective jets and JC sectors

### Complete jets

The projective-jet theorem makes inverse denominators intrinsic by keeping complete local Hasse banks and unit coordinate changes. The small-block theorem supplies an exact module of the same kind: maximum recovery loss 1, but r-1 independent deficient coordinates. This distinction matters when the question is fiber size or a simultaneous family of conditions rather than recovery of one coefficient. Eigenvalue orders (r-1)/r, determinant order r-1, and worst Smith loss 1 are three different observables.

The map is an abstract DVR-module transfer, not a claim that an arbitrary parameter specialization δ=p^e preserves all rational units. An arithmetic consumer must clear denominators and exclude primes dividing pivot units before asserting scaled exponents (0,e,...,e).

### JC split-memory conic: a proved local response formula

THM-4426 writes its unpaid bracket after restoring the two memory coordinates as

  k u v + αu + βv + c=0, k!=0.

Put U=u+β/k, V=v+α/k, and D=(αβ-kc)/k². The equation is UV=D. Over a DVR A, at a point with U,V in A, the relative linearized response to the restored memory pair is

  A² -> A, (a,b) -> V a+U b.

Its sole nonzero Smith exponent is exactly min(v(U),v(V)); the cokernel is A/(δ^min(v(U),v(V))). This follows because its image ideal is (U,V). The smoothing datum v(D)=v(U)+v(V) alone does not determine that response loss. For example U=1,V=δ⁴ and U=V=δ² have the same D=δ⁴ but losses 0 and 2. Thus a useful boundary sidecar is the pair of branch valuations, not merely the conic discriminant.

This is an exact local formula for the inherited split-conic normal form. It does not assert new JC source entry, polynomial termination, a new row-sixteen lift, or that every such DVR arc is realized by a global Keller pair. Its practical use is to stratify source-memory arcs by branch imbalance before trying to pay further row conditions. On a unit G_m chart, min(v(U),v(V))=0 and there is no linearized local precision debt; an apparent debt there must come from a restricted source packet or another constraint.

## 6. Closest next progress, ranked by actual remaining obligation

1. **Laurent boundary-module theorem: complete here.** It gives an all-height exact small-block partition and exact multiplicity at x=-h, plus a genuine interface/liftability bridge. This is a useful small theorem, not the sharp first-return conjecture.
2. **All-height complement invertibility: a focused next theorem beyond the new four strips.** At x=-r, divide p by t^r and q by t^(2r). The surviving factorial rows have indices shifted by r. Proving their coprimality would upgrade (6) to the full partition (0^(h-r+1),1^(r-1)) at every collision. This is not supplied by individual real-rootedness; it is a reduced cross-row separation problem. Section 8 supplies an unbounded-parameter proof for complementary degrees 0 through 3. General complementary degree is still open.
3. **Physical two-rung signs: still the important hard step.** Forced characteristic factors, even with exact local Smith modules, do not prove positive residual coefficients for x>=1. A possible usable target is a deformation-compatible Gram factorization or an inequality for the full Bezout form. The inherited coefficientwise-PSD hostile and phase-cone separator must remain controls. Do not equate our centralizer with a positivity cone.
4. **JC boundary arc stratification: inexpensive and concrete.** Transport the actual next-row residual along the split-conic coordinates and keep (v(U),v(V)), distinguishing source restriction from genuine singular response. The displayed conic alone cannot prove the next compatibility; the next-row residual is the missing input.
5. **Chan–Pak operation filtering: explicit candidate.** Test which proposed state rewrites land in or preserve the liftable three-dimensional subalgebra (8), versus using the two singular-only directions. This is a diagnostic for this specific bridge; all legal Chan–Pak rewrites need not preserve (8), so it is not an additional universal restriction on that problem.

## 7. Exact reproduction

Run `python3 probe.py` from this directory. The script constructs the actual polynomial remainder by exact lower-carry division, builds the multiplication matrices for h=2,...,6, and checks the full terminal-collision Smith partition through a unit pivot plus the first-order Schur complement. It also checks ranks for every 2<=r<=h<=6, the determinant-only and non-Eisenstein hostiles, and the ternary centralizer/liftable equations. It uses SymPy and exact rationals. The all-height claims are proved above, not inferred from that finite universe.

## 8. Native upgrade: exact full Smith partitions in four all-height boundary strips

This continuation closes part of the complementary-block obligation identified above. It is **PROVED BY EXACT SYMBOLIC CERTIFICATE**, uniformly in the parameter r, not a bounded scan in h or r.

Put ell=h-r and z=2r+1. At x=-r, after dividing out the zero-root powers, the monic complementary first polynomial is

  P_(ell,r)(t)=sum_(i=0)^ell [ell!/(i!(3ell-3i)!)]
                       product_(j=2i+1)^(2ell)(z+j) t^i.

Up to a fixed nonzero factorial scalar and t^(2r), the second polynomial is

  Q_(ell,r)(t)=sum_(i=0)^(2ell) [(2ell)!/(i!(6ell-3i)!)]
                       product_(j=2i+1)^(4ell)(2z+j) t^i.

Both formulas follow directly by setting x=-r in (3): the first surviving coefficients have indices r and 2r. The inverse carry vanishes at a nonzero complementary root. Thus the complementary multiplication block is invertible exactly when Res_t(P,Q) is nonzero.

For ell=1,2,3 the exact resultants factor as

  Res(P,Q)=(-1)^ell product_(j=1)^ell (r+j)^j(2r+2j+1)^j * B_ell(r),

where B_ell has strictly positive rational coefficients. The residual degrees are 2,10,24. Their complete coefficients and normalization are printed in `probe.out` under `COMPLEMENT_CERTIFICATE`; they are part of this certificate. In particular

  B_1=(54r²+163r+120)/45.

The symbolic resultants are independently checked against literal Sylvester determinants at 8ell²+1 integer values of r. This is an identity proof: every determinant term contains 2ell coefficients of P, each of degree at most 2ell, and ell coefficients of Q, each of degree at most 4ell, giving degree at most 8ell². Agreement at that many distinct points plus one proves equality, rather than estimating a sign from sampled points. Positivity uses the complete residual coefficient list, not those evaluations.

It follows for every integer r>=2 and ell in {0,1,2,3}, with h=r+ell, that the FULL h-dimensional response multiplication matrix at the parameter DVR x+r has exact Smith exponents

  (0 repeated ell+1 times, 1 repeated r-1 times).             (9)

Consequently the multiplicity of x+r in its determinant is exactly r-1, for each of these four strips, with h unbounded. Its full inverse has a simple pole, its special-fiber rank is ell+1, and its cokernel is (A/δ)^(r-1). This is new full-matrix information, not merely the small-block statement.

The next reasonable target is ell=4 with a factored positive resultant certificate, or a structural proof that all B_ell have the indicated sign and no positive roots. The stronger coefficientwise-positive pattern is a useful conjecture, not asserted for arbitrary ell. Even an all-ell complement theorem would concern x=-r boundary parameters and would not by itself prove physical x>=1 noncancellation.

### A bounded, literal next JC proof task

A separate near-term native target is the row-sixteen continuation of THM-4438's *specific* rational boundary family, not JC(2). Inherit the exact G_m parameterization and its ten-dimensional terminal tangent fiber. Add the complete prefix-preserving valuation-sixteen Hamiltonian packet

  R_16=sum_(b=0)^8 s_b p^(16-2b)y^b.

Completeness follows from the same distinct-leading-monomial argument as THM-4415: ord_t(p^a y^b)=a+2b and its leading polynomial is x^b. The missing computation is precise: form the row-sixteen bracket matrix and its full left nullspace on that inherited family; map these nine source coefficients to that quotient; then test the projected-depth remainder after solving the paid classes. One must retain any new unpaid class and the dependence on the G_m parameter instead of extrapolating the row-fifteen rank-one result. If an unpaid class survives, restore lower-valuation source memory before declaring a geometric obstruction, as THM-4426 demonstrated. This is a realistic finite-row proof target with a known input and an exact stopping criterion. No row-sixteen closure is claimed here.

## 9. Incoming module/switch connection, preserving the deformation boundary

Concurrent commit787a2f04 supplies
`chanpak_five_core_bridges_20261110/algebra/FIVE_MODULES.md`: the complete
five indecomposable modules of this classical algebra, and the inherited
B95 central unit switch. These are incoming results, not this lane's
discoveries. The natural three-dimensional carrier is its length-three
indecomposable; the smaller four modules kill the central E13 direction.

This sharpens our liftability statement: the central-switch direction IS
one of the liftable directions, because E13=T(0)^2 and T^2 commutes with
the entire carried response in the explicit model. The two nonliftable
directions are the independent diagonal split and independent off-diagonal
split, not the central defect itself. This is an algebraic lift of the
direction. It does not lift an arbitrary pair of actual tournament
operators with that difference, preserve their signs, or select an absolute
Hamilton target. The incoming source-sensitive H141 frontiers make that
physical-source obligation concrete.
