# A negative transitive-interface quartic at all large block sizes

## 1. Exact family and result

Fix the five-name quotient Q with outgoing bit rows

    Q=(4,13,8,17,7), on names A,X,Y,F,L.

Thus X->A->Y and X->Y; X,Y both beat F and lose to L; F->L; both F,L beat A. Inflate A,X,Y by directed triangles with their physical names respectively 0,1,2; 3,4,5; 6,7,8, oriented cyclically in that order. Inflate F and L by transitive tournaments TT_a and TT_b, with a,b>=1. Their names follow consecutively after8. Call the result T_(a,b), and let U_(a,b) converse Q while leaving all five child interiors fixed.

For the exact named response define

    Delta(a,b)=[a_0 a_3 c_1 c_6](F_T-F_U).

The coefficient is positive at Delta(1,1)=146 and Delta(2,2)=3716, but

    Delta(6,6)=-573454264.

In fact Delta(a,b)<0 for every a,b>=6. This is an all-order separator for this family and a counterexample to extending the positive cyclic-interface theorem unchanged to the transitive interface. It is not an exact-response collision.

The same negative quartic, together with one fixed exposed-anchor battery, separates T_(a,b) from every original-module XOR-family containing the root, for every a,b>=6. At most22 actual physical queries suffice against all such targets from this source.

## 2. The transitive interface kernel

Use the arbitrary-exterior notation in ordered_cyclic_module/PROOF.md, but orient the attachment as X->A->Y rather than A->X->Y->A. Let x,y count X,Y runs, h=x+y+xy and D=1-sigma*h. The path matrix N and its sums C_x,C_y,R_x,R_y,H are unchanged. The altered quantities are

    alpha_R=alpha+a(C_x+kappa H),
    beta_R=beta+b(R_y+mu H),
    tau_R=tau+N_yx+mu C_x+kappa R_y+mu kappa H,
    rho_x=N_yx+mu C_x, rho_y=N_yy+mu C_y,
    ell_x=N_xx+kappa R_x, ell_y=N_yx+kappa R_y.

The use of a,b in this paragraph denotes the interface polynomials of the cited proof, not the two block sizes. The complete four-term triangle kernel, with these substitutions, becomes

    K/(xy)=(2 alpha beta-alpha b-a beta)/D
      +[alpha b V+a beta W+ab((2 tau-kappa-mu)h-xy)]/D^2
      +2abVW/D^3,                                      (1)

where V=mu*h+y and W=kappa*h+x. This exact rational identity retains the full disjoint quartic; the negative terms cannot be dropped. For empty exterior it is K=xy(x+y+xy), already different from the cyclic-interface expression.

Here the actual exterior is the transitive join F->L. Put

    P=product_(v in F)(1+z_v)-1,
    R=product_(v in L)(1+z_v)-1.

The squarefree path decomposition gives exactly

    alpha=(1+P)(1+R), beta=1,
    a=1+R+PR, b=1+P+PR,
    sigma=PR, mu=tau=0, kappa=P+PR.                    (2)

For example every path starting in F and ending in L selects a nonempty subset of each chain, giving PR. Paths ending in L give R+PR. Every nonempty exterior path ends in A's in-neighborhood, giving P+R+PR. These facts account for every interface in (2).

Each physical marked triangle profile is b(1)=1,b(2)=2,b(3)=2 and zero otherwise. Apply x^s y^t -> b(s)b(t) to (1). This gives the following polynomial C(P,R), with all omitted coefficients zero:

    -120 P^6 R^6 -240 P^6 R^5 -120 P^6 R^4
    -120 P^5 R^6 -336 P^5 R^5 -288 P^5 R^4 -96 P^5 R^3
    +396 P^4 R^4 +480 P^4 R^3 +84 P^4 R^2
    +396 P^3 R^4 +1136 P^3 R^3 +720 P^3 R^2 +80 P^3 R
    +384 P^2 R^3 +701 P^2 R^2 +270 P^2 R +13 P^2
    +113 PR^2 +146 PR +28 P +8 R +8.                (3)

This coefficient extraction is finite: after the leading xy, powers (sigma*h)^r with r>4 cannot contribute to x^s y^t with s,t<=3. The checker independently expands the exact four-term kernel numerator, compares every coefficient of (1) that can contribute to this physical marked extraction, then performs the complete truncated extraction. It does not infer (3) from a fit to numerical data.

Let p_n(k)=k! S(n,k), with p_n(0)=0 for n>0. Since every nonempty subset of a transitive tournament has a unique directed path,

    [F]P^s=p_a(s),   [L]R^t=p_b(t).

Consequently, writing c_st for the23 coefficients in (3),

    Delta(a,b)=sum_(s,t) c_st p_a(s)p_b(t).           (4)

The exact physical marked-profile restoration from cyclic_child_disjoint/PROOF.md identifies the local triangle-flip coefficient with the displayed quotient-converse coefficient, since the other two marked children are balanced and the remaining children are unmarked. Thus (4) concerns the fixed physical names, not a substituted quotient response.

## 3. Strict negativity for every a,b>=6

The monotone Stirling ratios proved in covered_quotient/PROOF.md imply that for n>=6,

    p_n(s)/p_n(4) <= p_6(s)/p_6(4)  if 0<=s<=4,
    p_n(s)/p_n(4) >= p_6(s)/p_6(4)  if 4<=s<=6.

The zero coordinate s=0 is harmless. All positive terms of (3) have both exponents at most4. All its negative terms except -96 P^5 R^3 have both exponents at least4. Drop that one negative term when forming an upper bound. At n=6,

    (p_6(1),...,p_6(6))=(1,62,540,1560,1800,720).

The total positive contribution at6,6 is2080753736; the retained negative contribution is-2560896000. Therefore for all a,b>=6,

    Delta(a,b) <= -480142264 * p_a(4)p_b(4)/1560^2 <0.   (5)

At6,6 the omitted negative term is-93312000, yielding the exact value-573454264. This is a uniform all-order sign proof. It uses monotone ratios only for transitive children, not for arbitrary tournament profiles.

## 4. Original-module root families

Q is the concrete W5 quotient under the map (A,X,Y,F,L)->(3,2,4,0,1), where0 is the W hub. Thus Q is prime. Its hub child here is TT_a, not the marked C3 child A, so the cyclic witness theorem at a C3 hub does not supply this comparison.

Consider any family of original modules of T_(a,b) that includes the whole root. Every proper module lies within one child. Scalar Hamilton monotonicity first separates any strict increase. In the scalar-equal case, the hereditary equality theorem restores the complete q profiles of F,L. The marked X,Y triangles retain their balanced physical profiles under their only nontrivial original-module operations.

If A is selected, its orientation is restored after the root converse, and the exact marked-profile argument gives the unchanged negative coefficient (4). If A is unselected, uniform queries on every other whole child restore its target to the whole-converse response, and the exposed-anchor theorem separates it.

The sixteen quartic queries use a supported on{0,3} and c on{1,6}. The anchor splits A into singleton names and uses the six other compressed blocks {1},{2},X,Y,F,L. Physical0 beats exactly {1} and Y. There are six main anchor queries

    (a=B union {0} if 0->B, otherwise B; c={0}),

plus (a={0},c=empty). The parity correction queries are absent because both the compressed exterior count6 and marked outdegree2 are even. The last query already occurs in the quartic battery, giving16+7-1=22 distinct actual queries. The scalar query(0,0) is also already in the quartic list. This battery depends on the source's physical blocks, not on the unknown selected family.

## 5. Boundaries

The strict negative sign in (5) is proved for both transitive exterior orders at least6. Smaller cases are retained as controls, not asserted to follow the same sign. This is a concrete all-order root-family result, not arbitrary-prime quotient closure, full orbit injectivity for every changed interior, or unrestricted named reconstruction.

The exact local triangle-flip identity permits the displayed arbitrary exterior. The quotient-converse and root-family conclusions are whole-core statements, with justified ordered-SCC padding only. No arbitrary strong outer-context transfer of the negative quartic is claimed. All sixteen finite-difference query values are retained, so the comparison does not suppress cubic or quartic interactions.
