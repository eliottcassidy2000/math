# Independent residue-floor audit and the retained-correlation tail

**Status: PASS -- independently audited analytic floor and original-phase tail from 2500, with FINITE-EXACT identity certificates.** The proof gives

    e4>161875/888583>9/50

in the nonnegative-root, two-anchor, two-weak-interlacer model. Retaining its already-proved coefficient correlation strengthens the same original-root envelope to s>=2500. No model-wide finite-phase or general actual Laurent noncancellation theorem follows.

The2500 refinement was proposed during this audit, incorporated by the producer, and independently checked analytically and rationally by the root referee. It is distinguished from the independent reconstruction of the producer's original residue-floor and4100 proof. No producer source was edited by this audit lane.

## 1. Scope and inheritance

The frozen model uses five nonnegative roots with sum13 and square sum59. Write x=e3,y=e4,z=e5, and let the prescribed monic quartics C,D weakly interlace the prescribed quintic B. Both exact coefficient anchors and the actual C,D formulas remain unchanged. In particular these are two different residue measures on the same beta-root support; no unproved mixed C/D Gram matrix is used.

Read `long_frontier_sep06_anchor.md`, the incoming C-residue envelope, and `continuing3_20260906_laurent_finite_phase.md`, the audited complete carried response and original-root positive-term envelope. The former supplies the earlier moment viewpoint and C-only bounds; the new argument adds the fifth C moment, the D ordinary moment determinant, and the common support bound. It does not depend on the earlier e4>1/100 theorem.

The exact coefficient surrogate(104,50,37435088/3898125) is a hostile to truncated coefficient/moment relaxations, not to the model: its beta polynomial does not have five nonnegative real roots. Its positive full response at an original phase therefore cannot refute the retained model statement.

## 2. Weak interlacing really supplies nonnegative residue measures

For a monic degree-four interlacer A of the monic degree-five B, a beta root of multiplicity m is a root of A of multiplicity at least m-1. After canceling common factors, A/B therefore has only simple poles. At each surviving pole, the ordered interlacing signs make the residue positive. Nodes whose pole canceled can be assigned zero weight. There is no polynomial part, and the monic degree difference forces the residues to sum to one.

This proves the probability-measure representation for C/B and D/B at zero, repeated, or common roots as well as strictly interlacing simple ones. Every measure is supported in[0,M], M=max a_i. Hence all ordinary and shifted Hankel matrices are positive semidefinite: their quadratic forms are integrals of p(v)^2 and v p(v)^2. Also every moment sequence satisfies m_(j+1)<=M m_j for j>=0. These are necessary conditions, not an asserted converse to the full beta-root geometry.

The referee separately checks the weak boundary B=v^2(v-3)(v-5)^2, C=v(v-2)(v-5)^2. After cancellation C/B=2/(3v)+1/(3(v-3)), so the zero atom is retained and the measure need not have five positive weights. Literal partial fractions at the strict tuple(1,3,9,22,30)/5 give positive weights for both rows and reproduce all six moments independently.

## 3. The determinant floor is exact and noncircular

Set delta=x-75. The independently reconstructed C moments through order5 are

    (1,1,3,9+delta/3,
     27+16delta/3-4y/7,
     81+54delta-59y/7+z).

The D moments are

    (1,2,7,99/4+7delta/12,
     347/4+115delta/12-6y/7,
     583/2+199delta/2-92y/7+z).

Multiplying the generating series by the reversed B denominator recovers each prescribed numerator through order5. The C shifted order2 determinant is delta/3. Its shifted order3 determinant and the D ordinary order3 determinant are exactly

    -x delta^2/27+(15/7)delta y-(16/49)y^2+delta z/3,
    (-333+2334delta-49delta^2)/144-(18/7)y.

All three are nonnegative in the model. Thus delta>=0. If delta=0, the first order3 determinant forces y=0, and the other becomes-37/16, a contradiction. Consequently delta>0.

Cauchy on the other four beta roots gives5M^2-26M-67<=0. The quadratic increases above M=13/5 and is9/20>0 at71/10, so M<71/10. The identity

    My-5z=sum_i(M-a_i) product_(j!=i) a_j>=0

holds even if a root is zero. Therefore the weak inequality z<=71y/50 is valid before proving y>0. Substituting this in the C determinant and dividing by delta>0 gives

    x delta/27 <=(2747/1050)y.

The left side is positive, so y>0 follows. Because x=75+delta>75, the same inequality yields the strict slope

    y>(8750/8241)delta.                                 (1)

The D determinant gives2334delta>=333+49delta^2+(2592/7)y. Substituting delta<(8241/8750)y and checking the positive denominator proves

    y>333/[2334*(8241/8750)-2592/7]
      =161875/888583>9/50.                              (2)

No strict z bound, division by y, or earlier e4 floor was assumed in reaching y>0. Equality/zero/repeated-root limits have been explicitly discharged. Every displayed rational constant agrees with independent exact arithmetic.

## 4. Complete carried response and the2500 consequence

The referee independently multiplies the original O/E rows, constructs beta^2+2t C_raw D_raw with both Laurent shifts, and takes the coefficientwise product. It retains the exact coefficient28 at exponent-1. Only then does it substitute

    z=(12/7)y u-xu^2+10u^3-u^4/11, u=1/s,

which is the exact equality P(-s)=0. Rebuilding every term gives the same polynomial R=Q(-s)/s^7 and all eight positive monomials as the producer's certificate.

With the producer's positive polynomials A,B,C1,D1, the entire positive-term bound is

    R/y^2 <= -26075790/7
              +A(u) x/y+B(u) x/y^2+C1(u)/y+D1(u)/y^2.  (3)

The initial independent-cap proof uses x<123 and(2), giving the certified4100 cutoff. Its exact envelope value is

    -5961609014706655321683472522194343
      /99603957308055354000000000000<-59000.

The stronger proof retains the correlation in(1). Since x=75+delta and1/y<50/9,

    x/y <75*(50/9)+8241/8750 =10962223/26250,
    x/y^2 <(10962223/26250)*(50/9)=10962223/4725.        (4)

Insert(4),1/y<50/9 and1/y^2<(50/9)^2 into(3). Every nonconstant coefficient in the resulting envelope is positive, so it increases with u>0. Its exact value at u=1/2500 is

    -653528391359305169367452997041401
      /13323669433593750000000000000<-49000.

For all original positive phases s>=2500 this proves Q(-s)<0, indeed the stronger normalized bound below-49000. The2500 proof no longer needs the independent upper cap x<123; the retained slope replaces that lost correlation. The4100 calculation remains a useful separately reconstructed intermediate control.

The original-root condition is essential: the strict rational model gives Q(-2500)>0 and P(-2500)!=0. The same omission hostile is independently checked at4100. No derivative or discriminant division is used, so multiple first zeros are included. Combined with the previously proved smallest-phase interval, the model's remaining possible original phases lie in(1/80,2500), with at most three per shape counted with multiplicity. That continuous interval remains open.

## 5. The support hostile and independent artifact audit

For the coefficient surrogate(104,50,37435088/3898125), the referee checks every principal minor of the C,D ordinary3 and shifted2 matrices, rather than only their determinants. They are positive semidefinite. Nevertheless

    mu4-(71/10)mu3=2159/105>0,
    nu4-(71/10)nu3=1091/42>0.

Thus both sequences violate the necessary bounded-support inequality. Their shifted3 determinants also reject the surrogate. This demonstrates the exact missing coordinate: lower Hankel positivity alone does not restore the anchor-compatible support interval. The new fourth-coefficient floor itself is merely a necessary consequence and is not claimed sufficient for model membership or full response negativity.

The independent source uses SymPy exact polynomial/matrix arithmetic plus the standard library, importing no producer. This is separate from the producer's standard-library sparse polynomial and permutation-determinant engine. It checks the formal moments by numerator/denominator multiplication, determinant identities, literal positive residue measures, all hostile principal minors and support violations, the complete original carry and eliminated polynomial, both envelopes, and every saved JSON polynomial as data.

The revised 2500 producer proof and source have been read in full and accepted. All saved polynomial identities, both tail envelopes, and the coupled ratios agree with the independent reconstruction. No repair is requested. The original4100 primary freeze remains identifiable by source039c9b8cc68b2b61ac23d45b1b54355fe76892f61f226a44ecedde46eaa6eb37, output21f7d93c5d1775834ead322b786a6ee0a9c5779902d9f3ca7155f84ab9ea6a67, and report18339d3dbb3377bb403f002d1c9980e7e334626c16bb923009adc23fd4f4b64d. This is a strengthening, not a correction of a false earlier claim.


## 6. Final freeze, reproduction, and promotion basis

The [independent source](../../04-computation/continuing3_20260906_residue_floor_audit.py) passes **122 always-active exact gates** in both normal and optimized Python. Captured output bytes are identical and contain only LF newlines; neither source output was post-normalized. The producer's revised 94-gate certificate is pinned as data, and its module is never imported. An explicit certificate path keeps the audit portable when the files are filed in separate computation and result directories.

```text
python -B 04-computation/continuing3_20260906_residue_floor_audit.py --certificate 04-computation/continuing3_20260906_residue_floor_certificate.json
python -B -O 04-computation/continuing3_20260906_residue_floor_audit.py --certificate 04-computation/continuing3_20260906_residue_floor_certificate.json
```

Final revised producer pins:

```text
report      5454a60bd0071cd98c55e75b97256e3d0cd0bf29ed6817ba184073e2fc0c1189
source      8974668b0c36e129910142b03a66f95626b18861a009d037325f73f22f00f205
output      fa8b9936221b9f6cbea053fdc91c84f0eafe4faa25f19ffa8af5afe2d63033f0
certificate c2a20f53670198790529fa980b7d33913b9e5a1f64c3b206f9fc34aa42e9104e
```

Final independent audit pins:

```text
source f1c005bf70e2b67c24cdeaa3b012ecf5d9c94ad7d682c124a3ee49d01be06a61
output 95253fedc6f0df0e2e482e7ea739a77a57d59bd0a359cceb8ddf3bf6466c1a29
```

The promotion basis is the complete analytic proof, including the weak-interlacing residue lemma and strict boundary discharges, together with exact identity reconstruction by a separate computational method and independent root-referee checks of the coupled improvement. It supports a **PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED** scoped result, with the original-root restriction and unresolved interval intact. It does not promote truncated moment conditions to sufficient conditions, nor establish general Laurent noncancellation. The root integrator owns the filed status change and its candidate-hash lineage.

This audit made no repository, maintained navigation, theorem namespace, or Git mutation.
