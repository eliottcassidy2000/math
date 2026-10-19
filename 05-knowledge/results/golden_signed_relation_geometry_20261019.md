# Period quotients, the signed golden dual, and a faithful Golay bridge

**Status, 2026-10-19:** PROVED elementary lemmas below, with self-checked
FINITE-EXACT controls. This is a research note, not an independently
audited canon promotion. Classical coding and finite-field inputs are credited.
No theorem about termination of every positive integer is claimed.

## 1. Inheritance and the exact new objects

We inherit, rather than re-audit:

- [THM-4581, Haar coalescence](../../01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md): the affine relation group and its almost-sure absorption, with all integer/source/clock boundaries retained.
- [THM-4590, class equidistribution](../../01-canon/theorems/THM-4590-collatz-classes-are-equidistributed-residues-windows-slowly-varying-density.md): density-one consequences do not identify every integer's component.
- [THM-4592, eleven-cell dictionary](../../01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md): for the standard, unaccelerated map C, the period-dividing-five rational points map to F_11 by the golden parity reader κ, and κ(Cx)=4κ(x).
- [THM-4594, maximal sieve and sign barrier](../../01-canon/theorems/THM-4594-the-maximal-class-decided-collatz-sieve-and-the-sign-barrier.md): actual source guards, affine thresholds and positive sign remain indispensable.
- [Golden resonance synthesis](golden_collatz_resonance_20261007.md), especially §§4,7, and [the all-depth unit-lift clock](unit_lift_clock_20261015a.md): the original group/code relationship and the general valuation mechanism are already present. The new work makes the syndrome map and signed loss explicit, and compares **tripling versus golden lifting at the same prime**.

The anchor is a lawful finite relation quotient; the niche is the dual module
lost by reflection; the wildcard is the sphere decoder as an equivariant map.
The board is: relation presentation / golden time / prime-power multiplicity /
signed character / syndrome fibre / integer endpoint height.
The canonical hostile is the out-of-phase positive pair (1,2); the corrected
near miss is identifying distinct 121-element carriers by cardinality. The
least-used sidecar here is the reciprocal five-dimensional syndrome module.

## 2. Universal finite period quotients

Write B=BS(1,3)=<a,t | tat^-1=a^3>, with a:x↦x+1 and t:x↦3x.
For every integer k≥1,

    B / <<t^k>> ≅ (Z/(3^k−1)) ⋊ C_k.                    (2.1)

The C_k generator acts by multiplication by 3. In particular, its translation
modulus is **3^k−1**, not an arbitrary prime of order k.

**Proof.** Conjugating a by t^k shows a=a^(3^k), so its order divides 3^k−1.
Every word reduces to a^b t^j with 0≤b<3^k−1 and 0≤j<k. Conversely these
normal forms form the displayed semidirect product and satisfy the defining
relations; thus no further collapse is forced. The dilation has exact order
k modulo 3^k−1, since 0<3^j−1<3^k−1 for 0<j<k. ∎

For k=5 this is C_242⋊C_5. Reducing translations modulo 11 gives the inherited
55-element G_5. Reducing modulo 121 gives a different, 605-element quotient
which retains the repeated prime in 242=2·11². These are quotients of a relation
group, not rational Collatz cycle sets.

There is a signed version. Adjoin r:x↦−x with r²=1, rt=tr and rar=a^-1.
Then

    B^± / <<t^k>> ≅ C_(3^k−1) ⋊ (C_k×C_2),             (2.2)
    B^± / <<t^k r^-1>> ≅ C_(3^k+1) ⋊ C_(2k).          (2.3)

In (2.2), t and r act by 3 and −1; at k=1 the reflection's action on C_2 is
trivial, so the abstract quotient need not act faithfully on that translation
set. In (2.3), r=t^k, and t acts faithfully by 3 modulo 3^k+1. Its order is
2k: no exponent below k reaches −1 or 1, and for k≤j<2k write 3^j≡−3^(j−k).
The normal-form proof is identical to (2.1).

This sign is an affine-group reflection. It is **not** a claim that x↦−x
commutes with the plus-one Collatz map.

To use a finite quotient for both Terras branches, halving must be defined.
Pass to a modulus coprime to 6, in particular the odd part of 3^k±1, and
retain the actual parity bit. The residue alone does not choose the branch.

## 3. A universal common quotient with a golden period ring

Let O=Z[φ], φ²−φ−1=0, and R_k=O/(φ^k−1). Fix integers r,k≥1 and put

    D_(k,r)=gcd(3^(2r)−3^r−1, 3^(rk)−1).

Then there is a canonical ring isomorphism

    O/(φ^k−1, φ−3^r) ≅ Z/D_(k,r),                     (3.1)

under which φ maps to 3^r. This is the **largest cyclic quotient with this
specified generator identification**: every such ring map factors through it.
The gcd is coprime to 3.

**Proof.** First impose φ=3^r. The quadratic relation becomes the integer
3^(2r)−3^r−1, while the period relation becomes 3^(rk)−1. The ideal generated
by these two integers is their gcd ideal. ∎

The map preserves golden rotation as multiplication by 3^r; it does not
assert equality of the full <3> and <3^r> multiplier groups. That requires
an additional order check. Examples from exact gcds:

- (k,r)=(4,1): D=5.
- (5,4): D=11, recovering φ↦4=3^4 modulo 11.
- (9,4): D=19, with φ↦5 modulo 19.
- (15,4): D=341=11·31, a genuinely composite common quotient.
- (35,2): D=71.

At eleven the time orientation is exact:

    4·3=1 mod11,  so  ×3=(×4)^−1.                     (3.2)

Thus tripling is the **inverse** of the standard Collatz step on golden
labels, not just another unspecified generator of the same subgroup.
The cycle labels are the inherited two five-cycles H and −H, with
H={1,3,4,5,9}. The rational points' denominators and real signs are carried
by κ^-1, not by the abstract affine group. Pulling the label translation
x↦x+1 back through κ is a permutation of these eleven rational points;
it is not their ordinary arithmetic translation by 1.

## 4. The match separates at the first prime-power lift

For any odd prime p≠3 let d=ord_p(3) and h=v_p(3^d−1). Multiplication by 3
on Z/p^a has, in the stratum of additive order p^b (1≤b≤a),

    orbit length = d p^max(0,b−h),
    orbit count  = (p−1)p^(b−1) / [d p^max(0,b−h)].     (4.1)

The identical formula holds on additive characters of conductor p^b.
Indeed χ_j(x)=exp(2πijx/p^a) pulls back under tripling to χ_(3j). The order
formula is the standard LTE/binomial lift: v_p(3^(d m)−1)=h+v_p(m).
This is an application of the inherited unit-lift mechanism, not a new LTE.

At p=11, d=5 and h=2 because 3^5−1=242. For every a≥1 the number of
**nonzero** tripling orbits on Z/11^a is

    2+22(a−1)=22a−20.                                  (4.2)

Their conductor-one layer has two five-cycles; every additional conductor
layer has 22 cycles, with length 5·11^max(0,b−2).

Now let φ_a be the unique root of X²−X−1 modulo 11^a reducing to 4.
At a=2 the root is 37 and

    37^5=67 mod121,  ord_121(37)=55,
    ord_121(3)=5,    37^44=3 mod121.                    (4.3)

Therefore the golden root has lift valuation h=1, not 2. For every a≥1 it
has order 5·11^(a−1), and every nonzero conductor layer has exactly two
orbits. Its total is

    2a nonzero golden-prime-branch orbits.              (4.4)

Equivalently, at the prime ideal p=(11,φ−4),

    v_11(3^(5m)−1)=2+v_11(m),
    v_p(φ^(5m)−1)=1+v_11(m)       (m≥1).              (4.5)

These valuations explain the extra factor of eleven. The tripling
period-five quotient already contains 11²; the golden period-five ring
contains only one power of that prime. A golden period of 55 is needed to
acquire the next power on that branch.

The three 121-element objects must remain distinct:

| Carrier and chosen map | Additive exponent | Cycle partition |
|---|---:|---|
| Z/121, multiplication by 3 | 121 | 1+24·5 |
| Z/121, multiplication by the golden root 37 | 121 | 1+2·5+2·55 |
| R_10=O/(11)≅F_11×F_11, multiplication by (4,8) | 11 | 1+2·5+11·10 |

The final row is the inherited 121-cell cover. It is not the first row's
prime-power thickening: even their additive groups are nonisomorphic.

## 5. The actual Golay syndrome map

Use the classical cyclic ternary Golay code with generator

    g(X)=X^5+X^4+2X^3+X²+2 over F_3.

This construction and its [11,6,5] parameters are documented by
[GAP GUAVA, §5.4-3](https://gap-packages.github.io/guava/doc/chap5_mj.html).
The companion checker independently verifies the complete 729-word weight
distribution: 1,132,132,330,110,24 at weights 0,5,6,8,9,11 respectively.
No novelty is claimed for Golay's code or its perfectness.

Let K=F_3[X]/g≅F_243 and α=X mod g, of order 11. For a coordinate vector
v=(v_i) indexed by i∈F_11, define

    S_+(v)=Σ_i v_i α^i.                                 (5.1)

Let A_{b,j} move coordinate i to 3^j i+b, where b∈F_11 and j∈Z/5. Then

    S_+(A_{b,j}v)=α^b S_+(v)^(3^j).                    (5.2)

This is an exact equivariant map from the eleven-coordinate permutation
representation of G_5 to its five-dimensional F_3 syndrome representation.
It is not a comparison of dimensions or cardinalities. The generators
R(z)=αz and F(z)=z³ obey FRF^-1=R³, so this is precisely the inherited
55-element affine relation group in a different representation.

By (3.2), **forward golden Collatz time corresponds to F^-1(z)=z^81** on
syndromes when the coordinate labels follow κ. Frobenius z↦z³ corresponds
to backward golden time. Both assertions concern coordinate permutations
and the syndrome observable, not a new Collatz map on F_243.

Let E_2 be vectors of Hamming weight at most two. Its size is

    |E_2|=1+2·11+4·C(11,2)=243.

Distinct elements cannot have the same syndrome, because their difference
would be a nonzero codeword of weight at most four. Hence S_+ restricts to
a G_5-equivariant bijection E_2→K. The group is regular on the 55 unordered
two-element supports, so the exact orbit partition is

    1 + 11 + 11 + 55 + 55 + 55 + 55.                  (5.3)

There are two weight-one orbits, one for each coefficient sign. There are
four weight-two orbits, since every support has trivial stabilizer and
four nonzero coefficient assignments. The equivariant decoder transports
this **seven-orbit partition** to K.

Thus the ball-volume identity 3^5=1+2·11² is accompanied by an actual map
and its orbit partition. The numerical identity alone would establish
neither minimum distance nor injectivity, and its orbit partition is not
1+5+5, nor one of the 121-element partitions above.

## 6. Signed reflection requires the other five-dimensional module

Let J reflect coordinate labels i↦−i, and define the reciprocal syndrome

    S_−(v)=Σ_i v_i α^(−i).

Then S_+(Jv)=S_−(v), and conversely. The irreducible polynomial of α^-1 is

    g*(X)=X^5+2X³+X²+2X+2,

which differs from g because the Frobenius exponent orbits H and −H are
disjoint. Consequently **no F_3-linear operator on the single five-dimensional
syndrome space can conjugate R to R^-1**: conjugate matrices have the same
characteristic polynomial, whereas these two have g and g*. A putative
linear extension of this G_5 representation to the signed 110-element
coordinate group is therefore impossible in dimension five.

The information-preserving repair is exact:

    v ↦ (Σ_i v_i, S_+(v), S_−(v))
       : F_3^11 ≅ F_3 ⊕ K ⊕ K.                       (6.1)

The factorization (X−1)g(X)g*(X)=X^11−1 proves the isomorphism by the
Chinese remainder theorem. In these coordinates,

    translation:  (c,z_+,z_−) ↦ (c,αz_+,α^-1 z_−),
    tripling:     (c,z_+,z_−) ↦ (c,z_+³,z_−³),
    reflection:   (c,z_+,z_−) ↦ (c,z_−,z_+).          (6.2)

The ten-dimensional nontrivial part is the minimal representation containing
the given five-dimensional module and extending its action by reflection:
reflection forces the distinct reciprocal irreducible module as well.
This follows already by the distinct simple factors of X^11−1; no assertion
about arbitrary unrelated faithful representations is intended.

There is nevertheless a faithful **nonlinear permutation** repair on the same
243-element set. Let d=(S_+|E_2)^-1 be the canonical decoder and put

    J_dec(z)=S_+(J d(z)).

Then J_dec is an involution, J_dec R J_dec=R^-1, and J_dec commutes with
Frobenius. These follow because E_2 is stable under the signed coordinate
group. Thus the full signed group acts on this set; the obstruction above
is specifically to a five-dimensional **linear** extension. In general
S_+(Jv)≠J_dec(S_+(v)) for arbitrary v: only the selected representative d(z)
is transported. A concrete nonlinearity witness is

    z=−1+α²,  w=−α³,  z+w=−α⁴−α⁵.

The three displayed expressions have canonical error weights two, one and
two. Hence

    J_dec(z)+J_dec(w)−J_dec(z+w)=g(α^-1)≠0.

More intrinsically, the decoder carry c=d(z+w)−d(z)−d(w) is a codeword of
weight at most six, and the additivity defect is S_−(c). It vanishes iff
c=0, because the intersection of the reciprocal codes is the repetition
code of minimum nonzero weight eleven. A nonlinear section restores the
signed set action, while explicitly retaining its nonzero addition carry.

The two reciprocal six-dimensional Golay codes are exchanged by J. Their
intersection is precisely the one-dimensional repetition code, and only
3 of their 729 words are shared. Global **coefficient** negation v↦−v is
a different operation: it acts as −1 on both syndrome modules and preserves
each code. It must not be confused with coordinate-label reflection.

This is the signed character sidecar: complex conjugation on characters of
F_11 exchanges H and −H; over F_3 the same exchange is between the two
irreducible five-dimensional factors. The inverse Frobenius orientation
remains necessary in either coefficient field.

## 7. What this adds toward an integer proof, and the exact missing step

None of the quotients or syndrome maps is a certificate that a specified
positive integer enters a desired Haar cylinder. A finite list of odd moduli
cannot recover that missing fact: for every prescribed residue mod 2^L and
residue mod M, M odd, CRT supplies arbitrarily large positive **and** negative
integers satisfying the same finite data. The maps also discard the actual
archimedean endpoint and, if not retained separately, the branch parity and
elapsed clocks.

There is, however, a precise usable lifting gate. If two positive integer
orbits have checked clocks a,b, let u=T^a(n), v=T^b(m). For any integer M>0,

    u≡v (mod M) and |u−v|<M  ⇒  u=v.                 (7.1)

For a coprime modulus bank replace M by their product; for an arbitrary bank
use the lcm. This is a deterministic equality certificate. If 1≤m<n it is
an induction-safe descent-through-a-join certificate. Congruence without
the strict height inequality is insufficient, and equal-time clocks cannot
be imposed universally: (1,2) is the closed phase hostile.

**OPEN research target, not a theorem.** Use the signed/golden quotient
family as cheap validators inside a source-preserving search for

    (m,a,b,M,H),  1≤m<n,  |T^a(n)−T^b(m)|≤H<M,
    T^a(n)≡T^b(m) mod M,

where the height bound is independently proved from lawful affine words,
not inferred from the modular state or sampled endpoint sizes. A universal
existence theorem for these receipts would imply positive Collatz by strong
induction. Conversely Collatz supplies such receipts, but that direction
alone does not construct a useful search bound. The new finite algebra
makes both the signed and prime-power validator losses explicit; it does
not supply the required source-dependent height/clock theorem.

The immediate next research problem is an exceptional-set mechanism that
forces a controlled-height receipt for **each** integer outside a finite
verified base. Another average syndrome count or another Haar-a.e. statement
cannot replace this obligation. The cycle-period/interval dictionary and the
resonance frequencies 5,12,41,53 remain inherited; no code action is claimed
to force those integer periods or individual termination.

## 8. Reproduction, hostile controls, and scope

Run from the repository root:

    python3 -B 04-computation/experiments/golden_signed_relation_geometry_20261019.py
    python3 -B -O 04-computation/experiments/golden_signed_relation_geometry_20261019.py

The saved JSON outputs agree. There are 16,176 explicit checks: all 243
syndromes and all 55 group elements on each radius-two error, all 729 codewords,
the complete finite Fourier rank, signed quotients k≤8, common quotients
r≤8/k≤40, and exact tripling/golden orbits on Z/11^a for 1≤a≤4. These finite
universes are controls; (2.1)–(4.5) have the general proofs above.

Hostiles: a reflected Golay word whose original syndrome becomes nonzero;
nonisomorphic additive groups at cardinality 121; different lift valuations;
1≡12 modulo 11 despite inequality; the closed (1,2) phase pair. The literal
sources retain T^7(223)=T^15(233)=425 and T^6(334)=T^6(335)=425.
The finite checks concern these new maps and quotient distinctions; they do
not repeat a broad Collatz census. Method cards used: **The same representation is not the same
carrier** and **Compose witness fibres before taking their scalar shadow**.
