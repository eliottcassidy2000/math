# Rooted divisor products force automatic balanced decoder closure through height 28

**Status: PROVED elementary arithmetic lemma and conditional actual-entry
classes; FINITE-EXACT controls;
[independent proof and actual-entry audit accepted](third_20260906_decoder_root_audit.md).** The result
forces a suitable pair from every connected primitive smaller component,
then applies the inherited dual-pair mechanism. In particular, every actual
balanced `6+7` equality entry whose primitive seven-component has maximum
at most **28** is safe at clearance at least `1/14`, with no unit or special
star assumption. This is an automatic class, not a census or a universal
LRC(14) proof. LRC(14) remains **OPEN**.

## 1. Inheritance and the precise new connection

The closest proved consumer is the independently audited
[continuing1 dual-pair theorem](continuing1_20260906_lrc_dual_pair.md).
It turns a sufficiently small selected pair gcd in `V`, relative to the
whole maximum of `U`, into actual entry closure. Its old automatic bound
used `gcd(pair)<=355^(a-2)`; in the balanced case this gave `max U<=3`,
impossible for seven distinct positive labels.

The second supplier is the independently audited leaf-cofactor theorem in
[second_20260906_entry.md](second_20260906_entry.md),
with [independent review](second_20260906_entry_audit.md).
It bounds the **minimum** coordinate of a connected primitive shape.
A small minimum does not itself make the pair used in the dual theorem
coprime. The recovered sidecar is the entire product of pair gcds anchored
at that minimum, which is bounded by a strict fractional power of the
minimum under primitivity.

Actual decoder and relation typing are inherited from **THM-3818**,
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`,
and **THM-2052**,
`01-canon/theorems/THM-2052-finite-height-forces-high-rank-bounded-relation-code.md`.
The exact radius is from the audited
[overnight12 signed-box theorem](overnight12_20260906_lrc_gcd_semigroup.md).
The component-scale caps come from the completed, audited
[joint-shadow sieve](lrc14_joint_shadow_empty_core_next_sep06.md), not a
reserved theorem namespace. Proper-component LRC remains the same **CITED**
supplier as in the dual-pair theorem; no new literature claim is made.

The canonical hostile is the primorial/complementary-product shape, which
has collective gcd one but large pair gcds. It attains the rooted product
identity below, rather than contradicting it. The corrected near miss is
turning a whole six-component into a six-subset of a larger connected
component and importing the false `t<=31,950` bound. No bound on `t` is
assumed here. The recent full-atlas and scalar-versus-profile warnings in
`01-canon/MISTAKES.md` were read. Exact rooted-gcd/product statements,
constants, and synonyms were searched in current canon and result notes;
no matching automatic-class result was found. This is no external priority
claim.

Anchor: balanced entry. Niche: divisibility products at a distinguished
coordinate. Wildcard: use the minimum-cofactor theorem as a pair supplier.
The live board has six objects: primitive minimum, prime-valuation support,
rooted gcd product, actual internal pair height, distinguished outside
coefficient, and the full protected phase interval. The map is

    leaf cofactors -> primitive minimum -> actual pair gcd
      -> prohibited mixed relation -> coherent phase-grid safety.

The product step forgets the identity of the successful partner; inspecting
at most `a-1` ordinary gcds restores it. That partner need not be adjacent
to the minimum in the decoder graph. The full actual-entry theorem supplies
its internal height at most `Q`. No large-support Bezout identity is
converted into a support-three relation.

## 2. Exact rooted gcd-product divisibility

Let `w_1,...,w_n` be positive integers, `n>=2`, and let
`h=gcd(w_1,...,w_n)`. Distinctness is unnecessary for this lemma.

**Lemma 1.** For every distinguished index `i`,

    product_(j!=i) gcd(w_i,w_j) divides h*w_i^(n-2).      (1)

Consequently, if `h=1`, some partner `j!=i` satisfies

    gcd(w_i,w_j)^(n-1) <= w_i^(n-2).                    (2)

If additionally `n>=3` and `w_i>1`, the inequality for the minimum partner
gcd is strict.

**Proof.** Fix a prime. Write `alpha` for its valuation in `w_i` and `eta`
for its minimum valuation over all labels. If `alpha=eta`, the sum of the
`n-1` pair-gcd valuations is exactly `(n-1)alpha=eta+(n-2)alpha`.
Otherwise a different label has valuation `eta`; that pair contributes
`eta`, and each of the other `n-2` pairs contributes at most `alpha`.
Thus the valuation sum is at most `eta+(n-2)alpha` in every case, proving
(1). When `h=1`, the smallest of the `n-1` factors has its `(n-1)`st power
at most their product, proving (2).

If equality held in (2) for the minimum factor, every pair gcd would be
that same factor `D`, because their product is bounded above by the right
side. But the gcd of all pair gcds rooted at `i` is the collective gcd,
which is one. Thus `D=1`, and equality would force `w_i^(n-2)=1`, contrary
to `n>=3,w_i>1`. QED.

The common factor in (1) cannot be erased: at `(2,4,6)` and distinguished
label two, the gcd product is four while `2^(3-2)=2`. The displayed formula
with `h=2` repairs the failed normalization. For primitive complementary
products `w_i=P/p_i` of distinct primes, (1) is equality at every root.
So the divisibility inequality is sharp, while passing from a product to
its smallest factor can still lose information.

## 3. From leaf minima to guaranteed selected pairs

For a connected primitive actual component `V` of size `a`, the inherited
leaf-cofactor bound gives these integer upper bounds `M_a` on `min V`:

| `a` | `M_a` |
|---:|---:|
| 2 | `177` |
| 3 | `31684` |
| 4 | `6684150` |
| 5 | `1694040507` |
| 6 | `468424857663` |

For `a>=3`, the unrounded bound is
`356^(a-1)*(a-2)^(a-2)/(a-1)^(a-1)`. Its integer content is divided out
in the direction that lowers the primitive minimum. For `a=2`, the two
actual primitive coordinates have sum at most 356 and are distinct, giving
minimum at most 177.

**Corollary 2.** Some pair containing `min V` has gcd at most

    S_2=1,
    S_3=177,
    S_4=35483,
    S_5=8350122,
    S_6=2170260903.                                    (3)

For `a=2`, primitivity gives pair gcd one directly. For `a>=3`, a unit
minimum gives gcd one. Otherwise Lemma 1 gives a strict inequality
`D^(a-1)<M_a^(a-2)`, whose exact largest possible integer is (3).
In particular the three-vertex endpoint is 177, not 178: `31684=178^2`,
and the strict boundary removes 178. The other entries are checked by

    S_a^(a-1) < M_a^(a-2) <= (S_a+1)^(a-1).

These are sufficient pair bounds, not claims that every `M_a` or `S_a` is
sharp for distinct primitive actual shapes. The three-vertex pair bound
is sharp: the primitive shape

    (177*178,177*179,178*179)=(31506,31683,31862)

has pair gcds `177,178,179`, and its actual spanning edges have primitive
sums 355 and 356. It attains the asserted minimum pair-gcd ceiling.

## 4. Automatic unit-free entry classes

Keep the full inherited actual-entry hypotheses. Namely, the thirteen
physical labels are positive, distinct, primitive, and have sum at most
`Q^2`, where `Q=91^6`; their actual decoder graph has exactly the two
components

    tV union gU,  gcd(V)=gcd(U)=gcd(t,g)=1,
    |V|=a, |U|=b, a+b=13, 2<=a<=6.

Assume equality of the full bounded support-three relation span with the
decoder span, `W_(Q,3)=V_dec`. Let `L=max U`. Neither shape needs a unit.
The hypothetical-failure scale ceilings are

    (G_11,G_10,G_9,G_8,G_7)=(2,4,9,30,90).

**Theorem 3.** Every such entry is weakly `1/14` safe whenever `L` is at
most the corresponding ceiling below:

| Split | New automatic `max U` ceiling | Inherited automatic ceiling |
|---|---:|---:|
| `2+11` | `13,520,696,477` | `13,520,696,477` |
| `3+10` | `124,998,734` | `62,323,312` |
| `4+9` | `914,513` | `257,485` |
| `5+8` | `5,397` | `1,007` |
| `6+7` | **`28`** | `3` (vacuous) |

The first row is inherited, not a new result. The new rows require no
special spanning-tree type, coprime star numerators, unit, parity, or
maximum-endpoint coprimality. The balanced row is nonvacuous because seven
distinct positive integers can have maximum at most 28.

**Proof.** If `g>G_b`, the inherited gcd supplier gives weak safety.
Otherwise choose the pair `A<B` of Corollary 2 and put `D=gcd(A,B)<=S_a`,
`p=A/D,q=B/D`, and

    R=Q(p+q)-(p-1)(q-1),
    u=min U, delta=gcd(tD,gu), c=tD/delta, x=gu/delta.

Every internal primitive pair has height at most `Q` under actual entry,
so `q<=Q` and the signed-box theorem gives `R>Qq>Q`. Each ceiling in the
table equals

    min(floor(Q/G_b), floor(aQ/[7(b+1)S_a])).           (4)

Thus `x<=G_b L<=Q<R`. If `c<=Q`, the signed box would supply a literal
mixed relation `c(gu)=r(tA)+s(tB)` with all coefficients bounded by `Q`.
Actual equality forbids it. Hence `c>Q` and

    t>Q delta/D>=Q/D>=Q/S_a>=7(b+1)L/a.

The proper-component phase supplier and coprime grid now give strict
clearance above `1/14`, exactly as in the inherited dual-pair proof.
Together with the earlier `g>G_b` branch this proves the stated weak
conclusion. No uniform bound on the whole six-component scale `t` was used.
QED.

The small `L` ceilings pay the target-containment gate directly through
`G_b L<=Q`; the auxiliary lower bound on `B` in the earlier general
scale-free dual statement is unnecessary in this restricted class.
The result is a forced existence theorem for the selected pair and hence
an automatic entry class. It does not claim the underlying dual relation
mechanism or the safety of each explicit control as a new discovery.

A remaining unsafe actual balanced entry must consequently have
`max U>=29`. This is an exact additional obligation on the primitive
seven-component, not a claim that enumerating small integers settles the
remaining unbounded component shapes or physical entry.

## 5. A genuine boundary-height control with no inherited gcd rejection

For a nonvacuity control at `L=28`, let

    q=(215,251,257,263,273), P=product(q)=995780689995,
    V={P} union {83P/q_i:i=1,...,5},
    U=(2,3,4,8,16,19,28), t=3251, g=1.

Then

    V=(302746510145,314257784295,321594541905,
       329282060835,384417661719,995780689995).

The five denominators are pairwise coprime and coprime to 83. Every star
edge has reduced sum in `{298,334,340,346,356}`, all admissible in the
strict actual atlas. No two leaves are adjacent, because their reduced
sum exceeds 356. The unique star center has all five leaf numerators
83, so the incoming sufficient coprime-numerator star subclass does not
apply to this shape. The minimum `302,746,510,145` exceeds the old balanced
minimum cutoff `60,843,134,147`. Both shapes are unitless.

The seven-shape is connected through, for example,

    2--3, 2--8, 3--19, 19--4, 4--16, 16--28,

whose reduced sums are `5,5,22,23,5,11`. The full physical sum is

    8,608,905,638,154,474 < Q^2.

The complete actual graph has exactly the claimed components. Every
mixed support is excluded in both orientations. Because `t` is coprime
to every `U` label, the minimum cleared outside coefficient for the
105 two-`V`/one-`U` supports is

    t * min gcd(A,B)/gcd(gcd(A,B),u)
      =3251*174684705 > Q.

For the other 126 supports,
`t min V>Q*(28+19)` excludes a nonzero `V` contribution by amplitude.
The internal edge rows span the two weighted kernels, so full equality
`W_(Q,3)=V_dec` follows. The source also reconstructs all actual graph
edges, checks rational rank eleven, and independently tests all 231 mixed
supports by exact signed-box membership.

Every subset of seven through twelve physical labels has gcd one, and
all 4,095 full inherited joint-shadow profiles pass. Thus this example is
not already rejected by the retained hereditary gcd tests. The source pins
and uses the complete profile JSON, including every complement word.

At the minimum `v=302746510145`, the rooted gcd product is exactly `v^4`.
The smallest partner gcd is `1,151,127,415`, at partner `314257784295`.
It satisfies the inherited dual pair inequality with `L=28`. This individual
row is therefore a positive control for a proved native consumer, not a
claim that no prior sufficient method could close it.

For contrast, every cross-divisor score using a `U` pair and a `V` label
fails: its lcm is at least `min V>3Q/28`. All corresponding arbitrary-pair
native phase comparisons fail too, since `g=1`, `delta<=D`, and

    6 delta R<=6Q(u+w)<=6Q(28+19)<56*28*min V.

These failures make the dual orientation substantive; they never imply
unsafety. The literal full-row phase

    x=1301/6502

has clearance exactly `1289/6502>1/14`. It is a half-integer lift of the
all-odd smaller shape near the seven-shape's safe phase `1/5`. The new
claim is the automatic all-shape class, not first safety of this particular
row or a refutation of another valid phase method.

## 6. Reproduction and stopping boundary

[Standalone source](../../04-computation/third_20260906_decoder.py) and
[frozen output](third_20260906_decoder.out) use exact integer and rational
arithmetic without producer imports. The rooted-product universe is every
positive tuple of lengths 2 through 6 over `{1,...,6}`, with every
distinguished root, including repeated labels and arbitrary common scale.
Complementary-prime controls attain the product identity; `(2,4,6)` is
the missing-normalization hostile. The actual-entry universe is the one
fully declared boundary-height control, all 231 mixed supports, and all
4,095 retained profiles. The sharp three-shape is a component-level
arithmetic control, not an independently claimed thirteen-speed entry.

    python3 -B 04-computation/third_20260906_decoder.py
    python3 -B -O 04-computation/third_20260906_decoder.py

The all-parameter statements follow from the displayed proofs and named
inherited suppliers. No new theorem ID, external priority, exhaustive
physical-box census, or universal LRC(14) claim is made. Beyond the new
height ceilings, the exact dual gcd product and the full translated phase
packet remain separate live coordinates.

Normal and optimized output agree byte for byte. The frozen producer executes
911,473 always-active gates, including 324,720 distinguished-root cases on
55,980 tuples. Raw LF SHA256 values:

```text
source af70250f34588f838f50d1b3dcbf2f92a1f0a497451a7d7318066f26c5a5923a
output b3deb3b98c5ecc813c995b4ac9898c64c16252933b80bb3db46e3b60a0a29216
```
