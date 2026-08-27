---
id: THM-4280
title: "Integral three-channel fat-contact observer and sharp five-jet bound"
status: >
  PROVED RELATIVE TO THM-4241/4259/4272/4279 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On the actual O=Z[omega] lattice
  Hom(J(C_0),E_0), the formal-log channels in degrees 1,2,4 recover the
  complete Hom class; equivalently, with the target-value sidecar, restriction
  to 5Q_epsilon is faithful on global C_0->E_0 morphisms. The only other
  minimal inherited triple is 2,4,7. Actual local ramification indices are
  exactly 1,2,4, and the length-five bound is sharp. A direct corollary
  classifies ramification-four degree shells by Eisenstein norms, makes 3Q
  sharp on every degree 8s^2 shell, and gives exact hidden/visible-coset
  observers. Degree-34/42 maps are already nonconstant on 2Q_epsilon. This
  arithmetic compression is not stable under scalar extension. No new
  incidence, raw Keller descent, exact-M=12 entry, JC(2), or DC(2) is proved.
source: codex-creative-frontiers-20260827
depends_on:
  - THM-4241-w0-hidden-cm-saturation-and-visible-hidden-index-four-gluing
  - THM-4259-w0-explicit-hlambda-normalization-and-glue-dictionary
  - THM-4272-lambda-zero-a23-contact-and-e0-infinity-jet-obstruction
  - THM-4279-four-channel-formal-log-hasse-observer-for-e0-hom-at-fat-contact
related:
  - THM-4274-confined-confluent-observer-fibre-and-density-transport
  - THM-4275-opposite-parity-attachment-observers-and-confluent-sample-matroids
scripts:
  - 04-computation/jc23_integral_three_channel_fat_contact_observer_thm4280.py
  - 04-computation/jc23_integral_three_channel_fat_contact_observer_independent_audit_thm4280.py
  - 04-computation/jc23_integral_three_channel_contact_shell_corollary_thm4280.py
outputs:
  - 05-knowledge/results/jc23_integral_three_channel_fat_contact_observer_thm4280.out
  - 05-knowledge/results/jc23_integral_three_channel_fat_contact_observer_independent_audit_thm4280.out
  - 05-knowledge/results/jc23_integral_three_channel_contact_shell_corollary_thm4280.out
script_sha256:
  - 0db05f57a3aee190e71e44cee0a6e70ccf4ded5fc847ead7d5e7af4aacaaaf01
  - 86d5de3bb737a06343c030a30942e8fe269b515a5d418e6c10aeeba452292747
  - 400f8ca05a0798197716f9190e57fa5a53803cf6c2acb3b1325d1dc8fbdeebea
output_sha256:
  - 8a30a8f8880ab470f38331ffb1387252d504cc067c6264674ee128c4bb5daa0a
  - 726ed51dcd2c6acab3bd64d8b8f839904969eeb0db6eda4050292fbc9c35da26
  - 55e333c03f6d58049a3dba557d96b70c25b0472764f5e25c76fcece8b2a81bb6
hash_basis: raw LF bytes
audit: >
  PASS. Two independent downstream algebra audits agree after inheriting the
  normalized geometric channel matrix from THM-4259/4272/4279; neither script
  independently rederives that matrix. The primary SymPy and standard-library
  implementations verify every stated rank, kernel, hostile, residue, and
  length-twelve census. A third optimization-safe audit checks the direct
  norm-shell corollary's mod-two descent, examples, and bounded hostiles.
  Normal, optimized, and fixed-hash-seed outputs byte-match.
---

# THM-4280 -- integral three-channel observer and sharp `5Q` bound

**PROVED RELATIVE TO THM-4241/4259/4272/4279 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THE ARITHMETIC CLAIM IS NOT BASE-CHANGE STABLE. NO NEW
INCIDENCE IS EXCLUDED; `JC(2)` AND `DC(2)` REMAIN OPEN.**

## 1. Statement

Retain THM-4272's characteristic-zero special fibre

```text
C_0:x^6+y^4=1,                 E_0:Y^2=X^3+1,
Q=Q_epsilon,                   b=1/S at Q,             (1)
```

and put

```text
O=Z[omega],                    K=Q(omega),
omega^2+omega+1=0,             kappa=zeta_12^5,
M=Hom(J(C_0),E_0).                                      (2)
```

THM-4241 and THM-4259 give the complete `O`-basis

```text
M=O<u,f,g,h>,                  2h=v+omega^2 f+g,       (3)
```

where `u,v` are an orthogonal degree-four visible basis. For a morphism
`m:C_0->E_0`, translate its value at `Q` to the target origin and apply the
elliptic formal logarithm as in THM-4279:

```text
mbar=m-m(Q),                   L_m(b)=log_E(mbar(b)).  (4)
```

Let `c_j(m)` be the nonzero normalized scalar multiple of `[b^j]L_m`
specified in Section 2. Row normalization never changes equality or
vanishing of a labelled channel.

> **Theorem.** For all morphisms `m,n:C_0->E_0`,
>
> ```text
> m=n
> iff
> m(Q)=n(Q) and
> (c_1(m),c_2(m),c_4(m))=(c_1(n),c_2(n),c_4(n)).      (5)
> ```
>
> Equivalently, restriction to the length-five fat point is faithful:
>
> ```text
> m|_(5Q)=n|_(5Q) iff m=n.                             (6)
> ```

Among subsets of THM-4279's inherited active channels `{1,2,4,7}`, exactly
two inclusion-minimal subsets recover the actual `O`-lattice:

```text
{1,2,4},                       {2,4,7}.                (7)
```

Every two-channel subset fails. In contrast, after tensoring through the
chosen embedding to `C`, all four channels are still necessary, exactly as
THM-4279 asserts.

The actual local ramification spectrum is

```text
{e_Q(m):m nonconstant}={1,2,4}.                        (8)
```

Thus every nonconstant actual morphism is nonconstant on `5Q`. This is sharp:
the based map `v` is constant on `4Q` and nonconstant on `5Q`.

Finally, for either inherited attachment degree

```text
deg(m) in {34,42},                                     (9)
```

one has `c_1(m)!=0`. Hence every such map is unramified at `Q` and already
nonconstant on `2Q`. Here `c_1` is only a zero test on these two degree
shells; it is not a complete one-channel fibre observer.

More generally, Section 6A proves that ramification four occurs on a
degree-`d` shell exactly when `d=4N(e)` for a nonzero `e in O`. It follows in
particular that every degree-`8s^2` shell has sharp uniform contact `3Q`.
That section also gives exact one-channel observers on affine hidden cosets
and a sharp two-channel observer on affine visible cosets.

## 2. Exact normalized channel matrix

Write an actual Hom class in the basis `(3)` as

```text
m=a u+p f+q g+d h,             a,p,q,d in O.          (10)
```

THM-4259 gives the exact pullback differentials of `f,g` and the glue
identity `(3)`. THM-4272 identifies the four distinct differential orders
`0,1,3,6`, while THM-4279 integrates them to formal-log degrees `1,2,4,7`.
Comparing the four leading coefficients and dividing each row by its proved
nonzero common factor gives, in row order `1,2,4,7` and column order
`u,f,g,h`,

```text
       [ 0  1  -kappa  (omega^2-kappa)/2 ]
C  =  [ 1  0     0               0       ]
       [ 0  0     0              1/2      ]
       [ 0  1   kappa  (omega^2+kappa)/2 ].           (11)
```

The hidden leading factors are nonzero because THM-4272, using THM-4259's
explicit formulas, proves

```text
Res(q,A_0)=6,                    Res(q,B_0)=-2.        (12)
```

The visible factors are nonzero from the explicit maps `u,v`. Exact
expansion gives

```text
det(C)=kappa!=0.                                        (13)
```

Thus the four rows form the free complex matroid `U_(4,4)` on the
chosen-embedding complexification. Equation `(13)` independently recovers
the four-channel part of THM-4279.

For completeness, the remaining eight coefficient positions on the
length-twelve contact are loops. THM-4272 proves that the source action fixes
`Q` and sends `b` to `zeta_12^(-1)b`. On a differential eigenline,
`b^n db` has character `zeta_12^(-(n+1))`, so the leading orders
`0,1,3,6` force every term on that line into the same residue class modulo
twelve. After integration, the only possible formal-log exponents in
`0,...,11` are therefore `1,2,4,7`. Thus the complete length-twelve
coefficient matroid is `U_(4,4)` plus eight loops; its 4,096 subsets have rank
census

```text
rank 0:256,  rank 1:1024,  rank 2:1536,
rank 3:1024, rank 4:256.                               (14)
```

## 3. The arithmetic field separation

Define

```text
A=p+omega^2 d/2,                 B=q+d/2.              (15)
```

The first and seventh rows of `(11)` are

```text
c_1=A-kappa B,                   c_7=A+kappa B.        (16)
```

The key point is arithmetic, not complex-linear. One has

```text
kappa^2=-omega,                  kappa notin K.        (17)
```

Indeed `kappa` is a primitive twelfth root, so its minimal polynomial
`Phi_12(T)=T^4-T^2+1` has degree four, whereas every element of the quadratic
field `K=Q(omega)` has degree at most two over `Q`.

Since `A,B in K`, equations `(16)--(17)` imply

```text
c_1=0 iff A=B=0 iff c_7=0.                             (18)
```

For an actual integral vector, `B=0` gives `d=2e` with `e=-q in O`; then
`A=0` gives `p=-omega^2e`. Using `(3)`,

```text
m=a u+e(2h-omega^2 f-g)=a u+e v.                     (19)
```

Conversely every vector in `O u direct-sum O v` kills both hidden channels.
Therefore

```text
ker(c_1|M)=ker(c_7|M)=O u direct-sum O v.             (20)
```

On this visible kernel, `(11)` reads

```text
(c_2,c_4)(a u+e v)=(a,e).                             (21)
```

Equations `(20)--(21)` prove that both triples in `(7)` are injective on
`M`. This does not contradict complex dimension: the exact value
`A-kappa B` carries two `K`-coordinates on the arithmetic lattice because

```text
K(kappa)=K direct-sum K kappa.                        (22)
```

For every subset `S` of the four active labels, its arithmetic block rank
over `K` is exactly

```text
r_ar(S)=2*[S intersects {1,7}]
        +[2 in S]+[4 in S].                           (23)
```

This proves that `(7)` lists every minimal full observer. It also supplies
explicit necessity witnesses: any observer missing channel `2` kills `u`;
one missing channel `4` kills `v`; and one missing both hidden labels `1,7`
kills `f`. In particular all six two-channel subsets fail.

## 4. Why the fourth channel returns after base change

Over the chosen-embedding complexification, `kappa` is an allowed scalar.
The nonzero vector

```text
kappa f+g                                             (24)
```

has pure channel support `{7}` under `(11)`. It lies in the kernel of the
triple `{1,2,4}`. Likewise

```text
kappa f-g                                             (25)
```

has pure support `{1}` and lies in the kernel of `{2,4,7}`. Neither vector
is an actual `K`-rational, hence neither an `O`-integral, Hom class because
`kappa notin K`.

Thus the arithmetic observer becomes rank-deficient after adjoining
`kappa`. Four scalar channels are necessary on the four-dimensional complex
space, while three labelled exact-algebraic channels suffice on the actual
Hom lattice. Any use of the three-channel theorem after unrestricted base
change is a type error.

## 5. Fibres, ramification, and the sharp fat-point boundary

Apply `(20)--(21)` to the Hom class of the based difference of two morphisms.
Equality of the three channels forces that class to vanish. Equality of the
target values then kills its remaining translation, proving `(5)`.
Conversely equality of morphisms gives equality of every sidecar and channel.

The three coefficients have exponents below five, so they factor through
`5Q=Spec C[b]/(b^5)`. This proves `(6)`. Notice that the target-value sidecar
is built into a full restriction; the three coefficients alone recover only
the Hom class modulo target translation.

The formal logarithm has a nonzero linear term, so the first nonzero exponent
of `L_m` is the local ramification index. From `(20)--(21)`, every nonzero
actual class is in exactly one of three cases:

```text
c_1!=0;                           first exponent 1;
c_1=0, c_2!=0;                   first exponent 2;
c_1=c_2=0, c_4!=0;               first exponent 4.    (26)
```

The maps `f,u,v` attain the three cases, respectively, proving `(8)`. In
particular every nonconstant map is visible on `5Q`. For `v`, `(11)` gives
pure support `{4}`, so its based expansion begins with a nonzero multiple of
`b^4`. It agrees with the constant map on `4Q` but not on `5Q`, proving
sharpness.

## 6. Degree `34/42`: one channel is a zero test

THM-4241 proves that `u,v` are an orthogonal degree-four `O`-basis of the
visible lattice. Hence `(20)` gives

```text
c_1(m)=0
implies
deg(m)=4(N(a)+N(e)) in 4Z.                            (27)
```

But

```text
34=2 mod 4,                       42=2 mod 4.          (28)
```

Therefore no degree-`34` or degree-`42` class can have `c_1=0`. Its based
formal logarithm begins in degree one, proving the final assertion of Section
1. This is stronger than the general sharp `5Q` bound on these two shells,
but weaker than recovering their complete fibres from one channel.

## 6A. Per-consumer observers and exact norm-shell contact thresholds

The preceding observer has sharper forms once its consumer is specified. Put

```text
H=O f direct-sum O g,                  V=O u direct-sum O v.       (28a)
```

On every affine `H`-coset in `M`, either `c_1` alone or `c_7` alone is an
injective difference observer. Indeed, restricting the matrix `(11)` to `H`
gives

```text
c_1(pf+qg)=p-kappa q,                  c_7(pf+qg)=p+kappa q.       (28b)
```

Here `p,q in O subset K`, while `kappa notin K` by `(17)`, so either
vanishing equation forces `p=q=0`. The difference identity `(14)` of
THM-4279 then shows that, with the target-value sidecar, restriction to `2Q`
is faithful on global morphisms whose Hom classes lie in a fixed affine
`H`-coset. This is sharp for every nontrivial such consumer: `1Q` retains
only the common target value, whereas every nonzero hidden difference has
nonzero degree-one channel.

On every affine `V`-coset, the pair `(c_2,c_4)` is an injective difference
observer by `(21)`. With the target-value sidecar it factors through `5Q`.
The length-five bound is sharp for this consumer because `v` has pure support
`{4}`, as proved after `(26)`. These are distinct restricted consumers;
neither assertion promotes the indicated channels to an observer on all of
`M`.

There is also an exact degree-shell refinement. For a nonconstant based
actual Hom class,

```text
e_Q(m)=4    iff    m=e v for some e in O-{0}.                       (28c)
```

In fact, the third case of `(26)` has `c_1=c_2=0`. Equation `(20)` first
puts the class in `V`, and `(21)` then kills its `u`-coordinate. The converse
follows from the pure support of `v`. Since `u,v` are orthogonal of degree
four, as used in `(27)`, a degree-`d` shell contains a ramification-four map
if and only if

```text
d=4N(e) for some nonzero e in O.                                    (28d)
```

Consequently, for every nonempty actual degree-`d` shell:

- if `4` does not divide `d`, `(27)` forces `c_1!=0`; every map has
  `e_Q=1`, and `2Q` is the sharp uniform noncollapse contact;
- if `d=4n` and `n` is not an Eisenstein norm, `(28c)` excludes
  ramification four, so every map is nonconstant on `3Q`;
- if `n=N(e)`, the pure map `ev` shows that the inherited `5Q` bound is
  necessary as well as sufficient on that shell.

The middle bound is sharp whenever

```text
n=N(a)+N(e),                  a!=0,                  n notin N(O), (28e)
```

because `au+ev` has degree `4n`, kills `c_1` by `(20)`, and has
`c_2=a!=0` by `(21)`. Thus its first formal-log exponent is two by `(26)`.

In particular, for every integer `s>=1`, the degree-`8s^2` shell has sharp
uniform contact `3Q`. To see that `2s^2` is not an Eisenstein norm, write

```text
N(A+B omega)=A^2-AB+B^2.                                            (28f)
```

Modulo two, the three nonzero pairs `(A,B)` all have norm one. Hence an
even norm forces both `A` and `B` even; iterating shows that the `2`-adic
valuation of every nonzero Eisenstein norm is even. But
`v_2(2s^2)=1+2v_2(s)` is odd. Therefore `(28c)` bounds ramification by two
on degree `8s^2`, while

```text
m_s=s u+s v,
deg(m_s)=4(s^2+s^2)=8s^2,
(c_1,c_2,c_4)(m_s)=(0,s,s)                            (28g)
```

shows that `2Q` does not suffice.

Further concrete sharp-three shells are

```text
degree 8:   u+v;                         8=4(1+1),
degree 20:  u+2v;                       20=4(1+4),
degree 24:  (1-omega)(u+v);             24=4(3+3),
degree 32:  2u+2v;                      32=4(4+4).      (28h)
```

The quotients `2,5,6,8` are not Eisenstein norms: `2,6,8` have odd
`2`-adic valuation, while `5=2 mod 3` and
`N(A+B omega)=(A+B)^2 mod 3`. Thus `(28h)` supplies honest sharpness
witnesses, not only congruence obstructions.

The connection ledger for this addendum is

```text
source:              actual global C_0->E_0 Hom classes modulo translation,
                     in a fixed H/V coset or a fixed degree shell;
target:              selected formal-log coefficients, equivalently the
                     stated truncation 2Q, 3Q, or 5Q;
map:                 translate by m(Q), apply log_E, retain the labelled
                     b-coefficients;
preserved predicate: equality inside the named affine consumer, or
                     noncollapse on the named degree shell;
destroyed data:      target value, directions outside the named consumer,
                     and the complete Hom class for a shell zero test;
restoring sidecars:  m(Q), actual global-Hom membership, the fixed coset or
                     degree, Q, b, and omega_E;
sharp hostiles:      v on norm shells; s(u+v) on degree 8s^2;
local-map hostile:   exp_E(b^3), inherited from THM-4279.             (28i)
```

Everything remains confined to the exact characteristic-zero
`W=Lambda=0` fibre. The statements are not observers on arbitrary local
maps, are not stable under unrestricted scalar extension, and supply no
transverse `W`-jet, resolved-to-raw Keller descent, new incidence deletion,
physical exact-`M=12` entry, `JC(2)`, or `DC(2)` conclusion.

## 7. Exact audits and reproduction

The primary symbolic audit verifies `(11)--(13)`, all sixteen complex and
arithmetic subset ranks, `(19)--(23)`, the two complexification hostiles, and
the degree congruence. The independent standard-library referee implements
the two-step field extension and exact Gaussian elimination without SymPy. It
also verifies the full 4,096-subset length-twelve census and exhausts all 256
Eisenstein residue tuples modulo four as a hostile control of `(27)`.

```bash
python3 -B 04-computation/jc23_integral_three_channel_fat_contact_observer_thm4280.py
python3 -B -O 04-computation/jc23_integral_three_channel_fat_contact_observer_thm4280.py
PYTHONHASHSEED=0 python3 -B 04-computation/jc23_integral_three_channel_fat_contact_observer_thm4280.py

python3 -B 04-computation/jc23_integral_three_channel_fat_contact_observer_independent_audit_thm4280.py
python3 -B -O 04-computation/jc23_integral_three_channel_fat_contact_observer_independent_audit_thm4280.py
PYTHONHASHSEED=0 python3 -B 04-computation/jc23_integral_three_channel_fat_contact_observer_independent_audit_thm4280.py

python3 -B 04-computation/jc23_integral_three_channel_contact_shell_corollary_thm4280.py
python3 -B -O 04-computation/jc23_integral_three_channel_contact_shell_corollary_thm4280.py
PYTHONHASHSEED=0 python3 -B 04-computation/jc23_integral_three_channel_contact_shell_corollary_thm4280.py
```

All three modes byte-match the maintained output for each implementation.
The corollary audit is a finite hostile control; the mod-two descent and
all-shell statement are proved symbolically in Section 6A.

## 8. Scope and firewalls

This is an arithmetic theorem about the actual characteristic-zero
`O=Z[omega]` Hom lattice. The word "integral" does not assert an integral
model, positive-characteristic formal logarithm, or validity in residue
characteristics `2` or `7`.

The exact connection ledger is

```text
source:              global C_0->E_0 morphisms modulo translation;
target:              exact formal-log channels 1,2,4;
map:                 based formal logarithm at the fixed Q and b;
preserved predicate: equality of the actual O-valued Hom class;
destroyed data:      target value and arbitrary local-contact directions;
restoring sidecars:  m(Q), global-Hom membership, Q, b, and omega_E;
base-change hostile: kappa f+g;
local-map hostile:   exp_E(b^3), inherited from THM-4279.             (29)
```

The global-Hom sidecar is load-bearing: the three channels do not observe an
arbitrary map from `C[b]/(b^5)`, just as THM-4279's four channels do not
observe every local map from `C[b]/(b^8)`. The integral exhaustion and glue in
THM-4241 are also load-bearing; omitting them repeats MISTAKE-521.

Everything here is confined to the exact `W=0`, `Lambda=0` fibre with fixed
`Q` and `b`. These tangential `b`-coefficients are not transverse `W`-jets.
The theorem does not prove that the resolved Keller response is regular on,
or descends to, the raw nonreduced contact; MISTAKE-455 remains controlling.
It supplies no new incidence deletion, no physical exact-`M=12` entry, and no
proof of `JC(2)` or `DC(2)`. **QED.**
