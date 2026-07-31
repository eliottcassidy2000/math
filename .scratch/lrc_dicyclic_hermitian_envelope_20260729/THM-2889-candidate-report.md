# THM-2889 candidate -- two dicyclic shadows and their joint carrier

**Status:** `PROVED ABSTRACT ALGEBRA + VERIFIED-EXACT COEFFICIENT
SPECIALIZATION`, scratch candidate for the reserved THM-2889 namespace.
No physical \(Q_8\) action, row exclusion, or LRC(14) conclusion is claimed.

Suggested dependencies:

- THM-2868, for the split Prony pair and projective \(C_{13}\) orbit;
- THM-2882, for the event-twisted honest \(C_{169}\) coefficient lift;
- THM-2884, for the physical semantic \(V_4\) horn;
- THM-2886, for the transported origin-odd marked subcopy and the
  \(\operatorname{diag}(E,1)\) recombination boundary; and
- THM-2887, for the unique Arf-one \(Q_8\) lift and its physical stopping
  line.

## 1. Inheritance and the two-character obstruction

Write

```text
QA=(1,0), QB=(0,1), QAB=(1,1)
```

in the semantic quotient \(V_4=Q_8/\{\pm1\}\).  THM-2886 gives the
fixed-source selector parity

```text
(q3 Q0, q11 QA, q7 QAB)=(0,1,1).                       (1)
```

There is a unique character producing (1):

\[
\chi_s(h)=\det(Q_B,h).
\tag{2}
\]

It has

```text
chi_s(QA,QB,QAB)=(1,0,1),
ker(chi_s in Q8)=<QB>={1,QB,-1,-QB}.                    (3)
```

Thus the actual \(q11\to q7\) semantic edge has direction \(Q_B\) and lies
in the kernel of the selector character.  The selector character cannot
also reverse the carry/Hermitian coordinate on that edge.

The character which keeps the first edge \(Q_A\) flat and reverses the
second edge \(Q_B\) is instead

\[
\chi_c(h)=\det(Q_A,h),
\tag{4}
\]

with

```text
target values on (Q0,QA,QAB)=(0,0,1),
edge values on (QA,QB,QAB)=(0,1,1),
ker(chi_c in Q8)=<QA>.                                  (5)
```

No single \(Q_8\to C_2\) character satisfies both (1) and
\(\chi(Q_B)=1\).  More sharply, the literal carry-event edge pattern

```text
(q3->q11, q11->q7, q3->q7)=(0,1,0)                     (6)
```

is not a character because \(Q_A+Q_B=Q_{AB}\).  Equation (6) is
Bockstein curvature, not a vertex action.  This is the first failed
implication and the main negative result of the candidate.

## 2. The two labelled dicyclic shadows

The two characters separately produce isomorphic but differently labelled
semidirect products

\[
G_s=C_{169}\rtimes_{\chi_s}Q_8,\qquad
G_c=C_{169}\rtimes_{\chi_c}Q_8.                         \tag{7}
\]

Both are the dicyclic group \(\operatorname{Dic}_{338}\), using the
convention \(|\operatorname{Dic}_n|=4n\):

\[
\langle a,x\mid a^{676}=1,\ x^2=a^{338},\
                   xax^{-1}=a^{-1}\rangle.              \tag{8}
\]

For the selector shadow, use

```text
c=a^344, <c>=<a^4> ~= C169,
QA=x, QB=a^169, QAB_section=QA QB.                      (9)
```

Then \(Q_A,Q_{AB}\) invert \(c\), while \(Q_B\) fixes it.  Swapping the
labels \(Q_A,Q_B\) gives the carry shadow: \(Q_B,Q_{AB}\) invert \(c\),
while \(Q_A\) fixes it.

There is one central-section convention to keep explicit.  The script's
normal form uses

```text
QAB_section=QA QB.
```

THM-2887's displayed canonical lift uses

```text
QA QB=-QAB_canonical.                                   (10)
```

The two differ only by the central sign, but (10) is load-bearing in the
via/direct clutch below.

For either shadow,

```text
|<c>|=169, |Q8|=8, <c> intersect Q8={1}, |<c>Q8|=1352. (11)
```

The full element-order census is

```text
1^1, 2^1, 4^678, 13^12, 26^12, 52^24,
169^156, 338^156, 676^312.                              (12)
```

There are \(341\) conjugacy classes.  Grouped by
`(element order, class size)`, their census is

```text
(1,1)^1, (2,1)^1, (4,2)^1, (4,338)^2,
(13,2)^6, (26,2)^6, (52,2)^12,
(169,2)^78, (338,2)^78, (676,2)^156.                   (13)
```

Also

\[
Z(G)=\{\pm1\},\qquad
[G,G]=\langle a^2\rangle\cong C_{338},\qquad
G_{\rm ab}\cong V_4.                                   \tag{14}
\]

Every one-dimensional character kills the \(C_{169}\) subgroup.

## 3. The carry-shadow via/direct clutch

In the carry labelling, \(Q_A\) fixes \(c\), \(Q_B\) inverts it, and
THM-2887's section has \(Q_AQ_B=-Q_{AB}\).  Exact multiplication gives

\[
(c^8Q_A)(c^9Q_B)
 =(c^{13}(-1))(c^4Q_{AB,\mathrm{canonical}}).           \tag{15}
\]

The clutch

\[
K=c^{13}(-1)                                            \tag{16}
\]

has

\[
|K|=26,\qquad Q_BKQ_B^{-1}=K^{-1},\qquad K^2=c^{26}.
\tag{17}
\]

Thus the via/direct defect combines the vertical \(C_{13}\) carry and the
quaternionic central sign in one exact order-\(26\) element.  This is the
clean positive carry shadow.  It does not repair (1), because its target
values are \((0,0,1)\).

## 4. Minimal joint action

Retain the two characters on separate cyclic coordinates:

\[
B=(C_{169,s}\times C_{169,c})\rtimes Q_8,               \tag{18}
\]

where

```text
QA:(s,c)->(-s,c),
QB:(s,c)->(s,-c),
QAB:(s,c)->(-s,-c).                                    (19)
```

This is the fibre product of the two \(\operatorname{Dic}_{338}\) shadows
over their shared \(Q_8\):

\[
B\cong G_s\times_{Q_8}G_c,\qquad
|B|=\frac{1352^2}{8}=228488.                            \tag{20}
\]

Its exact anatomy is

\[
Z(B)=\{\pm1\},\qquad
[B,B]=C_{169}^2\times\{\pm1\},\quad |[B,B]|=57122,
\qquad B_{\rm ab}\cong V_4.                             \tag{21}
\]

Unambiguously,

\[
B/\langle-1\rangle
\cong (C_{169}\rtimes C_2)\times(C_{169}\rtimes C_2),
\tag{22}
\]

where each factor on the right has order \(338\).

The full order census is

```text
1^1, 2^1, 4^57798, 13^168, 26^168, 52^8112,
169^28392, 338^28392, 676^105456.                       (23)
```

There are \(14789\) conjugacy classes:

```text
(1,1)^1, (2,1)^1,
(4,338)^2, (4,57122)^1,
(13,2)^12, (13,4)^36,
(26,2)^12, (26,4)^36,
(52,338)^24,
(169,2)^156, (169,4)^7020,
(338,2)^156, (338,4)^7020,
(676,338)^312.                                         (24)
```

Here again `(order,class size)^number of classes` is the convention.

### 4.1 Why two cyclic coordinates are minimal

\(\operatorname{Aut}(C_{169})\) has a unique involution, so an action on one
cyclic coordinate can retain at most one nonzero \(V_4\to C_2\) character.
The two independent characters (2) and (4) require cyclic rank at least two.

The dual \(C_{169}^2\) weight orbits have census

```text
size 1: 1 orbit,
size 2: 168 orbits,
size 4: 7056 orbits.                                   (25)
```

A representation of degree below four can contain at most one nonzero
axis orbit and therefore cannot be faithful on both cyclic directions.
One generic four-orbit spans both axes, so degree four is sharp.

There is also a physical-address version of the same minimality.  The
inherited address supplies only the diagonal

\[
\Delta=\{(n,n):n\in C_{169}\}.                          \tag{26}
\]

It is not \(Q_8\)-stable.  The orbit of \((1,1)\) contains

```text
(1,1), (-1,1), (1,-1), (-1,-1).
```

Since the determinant of \((1,1),(-1,1)\) is \(2\), a unit modulo \(169\),
these points generate all \(C_{169}^2\).  Hence \(B\) is the minimal
\(Q_8\)-stable abelian envelope of the physical diagonal within this
two-character model.  Most off-diagonal states of \(B\) have no physical
realization.

## 5. Typed right-action horn law

The joint carrier types the horn only when the action side and the sign of
the second selector coordinate are retained.  Use **right multiplication**.
In the section gauge \(Q_AQ_B=Q_{AB,\mathrm{section}}\), put

```text
e8=( 8, 8,QA),
e9=(-9, 9,QB).                                         (27)
```

The minus sign is forced: at the intermediate \(Q_A\) state, the selector
coordinate is inverted while the carry coordinate is fixed.  For every
ancestry \(a\), starting from the diagonal q3 state gives

```text
q3  (13a+3, 13a+3,Q0)
 -> (13a+11,13a+11,QA)
 -> (13a+20,13a+20,QAB_section).                        (28)
```

The direct zero-borrow edge in the canonical section is

```text
e4=(4,4,QAB_canonical).                                 (29)
```

Equations (27)--(29) satisfy

\[
e_8e_9=(13,13,-1)e_4.                                  \tag{30}
\]

Thus the joint clutch is exactly the ancestry difference in both cyclic
coordinates plus the quaternionic central sign.  The companion checks all
thirteen ancestry starts.

Using \((+9,+9,Q_B)\), or leaving the action side unspecified, is false.
In the selector shadow alone the analogous right edge is
\(c^{-9}Q_B\), not \(c^{+9}Q_B\).

## 6. q7 filler separation and its boundary

At the quotient level THM-2884 cannot distinguish the two q7 fillers.  In
the joint carrier, their canonical lifts

\[
(13a+7,13a+7,Q_{AB}),\qquad
(13(a+1)+7,13(a+1)+7,Q_{AB})                            \tag{31}
\]

are distinct, square to the central sign, and are conjugate by
\((91,91)\), because \(2\cdot91=13\pmod{169}\).

This is pointwise separation, not invariant separation.  The two fillers
are conjugate; every class function still identifies them.  In fact the
joint \(Q_{AB}\) class has size \(57122\).  A physical consumer must retain
the marked lift or a non-class coordinate.

## 7. Four-channel Pauli realization

Both endpoint fields contain a primitive \(169\)-th root \(s\) and
\(i^2=-1\).  Let

\[
D(s)=\operatorname{diag}(s,s^{-1}),\quad
X=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
Z=\operatorname{diag}(1,-1).
\]

On \(K^2\otimes K^2\), define

\[
\begin{aligned}
c_s&=D(s)\otimes I,&
c_c&=I\otimes D(s),\\
Q_A&=iX\otimes I,&
Q_B&=iZ\otimes X.
\end{aligned}                                           \tag{32}
\]

Then

```text
QA^2=QB^2=-I, QA QB=-QB QA,
QA inverts cs and fixes cc,
QB fixes cs and inverts cc.                             (33)
```

The four torus weights are \((\pm1,\pm1)\).  Exact specialization proves
that the \(169^2\) torus matrices are distinct, their image does not contain
\(-I\), and the four \(V_4\) directions have distinct monomial permutation
supports.  Hence (32) is faithful.  Section 4.1 proves its semisimple degree
four is minimal.

The tensor factors suggest

```text
selector/origin channel tensor Prony/carry channel,
```

but no physical four-channel current with this action is constructed.

### 7.1 Balanced torus versus Hermitianization

Let \(E=\omega^4\), so the raw THM-2886 operator is
\(D=\operatorname{diag}(E,1)\).  In the tensor-weight basis of (32), exact
specialization in both endpoint fields gives

\[
c_s^{130}c_c^{130}
 =\operatorname{diag}(E,1,1,E^{-1}).                   \tag{33a}
\]

This diagonal has a precise Hermitian meaning.  On the ordered rank-one
coordinates

\[
\bigl(U\bar V,\ |U|^2,\ |V|^2,\ V\bar U\bigr),          \tag{33b}
\]

conjugation \(M\mapsto D M D^*\) has weights
\((E,1,1,E^{-1})\).  Thus (33a) identifies the **balanced torus
diagonal** with raw-\(D\) Hermitianization.

It does not identify the two four-dimensional representations.  In the
faithful module (32), the quaternionic center \(z=-1\) acts as \(-I_4\).
Hermitian conjugation kills \(z\), and indeed every common scalar phase,
so its action factors through \(B/\langle z\rangle\).  The former is a
faithful abstract Clifford tensor carrier; the latter is a centered
quadrature carrier.  Only their torus diagonals in (33a) have been
identified.

There is a second exact boundary.  The literal Prony plane
\(e_+\otimes K^2\), spanned by the first two tensor-weight vectors, is sent
by \(Q_A\) to the complementary plane \(e_-\otimes K^2\).  It is therefore
not \(Q_8\)-stable.  This is why the balanced four-coordinate identity
does not itself construct a physical action on the original two-channel
current.

## 8. Coefficient connection and scalar no-go

Put

```text
xi   = primitive 2366th root,
omega=xi^182,
rho  =xi^42.                                            (34)
```

Then \(\rho\) has order \(169\) and

\[
\rho^{13}=\omega^3.                                    \tag{35}
\]

The normalized THM-2868 projective coordinate is

\[
t_r=U_r/V_r=\xi^{955+546r},\qquad
t_{r+1}=\omega^3t_r.                                   \tag{36}
\]

Cyclotomic root inversion gives \(\bar t_r=t_r^{-1}\).  This bar is applied
before finite-field specialization; it is not a nontrivial automorphism of
the prime field.  In the **selector-shadow representation** (9), normalized
projective coherence obeys

```text
c   : t -> rho t,
QA  : t -> -bar(t),
QB  : t -> -t,
QAB : t -> bar(t),
c^13: t -> omega^3 t.                                  (37)
```

The carry-shadow table swaps \(Q_A,Q_B\):

```text
QA  : t -> -t,
QB  : t -> -bar(t),
QAB : t -> bar(t).
```

The THM-2882 lifted line admits the faithful state-dependent gauge

\[
\gamma(a,q)=\rho^{13a+q}\omega^{-3(a+q-3)}.             \tag{38}
\]

For every \(a,q,h\in C_{13}\),

\[
\gamma(L_h(a,q))\omega^{3(h+\kappa(q,h))}
 =\rho^h\gamma(a,q).                                   \tag{39}
\]

The companion verifies all \(2197\) rows in each endpoint field.  This is a
relative-channel gauge, not a common scalar gauge on \(U\oplus V\).

### 8.1 Raw wrapped frame versus normalized C169 frame

THM-2886's raw wrapped edge is

\[
\operatorname{diag}(\omega^4,1)
 =\omega^2\operatorname{diag}(\omega^2,\omega^{-2}).
\tag{40}
\]

The determinant-one factor is exactly

\[
\operatorname{diag}(\omega^2,\omega^{-2})
 =\rho_{\rm std}(a^{104})
 =\rho_{\rm std}(c^{130}).                              \tag{41}
\]

This is the raw Prony frame.  In the gauge-normalized honest
\(C_{169}\) frame, displacement \(+9\) is \(c^9\); it must not be
identified with the raw-frame \(c^{130}\).

The scalar \(\omega^2\) in (40) cannot be globalized as a one-dimensional
character.  In each dicyclic shadow \(c^{130}\in[G,G]\), and in the joint
carrier both \(C_{169}\) circles lie in \([B,B]\).  All one-dimensional
characters therefore kill the determinant-one carry rotation.  A global
linear lift with the same projectivization would have to differ by such a
character, so the raw scalar debt remains.

This gives the sharp positive/negative boundary:

- projective and normalized-coherence relative phases have exact dicyclic
  and four-channel realizations;
- the raw common scalar does not extend as a character; and
- recombination \(U+V\) still has no charged scalar descent.

## 9. Exact consequence and non-consequence

The candidate proves:

1. the selector and carry/Hermitian requirements are two distinct
   \(Q_8\to C_2\) characters;
2. no single character types both, and the literal carry pattern is
   Bockstein curvature rather than a character;
3. each character has an exact \(\operatorname{Dic}_{338}\) shadow;
4. the carry shadow has the order-\(26\) clutch (15)--(17);
5. their minimal joint action is the fibre product \(B\) in (18);
6. the complete group, conjugacy, derived, quotient, and weight anatomy;
7. the exact typed right-action horn law (27)--(30);
8. pointwise, but not class-invariant, q7 filler separation;
9. a faithful minimal four-dimensional Pauli representation; and
10. the exact balanced-torus/Hermitian diagonal (33a), together with its
    central-quotient and non-stable-two-plane boundary; and
11. the exact Prony/projective match and scalar-character no-go.

It does **not** prove:

- that either abstract \(Q_8\) action operates on physical horn packets;
- that the physical diagonal address realizes the off-diagonal states of
  \(B\);
- that the THM-2886 marked subcopy descends to the canonical unmarked
  current;
- a lawful physical Hermitian polarization or four-channel observable;
- positivity, a row exclusion, or LRC(14).

The cheapest next physical test is the right-action two-traversal law.
Construct the \(e_9=(-9,+9,Q_B)\) action on one common marked packet and
test whether the via endpoint differs from the direct endpoint by
\((13,13,-1)\).  This simultaneously tests the second character, ancestry
carry, and quaternionic sign before scalar recombination.

## 10. Reproduction

Run

```text
python3 .scratch/lrc_dicyclic_hermitian_envelope_20260729/lrc14_dicyclic_reverse_action_joint_carrier_thm2889.py
python3 -O .scratch/lrc_dicyclic_hermitian_envelope_20260729/lrc14_dicyclic_reverse_action_joint_carrier_thm2889.py
```

Normal and optimized stdout byte-match

```text
.scratch/lrc_dicyclic_hermitian_envelope_20260729/lrc14_dicyclic_reverse_action_joint_carrier_thm2889.out.
```

The script contains no executable Python `assert`; all load-bearing checks
use `require`.
