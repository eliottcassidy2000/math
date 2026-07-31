# Carry/selector product extension and recombination quotient obstruction

**Status:** FINITE-EXACT SCRATCH CANDIDATE; hostile audit PASS, with
scope qualifications incorporated below.

**Scope:** abstract residue-path, selector, and coefficient
representations inherited from THM-2882.  No physical origin coupling,
QA/QAB attachment, row exclusion, or LRC(14) conclusion is constructed.

## 1. Result

The order-thirteen carry and the Boolean origin-selector defect are
independent in a stronger sense than disjoint edge support.

- The carry is the nonsplit \(C_{169}\to C_{13}\) extension.
- Selector parity is an exact \(\mathbf F_2\)-valued coboundary.
- After the two-origin truth is retained, the distinguished positive
  selector is a constant mismatch section and needs no new state.
- If selector parity must be independently variable and faithful, the
  smallest joint carrier has

  \[
  13\text{ residues}\times13\text{ carry states}\times2
  =338
  \]

  states and is the direct-product torsor \(C_{169}\times C_2\).
- Under the same retained-residue, faithful-carry, and independent free
  commuting-selector hypotheses, transport of both selector bits and
  their signed orientation requires

  \[
  C_{169}\times V_4,\qquad 169\cdot4=676
  \]

  states.
- The exact 449-sheet horn carries a second, semantic copy

  \[
  H=\ker(e+u+v)\leq \mathbf F_2^3,\qquad H\cong V_4.
  \]

  Its three old word states are \(H\setminus\{0\}\), and the newly
  reconstructed all-out word supplies \(0\) at the stepped \(q3\)
  endpoint.  This is a complete support section, not yet a free physical
  \(V_4\)-action and not an identification with the selector \(V_4\).

The parity quotient preserves whether the signed coefficient is zero,
but loses its sign.  No quotient retaining the character-three vertical
translation \(T\) induces a one-dimensional linear action on the
recombined coefficient \(U+V\).  The split pair \((U,V)\), or the charged
channel \(U\) alone, is the minimal surviving coefficient object.

## 2. Truth, selector, and edge cocycles

Encode E3 by \(1\) and its full complement by \(0\).  The two-origin
truth is

\[
t(q)=(t_0(q),t_1(q)),
\]

with

```text
q0,q11: (1,1)
q3:     (1,0)
other:  (0,0).
```

Thus

\[
d(q)=t_0(q)\mathbin{\mathtt{XOR}}t_1(q)=\delta_{q,3}.
\]

The THM-2882 zero-origin-positive Boolean selector is

\[
s^+(q)=(t_0(q),\,1\mathbin{\mathtt{XOR}}t_1(q)).
\]

Its absolute selector parity is

\[
p(q)=s^+_0(q)\mathbin{\mathtt{XOR}}s^+_1(q)
    =1\mathbin{\mathtt{XOR}}\delta_{q,3}.               \tag{1}
\]

For a reduced positive edge \(q\xrightarrow{h}q+h\bmod13\), define

\[
\sigma(q,h)=p(q)\mathbin{\mathtt{XOR}}p(q+h\bmod13).
\]

This is the exact coboundary \(\delta p\), so all \(2,197\) triples obey

\[
\sigma(q,h)\mathbin{\mathtt{XOR}}\sigma(q+h,k)
=\sigma(q,h+k\bmod13).                                  \tag{2}
\]

Its unit-edge support is exactly

```text
q2 -> q3,   q3 -> q4.
```

The carry transition is

\[
\kappa(q,h)=\left\lfloor\frac{q+h}{13}\right\rfloor
\]

and satisfies

\[
\kappa(q,h)+\kappa(q+h\bmod13,k)
=\kappa(q,h+k\bmod13)+\left\lfloor\frac{h+k}{13}\right\rfloor .
                                                               \tag{3}
\]

Its unit-edge support is only \(q12\to q0\).  Two minimal witnesses show
that neither coordinate determines the other:

```text
q0 --3--> q3:  sigma=1, kappa=0;
q12--1--> q0:  sigma=0, kappa=1.
```

The first is exactly THM-2880's right-selector-provenance edge.  Its
selector changes while the \(u_3\) event count is zero.

## 3. The 338-state parity lift

Use the base-\(13\) set coordinates for \(C_{169}\), and adjoin an
independent parity bit:

\[
(a,b,q),\qquad a,q\in C_{13},\quad b\in\mathbf F_2.
\]

Define

\[
\widetilde L_h(a,b,q)
=\bigl(a+\kappa(q,h),\,
       b+\sigma(q,h),\,
       q+h\bmod13\bigr).                                \tag{4}
\]

The script checks all \(57,122\) initial-state composition laws.  The
central carry in (3) changes \(a\), while the exact selector cocycle in
(2) contributes no central defect.

The gauge

\[
c=b\mathbin{\mathtt{XOR}}p(q)                           \tag{5}
\]

is constant under every edge.  Hence the positive unit action consists
of two disjoint \(169\)-cycles, one for each \(c\).  Adjoining the free
selector flip

\[
J(a,b,q)=(a,b+1,q)
\]

connects the two cycles, commutes with the positive carry generator, and
acts regularly on all \(338\) states.  This is

\[
C_{169}\times C_2\cong C_{338}
\]

as a group, equipped with the residue-path gauge (5).

The state lower bound is conditional and exact: target residue is
retained, the carry \(T\)-orbit is faithful of size \(13\) in every
residue fibre, and an independent free parity flip needs two states over
each carry state.  Thus any deterministic faithful joint carrier of
these coordinates has at least \(13\cdot13\cdot2=338\) states.

If the particular external section \(s^+\) is frozen rather than made
independently variable, no new fibre is needed.  Indeed

\[
x(q)=s^+(q)\mathbin{\mathtt{XOR}}t(q)=(0,1)
\]

is constant.  This explains why THM-2882 itself has \(169\), not \(338\),
decorated states; it does not physically construct the missing selector
flip \(J\).

## 4. Why the product is direct

There are two abstract order-\(338\) semidirect possibilities once the
normal \(C_{169}\) and an order-two complement are fixed.  The
order-two action on \(C_{169}\) is either identity or inversion:

```text
u^2=1 mod169  iff  u in {1,-1}.
```

Inversion is inadmissible in the inherited typed gauge, where \(q0\) is
the fixed identity and the positive generator is preserved.  It reverses
that orientation and sends exceptional truth residue \(q3\) to ordinary
residue \(q10\):

```text
t(3)=(1,0),    t(10)=(0,0).
```

Conversely, \(C_{169}\) cannot act nontrivially on \(V_4\):

\[
\operatorname{Aut}(V_4)=GL_2(\mathbf F_2)\cong S_3
\]

has element orders \(1,2,3\), none divisible by \(13\).  A
\(V_4\)-action on \(C_{169}\) has four possible inversion characters,
and only the trivial character preserves the positive generator.
Therefore the gauge-fixed typed joint action is direct, not dihedral or
another semidirect regrading.  This statement does not classify arbitrary
affine reoriginings that move the chosen identity.

This is also visible cohomologically.  The normalized exact ranks over
\(\mathbf F_{13}\) are

```text
V4: dim C1=3, dim C2=9, rank B2=3, dim Z2=3, H2=0;
C13 carry: coboundary rank=11, carry-augmented rank=12.
```

The selector \(V_4/\mathbf F_2\) side cannot intrinsically host the
order-thirteen Bockstein class.  The class remains entirely on the
\(C_{169}\) factor.

## 5. Full two-origin selectors require \(676\) states

Let the selector itself vary through

\[
s=(s_0,s_1)\in V_4.
\]

Transport it by the truth boundary:

\[
s' = s\mathbin{\mathtt{XOR}}t(q)
         \mathbin{\mathtt{XOR}}t(q+h\bmod13).           \tag{6}
\]

Then the mismatch

\[
x=s\mathbin{\mathtt{XOR}}t(q)                           \tag{7}
\]

is constant.  The unit action has four disjoint \(C_{169}\)-cycles,
indexed by \(x\in V_4\); adjoining the two independent selector flips
gives the regular direct-product torsor

\[
C_{169}\times V_4
\]

on \(676\) states.

The signed endpoint amplitude attached to mismatch is

\[
A(x)=x_1-x_0\in\{-1,0,1\}.                              \tag{8}
\]

The parity quotient \(V_4\to C_2\) has classes

```text
even: {(0,0),(1,1)}   -> zero support;
odd:  {(0,1),(1,0)}   -> nonzero support.
```

It preserves support but identifies amplitudes \(+1\) and \(-1\).
Retaining the three values in (8) also fails to give an equivariant
quotient: \((0,0)\) and \((1,1)\) both have amplitude zero, but
translation by \((0,1)\) sends them to \((0,1)\) and \((1,0)\), whose
amplitudes are opposite.  Thus, under the retained-residue,
faithful-carry, and independent free selector-operation hypotheses:

- \(338\) states suffice for support parity;
- \(676\) are required for faithful selector operations and signed
  orientation;
- the frozen positive section uses \(169\) only by choosing one
  external \(V_4\) mismatch.

## 6. The 449-sheet semantic diagonal is a second \(V_4\)

Let

\[
(e,u,v)\in\mathbf F_2^3
\]

record the macro E3 bit and the `slot7,slot8` outer-word bits.  The exact
diagonal condition is the kernel of the ambient character

\[
\chi(e,u,v)=e+u+v,\qquad
H=\ker\chi.                                               \tag{9}
\]

Thus \(H\) is the graph \(e=u+v\), and projection to \((u,v)\) identifies
it with \(V_4\).  The exact semantic states are

```text
Qempty = 000,    QB = 101,    QA = 110,    QAB = 011.
```

The old horn states satisfy

\[
QB+QA=QAB
\]

and are exactly the three nonzero points of \(H\).  Hence they form the
punctured kernel \(H^\times=H\setminus\{0\}\), not a subgroup or a coset.
Equivalently, they are the three points of
\(\mathbf P^1(\mathbf F_2)\), permuted by
\(\operatorname{Aut}(H)\cong S_3\).

The independent reconstruction starts from the full circle and subtracts
the guard, five unit, owner, and two outer combs to build the previously
omitted all-out word.  It gives

```text
global Qempty counts q0..q12:
(0,66100,66101,66104,66107,66107,65653,
 0,66108,66106,66106,66103,66100);

on every one of the 449 horn labels:
q1,q2,q3,q4,q5 are Qempty;
q7 at a+1 is QAB.
```

At \(q3\), the stepped origin has complement truth and state \(000\), so
it completes \(H\).  The zero origin has E3 truth and state \(100\), in
the opposite coset.  On the four seam checkpoints the two-origin
character values are

```text
q0/QB:     (0,0)
q3/Qempty: (1,0)
q11/QA:    (0,0)
q7/QAB:    (0,0).
```

Therefore

\[
\chi_0(q)+\chi_1(q)=\delta_{q,3}.                         \tag{10}
\]

The diagonal gives a canonical origin-odd \(q3\) leg, but it is
origin-even at \(q11\) and \(q7\), where the signed coefficients still
cancel.

There is a useful seam-local coincidence:

\[
p(q)=1+\delta_{q,3}
    =1_{H^\times}(h(q))
    =u\mathbin{\mathrm{OR}}v
    =u+v+uv                                                \tag{11}
\]

on \(q=0,3,11,7\).  This makes the three live words a canonical
support-parity object.  It is nonlinear: \(1_{H^\times}\) is not a
character, since two nonzero elements can add to a third nonzero element.
Nor is `(11)` an all-\(q\) rule.  At \(q1,q2,q4,q5\), the same 449 sheets
are \(Qempty=0\), while the THM-2882 desired parity is one.

More strongly, the script exhausts all \(2^8=256\) Boolean observables on
\(\mathbf F_2^3\).  Exactly \(128\) can distinguish the two \(q3\)
states, but none can distinguish the origins at \(q11\) or \(q7\),
because those origins have literally the same \((e,u,v)\).  Thus no
semantic-only quotient or nonlinear Boolean rule can supply the missing
signed current there.

Adjoining the endpoint-origin bit \(o\) gives the smallest set-theoretic
repair.  For example,

\[
\ker(\chi+o)\leq\mathbf F_2^4                              \tag{12}
\]

selects one origin at \(q11\) and \(q7\).  Equation `(12)` is only a
precise target: it explicitly inserts the independent origin coordinate
and is not yet typed by a source/QA/QAB packet operation.

This sharpens the \(C_{169}\times V_4\) classification.

- \(H\cong V_4\) is a semantic word carrier.
- The \(V_4\) of Section 5 is the two-origin selector carrier.
- They are abstractly isomorphic but cannot be identified by the current
  data: \(H\) forgets origin exactly where sign is needed.
- A freely variable semantic \(H\) would give a \(676\)-state
  \(C_{169}\times H\) carrier, but the 449-sheet result supplies only a
  distinguished support section.  A conditional **replacement** model
  that uses an eight-state central-sign extension \(E\to H\) in place of
  the selector \(V_4\) has \(169\cdot8=1352\) states.  As a set this
  adjoins one origin/sign bit; the group law could be split
  \(H\times C_2\) or nonsplit.  Section 7 identifies the natural nonsplit
  law.  This replacement does not retain both copies of \(V_4\).
  Retaining semantic \(H\), selector \(V_4\), and the origin bit all
  faithfully and independently would instead cost
  \(169\cdot4\cdot4\cdot2=5408\) states.

The semantic increments also fit the direct-product bookkeeping:

```text
source QB(a) -> q11 QA(a):       increment QAB;
q11 QA(a)   -> q7 QAB(a+1):     increment QB plus one carry.
```

The second arrow is therefore a skew edge in
\(C_{169}\times H\), not a pure semantic translation.  Object-level
closure and 449-fold support realization are exact; complete-packet
intertwiners implementing a free \(H\)-action remain open.

## 7. The punctured \(V_4\) has a canonical quaternionic sign lift

The nonlinear seam parity in `(11)` is not an arbitrary failure of
linearity.  On \(H\cong\mathbf F_2^2\), put

\[
Q(u,v)=u+v+uv=1_{(u,v)\ne0}.                             \tag{13}
\]

Its polarization is

\[
Q(x)+Q(y)+Q(x+y)
 =x_1y_2+x_2y_1
 =\det(x,y).                                             \tag{14}
\]

Thus \(Q\) is the anisotropic, Arf-one quadratic refinement of the unique
symplectic form on \(H\).  This is the exact dual of THM-2779's
characteristic-two boundary

\[
Q_0(u,v)=uv,
\]

the Arf-zero refinement arising from \(D_8\).  Both have determinant
polarization, but their values and central extensions differ:

```text
Q_0 values on (00,01,10,11): (0,0,0,1) -> D8;
Q   values on (00,01,10,11): (0,1,1,1) -> Q8.
```

For one forward section gauge, the \(Q_8\) cocycle is

\[
c(x,y)=x_1y_1+x_2y_2+x_1y_2.                            \tag{15}
\]

On \((x,z)\in H\times C_2\), use

\[
(x,z)(y,w)=(x+y,z+w+c(x,y)).                             \tag{16}
\]

Then every nonzero semantic direction has square equal to the central
sign.  With \(QA=i\), \(QB=j\), and \(QAB=k\), the distinguished order is

```text
i*j=-k,       j*i=+k.
```

Changing the section by \(Q\) reverses which ordered product pays the
sign; the commutator and polarization `(14)` are gauge invariant.

The script exhausts all normalized \(C_2\)-valued cocycles on \(V_4\):

```text
normalized cocycles:                         16
with determinant commutator:                  8
quadratic refinements:     four, two gauges each
all-three-nonzero square law: one Arf-one class.
```

Thus \(Q_8\), rather than \(H\times C_2\) or \(D_8\), is the unique
central-sign extension with the THM-2884 support parity as its square
law.  Its exact symmetry anatomy is

\[
1\longrightarrow V_4=\operatorname{Inn}(Q_8)
\longrightarrow \operatorname{Aut}(Q_8)\cong S_4
\longrightarrow S_3\longrightarrow1.                   \tag{17}
\]

The enumeration gives automorphism order census

```text
1^1, 2^9, 3^8, 4^6,
```

kernel size four and quotient size six.  This is the literal
\(S_4/V_4\cong S_3\) object recurring in the repo, now attached to the
449-sheet horn rather than imported by analogy.  Unlike the THM-2779
\(D_8\) refinement \(uv\), \(Q\) is preserved by all of
\(GL_2(\mathbf F_2)\cong S_3\).

There is also an exact local carry match.  The two consecutive semantic
edge directions are \(QA,QB\), so

\[
Q(QA)+Q(QB)+Q(QAB)
=\det(QA,QB)=1
=\left\lfloor\frac{8+9}{13}\right\rfloor.                \tag{18}
\]

The Arf polarization therefore detects precisely the central event on
the \(q3\to q11\to q7\) triangle.

This match is sharp and local.  Modulo two, the full \(C_{13}\) carry
cocycle has the unique normalized primitive

\[
r(h)=h\bmod2,\qquad
\delta r(h,k)=\left\lfloor\frac{h+k}{13}\right\rfloor.    \tag{19}
\]

If a base-independent semantic label
\(\ell:C_{13}\to H\) globalized `(18)`, then
\(Q(\ell(h))=r(h)\).  Since \(Q\) vanishes only at zero, every even \(h\)
would have \(\ell(h)=0\), contradicting the observed

```text
ell(8)=QA,       ell(4)=QAB.
```

So the quaternionic sign cannot absorb or identify the order-thirteen
Bockstein globally.  Also, exhaustive \(S_3\)-equivariance shows there is
no map from the three semantic directions to a two-valued sign set when
\(S_3\) acts on signs by its sign character.  With trivial sign action,
only the two constant maps survive.  A physical endpoint orientation
still requires an extra lift/choice.

Every nontrivial subgroup of \(Q_8\) contains its center, so its minimal
faithful permutation degree is eight: the four-state semantic quotient
necessarily forgets the sign.  Since
\(\operatorname{Aut}(Q_8)\cong S_4\) has no order-thirteen element, the
positive gauge-fixed joint candidate is the direct product

\[
C_{169}\times Q_8,\qquad 169\cdot8=1352.                 \tag{20}
\]

This is a structured candidate, not a physical construction.  No horn
sheet has yet been lifted to a signed quaternion unit, and no \(q11/q7\)
coefficient follows from `(20)`.

## 8. Recombination quotient no-go

Let \(T\) be the vertical ancestry translation.  It fixes \(q\), truth,
and every selector coordinate.  On the two Prony channels it acts as

\[
D_T=\operatorname{diag}(\omega^3,1),\qquad \omega^3\ne1. \tag{21}
\]

The recombination map is

\[
R(U,V)=U+V.
\]

If a one-dimensional scalar action \(\lambda(T)\) descended through
recombination, then

\[
R D_T=\lambda(T)R
\]

would force simultaneously

\[
\lambda(T)=\omega^3,\qquad \lambda(T)=1,
\]

a contradiction.  Equivalently, the kernel of \(R\), spanned by
\((1,-1)\), is not \(D_T\)-invariant.

This obstruction survives every selector quotient.  An order-two
selector character cannot cancel an order-thirteen carry character:

\[
\mu_2\cap\mu_{13}=\{1\}.
\]

For any group quotient, the image of the prime-order subgroup
\(\langle T\rangle\) is either trivial or still order \(13\).

- If it remains nontrivial, recombination has no induced scalar action.
- If it is killed, the character-three Bockstein coordinate is lost.
- Projecting to \(U\), or retaining the split pair \(U\oplus V\), is
  lawful but is not the recombined physical scalar.

The exact edge census independently recovers THM-2882's boundary:

```text
scalar E-transport of U+V: 144 failures, 25 E=1 survivors.
```

The no-go is linear and one-dimensional.  It does not rule out a
nonlinear observable or a new physical operation that natively retains
the two-channel state.

## 9. Relation to the remaining LRC gate

THM-2880 closes the affine \(q0\to q3\) one-fibre route because it cannot
preserve right-selector provenance and source support together.  The
present product classification locates exactly what that theorem leaves
open:

```text
an independent-origin F2/V4 coupling, commuting with the C169 carry.
```

Such a coupling must do more than supply a phase.  It must:

1. realize a physical selector flip \(J\) while preserving the fixed
   source;
2. lift the exact semantic diagonal through an origin-coupled rule such
   as \(\chi+o\), rather than another observable of \((e,u,v)\).  The
   \(Q_8\) extension gives the canonical multiplication law for such a
   central sign, but not its physical sheetwise lift;
3. align the \(449\) semantic
   \(QA(q11,a)\to QAB(q7,a+1)\) sheets with the charged \(U\)-line; and
4. retain either the full split coefficient pair or a new lawful
   two-channel observable, since scalar recombination cannot carry \(T\).

The old unrestricted sheetwise selector search is now too broad.  The
cheapest decisive test is an origin-resolved audit of the tilted kernel
\(\chi+o\) on the complete source and endpoint packet.  Every rule that
factors through \((e,u,v)\) alone has already been excluded at \(q11/q7\).
A construction that records only the punctured-\(H\) support parity can
never certify signed orientation.

## 10. Reproduction

```bash
python3 .scratch/lrc_carry_selector_product_extension_20260729/audit.py
python3 -O .scratch/lrc_carry_selector_product_extension_20260729/audit.py
python3 .scratch/lrc_carry_selector_product_extension_20260729/semantic_diagonal_audit.py
python3 -O .scratch/lrc_carry_selector_product_extension_20260729/semantic_diagonal_audit.py
python3 .scratch/lrc_carry_selector_product_extension_20260729/quaternionic_selector_lift.py
python3 -O .scratch/lrc_carry_selector_product_extension_20260729/quaternionic_selector_lift.py
```

The scripts contain no executable Python `assert`.  Both normal and
optimized outputs byte-match their stored scratch outputs.

```text
audit.py                     1d029b86f4241e856abb5a8cbecd0d3a2f3ad589f5bc8d99d8f1983252db0a48
audit.out                    94da55063ae614a9d4cb274c53f3c1e4a04e9fe12a896877d1e1b4c0591a0b12
semantic_diagonal_audit.py   fe174cfd76229ab9c07b0db23e55c97c5818750c6740701aed292df238e8e6f4
semantic_diagonal_audit.out  e94f6b6a15fd67db782ad3975e76dd799a605b8dd7dacddb86d8dc798cc8fb8f
quaternionic_selector_lift.py  41b54fa12ea569da7b2d1f1d705879f76b09f467a1fadd22929923a53fe08fc9
quaternionic_selector_lift.out b8825a66c1050783b1609588cd32dc5130d8575695aa05bbef7ee6c2abb098ac
```

Hashes use LF-normalized bytes.  The exact claims received a hostile
scope audit; canonical promotion remains a separate decision.
