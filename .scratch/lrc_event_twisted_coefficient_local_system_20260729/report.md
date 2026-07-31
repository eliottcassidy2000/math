# Event-twisted all-\(q\) coefficient lift

**Status:** FINITE-EXACT SCRATCH CANDIDATE; no LRC(14) consequence.

**Closest proved mechanisms:** THM-2874 (truth-polarized endpoint
coefficient charts and their flat frequency seam), THM-2878 (the unique
positive \(u_3:D\to S\) exit computes the carry), and THM-2851 (the
abstract \(C_{169}\) lift and the \(449\)-sheet
\(QA(q11,a)\to QAB(q7,a+1)\) support incidence).  THM-2880 is the
nearest selector-provenance hostile: it closes the source-induced
one-fibre route but explicitly leaves broader independent-origin
couplings open.

**Hostile inherited boundary:** fixed-frequency target translation has
zero endpoint phase; q11 cancels between the two endpoint origins; q7 is
absent under E3.  A coefficient construction must therefore be tested
after restoring origins and the fixed source, not only on one endpoint
atom.

## 1. The thirteen truth-polarized atoms

Let \(J_q\) be the translate by \(qT/13\) of the selected target atom.
At the zero endpoint character, the literal nine-factor words are

```text
q0  SSSSSSSSD       q7  DSSSSDSSD
q1  SSDSSSSSD       q8  DSSSSSSSD
q2  SSDSSSSSD       q9  SSSSDSSSD
q3  SSSSSSSSD       q10 SSSSDSSSD
q4  SDSSSSSSD       q11 SSSSSSSSD
q5  DDSSSSSSD       q12 SSSDSSSSD
q6  DSSSSDSSD.
```

The full endpoint pattern `PAT_E3` contains all of \(J_q\) exactly for

```text
q in {0,3,11}.
```

Its full set-theoretic complement contains all of every other \(J_q\).
Thus choosing E3 at \(0,3,11\) and not-E3 elsewhere gives thirteen
full, translated, weight-one atoms.  This is a zero-chart Boolean
coefficient atlas.  It is not yet a common two-origin selector.

Every atom has the same endpoint parameters
\((\alpha_L,\alpha_R,\lambda_L,\lambda_R)\).  In the two certified
finite fields, the common Prony nodes are respectively

```text
(318667610960925913, 240826149255918765)
(330274131806661097, 567852948617574198).
```

For section \(r\), use the formal multiplier

\[
 m_r=1+42r
\]

and the existing unit repair \(e_4=1,e_8=2\), with all other \(e_r=0\).
All twenty-six \(m_r+e_r,m_r+e_r+1\) and all twenty-six corresponding
\(12+26m\) frequencies are units modulo \(91\).  Two-node Prony
splitting on every one of the \(13\cdot13\) atom/section pairs gives

\[
 U_r=U_0\omega^{3r},\qquad V_r=V_0,
\]

independently of \(q\).  Both channels are nonzero.
The nonzero THM-2868 common source scalar \(P_{\rm src}\) is retained
algebraically in both \(U_r\) and \(V_r\); no all-\(q\) fixed-source
packet operation is inferred from this reuse.

## 2. The lifted \(C_{169}\) coefficient system

Write each \(n\in C_{169}\) uniquely in the base-\(13\) set chart

\[
n=13a+q,\qquad 0\leq a,q<13.
\]

This is not an identification with the product group
\(C_{13}\times C_{13}\).  The carry in \(L_1(0,12)=(1,0)\) is exactly
what makes the group law nonsplit.  Assign the lifted state \(n\) the
frequency section

\[
 r(a,q)=a+q-3\pmod {13}.
\]

The decorated vertex is

\[
 \mathcal C_{a,q}=(J_q,\ U_{r(a,q)}).
\]

The support \(J_q\) distinguishes \(q\), and the thirteen nonzero,
distinct values \(U_r\) distinguish \(a\).  Hence the \(169\) decorated
vertices are pairwise distinct.

Let

\[
 \kappa(q,h)=\#\{u_3:D\to S\text{ exits on the positive }h
 \text{-step path from }q\}
              =\left\lfloor\frac{q+h}{13}\right\rfloor .
\]

The equality is independently reconstructed from the full literal word,
and \(u_3:D\to S\) is the unique one of the eighteen oriented
single-factor events with this property on the zero chart.  Define

\[
 L_h(a,q)=
 \left(a+\kappa(q,h),\,q+h\pmod {13}\right).
\]

If \(q'=q+h-13\kappa(q,h)\), then

\[
\begin{aligned}
r(L_h(a,q))-r(a,q)
 &=h-12\kappa(q,h)\\
 &\equiv h+\kappa(q,h)\pmod {13}.
\end{aligned}
\]

Therefore the exact coefficient ratio is

\[
 E(q,h)
 =\frac{U_{r(L_h(a,q))}}{U_{r(a,q)}}
 =\underbrace{\omega^{3h}}_{F(q,h)}
  \underbrace{\omega^{3\kappa(q,h)}}_{\tau(q,h)}.       \tag{1}
\]

Here \(F\) is the old flat frequency-section transition and \(\tau\) is
the event twist.  All \(169\) reduced edges and all \(2,197\) lifted
state transports were checked in both fields.

For \(s=h+k\bmod13\) and
\(\epsilon(h,k)=\lfloor(h+k)/13\rfloor\), the event identity is

\[
\kappa(q,h)+\kappa(q+h\bmod13,k)
=\kappa(q,s)+\epsilon(h,k).
\]

Consequently all \(2,197\) reduced compositions satisfy

\[
E(q+h\bmod13,k)E(q,h)
=\omega^{3\epsilon(h,k)}E(q,s).                         \tag{2}
\]

All \(28,561\) versions with an initial \(a\) were also checked.  Thus
\(E\) is a projective system on the reduced \(C_{13}\) residue
groupoid, with the character-three carry curvature, and an honest lifted
coefficient system after adjoining the \(C_{13}\) ancestry fibre.  It
must not be called a flat local system on the reduced \(C_{13}\) object.

The flat factor can be gauged away.  With

\[
\phi(q)=\omega^{-3(q-3)},
\]

all \(169\) edges satisfy

\[
\tau(q,h)\phi(q)=\phi(q+h\bmod13)E(q,h).                 \tag{3}
\]

Hence the lifted system is exactly gauge-equivalent to the
character-three ancestry line.

## 3. Reduced seam curvature and the vertical clutch

For

```text
q3 --8--> q11 --9--> q7
q3 --------4-------> q7,
```

the event counts are \((0,1,0)\).  The flat frequency factors have
reduced seam curvature \(1\), while both the event and event-twisted
factors have

\[
\frac{E(11,9)E(3,8)}{E(3,4)}=\omega^3.
\]

This is exactly the THM-2851 character-three Bockstein value.  It is not
the holonomy of a closed path in the honest lifted system: the via path
ends at \((a+1,7)\), while the direct path ends at \((a,7)\).  The value
\(\omega^3\) is the vertical \(T\)-clutch between those endpoints, or
equivalently the projective curvature after reducing away \(a\).
Appending \(T^{-1}\) closes the lifted loop and contributes
\(\omega^{-3}\), so the honest lifted loop has holonomy \(1\).

The abstract, event-decorated rank-one coefficient transition invoice is
therefore closed.  The physical clutch invoice is not.

## 4. Channel typing

The rank-one system above lives on the charged split-left channel
\(U_r\), equivalently the centered line \(c_r-A\).  The other Prony
channel is trivial:

\[
(U,V)\longmapsto(EU,V).
\]

It is not scalar transport of the recombined coefficient \(S=U+V\).
Exact exhaustion gives failure of \(S'\!=ES\) on \(144\) of the \(169\)
reduced edges (the \(25\) surviving cases are precisely the edges with
\(E=1\)).  Thus the result is not a scalar action on the raw
two-node/full-current coefficient.

## 5. The two-origin selector dichotomy

The zero endpoint origin has E3 support

```text
{0,3,11},
```

whereas the stepped origin has E3 support

```text
{0,11}.
```

All Boolean selector choices admit an exact classification.  Encode
E3 by \(1\), complement by \(0\).  Let \(t_o(q)\) be the block containing
\(J_q\) at origin \(o\), and let \(s_o(q)\) be the selected block.  The
origin contributes precisely when \(s_o=t_o\).  Therefore the signed
coefficient is nonzero exactly when

\[
s_0\mathbin{\mathtt{XOR}}s_1
\mathbin{\mathtt{XOR}}t_0\mathbin{\mathtt{XOR}}t_1=1.    \tag{4}
\]

Here the truth-origin parity is

\[
t_0(q)\mathbin{\mathtt{XOR}}t_1(q)=\delta_{q,3}.          \tag{5}
\]

Equations (4)--(5) sharpen the gate:

- Every origin-independent same-label rule \(s_0=s_1\) has signed
  support exactly \(\{3\}\).  In particular, q11 E3 coefficients and
  q7 complement coefficients are individually nonzero and equal at the
  two origins, hence cancel.  This was replayed on all twenty-six lawful
  samples.
- The everywhere-full rule \(s_o=t_o\) is origin-dependent at q3 and
  cancels identically.
- Any nonzero q7 or q11 coefficient requires an origin-odd selector
  \(s_0\mathbin{\mathtt{XOR}}s_1=1\) there.

Within the Boolean E3/complement selector class, there is a unique
zero-origin-positive selector that retains the full zero-origin atom and
the empty stepped-origin block at every q.  On the seam it is

```text
q3:  (E3,E3)
q11: (E3,not-E3)
q7:  (not-E3,E3).
```

Its signed coefficient is exactly the all-\(q\) zero-chart coefficient
used in Sections 1--3.  Thus the all-\(q\) atlas can be recast as a
precise origin-polarized coefficient difference, and this removes q11
cancellation and supplies q7 coefficient support.  It does so by
projecting away the stepped origin with a target-dependent truth rule.
No source/QA-to-QAB operation currently types that origin-odd selector.

The selector parity is independent of the carry event.  For the unique
zero-origin-positive rule,

\[
p(q)=s_0(q)\mathbin{\mathtt{XOR}}s_1(q)
    =1\mathbin{\mathtt{XOR}}\delta_{q,3}.
\]

It toggles only on \(q2\to q3\) and \(q3\to q4\).  At the unique carry
edge \(q12\to q0\), \(p\) remains \(1\).  Conversely, \(q0\to q3\)
changes \(p\), while \(\kappa(0,3)=0\), exactly as pinned in THM-2880.
Thus neither coordinate determines the other: the \(u_3\) event twist
cannot manufacture or type the origin-odd selector.  The polarization
belongs to THM-2880's still-open broader signed/independent-origin
coupling class, not its closed source-induced fibre route.

Multiplying the inherited origin-independent current by the event factor
cannot change its q11/q7 zeros.  The new exact missing sidecar is not
another phase: it is a lawful origin-odd selector parity at q7 and q11.

## 6. Why the twist is external to the physical current

At one fixed frequency section, raw target translation has coefficient
ratio \(1\) on all \(169\) \((q,h)\) edges: both the \(12\)- and
\(26\)-endpoint translation phases vanish exactly.  The nontrivial
factor \(E\) arises by reindexing the external frequency measurement
according to \(r(a,q)\), then adjoining the directed event clutch.  It
is not observed as a raw translation phase of one current.

THM-2851/THM-2835 still provide the exact support incidence

\[
449:\quad QA(q11,a)\longrightarrow QAB(q7,a+1).
\]

But the present \(a\in C_{13}\) is an externally chosen first-carry
coordinate.  No coefficient-to-\(\Lambda\) basepoint or sheetwise
alignment is constructed.  The pinned marked-gauge hostile says that
retaining the q7 endpoint copy moves the fixed source into the all-safe,
no-word chamber.  Together with the common-selector q11/q7
cancellations—and the untyped origin-odd repair above—this leaves the
physical \(QA/QAB\) attachment open.

## 7. Connection contract

```text
source:
  thirteen zero-chart E3/complement atoms, THM-2868 split frequencies,
  and THM-2878's positive u3 exit event;

target:
  the first-carry C169 coefficient lift;

map:
  (a,q) -> frequency section r=a+q-3, with L_h using the event count;

preserved:
  target residue, full endpoint atom, common Prony nodes, lawful
  91-unit samples, positive orientation, character-three carry, and the
  common nonzero source scalar algebraically;

destroyed/not supplied:
  common two-origin selector, a fixed-source packet operation/intertwiner,
  recombined scalar current, Lambda-sheet basepoint, and physical QA/QAB
  action;

needed sidecar:
  a lawful origin-odd selector at q7/q11 and an origin-resolved,
  source-preserving map from the 449 semantic sheets to the charged
  coefficient fibre;

decisive next test:
  find a typed origin-asymmetric q7 selector justified by the same
  source/QA-to-QAB operation, or prove none exists.  Merely inserting
  another scalar phase cannot work because the common-selector q7 and
  q11 coefficients are already zero.
```

## 8. Reproduction

```bash
python3 04-computation/lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.py
python3 -O 04-computation/lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.py
```

Normal and optimized outputs byte-match.  The script contains no
executable Python `assert`.

```text
script 3ed346e0c631b34bd61f0c4d27d7f161e8d35b70decfb95f5207c5f57893d005
output 0faa0a24f6ba8b6c88b6bbfc4f225e38667097b1a937d977741453499884901c
```

Hashes use LF-normalized bytes.
