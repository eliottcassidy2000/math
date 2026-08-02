---
id: THM-3156
title: "Uniform two-star law on the cap-two boundary through width 100"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For H=(1,2,3,4,6,12), reflected cell 90, and every injection of its six
  labels into [D,2D] with 6<=D<=100, one fixed uniform law on the nine edges
  of the two-star has expected overlap strictly above the full singleton
  debt.  The exact global margin is
  2942257539652583/320491868062929255.  Hence some two-star edge closes this
  fixed cell by Bonferroni.  This is not an all-width cone or LRC(14).
audit: >
  The canonical companion imports and hash-pins only THM-2941's proved exact
  reflected-arc engine.  An independent implementation reconstructs the
  modular arcs, brute-forces the complete D=6,7 universes, and checks a
  separate full 16-state matching DP at the exact/rounded seam and selected
  tail widths.  Directed floor/ceiling signs, the factor 9, injective
  coverage, all 95 witness rows, both global minima, the digital-strip and
  periodic-ray formulas, and the Bonferroni consequence were independently
  checked.  Normal and optimized transcripts byte-match.
source: root/frontier-synthesis/uniform-two-star-beach-2026-08-02
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - THM-3150-cap-two-two-star-skeleton-and-common-dilation-tail
  - THM-3135-directed-cycle-weak-order-lane-cover-and-reflected-h-boundary
script: 04-computation/lrc_H_upper_median_uniform_two_star_boundary_thm3156.py
output: 05-knowledge/results/lrc_H_upper_median_uniform_two_star_boundary_thm3156.out
script_sha256: 0f8e87a4f00b1af78aeee2a8d2160e31830a9e464b01e6abef28a7be02395454
output_sha256: c477a863d8cd211984343e891a4924d1173040b40e81c7db2f75820a443b788d
independent_script: 04-computation/lrc_H_upper_median_uniform_two_star_boundary_independent_audit_thm3156.py
independent_output: 05-knowledge/results/lrc_H_upper_median_uniform_two_star_boundary_independent_audit_thm3156.out
independent_script_sha256: 066885479414ba9dce6bad0fb39d0567a3592502c2bf78f151bbad1558dd7a6a
independent_output_sha256: 148a323dad236f0d795122be8aa9286f188941dbb9e1e9e7c0411f649d74a7d0
hash_basis: LF-normalized bytes
---

# THM-3156 -- uniform two-star law on the cap-two boundary through width 100

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Work inside the sufficient reflected \(k=1\) family of THM-2941.  Fix

~~~text
H=(1,2,3,4,6,12),              L=168,              cell j=90. (1)
~~~

Labels \(0,1\), carrying body values \(1,2\), are pivots.  Let

~~~text
E_*={01,02,03,04,05,12,13,14,15}                            (2)
~~~

be their nine-edge two-star.

For label \(i\), level \(q\), and normalized cell coordinate
\(0\le u\le1\), let \(A_i(q)\) be the exact THM-2941 reflected arc set

\[
 A_i(q)=\left\{u:
 \left\|\frac{(qL-H_i)(j+u)}L\right\|<\frac1{14}\right\}.     \tag{3}
\]

Its exact singleton mass is

\[
 \mu(A_i(q))=\frac{qL}{7(qL-H_i)}
             =\frac17+\delta_i(q),\qquad
 \delta_i(q)=\frac{H_i}{7(qL-H_i)}.                          \tag{4}
\]

For every integer \(6\le D\le100\) and every injection

~~~text
q:{0,1,2,3,4,5} -> {D,D+1,...,2D},                          (5)
~~~

put

\[
 M_D(q)=\frac19\sum_{ij\in E_*}
       \mu(A_i(q_i)\cap A_j(q_j))
       -\sum_{i=0}^5\delta_i(q_i).                           \tag{6}
\]

Then

\[
 \boxed{M_D(q)\ge
 \frac{2942257539652583}{320491868062929255}>0.}             \tag{7}
\]

The lower bound is the exact global minimum over the entire declared finite
universe.  It is attained at

~~~text
D=6,                 q=(6,9,12,11,8,7).                     (8)
~~~

Equation (7) is a theorem about one fixed probability law—weight \(1/9\) on
every edge in (2)—not a separately optimized winning edge for each
assignment.

## 2. Bonferroni consequence

Write

\[
 \Delta(q)=\sum_i\delta_i(q_i),\qquad
 I_{ij}(q)=\mu(A_i(q_i)\cap A_j(q_j)).                        \tag{9}
\]

By (4),

\[
 \sum_i\mu(A_i(q_i))=\frac67+\Delta(q).                      \tag{10}
\]

Equation (7) says that the average of the nine \(I_{ij}\) is larger than
\(\Delta(q)\).  Therefore at least one edge \(ij\in E_*\) satisfies

\[
 I_{ij}(q)>\Delta(q).                                        \tag{11}
\]

Pairwise Bonferroni gives

\[
 \mu\left(\bigcup_i A_i(q_i)\right)
 \le\sum_i\mu(A_i(q_i))-I_{ij}(q)<\frac67.                   \tag{12}
\]

The cell \(90\) is body-safe for \(H\).  THM-2941's reflected completion
criterion says that one body-safe cell with union mass below \(6/7\)
contradicts aligned completion.  Thus every assignment in (5) closes this
fixed reflected cell.

This conclusion stays inside the inherited sufficient family.  It is not a
classification of physical LRC survivors and does not prove \(LRC(14)\).

## 3. Exact matching decomposition

Fix distinct pivot levels \(q_0,q_1\).  The objective (6) separates as

\[
 M_D(q)=B(q_0,q_1)+\sum_{\ell=2}^5 C_\ell(q_\ell),            \tag{13}
\]

where

\[
 B=\frac19I_{01}-\delta_0-\delta_1,\qquad
 C_\ell(x)=\frac19(I_{0\ell}+I_{1\ell})-\delta_\ell(x).       \tag{14}
\]

The four leaves must receive four distinct levels different from the
pivots.  This is an exact four-left-vertex assignment problem.

For any fixed leaf, retain its four cheapest available levels.  If an
optimal assignment used a level outside those four, at most three of the
four cheaper levels could be occupied by the other leaves.  Moving that leaf
to an unused cheaper level preserves injectivity and does not increase the
objective.  Hence an optimum exists inside the at-most \(4^4\) top-four
product.

The companion checks every one of the \(D(D+1)\) ordered pivot pairs and
exhausts this top-four product.  This covers every injection in (5), without
enumerating all \((D+1)P6\) six-tuples.  For \(6\le D\le12\), every cost and
minimum is computed in exact rational arithmetic.

## 4. Directed rounding and the sharp finite minimum

For \(13\le D\le100\), put \(S=10^{15}\) and define the integer score

\[
 C_D(q)=
 \sum_{ij\in E_*}\lfloor S I_{ij}(q)\rfloor
 -9\sum_i\lceil S\delta_i(q_i)\rceil.                        \tag{15}
\]

The directions are load-bearing:

\[
 C_D(q)\le
 S\sum_{ij\in E_*}I_{ij}(q)-9S\Delta(q)
 =9S M_D(q).                                                 \tag{16}
\]

The same exact matching decomposition minimizes \(C_D\).  All 88 rounded
minima are positive.  Their smallest value is

\[
 \min_{\substack{13\le D\le100\\q}} C_D(q)
 =85295963189446                                             \tag{17}
\]

at

~~~text
D=24,                 q=(24,36,48,27,32,45).                (18)
~~~

Consequently every rounded-range assignment satisfies

\[
 M_D(q)\ge
 \frac{85295963189446}{9\cdot10^{15}}
 =\frac{42647981594723}{4500000000000000}.                  \tag{19}
\]

Exact comparison gives

\[
 \frac{42647981594723}{4500000000000000}
 >
 \frac{2942257539652583}{320491868062929255}.               \tag{20}
\]

The exact minima for \(D=7,\ldots,12\) also exceed the \(D=6\) value.
Equations (8), (17), and (20) prove the sharp global statement (7).

## 5. Exact digital-strip form

The cell geometry also has a direct integer-tent form useful beyond this
finite theorem.  For labels \(e,f\), levels \(p\le q\le2p\), put

\[
 z=168p-e,\quad w=168q-f,\quad
 r=90e\bmod168,\quad s=90f\bmod168,\quad
 B=\frac{rw-sz}{168}.                                       \tag{21}
\]

The last quantity is an integer because \(z\equiv-e\), \(w\equiv-f\),
\(r\equiv90e\), and \(s\equiv90f\pmod {168}\).  For the six residues of
\(H\), direct endpoint checking shows that the \(p\) arcs for \((e,p)\) are
the full
intervals

\[
 \frac{[r+168k-12,\ r+168k+12]}{z},\qquad 0\le k<p,           \tag{22}
\]

and similarly for \((f,q)\).  For a fixed \(k\), only

\[
 \ell=\left\lfloor\frac{B+kw}{z}\right\rfloor
 \quad\hbox{or}\quad \ell+1                                 \tag{23}
\]

can contribute.  Indeed, the length of the possible \(\ell\)-strip is
\((z+w)/(7z)<1\), using \(q\le2p\), \(e,f\le12\), and \(z\ge156p\).
If

\[
 d_{k,\ell}=|B+kw-\ell z|,                                  \tag{24}
\]

then the exact pair mass is

\[
 I_{e,f}(p,q)=\frac1{zw}\sum_{k=0}^{p-1}\sum_{\ell\ {\rm in}\ (23)}
 \left(
 [12(z+w)-168d_{k,\ell}]_+
 -
 [12|z-w|-168d_{k,\ell}]_+
 \right),                                                   \tag{25}
\]

with \(0\le\ell<q\).

Indeed, (22) follows by solving the modular inequality (3).  The overlap of
two intervals with radii \(a,b\) and centre distance \(d\) is

\[
 [a+b-d]_+-[|a-b|-d]_+.                                     \tag{26}
\]

Clearing the common denominator \(zw\) gives (25).  Thus one pair is
computable by an \(O(\min(p,q))\) integer tent sum rather than by merging two
long rational interval lists.  Formula (25) is an exact closed form for the
sequence of located overlaps; it does not by itself supply an all-\(D\)
lower bound.

As one proved periodic ray, swap the order and take
\((e,p)=(6,3g)\), \((f,q)=(1,5g)\).  Then

\[
 z=504g-6,\qquad w=840g-1,\qquad B=-90g+3.                 \tag{27}
\]

Write \(k=3t+\rho\).  The class \(\rho=0\) uses \(\ell=5t\) and has
\(d_{k,\ell}=|-90g+3+27t|\); the class \(\rho=1\) contributes nothing; and
\(\rho=2\) uses \(\ell=5t+3\), with
\(d_{k,\ell}=|78g+19+27t|\), for exactly
\(n=\lfloor2g/3\rfloor\) contributing values.  Summing (25) gives

\[
 N_g=3276g^2-1848g+n(3024g-1008-2268n)
     =4284g^2-2520g+C_{g\bmod3},                            \tag{28}
\]

\[
 I_g=\frac{N_g}{(840g-1)(504g-6)},\qquad
 (C_0,C_1,C_2)=(0,-336,84).                                \tag{29}
\]

After clearing the positive denominators in \(I_{g+1}-I_g\), the three
residue transitions have numerators

\[
\begin{array}{ll}
0\mathbin{\to}1:&504(1787436g^2+2073474g+17),\\
1\mathbin{\to}2:&1008(1211238g^2+1314819g+139285),\\
2\mathbin{\to}0:&1008(999558g^2+964791g-34808).
\end{array}                                                  \tag{30}
\]

They are positive for \(g\ge1\), so this sequence is strictly increasing.
At the admissible floor \(g=2\),

\[
 I_g\ge I_2=\frac{2030}{280393}>\frac1{140}.                 \tag{31}
\]

This ray is not an obstruction to a tail theorem.

## 6. Independent replay

The canonical companion imports and hash-pins only

~~~text
04-computation/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py
SHA256 2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31
~~~

which is THM-2941's proved exact reflected-arc engine.  It does not import the
unpromoted one-cone proof candidate.

The independent audit reconstructs the arcs directly from (3), rather than
importing the main interval or matching routines.  It:

- reconstructs all 95 reported witness assignments and their exact margins;
- brute-forces every injection at \(D=6,7\);
- uses a separate full 16-state matching DP at
  \(D=6,7,12,13,37,62,100\);
- checks singleton masses, injectivity, rounding directions, the factor
  \(9S\), both global minima, and the semantic digest; and
- verifies (25) in 6,900 direct-versus-strip cases, the periodic ray through
  \(g=100\), its three symbolic difference polynomials, and the phase-blind
  zero-floor witness; and
- replays ordinary and optimized Python byte-for-byte.

Run

~~~text
python3 04-computation/lrc_H_upper_median_uniform_two_star_boundary_thm3156.py
python3 -O 04-computation/lrc_H_upper_median_uniform_two_star_boundary_thm3156.py
python3 04-computation/lrc_H_upper_median_uniform_two_star_boundary_independent_audit_thm3156.py
~~~

and compare with the declared outputs.

## 7. Coherent selector and boundary

For the finite zero-sum game with adversary \(q\) and prover edge
\(ij\in E_*\), use payoff

\[
 \Pi(ij,q)=I_{ij}(q)-\Delta(q).                              \tag{32}
\]

An existential pair certificate only proves
\(\max_{ij}\Pi(ij,q)>0\) separately for each \(q\).  Equation (7) is stronger:
the single mixed strategy uniform on \(E_*\) has

\[
 \mathbb E_{ij\sim{\rm Unif}(E_*)}\Pi(ij,q)
 =M_D(q)>0                                                   \tag{33}
\]

for every adversary in the finite universe.  This is a coherent mixed
section, not a moving fibrewise witness.

The graph in (2) is undirected and the payoff retains overlap magnitudes.
No intrinsic binary orientation appears, so forcing it into a tournament
would destroy the coordinate that makes (33) true.

Two boundaries remain exact.

1. The theorem stops at \(D=100\).  Formula (25) exposes the next target:
   prove a uniform lower bound for every high two-star channel, with the
   periodic ray (27) treated separately.
2. The standard phase-blind one-edge transport floors cannot prove the
   result.  To state this exactly, reduce \((p,q)\) to the unordered coprime
   channel \(P\le Q\), and put

   \[
   \phi(P,Q)=
   \begin{cases}
   0,&P+Q\le7,\\
   \max\{1/105,\ 1/49-1/(2PQ)\},&P+Q\ge8.
   \end{cases}                                               \tag{34}
   \]

   For one orientation define

   \[
   \eta_{a,p;b,q}=\frac{qa-pb}{pL-a},\qquad
   T_{a,p;b,q}=\left[
   \frac{\phi(P,Q)-2|\eta_{a,p;b,q}|}{(pL-a)/(pL)}
   \right]_+,                                                \tag{35}
   \]

   and use the larger of \(T_{a,p;b,q}\) and \(T_{b,q;a,p}\).
   At \(D=100\), the assignment

   \[
   (136,119,102,103,101,100)                                \tag{36}
   \]

   makes all nine resulting edge floors exactly zero, while its exact debt is

   \[
   \frac{656518999793763784}{2839155148366624716475}>0.      \tag{37}
   \]

   Thus this phase-blind route returns only the negative debt.
   Location-specific cancellation in (25), or an equivalent sidecar, is
   indispensable for repairing that route.

**End of proof.**
