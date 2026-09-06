---
id: THM-4442
title: "LRC14 bounded ten-body parity-free scale-three completion"
status: >
  PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED. For every
  ten-subset C of {1,...,13} and every three distinct positive integers not
  divisible by three, G_(3C union T) is nonempty. The proof is a parity-free
  longest-component lemma plus an exhaustive 174,045-row exact residual.
  This is a bounded local completion theorem, not arbitrary chart entry;
  LRC(14) remains open.
source: root component-address continuation + independent referee, 2026-09-06
depends_on: []
related:
  - THM-737-pack-clock-sampling-measure-dispatch
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits
  - THM-4446-lrc14-primitive-ten-pack-descent-and-dilation-rays
geometry_script: 04-computation/lrc14_bounded_ten_body_geometry_thm4442.py
geometry_script_sha256: b03fd618d24faaf7edc1ec6e901cabcd4739886be6f472c61ba4a2f2baf7d0a5
primary_script: 04-computation/lrc14_bounded_ten_body_parity_free_thm4442.py
primary_output: 05-knowledge/results/lrc14_bounded_ten_body_parity_free_thm4442.out
primary_script_sha256: 25f6d300e66fe114f6c6a23b0b8c2b8eb4274e475eb112c1e3d5600cdae8761f
primary_output_sha256: 50e52182b2e2d2f42e59b27956918a417f7faeaefb60ff91bd64bc075ed37ee7
independent_script: 04-computation/lrc14_bounded_ten_body_parity_free_thm4442_independent.py
independent_output: 05-knowledge/results/lrc14_bounded_ten_body_parity_free_thm4442_independent.out
independent_script_sha256: 0f3ee18bbdc7e8aa656afa2e69d97f018489f0b33ce4623014f26437ac4a815b
independent_output_sha256: 7f8f1dc3e1cb427d9657c046f41ea5be6aa98e08ec17fc6d67713316a82c76f7
report: 05-knowledge/results/lrc14_bounded_ten_body_parity_free_thm4442.md
audit: 05-knowledge/results/lrc14_bounded_ten_body_parity_free_thm4442_independent_audit.md
hash_basis: raw LF repository bytes
---

# THM-4442 -- LRC14 bounded ten-body parity-free scale-three completion

**PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED.** This closes a
complete bounded clock-three branch. It does not supply entry into that branch,
and `LRC(14)` remains **OPEN**.

## 1. Statement

For a finite positive integer set \(A\), put
\[
 G_A=\{x\in\mathbb R/\mathbb Z:\|ax\|\ge 1/14\text{ for every }a\in A\}.
\]
If \(C\) is any ten-element subset of \(\{1,\ldots,13\}\), and
\[
 T=\{a,b,c\},\quad 1\le a<b<c,\quad 3\nmid abc,
\]
then
\[
                    G_{\,3C\cup T}\ne\varnothing.                 \tag{1}
\]
There is no parity or tail-primitivity condition. The thirteen speeds are
distinct and the full row is primitive: the two parts occupy different
classes modulo three, and the common gcd divides three but not any tail.

## 2. Parity-free longest-component lemma

**Lemma.** Suppose \(G_C\) has a closed component \(J\) of positive length
\(L\). For \(T\) as above, \(7cL\ge3\) implies
\(G_{3C\cup T}\ne\varnothing\).

**Proof.** For \(y\in J\), all three lifts \(x_j=(y+j)/3\), \(j=0,1,2\),
are safe for \(3C\). For a tail \(t\), define
\[
 K_t(y)=\{j:\|t(y+j)/3\|<1/14\}.
\]
Because \(3\nmid t\), its three phases are separated by \(1/3\), up to
permutation. The strict danger arc has length \(1/7<1/3\), so
\(|K_t(y)|\le1\).

Assume the full good set is empty. At every \(y\), the three sets
\(K_a(y),K_b(y),K_c(y)\) cover the sheet labels. Hence all are singletons.
The relatively open sets
\[
 U_j=\{y\in J:j\in K_c(y)\}
\]
are pairwise disjoint and cover connected \(J\), so one fixed \(U_j\) equals
\(J\). The closed lift interval \(I_j=(J+j)/3\), of length \(L/3\), lies in
the open danger set for speed \(c\). Its disjoint teeth have length
\(1/(7c)\). Connectedness puts \(I_j\) in one tooth, and closed containment
in an open interval forces \(L/3<1/(7c)\), a contradiction. \(\square\)

The fixed-sheet connectedness step is essential: sheet counting alone does
not produce the strict width inequality, and strictness makes equality safe.

## 3. Exact finite reduction

For all \(286=\binom{13}{10}\) bodies, exact endpoint geometry constructs
every positive component of \(G_C\). Complements of merged open danger teeth
agree with exact endpoint-arrangement cells. The global longest-component
floor is
\[
 \min_{C\in\binom{[13]}{10}}\max_{J\subset G_C}|J|={9\over1232},    \tag{2}
\]
attained at \(C=(1,2,3,5,6,7,8,9,11,13)\). Therefore an unresolved largest
tail satisfies
\[
 c<{3\over7L(C)}\le {176\over3},
\]
so \(c\le58\). The executable derives each body-specific ceiling from
\(L(C)\) and asserts both (2) and the global cutoff.

The exact residual has \(174,045\) rows. Every row has positive safe measure,
with no counterexample and no isolated-only completion. Two physical
representations agree on each row: quotient-component cells carrying
three-bit lift masks, and the direct merged danger union for \(3C\cup T\).

```text
residual rows                         174,045
positive-measure rows                 174,045
isolated-only rows                          0
literal positive components         5,171,992
quotient mask pieces                3,643,548
quotient pieces with a safe sheet   3,559,274
```

All eight parity masks occur. Overlapping occurrences of the signed circuits
\((1,1,1),(1,1,2),(1,2,2)\) number \(7,092,14,249,11,765\), respectively.
The smallest positive completed-row measure is
\[
 {10517879\over643242600}
\]
at
\[
 C=(1,2,3,5,6,7,8,9,11,13),\qquad T=(8,34,50),
\]
with twenty positive components. The tail gcd is two, explicitly controlling
against a hidden primitive-tail assumption.

## 4. Independent audit

The clean-room implementation imports no primary code. It constructs every
literal row on denominator \(14\operatorname{lcm}(3C\cup T)\), classifies
cells by integer modular arithmetic, checks every weak endpoint, and
reverse-lifts a physical midpoint to a quotient body point and labelled
sheet. It independently reproduces (2), the cutoff, all \(174,045\) rows,
the \(5,171,992\) components, every parity-mask count, and the exact minimum.
The [audit note](../../05-knowledge/results/lrc14_bounded_ten_body_parity_free_thm4442_independent_audit.md)
records the engine separation.

## 5. Scope and proof-map effect

This strengthens [THM-737](THM-737-pack-clock-sampling-measure-dispatch.md)
at its zero-margin case: clock three, ten body speeds, and three detuned
speeds. Component width plus a finite exact residue settles all 286 bounded
bodies at the clock where the uniform counting bound has no margin.

It does not prove that an arbitrary unresolved LRC(14) row enters a form
\(3C\cup T\) with \(C\subseteq[13]\), nor resolve the owner/arrival and
crossing-height obligations retained by
[THM-3818](THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md).
This is a local consumer after bounded entry, not entry itself.

## 6. Reproduction

```powershell
python -B 04-computation/lrc14_bounded_ten_body_parity_free_thm4442.py
python -B -O 04-computation/lrc14_bounded_ten_body_parity_free_thm4442.py
python -B 04-computation/lrc14_bounded_ten_body_parity_free_thm4442_independent.py
python -B -O 04-computation/lrc14_bounded_ten_body_parity_free_thm4442_independent.py
```

All theorem decisions use integer or `Fraction` arithmetic. Canonical raw-LF
hashes are recorded in the frontmatter.
