---
id: THM-1192
title: THE ALL-PAIRS BEAT-PUNCTURE LAW — a six-comb cover of one seven-wall slow gap must cover every safe sum/difference beat of every pair by the complementary four or five combs; equivalently an exact phase-aware puncture deficit vanishes, and reduction modulo q/gcd(c,q) gives the explicit mixed-gcd count A(Q)=2⌈Q/14⌉−1
status: the all-range inclusion, exact count inequality, and A(Q) reduction are PROVED below; the reported finite eliminations are explicitly telemetry, not a proof that every harmonic-crowded packet violates the law
source: codex-2026-07-18 beat-puncture bridge audit
depends_on: [THM-1170 (critical beat points), THM-1175 (sum/difference beat pairing), THM-1176 (six-on-seven slow-gap cone), THM-1180 (witness-point branching)]
related: [THM-1181 (cyclic-gap events as faithful carriers), THM-1190 (global beat certificate), THM-1191 (thirteen-adic toothpick ladder)]
scripts: 04-computation/beat_puncture_mixed_gcd_codex_20260718.py -> 05-knowledge/results/beat_puncture_mixed_gcd_codex_20260718.out
---

# THM-1192 — the all-pairs beat-puncture law

Put

\[
 D_w=\{t\in\mathbb R/\mathbb Z:\|wt\|<1/14\}
\]

and let

\[
 G_k(a)=\left[\frac{14k+1}{14a},\frac{14k+13}{14a}\right]
 \qquad(0\le k<a)
\]

be a closed safe gap of the slow carrier $a$.  Let

\[
 a<b_1<\cdots<b_6
 \quad\text{and suppose}\quad
 G_k(a)\subseteq\bigcup_{i=1}^6D_{b_i}.                 \tag{1}
\]

The point of the theorem is that a beat does more than identify two equal
kill sets globally.  **Inside a safe carrier gap, every safe beat deletes its
defining pair from the list of possible coverers.**  A fast-fast beat therefore
turns the six-comb cover into an exact four-comb witness obligation.  This is
the discrete $j=4$ skeleton that was missing from the continuous flood-tail
view.

## The exact inclusion

Choose two distinct speeds $u,v$ from

\[
 S=\{a,b_1,\ldots,b_6\}
\]

and choose either

\[
 q=u+v
 \quad\text{or}\quad
 q=|u-v|>0.                                             \tag{2}
\]

Define the consecutive beat block

\[
 P_q(k)=
 \left\{p\in\mathbb Z:
 \left\lceil\frac{q(14k+1)}{14a}\right\rceil
 \le p\le
 \left\lfloor\frac{q(14k+13)}{14a}\right\rfloor
 \right\},                                             \tag{3}
\]

and, for any speed $c$, define its modular kill set

\[
 K_q(c)=\{p\in\mathbb Z:14\,d_q(cp)<q\},
 \qquad
 d_q(x)=\min(x\bmod q,(-x)\bmod q).                    \tag{4}
\]

Let

\[
 \mathcal C_{u,v}=\{b_1,\ldots,b_6\}\setminus\{u,v\}; \tag{5}
\]

if one of $u,v$ is $a$, this has five elements, and if both are fast it
has four.  Then every cover (1) satisfies the **all-pairs beat-puncture
inclusion**

\[
 \boxed{
 P_q(k)\setminus K_q(u)
 \ \subseteq\
 \bigcup_{c\in\mathcal C_{u,v}}K_q(c).
 }                                                       \tag{6}
\]

Equivalently, the exact phase-aware puncture deficit

\[
 \Delta^{q}_{u,v}(k)=
 \#\left(
 P_q(k)\setminus
 \left(K_q(u)\cup\bigcup_{c\in\mathcal C_{u,v}}K_q(c)\right)
 \right)                                                \tag{7}
\]

must vanish for every pair, both signs in (2), and every covered gap:

\[
 \boxed{\Delta^{q}_{u,v}(k)=0.}                         \tag{8}
\]

This is the sharpest exact statement supplied by the beat witnesses: no
union bound or discarded phase remains in (6)--(8).

### Proof

Take $p\in P_q(k)\setminus K_q(u)$ and put $t=p/q$.  By (3),
$t\in G_k(a)$, so $t\notin D_a$.  If $q=u+v$, then

\[
 vt=p-ut,
\]

while, after ordering $v>u$, if $q=v-u$, then

\[
 vt=p+ut.
\]

In either case

\[
 \|vt\|=\|ut\|.                                        \tag{9}
\]

The assumption $p\notin K_q(u)$ says $\|ut\|\ge1/14$; hence (9) says
that neither defining speed lies in its dangerous comb at $t$.  The cover
(1) must therefore be supplied by a faster speed other than $u,v$, which
is exactly membership in the right side of (6).  This proves (6), and (7)--
(8) are its cardinality form.  Notice that equality at $1/14$, including
an endpoint of $G_k(a)$, is handled correctly because every dangerous
window is strict.  ∎

## Exact and phase-free count inequalities

For a consecutive integer block $I$, write

\[
 C(c;q,I)=\#(I\cap K_q(c)),\qquad N=\#P_q(k).           \tag{10}
\]

Counting (6), before discarding overlaps, gives the exact necessary scalar
inequality

\[
 \boxed{
 N-C(u;q,P_q(k))
 \le
 \sum_{c\in\mathcal C_{u,v}}
 \#\big((P_q(k)\setminus K_q(u))\cap K_q(c)\big)
 \le
 \sum_{c\in\mathcal C_{u,v}}C(c;q,P_q(k)).
 }                                                       \tag{11}
\]

The first inequality still remembers the defining pair's phase; (8) is
stronger because it also remembers overlaps among complementary kill sets.

There is also a closed mixed-gcd relaxation.  Put

\[
 g_c=\gcd(c,q),\qquad Q_c=q/g_c.                        \tag{12}
\]

Writing $c=g_cc'$, multiplication by $c'$ permutes the residues modulo
$Q_c$, and

\[
 p\in K_q(c)
 \iff
 14\min(c'p\bmod Q_c,-c'p\bmod Q_c)<Q_c.               \tag{13}
\]

Thus the kill pattern has exact period $Q_c$.  In one full period its
cardinality is

\[
 \boxed{A(Q)=2\left\lceil\frac Q{14}\right\rceil-1.}   \tag{14}
\]

Indeed, the allowed residues are $0$, the positive residues

\[
 1,\ldots,\left\lceil Q/14\right\rceil-1,
\]

and their negatives.  They are disjoint because $1/14<1/2$, giving
$1+2(\lceil Q/14\rceil-1)$.  The ceiling, rather than a floor, is forced by
the strict inequality when $14\mid Q$.

For a block of $N$ consecutive integers, split off full periods.  The
remaining $R=N\bmod Q$ positions are distinct modulo $Q$, so at most
$\min(R,A(Q))$ are dangerous.  Hence

\[
 C(c;q,I)\le
 U(N,Q_c):=
 \left\lfloor\frac N{Q_c}\right\rfloor A(Q_c)
 +\min(N\bmod Q_c,A(Q_c)).                              \tag{15}
\]

Combining (11) with (15) on the defining and complementary speeds yields
the purely mixed-gcd necessary law

\[
 \boxed{
 N-U(N,Q_u)
 \le
 \sum_{c\in\mathcal C_{u,v}}U(N,Q_c).
 }                                                       \tag{16}
\]

For either choice in (2),
$\gcd(u,q)=\gcd(u,v)$.  The left side of (16) therefore measures how many
beat points survive the defining pair gcd, while the right side is the
maximum capacity of the complementary gcd periods.  Formula (16) is weaker
than the phase-aware law (8), but it is the promised direct bridge from beat
pairing to the mixed-gcd cone.

## What the computation says — and does not say

The exact companion script checked (8), using integers only.

* In the exhaustive cone
  $2\le a\le10$, $a<b_1<\cdots<b_6\le3a$,
  $a\sum_i1/b_i>1$, and
  $\gcd(a,b_1,\ldots,b_6)=1$, there were **69,258** packets.  For every
  packet and every $k$, some distinct-pair sum/difference beat violated
  (6).  Survivors: **0**.
* A seeded audit of **50,000** accepted harmonic-crowded mixed-gcd draws
  with $a\le60$ and $b_6\le6a-3$ likewise found **0** survivors.
* Carrier-fast sum beats alone are not enough.  In the smaller exhaustive
  cone $a\le7$, $b_6\le3a$, all **4,165** packets pass the phase-free count
  for some gap, **320** pass the exact sum count, and **107** pass the exact
  carrier-beat union.  The fast-fast obligations are doing real work.

A concrete carrier-only survivor is

\[
 (a;b_1,\ldots,b_6;k)=(5;6,7,9,11,12,15;2).             \tag{17}
\]

It passes every carrier-fast sum-beat inclusion, but the fast pair $6,7$
has $q=13$ and $p=6$.  The point $t=6/13$ lies in $G_2(5)$, and the
seven distances have least residues

\[
 4,3,3,2,1,6,1\pmod {13},                               \tag{18}
\]

all at least $1/13>1/14$.  Thus $t$ is an exact lonely witness.  This is
the smallest useful adversarial lesson from the computation: restricting
the defining pair to $(a,b_i)$ destroys precisely the fast-fast $j=4$
information.

The two zero-survivor statements are telemetry, not an all-$a$ proof.
What is missing is an analytic reason that some deficit (7) must be positive
for every packet in the harmonic-crowded cone.  The law samples finitely many
rational points; without such an argument it remains a necessary condition,
not the final contradiction.

## Structural reading: H-drift, toothpicks, Fano, and tournaments

The harmonic drift

\[
 H_j=a\sum_{i\le j}\frac1{b_i}
\]

records capacity but forgets every phase.  The natural missing state is not
one more scalar correction to $H_j$, but the vector of deficits

\[
 \big(\Delta^{u+v}_{u,v}(k),
       \Delta^{|u-v|}_{u,v}(k)\big)_{\{u,v\}\subset S}. \tag{19}
\]

The harmonic inequality selects the crowded cone; (19) tests whether its
claimed capacity can be transported to the exact pair-beat obligations.

The reduction $q\mapsto Q_c=q/\gcd(c,q)$ is also a precise form of the
ladder's toothpick self-similarity.  After quotienting by the gcd, every
complementary comb is the *same* strict $1/14$ window on a smaller cyclic
clock.  What changes from toothpick to toothpick is only the unit multiplier
and the truncated-block phase.  Formula (14) is the scale-invariant tooth
count; the residual phase in (8) is the part that cannot be replaced by a
density.

This also separates the result cleanly from the two newest beat/toothpick
threads.  THM-1190's global beat certificate counts blockers on the unit
numerators modulo $q$ and is still conjectural.  Here $P_q(k)$ is a local
consecutive block, numerators need not be units, and (6) is a proved pointwise
inclusion; the gain comes from deleting the defining pair, not from a
$\varphi(q)$ union bound.  Conversely, THM-1191 shows why a scalar overlap
floor can be fooled by a thirteen-adic scale stalk.  The exact deficit (7)
retains the unit multiplier and truncated-block phase after the gcd quotient;
replacing it by $A(Q)/Q$ alone would repeat precisely that loss of structure.

For the Fano/$\chi_7$ probe, (6) supplies an obligation on every edge of the
labelled seven-speed complete graph: a safe point created by an edge must be
covered from its four- or five-vertex complement.  This is compatible with
Fano line aggregation, but it exposes the missing transport identity.  The
periods here are

\[
 (u\pm v)/\gcd(c,u\pm v),                               \tag{20}
\]

whereas the existing Fano ledger uses line gcds such as
$\gcd(s_i,s_j,s_k)$.  Relating (20) to those line gcds without losing the
truncated phase is the still-open bridge.

Finally, the natural pair observable
$\Delta^q_{u,v}(k)$ is symmetric in $u,v$.  A switch or gauge that orients
it as a tournament edge is therefore arbitrary; the increasing-speed tie
Hamiltonian path contains no new information.  Vertices were challenged as
required: speeds, gaps, beat denominators, beat points, and proof obligations
were all considered.  The faithful object is the three-level incidence

\[
 \text{defining pair}\longrightarrow
 \text{safe beat point}\longrightarrow
 \text{complementary covering combs},                  \tag{21}
\]

not a binary tournament.  A runner tournament preserves only which pair was
chosen and destroys both the numerator $p$ and the complementary-cover
relation.  Consequently score histograms, cycles, SCCs, edge flips, and
Hamiltonian-path counts are not meaningful fingerprints for this script.

## Formalization target

The proof needs only integer modular arithmetic:

1. define $d_q(x)=\min(x\bmod q,(-x)\bmod q)$;
2. prove $d_q(up)=d_q(vp)$ from $q=u+v$ or $q=|u-v|$;
3. translate $p/q\in G_k(a)$ to the two integer inequalities in (3);
4. apply the assumed finite cover to obtain (6);
5. reduce by $g=\gcd(c,q)$, use that $c/g$ is a unit modulo $q/g$,
   and count the residues in (14).

No compactness, measure theory, floating arithmetic, or unproved critical-
point classification is required for the theorem itself.

## Exact replay

The checked-in output records these two restartable commands:

```text
python3 -u 04-computation/beat_puncture_mixed_gcd_codex_20260718.py --random 0
python3 -u 04-computation/beat_puncture_mixed_gcd_codex_20260718.py --skip-carrier --skip-exhaustive --random 50000
```

At close-out the replay artifacts had SHA-256 hashes

```text
script d45b9782e7efc0a0db8f388ccd9271b6c65237f5824077b9c426127251d20736
output e7f99276fd1628f331de9146960732c32ae82e2a95904c592014fec7ae519a10
```
