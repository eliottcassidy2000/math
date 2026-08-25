---
id: THM-4026
title: "Sun 2-4-6-8 binomial counterexample"
status: >
  REFUTED CONJECTURE + PROVED COUNTEREXAMPLE + VERIFIED-EXACT +
  INDEPENDENTLY HOSTILE-AUDITED. The 2019 assertion that every positive
  integer is C(w,2)+C(x,4)+C(y,6)+C(z,8), with w,x,y,z>=2, is false:
  N=896315812331399 has no representation. Three exact routes in two
  implementation languages exhaust the height-forced universe; two disjoint
  residue banks end in independent integer-square tests and a redundant
  bounded modular cover ends with zero survivors. Both N-1 and N+1 are
  represented. No least-counterexample, priority, finiteness, density, or
  pointwise asymptotic claim follows.
source: user-supplied candidate + T. Adamczewski public gist + root / independent exact audits, 2026-08-24
audit: >
  PASS. The primary Python route enumerates 2,755,643,831 admissible canonical
  (x,y,z) triples, applies integer-formed binomial residue filters, and rejects
  all 263,434 in-range survivors by exact isqrt. A dependency-free Python
  reconstruction covers the 3,090,472,000 rectangular triples: primes through
  89 leave 324 terminal candidates, of which 31 fail the height test and 293
  fail exact isqrt; primes through 137 independently leave zero modular
  survivors. An independently written C++20/OpenMP verifier uses
  two disjoint prime banks, obtains 67,181 and 99,556 nonnegative square-test
  candidates with zero representations, and separately obtains a zero-survivor
  bounded modular cover. Small, high-scale, immediate-neighbor, normal/-O, and
  UBSan controls pass.
depends_on: []
related:
  - THM-4027-sun-two-four-six-eight-universal-modular-solubility
  - THM-4028-sun-two-four-six-eight-average-order-criticality
  - THM-2059-crt-fiber-product-phase-packet
  - THM-2043-period14-parity-hasse-jet-completeness
  - THM-2050-period14-top-germs-do-not-determine-global-loneliness
  - THM-2412-delta-exponential-and-central-newton-layer-split
  - HYP-1953-additive-basis-normal-form-spectrum
  - MISTAKE-209
  - MISTAKE-219
  - MISTAKE-222
  - MISTAKE-363
references:
  - "Z.-W. Sun, MathOverflow question 323541 (2019): https://mathoverflow.net/questions/323541"
  - "OEIS A306477: https://oeis.org/A306477"
  - "T. Adamczewski, public counterexample certificate gist (2026-08-24): https://gist.github.com/tadamcz/0c578c8b2b3fb92fe8584bc0725187e3"
source_note: 05-knowledge/reference/CORE-PAPERS-SUN-2468-2026-08-24.md
script: 04-computation/sun_2468_counterexample_thm4026.py
output: 05-knowledge/results/sun_2468_counterexample_thm4026.out
script_sha256: 5cb77eab853e40e7698c2a20fbdf4645e7162bc55beffbc11481502eee510a2d
output_sha256: 5261ac45a367c5f7f61a2bafd98e4f6eb04fcc97bff3ac304575ca3ede7d7aeb
semantic_sha256: 7fe93b58bc221c6d9859c411de7d5215743db51a211609100cc87015817a6b7e
independent_python_script: 04-computation/sun_two_four_six_eight_counterexample_independent_audit_thm4026.py
independent_python_output: 05-knowledge/results/sun_two_four_six_eight_counterexample_independent_audit_thm4026.out
independent_python_script_sha256: 39febe3affe818f03c1b3d83161aad688fc1e6da271495617757eb097dd022cc
independent_python_output_sha256: df11680c1a63266911a4f12bf1d3a71d27101dfe7dd6318ecffa34ddba05cf25
independent_cpp_script: 04-computation/sun_2468_counterexample_thm4026_independent.cpp
independent_cpp_output: 05-knowledge/results/sun_2468_counterexample_thm4026_independent.out
independent_cpp_script_sha256: 6fc67f593dc77fbef2f6bcee1c9fe917548c6144fa33790d459185f98c5dae25
independent_cpp_output_sha256: 990f637c740f9330b993c699cd0307c04912bf6280642f7b7bffb72d6df615db
hash_basis: raw LF bytes
---

# THM-4026 -- Sun's 2-4-6-8 conjecture is false

**REFUTED CONJECTURE + PROVED COUNTEREXAMPLE + VERIFIED-EXACT +
INDEPENDENTLY HOSTILE-AUDITED.** Put

\[
N=896315812331399.
\]

There do not exist `w,x,y,z in {2,3,...}` such that

\[
N={w\choose2}+{x\choose4}+{y\choose6}+{z\choose8}.       \tag{1}
\]

Thus the universal conjecture posed by Zhi-Wei Sun in 2019 is false. The
candidate was supplied to this session by the user and also appears in the
public Adamczewski certificate linked above. This theorem certifies the
integer claim; it makes no discovery-priority claim and does not assert that
`N` is the least counterexample.

## 1. Exact finite reduction

For `k=4,6,8`, all legal indices `2<=t<k` give the same zero value. They may
therefore be represented canonically by `t=k-1`. Exact monotonicity and binary
search give

```text
x=3..12112,       y=5..932,       z=7..281,             (2)
```

with support sizes `12,110`, `928`, and `275`. The bracketing values are

```text
C(12112,4)=896266258399820 < N < C(12113,4)=896562324551740,
C(932,6) =895693597430352 < N < C(933,6) =901490966993008,
C(281,8) =871896500955975 < N < C(282,8) =897353333100675. (3)
```

For a fixed canonical triple define

\[
R=N-{x\choose4}-{y\choose6}-{z\choose8}.
\]

Since `w>=2`, equation (1) holds exactly when

\[
R\ge1,\qquad 8R+1=q^2\quad(q\text{ odd},\ q\ge3),        \tag{4}
\]

in which case `w=(q+1)/2`. This is just

\[
8{w\choose2}+1=(2w-1)^2.                               \tag{5}
\]

No floating-point arithmetic or modular division is used. In particular,
every `C(t,k)` is formed over the integers before reduction, avoiding the
factorial-denominator error in MISTAKE-363.

## 2. Three exhaustive certificate routes

The primary script enumerates the `248,160` admissible `(y,z)` pairs and the
exact `2,755,643,831` triples satisfying
`C(x,4)<=N-C(y,6)-C(z,8)`. Necessary odd-prime residue masks leave `287,120` triples;
`263,434` also satisfy the nonnegative height test. Exact `isqrt` rejects all
of them. The survivor stream and exact-test stream have SHA-256

```text
7e9021758b26c124f24ba1d34e1e175f321b93326a1613507528a8d47b3ad66a
6339b50611656ce14a6ec748815863b7b5f52ab393483c7623f0f31811aa697c. (6)
```

The dependency-free Python audit reconstructs the full rectangular universe
of `3,090,472,000` canonical triples. Quadratic-residue filters through prime
`89` leave `324` candidates; `31` have residual below one and exact `isqrt`
rejects the remaining `293`. Continuing the same necessary filters through
`137` leaves zero candidates, giving a redundant pure bounded modular-cover
certificate.

The independent C++20 verifier works over the smaller admissible universe and
uses two disjoint prime banks. They leave respectively `67,181` and `99,556`
nonnegative candidates; independent exact square tests again return zero. Its
third bank, containing primes `11` through `137`, leaves no survivor. After
collapsing each higher-degree zero fibre, the OEIS control returns exactly

```text
4655=C(85,2)+C(14,4)+C(9,6)+C(7,8)
    =C(94,2)+C(7,4)+C(9,6)+C(11,8),                   (7)
```

The literal domain has extra index multiplicity when a zero summand is used;
this does not change existence. Two planted targets on either side of `N`
retain `67` and `143` canonical representations. Normal and optimized Python
agree byte-for-byte; the independent routes agree on the zero verdict and
support bounds, and a serial UBSan C++ replay reports no diagnostic.

The modular-cover route is **height-bounded**. Its zero mask is not a fixed
congruence obstruction to the unbounded problem; THM-4027 proves the opposite.

## 3. The hole is locally isolated

Both immediate neighbors are represented:

```text
N-1 = C(33663667,2)+C(9433,4)+C(16,6)+C(9,8),
N+1 = C(40920205,2)+C(6138,4)+C(22,6)+C(13,8).        (8)
```

This proves local isolation only at distance one; it does not claim that the
gap is globally minimal or statistically typical.

## 4. Two exact quadratic-form lenses

Put

```text
A=x^2-3x+1,                  B=2w-1.
```

The identities

\[
24{x\choose4}+1=A^2,
\qquad 8{w\choose2}+1=B^2                                  \tag{9}
\]

turn (1), for `S=C(y,6)+C(z,8)`, into

\[
A^2+3B^2=24(N-S)+4.                                      \tag{10}
\]

This is a `Q(sqrt(-3))`-norm compatibility fibre. It becomes equivalent to the
original low-pair equation only with the thin sidecars

```text
A=x^2-3x+1 for x>=3,              B odd and B>=3.          (11)
```

Forgetting (11) enlarges the search and is not a lawful converse.

There is also an infinite exact family of low-pair residue holes:

\[
384\left({w\choose2}+{x\choose4}+{1\over6}\right)
=\left((2x-3)^2-5\right)^2+48(2w-1)^2.                  \tag{12}
\]

If a prime `p` is `17` or `23 mod 30`, then `(5/p)=(-48/p)=-1`. Hence the
right side of (12) cannot vanish modulo `p`: vanishing would first force both
squares to vanish and then make `5` a square. Therefore

\[
{w\choose2}+{x\choose4}\not\equiv-1/6\pmod p.          \tag{13}
\]

The target occupies this low-pair hole at both `p=17` and `p=23`. Thus every
hypothetical representation would require the `(6,8)` pair to repair both
characters. This is useful structure, but it is not a full local obstruction.

## 5. Local rarity without local failure

The target is `20 mod 33`, the unique minimum-density class for the true
binomial periods, with exact probability `16/1089` instead of the uniform
`1/33`. Exact normalized target-density diagnostics are

```text
2^4:1       3^2:68/81       5^2:566/625       7^2:310/343
p=11:72/121 p=13:154/169    p=17:240/289      p=19:316/361
p=23:472/529.                                             (14)
```

The four small-prime entries are **FINITE-EXACT at the displayed levels**;
no all-level stabilization claim is made. For `p=11,13,17,19,23`, the
displayed first-level factors are stable because the exact critical-tuple
audit finds no fully critical target solution.

These factors diagnose a thin lower tail. They do not prove a zero count:
THM-4027 proves that every residue modulo every positive modulus is attained.
The missing coordinate is global height-compatible intersection, not residue
support.

## 6. Combinadic carry anatomy

The rank-eight Macaulay/combinatorial-number-system normal form is

\[
\begin{aligned}
N={}&{281\choose8}+{279\choose7}+{234\choose6}+{212\choose5}\\
   &+{188\choose4}+{136\choose3}+{43\choose2}+{15\choose1}. \tag{15}
\end{aligned}
\]

THM-2412 makes the adjacent-rank operation exact through
`Delta C(t,k)=C(t,k-1)`. This suggests a precise Pascal-carry reachability
program: can the unique decreasing normal form (15) be reached from at most
one legal atom in each even rank `8,6,4,2`, allowing the canonical zero
fibres? Value is preserved, but an unlabelled normal form forgets the carry
path and atom budget. A confluent and complete rewrite system has not been
proved, so this is not yet an equivalent reformulation. A finite-state
obstruction is OPEN; shared Pascal syntax alone supplies no bridge
(MISTAKE-222).

## 7. Typed connection ledger

| Source -> target | Map and preserved predicate | Lost information | Required sidecar / decisive test |
|---|---|---|---|
| bounded `(x,y,z)` box -> residue masks | reduce the square discriminant modulo odd primes; every exact square survives | sign, magnitude, and a common unbounded lift | exact index bounds plus terminal `isqrt`, or a zero mask on the bounded box |
| low/high summand split -> norm fibre | (9)--(10) preserve exact equality | the thin image of `A` and positivity/parity of `B` | constraints (11); exact intersection with the `(6,8)` pair |
| binomial coordinates -> falling-factorial/Hasse observers | finite differences preserve adjacent-rank structure and exact periods | reduction modulo a period forgets height | THM-4027 period classes plus a global height clock |
| normal form -> carry graph | Pascal moves preserve integer value | path, labels, and four-atom budget | rank-labelled reachability or an invariant separating (15) |
| local histograms -> exact support | CRT joins preserve compatible residue fibres | pointwise global equality | THM-2059-style fibre join plus the bounded height/square sidecar |

This resembles THM-2043/2050's complete-local-data/global-exit failure, but
there is no map from (1) to LRC(14); it is a methodological hostile, not an LRC
reduction. The natural carrier here is a symmetric bipartite compatibility
graph between the `(2,4)` and reflected `(6,8)` pair values. There is no
intrinsic orientation, so Tournament Analysis is not invoked.

## 8. Repaired frontiers

The universal statement is refuted. The following replacements remain
pairwise distinct and OPEN unless explicitly scoped otherwise:

1. determine the least counterexample, or certify a complete prefix below `N`;
2. decide whether the exceptional set is finite, infinite, or density zero;
3. prove an almost-all or cofinite representation theorem with a uniform
   lower-tail estimate, rather than an average count;
4. search low-density CRT progressions while retaining the Archimedean height
   box, and decide whether they contain infinitely many holes;
5. classify the constrained norm intersections in (11), using the character
holes (13) as a sidecar rather than as a false full obstruction;
6. decide the rank-labelled Pascal-carry reachability problem attached to
   (15); and
7. minimize and formally verify the bounded prime cover. A fixed prime bank
   without its height sidecar cannot work, by THM-4027.

THM-4028 proves only the summatory mean-order mechanism. It explains why the
tuple count is barely supercritical, but it supplies no pointwise positivity.

## 9. Reproduction

```bash
python3 04-computation/sun_2468_counterexample_thm4026.py
python3 -O 04-computation/sun_2468_counterexample_thm4026.py
python3 04-computation/sun_two_four_six_eight_counterexample_independent_audit_thm4026.py
python3 -O 04-computation/sun_two_four_six_eight_counterexample_independent_audit_thm4026.py
g++ -std=c++20 -O3 -fopenmp -Wall -Wextra -Wconversion \
  04-computation/sun_2468_counterexample_thm4026_independent.cpp \
  -o /tmp/sun_2468_counterexample_thm4026_independent
OMP_NUM_THREADS=4 /tmp/sun_2468_counterexample_thm4026_independent \
  896315812331399 1
```

The frozen outputs and all raw-LF hashes are recorded in the frontmatter.
