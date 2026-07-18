---
id: THM-1123
title: Bounded sharp four-comb component bank and the failure of scalar r=5 gates
status: PROVED BOUNDED FINITE-EXACT — all 45,238,050 rows in the first forty legal speeds over all 495 eight-speed cores satisfy the sharp component target; an independent Fraction arrangement replay verifies the complete core atlas, the reported extremal row, and exact counterrows to the scalar half-coverage and MST gates. This is not a uniform r=5 theorem
source: codex-2026-07-18-S73 r5 sharp-component continuation
depends_on: [THM-1097, THM-1126]
related: [THM-1101, THM-1111, THM-1127, MISTAKE-164]
script:
  - 04-computation/lrc14_r5_four_comb_component_bounded_exact_codex_S73.cpp
  - 04-computation/lrc14_r5_four_comb_bounded_fraction_replay_codex_S73.py
result:
  - 05-knowledge/results/lrc14_r5_four_comb_component_bounded_exact_codex_S73.out
  - 05-knowledge/results/lrc14_r5_four_comb_bounded_fraction_replay_codex_S73.out
---

# THM-1123 — bounded sharp four-comb component bank

For `D_k={t:||kt||<1/14}` and an eight-speed core
`P subset {1,...,12}`, put

```text
S(P)=[0,1] minus union_{p in P}D_p,
lo(P)=13 max(P)+1.
```

## 1. Bounded theorem

For every one of the `C(12,8)=495` cores and every ordered quadruple

```text
lo(P) <= k1 < k2 < k3 < k4 <= lo(P)+39,              (1)
```

the remainder `S(P) minus union_i D_ki` has a component of length `L`
satisfying

```text
7 k4 L > 1.                                            (2)
```

There are exactly

```text
495*C(40,4)=45,238,050                                 (3)
```

rows, and the exact referee finds zero failures.  This theorem has precisely
the bounded scope (1).  In particular, (3) supplies no reduction of larger
quadruples and uniform `r=5` remains open.

## 2. The 495-core atlas

The complete exact atlas of largest core-safe intervals gives

```text
ell(P) >= 1/70.                                        (4)
```

Equality is unique at

```text
P={1,2,6,7,8,9,10,11},                                (5)
```

whose longest intervals are

```text
[43/140,9/28], [19/28,97/140].                         (6)
```

The C++ certificate obtains the atlas by successive exact tooth subtraction.
The Python replay independently builds every breakpoint
`(14j+-1)/(14p)`, classifies the cells by rational midpoints, and reproduces
(4)--(6) over all 495 cores.

## 3. Exact bounded extremal

The global metric minimum in (1) is attained at

```text
P={1,2,4,5,7,9,11,12},
(k1,k2,k3,k4)=(158,160,162,164),
L=41/25920,
7 k4 L=11767/6480=1.815895... .                        (7)
```

The two longest surviving intervals are the reflected pair

```text
[673/2240,685/2268], [1583/2268,1567/2240].           (8)
```

All endpoint comparisons and the final strict test use signed `__int128`
cross-products; no floating-point value enters a certificate decision.  The
independent arrangement replay recovers (7)--(8) without calling the C++
subtraction logic.

Frozen SHA-256 values are

```text
C++ source     f1833ba222f3964141de97de1f1a27f4e444b0fce9393ca37db08b7abb40599f
C++ output     aac2577804199c4adb66cf627845e13cb8c2c7739036001a2bcf262a352c600f
Python source  df9c1de72126a96b72dd63772270e1223285d7829b2dfda6a02502c17e03789c
Python output  c61f72a0a7bbd9c0af342795dde7fc324a3320914c4b5d223a31c698e6802fe2
```

Fresh `-O0` and `-O2` C++ runs byte-match the frozen C++ output.  Normal and
optimized Python runs pass the assertions, and the normal run byte-matches the
frozen Fraction output.

## 4. Scalar sufficient gates are not the missing theorem

THM-1126's local half-coverage gate, its whole-core variant, and its
maximum-spanning-tree overlap criterion all remain useful sufficient tests.
They are not uniform explanations of (2), even inside the bounded bank.
Exact counterrows are:

```text
gate                         P                         (k1,k2,k3,k4)       exact slack
local half-coverage          {2,3,4,5,6,7,8,9}        (128,131,134,137)  -104610971/19393097472
maximum overlap tree         {2,3,4,5,6,7,8,9}        (121,124,127,130)  -466976387/54621386820
whole-core half-coverage     {1,2,3,6,7,9,10,11}      (153,156,159,162)  -4926398/150405255
```

The corresponding sharp metrics are respectively `1107/256`, `93/22`, and
`7/2`, so every row is safely on the positive side of (2).  Their negative
gate slacks prove that union mass and pairwise overlap do not encode the
load-bearing endpoint-spacing dispersion.

## 5. A better dispersion observable

Let the final safe gaps in a chosen core interval have lengths `l_j`, and put
`tau=1/(7k4)`.  Consecutive safe gaps are separated by a danger component of
length at least `tau`.  Therefore the short-shift autocovariogram satisfies

```text
C_tau = integral_0^tau |A intersect (A-h)| dh
      = sum_j phi_tau(l_j),

phi_tau(l)=l^2/2                    if l<=tau,
           tau*l-tau^2/2           if l>=tau.          (9)
```

If all `l_j<=tau`, then (9) gives `C_tau<=tau*|A|/2`.  Hence

```text
C_tau > tau*|A|/2                                  (10)
```

is another rigorous sufficient certificate for (2).  On the difficult row
(7), inside either longest core interval, the exact ratio is

```text
C_tau/(tau*|A|/2)=336218180863/317189087460
                 =1.0599929... >1.                    (11)
```

Thus the autocovariogram sees the hard row that the scalar half/MST gates
miss.  It is still only a sufficient observable, not a uniform tail theorem;
the promising next object is the full short-shift profile (or an equivalent
endpoint-owner word), rather than one scalar sample of it.

## 6. Carrier and tournament audit

Runner, comb, core-component, tooth, endpoint, residue, and proof-obligation
vertices were all challenged.  On endpoints, exact rational order gives a
transitive tournament: after coalescing exact ties its sorted word is the
unique Hamiltonian path, with no directed cycles and singleton SCCs.  Plain
order destroys lengths, owners, removal stages, and the dispersion in (9).
The faithful carrier for this bank is therefore the exact coordinate-plus-
owner endpoint word.  Comb vertices with intersection weights retain the MST
observable but lose higher overlap and endpoint order; scalar mass loses both.

## 7. Frontier

The new information is deliberately two-sided:

1. the sharp target survives a large exact bank with a margin of more than
   `0.815` at its hardest row; and
2. the natural scalar gates fail on exact rows where the target remains easy.

The feasible all-scale strategy is to combine THM-1126's analytic gates with
an endpoint-dispersion or toothpick self-similarity theorem (THM-1127) and use
a guarded exact bank only for the true finite complement.  Enlarging an
unguarded window would repeat MISTAKE-164.
