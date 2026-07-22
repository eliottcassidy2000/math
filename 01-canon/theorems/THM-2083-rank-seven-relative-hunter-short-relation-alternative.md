---
id: THM-2083
title: "Every rank-seven terminal guard obstruction has a uniformly bounded three-term relation"
status: >
  PROVED by character convergence and THM-2081. There is an absolute H such
  that if a seven-speed safe set is contained in one guard, then some triple
  (h,q_i,q_j) has a nonzero integer relation of coefficient height at most H.
  Equivalently, packets whose every guard/two-speed relation height tends to
  infinity have relative-Hunter margin tending to 30/343 and cannot obstruct.
  This compactness proof is uniform but qualitative. THM-2085 subsequently
  makes it effective with H=57; classification of those finite coefficient
  templates remains open.
source: codex-2026-07-22-LRC14-relative-Hunter-compactness
depends_on:
  - THM-2081
related:
  - THM-685
  - THM-2051
  - THM-2054
  - THM-2082
  - THM-2085
---

# THM-2083 -- the rank-seven relative-Hunter short-relation alternative

For positive integers `u_1,...,u_k`, define their integer relation height by

```text
lambda(u_1,...,u_k)
 =min{||m||_infinity:m in Z^k minus {0}, sum m_i u_i=0}. (1)
```

The minimum is finite for `k>=2`. For a guard `h` and a labelled seven-set
`Q={q_1,...,q_7}`, put

```text
Lambda(h,Q)=min_(i<j) lambda(h,q_i,q_j).                (2)
```

## 1. Statement

There is an absolute positive integer `H_7` such that

```text
G_Q subset E_h     implies     Lambda(h,Q)<=H_7.        (3)
```

No parity, divisor-completeness, hereditary-primitivity, or height bound is
needed for (3). In the dyadic-tower application `h` is odd and those extra
properties remain available after the reduction.

Equivalently, if `(h_n,Q_n)` is any sequence with

```text
Lambda(h_n,Q_n)->infinity,                              (4)
```

then the THM-2081 relative-Hunter quantities satisfy

```text
sum_i I_i(h_n,Q_n) ->2/7,
Delta_(h_n)(Q_n)   ->0,
w_ij(h_n,Q_n)      ->5/343       for every i<j,
tau_(h_n)(Q_n)     ->30/343,                              (5)
```

and hence

```text
tau_(h_n)(Q_n)-Delta_(h_n)(Q_n)->30/343>0.              (6)
```

Thus every sufficiently far member of the sequence has positive-measure
safe time outside the guard.

## 2. Character convergence lemma

Let `v_n=(v_(n,1),...,v_(n,k))` be integer vectors and let `mu_n` be the
pushforward of Haar measure on `R/Z` under

```text
t ->(v_(n,1)t,...,v_(n,k)t) in (R/Z)^k.                (7)
```

If `lambda(v_n)->infinity`, then `mu_n` converges weakly to Haar measure on
the full `k`-torus.

Indeed, for a character `m in Z^k`,

```text
integral exp(2pi i m.(v_n t)) dt
 =1 if m.v_n=0,
 =0 otherwise.                                         (8)
```

For every fixed nonzero `m`, the first case is eventually excluded by the
relation-height hypothesis. Hence every nonconstant Fourier coefficient of
`mu_n` eventually vanishes. Trigonometric polynomials are uniformly dense in
the continuous functions on the torus, proving weak convergence. For boxes
whose boundary has Haar measure zero, their indicator integrals converge as
well by the standard inner/outer continuous approximation. QED.

This is the varying-direction form of the Kronecker/relation-lattice split in
THM-685. The relation-height sidecar is essential: without it, the image of
(7) can remain trapped in a proper subtorus.

## 3. Proof of the asymptotic relative-Hunter margin

Assume (4). For each `i<j`, relation vectors supported on only two of the
three coordinates are allowed in (1). Therefore

```text
lambda(h_n,q_(n,i))->infinity,
lambda(h_n,q_(n,i),q_(n,j))->infinity.                 (9)
```

Apply the character lemma with `k=2` to the guard/danger box. Its two side
measures are `2/7` and `1/7`, so

```text
I_i=measure(E_(h_n) intersect D_(q_(n,i)))->2/49.       (10)
```

Apply it with `k=3` to

```text
E_(h_n)^c times D_(q_(n,i)) times D_(q_(n,j)).          (11)
```

The side measures are `5/7,1/7,1/7`, hence

```text
w_ij->5/343.                                           (12)
```

There are seven vertices and finitely many labelled spanning trees. Every
tree has six edges. Since all 21 edge weights converge to the same number,
the maximum of their tree sums converges to

```text
6*(5/343)=30/343.                                      (13)
```

Equations (10), (12), and the definition
`Delta=2/7-sum I_i` prove (5)--(6). THM-2081 then supplies safe mass outside
the guard for all sufficiently large `n`. QED.

## 4. Uniform bounded relation

If no `H_7` in (3) existed, then for every positive integer `n` one could
choose a contained packet with `Lambda(h_n,Q_n)>n`. This sequence would
satisfy (4), while containment and THM-2081 would force

```text
tau_(h_n)(Q_n)<=Delta_(h_n)(Q_n)                       (14)
```

for every `n`. Equation (6) contradicts (14). Therefore some absolute `H_7`
exists, proving (3). QED.

The proof is qualitative but effective in principle: replace the continuous
boxes by one-sided trigonometric polynomials and retain all characters through
a chosen degree. THM-2085 carries this out with signed Selberg box minorants
and proves the explicit value `H_7=57` is sufficient.

## 5. Frontier effect

Every depth-four THM-2073 terminal obstruction now lies on one of finitely
many relation templates

```text
a h+b q_i+c q_j=0,
max(|a|,|b|,|c|)<=H_7.                                 (15)
```

This is a support-at-most-three relation involving the **guard**, not merely
a relation internal to the terminal core. It is exactly the coordinate that
the global THM-2051 relation alternative does not retain after forgetting the
last safe-child guard.

The result does not make the terminal speeds finite. Each template (15) is an
unbounded rational plane and may still carry divisor-complete hereditary
families. THM-2085 makes `H_7=57` explicit. The remaining decisive task is a
CRT/endpoint classification on those finitely many planes, or a direct proof
that each plane violates the relative-tree inequality.

## 6. Assumption challenge and Tournament Analysis

The challenged assumption is that genericity should mean large or separated
speeds. The correct generic coordinate is absence of bounded characters on
each guard/two-speed torus. Widely separated speeds can still have a tiny
relation; comparable speeds can be relation-free to any fixed height.

No tournament quotient carries the proof. In the relation-free limit all
restricted edge weights become equal, so every tournament orientation is a
tie-breaking artifact while every spanning tree has the same limiting weight.
The faithful carrier is the character-relation lattice together with the
weighted graphic matroid of THM-2081. QED.
