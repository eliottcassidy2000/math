---
id: THM-2167
title: "Rank-two adaptive small-base carry descent"
status: >
  PROVED from THM-2163 and THM-2164. Every zero-safe distinct thirteen-speed
  row has two relations of height at most 105 which remain independent modulo some
  prime q in {2,3,5,7,11,13}. In that base every digit layer lies in an
  affine codimension-two fibre of size q^11. The simultaneous carry has
  at most 2729^2 states; after sorting, the owner mask has fourteen
  possible suffix values. Repeated carry-owner states may be deleted while
  preserving positivity and both relations, but not distinctness or the
  lonely-runner target. Thus the theorem gives a finite adaptive automaton
  and an exact pumping boundary, not a proof of LRC(14).
source: codex-2026-07-24-rank-carry-synthesis
depends_on:
  - THM-2163
  - THM-2164
related:
  - THM-2144
  - THM-2161
  - THM-2163
  - THM-2164
script: 04-computation/lrc14_rank_two_adaptive_carry_referee_thm2167.py
output: 05-knowledge/results/lrc14_rank_two_adaptive_carry_referee_thm2167.out
script_sha256: 40ab89f2679fe374791eaee0aa4e89df1ec11f9f981881c8983485478bf47ae3
output_sha256: e066ab1d834140f3620bc7e2b93ea3915fa116f4d5f4af5984f81000ffa11147
hash_basis: working-tree bytes (LF)
---

# THM-2167 -- rank-two adaptive small-base carry descent

For a positive integer row `V=(V_1,...,V_13)`, put

```text
Lambda(V)={m in Z^13:m.V=0},
W_H(V)=span_Q{m in Lambda(V):||m||_infinity<=H}.       (1)
```

The theorem turns the analytic rank conclusion of THM-2164 into a bounded
radix carrier.  Its final pumping statement also identifies why bounded
relations alone do not make the speed search finite.

## 1. A bounded rank-two lattice has a small good prime

Let `r,s in Z^13` be linearly independent over `Q`, with

```text
||r||_infinity,||s||_infinity<=H.                     (2)
```

Some two-coordinate minor

```text
Delta_(i,j)=r_i s_j-r_j s_i                           (3)
```

is nonzero.  Every such minor obeys

```text
|Delta_(i,j)|<=2H^2.                                  (4)
```

If the reductions of `r,s` were dependent over `F_p` for every prime

```text
p in P={2,3,5,7,11,13},                               (5)
```

then every prime in `P` would divide (3).  Since the primes are distinct,

```text
30030=product_(p in P)p
```

would divide its nonzero value.  At `H=105`, however,

```text
2H^2=22050<30030.                                     (6)
```

This is impossible.  We have proved:

> **Small-prime lemma.** Two independent integer vectors of coefficient
> height at most `105` remain independent modulo at least one prime
> `q in {2,3,5,7,11,13}`.

The endpoint `13` cannot simply be removed from this abstract argument:
the two vectors `(105,0)` and `(0,22)` have determinant
`2310=2*3*5*7*11`, so they are dependent modulo every earlier prime and
independent modulo `13`.

## 2. Application to a zero-safe LRC row

Let `V_1,...,V_13` be pairwise distinct positive integers and suppose

```text
mu{t:||V_i t||>=1/14 for every i}=0.                  (7)
```

THM-2164 proves

```text
dim_Q W_105(V)>=2.                                    (8)
```

Choose two independent members

```text
r,s in Lambda(V),
||r||_infinity,||s||_infinity<=105.                   (9)
```

The small-prime lemma supplies a prime `q<=13` for which their reductions
are independent in `F_q^13`.  The base is chosen from the relation lattice
of the row; it is not a fixed universal residue bank of the kind ruled out
by THM-2161.

## 3. The simultaneous carry automaton

Write the common base-`q` expansion as in THM-2163:

```text
V=q^j Z_j+R_j,             0<=R_(j,i)<q^j,
D_j=Z_j mod q in {0,...,q-1}^13.                      (10)
```

For `a in {r,s}`, define

```text
kappa_j^a=(a.R_j)/q^j=-a.Z_j.                         (11)
```

THM-2163 gives, coordinatewise for both relations,

```text
kappa_0^a=0,
q kappa_(j+1)^a=kappa_j^a+a.D_j,
|kappa_j^a|<||a||_1,
kappa_J^a=0 once q^J>max_i V_i.                       (12)
```

At every layer the digit vector therefore satisfies the two affine
congruences

```text
r.D_j=-kappa_j^r mod q,
s.D_j=-kappa_j^s mod q.                               (13)
```

Because `r,s` are independent modulo `q`, the linear map

```text
F_q^13 -> F_q^2,          D |->(r.D,s.D)              (14)
```

has rank two.  Each fibre of (14), and hence each unrestricted digit fibre
in (13), has exactly

```text
q^(13-2)=q^11                                           (15)
```

elements.  Owner restrictions may delete elements of this fibre, but can
never enlarge it.

The strict bounds in (12) give at most

```text
(2||r||_1-1)(2||s||_1-1)
 <=2729^2
 =7,447,441                                             (16)
```

carry pairs.  Conversely, any finite digit word and carry-pair path obeying
(12), with both endpoint carries zero, reconstructs by THM-2163 a
nonnegative integer row satisfying both relations.  Thus the construction
is exact in both directions; it is not merely a necessary congruence test.

## 4. Owner masks and the pumping boundary

After relabelling, assume

```text
0<V_1<...<V_13.
```

Put

```text
O_j={i:V_i>=q^j}.                                     (17)
```

The masks are nested suffixes, so only fourteen values are possible,
including the empty suffix.  The combined carry-owner carrier has at most

```text
14*7,447,441=104,264,174                              (18)
```

states.  It records that

```text
support(D_j) subset O_j,       O_J=empty.              (19)
```

A zero digit at an active coordinate still does not mean termination:
higher nonzero digits may remain.

There is an exact pumping statement.  Suppose `0<=j<k<=J` and

```text
(kappa_j^r,kappa_j^s,O_j)
 =(kappa_k^r,kappa_k^s,O_k).                          (20)
```

Delete the digit block

```text
D_j,D_(j+1),...,D_(k-1)                               (21)
```

and concatenate the remaining lower and upper blocks.  Equality of the
carry pairs makes the two recurrences splice, so the reconstructed row
still satisfies

```text
r.V'=s.V'=0.                                          (22)
```

It is also positive.  A coordinate outside `O_j` had already terminated
below level `j` and is unchanged.  A coordinate in
`O_j=O_k` has a nonzero digit at some level at least `k`, which survives
after the upper block is shifted down.

What is *not* preserved is load-bearing:

1. two distinct coordinates can differ only inside the deleted block and
   become equal;
2. order and primitive normalization can change;
3. residues at later moduli and the phase geometry of (7) can change.

Therefore the finite carrier is pumpable as a relation language but is not
a target-preserving quotient for LRC.  A successful finite descent needs at
least a distinctness/order sidecar and a phase or modulus-owner certificate.
THM-2163's explicit target-mixed pair shows that even a full bounded relation
box plus a fixed residue prefix and scalar maximum does not supply that
missing coordinate.

This is the precise gain and limitation:

```text
zero-safe row
 -> two bounded independent relations
 -> an adaptive q<=13 with rank-two digit constraints
 -> finite carry-owner automaton
 -/-> finite LRC search without a target-preserving pump.             (23)
```

## 5. Exact hostile controls

The companion script checks (6), the endpoint example after the small-prime
lemma, an independent relation pair on the row `(1,...,13)`, every one of
the `2^13` binary digit vectors for the resulting rank-two map, the exact
`2^11` fibre count, both carry recurrences, and a block deletion with equal
carry-owner endpoints.  It also exhibits the promised loss: that deletion
preserves positivity and both relations but merges two coordinates.

The finite checks are hostile controls for the constants and the pumping
boundary.  The all-row theorem is the algebraic proof above.

QED.
