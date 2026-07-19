---
id: THM-1211
title: Exact four-torsion centre-hit congruence and refutation of the proposed uniform standoff
status: PROVED (elementary algebra; exact referee replay)
source: codex-2026-07-18-S79
depends_on: [THM-1181]
related: [THM-1174, THM-1177, THM-1210, THM-1214, MISTAKE-180]
script: 04-computation/lrc14_center_hit_congruence_referee_codex_S79.py
output: 05-knowledge/results/lrc14_center_hit_congruence_referee_codex_S79.out
script_sha256: 743c485510642732f10302693b49da289e03dee571b2d94a7778e0dcef95240f
output_sha256: dbdd30ac4508a1554f958591d0abe67eee4bef7a9c828df8fa22e40fa881af8a
---

# THM-1211 — the centre-hit condition is congruence, not proportionality

## Statement

Let `d=(d1,d2,d3)` be a nonzero integer direction, put

```text
g = gcd(d1,d2,d3),       e = d/g,
gamma_d(u) = ({-d1 u}, {-d2 u}, {-d3 u}) in T^3,
```

and let `r` be any permutation of `(1,2,3)`.  Then

```text
gamma_d(u) = r/4 for some u
    iff
e = r (mod 4) or e = -r (mod 4), coordinatewise.              (1)
```

More precisely, the two congruence classes have canonical witnesses:

```text
e =  r (mod 4)  =>  u = 3/(4g),
e = -r (mod 4)  =>  u = 1/(4g).
```

Consequently there are infinitely many nonproportional directions hitting
each centre.  For the ordered centre `(1/4,1/2,3/4)`, for example,

```text
d_m=(1,2,4m+3),       u=3/4,       m>=1                         (2)
```

gives `gamma_d_m(u)=(1/4,1/2,3/4)` exactly.  The first row is the small
counterexample `d=(1,2,7)`.

Thus the superseded six-box draft's assertion “a centre is hit iff the
direction is proportional to `(1,2,3)`” and its proposed positive standoff for
every nonproportional integer direction are false.  This does **not** refute
the possible measure ceiling `bad(d)<=2/21`: an isolated centre hit need not
give positive, let alone maximal, sojourn time.  Uniform clustered `r=5` is
instead closed independently by THM-1214's carrier-owner window argument.

## Proof

Suppose first that `gamma_d(u)=r/4`.  Then `4 d_i u` is an integer for every
coordinate.  Since the primitive coordinates `e_i` have gcd one, Bezout
coefficients `b_i` exist with `sum b_i e_i=1`.  Therefore

```text
a := 4 g u = sum_i b_i (4 d_i u)
```

is an integer.  Multiplying the centre congruences by four gives

```text
-a e_i = r_i (mod 4)                                                (3)
```

for all three coordinates.  Two entries of `r` are odd, so `a` is odd.
Hence `a=1` or `a=3` modulo four.  Equation (3) gives respectively
`e=-r (mod 4)` or `e=r (mod 4)`.  This proves necessity.

Conversely, if `e=r (mod 4)`, take `u=3/(4g)`.  Then
`-d_i u=-3e_i/4=r_i/4 (mod 1)`, because `-3=1 (mod 4)`.  If
`e=-r (mod 4)`, take `u=1/(4g)` and obtain
`-d_i u=-e_i/4=r_i/4 (mod 1)`.  This proves sufficiency and (1).
Equation (2) is the first congruence class with primitive residue vector
`(1,2,3)` modulo four, and is visibly nonproportional for every `m>=1`.

## Exact replay

The referee checks all `249,984` labelled rows consisting of a sorted triple
`1<=d1<d2<d3<=64` and one of the six labelled centres.  It compares the
common-gauge condition (3), the signed residue test (1), and the canonical
rational witness.  All agree.  There are `9,362` exact centre hits, of which
`9,320` lie outside the proportional `(1,2,3)` ray.  It also checks the first
64 members of (2).  The computation is a replay, not a premise of the proof.

## Structural viewpoint and tournament challenge

The faithful quotient here is not the ordering tournament of the three
frequencies.  It is the primitive residue vector in `(Z/4Z)^3`, modulo the
two units `+1` and `-1`.  What must be preserved is one *common gauge* `a`
in (3).  Runner vertices, gap vertices, centre-slot vertices, pairwise residue
comparisons, and proof-obligation vertices were all considered.  Any binary
tournament on them forgets whether the same gauge works simultaneously in
all three coordinates.  Consequently no natural tie Hamiltonian path or
tournament fingerprint controls centre membership; imposing one would add
telemetry while destroying the predicate.  The useful “switch” is instead
the unit action `e -> +/-e` on the four-torsion residue vector.

## What this sharpens

The six-centre picture remains a useful geometric visualization of the bad
region, but “hits a centre” is far too weak to characterize the maximising
flow.  A viable maximizer argument must distinguish **contact order or
sojourn length** from incidence: proportional directions may track a balanced
chamber for an interval, whereas the modular families above can pass through
its centre at a single phase.  This is the missing structural coordinate in
the centre-only view.
