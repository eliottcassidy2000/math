# LRC(14) post-full-pool two-newcomer frontier — 2026-08-26

**Status:** **FINITE-EXACT SINGLE IMPLEMENTATION + REPLAY**, with the
scale-separated consequence **CONDITIONAL** on accepting that finite census.
THM-4191 supplies the independently audited lower bound used in the exact
transversal conclusion. This note is not a theorem promotion and does not
prove LRC(14).

## 1. Exact base transversal number

Let `P` be THM-4191's thirty-label pool and let `E_7(50)` be its lawful global
depth-seven repair hypergraph. THM-4191 proves, by two independent exact paths,
that `E_7(50)` has no ten-vertex transversal. The first probe checks against
all `821,737` edges that

```text
{8,80,85,88,95,145,168,193,240,252,286}
```

is an eleven-vertex transversal. In fact, deleting `8` leaves the sharp
depth-six hostile ten-body, which misses exactly `278` depth-seven edges, and
all `278` contain `8`. Therefore

```text
tau(E_7(50))=11.                                       (1)
```

The probe also extends THM-4191's closest ten-body by the same vertex `8`.

## 2. First exact two-newcomer staircase

For distinct positive `q_1,q_2` outside `P`, define

```text
E_d(q_1,q_2)={R in binom(P,d):
  mu(G_((P union {q_1,q_2})\R))>=4/63}.                (2)
```

The exact joint-wall program exhausts every one of
`binom(30,9)=14,307,150` pool nine-bodies. For `{q_1,q_2}={49,50}`:

```text
d=7: 145,462 edges, 0 equalities, 120 nine-covers;
d=8: 1,536,023 edges, 0 equalities, 0 nine-covers.      (3)
```

As in THM-4191, a transversal of size below nine would extend inside `P` to a
nine-transversal. Thus the second row gives `tau(E_8(49,50))>9`. Every
`B in binom(P,9)` consequently misses a lawful repair `R`, so safe-set
monotonicity gives

```text
mu(G_(B union {49,50}))>=4/63.                          (4)
```

The depth-seven row is sharp only for this certificate. Its first nine-cover
`{16,88,95,126,143,145,168,193,240}` has actual positive scaled body surplus
`63m-4D=4,150,835,328,292,128`; it is not an unsafe-body witness.

Uniform depth eight is false. At `{6,50}`, `E_8` has `13,497` edges and exactly
`472,050` nine-covers. The first is

```text
{8,63,80,84,85,88,120,143,145},                        (5)
```

whose actual scaled body surplus is again positive,
`601,044,495,065,784`. This is the hostile control separating repair-deck
failure from danger.

## 3. Strict limit deck and a cofinal ordered wedge

Write `D=18,241,159,416,480`. For `R in binom(P,8)`, let
`m_R/D=mu(G_(P\R))`, and define the strict double-limit deck

```text
L_8={R:81m_R>7D}.                                      (6)
```

The exact scan finds `5,267,460` edges, zero equalities, and no nine-cover.
Its minimum strict numerator is

```text
min_(R in L_8)(81m_R-7D)=944,928,                      (7)
```

at `R={15,16,20,40,170,190,193,240}`. Hence the limiting excess over `4/63`
is at least

```text
epsilon=4/8,513,189,685.                               (8)
```

THM-4170, equation (9), says that for a union `U` of `c` intervals,

```text
|mu(U intersect G_q)-(6/7)mu(U)|<=6c/(49q).            (9)
```

Apply `(9)` first to `U_R=G_(P\R)` and then to
`U_R intersect G_(q_1)`. Since the latter has at most `c_R+q_1` components,

```text
mu(U_R intersect G_(q_1) intersect G_(q_2))
 >=36M_R/49-36c_R/(343q_1)-6(c_R+q_1)/(49q_2).        (10)
```

The fixed pool has `7,133` wall cells, so `c_R<=7,133`. Splitting `(8)` evenly
between the two errors proves the explicit sufficient wedge

```text
q_1>=3,186,712,759,230,
q_2>=ceil(521,215,695*(7,133+q_1)),                     (11)
```

and the wedge with `q_1,q_2` swapped. Conditional on the single-path census of
`L_8`, every `B in binom(P,9)` obeys `(4)` throughout `(11)`. THM-4150 then
accepts every distinct positive odd-tail pair after doubling and arbitrary
positive common content.

The constants in `(11)` are intentionally crude. Their role is structural:
the separated two-newcomer tail is finite-deck controlled, leaving the
comparable-scale two-dimensional resonance strip as the sharp next lane.

## 4. Connection contract and validity firewall

```text
source:       lawful exact repair edges for two labelled newcomers
target:       every B in binom(P,9)
map:          R disjoint B gives
              B union {q_1,q_2} subset (P union {q_1,q_2})\R
preserved:    both newcomer labels, full deletion labels, exact 4/63
              threshold, Haar measure, content, and odd-tail transfer
destroyed:    safe-component address and witness phase after thresholding
sidecar:      complete R, joint-wall mass/equality ledger, and (m_R,c_R)
positive:     {49,50} closes at depth eight
hostile:      {49,50} at depth seven and {6,50} at depth eight
decisive test: exact nine-cover exhaustion on binom(P,9).             (12)
```

This does not control the comparable-scale strip, arbitrary two outsiders,
physical entry into the fixed pool, or LRC(14). The pair program is one exact
implementation; an independent joint-wall or fixed-prefix path is required
before theorem promotion.

## 5. Artifacts and replay

```text
04-computation/lrc14_full_pool_e7_exact_transversal_probe_20260826.cpp
sha256 ef8130b2cf24fe61ab51376c7ea452479f3afc45ed00dbec6a01e8788a5f5516

05-knowledge/results/lrc14_full_pool_e7_exact_transversal_probe_20260826.out
sha256 a5197e48c17d529752a622acd2e9d8b30c0a58ffd27bc92f23306cfbb67425c6

04-computation/lrc14_two_newcomer_deletion_staircase_probe_20260826.cpp
sha256 53ffe1642457241153a858c3da3c20cf31e8cc79ab6009f9321895c2f20b06c2

05-knowledge/results/lrc14_two_newcomer_deletion_staircase_probe_20260826.out
sha256 850f67053a7bfcbe73a0e117d5a2613836e32555e60cc8b49df66eb18ba36f5c

hash basis: raw LF bytes
```

Replay:

```bash
g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_full_pool_e7_exact_transversal_probe_20260826.cpp \
  -o /tmp/lrc-post-tau-e7
/tmp/lrc-post-tau-e7 | diff -u \
  05-knowledge/results/lrc14_full_pool_e7_exact_transversal_probe_20260826.out -

g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_two_newcomer_deletion_staircase_probe_20260826.cpp \
  -o /tmp/lrc-two-newcomers
/tmp/lrc-two-newcomers all | diff -u \
  05-knowledge/results/lrc14_two_newcomer_deletion_staircase_probe_20260826.out -
```
