---
id: HYP-9080
title: "Tournament deletion-slack local unavoidability"
status: >
  OPEN CONJECTURE + PROVED REDUCTION IDENTITY + FINITE-EXACT EVIDENCE, NOT
  PROVED. For a tournament C and vertex x, put S=H-disc and
  chi_x(C)=2(S(C)-S(C-x)). THM-3729 and THM-4094 identify this exactly with
  twice the Hamiltonian insertion gain minus the actual rooted odd-Pfaffian
  boundary energy plus the deletion discriminant. The conjecture says every
  tournament of order at least two has some x with chi_x>=0. It would prove
  H>=disc by deletion induction. Pointwise monotonicity is false: a source
  over C3 has one chi=-2. Exhaustive labelled orders 2..6 and all 456 order-7
  isomorphism representatives have no all-negative row; the stronger summed
  charge is also nonnegative on that universe. Neither statement is proved
  in general.
source: codex-snark-apex-260822870-20260825
depends_on:
  - THM-3729-rooted-pfaffian-response-and-sign-root-deletion-average
  - THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness
related:
  - THM-1950-h-ge-disc-reduced-to-strongly-connected
  - THM-2290-context-selected-colored-pair-kernel-is-hafnian-complete
  - THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary
  - THM-4104-selected-order-eleven-strong-ear-solid-interval
  - THM-4111-uniform-ear-average-and-recursive-selected-bank-growth
  - THM-4113-maximal-noncrossing-half-kempe-atlas
script: 04-computation/tournament_deletion_slack_unavoidability_hyp9080.py
output: 05-knowledge/results/tournament_deletion_slack_unavoidability_hyp9080.out
script_sha256: 6ff91dacc06b63c7e7fcd34163ca647e35a814f683eb0aa9fdb5a8ebd54000e6
output_sha256: 0c41b86d3a445df8fcfcf2e2a49e7363f9410766de7c495a6ac02ac83610549b
hash_basis: raw working-tree bytes; normalized-LF hash for the inherited n=7 representative bank
audit: >
  FINITE-EXACT PASS. The script exhausts all 33,866 labelled tournaments of
  orders 2 through 6 and the inherited 456-representative order-7 bank
  (353 strong classes). It computes H for every deletion by subset DP and
  separately evaluates disc and the rooted response with exact Bareiss
  determinants; these are formula cross-checks sharing one determinant
  kernel, not independent implementations. It also reconstructs full
  deletion-incidence fibers through order five. It checks
  410,408 rooted-response gates and 5,404 incidence gates. Normal and
  optimized runs reproduce the frozen transcript. Local unavoidability,
  its averaged strengthening, and H>=disc remain OPEN.
---

# HYP-9080 -- tournament deletion-slack local unavoidability

**OPEN + FINITE-EXACT.** The apex-cubic proof has two logically separate
halves: local reducibility and global unavoidability. The same shape gives a
precise proposal for the repo's open tournament inequality `H>=disc`, but no
planar theorem transfers.

## 1. Exact deletion charge

Let `C` be a finite tournament of order at least two, let `x` be a vertex,
and put

```text
T=C-x,                    S(W)=H(W)-disc(W).                 (1)
```

Let `u_x in {+1,-1}^V(T)` be the **actual incident sign root** of `x`; in
skew-matrix notation its coordinate at `v` is `K_C(v,x)`. Write

```text
Delta H_x=H(C)-H(T).                                        (2)
```

THM-4094 proves `Delta H_x>=0` and computes it from all insertion fibers and
orphan paths. THM-3729, equation (8), gives

```text
E_odd(T;u_x)=2 disc(C)-disc(T).                              (3)
```

Define the local deletion charge

```text
chi_x(C)=2 Delta H_x-E_odd(T;u_x)+disc(T).                  (4)
```

Substituting `(2)--(3)` yields the exact identity

```text
boxed:
chi_x(C)=2[(H-disc)(C)-(H-disc)(C-x)]
        =2[S(C)-S(C-x)].                                    (5)
```

Thus `chi_x` is an even integer. The sign condition

```text
chi_x(C)>=0  iff  S(C-x)<=S(C)                              (6)
```

is a genuine size-reducing certificate: deleting `x` does not increase the
slack which the open inequality seeks to keep nonnegative.

The rooted sign is load-bearing. Resetting it to the all-one vector changes
the energy in general by THM-3729. Likewise, a selected injection of old
Hamiltonian paths into new paths retains only `Delta H_x>=0`; THM-4094's
transitive-triple/`C_3` hostile proves that it loses the exact increment.

## 2. Local-unavoidability conjecture

> **Conjecture LU.** Every finite tournament `C` with `|C|>=2` has a vertex
> `x` such that
>
> ```text
> chi_x(C)>=0.                                               (7)
> ```

Conjecture LU implies the open inequality `H(C)>=disc(C)` for all finite
tournaments. The one-vertex base has `S=0`. Given `C`, choose `x` by `(7)`;
induction gives `S(C-x)>=0`, while `(6)` gives

```text
S(C)>=S(C-x)>=0.                                             (8)
```

Equivalently, a smallest counterexample to `H>=disc` would have

```text
chi_x(C)<0                 for every vertex x,               (9)
```

because all proper deletions would have nonnegative slack while `S(C)<0`.

LU is stronger than the desired inequality; no converse is known. A
nonnegative object may conceivably have every deletion with still larger
slack. THM-1950 reduces `H>=disc` to strong tournaments, but that theorem does
not prove `(7)` even on the strong stratum.

## 3. Averaged global-charge strengthening

Summing `(5)` over all vertices gives another exact identity:

```text
1/2 sum_x chi_x(C)
 = |C| S(C)-sum_x S(C-x).                                  (10)
```

This motivates the stronger conjecture

> **Conjecture ALU.** For every finite tournament `C`,
>
> ```text
> |C|S(C)-sum_xS(C-x)>=0.                                  (11)
> ```

ALU implies LU by averaging. It is the closest exact analogue to a
discharging total: the left side is a global charge, while one nonnegative
summand supplies a reducible deletion. It remains open. In particular, the
positive total charge of the apex-cubic paper supplies no inequality for
`(10)`.

## 4. Pointwise monotonicity fails at order four

Use upper-pair bit order

```text
(01,02,03,12,13,23),       bit 1 meaning i->j.               (12)
```

At mask `2`, vertex `3` is a source over the directed triangle

```text
0->2->1->0.                                                   (13)
```

The full tournament and its source deletion have

```text
(H,disc,S)(C)   =(3,2,1),
(H,disc,S)(C-3) =(3,1,2),
chi(C)          =(2,2,2,-2).                                (14)
```

Hence the tempting statement `chi_x>=0` for every vertex is false. The
strongest survivor is existential unavoidability: the other three vertices
in `(14)` still have positive charge.

## 5. Finite-exact frontier

The exact scout finds:

```text
universe                         all-negative rows   minimum individual chi
labelled n=2                               0                    0
labelled n=3                               0                    0
labelled n=4                               0                   -2
labelled n=5                               0                   -2
labelled n=6                               0                   -6
all 456 n=7 isomorphism reps               0                   -8. (15)
```

The minimum averaged half-charge in every row of `(15)` is zero, attained on
transitive examples. Among the 353 strong order-seven classes, the smallest
best local certificate is already positive:

```text
mask 171: (H,disc,S)=(27,4,23),
chi=(44,34,20,44,20,20,34),       max chi=44.               (16)
```

These are `FINITE-EXACT` facts only. The order-seven bank is checked against
its normalized-LF SHA-256 before use. No extrapolation in `n` is licensed.

## 6. Typed relation to the snark proof

The source/target ledger is

```text
source:      a local island/configuration in an apex counterexample;
target:      the deletion packet (C-x,u_x,full insertion fibers);
source test: local configuration is reducible;
target test: chi_x(C)>=0;
source global step: discharging makes one configuration unavoidable;
target open step: prove LU or ALU;
destroyed by scalarization:
             the covariant root, insertion multiplicities, and orphans;
required sidecar:
             u_x plus THM-4094's full deletion-incidence response.          (17)
```

This is an architectural analogy with an exact target identity, not a graph
transformation between cubic graphs and tournaments. The apex theorem's
Euler charge, cyclic embeddings, Kempe moves, and computer-checked
configuration atlas do not exist on the tournament side.

THM-4104 and THM-4099 suggest the next native operation: expand `(4)` into
the full Boolean ear-cut response rather than its mean. THM-4111 is a
firewall: two parents can have the same aggregate ear mean and different ear
images, so average growth alone cannot prove `(7)`.

## 7. Reproduction and next decisive tests

Reproduce the current frontier with

```bash
python3 -B 04-computation/tournament_deletion_slack_unavoidability_hyp9080.py
python3 -B -O 04-computation/tournament_deletion_slack_unavoidability_hyp9080.py
```

High-value next tests are:

1. derive a local sign decomposition of `chi_x` from the complete
   `Start/End/Q` ear tensor of THM-4104;
2. seek a double-counting proof of ALU, retaining orphan paths and the actual
   rooted Pfaffian energy rather than separate upper bounds;
3. restrict first to strong tournaments and classify equality or near-zero
   rows under deletion;
4. enumerate order eight by isomorphism classes before using random evidence;
5. construct a hostile with all `chi_x<0`, which would refute LU without
   refuting `H>=disc`, or prove a structural reason one cannot occur.

Neither LU nor ALU is a proved dependency of any current theorem. In
particular, **`H>=disc` remains OPEN.**
