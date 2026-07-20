# The Shear-Delay Law and the Relational Triad

**boxeph-2026-07-20-S148** · companion scripts: `04-computation/triangular_continuations_wide_birth_boxeph_S148.py`, `..._laws_boxeph_S148b.py` (HYP-8175)

The owner's directive: the triangular numbers "arise from relation itself, the
edges of complete graphs; the nth triangular number is n in the relational
sense." Take that literally — T(n) = C(n,2) counts the relations among n
things — and every "continuation of the triangular numbers" becomes a
**deformation axis of the concept of relation**. This session built eleven such
families, sheared them down 0/1/2/3, summed and multiplied the columns, and ran
every derived sequence through an exact recognizer. Five laws came out, three
with proofs. But the frame matters more than the laws.

## 1. The six deformation axes

| axis | family | triangular at | what "relation" becomes |
|---|---|---|---|
| dimension | C(n+d−1,d) simplicial | d=2 | d-ary incidence (polyhedral view) |
| step | k-gonal P_k(n) | k=3 | gnomon growth (polygonal view) |
| exponent | Faulhaber S_p(n) | p=1 | weighted relations (S147's third view) |
| ground field | [n,2]_q | q=1 | F_q-linear relations (2-subspaces) |
| truncation | R_d(n)=Σ_{i≤d}C(n,i) | d=2 | bounded-arity relations (Moser/lazy-caterer) |
| group-meaning | C(n,2) / S(n,2) / c(n,2) | — | subset / partition / cycle sense of "two groups" |

The last axis is the deepest and is the **relational triad**: the three
fundamental ways to organize n things into two groups are 2-subsets (C(n,2) =
triangular), 2-blocks (S(n,2) = 2^{n−1}−1 = **Mersenne**), and 2-cycles
(c(n,2) = (n−1)!·H_{n−1} = the **harmonic numerators**). All three surfaced
independently in the owner's S147 riddle — the triangular column, the
Mersenne–Moser break 1,3,7,15,29, and the harmonic series H_n — because the
riddle was walking the triad without naming it. At the top of the three
triangles sit the three total-organization counts: shear-0 row sums = **2^n,
Bell, n!** — the exponential, the intermediate, the factorial world.

## 2. The Shear-Delay Law (proved)

Moving column x down t·x rows and summing:

**A_t(m) = Σ_x C(m−tx, x) obeys a(m) = a(m−1) + a(m−t−1).**

Proof is one Pascal identity resummed. The tower: t=0 → 2^m, t=1 → Fibonacci,
t=2 → Narayana's cows, t=3 → A003269, … (verified through t=6). **Shear trades
growth for memory**: each unit of downward shear converts one factor of
"branch now" into one step of "wait t steps, then branch." The owner's
observation that shear-2 "makes it Fibonacci-analogous" is the t=1 rung seen
from below; the whole ladder is delay-(t+1) generalized Fibonacci. On the
square Pascal grid the same shears appear one step later (shear-t of the square
= shear-(t−1) of the triangle) — the cross-family coincidence engine caught
this automatically.

In the other two relational worlds the same shear-1 produces the Fibonacci
analogues Σ_x S(m−x,x) = 1,0,1,1,2,4,9,22,… and Σ_x c(m−x,x) =
1,0,1,1,3,9,36,… — and the detector proves neither is C-finite (no linear
recurrence of order ≤ 6). **C-finiteness dies as relational richness grows**:
the set world is rational, the partition and permutation worlds are not. This
is the same signature as our tournament metagraph story — the objects that
count subsets stay linear-recursive; the objects that count structures do not.

## 3. Proth laws: the owner's grid is a Mersenne machine (proved)

The owner's n·2^x+1 grid, sheared:

- **shear-1 row sums = 2^{m+1}−1, exactly Mersenne.** (Telescoping; asserted
  in-script.) The grid whose rays are our two key prime families (2n+1
  pigeonhole moduli, 2^x+1 Fermat/gate primes) has Mersenne as its *sheared
  shadow* — the third prime family appears as soon as you triangulate.
- **shear-2 row sums** = 1,2,4,7,11,18,26,41,… with closed form 2^{K+2}−K−3 /
  3·2^{K+1}−K−4 (m = 2K / 2K+1) and characteristic polynomial
  **(x−1)²(x+1)(x²−2)** — a √2-geometric strand (the "half-tower"), a parity
  strand, and a linear drift. The detector found order 5; the factorization
  confirms the hand derivation. Note the Lucas overture 1,2,4,7,11,18 before
  the √2 strand takes over — a genuine incidental overlap, not a law.

## 4. q-confluence (proved): triangular numbers are a collapsed chord

[n,2]_q = a + b·q^n + c·q^{2n} with a = 1/((q²−1)(q−1)), b = −1/(q(q−1)²),
c = 1/(q(q²−1)(q−1)); characteristic polynomial **(x−1)(x−q)(x−q²)**. As q→1
the three geometric strands collide into a triple root and the mixture
degenerates into the polynomial n(n−1)/2. **The triangular numbers are the
confluent limit of a three-note geometric chord {1, q, q²}.** The q-axis
resolves the degeneracy — exactly the pattern we know from the LRC certificate
ladder, where the rational rungs 2/p resolve into distinct mod-p strands, and
from the JC work, where the caustic is the confluence locus of fiber sheets.
Confluence = where counting becomes polynomial = where distinct exponential
mechanisms become indistinguishable. The interesting arithmetic (Mersenne at
q=2 in S(n,2); the 2-subspace counts 1,7,35,155 at q=2 here) lives just off
the confluence point.

## 5. Moser strand (proved) and the owner triangle's ledger

- Truncation grid sheared 1: Σ_x R_x(m−x) = **F(m+3)−1** (hockey-stick +
  shallow diagonals). The Moser family, triangulated, is Fibonacci-minus-one.
- Owner triangle row sums = Faulhaber-shear-1 sums + 2^{m−6} for m ≥ 6
  (106 = 105+1, 317 = 315+2), and its shear-1 sums equal Faulhaber-shear-2
  sums until the Pascal cascade joins at m=9. The riddle triangle is
  *literally* a sheared Faulhaber grid plus a delayed Pascal echo — the S147
  closed form restated as a sum identity. (Caveat: rows ≥ 8 use the closed
  form, verified only on the 28 given entries.)

## 6. Where this points

The shear-delay law is a **one-parameter interpolation between growth and
memory**, and the delay parameter is the same object as our tower heights
(2^a·7^b·p lifts; the x-direction of the Proth grid). The confluence principle
says polynomial laws sit at collision points of exponential strands — which is
how we should read every "polynomial in n" invariant on G_n: look for the
q-deformation whose strands collide there. Candidate next steps: (i) q-deform
the tournament count via [n,2]_q-weighted arc counting and see which metagraph
invariants become geometric mixtures; (ii) the partition/cycle Fibonacci
analogues as tournament objects (what counts "Hamiltonian paths" in the
partition world?); (iii) primes in sheared rows (Proth-row products) as a
testbed for the gate-prime heuristics.
