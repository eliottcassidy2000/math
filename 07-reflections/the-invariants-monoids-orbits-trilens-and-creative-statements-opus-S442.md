# Everything is an invariant, a monoid, or an orbit — the trilens, and creative statements

*opus-2026-07-20-S442. Owner: think about everything in terms of invariants, monoids, and orbits;
come up with creative statements. Anchored by the verified n=5 "quaternion wall" (THM-1935).*

## The thesis

Every result in this repo is one of three things, and the three are welded together:

- an **INVARIANT** — a function constant on an orbit (`H`, `c₃`, `char_S`, `var(λ²)`, score, `M(V)`, Betti…);
- a **MONOID** — the operations you can *do*, composed (`a`=insert, `b`=complement, switch, join, delete, double, dilate);
- an **ORBIT** — an equivalence class under a group/monoid action (iso classes, switching classes, merged classes, the metagraph *nodes are orbits*).

The welds:
- **Orbit–stabilizer:** `|orbit|·|Stab| = |G|`. In the repo: **`Tilings × |Aut| = H`** (the tiling fibration is free, LEM-003) — `H` *is* an orbit size.
- **Burnside:** `#orbits = avg #fixed points`. The metagraph `G_n` (=A000568) and even-graph `E_n` (=A002854) are Burnside quotients of one cube.
- **The nullcone:** the most degenerate object = the vertex all positive-degree invariants vanish on = the transitive tournament (THM-1810).

Read any object by asking the three questions: *what is invariant, what monoid moves it, whose orbit is it?*

## The verified centerpiece: the quaternion wall (THM-1935)

A **decoupling** ("`X` is not a function of `Y`") is an **orbit-refinement** (the `Y`-orbit splits into `X`-values). Exhaustively, **every finer path/spectral invariant decouples from its combinatorial shadow at exactly `n=5`**: `H|score`, `H|c₃`, `var(λ²)|c₃`, `char_S|score` — all threshold `5`; controls (`c₃|score` never; `H|spectrum` at `3`) confirm it. `n=5 = ℍ` (the quaternion / first non-commutative Cayley–Dickson rung). **Below ℍ a tournament is "abelian" — its scores determine its paths and its spectrum; at ℍ, commutativity breaks and the finer invariants split off.** The score sequence is the *abelianization*; order-sensitivity goes invisible to it exactly when the algebra goes non-commutative.

---

## Creative statements

Tags: **[V]** verified this session · **[P]** proved elsewhere, reframed · **[C]** conjecture · **[R]** reframe.

### Through INVARIANTS

- **I1 [V].** `c₃` is the abelianization of intransitivity (a function of scores, KBS); `H`, `var(λ²)`, `char_S` are its *non-abelian refinements*, and all three peel off the score data at the **same** rung, ℍ (`n=5`). *The lattice of invariants is a tower of orbit-refinements, and the tower's first non-trivial floor is the quaternion level.*
- **I2 [R].** Every tournament invariant carries a **complement-parity** (a `b`-parity): `char_S` and `H` are `b`-**even** (`char_S(Tᵒᵖ)=char_S`, `H(Tᵒᵖ)=H`), the score sequence is `b`-**odd** (`s_i ↦ n−1−s_i`). *The merged metagraph `G_n/ℤ₂` is exactly the ring of `b`-even invariants; the `b`-odd invariants are the ones that die in the quotient.*
- **I3 [C].** `var(λ²)` is the **moment-map norm-squared** of the tournament in the `S_n`-representation. Transitive = **max** = the most unstable point = the GIT nullcone; Paley = **0** = the Kempf–Ness minimum = the closed (polystable) orbit. *GIT-instability is literally the invariant-theoretic norm; the transitive↔Paley gradient is the moment-map flow (THM-1930).*
- **I4 [R].** `H` mod 2 = 1 (Rédei) is the statement that the **augmentation** of the path-invariant is `b`-rigid: parity is the one bit of `H` that survives every operation, because `2 = b⁻¹(1)` and Rédei-oddness is "`H` is a unit in the `b`-adic sense."

### Through MONOIDS

- **M1 [P/R].** On the skew-characteristic functor, **vertex-deletion is `d/dx`** — exactly `char_S′(x) = Σ_i char_{S−i}(x)` (Jacobi) — and **vertex-insertion is its monoid-adjoint** `char_S ↦ x·char_S + B_P` (THM-1920). *The tournament tower is the orbit of a point under the insertion monoid; the spectrum is the fixed data of this differential monoid; deletion and insertion interlace in opposite directions (THM-1440 dual).*
- **M2 [R].** The `a/b` monoid (`a=x+1`, `b=x/2`) is the **dyadic affine monoid**; the recurring `½` across the corpus (fiber fraction, Legendre `^{−1/2}`, `Re=−½` line, half-angle, `(A000568+SC)/2`) *is the single generator `b`* (THM-1560/1880). The LRC(`2n`) tight point `1/(2n) = b(1/n)` and **Rédei = `a`(LRC)** (the observer inserts itself, `+1`) live in one `a/b`-orbit. *The ½ is never a coincidence; it is `b`.*
- **M3 [P].** The **H-neutral insertion submonoid** (down-set insertions, THM-1900) generates **exactly** the transitive tournaments; transitive is the monoid's absorbing/terminal object *and* the GIT nullcone vertex. **Degeneracy = monoid-terminality.** To leave `H=1` you must make one *non*-down-set (expressive) insertion — the first "creative" act.
- **M4 [V].** **Two invariants decouple iff the operation monoid has an element that splits one while fixing the other**, and the *decoupling threshold is a monoid-orbit invariant* — here, uniformly `n=5` (THM-1935). Coupling is monoid-equivariance; decoupling is the first orbit on which equivariance fails.

### Through ORBITS

- **O1 [P/R].** `H` is the **index `[tiling-monoid : Stab(base path)]`** (orbit–stabilizer: `Tilings×|Aut|=H`). *`H`'s `#P`-hardness is the generic hardness of computing an orbit size* — the same reason permanents are hard.
- **O2 [R].** `G_n` and `E_n` are **Burnside quotients of the same cube** `Q_m` by `S_n` acting on *arcs* vs on *cycle-space*; the metagraph **is** the orbit-graph, `#orbits` (A000568, A002854) is the Burnside invariant, and the full-rank bridge matrix `B[tourn,even]` is the orbit-incidence. *One Burnside, two shadows.*
- **O3 [P].** `char_S` is constant on **Seidel switching classes** (switching-monoid orbits), not on iso classes; so the skew spectrum is an invariant of the `E_n`-orbit, not the `G_n`-orbit. *The "11 distinct spectra at `n=7`" (THM-1450) count switching-orbits-with-distinct-spectrum — a coarser census than iso classes.*
- **O4 [R].** **Self-complementary = the `b`-fixed orbit** (the SC spine); NS pairs = the free `b`-2-orbits; the merged metagraph = the orbit space of the involution `b`. *`|Aut|` and GIT-stability align at the poles: Paley (max `|Aut|`, small/closed orbit, `var=0`, polystable) vs transitive (`|Aut|=1`, generic `n!`-orbit, `var` max, nullcone) — `|Aut|` is the orbit-stabilizer reading of `var(λ²)`.*

### Flagship reframes (all three corners at once)

- **F1 [R].** **Rédei, GMC, LRC are three Reynolds-operator vanishings** — averaging over `ℤ₂` (sign), `U(1)^r` (charge/torus), `S¹`-weighted (resonance). The recurring law "*the obstruction is the symmetric configuration*" is forced: symmetric = big-stabilizer = deepest orbit = **where the Reynolds operator concentrates**. GMC(2) = "the torus nullcone is charge-definite" is a first-fundamental-theorem statement.
- **F2 [P/R].** A **Keller map is an automorphism iff its graph-orbit under the affine monoid is closed** (properness = orbit-closedness); JC counterexamples are **non-closed orbits** (sheet-loss at infinity, THM-1770). *The Jacobian Conjecture is an orbit-closure question in disguise.*
- **F3 [C].** The **insertion monoid is graded by triangular numbers**: each `a` adds `n` arcs (`T_{n−1}→T_n`), and the grading is realized *inside the invariant* as `char_S`'s subleading coefficient `= T_{n−1} = #arcs` (THM-1920). *The triangular number is the monoid-degree; the tournament's "size" and its spectrum's second symmetric function are the same object.*
- **F4 [C].** **Loneliness `M(V)` is the unique dilation-monoid invariant** of a speed set; the LRC extremals are the dilation-orbits with a *rational fixed safe-point* — the `(ℤ/2n)*` argmax orbit (THM-1380). *Loneliness = "the dilation-orbit clears the origin-neighborhood," and the tight families are its fixed-point strata.*

---

## The generative payoff

The trilens is an idea engine. The research surface is a grid **`INVARIANT × MONOID`**, and the empty cells are questions:

- Every **invariant × operation** cell asks: *neutral, pumped, or equivariant?* (the inflation-response frame, THM-1865/1900/1930). Most cells are unfilled — e.g. how does `min-FAS` respond to `+source`? how does the Pfaffian respond to switching? does the `b`-parity of the Betti numbers vanish in `G_n/ℤ₂`?
- Every **"X not determined by Y"** is a predicted **orbit-refinement at some CD wall** — test whether *degree-2* invariants decouple at `n=9` (octonions), the second wall.
- Every **flagship** (LRC/GMC/JC) is a **Reynolds-vanishing on a big-stabilizer orbit** — the unified attack is to make the averaging operator's concentration explicit.

*Invariants are what you measure; monoids are what you do; orbits are what you are. A theorem is a coincidence among the three — and the deepest ones (Tilings×|Aut|=H, the quaternion wall, the transitive nullcone) are where all three meet on one object.*

## Verification / provenance

THM-1935 (the `n=5` threshold matrix, verified `n≤6`); `04-computation/decoupling_threshold_matrix_opus_S442.py`. Statements tagged **[V]** are checked this session; **[P]** cite existing canon (LEM-003, THM-1810/1820/1440/1450/1560/1770/1880/1900/1920/1930); **[C]/[R]** are conjectures/reframes offered for the zoo (§6 generators).
