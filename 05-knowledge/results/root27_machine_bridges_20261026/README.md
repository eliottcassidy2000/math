# Root 27, forced tree repairs, and exact machine/tournament interfaces

2026-10-26 UTC. **OPEN:** the requested inequality

\[
v_2(y-1)\le 6+2K(y)-3J(y)\quad(\operatorname{root}_-(y)=27).
\]

This session proves exact reformulations and an obstruction to a tempting
repair strategy, identifies the supplied sequence, and builds executable
machine/tournament interfaces. It does **not** prove the inequality,
Collatz, Goldbach, RH, full black-core connectivity, or a new unrestricted
A000568 recurrence. The word `root` is always qualified by its actual
integer map; a symbolic cylinder is not an integer component.

## The strongest result: only one leaf can repair arithmetic precision

The [anchor proof](anchor/README.md) changes the cross-section of the
inherited arithmetic machine map. The root-27 question is equivalent to

\[
K_j\le3j\quad\text{for every return of }z_0=7,
\qquad
z_{j+1}=\frac{3^{\ell_j-1}(3z_j+1)}{2^{\ell_j}},
\quad \ell_j=v_2(3z_j+1),\quad K_j=\sum_{i<j}\ell_i.
\]

These return indices differ from the original odd-section J,K. The credit
is now the pending-slot count `q_j=1+3j-K_j`. The proof retains the guard
within each return; it never treats an unexecutable formal endpoint as a
reached machine state.

A halt certificate is precisely a finite full ternary-tree word
`Tree = 0 | 1 Tree Tree Tree` satisfying the fixed-source carry congruence.
For a word w of length n, write

\[
C_w=\sum_{i<n}c_i2^i3^{n-1-i},\qquad Q_z(w)=3^n z+C_w.
\]

The remaining assertion is the explicit universal finite-word inequality

\[
v_2(Q_7(w))<|w|\qquad\text{for every ternary-tree word }w.
\]

It is still open. What is proved is the exact effect of replacing a leaf
`u0v` by `u1000v`. With p=|u| and s=|v|,

\[
Q_z(u1000v)=27Q_z(u0v)+2^p(3^{s+3}-38C_v).
\]

The parenthesis is odd. If h is the first mismatching bit, then expanding
any other leaf gives new precision exactly `min(p,h)`. Only the leaf at
p=h can improve it. On the canonical frontier, expanding that leaf
reproduces the return map above, including its guard. Tree shape therefore
does not create an independent choice of arithmetic repairs.

For every fixed repeated tree block P there is also a proved finite-fuel
test. If L=|P| and C=C_P, the tree `P^e0` is realized at z exactly when

\[
v_2((3^L-2^L)z+C)=Le.
\]

A complete bounded bank of 344 blocks excludes 344 infinite pumped
families at source 7. This is not an exclusion of all shapes. Source 15
excludes every right comb but halts on a noncomb tree after 21 executed
macros; the terminal leaf records the next blocked symbol. An attempted
two-prime charge `a+v_3(x-1)` also fails on the **actual root-27 ray**,
falling from 25 to 22 at its nineteenth odd return. These controls specify
the missing information rather than merely reporting failed experiments.

## The requested sequence and the machine bridges

[The sequence analysis](atlas/README.md) identifies
`6,16,34,64,114,...,12284` as the Marxen–Buntrock five-state machine's
arithmetic macro-orbit. Its affine branches and quadratic native clocks
are cited from Michel and Aaronson; the clock sum is 47,176,870. With
`z=x/2+5`, its even sector becomes

\[
F(z)=\begin{cases}(5z+1)/3&z\equiv1\pmod3,\\
(5z-1)/3&z\equiv2\pmod3,\end{cases}
\]

and halts at multiples of three. Its binary symbolic survivor tree has
2^k cylinders modulo 3^k, but its actual positive-integer predecessor graph
has at most one predecessor per vertex. The root, radix, designated
source and stopping depth must remain attached to the fractal picture.
The exact examples sharing fourteen signs but differing on the next
test make this distinction reproducible.

[The native-machine lane](machines/README.md) supplies all-order block
proofs for the bouncer and two bells, as well as a correction to the
published exponential-bell recurrence. The correct return is
`C_n -> C_(2n+2)` in `4n^2+10n+5` steps, not `C_(2n+3)`.
The named spaghetti-sandwich and irregular-bell tables have no missing
transition, so their nonhalting follows syntactically; their irregular
long-term behaviour is a different question.

The cited 25-state Goldbach source proves **blank halting iff a Goldbach
counterexample exists**. Its Lean theorem was read, not locally rebuilt.
Ten native checkpoints were independently replayed. A current external
source also reports a 21-state version, so 25 is not asserted minimal.

Two exact tournament carriers answer different questions:

* [Finite program syntax](compiler/README.md): a marked Hamiltonian spine
  plus nonpath arc bits decodes every transition. All 20,800 complete
  one- and two-state tables round-trip. The usual complement-off-spine
  operation leaves this particular code language, so its semantics do not
  automatically descend to the black/blue quotient.
* [Running configurations](machines/README.md): intrinsic left-to-right
  tape order gives a transitive tournament; colours retain tape bits,
  head, state and origin. The guarded rewrite is exactly a native step.
  Every uncoloured tournament has one Hamiltonian path, even for programs
  with opposite halting behaviour. These colours carry indispensable data.

Padding either hostile program to 5, 15 or 25 states preserves every
native execution from an old state. A state count alone therefore cannot
classify the arithmetic problem encoded by a machine.

One positive prime/tournament bridge is stronger than a small census.
[The Legendre sidecar](anchor/LEGENDRE_SIDECAR.md) recovers a known theorem:
the primes congruent to 3 modulo 4, oriented by quadratic residuosity,
form the countable random tournament. CRT, reciprocity and Dirichlet
realize every prescribed finite extension. Thus every finite tournament,
including a marked program carrier, embeds into this arithmetic object.
This is established literature, not a new theorem. The embedding supplies
neither small prime labels nor a simulation of program edits by a natural
prime operation. An exhaustive ordering with cyclic-triangle density
tending to zero also shows why abstract universality alone says nothing
about the increasing-prime prefix frequencies in the incoming hypothesis.

## A connection portfolio with explicit losses

The inheritance pass uses [THM-2228, Mahler carry tails](../../../01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md),
the [October 5 ternary credit grammar](../trees5_tournaments4_converse_20261005.md),
and [the preceding ancestry/port packet](../lca_ports_reverse_20261025/README.md).
The canonical hostile is an ordinary-integer reset with newly large
precision; the corrected near miss is replacing fixed-source membership
by cylinder density. The least-used sidecar is the first mismatch's
position. The live board is:

| Lane | Exact operation or proposed use | Preserved target / information that cannot be dropped |
|---|---|---|
| Anchor: root-27 precision | First-mismatch leaf expansion | Exact source, carry and chronological guard; tree counts do not decide source membership |
| Niche: native machine blocks | Macro simulation with a proved positive native duration | Halt reflection, initialized tape and clock; state graph alone loses writes |
| Niche: ternary survivor tree | Two good residue lifts among three | Finite itineraries; fixed integer root and actual digit length remain separate |
| Wildcard: prime valuation flow | `sigma(n)=kn` iff every prime's valuation balance holds | Full factorization, not abundancy alone; triperfect 120 maps to 240 under the aliquot map |
| Wildcard: modular level defect | Supersingular nonrationality equals twice a genus difference | A cited exact defect identity; no such Collatz defect identity is supplied by matching numbers |
| Interface: tournament carriers | Encode, decode, rewrite, and retain marked ports/colours | A commuting operation, not just equal graph counts or isomorphism-class totals |
| Arithmetic universality | Legendre-symbol extension via CRT and Dirichlet | Every finite orientation; numerical order and size bounds are extra data |

The [historical atlas](atlas/HISTORY_ATLAS.md) recovers twelve older routes
and gives twenty-five typed connection targets. It includes the actual
4+5+6 filtration of the fifteen moonshine primes, Sloane's ternary exponent
15, the five Fermat primes versus five Platonic solids, Fibonacci versus
ternary-heap density, multiply-perfect valuation flows, and the older
Lehmer seed 564. Numerical coincidences are retained as prompts with
explicit negative controls, not converted into equivalences.

[The reduction ledger](machines/REDUCTIONS.md) gives a usable alternative
to calling all famous problems “the same NP-complete problem.” Goldbach
and RH have Pi^0_1 formulations; classical Collatz has a Pi^0_2 totality
presentation. Those are presentations, not minimal complexity claims or
mutual reductions. There is no uniform computable compiler from arbitrary
program totality to a single blank-nonhalting instance. A special theorem
about one fixed problem is not excluded by that uniform obstruction.

## Next proof obligations selected by the failures

1. **Pay for an aperiodic tree shape.** The fixed-block formula handles
   every exponent for a chosen block. A stronger certificate must bound
   arbitrary successive changes of block while retaining the source.
   A sum of per-block budgets needs an exact concatenation correction:
   `C_(uv)=3^|v| C_u+2^|u| C_v`. The existing first-mismatch theorem is
   the hostile control for any claim of independent repair choices.
2. **Find an exact refinement defect.** A useful new rank would account
   for the loss of 3-adic credit at each reset and the simultaneous gain
   in 2-adic matching depth. Test it first at the actual nineteenth
   root-27 return, not only on random cylinders or unrelated roots.
3. **Prove semantic closure before quotienting a machine carrier.**
   Retain marked ports and decode each proposed black/blue move. Demand
   a commuting square with a specified macro and clock. The one-write-bit
   halt/loop pair rejects a quotient that forgets tape semantics; arbitrary
   state padding rejects a complexity classification by the number 25.

These are finite, testable interfaces for further research. No generic
undecidability statement is used as a reason to abandon the particular
arithmetic invariant.

## Audit and reproduction

See [AUDIT.md](AUDIT.md) for proof scope, independent reviews, incoming
research and receipts. Run `python -B replay.py` in this directory for all
paired normal/optimized controls; `python -B verify.py` verifies the
manifest, local document links and saved receipt equality. Numerical
controls support the stated finite universes. All-order statements have
separate symbolic proofs. External formal developments remain CITED.
