# Independent audit of the Session 188 fixed-quotient obstruction

**Verdict: the scoped all-order obstruction checks out after one prose
definition correction.** The circuit-branch worker performed this audit
independently of the producer. No global spectrum or arbitrary-replacement
claim is supported or inferred.

Sections 1-5 audit the original transitive-family theorem. Section 6 audits
the subsequent arbitrary-tournament-module extension, which supersedes the
transitivity restriction discussed in Section 3.

## 1. Correction found and verified

The first version of `returns/PROOF.md` described F as ending the "first
nonempty slot" in W, and G similarly. That is not the observable used by
the program or the inherited theorem: if slot one is empty, slot two must
not contribute to F. The producer corrected the definitions to:

- F: the first slot is nonempty and ends in W;
- G: the second slot is nonempty and starts in W.

The corrected text was reread. All computations already used these correct
definitions, so the source verifier and its receipts needed no change.

## 2. All-order symbolic mechanism

I read `carrier.py`, `skeleton.py`, `family_coeffs.py`, `inherited.py`,
`verify.py`, and the complete receipt. The finite skeleton is a genuine
all-order certificate rather than interpolation:

1. Delete the clone module from a two-slot cover. This leaves one of 5!
   exterior permutations and a cut among its six positions.
2. Clone runs may occupy exactly the gaps whose two exterior incidences
   permit passage through the homogeneous module. A backwards exterior
   adjacency must receive at least one clone.
3. A transitive clone run has a unique internal order. Therefore assigning
   each labelled clone to an allowed gap is a bijection with the missing
   path pieces. Requiring r specified gaps nonempty gives the stated
   inclusion-exclusion expression.
4. For clone x_i at a first-slot end, only the i smaller clones may share
   its forced terminal gap; none of the j larger clones may follow it.
   Thus the two capacities are k+1 and k before mandatory-gap exclusions.
   For clone x_i at a second-slot start the capacities reverse. These are
   exactly the base pairs used in `skeleton.py`.
5. The both-external, clone-end/external-start, and external-end/clone-start
   cases are disjoint and exhaustive for a singleton clone port. A single
   clone cannot be in both disjoint path slots.

The resulting singleton formula and all seven first-clone coordinates are
consistent with this bijection. In particular, the 0^0 terms at m=1 are
retained in the symbolic certificate; the simplified m>=2 family is not
incorrectly applied at the contraction endpoint.

For i>=1, the displayed bracket obeys

    18*3^i + 5 - 23*2^i >= 4*2^i + 5 > 0.

The remaining singleton-difference terms are positive. For i=0 the singleton
increment over the external-only J is `23*2^(m-1)+4>0` when m>=2. Because
the full J is a quadratic sum of nonnegative named endpoint-pair counts,
adding a second clone port is strictly increasing. This proves both the
all-markings lower bound and its unique equality mask `{x_0}`.

## 3. Exact mask and physical scope

The five exterior vertices are fixed as original names (1,2,4,5,6), and
their ports remain exactly {1,5}. Only the clone-port subset S varies.
The no-return theorem requires m>=2 and S nonempty. It does not cover
changing exterior ports, an empty clone-port set, a nontransitive replacement,
or another quotient adjacency. "Arbitrary masks" must be read in this
specified clone-only sense.

The theorem actually tests more masks than those preserving the original
old3-marked/old0-unmarked assignment, so it includes every such preserving
inflation. It does not assert that every tested mask preserves that old
marking. The first-clone family has the literal old-name embedding described
in the proof. Each added clone has homogeneous exterior incidence and is
inserted immediately before old0; all old arcs and old marks survive.

The quotient's directed cycle certifies strongness, and inflating a vertex
of a strong quotient by a nonempty module preserves it. The fixed external
ports are nonempty. No native Hamilton-complement legality is established.

## 4. Ordered contexts and clock shifts

I checked the actual inherited seven-coordinate law in the pinned Session
182 proof, rather than assuming the contextual inequality. For X -> Y it
gives exactly

    c_(X->Y) = h_X*c_Y - J_X*h_Y - e_X*s_Y - G_X*F_Y.

All subtracted terms are nonnegative. Hence arbitrary contexts X -> B -> Y
give `c_D <= h_X*h_Y*c_B`. Because c_B<0, contexts cannot increase it.
All m>=4 are therefore excluded by the sharp bound `c_B<=-137`.

The three remaining small signatures have c=-35,-59,-69. For either target
C=-103 or -101, each satisfies `3c<C<c` and `c-min(F,G)<C`.

- A nontransitive context has at least three Hamilton paths, so the first
  inequality excludes it. A directed triangle supplies three paths, and
  deterministic vertex insertion injects them into full paths; transitive
  tournaments have exactly one.
- A marked transitive left context has G_X>=1 and incurs a subtraction
  of at least F_B.
- If the left context is unmarked and the right is marked transitive,
  G_(X->B)=G_B and F_Y>=1, incurring a subtraction of at least G_B.
- Unmarked transitive contexts leave c_B unchanged.

This covers arbitrary sizes and arbitrary context markings. The full
attachment need not be assumed strong for this obstruction, since its
necessary count equation already fails.

For the clock statement, the target carrier, its marks, its ordered contexts,
and its integer shift q are fixed independently of t. If

    c_D + 2^(t+q)*F_D = 3*(-35+2^t*252) + d

holds for all sufficiently large legal clocks t, subtracting two consecutive
instances forces `2^q*F_D=756`, then `c_D=-105+d`. The constant condition is
therefore independent of every fixed integer q, including negative shifts
once t+q>=1. The result does not exclude a one-clock numerical coincidence,
an arbitrary t-dependent change of target carrier, or a non-shift clock map.

## 5. Independent exact checks

`verify_independent.py` imports no producer module. It builds the inflated
tournaments directly from the original seven-row adjacency data, enumerates
every literal directed simple path recursively, and assembles the named
endpoint counts from complementary visited sets. It uses no subset-DP
aggregation and no symbolic gap coefficients.

Normal and optimized runs match. The independent checks cover:

- m=1,2,3,4, including the newly checked nine-vertex m=4 family;
- all 26 nonempty clone masks in those sizes;
- the complete seven-coordinate first-clone formulas and every singleton
  position formula at those sizes;
- every literal row and signature in the producer's m=2,3 receipt;
- all 75 marked tournament context presentations of order at most three,
  plus absence, on each side of the three analytic survivors: 16,875 pairs.

The numbers of literal directed nonempty simple paths enumerated at these
four sizes are 139, 385, 1,063, and 2,977. In addition, every producer-source
hash matches its receipt, and its stored normal/optimized outputs agree.
`producer_replay.json` records a fresh normal and optimized execution in an
isolated temporary directory, checking that its regenerated receipt agrees
byte-for-byte with the supplied original. No original source or receipt is
modified by that replay.

The independent finite checks corroborate the theorem. The all-order proof
rests on the bijection, positivity, contextual law, and fixed-clock coefficient
comparison audited above. No other mathematical correction was found.

## 6. Stronger arbitrary-module extension: independently verified

The producer subsequently added `PROOF.md` Section 7. This strengthened
claim also passes: any tournament module A on m>=2 vertices, with any
nonempty internal port subset U, is excluded in the same fixed quotient
with the same fixed external ports {1,5}, even after arbitrary marked ordered
contexts and a fixed integer clock shift.

I checked `module_coeffs.py`, `module_extension.py`, the new proof, and the
new receipt. The nonnegative skeleton coefficients are consistent with an
independent conversion of the earlier empty-slot formula. Let H, P_2 be
the module's one-path and two-possibly-empty-slot counts, and let tilde E_k
and tilde S_k permit empty slots except the distinguished endpoint slot.
Writing the new nonempty-slot quantities as E_k and S_k gives

    P_2 = 2h + q2,
    tilde E_2 = E_1 + E_2,
    tilde E_3 = E_1 + 2E_2 + E_3,
    tilde S_2 = S_1 + S_2.

Substitution in the earlier symbolic endpoint identity yields exactly

    c_single = 9h+3q2-25E_1-44E_2-18E_3-S_1-5S_2.

The nonnegative occupied-gap construction supplies the same identity
directly for arbitrary modules: each occupied gap receives a nonempty
directed path, the forced distinguished gap is put first in the module cover,
and the other paths retain exterior-gap order. This is reversible and does
not assume the module is transitive or that its paths are unique.

Both injections are valid with their recorded data. Splitting a Hamilton
path after u yields an E_1 or E_2 object and is reversed by concatenation.
Splitting u's path in a two-nonempty-path cover yields an E_2 or E_3 object;
retaining its original slot bit is necessary and sufficient for inversion.
Thus `h<=E_1+E_2` and `q2<=2(E_2+E_3)`. The coefficient arithmetic is

    c_single <= -16E_1-29E_2-12E_3-S_1-5S_2.

The singleton-u construction gives E_2,S_2>=1, and gives
`E_3>=q2(A-u)>=2^(m-1)-2`. The last inequality counts each nonempty proper
subset and its complement, selecting a Hamilton path in each induced
tournament. Consequently

    c(B[A,U]) <= -10-12*2^(m-1) <= -106   (m>=4).

For m=3 the only additional isomorphism type is C3. Its quantities are
`(h,q2,E_1,E_2,E_3,S_1,S_2)=(3,6,1,2,2,1,2)`, so c=-115. Adding ports
only lowers c. The remaining m=2,3 transitive masks are exactly the cases
already excluded by the ordered-context argument. No larger-module census
is needed for this proof.

`verify_modules_independent.py` imports no producer code. It enumerates
literal module permutations and every nonempty cut pattern, explicitly
checks the two injections for collisions and valid targets, and independently
enumerates literal paths in every substituted tournament. Its normal and
optimized runs agree on all 75 labelled modules through order four, all
1,023 nonempty masks, 285 marked singleton vertices, and 5,569 individual
injection images. It checks the module receipt's source and proof hashes.
The exact receipt is `modules_receipt.json`.

The strengthened theorem does not retain the old internal arcs for every
arbitrary A; it covers replacements, so this is a deliberate larger universe.
It retains the exterior quotient incidences and external port set. Those
restrictions, nonempty U, m>=2, and a target independent of the clock remain
essential scope conditions. A nontransitive replacement of this pair alone
is no longer a viable escape from the proved obstruction.
