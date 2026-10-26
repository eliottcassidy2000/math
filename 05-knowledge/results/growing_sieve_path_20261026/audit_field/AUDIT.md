# Independent audit of the sixteen-state tournament observer

**Verdict: PASS for the stated elementary deductions, counterexamples and
finite universes.** No correction is required. This is an independent review
of the field lane, not certification of the unrelated global conjectures
mentioned for context. Producer files were not modified.

Reviewed producer bytes:

| File | SHA256 |
|---|---|
| `../field/README.md` | `5dd8eaf5fb47c0029cb291d1aa53da843d0719714bb4aacccc68a539b9cef898` |
| `../field/check.py` | `032410f418a4a3303220e0f31fc401c5e35b21492071e5ddc95897ee7b60325c` |
| `../field/normal.json`, `optimized.json` | `4a770d4f198dd5509107530c0c8c458b1050bac043672d90c0d639b45710aaa3` |
| `../field/SOURCES.json` | `d71df0a0fa76fae87539c0e75abafeb35f1e34cc206c198def2ced141b0d2230` |

## Proof review

1. **Bracket transport and the root gauge.** Each representative occurs in
   two brackets, so its rescaling changes K by a square. A projectivity
   multiplies the triple product by its determinant cubed; this induces
   precisely the determinant square class. Reversing an arc changes K by
   −1, a nonsquare modulo 7. At infinity K is v−u. Thus the stabilizer of
   the marked orientation is exactly the order-21 square-affine subgroup,
   inside the order-42 pole stabilizer. The coset description uses left
   actions and left cosets consistently. The determinant-one U maps
   infinity to t=2, so conjugation is necessary for the actual root, whose
   reciprocal common slope is r=4. Distinct poles give different labelled
   vertex domains; the two sheets at one pole give converse arc sets.

2. **Literal tournaments and scope.** The affine subgroup acts regularly
   on the 21 arcs, and on the 21 reversed pairs. The two invariant
   tournaments are exactly the two orbitals. Direct permutation enumeration
   independently gives 189 directed Hamilton **paths** and 21 automorphisms,
   hence tc=9, for this transported family. The assertion that a finite
   tournament automorphism group has odd order follows from Cauchy's theorem
   and the impossibility of a transposed vertex pair. It does not identify
   the nonsquare anti-isomorphisms with automorphisms. The joint Frobenius
   stabilizer statement over F8 uses its prime extension degree three:
   a coordinate outside F2 has orbit length three, whereas the determinant
   character is trivial on every power of t↦2t.

3. **The wheel refinement has the declared carrier.** Integer branch
   matrices are invertible modulo every modulus; observing them gives a
   closed transition system. Primitive residue vectors alone need not be
   reachable at an arbitrary modulus, as the note now explicitly states.
   Opposite parity and the mod-7 consistency gate are retained as necessary.
   The complete 3072-state assertion is restricted to M=15. The independent
   calculation below proves it by a different, larger carrier. On the
   D-ray, period 105 and 21 surviving states follow literally from the
   residues, without identifying these 21 states with the 21 tournament arcs.
   Prime self-hit exceptions are explicitly retained.

4. **Arithmetic witnesses.** Both D^5 and D^215 have genuine twin neighbours,
   identical (t,e,m,n) observations modulo (7,2,15,15), and opposite strict
   inequalities c²>rad(abc). Their radicals and every factor were checked.
   The words UUUDUDDDU and DUUDDDDDUA have the stated primitive Euclid pairs,
   opposite A-parities, and identical full prime support and radical.
   These disprove the specified functional-descent claims. They neither
   contradict abc nor prove a universal independence or distribution claim.

5. **Monoid maps.** Commutativization retains A-count parity but loses the
   endpoint, witnessed by UA and AU. The contracted algebra of the pointed
   free commutative monoid is correctly distinguished from an ordinary
   monoid algebra that would retain an extra basis element for the absorbing
   zero. Augmentation, sign evaluation and prime-weight evaluation are
   distinct maps of that algebra. The proper-divisor polynomial omits only
   the full-support monomial, so its three evaluations are exactly as
   stated for squarefree N and k≥1. Every homomorphism from an idempotent
   union monoid to a group is trivial by cancellation of y²=y. The actual
   radical witness separately obstructs arbitrary functional descent of
   this particular bit; the text does not confuse the two obstructions.

## Independent exact verification

The producer replay passed all **18,007 gates**, normally and with `-O`, and
both replays match its frozen receipt byte for byte.

`check_independent.py` imports no producer code and passes **6,499 explicit
gates** in both modes. Its principal independent routes are:

- Build the 16 tournaments by pulling Paley7 back through determinant-one
  pole charts, then compare with the bracket formula. Reconstruct all 336
  projectivities from the ordered images of infinity, zero and one, rather
  than enumerating all field matrices. Check every tournament transport.
- Count Hamilton paths by literal 7! enumeration, rather than subset dynamic
  programming; separately exhaust the full automorphism group.
- Reach the full signed primitive-vector carrier modulo **105**, of size
  **18,432**. Its projection has exactly six preimages over each of the
  **3072** claimed projective-plus-mod-15 states. This supplies an independent
  finite proof of reachability at M=15 and makes no arbitrary-M assertion.
- Recompute both arithmetic witnesses with independently generated primes,
  full factor reconstruction, direct integer branch formulas, and exact
  squared-height comparisons. Check the commutativization and divisor-sign
  controls independently.

Independent checker SHA256:
`eca7a85a519f59ba2cc494e44f5d1f2c50968bd5cb074287b45a49ad634d894f`.
Independent normal and optimized receipt SHA256:
`b10502455c41cc76fb7619adf7cf1fbb7900342dbec269632787ff5a12f48a1a`.

Reproduce from this directory:

```sh
python3 -I -B check_independent.py > replay.json
python3 -I -B -O check_independent.py > replay.optimized.json
cmp replay.json normal.json
cmp replay.json replay.optimized.json
python3 -I -B ../field/check.py > producer.replay.json
cmp producer.replay.json ../field/normal.json
```

No twin-prime infinitude, global abc bound, unrestricted normalized-Hamilton
exclusion, arbitrary-modulus reachability theorem, or literal field of one
element is certified or needed. A subsequent status/link-only promotion of
the reviewed note is covered by this audit; any mathematical edit requires
renewed review.
