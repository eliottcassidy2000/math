# Chan–Pak: physical carries from pointed supports

The field-with-one-element analogy is used literally at the level of
finite choices: supports, directed marks and Cartesian products, not
formal substitution of q=1 for a finite field. It yields genuine
tournament carry digits and also exposes a structural obstruction.

## What is proved

- [The smallest possible puncture count for a unit switch](MINIMAL_SWITCH.md)
  is 5, with literal augmented counts 579 and 581 and identical
  six-coordinate boundary data. The motif moves three real contacts.
  External parities make the switch neutral or change its count by two.
- [Pointed support decoder](SUPPORT_DECODER.md): for prescribed
  (H,F2,P) it builds or rejects physical punctures without an arbitrary
  order cutoff. In the triangle-only family, the forward word and
  backward signed digits recover every mixed triangle slot. Fixed
  (H,F2,P,Q) then fixes K, so those triangle assemblies cannot supply
  an independent carry choice.
- [Physical carry bank](f1.md): matched 5- and 11-path blocks have
  consecutive normalized gains 11 and 10. With a≥9 and b≥10, every
  normalized carry in [90,11a+10b−90] has literal subset choices, and
  moving the reader realizes every baseline residue modulo 5^a·11^b.
- [Absolute compiler](ABSOLUTE_COMPILER.md): given a prescribed N,
  find a witness or certify rejection across all widths of this
  entire ordered 5/11 bank. The proof keeps an actual baseline and
  an actual pair of subset ranks, and checks their sum, not just a
  modular condition.
- [Induced enlargement barrier](bank.md): if an induced extension
  raises H by Δ>0, the old count is at most 2Δ (sharply 2Δ−3
  for Δ≥6, with infinitely many equality cases). Adding any number
  of vertices without changing H only adds singleton strong
  components. Therefore an exact small carry at large H really
  must change old material, not merely add a longer reader.

The open step is **absolute alignment**: surjective baseline residues
and a physical window at each selected baseline do not imply that some
window contains a given N. The compiler is complete within its bank;
completeness of the bank, or a replacement construction covering the
rejected admissible values, is not proved here. Chan–Pak remains open.

## Follow-up

[Right-only switches](RIGHT_ONLY_SWITCH.md) have minimum puncture
count 3 rather than 5. A single old arc gives the literal pair
40,241 and 40,243, and finite probes characterize the precise
boundary data seen by right-only and two-sided readers.

[Positional five-way carry](POSITIONAL_FIVE.md) gives a physical digit
word at each residue modulo 5^k, with minimum eleven-vertex blocks
in this model. Both wings may move; a final quotient test decides
absolute attainment for the selected wings.

[Joint positional five-way carry](JOINT_FIVE.md) selects *both*
ordinary and rooted responses modulo 5^k, including nonunit reader
residues. Separate first and later alphabets reduce the order, and
actual valuations come from a physical appended reader. Absolute
coverage remains open.

## Reproduction

Run **python3 replay_support.py**, **python3 replay_f1.py**,
**python3 replay_carry_bank.py**, and **python3 compile_absolute.py**.
All use the Python standard library. The small tournament profiles
are checked by literal subset counting, and large-address readouts
use the exact transfer identities.
