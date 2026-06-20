# LRC14 AP-Tail Root Packets

The useful old artifact this session was HYP-2537: low-energy FKN shells are
root-packeted, not Hamming-uniform.  That note was about tournaments near the
transitive ground state, but it says exactly what the LRC AP-tail layer needs:
do not count perturbations; keep their addresses.

The AP-tail rows are not just "one replacement", "two replacements", "three
replacements".  They are packets:

```text
holes in the AP collar
tail insertions
Glaisher odd-shell carry
drop-6 mouth survival/damage
endpoint owner geometry
```

That packet viewpoint explains the recent exact results.  THM-543 has one
below-second exception because the packet `(6,10)->20` is pure existing-shell
carry and keeps all four drop-6 mouth intervals.  THM-544 clears the two-tail
layer because the next packets either damage the mouth or open new carry debt
before they can stay below `426/35035`.

KPS HYP-2661 sharpened this while I was working: the mouth seems owned by the
shell-1 dyadic tower `{1,2,4,8}`.  Sub-second AP-tail rows should preserve
odd-shell carry `{1:15}`; spending any shell-1 bit is the clean algebraic form
of mouth damage.

The S38 scan pushes this pattern one layer further.  Among `1,076,482` exact
three-tail rows with tails up to `35`, none lies below the AP-second threshold.
The best row is `(3,4,10,12)->(17,19,20)` with measure `4309/255255`, already
above threshold by `59/12495`.  It has no old-mouth survivor and carry damage
`{1:-4}`.  In other words: once the packet has left the retained-mouth/shell-1
template, it is no longer close to the dangerous floor.

The proof target now looks more precise than a global induction:

```text
full shell-1/mouth survival -> finite template ledger
shell-1/mouth damage or new-shell carry -> AP-second threshold paid
large tails -> HYP-2655/HYP-2658 recursion
```

This is also the right way to read the Glaisher/Witt bridge.  Glaisher carry
does not totally order the rows, and Witt/ghost closure does not replace
geometry.  The carrier is the closed packet: carry plus CRT cell plus mouth
ownership.  Scalar measure should be the last line of the proof, not the first.
