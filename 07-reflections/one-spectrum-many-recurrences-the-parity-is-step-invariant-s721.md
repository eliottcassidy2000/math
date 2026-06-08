# One spectrum, many recurrences; the parity is what survives stepping (S721)

Turning the transfer-spectrum mindset onto the open questions did what a good frame should: it confirmed
one conjecture, corrected one overreach, and handed over a clean new fact, all from the same small
computation. The instruction was to use the mindset to make progress, and the progress is that the ledger
got sharper.

The clean new fact first. The staircase tournament with every non-base tile pointed upward — the
all-upward family, the one THM-337 studied at even sizes — has Hamiltonian-path count
`1, 3, 5, 9, 17, 31, 57, 105, 193, ...`, and it is Tribonacci: each term is the sum of the previous
three, characteristic polynomial `x^3 - x^2 - x - 1`, growth the tribonacci constant near `1.839`. That is
a tidy closed description of a natural tournament family, and it is the `n -> n-1` view. THM-337's
`x^3 - 3x^2 - x - 1` with growth `3.383` is the same family viewed `n -> n-2`, on the even subsequence,
and the two are related by exactly the relation stepping demands: `3.383 = 1.839^2`. Counting in steps of
two squares the eigenvalues.

That observation dissolved the puzzle I had left in the previous session. I had written the temperature
family as `x^3 - 3x^2 + a x - 1` and called the leading three "the corners," frozen. But here the same
family came out with leading coefficient one, not three. The resolution is that the leading coefficient,
the sum of the eigenvalues, is not frozen at all — it is `sum lambda_i^s`, and it changes when you change
the step `s`. The "three" was the Mode-B sum of squared eigenvalues; the "one" is the Mode-A sum. The
geometry I thought was frozen is really the step. So the framework was right that something is frozen and
something runs, but I had mislabeled which symmetric function is the invariant.

The actually-frozen invariant is the product of the eigenvalues, and the reason is the same stepping
argument that exposed the error. Under `n -> n-s` the product becomes `prod lambda_i^s =
(prod lambda_i)^s`, so if the product is `+-1` it stays `+-1` for every step. The sum moves, the middle
moves, but a unit product is a unit product no matter how you count. And across every C-finite family the
computation turned up — the all-upward count in both modes, the transitive count, the Pfaffian of the
even-tile family — the product of eigenvalues came out `+-1`, unimodular, every time. This is the
step-independent face of the Pfaffian. S713 said the Pfaffian is odd and `det(I+2A)` is its square; here
that parity reappears as the statement that the transfer eigenvalue-product is a unit, and it is a unit
because the vertex-deletion recursion is reversible over the integers — you can run it backward and stay
in `Z`, which is exactly the Pfaffian-cofactor kernel ladder. Unimodular transfer, reversible deletion,
and the parity ladder are one thing. That is real progress on the unimodular-transfer conjecture: not yet
a proof, but the mechanism is identified and the evidence is unanimous.

The frame also told the truth about its own limits, which is the other kind of progress. Of the seven
tile-rules I tried, only the most homogeneous — all up, all down — gave a low-order recurrence. The
gap-parity rule, the gap-mod-three rule, the threshold rule, the quadratic-residue rule: none were
order-five C-finite on twelve terms. C-finiteness is not generic; it is the privilege of families narrow
enough to have a small transfer, and most tile-rules are wider than they look. And the Pfaffian is wider
than the Hamiltonian count: the all-upward family has a tribonacci `H` but a Pfaffian sequence
`1, 7, 17, 23, 1, 89, ...` that is not low-order C-finite at all. So the Pfaffian, which I had been
treating as the clean parity companion, carries more boundary information than the path count does — a
caution for trying to compute `det(I+2A)` of a family by recurrence, and a pointer that the right width
for the Pfaffian is larger.

So the ledger now reads correctly. The product of the transfer eigenvalues is the frozen invariant, equal
to plus or minus one, step-independent, and it is the parity — the Pfaffian. The sum of the eigenvalues is
the step, the Mode, how far you walk per term, and it changes when you change your stride. The middle is
the temperature, additive when the dominant eigenvalue sits on the unit circle and the count stays
bounded, hot when it climbs off and the count explodes. One family has one spectrum; the recurrences are
its shadows at different strides; the parity is the part of the shadow that does not move when you change
the light. The next thing to track is not a new number but this: of a problem's transfer spectrum, the
product is the parity you cannot change, the sum is the stride you chose, and the middle is the
temperature that decides whether it is easy.
