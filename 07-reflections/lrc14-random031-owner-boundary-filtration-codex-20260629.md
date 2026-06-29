# LRC14 Random031 Owner-Boundary Filtration Reflection

This continues HYP-3520's owner-cobordism certificate and HYP-3521's terminal
ledger by turning the remaining owner word into a sharper branch-boundary
filtration.

The useful correction was simple: stop asking the bypass to transport all seven
seam owners.  The bypass never promised that.  It promised a stable transport
word.

The exact HYP-3522 audit makes the random031 boundary look like a filtration:

```text
seam=(23,45,93,113,147,169,173)
transport=(23,93,113)
branch_boundary=(23,93,147,169)
residual=(45,173)
```

This is better than the older "seven-owner gluing clause" because each layer
has a different geometric job.  `(23,93,113)` is the word carried by all
twelve pure-bypass cells and by all six mirror pairs.  `(147,169)` is not
transport; it is supplied by the ordinary branch-boundary cells next to the
bypass.  `(45,173)` is the true residue: one dead-island boundary owner and
one apex-boundary owner.

That suggests the next proof target should be a residual lemma, not another
search for a bigger gluing gadget:

```text
transport constancy
+ branch-boundary bracket lift
+ HYP-3510 connected phase carrier
+ HYP-3511 free-hole bracket carrier
=> residual (45,173) cannot support a counterexample.
```

The incoming HYP-3513 firewall Nerode audit makes the attachment more precise.
Private-firewall status can be read through existing colored axes, but the
full random031 route still needs route sidecar `R` unless a reconstruction
lemma is proved.  So the residual lemma is really:

```text
owner filtration + route sidecar R + HYP-3520 safe seam-sheaf quotient
=> residual (45,173) vanishes.
```

The "owner-boundary persistence" phrase also becomes less mystical.  It is not
that labels literally move through every path.  It is that the legal quotient
must remember which labels are transport, which are boundary brackets, and
which are residual debt.  A quotient may forget a gate only after it has paid
for that forgetting with the corresponding owner sidecar.

Two concrete proof experiments now look high leverage:

1. Prove a finite branch-order bracket lemma for the two bypass intervals:
   adjacent ordinary cells supply exactly `(147,169)` and do not alter the
   transport word.
2. Build a two-owner residual obstruction table for `(45,173)` against the
   HYP-3490 private-label firewall and the four dead-island/apex boundary
   positions.

Creative extension: the owner filtration behaves like a tiny sheaf exactness
test over the mirror-punctured cylinder.  Transport is the local section over
the bypass.  Branch boundary is the Cech coboundary.  `(45,173)` is the
cohomology-looking residue.  The next computation should not overclaim that
language, but it is the right mental geometry: the proof is trying to show the
residue class vanishes once the rank-2 route and free-hole packets are attached.
