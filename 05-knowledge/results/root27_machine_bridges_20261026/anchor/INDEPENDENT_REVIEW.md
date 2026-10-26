# Independent source, guard, tree and pump audit

2026-10-26. **PASS after the two explicit wording repairs below.** The
frozen README, complete checker, saved receipts, and inherited
`../../lca_ports_reverse_20261025/arithmetic/ANCESTRY.md` and arithmetic
README were read independently. No remaining mathematical gap was found.
Producer normal and optimized runs replay byte-for-byte with saved receipt
SHA-256 `90c7c8c4751871104028ab613b2d57f7627479586dd61bf36688c09db805baba`.
No producer file was edited by this reviewer.

## Repairs and hostile boundaries

The original wording “V is strictly increasing” would ordinarily assert
order preservation, which is false: V(5)=27>V(7)=11. Both root and this
reviewer independently found that witness. The producer now states the
needed property exactly: V(z)>z for every positive odd z. The inverse
termination argument uses this pointwise growth, never order preservation.

The cited physical parameter sector has b>=0, hence z=b+3>=3. The producer
now explicitly labels z=1,2 as arithmetic extensions. Its blank trajectory
and roots 7/39/27 stay inside the cited sector. This prevents valid extended
arithmetic controls from becoming an unsupported native-machine claim.
The machine association remains CITED; neither this audit nor the producer
builds a native Turing-machine formalization or identifies its transition
clock with the arithmetic macro clock.

## Proof audit

**Two sections, one guarded process.** The literal conjugacy is between
the full (a,b) macro process and `(q,z)=(a+1,b+3)`, not between its two
odd-integer root functions. On the physical blank run, z=7 starts at
macro time zero; x=27 is the first odd-x state at macro time three, with
a=6. It corresponds to even z=26, outside the odd-z section. Thus no
subtraction-by-one identification of root_- and root_+ is being used.

From odd z, `3z+1` has valuation ell; its first macro is a carry 1,
followed by exactly ell-1 carry zeros before the next odd z. The return
credit is q+3-ell. The first carry adds two slots; later zeros consume
one each. Therefore the return is reached precisely when ell<=q+2.
If ell>=q+3, exactly q+2 macros execute before the next zero is blocked.
This includes the equality boundary: a *formal* return with q'=0 was
not reached by executing the blocked symbol. This is checked independently
against the original (a,b) rules.

The old ancestry theorem identifies (A) with guarded survival from (6,27).
The certified three-macro prefix identifies that survival with survival
from blank (0,4). The new section identifies the latter with positivity
of q_J=1+3J-K_J at every odd-z return, equivalently K_J<=3J for every J.
If the process halts within a return, the first unexecutable formal return
already violates this inequality. Later unguarded states do not repair
that failure. This proves the claimed equivalence without swapping the
old and new meanings of J,K or asserting a native transition clock.

**Tree language and carry.** The slot recursion exactly recognizes ordered
full ternary forests. The last leaf is blocked; a realizing word of length
n records n-1 executed macros. Component-tree completion is not global
halt while other slots remain. The affine carry congruence determines every
prefix parity inductively, so its converse as well as forward direction
is valid. Carry concatenation preserves ordered composition. Rotation
changes the cylinder residue; it is not a permitted source quotient.
For source 7 all Q values are positive, so every displayed valuation is
well defined and finite.

**Leaf update and first mismatch.** Expanding the three affine factors
independently gives the coefficient 16-54=-38 in (L1). Its correction
factor is odd. Unequal valuations of the two summands therefore give
exactly min(p,h), while equal valuations give at least h+1. The only
improving position must equal the first mismatch and must actually be a
leaf. The canonical all-leaf suffix guarantees that latter condition;
it is not claimed for every arbitrary tree. Expanding the active leaf
gives exactly `2^h 3^(q+2)(3z+1)`, so its gained precision is ell and
remaining frontier length is q+3-ell. Nonpositive remainder means a
matching finite forest and a halt certificate, not an additional live
state. This is exactly the guarded odd-z return, not extra search freedom.

**Pumped criterion and quantifiers.** The block P is nonempty, begins in
1, has length L a positive multiple of three, and has balanced credit
with nonnegative prefix changes. Hence P^e0 is a tree for every e>=0.
Its carry C is odd; A-B is odd; N_z=(A-B)z+C is positive at positive z.
After removing the odd denominator, the two summands in the displayed
Q expression have valuations v2(N_z) and Le. If unequal their minimum
is too small; if equal the odd parts cancel once and meet the full
Le+1-bit requirement. This proves both directions of (P1), including
e=0 (the single leaf, realized exactly at even sources).

Thus the finite bank of 344 blocks gives 344 all-exponent exclusions by
the proved criterion. The bank itself remains bounded. There is no
interchange of “every repeat of each listed block” with “every tree”.
The explicit noncomb source 15 halts despite excluding every right comb,
so even that positive finite-family result does not justify universal
survival. Likewise the local two-prime charge is actually refuted on the
marked ray, not merely on an unrelated starting value. The Haar/Kraft
calculation is a cylinder measure, not a probability for any fixed integer.

## Independent executable controls

`review_check.py` imports no producer or earlier research code. Its main
independent route uses the original (a,b) rules directly, whereas the
producer mainly uses the translated arithmetic map. It checks:

| Control | Finite universe / count |
|---|---|
| Full translated macro trajectories | z=3..129, initial a=0..5, at most 40 macros; 23,023 executed comparisons |
| Guarded odd-z returns | odd z=3..513, q=1..9; 2,304 cases |
| Both blank cross-section clocks | 32 checkpoints in each section, 62 return edges |
| First mismatch versus literal bits | all binary words of lengths 1..7 at sources 1,2,7,15,39; 1,270 cases |
| Arbitrary-word leaf expansions | all zero positions in those words; 3,845 exact identities and valuation cases |
| Pumped iff against literal itineraries | independent ballot-word enumeration through four internal vertices; 71 blocks, sources 1..64, repeats 0..5; 27,264 cases |

The source-1 and source-2 rows are explicitly arithmetic extensions.
The pump test independently enumerates words by positions of their 1s and
the ballot predicate, rather than importing the recursive producer tree
generator. The all-order conclusions come from the proofs, not these caps.

Run from the repository root:

    python -B 05-knowledge/results/root27_machine_bridges_20261026/anchor/review_check.py
    python -O -B 05-knowledge/results/root27_machine_bridges_20261026/anchor/review_check.py

Both modes pass identically. Saved `review-normal.json` and
`review-optimized.json` have SHA-256
`70b0cda47a78371a9ab12caa2a377e8d8b29cd40fc145e327a675144666f4822`.
The root-7 tree exclusion, root-27 debt inequality and blank-state
nonhalting remain OPEN.
