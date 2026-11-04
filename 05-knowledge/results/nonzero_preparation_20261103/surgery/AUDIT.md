# Independent audit of the two actual H1429 rearmings

**PASS — analytic proof, complete source read, independent normal/optimized
replays, and independent named-graph reconstruction.** This audits the frozen
first-round surgery bundle. Producer status promotion and Git remain with the
parent. The conclusion is two actual next nonzero branches at the stated
source, not a reusable two-branch closure.

## 1. Actual objects and joint response

I independently checked the eleven-vertex input, physical root8, both literal
six-pair schedules, and the gap1 permanent birth. Every old physical name is
retained. The gap2 target has `(H,F,S,L,K,G)=(1429,1493,1431,469,1024,1493)`;
the gap1 target has `(1429,1923,1430,899,1024,1923)`. In both targets the
designated auxiliary has root8 as its sole outgoing neighbour.

For such an auxiliary a, the path from a to r must be exactly ar. In the
L-cover the a-starting path avoids r and is therefore the singleton a.
Consequently `L=F_r(T-a)` and `K=P_(T-{a,r})(2)`. Splitting an r-ending
two-slot cover according to the slot containing a proves `F=L+K`, including
the empty-slot cases: if a lies on the r-ending slot, delete the final ar;
otherwise delete the terminal a in the other slot. All remaining vertices
beat a, so the latter deletion is reversible. This establishes `G=F` from
actual paths, without inferring a joint response from H alone.

The inherited C2 interface then gives the literal children

|Digit|Order|H|F at retained root8|Defect|
|---|---:|---:|---:|---:|
|4|13|4291|4479|188|
|2|14|4289|5769|1480|

The preparations preserve H only as complete programs. Their F-values change,
and the gap1 preparatory birth itself has H6209 and F6939. The proof records
these changes and makes no individual-flip neutrality assertion.

## 2. Native schedules and an independent reconstruction

I read the complete `surgery_native.py` implementation, including pair-block
word construction, orientation reversal, helper cuts, both Hamilton
complements, named-set interpreter, helper deletion, observers, and invalid
guard controls. The compiled identity changes exactly the designated old
pair. A pair reversal is never admitted as a native primitive. Strongness is
explicitly checked at every helper birth and each complement, not merely
stored in the receipt. The root remains old name r at index8 throughout;
an appended helper is never silently substituted as root.

In addition to the producer replays, a separate verifier with no producer
imports reconstructed every native command from the receipt. Its interpreter
first complements all arcs and then restores the consecutive guard arcs,
independently of the producer's implementation. It checked exact guard
permutations and directions, all helper cuts, each intermediate adjacency,
helper-only deletion, and that each completed gadget changes precisely one
old pair. A separate forward endpoint recurrence recomputed H, F_r and the
root defect at all36 native states, and all completed stages and children.
Independent reachability checks confirmed every one of these36 states is
strong. This covers all24 complements and12 helper deletions in the two
fixed programs. The durable independent verifier `surgery_audit.py` returned207 always-active gates; its frozen stdout is `AUDIT.out`.

I also constructed the canonical funnel targets directly and verified the
following canonical-to-physical maps against the literal surgery targets:

```
gap2: [2,1,9,10,5,4,6,3,0,8,7]
gap1: [10,7,9,2,1,4,5,6,3,0,8,11]
```

They carry the canonical retained root to physical8 and the newborn to the
stated auxiliary. Thus the canonical positive witnesses and the physical
surgery programs refer to the same marked objects.

## 3. Unbounded formulas and their failure boundary

The all-clock adjacency proof is valid. Every flipped pair is outside the
TT_k prefix, and the permanent gap1 auxiliary also has uniform incidence to
it. The prefix remains a module. The contracted graph is the checked k1
stage, and its outside return paths connect prefix vertices in either order.
Thus all adjacency-stage expansions are strong for every k>=1. This does
not assert that the particular finite native guard programs have already
been lifted to every k; the primary keeps that distinction explicit.

I checked both extracted suffix graphs, port masks and seven-coordinate
signatures. Ordered-prefix assembly gives

```
gap2: H=744*2^k-59, F=744*2^k+5;
gap1: H=961*2^k-493, F=961*2^k+1.
```

After deleting the auxiliary, a Hamilton path ending at r must traverse the
entire prefix first, since the only outside vertex entering it is r itself.
The remaining fixed suffix ending count is2 or1. Together with the
sole-exit identity this pays both the start gap and `G=F` for every clock.

Subtracting the original `756*2^k-83` gives respectively
`-12(2^k-2)` and `205(2^k-2)`. Neutrality is therefore exactly k1. I also
checked the independent-clock equation
`31*2^(q-1)-63*2^t=-8`: equal exponents are impossible; distinct exponents
force the smaller to be3, leaving uniquely `(q,t)=(5,3)` in the declared
domain. Neither a finite coincidence nor an equality of scalar states is
being promoted to uniform source equality.

## 4. Complete finite-bank inverse

The34 entry alternatives form an explicitly declared bank, not a claimed
classification of every ready tournament. I read all source logic for core
response data, the tail automaton, target construction, short controls and
inverse search. The ready count is supplied by the proved
[polarized FUNNELS §§F3–F4](../../polarized_path_calculus_20261103/FUNNELS.md):

```
H_parent=L2*G-C+B1*h-D1*e,
F_parent=L2*p2-F,
s_parent=b.
```

An indegree-one birth has count `H_parent+F_parent-s_parent`, and converse
preserves H. These identities give the displayed all-word formula. I
recovered this supplier directly and requested the explicit link/equations;
the producer added them without changing either mathematical source or
receipt.

For an admissible nonempty length-n tail, `B1>=0`, `D1` is0 or1 and
`L2>=1+2^(n-1)`. Since `G+p2>0`, the displayed lower bound excludes that
length and all larger lengths once it exceeds the maximum requested target.
Thus the21,247-word enumeration is complete for all lengths in this bank.
Empty tails are handled exactly when the entry set is proper. The271 literal
short-word checks test the formula, actual gap and target strongness; the28
hits are actual constructed graphs. These finite checks support the proved
automaton formula and do not replace its unbounded argument.

The primary correctly leaves the next round uncompiled. In particular a
thirteen-name scalar target cannot receive the fourteen-name +2 child
without deleting an old name. This is a real lost predicate, not a harmless
relabeling. No closure or arbitrary ready-core classification follows here.

## 5. Reproduction and frozen pins

Run from the session worktree, with temporary output paths:

```
python3 -B 04-computation/nonzero_preparation_20261103/surgery_ready.py --output /tmp/nonzero_surgery_ready_audit.json
python3 -B -O 04-computation/nonzero_preparation_20261103/surgery_ready.py --output /tmp/nonzero_surgery_ready_audit_optimized.json
python3 -B 04-computation/nonzero_preparation_20261103/surgery_native.py --output /tmp/nonzero_surgery_native_audit.json
python3 -B -O 04-computation/nonzero_preparation_20261103/surgery_native.py --output /tmp/nonzero_surgery_native_audit_optimized.json
```

Both mode pairs reproduce the frozen receipts and stdout byte for byte:
ready has1,728 gates; native has430 gates. The independent207-gate verifier is also preserved and reproducible:

```
python3 -B 04-computation/nonzero_preparation_20261103/surgery_audit.py
```

Its stdout must match `AUDIT.out`. It reads the frozen native receipt and
imports no producer implementation. The source is4,336 bytes, SHA256
`128e0b033aed223ea7fc691b91cbb35dd011229bdb5e24d305c231f397863cbb`; stdout is103 bytes, SHA256
`18d264abac4e6bdc1cb31c0fb366081326573472e30af1a3d3017c08b8facc11`.

|Artifact|Bytes|SHA256|
|---|---:|---|
|surgery_ready.py|8410|024ac0b483d832865c6f12eabf0fe0f8676ac39188631e44910a0a4c9e66ca2f|
|surgery_native.py|16235|66eb9d682e0128bdf087c79b24a44dcfe127daf810ec560de85daa9689e220e2|
|ready.json|15717|4cc5a76fb0bd3ddb1276b86b3118c3dd7f51849a609b5cbad576cd1fbe128f44|
|native.json|42622|5d3af77198b041ab4280ad250252fbd153ea6853869ce2b5caba87d4511680d1|
|ready.out|191|6f6312f69ee3ee9944c39825aadcd98cf96e8b18423c576021561e2d5c3fc2a4|
|native.out|191|a200982060d1b2a41efbb91692617b05413107195828bbbeafff4a03ee916252|
|PROOF.md, audited prepromotion version|12294|6d3ca284a31a1d38cc324ba434cda8030e332898afd6b65d3ad423a511d3935e|

All pins match the producer manifest. No mathematical correction remains.
The accepted scope is two actual H-neutral preparations with full joint
root response and strong finite native implementation, plus the exact
all-clock response and unique-neutrality boundary. Unbounded reusable
nonzero branching remains open.
