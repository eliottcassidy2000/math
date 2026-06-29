# LRC14 Minimal Bad-Core Cover Extractor

HYP-3436 is the inverse view of HYP-3435.  HYP-3435 says every checked covering
row has a surviving two-adic branch cell.  HYP-3436 asks what the failed cells
look like locally.  The object is no longer a scalar branch-union measure but a
finite interval cover:

```text
E_safe cap B0_odd cap B1_odd.
```

That change matters.  A counterexample to the branch lemma must not merely make
the good measure small.  It must concatenate local two-color bad covers across
every component of `E_safe`.  The new extractor records those local covers by
minimal odd-owner sets in each branch.

The stress result is encouraging: `135/135` rows still have a survivor, and
the measure identity with HYP-3435 holds `135/135`.  But the stronger signal is
structural.  Across `11670` bad-core intervals, `10288` are `(1,1)` covers:
one odd owner covers branch `0`, one odd owner covers branch `1`.  The maximum
minimal two-color owner count is only `6`.  Endpoint support is almost always
one labelled wall on each side: `(1,1)` for `11634/11670` bad-core components.

This makes the next proof target feel much less foggy.  We do not need to
prove a uniform scalar density floor first.  We can try to prove that these
small local owner covers cannot tile all even-safe components without emitting
a named debt.  HYP-3434's overlap tax becomes the one-branch shadow of the same
phenomenon: local owner overlaps rescue deficits that raw sums think are fatal.
HYP-3429's endpoint spines are the compressed survivor side.  HYP-3427's wall
atlas is the alphabet.  HYP-3423 says which quotients are not allowed to forget
the branch/magnitude data.

The canonical family gives the cleanest new clue.  For

```text
{1,2,3,4,5,6,7,8,9,10,11,13,84m}
```

through `m=30`, every bad-core component is singleton/singleton from `m=3`
onward.  That is the local-cover form of the HYP-3431 corridor-fence theorem:
the moving high even wall slices the obstruction into small cells, but the
survivor gaps remain.  The first two cases have pair-owner overlaps; from
`m=3`, the bad cells are simple enough that an all-`m` singleton-tail lemma
looks plausible.

Candidate exits are HYP-3431 corridor-fence structure, HYP-3429 endpoint
spines, HYP-3428 two-adic loss debt, HYP-3427 wall words, HYP-3426 mirror
reduction, HYP-3417/HYP-3420 owner-current, HYP-3423 topology legality, and
HYP-3421/HYP-3129 signed-SPEC/Rprime.

After rebasing over HYP-3437, the handoff is clearer: HYP-3436 supplies the
bad-core atoms and owner-cover signatures, while HYP-3437 should build the
incidence/Menger graph that explains negative-slack overlap tax and multi-owner
no-gluing.

The most useful next agent task is therefore not another random stress bank.
It is a labelled packet theorem:

```text
local minimal bad-cover packets cannot concatenate into a global cover of E_safe
unless owner-current, endpoint-spine, overlap-tax, or two-adic descent debt is
explicitly emitted.
```

There are two immediate proof directions.

1. Prove the canonical singleton-tail law for all `m >= 3`, as a local-cover
   strengthening of HYP-3431.
2. For general rows, build a graph whose vertices are `E_safe` components and
   whose edges are shared odd-owner cover tokens.  Try to show every connected
   all-covered component chain either has a dischargeable endpoint spine or a
   two-adic child row with a survivor.

Tournament Analysis helped because the right vertices were not runners.  The
selected vertices were proof obligations and cover certificates; alternate
vertices considered were safe components, bad-core intervals, endpoint events,
odd-owner tokens, Fourier modes, and runners.  This quotient preserves the
exact two-branch survivor predicate and forgets only raw runner ordering plus
open-endpoint conventions.  That is a legal forgetting move.
