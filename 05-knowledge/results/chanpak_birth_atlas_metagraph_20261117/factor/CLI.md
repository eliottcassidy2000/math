# Standalone exact first-birth CLI

Python3 standard library only. No network or repository imports.

```sh
python solve.py 9
python solve.py 21 --max-labelled 0
python solve.py 33 --max-labelled 10000
python solve.py 15 --atlas /path/to/births.tsv
python solve.py 21 --max-labelled 0 --output result.json
```

The default budget is1,000,000 labeled tournament tests. `--max-labelled 0` removes the budget; large inputs may be completely impractical. Factoring the input by trial divisors is also intentionally elementary, not claimed efficient. Budget exhaustion returns status`incomplete` and exit code2, never`unattainable`. A fully decided result exits0.

## Outputs

- `minimum`: exact minimum nonempty order and literal outgoing adjacency masks. The final witness is independently recounted.
- `unattainable`: the parity theorem excludes an even target, or every relevant strong search was exhausted through its exact Busch cap and no SCC factorization exists.
- `incomplete`: a budget stopped a search. Any accompanying witness is labeled only as an upper bound.

The certificate lists each target divisor, its Busch strong cap, fully exhausted orders, mathematically excluded orders(H<=n!), and any partial/witness search. It is a replayable exhaustive-search receipt, not a succinct independently checkable certificate of every omitted graph. Nonexistence depends on the cited Busch theorem and completeness of the enumerator.

## Exact caps, without floating point

For n>=3 the inherited Busch minimum is

    f(3q)=3*5^(q-1), f(3q+1)=5^q, f(3q+2)=9*5^(q-1).

`strong_cap(d)` is the largest n with f(n)<=d, evaluated using integers. Only strong targets dividing h are searched. A proper-factor witness or an atlas witness provides an upper bound U; strong search then stops below U, since a witness at U is already available and larger strong classes cannot improve the minimum. No canonicalization is needed for this existence/minimum task: all labeled orientations are enumerated, and one strong witness suffices.

The output also gives `Busch_factor_core_cap`, the maximum sum of strong caps across all divisor factorizations of h. This is the exact core bound implied by the Busch inequalities, not a claim that the bound itself is attained. Earlier repository work supplies the convenient coarser bound floor(3log_3 h); the implementation avoids approximate logarithms altogether.

## Optional atlas

Accepted formats:

```json
{"9": [14,28,24,16,1]}
```

or TSV columns `h`, `order`, `outgoing_masks`, with masks comma-separated. The existing `birth_program/births.tsv` format is supported. Multiple witnesses per h are allowed; the shortest is used. Only target divisors are loaded. Every relevant graph is checked to be a tournament and its H is recomputed. An atlas order is **not trusted as minimal**: the lower-order proof is still performed. A valid nonminimal witness remains a useful upper bound. Wrong-count or malformed relevant witnesses are rejected.

## Inheritance

The exact Busch formula is cited in `01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`, citing A.H.Busch, EJC13(2006), N3, Corollary1. The earlier `birth-atlas-now/PROOF.md` and `birth_program/README.md` already state first-birth factorization and logarithmic finite-core decidability. This CLI implements that reduction, adds safe budget/atlas handling and literal witnesses, and does not claim those earlier theorems as new.

Run `python test_cli.py` for exact minima at1,3,5,9,11,13,15, nonexistence at7,21, an even-target test, budget exhaustion, valid/minimal/nonminimal atlas controls, and rejection of a wrong-count atlas.
