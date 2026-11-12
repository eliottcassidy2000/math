# A000568: adaptive exact backbones and finite moment jets

[Proof and precise scope](PROOF.md).

## Main results

- **Proved:** condition on a finite set of exactly handled odd cycle lengths,
  then compute the envelope by a size-only recurrence. A finite triangular
  hierarchy of marked multiplicity moments computes every fixed correction
  order without enumerating unrestricted partitions.
- **Implemented:** C0,C1 first jet through96 and C0,C1,C2 second jet through40,
  controlled independently against partition sums. General fixed-order closure
  is proved, not an optimized arbitrary-order implementation.
- **Certified interval:** at n64, adding exact backbone lengths changes the
  sufficient residue budget from738 bits ({3}) to493 ({3,5}),310 ({3,5,7}),187
  ({3,5,7,9}),121 ({3,5,7,9,11}). The table uses partition-backed degree-four
  coefficients. Sibling ground-state corrections improve these further.
- **Exactness boundary:** for a prefix backbone whose first residual length is
  a>=5, S3 is exact below4a+4; the two-sided C4 interval is exact through4a+5.
  Fixed backbone and fixed correction depth still eventually lose absolute
  integer accuracy.

## Reproduction

From any working directory:

    python3 /path/to/enumeration/backbone.py > results.json
    python3 /path/to/enumeration/second_jet.py > second_jet_results.json
    python3 /path/to/enumeration/benchmark.py > benchmark.json
    python3 /path/to/enumeration/integer_jet.py > integer_jet_results.json

`backbone.py` imports Session278's partition generator and elementary-symmetric
routine by a path anchored to its own file, not the current working directory.
All certificates are rational; timings are descriptive wall-clock samples.
The `dp_results.json` supplementary first-jet controls are redundant with the
main controls/benchmark and contain exact C0/C1 outputs through64.

## Matched-workload honesty

For B={3,5}, both methods compute exactly C0 and C1. In the recorded run:

| n | marked DP seconds | direct partition seconds |
|---:|---:|---:|
|32|0.060|0.006|
|48|0.495|0.076|
|64|2.242|0.613|
|80|7.502|3.958|

The initial rational DP is **slower at every tested size**. Compiling factorial
and dyadic denominators into integer transitions changes this:

| n | integer jet seconds | partition C0/C1 seconds |
|---:|---:|---:|
|32|0.00695|0.00609|
|48|0.0437|0.0684|
|64|0.1748|0.5343|
|80|0.5561|3.4949|

At96 the integer jet agrees with the rational jet and takes1.767s versus16.470s.
These local measurements concern exactly C0/C1, not full A000568 enumeration;
they are not compared to the optimized GMP counter. The polynomial operation
exponent worsens as correction depth/backbone cardinality grows. General
fixed-order closure is a theorem; arbitrary-order optimized software is not
claimed.

## Structural connection

The finite state is not just an approximate count: it includes precisely the
multiplicity moments required by future interactions. This is the enumeration
analogue of retaining endpoint response in a tournament construction. Neither
an ordinary-class total nor an accurate scalar forecast produces a physical
Chan–Pak target selector on its own.
