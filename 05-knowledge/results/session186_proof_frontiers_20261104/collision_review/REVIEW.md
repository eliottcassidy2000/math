# Independent review: exact gate-edit response separation

**ACCEPT** the finite theorem and the three explicit prime scalar hostiles in `collision/PROOF.md`. No corrections requested. Source hashes are frozen in `inspected.sha256`.

## Independently rebuilt complete search

`check.cpp` is separately written and does not import author code. It builds each inflation directly by mapping vertex names to quotient seed names. Its enumeration uses a three-index representation with repeated indices encoding radius one and two, independently of the author's nested incremental flip/restore loop. Every target is rebuilt from the untouched base. Homogeneous pairs are checked by comparing every exterior vertex's *incoming* incidences to the candidate pair, rather than XOR of outgoing rows. Exact Hamilton counts use predecessor subset sums, rather than the author's forward propagation.

The full independently executed census agrees with all 27 published rows, entry by entry:

- 464,709 source/edit-set cases;
- 28,293 pair-free targets;
- exactly three scalar matches, with exactly the published rows;
- every radius-k universe size equals `binomial(binomial(n,2),k)` for each source.

`census.out` retains all rows and the entire match list. This is a complete finite search over the declared nine sources, not a sample. Source/edit pairs may refer to targets reached from other sources; the theorem does not require unique targets. Pair-free is only a necessary condition for prime, as the author explicitly states.

## Independent witness validation

`check.py` independently evaluates Hamilton counts by **uncached DFS enumerating every valid directed Hamilton word**. It does not use subset DP or the author's counting routines. It confirms:

- the source and all three targets have H=3489;
- the literal named attachment `a={0},c=empty`, with x->y, gives 17055 on the source and 15579 on each target;
- the three exact reversal sets match the receipt;
- all three targets are valid tournaments, strong by all-pairs Floyd reachability, and prime by checking every proper nontrivial vertex subset with an outgoing-incidence criterion.

The source/query counts imply linear coefficients 13566 and 12090 and difference 1476. The identity F(0,0)=H is valid: x then y must precede all old vertices. Thus the scalar mismatch excludes 28,290 targets and the actual singleton query excludes the remaining three. No sufficiency claim about H,S,E,C or a degree-three battery is needed.

The independent Python normal and optimized runs agree byte-for-byte. The author's `verify.py` was also replayed normally and with `-O`; both outputs agree with its frozen receipt.

## Proof and scope inspection

The finite proof is sound. The statements about fixed names, radius at most three, orders 11--13, and pair-free targets are all essential and retained. The small-universe coefficient search is a supplementary author control, not load-bearing in this finite neighborhood theorem. Its `search.cpp` and recorded output were inspected: it retains all bidegree-(2,2) coefficients and compares complete vectors. The order-six coefficient census was not independently recomputed in this review. Its degree-two diagnostic is explicitly bounded and does not justify truncation elsewhere.

The main C++ counts fit uint64: old graphs have at most 13 vertices, so H<=13!. Literal witness attachments also have 13 vertices. No floating-point, probabilistic equality, unsafe coefficient truncation, or relabelling occurs.

No universal reconstruction, four-edge-neighborhood, hidden-module-converse, or literature-novelty conclusion follows. The packet's stated boundaries accurately retain these limitations.

## Replay

    g++ -std=c++17 -O3 -Wall -Wextra -Werror check.cpp -o check
    ./check > census.out
    python check.py > normal.json
    python -O check.py > optimized.json

The checker uses explicit exceptions, not removable assertions. The checker resolves its author packet through the sibling `collision` directory and its census through its own directory, so the packet is portable.
