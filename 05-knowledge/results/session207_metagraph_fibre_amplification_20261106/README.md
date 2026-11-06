# Fixed tournament fibres can have arbitrarily large neighbours

For the actual path-complement metagraph, not a single-arc-flip graph, an explicit six-vertex padding keeps one endpoint's ordinary Hamilton-path count and its original/converse-merged tiling-fibre degree constant. The opposite endpoint's count grows at least fourfold per step, and its automorphism group cannot grow. Thus its fibre degree is unbounded too.

At order n+6r the constructed pair satisfies H(T_r)=H(T) and H(U_r)>=4^r H(U)+(4^r-1)/3 for r>=1. This holds starting with any tournament and chosen Hamilton path. The low-count endpoint is nonstrong and the opposite endpoint is strong.

The full proof and exact automorphism identity are in RESULT.md. This does not classify all attainable counts, prove general metagraph connectivity, or bound degrees after parallel edges are discarded.

Reproduce from this directory:

    python3 check.py > replay.json
    cmp replay.json RESULT.json
    python3 -O check.py > replay-opt.json
    cmp replay-opt.json RESULT.json
    python3 audit/check_padding.py
    sha256sum -c SHA256SUMS

The main verifier has explicit checks that survive Python optimization. It covers all fixed-path tilings of orders1..4 with three padding lengths (36 count checks), plus orders1..5 with six-vertex padding (76 automorphism checks). The independent checker uses a separate endpoint recurrence for the36 counting controls. It uses ordinary Python assertions and must run without -O. The all-order statements rely on the audited proof, not these finite tests.

The independent mathematical review is in audit/REVIEW.md. Its checker was relocated into this package by changing only the input path to the package-relative RESULT.json; mathematics and checks were unchanged. Historical paths in the proof/review are source provenance, not replay dependencies. No external dependency beyond Python3 is needed.
