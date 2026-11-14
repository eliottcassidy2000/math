# Changing physical sources and readers for Chan–Pak

Continues the [67/69 source–reader pair](../chanpak_pair_67_69_20261114/README.md)
and [single-count cover readout](../chanpak_pair_cover_readout_20261114/PROOF.md).
This session changes actual marked tournaments. It does not extend a numerical
coverage threshold. Chan–Pak and universal absolute placement remain open.

## What changed

| Result | Exact scope |
|---|---|
| [Four-odd two-height window and neutral return](PRODUCT_RETURN.md) | For every k,m≥1 the four counts are 73(2^k−1)(2^(m+1)−1)+{200,202,204,206}. Exchanging the Mersenne factors preserves the same absolute count, order and both full path endpoint counts at the retained root, but changes its weighted continuation budget. A literal next seam detects the change. |
| [Physical reader gauges](READER_GAUGE.md) | One-axis windows of lengths ten, nine and seven, including a seven-odd window that changes the source class as well as the reader. The maximum consecutive width in a fixed normalized reader fiber is exactly three through order seven and four through order eight; not a global bound. |
| [Cyclic source switch and determinant compensation](SOURCE_SWITCH.md) | A smaller two-height triple starts at 97/99/101. Changing the reader and swapping two actual nontransitive blocks gives 88021/88023; neither numerical output is presented as a new spectrum record. A separate compensation pair 227125165/227125167 has two strong outputs, recounted by standalone C++ subset DP. With one block fixed, this exact compensation cannot grow the other block without bound. |
| [Apex-address reader change](APEX_READER.md) | Exact inverse on guarded rays and new physical pair/triple families; holding H+P(2) fixed allows a sharp gap of four but not two under the same reader. |
| [Factorial lock](FACTORIAL_LOCK.md) and [one-contact-stage bypass](SELECTIVE_BYPASS.md) | Matching the first r evaluations forces a difference divisible by 2r!(2r+1)! in every common ordinal/singleton-cycle context, with an extra r+1 for even source order. The r=2 bound is physically sharp in both parities. Two selective contacts from a new vertex break the restriction: 8931/8933 and 37689/37691 are literal examples even though the unmarked source H,P(2) agree. |
| [Bounded-source sieve](BOUNDED_SOURCE_SIEVE.md) and [composite subgrid](COMPOSITE_SUBGRID.md) | The apex interface with bounded source H misses a positive density of odds even with unbounded source order and freely chosen reader addresses. Any finite polynomial window bank on a dyadic grid also has a common parameter subgrid with no primes or 7p outputs past its corner. |

## A second reason both source cores must change

Use the notation of [the source-swap identity](SOURCE_SWITCH.md).
For two equal-(a,s) readers let C=s−a>0, let K=h_A h_B,
and D=Q_Ah_B−Q_Bh_A. A gap of ±2 requires

    C D + (c−c') K = ±2.

Each of gcd(h_A,Q_A), gcd(h_B,Q_B) and gcd(h_A,h_B) divides both
terms, and is odd, so must equal 1. If both source orders are at
least two then 4 divides D, and c−c' must be 2 modulo 4.

More strongly, if K>1 and both h_A,h_B are bounded, **both source
orders are bounded for these unit-gap swaps**, even with changing
readers of unrestricted order. As in the sieve, under bounded H each
Q/h is on one of finitely many rays α·2^u prescribed by a finite
strong-component word and an arbitrary number u of singleton pads.
The swap inequality forces the two ray values to be within 4.
They cannot be equal: D=0 would make the gap an even multiple of K.
For any pair of slopes, factor out 2^min(u,v). The exponent difference
can take only finitely many values when the distance is at most 4.
Among those values the nonzero difference of normalized slopes has
a positive minimum, so min(u,v), and then max(u,v), is bounded.
This closes the bounded-H source-padding escape for this swap route.
It is independent of the apex sieve and says nothing about other
reader geometries.

## Why the old gadgets helped

The planar Jacobian response fibers suggested keeping the reader
profile fixed while changing a marked realization. Lonely Runner's
two clock carriers suggested matching actual factors by exchanging
the clocks, not merely matching their periods; the retained-root
budget shows exactly which coordinate that exchange fails to preserve.
The June Green-conductance note supplied the sidecar heuristic, while
the August Rule 30 and factorial–Jacobian notes suggested separating
observability from divisibility. The October Fermat/repunit analysis
supplied the right way to group the unbounded singleton directions
before constructing an absolute-target sieve. The linked notes in
each proof retain the original scopes; these are transfers of a method
or a question, not claims to solve the unrelated conjectures.

## Reproduce

Run from this directory:

    python replay_source_switch.py
    python replay_apex_reader.py
    python literal_windows.py
    python replay_product_return.py
    python replay_factorial.py
    python replay_sieve.py
    g++ -O3 canonical_readers.cpp -o /tmp/chanpak-reader-census
    /tmp/chanpak-reader-census 8 > /tmp/chanpak-reader-census.txt
    python reader_census.py /tmp/chanpak-reader-census.txt

Optional strong-pair check (the 22-vertex subset DP needs roughly 370 MB):

    g++ -std=c++20 -O3 strong_recount.cpp -o /tmp/chanpak-strong-recount
    /tmp/chanpak-strong-recount

The all-height statements follow from the physical composition equations,
not from the numerical tests. The remaining constructive target is an
unbounded rule choosing a source and marked reader from the prescribed
integer, across profile fibers. A promising next question is to retain
the exact neutral return while varying the *selective* endpoint cut: the
return changes a real port budget, while the cut is the one tested
operation that changes a locked profile by the required unit step.
