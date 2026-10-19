# Complete three-color obstruction classification

Read [the analytic proof and release results](CLASSIFICATION_AND_TWO_RELEASES.md).

Main results:

- In any representable rank-three matroid, every disconnected nonempty three-color independent-transversal graph is exactly a cyclic assignment of three independent directions. It has precisely two components.
- In a graceful fixed-edge fiber, this is exactly the tree-supported chordless-six-cycle obstruction.
- Releasing only bridges of the original fixed-plus-free-support graph gives a product and cannot remove an old obstruction.
- A nonbridge edge whose color crosses the failed color-1 bridge cut supplies a two-release certificate.
- An explicit 17-label arithmetic example remains disconnected after releasing colors 1 and 2 together, even though both are full-rank buffers. Another single edge repairs it, so existential one- and two-release statements remain open.

Dependency: [the preceding tight-kernel checkpoint](../graceful_three_color_release_20261019/THREE_COLOR_RELEASE_DICHOTOMY.md).

Exact replay, Python 3 with standard library only:

    python verify_classification_and_releases.py
    python -O verify_classification_and_releases.py

Both outputs agree byte for byte with the included receipt. The verifier checks the explicit arithmetic fiber, selected one- and two-edge releases, the bridge product counts, the tight-kernel certificate, and all 64 zero-diagonal binary matrices used in the graphic specialization. The general-field classification is proved analytically, not inferred from those 64 matrices. No search transcript or exploratory support census is included.
