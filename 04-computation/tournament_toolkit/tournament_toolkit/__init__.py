"""
tournament-toolkit: Mathematical tools for ranking, fraud detection, and AI diagnostics.

Built on tournament theory — the mathematics of pairwise comparison.

Core tools:
  FormalRank       - One-pass ranking from pairwise comparisons (O(m), no iteration)
  CycleDetector    - Streaming intransitivity/fraud detection
  CartanProbe      - AI confidence diagnostics via Cartan decomposition
  SpectralAnalyzer - Tournament quality and ambiguity measurement

The formal group F(x,y) = (x+y)/(1+xy) unifies all tools:
  - FormalRank uses arctanh(evidence) for additive aggregation
  - CycleDetector uses F(x,-x) = 0 for cancellation detection
  - CartanProbe decomposes matrices into competition + cooperation
  - SpectralAnalyzer uses eigenvalue flatness for quality assessment

Zero parameters. Zero training. Zero calibration.
Works on any pairwise comparison data.

Quick start:
    from tournament_toolkit import FormalRank

    ranker = FormalRank()
    ranker.add("Alice", "Bob")
    ranker.add("Bob", "Carol")
    ranker.add("Carol", "Alice")  # A cycle!
    print(ranker.summary())
"""

__version__ = "0.1.0"

from tournament_toolkit.ranking import FormalRank
from tournament_toolkit.cycles import CycleDetector
from tournament_toolkit.cartan import CartanProbe
from tournament_toolkit.spectral import SpectralAnalyzer
from tournament_toolkit.quaternion import QuaternionLinear, QuaternionAttention
from tournament_toolkit.shadow_cache import ShadowCache
