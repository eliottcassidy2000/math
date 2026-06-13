# tournament-toolkit

Mathematical tools for ranking, fraud detection, and AI diagnostics — all based on tournament theory.

## Install

```bash
pip install -e .
```

## Quick Start

### Rank anything from pairwise comparisons

```python
from tournament_toolkit import FormalRank

ranker = FormalRank()
ranker.add("Claude", "GPT-4o", wins=156, losses=134)
ranker.add("Claude", "Gemini", wins=178, losses=112)
ranker.add("GPT-4o", "Gemini", wins=162, losses=128)
print(ranker.summary())

# Find the most informative next comparison to run
print(ranker.most_informative_comparison())
```

### Detect fraud / intransitivity in real-time

```python
from tournament_toolkit import CycleDetector

detector = CycleDetector()
detector.observe("TraderA", "TraderB")  # A sells to B
detector.observe("TraderB", "TraderC")  # B sells to C
detector.observe("TraderC", "TraderA")  # C sells to A — cycle!

print(detector.summary())
print(detector.suspects())  # Who's in the wash trading ring?
```

### Diagnose AI confidence

```python
from tournament_toolkit import CartanProbe
import numpy as np

probe = CartanProbe()
attention = np.random.randn(8, 8)  # An attention matrix
result = probe.analyze(attention)
print(result['verdict'])           # CONFIDENT / UNCERTAIN / RISKY
print(result['tournament_fraction'])  # How decisive is this head?
print(result['entanglement'])      # How coupled are competition and self-knowledge?
```

### Analyze tournament quality

```python
from tournament_toolkit import SpectralAnalyzer

analyzer = SpectralAnalyzer()
# Paley tournament on 7 vertices
qr = {1, 2, 4}  # quadratic residues mod 7
adj = [[1 if (j-i) % 7 in qr else 0 for j in range(7)] for i in range(7)]
for i in range(7): adj[i][i] = 0
print(analyzer.quality_report(adj))
```

## The Mathematics

All tools are built on the **formal group law** F(x,y) = (x+y)/(1+xy):

- **FormalRank** uses arctanh (the formal group logarithm) to convert multiplicative evidence into additive scores
- **CycleDetector** uses F(x,-x) = 0 (evidence cancellation detects circular patterns)
- **CartanProbe** decomposes matrices into competition (antisymmetric) + cooperation (symmetric) via the Cartan decomposition of gl(n,R)
- **SpectralAnalyzer** uses eigenvalue flatness — Paley tournaments (flattest spectrum) maximize Hamiltonian paths

**Zero parameters. Zero training. Zero calibration.**

## Key Theorems Behind the Tools

- **Redei's theorem**: Every tournament has an odd number of Hamiltonian paths
- **OCF (Grinberg-Stanley)**: H(T) = I(Omega(T), 2) — path count = independence polynomial at fugacity 2
- **Spectral flatness principle**: Flat eigenvalue spectrum ↔ doubly regular ↔ maximum H
- **Cartan commutator formula**: ||[A,S]||^2 = n * S_2 / 2 (Landau irregularity controls sector coupling)
- **Ghost 13**: [competition, self-knowledge] = self-knowledge (the Cartan bracket relation)
