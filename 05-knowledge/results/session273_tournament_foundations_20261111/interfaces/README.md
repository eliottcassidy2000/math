# Integer gauge packets

[Proof and scope](PROOF.md) · [structural controls](RESULTS.json)

Run from the repository root:

```sh
python3 05-knowledge/results/session273_tournament_foundations_20261111/interfaces/check.py
```

Main results:
- A finite graph-coboundary certificate classifies positive multiplicative local changes that preserve every closed-word count.
- An explicit physical packet at order N(10+3K) has2+(K+1)(2^N−2) ordinary tournament classes, exactly floor(N²/4)+2 H values, and exact fiber/automorphism formulas.
- It preserves the complete unmarked path-cover profile, not just H. A single rooted seam can expose the hidden placement: the26-vertex equal-profile pair with H131565 yields27-vertex counts131565 and219275.

Sources: Session268's interacting ten-vertex cells and seam law; THM-3166's order-join path-colour transform; the old lonely-runner signed-telescoping reflection. The graph-cohomology principle is standard; the explicit tournament lift is the new application. No universal coverage, new Chan–Pak range, or reduction of another open conjecture is claimed.

Controls use actual adjacency matrices and independent subset DP on eight13-vertex cells and six endpoint probes; complete small word families use the proved SCC encoding rather than pretending to run whole-graph DP at larger orders. The all-order assertions are established in the proof, not extrapolated from the controls.
