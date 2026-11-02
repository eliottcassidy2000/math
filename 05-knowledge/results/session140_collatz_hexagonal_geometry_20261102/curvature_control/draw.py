#!/usr/bin/env python3
"""Small scientific figure; exact dyadic sources, sampled circles for display."""
from fractions import Fraction as F
from pathlib import Path
import math
import os

os.environ["MPLCONFIGDIR"] = str(Path(__file__).resolve().parent / ".mplconfig")
os.environ["XDG_CACHE_HOME"] = os.environ["MPLCONFIGDIR"]
Path(os.environ["XDG_CACHE_HOME"], "fontconfig").mkdir(parents=True, exist_ok=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({"font.size": 11, "svg.fonttype": "none"})
OUT = Path(__file__).resolve().parent
DEPTH = 4


def point(k, j):
    h = F(1, 2**k)
    return (j+F(1, 2))*h, h


fig = plt.figure(figsize=(10, 6), facecolor="white")
ax = fig.add_axes([0.07, 0.17, 0.53, 0.66])
ax.set_aspect("equal")
for k in range(DEPTH):
    for j in range(2**k):
        a, h = point(k, j)
        for s in (-1, 1):
            child_x, child_y = a+s*h/4, h/2
            center = a-s*11*h/8
            radius = math.sqrt(float(F(185, 64)*h*h))
            # The endpoints and center are exact rationals. Float conversion
            # and linear angular sampling are only for this illustration.
            start = math.atan2(float(h), float(a-center))
            end = math.atan2(float(child_y), float(child_x-center))
            angles = [start+(end-start)*i/64 for i in range(65)]
            xx = [float(center)+radius*math.cos(t) for t in angles]
            yy = [radius*math.sin(t) for t in angles]
            ax.plot(xx, yy, color="#2463a6", linewidth=1.6, zorder=2)

for k in range(DEPTH+1):
    xy = [point(k, j) for j in range(2**k)]
    ax.scatter([float(p[0]) for p in xy], [float(p[1]) for p in xy],
               s=30 if k < 3 else 18, color="#123957", zorder=3)

ax.axhline(0, color="#bd4f36", linestyle="--", linewidth=1.4)
ax.text(0.5, 1.035, "root", ha="center", color="#123957")
ax.set_xlim(-0.04, 1.04)
ax.set_ylim(-0.03, 1.09)
ax.set_xticks([0, .5, 1])
ax.set_yticks([0, .25, .5, 1])
ax.set_xlabel("x (Euclidean chart coordinate)")
ax.set_ylabel("y (Euclidean chart coordinate)")
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.spines["left"].set_color("#a8afb7")
ax.spines["bottom"].set_color("#a8afb7")
ax.tick_params(colors="#4c5968")

fig.text(.07, .93, "Binary branching with bounded hyperbolic routes", fontsize=18,
         weight="bold", color="#123957")
fig.text(.07, .875, "An abstract tree control — not a Collatz graph or a hexagonal mesh",
         fontsize=11, color="#596575")
fig.text(.64, .77, r"$ds^2=(dx^2+dy^2)/y^2$", fontsize=16, color="#123957")
fig.text(.64, .69, r"$z_{k,j}=((j+\frac{1}{2})2^{-k},\;2^{-k})$", fontsize=14)
fig.text(.64, .59, "Every blue edge has length", fontsize=12)
fig.text(.64, .54, r"$L=\operatorname{arcosh}(21/16)$", fontsize=15, color="#2463a6")
fig.text(.64, .45, "Distinct vertices are at least", fontsize=12)
fig.text(.64, .40, r"$\log 2$ apart in the hyperbolic metric.", fontsize=12)
fig.text(.64, .32, "Levels 0–4: 31 vertices, 30 edges.\n"
         "The Euclidean view crowds near y = 0;\n"
         "that boundary is outside the space.", linespacing=1.5,
         fontsize=11, color="#596575", va="top")
fig.text(.07, .065, "Dashed line: excluded boundary y = 0. Circle-arc routes do not cross. "
         "The tree metric is not preserved quasi-isometrically.", fontsize=10, color="#596575")
fig.savefig(OUT / "hyperbolic_binary_tree.svg", metadata={"Description":
    "Exact coordinate and circle-center source is draw.py; all-level proof is PROOF.md."})
fig.savefig(OUT / "hyperbolic_binary_tree.png", dpi=110)
plt.close(fig)
