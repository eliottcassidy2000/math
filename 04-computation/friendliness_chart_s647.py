#!/usr/bin/env python3
"""
S647 — Friendliness instead of loneliness (the user's flip), as a chart.

LONELINESS (LRC, the usual view): the watched runner is LONELY at time t if it is
far (>= 1/n) from EVERY other runner. Lonely measure p0 = fraction of time lonely.

FRIENDLINESS (the flip): friendliness(t) = how many other runners the watched runner
is CLOSE to (clock-distance < 1/n) right now = the COVERING DEPTH (CoveringDepth.lean).
loneliness is just friendliness == 0. The whole covering-depth distribution {p_k} IS
the friendliness data; p0 (lonely) is only its bottom bar.

Renders a 2-panel SVG: (A) friendliness over one lap, (B) the friendliness distribution.
Pure python (no numpy/matplotlib) -> hand-written SVG.
"""
from fractions import Fraction

# ---- config: the n=14 "wall" (the arc's signature LRC case) ----
SPEEDS = [1,2,3,4,5,6,7,8,9,10,11,13,14]      # the 13 other runners
N = 14                                         # total runners; gap delta = 1/N
DELTA = 1.0 / N

def clock(x):                                  # distance to nearest integer in [0,1/2]
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def friendliness(t):                           # how many others are within the gap
    return sum(1 for v in SPEEDS if clock(v * t) < DELTA)

# ---- sample one lap finely ----
M = 200000
depth_hist = [0]*(len(SPEEDS)+1)
series = []                                    # (t, friendliness) subsampled for drawing
for i in range(M):
    t = i / M
    f = friendliness(t)
    depth_hist[f] += 1
    if i % (M//1000) == 0:
        series.append((t, f))
p = [c / M for c in depth_hist]
p0 = p[0]
mean_friend = sum(k*p[k] for k in range(len(p)))
maxk = max(k for k in range(len(p)) if p[k] > 0)

print(f"n={N}, speeds={SPEEDS}, gap=1/{N}")
print(f"LONELY (friendliness=0): p0 = {p0:.4f}  ({100*p0:.2f}% of the lap)")
print(f"FRIENDLY (>=1 neighbor):       {100*(1-p0):.2f}% of the lap")
print(f"mean friendliness = {mean_friend:.3f} neighbors;  max simultaneous = {maxk}")
print("distribution p_k (friendly to exactly k):")
for k in range(maxk+1):
    bar = '#'*int(round(p[k]*100))
    print(f"  k={k:2d}: {p[k]*100:5.2f}%  {bar}")

# ===================== SVG =====================
W, H = 1100, 520
def esc(s): return s
svg = []
svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="Helvetica,Arial,sans-serif">')
svg.append(f'<rect width="{W}" height="{H}" fill="#0f1424"/>')
svg.append(f'<text x="{W/2}" y="34" fill="#eaf0ff" font-size="22" font-weight="bold" text-anchor="middle">'
           f'Friendliness instead of loneliness &#8212; Lonely Runner, n=14 (the &#8220;wall&#8221;)</text>')
svg.append(f'<text x="{W/2}" y="58" fill="#9fb0d8" font-size="13" text-anchor="middle">'
           f'friendliness(t) = how many of the 13 runners you are within the gap 1/14 of '
           f'&#183; lonely = friendliness 0</text>')

# ---------- Panel A: friendliness over one lap ----------
ax0, ay0, aw, ah = 70, 90, 580, 360
svg.append(f'<rect x="{ax0}" y="{ay0}" width="{aw}" height="{ah}" fill="#161d33" stroke="#2b3759"/>')
svg.append(f'<text x="{ax0}" y="{ay0-10}" fill="#cfe" font-size="14">A &#183; one lap of time t &#8594;</text>')
ymax = maxk + 1
def AX(t): return ax0 + t*aw
def AY(f): return ay0 + ah - (f/ymax)*ah
# gridlines + y labels (friendliness levels)
for k in range(0, ymax+1, 1):
    yy = AY(k)
    svg.append(f'<line x1="{ax0}" y1="{yy:.1f}" x2="{ax0+aw}" y2="{yy:.1f}" stroke="#222c4a"/>')
    if k % 2 == 0:
        svg.append(f'<text x="{ax0-8}" y="{yy+4:.1f}" fill="#7e8cb5" font-size="11" text-anchor="end">{k}</text>')
svg.append(f'<text x="{ax0-46}" y="{ay0+ah/2}" fill="#9fb0d8" font-size="12" '
           f'transform="rotate(-90 {ax0-46} {ay0+ah/2})" text-anchor="middle">friendliness (neighbors)</text>')
# filled area (step)
pts = [f"{AX(0):.1f},{AY(0):.1f}"]
prev_t = 0.0
for (t, f) in series:
    pts.append(f"{AX(t):.1f},{AY(f):.1f}")
pts.append(f"{AX(1):.1f},{AY(series[-1][1]):.1f}")
pts.append(f"{AX(1):.1f},{AY(0):.1f}")
svg.append(f'<polygon points="{" ".join(pts)}" fill="#5aa9ff" fill-opacity="0.35" stroke="#5aa9ff" stroke-width="1"/>')
# highlight the lonely band (friendliness 0): a red strip at the baseline
svg.append(f'<line x1="{ax0}" y1="{AY(0):.1f}" x2="{ax0+aw}" y2="{AY(0):.1f}" stroke="#ff5a7a" stroke-width="2"/>')
svg.append(f'<text x="{ax0+aw-6}" y="{AY(0)-6:.1f}" fill="#ff8aa3" font-size="11" text-anchor="end">'
           f'touches 0 = LONELY ({100*p0:.1f}% of the lap)</text>')
svg.append(f'<text x="{ax0+aw/2}" y="{ay0+ah+22}" fill="#7e8cb5" font-size="11" text-anchor="middle">t (one full lap, 0 &#8594; 1)</text>')

# ---------- Panel B: friendliness distribution ----------
bx0, by0, bw, bh = 720, 90, 320, 360
svg.append(f'<rect x="{bx0}" y="{by0}" width="{bw}" height="{bh}" fill="#161d33" stroke="#2b3759"/>')
svg.append(f'<text x="{bx0}" y="{by0-10}" fill="#cfe" font-size="14">B &#183; how often friendly to k runners</text>')
nb = maxk + 1
pmax = max(p[:nb]) if nb else 1
gap = 6
bwid = (bw - gap*(nb+1)) / nb
for k in range(nb):
    frac = p[k]
    hh = (frac/pmax)*(bh-30)
    xx = bx0 + gap + k*(bwid+gap)
    yy = by0 + bh - hh
    color = "#ff5a7a" if k == 0 else "#5aa9ff"
    svg.append(f'<rect x="{xx:.1f}" y="{yy:.1f}" width="{bwid:.1f}" height="{hh:.1f}" fill="{color}" rx="2"/>')
    svg.append(f'<text x="{xx+bwid/2:.1f}" y="{by0+bh+13}" fill="#7e8cb5" font-size="10" text-anchor="middle">{k}</text>')
    if frac > 0.02:
        svg.append(f'<text x="{xx+bwid/2:.1f}" y="{yy-3:.1f}" fill="#cfe" font-size="9" text-anchor="middle">{frac*100:.0f}%</text>')
svg.append(f'<text x="{bx0+bw/2}" y="{by0+bh+30}" fill="#7e8cb5" font-size="11" text-anchor="middle">k = number of runners you are close to</text>')
svg.append(f'<text x="{bx0+gap+bwid/2:.1f}" y="{by0+18}" fill="#ff8aa3" font-size="10">&#8592; lonely</text>')

# ---------- caption ----------
svg.append(f'<text x="{W/2}" y="{H-12}" fill="#9fb0d8" font-size="12" text-anchor="middle">'
           f'You are FRIENDLY {100*(1-p0):.1f}% of the lap (mean {mean_friend:.1f} neighbors), '
           f'LONELY only {100*p0:.1f}%. Loneliness is the rare event; the bottom red bar is the LRC&#8217;s p&#8320;.</text>')
svg.append('</svg>')

with open('friendliness_chart_s647.svg','w') as fh:
    fh.write("\n".join(svg))
print("\nwrote friendliness_chart_s647.svg")
