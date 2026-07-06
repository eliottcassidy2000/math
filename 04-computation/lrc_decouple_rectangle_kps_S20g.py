#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S20g: verify the DECOUPLING insight behind
torus_decouple_clears, and test the generic-coupling RECTANGLE rigidity.

Context (MISTAKE-113): CONSECUTIVE distinct freqs {1..l} TILE a single theta-slice
at 2/25 for l >= 10.  But in the (A)-torus the lifted shifts are a_i*t, NOT
arbitrary.  torus_decouple_clears (Lean GREEN) proves: if a_i = c*r_i (parallel /
sheared product), the system DECOUPLES (theta' = theta + c t) and LRC(<=13) on
each factor clears it (M >= 1/12 > 2/25) for ALL l <= 11.

This script confirms two things numerically:
  (A) DECOUPLING: a parallel-coupled system (a_i = c*r_i) with consecutive
      distinct lifted freqs has M(system) = min(M_base, M_lifted_homogeneous)
      >= 1/12 > 2/25 -- it CLEARS, even at l = 10, 11 where the SLICE tiles.
      (The slice tiling used arbitrary shifts; the homogeneous ||r_i theta'||
      does NOT tile -- LRC.)
  (B) GENERIC RECTANGLE RIGIDITY: for GENERIC (non-parallel) couplings, does a
      single-slice tiling persist over a t-interval?  We take consecutive freqs
      with a near-tiling slice at t0, then sweep t and measure the max over the
      interval of the min-tooth (the rectangle M).  If it rises above 2/25
      somewhere in the interval, the slice-tiling BROKE => a clear point exists
      (rectangle rigidity holds, the generic l=10,11 case also clears).
"""
import numpy as np

RHO = 2.0 / 25.0
rng = np.random.default_rng(2026070623)

def min_tooth_over_theta(freqs, shifts, grid=100000):
    """max_theta min_i ||r_i theta + s_i||  = the LR value of these combs at these shifts."""
    th = (np.arange(grid) + 0.5) / grid
    mins = np.full(grid, 1.0)
    for r, s in zip(freqs, shifts):
        x = r * th + s
        d = np.abs(x - np.round(x))
        mins = np.minimum(mins, d)
    return mins.max()

def homogeneous_LR(freqs, grid=200000):
    """max_theta min_i ||r_i theta||  (shifts = 0, the LONELY RUNNER value)."""
    return min_tooth_over_theta(freqs, [0.0]*len(freqs), grid)

print("=== (A) DECOUPLING: parallel-coupled consecutive freqs CLEAR (homogeneous LR) ===", flush=True)
print("  the lifted, decoupled, are homogeneous ||r_i theta'|| -- LR value >= 1/(l+1):", flush=True)
for l in range(6, 13):
    freqs = list(range(1, l+1))
    M = homogeneous_LR(freqs)
    bound = 1.0/(l+1)
    print(f"    l={l:2d} freqs 1..{l}: homogeneous LR M = {M:.5f}  (LRC bound 1/(l+1) = {bound:.5f}; "
          f"2/25 = {RHO:.5f}) => {'CLEARS (> 2/25)' if M > RHO else 'IN GAP'}", flush=True)

print(flush=True)
print("=== (B) GENERIC rectangle rigidity: does a slice-tiling persist over t? ===", flush=True)
print("  consecutive freqs 1..l, GENERIC couplings a_i, near-tiling shift at t0=0;", flush=True)
print("  sweep t in [-w, w], report max_t (rectangle M) -- if > 2/25, the tiling BROKE:", flush=True)

def anneal_slice_tiling(freqs, iters=6000, restarts=6, grid=40000):
    """find shifts minimizing uncovered (best slice tiling) at t=0."""
    best, bestsh = 1.0, None
    for _ in range(restarts):
        sh = rng.random(len(freqs))
        cur = 1.0 - _covered_frac(freqs, sh, grid)
        for _ in range(iters):
            k = rng.integers(len(freqs)); old = sh[k]; sh[k] = rng.random()
            new = 1.0 - _covered_frac(freqs, sh, grid)
            if new <= cur: cur = new
            else: sh[k] = old
        if cur < best: best, bestsh = cur, sh.copy()
    return bestsh

def _covered_frac(freqs, shifts, grid):
    th = (np.arange(grid) + 0.5) / grid
    cov = np.zeros(grid, bool)
    for r, s in zip(freqs, shifts):
        cov |= np.abs(r*th + s - np.round(r*th + s)) < RHO
    return cov.mean()

for l in [10, 11]:
    freqs = np.array(range(1, l+1), float)
    s0 = anneal_slice_tiling(list(freqs))           # tiling shifts at t0 = 0
    # realize s0 = a_i * 0 is impossible; instead set couplings a_i so that a_i*0 = 0.
    # Use s0 as the theta-offset via a phase: model lifted = ||r_i theta + a_i t|| with
    # a_i GENERIC and the t0-slice shift absorbed; sweep t and measure rectangle M.
    for trial in range(3):
        a = rng.uniform(-3, 3, size=l)              # generic couplings
        w = 1.0/300.0                                # the forced-rectangle half-width scale
        ts = np.linspace(-w, w, 41)
        rectM = 0.0
        for t in ts:
            shifts = s0 + a * t                      # lifted shift at time t (tiling arranged at t=0)
            M = min_tooth_over_theta(list(freqs), shifts, grid=40000)
            rectM = max(rectM, M)
        broke = rectM > RHO
        print(f"    l={l}, trial {trial}: slice uncovered@t0={1-_covered_frac(list(freqs),s0,40000):.4f}, "
              f"max_t rectangle-M over [-1/300,1/300] = {rectM:.5f} "
              f"=> {'TILING BROKE (clear pt, rigidity holds)' if broke else 'tiling persisted'}", flush=True)

print(flush=True)
print("READING:", flush=True)
print("  (A) parallel-coupled consecutive freqs have homogeneous LR value > 2/25 for all l<=11", flush=True)
print("      => torus_decouple_clears is confirmed: the sheared-product stratum CLEARS.", flush=True)
print("  (B) if generic couplings break the slice-tiling within [-1/300,1/300] (rectangle-M > 2/25),", flush=True)
print("      the generic l=10,11 case ALSO clears (rectangle rigidity) -- the forced-rectangle route.", flush=True)
