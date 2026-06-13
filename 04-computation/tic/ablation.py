#!/usr/bin/env python3
"""
Ablation study: measure the exact impact of each GR4 feature.

Encode a real photo plane with:
  A. Pure GR2 baseline (64 ctx, EMA, no sign, no run)
  B. + sign correction
  C. + run mode
  D. + adaptive escape
  E. All (= GR4)

Python implementation for precise measurement.
opus-2026-04-02-S2
"""
import numpy as np
from PIL import Image
import sys, time, math

sys.stdout.reconfigure(line_buffering=True)

IMG = "../../inbox/processed/2026-04-02/new/IMG_2811.jpg"
pil = Image.open(IMG).convert('RGB')
cx, cy = pil.size[0]//2, pil.size[1]//2

# Test at multiple sizes
for SZ in [256, 512, 1024, 2048]:
    crop = pil.crop((cx-SZ//2, cy-SZ//2, cx+SZ//2, cy+SZ//2))
    arr = np.array(crop)
    h, w = arr.shape[:2]
    G = arr[:,:,1].astype(np.uint8)  # Just test luma (G) plane

    def med(a, b, c):
        mn, mx = min(a,b), max(a,b)
        if c >= mx: return mn
        if c <= mn: return mx
        return a + b - c

    def get_mag_ctx(a, b, c, d):
        g = abs(a-c) + abs(b-c) + abs(d-b)
        if g<3: return 0
        if g<8: return 1
        if g<15: return 2
        if g<25: return 3
        if g<40: return 4+(g-25)//5
        if g<100: return 7+(g-40)//10
        r = 13+(g-100)//20
        return min(r, 63)

    def get_sign_flip(a, b, c, d):
        D1, D2, D3 = d-b, b-c, c-a
        if D1<0 or (D1==0 and D2<0) or (D1==0 and D2==0 and D3<0):
            return -1
        return 1

    def is_flat(a, b, c, d):
        return abs(d-b)<=2 and abs(b-c)<=2 and abs(c-a)<=2

    def k_from(avg):
        if avg < 0.5: return 0
        k = int(math.log2(avg / 0.693) + 0.5)
        return max(0, min(k, 12))

    def zigzag(v):
        return 2*v if v >= 0 else 2*(-v)-1

    def gr_bits(val, k, esc=30):
        """Count bits for GR encoding."""
        v = zigzag(val)
        q = v >> k
        if q >= esc:
            return esc + 1 + 16  # escape
        return q + 1 + k

    def encode_plane(plane, use_sign=False, use_run=False, adaptive_esc=False):
        """Simulate encoding and count total bits."""
        total_bits = 0
        ctx_avg = [4.0] * 64
        run_ct = 0
        i = 0
        while i < h * w:
            r, c = i // w, i % w
            a = int(plane[r, c-1]) if c > 0 else 0
            b = int(plane[r-1, c]) if r > 0 else 0
            cv = int(plane[r-1, c-1]) if (r>0 and c>0) else 0
            d = int(plane[r-1, c+1]) if (r>0 and c+1<w) else b

            if use_run and is_flat(a, b, cv, d):
                # Count run
                run = 0
                j = i
                while j < h*w:
                    rj, cj = j//w, j%w
                    aj = int(plane[rj, cj-1]) if cj>0 else 0
                    bj = int(plane[rj-1, cj]) if rj>0 else 0
                    cvj = int(plane[rj-1, cj-1]) if (rj>0 and cj>0) else 0
                    dj = int(plane[rj-1, cj+1]) if (rj>0 and cj+1<w) else bj
                    if not is_flat(aj, bj, cvj, dj): break
                    pj = med(aj, bj, cvj)
                    if plane[rj, cj] != (pj & 0xFF): break
                    run += 1; j += 1
                rk = 0
                if run_ct>2: rk=1
                if run_ct>5: rk=2
                if run_ct>10: rk=3
                if run_ct>20: rk=4
                total_bits += gr_bits(run, rk)
                run_ct = (run_ct * 7 + run) // 8
                i += run
                # Interruption
                if i < h*w:
                    ri, ci = i//w, i%w
                    ai = int(plane[ri, ci-1]) if ci>0 else 0
                    bi = int(plane[ri-1, ci]) if ri>0 else 0
                    cvi = int(plane[ri-1, ci-1]) if (ri>0 and ci>0) else 0
                    di = int(plane[ri-1, ci+1]) if (ri>0 and ci+1<w) else bi
                    if is_flat(ai, bi, cvi, di):
                        pred = med(ai, bi, cvi)
                        res = (int(plane[ri, ci]) - pred)
                        if res > 128: res -= 256
                        elif res < -128: res += 256
                        k = k_from(ctx_avg[0])
                        esc = (32-k) if adaptive_esc else 30
                        total_bits += gr_bits(res, k, esc)
                        ctx_avg[0] = 0.9*ctx_avg[0] + 0.1*abs(res)
                        i += 1
                continue

            pred = med(a, b, cv)
            res = int(plane[r, c]) - pred
            if res > 128: res -= 256
            elif res < -128: res += 256

            cx = get_mag_ctx(a, b, cv, d)
            if use_sign:
                sf = get_sign_flip(a, b, cv, d)
                coded = res * sf
            else:
                coded = res

            k = k_from(ctx_avg[cx])
            esc = (32-k) if adaptive_esc else 30
            total_bits += gr_bits(coded, k, esc)
            ctx_avg[cx] = 0.9*ctx_avg[cx] + 0.1*abs(coded)
            i += 1

        return total_bits

    print(f"\n{'='*70}")
    print(f"  ABLATION STUDY: {SZ}x{SZ} crop of IMG_2811 (G plane)")
    print(f"{'='*70}")

    configs = [
        ("A. GR2 baseline",      False, False, False),
        ("B. + sign correction",  True,  False, False),
        ("C. + run mode",         False, True,  False),
        ("D. + adaptive escape",  False, False, True),
        ("E. sign + run",         True,  True,  False),
        ("F. ALL (GR4)",          True,  True,  True),
    ]

    baseline_bits = None
    for label, sign, run, esc in configs:
        t0 = time.time()
        bits = encode_plane(G, use_sign=sign, use_run=run, adaptive_esc=esc)
        dt = time.time() - t0
        bpp = bits / (w * h)
        if baseline_bits is None:
            baseline_bits = bits
        diff = (bits - baseline_bits) / (w * h)
        print(f"  {label:<26} {bits:>12,} bits  {bpp:.3f} bpp  {diff:+.3f} vs baseline  ({dt:.1f}s)")
