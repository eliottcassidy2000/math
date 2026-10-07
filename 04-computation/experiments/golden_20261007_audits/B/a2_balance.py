# Balance 1/2((m/2)^th / s + s / 2^th) = 1 for (mx+1)/2 maps, 0 < th < 1.
# Roots: s = 2^th (1 +- sqrt(1 - (m/4)^th)); real iff m <= 4; smaller root s_- < 1 needed for THM-4581's weight.
import math
for m in (1, 2, 3, 3.5, 3.9, 3.99, 4, 4.01, 5):
    best = None; real_any = False
    for i in range(1, 1000):
        th = i / 1000
        disc = 1 - (m / 4) ** th
        if disc < -1e-15: continue
        real_any = True
        s = 2 ** th * (1 - math.sqrt(max(disc, 0)))
        lhs = 0.5 * ((m / 2) ** th / s + s / 2 ** th)
        assert abs(lhs - 1) < 1e-9
        if best is None or s < best[0]: best = (s, th)
    print(f"m={m:5}: real root for some th in (0,1): {real_any}; min over th of s_- = "
          f"{best[0] if best else None!s:.24} at th={best[1] if best else None}; s<1 possible: {bool(best and best[0] < 1)}")
# m = 3: s(th) < 1 on all of (0,1)?
print("m=3: max over th in (0,1) of s(th) =", max(2**t*(1-math.sqrt(1-0.75**t)) for t in [i/1000 for i in range(1,1000)]))
