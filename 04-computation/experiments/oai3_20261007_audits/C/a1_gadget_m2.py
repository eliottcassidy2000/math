# m=2 minimal-shift gadget in the region B=2h+c < a-1 (the only region not covered by the tangent-line argument)
import math
NA = 6000
ls = [0.0, 0.0]  # log sigma_m, index m
for m in range(2, NA + 5): ls.append(ls[-1] + math.log1p(1/(3*(m-1))))
worst = (9, None)
for a in range(3, NA):
    c = -(-(a-1)//2)
    for h in range(1, a):
        B = 2*h + c
        if B >= a - 1: break
        # log of P_min(a,B) / (2 P_min(a,h)), both below the diagonal
        r = ls[B] + math.log(2*a + B - 1) - ls[h] - math.log(2*a + h - 1) - math.log(2)
        if r < worst[0]: worst = (r, (a, h, B))
print("min log-ratio over a<6000 (should be > 0):", worst, "ratio", math.exp(worst[0]))
