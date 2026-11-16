# Exhaustive: quadratic f(y) = sum D_i y_i + sum_{i<j} s_ij y_i y_j on {0,1}^3, bijective onto 8 consecutive integers.
# Membership coordinates of a tournament bank: s_ij = -(W(a,b) + W(b,a)) <= -1 for every pair (cut functions are
# strictly submodular).  How many bijections have all three s_ij negative?  (f(000) = 0 by translation.)
import itertools
k = 3; N = 8; total = 0; neg = []; signs = {}
for c in range(-7, 1):
    rng = [v for v in range(c, c + N) if v != 0]
    for vals in itertools.permutations(rng, 6):          # values at 100, 010, 001, 110, 101, 011
        D = vals[:3]
        s12 = vals[3] - D[0] - D[1]; s13 = vals[4] - D[0] - D[2]; s23 = vals[5] - D[1] - D[2]
        f111 = sum(D) + s12 + s13 + s23
        if f111 == 0 or f111 in vals or not (c <= f111 < c + N): continue
        total += 1
        key = tuple(sorted(( (s12 > 0) - (s12 < 0), (s13 > 0) - (s13 < 0), (s23 > 0) - (s23 < 0) )))
        signs[key] = signs.get(key, 0) + 1
        if s12 < 0 and s13 < 0 and s23 < 0: neg.append((D, (s12, s13, s23)))
print("quadratic bijections of the 3-cube onto 8 consecutive integers:", total)
print("by sign pattern of the three interactions (sorted signs: count):", signs)
print("with all three interactions negative (strictly submodular):", len(neg))
for D, s in neg[:12]: print("   D =", D, " s =", s)
