# #158 Prop 7.2 (five-cycle exclusion), combinatorial core: six positions j in Z/6 carry palette indices p_j in Z/5
# (palette D_p = {p, p+2}); R-adjacent palettes disjoint  <=>  p_{j+1}-p_j = +-1 mod 5;
# antipodal positions (j, j+3) use the same two-element option set {D_a, D_{a+2}} or {D_a,D_{a-2}}:
# difference p_{j+3}-p_j in {0, +2, -2} mod 5.  Claim: no such cyclic word exists.
import itertools
sols = []
for p0 in range(5):
    for steps in itertools.product((1, -1), repeat=6):
        if sum(steps) % 5: continue
        p = [p0]
        for s in steps[:-1]: p.append((p[-1] + s) % 5)
        if all(((p[(j+3) % 6] - p[j]) % 5) in (0, 2, 3) for j in range(6)):
            sols.append(p)
print("closed six-position words with +-1 steps and antipodal difference in {0,+-2}:", len(sols))
# Also the stronger version actually used (option set fixed to {D_1,D_4} up to relabelling):
# palettes D_1={1,3}, D_4={4,1}: indices 1 and 4 -> antipodal difference in {0, 3} = {0,-2}
