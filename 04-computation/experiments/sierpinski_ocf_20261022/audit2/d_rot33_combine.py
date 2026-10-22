# Audit2 (D): combine the four partial inclusion-exclusion sums for H(R_33) (hcount.c "iecircpart k 4 33 1..16",
# masks X = k mod 4) and compare with the session's R33.out.
import re, glob
parts = {}
for f in sorted(glob.glob('d_rot33_parts/part_*.out')):
    m = re.match(r'part (\d) of 4: (-?\d+) \((\d+)s\)', open(f).read().strip())
    parts[int(m.group(1))] = (int(m.group(2)), int(m.group(3)))
assert sorted(parts) == [0, 1, 2, 3], parts
H = sum(v for v, t in parts.values())
repo = int(re.search(r'H=(\d+)', open('../R33.out').read()).group(1))
print(f"H(R_33) = {H}  (parts: {[parts[k][0] for k in range(4)]}, times {[parts[k][1] for k in range(4)]} s)")
print(f"session R33.out: {repo}  -> {'MATCH' if H == repo else 'MISMATCH'}")
