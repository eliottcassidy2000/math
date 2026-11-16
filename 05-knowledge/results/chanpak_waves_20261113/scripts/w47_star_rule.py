"""Star rule: in one-parent growth, a star (root 0, leaves 1..m) is the only shape with H = 1; attaching the next
vertex to leaf i gives H = 1 + 2^(m-i+1), and to the root gives H = 1.  Needs w46_treeH compiled in this folder."""
import subprocess
lines, keys = [], []
for m in range(1, 9):
    for i in range(0, m + 1):
        lines.append(f"{m + 2} " + " ".join(map(str, [0] * m + [i]))); keys.append((m, i))
out = subprocess.run(['./w46_treeH'], input="\n".join(lines) + "\n", capture_output=True, text=True).stdout.split()
bad = sum(int(h) != (1 if i == 0 else 1 + 2 ** (m - i + 1)) for (m, i), h in zip(keys, out))
print("star rule mismatches:", bad, "of", len(keys))
