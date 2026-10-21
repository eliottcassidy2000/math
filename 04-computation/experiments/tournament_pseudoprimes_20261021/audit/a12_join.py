# AUDIT a12_join: brute force (all permutations, a12_circ_brute) vs subset DP (a11_circ) for every connection set n <= 11.
import subprocess, re
bad = 0; tot = 0
for f in ["a12_circ_brute.out", "a12_circ_brute_11.out"]:
    for line in open(f):
        m = re.match(r"n=(\d+) S=\{([\d,]+)\} H=(\d+) C=(\d+) rho0=(\d+)", line)
        n, S, H, C, R = m.groups()
        out = subprocess.run(["./a11_circ", n, "set", S], capture_output=True, text=True).stdout
        mm = re.search(r"H=(\d+) Hfull=(\d+) C=(\d+) rho0=(\d+)", out)
        tot += 1
        if (mm.group(1), mm.group(3), mm.group(4)) != (H, C, R) or mm.group(2) != H:
            bad += 1; print("MISMATCH", line.strip(), out.strip())
print(f"brute force (all permutations) vs subset DP: {tot} connection sets, mismatches {bad}")
