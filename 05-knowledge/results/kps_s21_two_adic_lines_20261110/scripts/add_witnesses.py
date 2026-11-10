"""add_witnesses.py - kind-pasteur-2026-11-10-S21 (after the audit, minor 15).
Rewrites ../below_rays.txt so that every certificate is self-contained: each value is given with a witness
tournament (order, upper-triangle code) taken from the independent class census of the audit
(../audit/out/values_witness.txt: value, least order, code).
  c s N census ORDER CODE
  c s N product A B ORDER_A CODE_A ORDER_B CODE_B        (order join of the two witnesses; H is multiplicative)
  c s N GAP (THM-1370)
verify_lines.py recounts every witness literally."""
import os
here = os.path.dirname(os.path.abspath(__file__))
wit = {}
for line in open(os.path.join(here, "..", "audit", "out", "values_witness.txt")):
    v, m, code = line.split()[:3]
    wit[int(v)] = (int(m), int(code))
src = os.path.join(here, "..", "below_rays.txt")
out = []
for line in open(src):
    if line.startswith("#"):
        out.append("# c s N certificate   (points 1 + c 2^s below the first ray, odd c <= 129; witnesses: order and upper-triangle code)\n")
        continue
    p = line.split()
    c, s, N, kind = p[0], p[1], int(p[2]), p[3]
    if kind == "census":
        m, code = wit[N]
        out.append("%s %s %d census %d %d\n" % (c, s, N, m, code))
    elif kind == "product":
        a, b = int(p[4]), int(p[5])
        (ma, ca), (mb, cb) = wit[a], wit[b]
        out.append("%s %s %d product %d %d %d %d %d %d\n" % (c, s, N, a, b, ma, ca, mb, cb))
    else:
        out.append(line)
open(src, "w", newline="\n").write("".join(out))
print("rewrote", src, "with", len(out) - 1, "certificates")
