"""Brute-force check of the fixed-port tail transfer (THM-5085 Lemma 1): build the actual tournament and compare I(Omega,2), I(Omega,3)
with the predicted h.  Uses packing_count (C) from the S18 packet for the independent count."""
import random, subprocess, sys
import os
# FP_DIR: directory holding classgen's cls_<n>.txt and ports_le6_n.txt (port_stats output with an order prefix);
# FP_PACKING_COUNT: the compiled packing_count program from the S18 packet (kps_s18_.../packing_count.c).
S = os.environ.get('FP_DIR', '.')
PC = os.environ.get('FP_PACKING_COUNT', './packing_count')
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
classes = {}
for n in range(1, 7):
    for line in open(f'{S}/cls_{n}.txt'):
        classes.setdefault(n, []).append(int(line.split()[0]))
stats = {}
for line in open(f'{S}/ports_le6_n.txt'):
    nn, code, u, h2, E2, O2, hu2, h3, E3, O3, hu3 = map(int, line.split())
    stats[(nn, code, u)] = (h2, E2, O2, hu2, h3, E3, O3, hu3)

def build_core(n, code):
    A = [[0] * n for _ in range(n)]
    b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code >> b & 1: A[i][j] = 1
            else: A[j][i] = 1
            b += 1
    return A

def add_vertex(A, beats_list, beaten_by_all_except):
    """append a vertex z: z -> v for v in beats_list, v -> z otherwise"""
    n = len(A)
    for r in A: r.append(0)
    A.append([0] * (n + 1))
    for v in range(n):
        if v in beats_list: A[n][v] = 1
        else: A[v][n] = 1
    return n

lines, preds = [], []
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 300):
    n = rng.randint(1, 6)
    code = rng.choice(classes[n])
    u = rng.randrange(n)
    A = build_core(n, code)
    h2, E2, O2, hu2, h3, E3, O3, hu3 = stats[(n, code, u)]
    st = {2: [h2, E2 + O2, E2 - O2, hu2], 3: [h3, E3 + O3, E3 - O3, hu3]}
    total = n
    while total < 15:
        letter = rng.choice(['S', 'I0', 'I1', 'I2'])
        t = 0 if letter == 'S' else int(letter[1])
        need = 1 if letter == 'S' else t + 1
        if total + need > 15: break
        if letter == 'S':
            add_vertex(A, [], None)          # beaten by everyone (incl. u)
            for x in (2, 3):
                h, F, D, hu = st[x]; st[x] = [h, 2 * F, 0, hu]
        else:
            for k in range(t):
                add_vertex(A, [], None)      # TT_t vertices: beaten by all earlier (increasing order)
            add_vertex(A, [u], None)         # w: beats only u
            for x in (2, 3):
                h, F, D, hu = st[x]
                off = x * (F + D) // 2 if t == 0 else x * (2 ** (t - 1)) * F
                st[x] = [h + off, 2 ** (t + 1) * F - hu, -hu, hu]
        total = len(A)
    lines.append(f"{len(A)} " + "".join(str(A[i][j]) for i in range(len(A)) for j in range(len(A))))
    preds.append((st[2][0], st[3][0]))
out = subprocess.run([PC], input="\n".join(lines) + "\n", capture_output=True, text=True).stdout.split("\n")
bad = 0
for (p2, p3), o in zip(preds, out):
    ham, i2, i3 = map(int, o.split())
    if not (ham == i2 == p2 and i3 == p3):
        bad += 1
        if bad < 5: print("MISMATCH", p2, p3, ham, i2, i3)
print(len(preds), "random fixed-port tails checked (orders <= 15):", bad, "mismatches")
