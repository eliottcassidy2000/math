import subprocess, sys
LABELG = '/opt/homebrew/bin/labelg'
def rows_from_d6(s):
    n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        for k in range(5, -1, -1): bits.append((v >> k) & 1)
    rows = [0] * n
    for i in range(n):
        for j in range(n):
            if bits[i * n + j]: rows[i] |= 1 << j
    return rows
start, stop = int(sys.argv[1]), int(sys.argv[2])
# parents: all ready classes (digit 1 or 2) at the start order from extgen
cur = [l.split()[3] for l in open('extgen_order%d.txt' % start)]
print('order %d: %d ready parents' % (start, len(cur)), flush=True)
targets = [4289, 4291, 12869, 12871, 12875, 12877, 38609, 38611, 38615, 38617, 38627, 38629, 38633, 38635]
for N in range(start + 1, stop + 1):
    inp = '\n'.join(' '.join(map(str, rows_from_d6(c))) for c in cur) + '\n'
    out = subprocess.run(['./gen2', str(N - 1)], input=inp, capture_output=True, text=True).stdout.split('\n')
    recs = [l.split() for l in out if l.strip()]
    can = subprocess.run([LABELG, '-q'], input='\n'.join(r[1] for r in recs) + '\n', capture_output=True, text=True).stdout.split()
    cls = {}
    for r, c in zip(recs, can): cls[c] = int(r[0])
    vals = sorted(set(cls.values()))
    print('order %d: %d digit-2 classes; R2(%d): %d values in [%d, %d]; tree values present: %s' % (N, len(cls), N, len(vals), vals[0], vals[-1], [t for t in targets if t in set(vals)]), flush=True)
    with open('gen2_order%d.txt' % N, 'w') as fh:
        for c, h in cls.items(): fh.write('%d %s\n' % (h, c))
    # next parents: digit-2 classes just generated plus digit-1 classes of the same order from gen1 (if available)
    try: d1 = [l.split()[1] for l in open('gen1_order%d.txt' % N)]
    except FileNotFoundError: d1 = []
    cur = list(cls) + d1
