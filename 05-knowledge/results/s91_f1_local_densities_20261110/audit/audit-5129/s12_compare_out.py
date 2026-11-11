# audit-5129, step 12: comparison of the auditor's lists with the committed outputs of the session (run only
# after the auditor's numbers were final).  Reads the .out files as text; nothing of the session is imported.
import re, sys

base = sys.argv[1]


def parse_lists(path):
    txt = open(path).read().split("\n")
    out = {}
    cur = None
    for line in txt:
        m = re.search(r"D = (\d+): (\d+) odd parts", line) or re.search(r"offset (\d+) \(level s = \d+\): (\d+) odd parts", line)
        if m:
            cur = int(m.group(1))
            out[cur] = {}
            continue
        if cur is not None and re.match(r"^\s+\d+:\d+", line):
            for tok in line.split():
                tok = tok.rstrip("*")
                d, j = tok.split(":")
                out[cur][int(d)] = int(j)
        elif cur is not None and line.strip() and not re.match(r"^\s+\d+:\d+", line):
            cur = None
    return out


mine = {}
for line in open("s3_slope_lists.txt"):
    D, rest = line.split("|")
    mine[int(D)] = {int(t.split(":")[0]): int(t.split(":")[1]) for t in rest.split()}
for D in (129, 257):
    mine[D] = {int(t.split(":")[0]): int(t.split(":")[1]) for t in open("s11_slopes_%d.txt" % D).read().split()}

theirs = parse_lists(base + "/source_supply.out")
theirs.update(parse_lists(base + "/source_supply_levels78.out"))
for D in sorted(theirs):
    print("offset %d: session %d odd parts, auditor %d, identical (odd parts and least exponents): %s" %
          (D, len(theirs[D]), len(mine[D]), theirs[D] == mine[D]))
