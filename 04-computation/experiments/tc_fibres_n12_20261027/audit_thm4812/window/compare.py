"""compare.py -- compare the independent window enumeration (window_enum.out) with the repository's window_scan.out.
Name mapping: mine a_n{n}_z{z}_c{c}_r{r}_p{p}_{name}, characters ordered (C3_1..C3_c, P7); the repository's
n{n}_K{z}_C{c}_R{r}_P{p}_{name}, ordered (P7, C3_1..C3_c).  My C3 copy is the converse of theirs, so the relabelling
x -> -x identifying them conjugates a rotation by e into a rotation by -e: C3 character coordinates change sign;
the P7 copies (non-residues vs residues) are identified by x -> -x, which keeps the multiplier, so the P7 coordinate
is unchanged.  Normalisation: first nonzero coordinate = 1."""
import re, sys

def norm(v):
    lead = next(x for x in v if x)
    inv = pow(lead, -1, 3)
    return tuple((x * inv) % 3 for x in v)

mine = {}
for line in open(sys.argv[1]):
    m = re.match(r'jobs/a_n(\d+)_z(\d+)_c(\d+)_r(\d+)_p(\d+)_(\w+)\.job: .*strong=(\d+) minH=(\d+) repo-style-autcalls=(\d+) hits=(\d+) minTc=(-?\d+)', line)
    if not m: continue
    n, z, c, r, p, name, st, mh, ac, h, mt = m.groups()
    n, z, c, r, p = map(int, (n, z, c, r, p))
    if name == 'P':
        theirs = 'P'
    else:
        v = tuple(int(x) for x in name[3:])
        cs, pv = v[:c], v[c:]
        w = tuple(pv) + tuple((-x) % 3 for x in cs)
        theirs = 'ker' + ''.join(map(str, norm(w)))
    key = (n, z, c, r, p, theirs)
    mine[key] = (int(st), int(mh), int(ac), int(h), int(mt), line.strip())

repo = {}
for line in open(sys.argv[2]):
    m = re.match(r'window_jobs/n(\d+)_K(\d+)_C(\d+)_R(\d+)_P(\d+)_(\w+)\.job: .*strong=(\d+) minH=(\d+) autcalls=(\d+) hits=(\d+)', line)
    if not m: continue
    n, z, c, r, p, name, st, mh, ac, h = m.groups()
    repo[(int(n), int(z), int(c), int(r), int(p), name)] = (int(st), int(mh), int(ac), int(h))

ok = 0; bad = []
for k, v in repo.items():
    if k not in mine:
        bad.append(('missing in mine', k)); continue
    if mine[k][:4] != v:
        bad.append(('differ', k, v, mine[k][:4]))
    else:
        ok += 1
extra = [k for k in mine if k not in repo]
print(f"repository jobs {len(repo)}, mine {len(mine)}, identical (strong, minH, autcalls, hits): {ok}")
print(f"differences: {bad}")
print(f"jobs only in mine: {extra}")
print(f"total strong (mine) {sum(v[0] for v in mine.values())}, total hits (mine) {sum(v[3] for v in mine.values())}, "
      f"min tc over all window tournaments {min(v[4] for v in mine.values() if v[4] >= 0)}")
