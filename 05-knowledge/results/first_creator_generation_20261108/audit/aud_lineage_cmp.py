"""aud_lineage_cmp.py - compares the packet's lineage_le10.txt with the auditor's aud_lineage_le10.txt.
The packet abbreviates long parent lists as 'p1,...,pk,...,q1,..,qm (N values)'; those are compared by head, tail and count."""
import re
def load(fn):
    d = {}
    for line in open(fn):
        if line.startswith('#') or not line.strip(): continue
        p = line.split(None, 4)
        d[int(p[0])] = (int(p[1]), int(p[2]), int(p[3]), p[4].strip() if len(p) > 4 else '')
    return d
P, A = load('../lineage_le10.txt'), load('aud_lineage_le10.txt')
bad, abbrev = [], 0
for h in sorted(set(P) | set(A)):
    if h not in P or h not in A: bad.append((h, 'missing')); continue
    (pn, pc, pp, pl), (an, ac, ap, al) = P[h], A[h]
    if (pn, pc, pp) != (an, ac, ap): bad.append((h, 'counts', P[h], A[h])); continue
    if al == 'ORPHAN' or pl == 'ORPHAN':
        if al != pl: bad.append((h, 'orphan', pl, al))
        continue
    full = [int(x) for x in al.split(',')]
    m = re.match(r'^(.*)\s+\((\d+) values\)$', pl)
    if m:
        abbrev += 1
        body, cnt = m.group(1), int(m.group(2))
        head, tail = body.split(',...,')
        head = [int(x) for x in head.split(',')]; tail = [int(x) for x in tail.split(',')]
        if cnt != len(full) or full[:len(head)] != head or full[-len(tail):] != tail: bad.append((h, 'abbrev', pl, al))
    else:
        if [int(x) for x in pl.split(',')] != full: bad.append((h, 'list', pl, al))
print("values: packet %d, auditor %d; abbreviated lines compared by head/tail/count: %d; disagreements: %d" % (len(P), len(A), abbrev, len(bad)))
for b in bad[:10]: print("  ", b)
orph = [h for h in A if A[h][3] == 'ORPHAN']
print("orphan values (no first creator has a first-creator deletion):", orph)
print("values whose first creators are not all grown from first creators:", [h for h in A if A[h][1] != A[h][2]])
