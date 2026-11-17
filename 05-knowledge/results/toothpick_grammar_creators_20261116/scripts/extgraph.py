#!/usr/bin/env python3
"""extgraph.py (R-3, R-4): the one-vertex extension graph between tournament classes of orders n-1 and n (n = 3..9):
canonical forms (nauty labelg, digraph6) of every class and of every one-vertex deletion; distinct parents per class
(contacts), distinct children per class, childless classes, and the canonical-parent tree (parent = the least canonical
deletion string).  Writes cparent_n.txt (one parent index per class, gentourng order) and parents_n.txt (distinct parent
indices per class) for n <= 9, and nparents_n.txt."""
import subprocess, sys
from collections import Counter
def d6_decode(s):
    s = s.strip(); n = ord(s[1]) - 63; bits = ''.join(format(ord(c) - 63, '06b') for c in s[2:])
    return n, [[int(bits[i*n+j]) for j in range(n)] for i in range(n)]
def d6_encode(n, M):
    bits = ''.join(str(M[i][j]) for i in range(n) for j in range(n)); bits += '0' * ((-len(bits)) % 6)
    return '&' + chr(n + 63) + ''.join(chr(int(bits[i:i+6], 2) + 63) for i in range(0, len(bits), 6))
def canon(lines):
    out = subprocess.run(['labelg', '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True).stdout.split()
    assert len(out) == len(lines), (len(out), len(lines)); return out
prev_index = None; prev_canon = None
print("n  classes  contacts  parents/class(min,mean,max)  children/class(min,mean,max; max attained by transitive?)  childless  tree-branching(max, childless, mean)")
for n in range(2, 10):
    raw = subprocess.run(['gentourng', '-q', '-z', str(n)], capture_output=True, text=True).stdout.split() if n >= 3 else ['&A_']
    if n == 2: raw = [d6_encode(2, [[0, 1], [0, 0]])]
    cn = canon(raw); index = {s: i for i, s in enumerate(cn)}
    assert len(index) == len(raw)
    if prev_index is not None:
        dels = []; owner = []
        for i, s in enumerate(raw):
            m, A = d6_decode(s)
            for v in range(m):
                keep = [u for u in range(m) if u != v]; dels.append(d6_encode(m - 1, [[A[a][b] for b in keep] for a in keep])); owner.append(i)
        cd = canon(dels)
        parents = [set() for _ in raw]; cparent = [None] * len(raw)
        for s, i in zip(cd, owner):
            p = prev_index[s]; parents[i].add(p)
            if cparent[i] is None or prev_canon[p] < prev_canon[cparent[i]]: cparent[i] = p
        children = [set() for _ in prev_canon]
        for i, ps in enumerate(parents):
            for p in ps: children[p].add(i)
        np_ = [len(p) for p in parents]; nc = [len(c) for c in children]; contacts = sum(np_)
        trans_prev = prev_index[canon([d6_encode(n - 1, [[1 if a < b else 0 for b in range(n - 1)] for a in range(n - 1)])])[0]]
        branch = Counter(cparent); bl = [branch.get(p, 0) for p in range(len(prev_canon))]
        print(f"{n}  {len(raw)}  {contacts}  ({min(np_)}, {contacts/len(raw):.3f}, {max(np_)})  ({min(nc)}, {sum(nc)/len(nc):.3f}, {max(nc)}; transitive has {nc[trans_prev]}, {'max' if nc[trans_prev] == max(nc) else 'not max'})  {sum(1 for c in nc if c == 0)}  ({max(bl)}, {sum(1 for b in bl if b == 0)}, {sum(bl)/len(bl):.3f})")
        with open(f'cparent_{n}.txt', 'w') as f: f.write('\n'.join(str(p) for p in cparent) + '\n')
        with open(f'nparents_{n}.txt', 'w') as f: f.write('\n'.join(str(k) for k in np_) + '\n')
        if n <= 9:
            with open(f'parents_{n}.txt', 'w') as f: f.write('\n'.join(' '.join(map(str, sorted(ps))) for ps in parents) + '\n')
        sys.stdout.flush()
    prev_index = index; prev_canon = cn
