# Batch canonical forms of tournaments (as digraphs) with nauty's labelg, in digraph6 format.
import subprocess
def to_d6(adj):
    n = len(adj)
    bits = []
    for i in range(n):
        for j in range(n):
            bits.append((adj[i] >> j) & 1)
    while len(bits) % 6: bits.append(0)
    s = '&' + (chr(63 + n) if n <= 62 else None)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]: v = (v << 1) | b
        s += chr(63 + v)
    return s
def canon_batch(adjs):
    """return a list of canonical digraph6 strings (one per input), computed by a single labelg call"""
    if not adjs: return []
    inp = '\n'.join(to_d6(a) for a in adjs) + '\n'
    out = subprocess.run(['labelg', '-q'], input=inp, capture_output=True, text=True).stdout.split('\n')
    out = [l for l in out if l]
    assert len(out) == len(adjs), (len(out), len(adjs))
    return out
if __name__ == '__main__':
    # self-test: isomorphic relabellings get the same canonical form; non-isomorphic classes differ
    import itertools, random
    from rooted import tournaments, relabel
    for n in (4, 5, 6, 7):
        T = tournaments(n)
        c = canon_batch(T)
        assert len(set(c)) == len(T)
        rnd = random.Random(n)
        rel = [relabel(a, rnd.sample(range(n), n)) for a in T]
        assert canon_batch(rel) == c
        print(f"n = {n}: {len(T)} classes, labelg canonical forms distinct and relabelling-invariant")
