#!/usr/bin/env python3
"""values.py -- the Hamiltonian-path value set of span-w ribbon words up to a bound N.

BFS over distinct profile state vectors (prefix H is monotone under appending a vertex, by the
insertion lemma, so prefixes with H > N are pruned).  Tracks, per vector, the uncrossed-cut mask
(bit j set iff the cut between vertices q-j-1 and q-j has no back arc across it yet; j < w) and a
'broken' flag (an uncrossed cut left the window: the word can never be strong again), so that the
strong-word value set is obtained as well.  Records the shortest witness word for every value.
Usage: values.py w N [maxvec]
"""
import sys, time
from ribbon import transfer_system, letters, word_str, tournament_from_word, ham_count_bruteforce, is_strong

def run(w, N, maxvec=2_000_000):
    order, states, mats, final = transfer_system(w)
    L = letters(w)
    finalset = set(final)
    # key: (vec_items, uncrossed_mask, broken) ; value: word
    init = (tuple(), 0, False)
    # vec for the empty word: {0: 1}
    cur = {((0, 1),): {(0, False): ()}}
    attained = {}        # H -> shortest word
    attained_strong = {} # H -> shortest strong word
    t0 = time.time()
    depth = 0
    nvec_total = 0
    while cur:
        depth += 1
        nxt = {}
        for vec, maskdict in cur.items():
            vd = dict(vec)
            for a in L:
                M = mats[a]
                nv = {}
                for (j, i), mult in M.items():
                    if i in vd:
                        nv[j] = nv.get(j, 0) + mult * vd[i]
                H = sum(c for j, c in nv.items() if j in finalset)
                if H > N:
                    continue
                nvt = tuple(sorted(nv.items()))
                # cut bookkeeping: new vertex q; cuts j=0..w-1 relative: cut between q-1-j and q-j.
                # back arcs of span s (1..w) cross cuts 0..s-1 (relative to the new vertex).
                smax = max(a) if a else 0
                crossed = (1 << smax) - 1
                for (mask, broken), word in maskdict.items():
                    # shift: the new cut (between q-1 and q) is bit 0, uncrossed unless a back arc exists
                    newmask = ((mask << 1) | 1) & ((1 << w) - 1)
                    dropped = (mask >> (w - 1)) & 1  # the oldest cut leaves the window
                    newbroken = broken or bool(dropped)
                    newmask &= ~crossed
                    nw = word + (a,)
                    if H not in attained or len(attained[H]) > len(nw):
                        attained[H] = nw
                    if newmask == 0 and not newbroken and depth >= 1:
                        if H not in attained_strong or len(attained_strong[H]) > len(nw):
                            attained_strong[H] = nw
                    d = nxt.setdefault(nvt, {})
                    key = (newmask, newbroken)
                    if key not in d or len(d[key]) > len(nw):
                        d[key] = nw
        cur = nxt
        nvec_total += len(cur)
        print(f"  depth {depth}: {len(cur)} vectors, {len(attained)} values, {len(attained_strong)} strong values, {time.time()-t0:.1f}s", flush=True)
        if len(cur) > maxvec:
            print("  (vector cap reached; truncating)")
            break
    return attained, attained_strong, depth

if __name__ == "__main__":
    w = int(sys.argv[1]); N = int(sys.argv[2])
    maxvec = int(sys.argv[3]) if len(sys.argv) > 3 else 2_000_000
    attained, strong, depth = run(w, N, maxvec)
    odd = [h for h in range(1, N + 1, 2)]
    miss = [h for h in odd if h not in attained]
    missS = [h for h in odd if h not in strong]
    print(f"w={w} N={N} depth={depth}: attained {len(attained)} values; odd misses <= N: {len(miss)}")
    print("misses:", miss[:80], "..." if len(miss) > 80 else "")
    print(f"strong-word values: {len(strong)}; strong misses <= N: {len(missS)}")
    print("strong misses:", missS[:80], "..." if len(missS) > 80 else "")
    # sanity: re-verify 30 witnesses by brute force
    import random
    rng = random.Random(1)
    items = list(attained.items())
    rng.shuffle(items)
    bad = 0
    for h, word in items[:30]:
        adj = tournament_from_word(w, list(word))
        if ham_count_bruteforce(adj) != h:
            bad += 1
    print("witness re-verification failures:", bad)
    items = list(strong.items()); rng.shuffle(items)
    badS = 0
    for h, word in items[:30]:
        adj = tournament_from_word(w, list(word))
        if ham_count_bruteforce(adj) != h or not is_strong(adj):
            badS += 1
    print("strong witness re-verification failures:", badS)
    with open(f"values_w{w}_N{N}.txt", "w") as f:
        for h in sorted(attained):
            f.write(f"{h}\t{word_str(list(attained[h]), w)}\t{'S' if h in strong else '-'}\t{word_str(list(strong[h]), w) if h in strong else ''}\n")
