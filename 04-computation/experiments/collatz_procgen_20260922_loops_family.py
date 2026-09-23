#!/usr/bin/env python3
"""collatz_procgen_20260922_loops_family.py -- the insertion closure of E-loops through 1.

Insertion lemma (proved in the note): if a loop through 1 of length s (s multiplications, K0(s) halvings)
passes through the hub h, and C is a cycle through h with q multiplications and p halvings, then splicing C
in at h gives a loop of length s+q with K0(s)+p halvings; this equals K0(s+q) iff 2^(K0(s)+p) < 3^(s+q+1),
i.e. iff the new ratio stays < 3.

Hub cycles used (reverse moves from the hub, verified exactly below):
   C1  hub 1   (q,p)=(1,2)   k: 2              1 -> 4 -> 2 -> 1
   C3  hub 4   (q,p)=(3,5)   k: 2,3,0          4 -> 13 -> 40 -> 20 -> 10 -> 5 -> 16 -> 8 -> 4
   C5  hub 13  (q,p)=(5,8)   k: 6,0,1,0,1      13 -> 40 -> 20 -> 61 -> 184 -> 92 -> 277 -> 832 -> ... -> 13
Every loop with ratio <= 128/81 has climb >= 2 (proved), hence passes through 1, 4, 13, and splicing
preserves passing through 1, 4, 13.  So all loops generated from base loops that pass through 13 do too.

s is IRREDUCIBLE if no q in {1,3,5} with s-q >= 2 has K0(s) = K0(s-q) + p_q  (equivalently, for s >= 7,
2^eps(s+1) <= 256/243, i.e. c(s) <= 128/81).  Base loops for irreducible s are read from files
('s=.. k:..' lines, e.g. produced by loops_recon1 with XMAX = H(s)); every other s is built by splicing.

Usage: python3 ..._loops_family.py SMAX baseglob OUTFILE [hubcycles]  (then ..._loops_verify.py OUTFILE)
       python3 ..._loops_family.py SMAX --list                (print the irreducible s <= SMAX)
With the optional file of record hub cycles ('q=<q> ... minK=<p> ... k:..' lines, reverse moves from the hub
1093 = y_6, produced by loops_recon1 q XMAX B 1093), s is spliced as C_q + loop(s-q) whenever
K0(s) = K0(s-q)+p_q and loop(s-q) visits 1093 (record cycles are preferred); the base set is then determined
dynamically: s needs a base loop iff no splice applies to the loops already built.
"""
import sys, re, glob

def K0(s):
    return (3 ** (s + 1)).bit_length() - 1      # largest K with 2^K < 3^(s+1)  (3^(s+1) is not a power of 2)

CYC = {1: (1, 2, [2]), 3: (4, 5, [2, 3, 0]), 5: (13, 8, [6, 0, 1, 0, 1])}   # q: (hub, p, reverse moves)

def check_cycle(h, p, ks):
    x = h
    for k in ks:
        t = (1 << k) * x - 1; assert t % 3 == 0; x = t // 3; assert x >= 1 and x % 3
    assert x == h and sum(ks) == p

def irreducible(s):
    # s = 2 (loop 1->4->2->1 twice) is taken as a base loop; for s >= 3 the splice C1 at hub 1 applies to
    # loop(s-1) whenever K0(s) = K0(s-1)+2, and C3/C5 need hubs 4/13 (all loops with s >= 4 pass through 13).
    if s == 2:
        return True
    for q, (h, p, ks) in CYC.items():
        if s - q >= 2 and K0(s) == K0(s - q) + p and not (q == 5 and s - q < 4):
            return False
    return True

def path_values(ks):
    xs = [1]
    for k in ks:
        xs.append(((1 << k) * xs[-1] - 1) // 3)
    return xs

def splice(ks, q):
    h, p, cks = CYC[q]
    xs = path_values(ks)
    i = xs.index(h)                      # first visit of the hub on the reverse path
    return ks[:i] + cks + ks[i:]

def main():
    SMAX = int(sys.argv[1])
    for q, (h, p, cks) in CYC.items():
        check_cycle(h, p, cks)
    irr = [s for s in range(2, SMAX + 1) if irreducible(s)]
    if sys.argv[2] == '--list':
        print(' '.join(map(str, irr))); return
    base = {}
    for fn in glob.glob(sys.argv[2]):
        for line in open(fn):
            m = re.search(r's=\s*(\d+).*?k:(\S+)', line)
            if m:
                base[int(m.group(1))] = [int(t) for t in m.group(2).split(',')]
    rec = {}
    if len(sys.argv) > 4:
        for line in open(sys.argv[4]):
            m = re.search(r'q=(\d+).*?minK=(\d+).*?k:(\S+)', line)
            if m:
                q, p = int(m.group(1)), int(m.group(2)); cks = [int(t) for t in m.group(3).split(',')]
                check_cycle(1093, p, cks); assert p == (3 ** q).bit_length(); rec[q] = (1093, p, cks)
    if rec:
        CYC.update(rec)
        loops = {}; how = {}
        for s in range(2, SMAX + 1):
            done = False
            for q in sorted(rec, reverse=True) + [5, 3, 1]:
                h, p, cks = CYC[q]
                if s - q >= 2 and K0(s) == K0(s - q) + p and h in path_values(loops[s - q]):
                    loops[s] = splice(loops[s - q], q); how[s] = 'C%d+loop(%d)' % (q, s - q); done = True; break
            if not done:
                if s not in base:
                    raise SystemExit('missing base loop for s=%d' % s)
                loops[s] = base[s]; how[s] = 'base'
        with open(sys.argv[3], 'w') as f:
            for s in range(2, SMAX + 1):
                f.write('s=%d %s k:%s\n' % (s, how[s], ','.join(map(str, loops[s]))))
        bs = [s for s in how if how[s] == 'base']
        print('SMAX=%d with %d record hub cycles at 1093: %d base loops, %d spliced; base s: %s' % (
            SMAX, len(rec), len(bs), len(how) - len(bs), bs))
        return
    loops = {}; how = {}
    for s in range(2, SMAX + 1):
        if s in irr:
            if s not in base:
                raise SystemExit('missing base loop for irreducible s=%d' % s)
            loops[s] = base[s]; how[s] = 'base'
            assert s < 4 or 13 in path_values(loops[s]), 'base loop %d misses hub 13' % s
            continue
        for q in (5, 3, 1):
            h, p, cks = CYC[q]
            if s - q >= 2 and K0(s) == K0(s - q) + p and h in path_values(loops[s - q]):
                loops[s] = splice(loops[s - q], q); how[s] = 'C%d+loop(%d)' % (q, s - q); break
        else:
            raise SystemExit('no splice available for s=%d' % s)
    with open(sys.argv[3], 'w') as f:
        for s in range(2, SMAX + 1):
            f.write('s=%d %s k:%s\n' % (s, how[s], ','.join(map(str, loops[s]))))
    nb = sum(1 for s in how if how[s] == 'base')
    print('SMAX=%d: %d irreducible (base) loops, %d built by splicing; irreducible s: %s' % (
        SMAX, nb, len(how) - nb, irr))

if __name__ == '__main__':
    main()
