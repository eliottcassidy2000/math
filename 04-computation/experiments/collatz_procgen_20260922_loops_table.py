#!/usr/bin/env python3
"""collatz_procgen_20260922_loops_table.py -- print the loop table (s <= SMAX) for the note/.out.

Input: MINH  file of minimal-height loops ('s=.. k:..', from loops_recon1 s H(s)),
       FAM   file of the splicing family ('s=.. <how> k:..', from loops_family.py with record cycles).
For each s: K0(s), ratio c(s) = 2^K0/3^s, eps(s+1) = log2(c/1.5), minimal height H(s), climb length,
the minimal-height loop as a forward E-word (run-length: M^a = a multiplications n->3n+1 in a row,
H^b = b halvings) and its reverse move list, and the construction of the family loop.
Usage: python3 ..._loops_table.py MINH FAM SMAX
"""
import sys, re, math

def K0(s):
    return (3 ** (s + 1)).bit_length() - 1

def read(fn):
    d = {}
    for line in open(fn):
        m = re.search(r's=\s*(\d+)(.*?)k:(\S+)', line)
        if m:
            d[int(m.group(1))] = (m.group(2).strip(), [int(t) for t in m.group(3).split(',')])
    return d

def fwd_word(ks):
    seq = []
    for i in range(len(ks), 0, -1):
        seq.append('M'); seq.extend('H' * ks[i - 1])
    runs = []
    for ch in seq:
        if runs and runs[-1][0] == ch: runs[-1][1] += 1
        else: runs.append([ch, 1])
    return ' '.join(c if n == 1 else f'{c}{n}' for c, n in runs)

def values(ks):
    xs = [1]
    for k in ks:
        xs.append(((1 << k) * xs[-1] - 1) // 3)
    return xs

def main():
    minh = read(sys.argv[1]); fam = read(sys.argv[2]); S = int(sys.argv[3])
    print("s   K0   c(s)=2^K0/3^s  eps(s+1)  H(s)   climb  family construction   | minimal-height loop: forward word (from 1) ; reverse moves k_1..k_s")
    for s in range(1, S + 1):
        ks = minh[s][1]; K = sum(ks); c = 2 ** K / 3 ** s
        eps = math.log2(c / 1.5)
        xs = values(ks); H = max(xs)
        climb = 1
        for i in range(s, 1, -1):
            if ks[i - 1] == 0: climb += 1
            else: break
        how = fam[s][0].split()[0] if s in fam else ('trivial' if s == 1 else '')
        print(f"{s:<3d} {K:<4d} {c:.6f}      {eps:+.5f}  {H:<6d} {climb:<5d}  {how:<20s}  | {fwd_word(ks)} ; {','.join(map(str, ks))}")

if __name__ == '__main__':
    main()
