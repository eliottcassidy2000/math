#!/usr/bin/env python3
"""Decode a prescribed count in the entire one-apex/transitive-spine family."""
import argparse

def h(a: int) -> int:
    if a < 0:
        raise ValueError('nonnegative input required')
    total = 1
    while a:
        p = a >> 1
        if not a & 1:
            total += p + (p & 1)
        a = p
    return total

def cut(x: int, y: int, hx=None) -> int:
    return x*y + (h(x) if hx is None else hx) + h(y) - 1 + (x & y & 1)

def decode(n: int):
    if n < 1 or n % 2 == 0:
        return None
    if n == 1:
        return '0'
    imax = (n - 1).bit_length() + 1
    for width in range(2, imax+1):
        right = width//2
        left = width-right
        y_min, y_stop = 1 << (right-1), 1 << right
        for x in range(1 << (left-1), min(1 << left, (n-1)//y_min+1)):
            hx = h(x)
            lo, hi = y_min, y_stop
            while lo < hi:
                mid = (lo+hi)//2
                if cut(x, mid, hx) >= n:
                    hi = mid
                else:
                    lo = mid+1
            if lo < y_stop and cut(x, lo, hx) == n:
                tail = format(lo, 'b')[::-1].translate(str.maketrans('01', '10'))
                return format(x, 'b') + tail
    return None

def path_dp(word: str) -> int:
    size = len(word)+1
    all_rows = [(1 << len(word)) - (1 << (i+1)) for i in range(len(word))]
    all_rows = [r | (1 << len(word) if word[i]=='0' else 0)
                for i, r in enumerate(all_rows)]
    all_rows.append(sum(1 << i for i,b in enumerate(word) if b=='1'))
    states = [[0]*size for _ in range(1 << size)]
    for i in range(size):
        states[1<<i][i] = 1
    for mask, row in enumerate(states):
        for i, v in enumerate(row):
            if not v:
                continue
            avail = all_rows[i] & ~mask
            while avail:
                candidate = avail & -avail
                avail -= candidate
                states[mask | candidate][candidate.bit_length()-1] += v
    return sum(states[-1])

def self_test() -> None:
    found = {1}
    for l in range(2,9):
        for mid in range(1 << (l-2)):
            w = format((1 << (l-1)) | (mid << 1),'b')
            found.add(h(int(w,2)))
            assert path_dp(w) == h(int(w,2))
    # Width <= 8 covers all normalized solutions for n <= 127.
    for n in range(1,128,2):
        w = decode(n)
        assert (w is not None) == (n in found)
        assert w is None or path_dp(w)==n
    assert cut(4,6)==cut(4,7)==33
    print('All path-DP and inversion tests passed.')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('target', type=int, nargs='?')
    p.add_argument('--test',action='store_true')
    a=p.parse_args()
    if a.test:
        self_test()
    if a.target is not None:
        print(decode(a.target))
