"""Small library for one-apex counts g (S84 exploration)."""
import numpy as np


def gtable(N):
    """g on [0, 2^N): g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod 2)."""
    g = np.zeros(1 << N, dtype=np.int64)
    g[0] = 1
    size = 1
    for _ in range(N):
        v = np.arange(size, dtype=np.int64)
        g[2 * v + 1] = g[v]
        g[2 * v] = g[v] + v + (v & 1)
        size *= 2
    g[0] = 1
    return g


def rc(x, n):
    """reversed complement of the n-bit word x."""
    r = 0
    for i in range(n):
        r |= (1 - ((x >> i) & 1)) << (n - 1 - i)
    return r


def rc_table(n):
    x = np.arange(1 << n, dtype=np.int64)
    r = np.zeros_like(x)
    for i in range(n):
        r |= (1 - ((x >> i) & 1)) << (n - 1 - i)
    return r
