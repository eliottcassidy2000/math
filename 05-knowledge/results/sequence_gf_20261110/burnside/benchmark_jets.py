#!/usr/bin/env python3
"""Reproduce larger fixed-order jets, not a full-count timing comparison."""
from time import perf_counter
from cluster_compiler import ClusterCompiler, reconstruct, burnside_control

for n, r in [(80, 3), (200, 2)]:
    start = perf_counter()
    compiler = ClusterCompiler(n)
    layers = [compiler.log_layer(j) for j in range(r + 1)]
    b = reconstruct(layers)
    elapsed = perf_counter() - start
    reference, _ = burnside_control(18, r)
    assert [row[:19] for row in b] == reference
    assert all(x >= 0 for row in layers for x in row)
    print(f'N={n} order={r}: compiler+recurrence seconds={elapsed:.3f}; '
          'all layers agree with independent partitions through n=18.')
    print(f'Final coefficient [x^{n} s^{r}] A(x,1+s) = {b[r][n]}')
