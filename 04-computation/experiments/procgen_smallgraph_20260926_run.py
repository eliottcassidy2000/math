"""procgen_smallgraph_20260926_run.py -- single runner for the smallgraph lane.

Session collatz-procgen-20260922, lane "smallgraph" (2026-09-26).

    cd <worktree>
    python3 -u 04-computation/experiments/procgen_smallgraph_20260926_run.py \
        > 05-knowledge/results/procgen_smallgraph_20260926.out

Every printed claim is a check(...) that raises on failure (lines starting with [OK]).
The C helper is compiled into scratch/procgen_smallgraph/bin/ on first use.
Dependencies: python3 (numpy, scipy, ortools), a C compiler.
"""
import os
import resource
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import procgen_smallgraph_20260926_lib as lib  # noqa: E402
from procgen_smallgraph_20260926_lib import check, sha256_file  # noqa: E402

FILES = ['procgen_smallgraph_20260926_lib.py', 'procgen_smallgraph_20260926_ham.c',
         'procgen_smallgraph_20260926_t1.py', 'procgen_smallgraph_20260926_t2.py',
         'procgen_smallgraph_20260926_t3.py', 'procgen_smallgraph_20260926_fence6.py',
         'procgen_smallgraph_20260926_run.py']


def main():
    t0 = time.time()
    print('# procgen_smallgraph_20260926 -- small graphs encoding arithmetic (square sums, fences, C_n)')
    for f in FILES:
        p = os.path.join(HERE, f)
        check(f'SRC {f}', os.path.exists(p), f'sha256 {sha256_file(p)}')
    lib.ensure_ham()
    import procgen_smallgraph_20260926_t1 as t1
    import procgen_smallgraph_20260926_t2 as t2
    import procgen_smallgraph_20260926_t3 as t3
    import procgen_smallgraph_20260926_fence6 as f6
    for name, fn in [('T0', t1.t0), ('T1.A', t1.t1a), ('T1.B', t1.t1b), ('T1.C', t1.t1c), ('T1.D', t1.t1d), ('T1.E', t1.t1e),
                     ('T1.F', t1.t1f), ('T1.G', t1.t1g), ('T2', t2.t2), ('T2.H', f6.t2h),
                     ('T3.A', t3.t3a), ('T3.BC', t3.t3bc), ('T3.D', t3.t3d), ('T3.E', t3.t3e),
                     ('T3.F', t3.t3f), ('T3.G', t3.t3g)]:
        ts = time.time()
        fn()
        print(f'   ({name} section time {time.time() - ts:.1f} s)', flush=True)
    wall = time.time() - t0
    self_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child_rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    scale = 1 if sys.platform == 'darwin' else 1024      # macOS reports bytes, Linux KiB
    mb_self, mb_child = self_rss * scale / 2 ** 20, child_rss * scale / 2 ** 20
    check('RUN resources', mb_self < 700 and mb_child < 700,
          f'wall {wall:.0f} s; peak RSS {mb_self:.0f} MB (python, incl. CP-SAT), {mb_child:.0f} MB (largest C child); '
          f'{lib.N_CHECKS[0] + 1} checks passed')


if __name__ == '__main__':
    main()
