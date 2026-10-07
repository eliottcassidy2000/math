#!/usr/bin/env python3
"""Re-solve a DIMACS CNF with a pysat solver under a time limit (mac-mini-2026-10-06-oaimath2).
Used on the certified Q_gap(3,7) clause set (oai2_20261006_readers/hindman_ramsey/finv_unsat_n3_t7.cnf.gz, gunzipped):
  python3 oai2_20261006_finv37_resolve.py cadical153 7000 finv_unsat_n3_t7.cnf
  -> solver=cadical153 vars=1098 clauses=669086 result=UNSAT time=1379s   (see oai2_20261006_finv37_resolve.out)
(Lingeling cannot be run this way: pysat does not support limited solves for it.)
"""
import sys, time, threading
from pysat.formula import CNF
from pysat.solvers import Solver
name, limit = sys.argv[1], float(sys.argv[2])
cnf = CNF(from_file=sys.argv[3])
t0 = time.time()
with Solver(name=name, bootstrap_with=cnf.clauses) as s:
    timer = threading.Timer(limit, s.interrupt); timer.start()
    r = s.solve_limited(expect_interrupt=True)
    timer.cancel()
    print(f"solver={name} vars={cnf.nv} clauses={len(cnf.clauses)} result={'UNSAT' if r is False else ('SAT' if r else 'UNKNOWN (interrupted)')} time={time.time()-t0:.0f}s", flush=True)
