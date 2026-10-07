import sys, time
from pysat.formula import CNF
from pysat.solvers import Solver
cnf=CNF(from_file=sys.argv[1]); name=sys.argv[2]; t0=time.time()
with Solver(name=name,bootstrap_with=cnf.clauses) as s:
    r=s.solve()
print(f"{sys.argv[1]} | solver {name}: {'SAT' if r else 'UNSAT'} in {time.time()-t0:.0f}s ({len(cnf.clauses)} clauses)",flush=True)
