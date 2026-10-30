"""Two-port SAT experiment for a saturated (3,3,5;56) coloring.

Default: necessary red-row and green-column counts. --equitable adds
extra row/column assumptions which are NOT necessary for every counterexample.

SAT is a local compatibility witness only, never a 56-vertex coloring.
UNSAT would require a separately checked proof before a Ramsey claim.
"""
from itertools import combinations
from pathlib import Path
from threading import Timer
import argparse,json,hashlib
from pysat.card import CardEnc,EncType
from pysat.formula import CNF,IDPool
from pysat.solvers import Solver

HERE=Path(__file__).resolve().parent


def h13(i,j): return (i-j)%13 in {1,5,8,12}


def build(equitable=False, force_defect=False):
    pool=IDPool(); cnf=CNF()
    def x(i,j,c): return pool.id((i,j,c))
    for i in range(13):
        for j in range(13):
            lits=[x(i,j,c) for c in range(3)]
            cnf.append(lits)
            for a,b in combinations(lits,2): cnf.append([-a,-b])
    for c in [0,1]:
        for i in range(13):
            for lits in [[x(i,j,c) for j in range(13)],[x(j,i,c) for j in range(13)]]:
                is_row = lits[0] == x(i,0,c) and lits[-1] == x(i,12,c)
                if equitable or (c == 0 and is_row) or (c == 1 and not is_row):
                    cnf.extend(CardEnc.equals(lits=lits,bound=4,vpool=pool,encoding=EncType.seqcounter).clauses)
    if force_defect:
        if equitable: raise ValueError('A defect is incompatible with the equitable ansatz')
        cnf.extend(CardEnc.atmost(lits=[x(0,j,1) for j in range(13)],bound=3,vpool=pool,encoding=EncType.seqcounter).clauses)
    # Red triangles use two G vertices, green triangles two R vertices.
    for a,b in combinations(range(13),2):
        if h13(a,b):
            for i in range(13):
                cnf.append([-x(i,a,0),-x(i,b,0)])
                cnf.append([-x(a,i,1),-x(b,i,1)])
    independent={k:[S for S in combinations(range(13),k)
                    if all(not h13(a,b) for a,b in combinations(S,2))] for k in range(1,5)}
    for k in range(1,5):
        for left in independent[k]:
            for right in independent[5-k]:
                cnf.append([-x(i,j,2) for i in left for j in right])
    return cnf,x,independent


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--seconds',type=int,default=45)
    ap.add_argument('--solver',default='g4'); ap.add_argument('--equitable',action='store_true')
    ap.add_argument('--prefix',default='port_probe'); ap.add_argument('--force-defect',action='store_true'); args=ap.parse_args()
    cnf,x,independent=build(args.equitable,args.force_defect)
    dimacs=cnf.to_dimacs().encode()
    (HERE/(args.prefix+'.cnf')).write_bytes(dimacs)
    with Solver(name=args.solver,bootstrap_with=cnf,with_proof=False) as solver:
        timer=Timer(args.seconds,solver.interrupt); timer.start()
        try: result=solver.solve_limited(expect_interrupt=True)
        finally: timer.cancel()
        data={'meaning':'root plus two ports; NOT a full 56-vertex coloring',
              'equitable_extra_assumption':args.equitable,
              'forced_green_row_defect':args.force_defect,
              'status':{True:'SAT',False:'UNSAT',None:'UNKNOWN'}[result],
              'solver':args.solver,'variables':cnf.nv,'clauses':len(cnf.clauses),
              'dimacs_sha256':hashlib.sha256(dimacs).hexdigest(),
              'independent_set_counts':{k:len(v) for k,v in independent.items()}}
        if result:
            model=set(solver.get_model())
            data['cross_colors']=[[next(c for c in range(3) if x(i,j,c) in model) for j in range(13)] for i in range(13)]
        elif result is False:
            data['warning']='UNSAT is not promoted without an independently checked proof.'
    (HERE/(args.prefix+'.json')).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:v for k,v in data.items() if k!='cross_colors'}))


if __name__=='__main__':main()
