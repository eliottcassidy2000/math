#!/usr/bin/env python3
"""Test each of 375 local branch assumptions against the WHOLE selected
six-coordinate additive-colouring CNF. All original vertices/ports remain.
A SAT witness defeats that branch's elimination in this selected geometry;
a UNSAT verdict only eliminates an additive branch, not all colourings.
"""
import argparse,collections,gzip,hashlib,itertools,json
from pathlib import Path
from pysat.formula import CNF
from pysat.solvers import Solver


def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cnf',type=Path,required=True);ap.add_argument('--geometry',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    geom=json.loads(args.geometry.read_text());groups=geom['support_groups'];six=[g for g in groups if g['family']=='u_times_v_square']
    sizes=geom['coordinate_sizes'];offset=[0]
    for n in sizes:offset.append(offset[-1]+n)
    end=6*offset[-1]
    def var(axis,i,c):return 6*(offset[axis]+i)+c+1
    formula=CNF(from_file=str(args.cnf));need(formula.nv==geom['variables'],'CNF variable universe')
    need(len(formula.clauses)==geom['clauses'],'CNF clause count')
    axes=set()
    for C in formula.clauses:
        if len(C)==2 and all(-end<=l<0 for l in C):
            i,j=(-l-1 for l in C)
            if i%6==j%6:axes.add(tuple(sorted((i//6,j//6))))
    selectors={};top=formula.nv
    for group in six:
        for axis,ps in group['factors']:
            for p in (2,3):
                for a in range(p):
                    key=(axis,tuple(map(tuple,ps)),p,a)
                    if key in selectors:continue
                    top+=1;selectors[key]=top
                    for i,j in ps:
                        for x in range(6):
                            for y in range(6):
                                if (y-x)%p!=a:formula.append([-top,-var(axis,i,x),-var(axis,j,y)])
    results=[];witnesses={};checks=0
    with Solver(name='g3',bootstrap_with=formula.clauses) as sol:
        for gi,group in enumerate(six):
            rows=[]
            for i,j in itertools.combinations(range(6),2):
                for p in (2,3):
                    for a in range(p):
                        for q in (2,3):
                            for b in range(q):
                                branch=[i,p,a,j,q,b];ass=[]
                                for pos,mod,res in ((i,p,a),(j,q,b)):
                                    axis,ps=group['factors'][pos];ass.append(selectors[(axis,tuple(map(tuple,ps)),mod,res)])
                                if not sol.solve(assumptions=ass):
                                    rows.append({'branch':branch,'status':'UNSAT','core':sol.get_core()});continue
                                model={l for l in sol.get_model() if 0<l<=end}
                                colour=[next(c for c in range(6) if 6*k+c+1 in model) for k in range(offset[-1])]
                                digest=hashlib.sha256(bytes(colour)).hexdigest()
                                if digest not in witnesses:
                                    need(all(colour[u]!=colour[v] for u,v in axes),'axis edge failure')
                                    col=[colour[offset[k]:offset[k+1]] for k in range(6)]
                                    for g in groups:
                                        S={0}
                                        for axis,ps in g['factors']:
                                            D={(col[axis][v]-col[axis][u])%6 for u,v in ps}
                                            S={(x+y)%6 for x in S for y in D}
                                        need(0 not in S,'contact failure')
                                    witnesses[digest]=col;checks+=1
                                col=witnesses[digest]
                                for pos,mod,res in ((i,p,a),(j,q,b)):
                                    axis,ps=group['factors'][pos]
                                    need(all((col[axis][v]-col[axis][u])%mod==res for u,v in ps),'branch failure')
                                rows.append({'branch':branch,'status':'SAT','witness':digest})
            need(len(rows)==375,'branch count')
            counts=dict(collections.Counter(x['status'] for x in rows));print(gi,counts,flush=True)
            results.append({'group':gi,'counts':counts,'branches':rows})
            partial={'status':'FINITE-EXACT solver verdicts; SAT witnesses directly checked','cnf_sha256':hashlib.sha256(args.cnf.read_bytes()).hexdigest(),'coordinate_sizes':sizes,'original_axis_edges':len(axes),'groups':results,'witness_checks':checks,'witnesses':witnesses,'scope':'Each local branch tested against all selected-family and axis constraints; original ports unquotiented. No complete ambient unit-distance census.'}
            payload=(json.dumps(partial,separators=(',',':'))+'\n').encode()
            args.out.write_bytes(gzip.compress(payload,mtime=0) if args.out.suffix=='.gz' else payload)
if __name__=='__main__':main()
