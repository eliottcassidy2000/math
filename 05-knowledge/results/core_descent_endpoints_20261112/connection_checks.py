"""Exact unit reduction, orientation gauge, and inherited small-number controls."""
from context import *
from itertools import permutations
from math import isqrt,comb
from fractions import Fraction
import importlib.util

def mul(x,y,p,trace,norm):
    a,b=x;c,d=y
    return ((a*c-norm*b*d)%p,(a*d+b*c+trace*b*d)%p)

def power(x,n,p,trace,norm):
    result=(1,0)
    while n:
        if n&1:result=mul(result,x,p,trace,norm)
        x=mul(x,x,p,trace,norm);n//=2
    return result

def run():
    unit=(0,1);powers=[power(unit,k,13,127,1) for k in range(8)]
    need(powers[7]==(1,0) and len(set(powers[:7]))==7,'unit order exactly7')
    squares={x*x%13 for x in range(13)};need(5 not in squares,'irreducible minimal polynomial modulo13')
    need(power(unit,13,13,127,1)==powers[6],'Frobenius equals conjugate inverse')
    residues={1,2,4};rows=[sum(1<<b for b in range(7) if (b-a)%7 in residues) for a in range(7)]
    autos=[];antis=[];paths=0
    for perm in permutations(range(7)):
        signs={bool(rows[a]>>b&1)==bool(rows[perm[a]]>>perm[b]&1) for a in range(7) for b in range(a+1,7)}
        if signs=={True}:autos.append(perm)
        if signs=={False}:antis.append(perm)
        if all(rows[a]>>b&1 for a,b in zip(perm,perm[1:])):paths+=1
    need(len(autos)==21 and len(antis)==21,'complete auto/anti action census')
    need(tuple((-a)%7 for a in range(7)) in antis,'quadratic conjugation reverses tournament')
    affine={tuple((s*a+t)%7 for a in range(7)) for s in residues for t in range(7)}
    need(set(autos)==affine,'exact7:3 automorphism group')
    # The older norm-minus-one clock has a different behavior modulo13.
    old_order=next(k for k in range(1,25) if power(unit,k,13,5,-1)==(1,0))
    traces={(2*a+5*b)%13 for a in range(13) for b in range(13)}
    paired={((2*a+5*b)%13,(5*a+27*b)%13) for a in range(13) for b in range(13)}
    need(old_order==12 and len(traces)==13 and len(paired)==169,'inherited12/13 trace sidecar')
    divisors=[d for d in range(1,13) if 12%d==0];shifted=sorted(d+1 for d in divisors)
    bern=[Fraction(1)]
    for n in range(1,13):bern.append(-sum(Fraction(comb(n+1,k))*bern[k] for k in range(n))/(n+1))
    need(shifted==[2,3,4,5,7,13] and bern[12]==Fraction(-691,2730),'divisor-shift and Bernoulli arithmetic')
    sp=importlib.util.spec_from_file_location('old_paths',HERE.parent/'joint_incidence_clocked_frames_20261102/path_kernel.py');mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)
    square=[]
    for n in range(2,26):
        edges,pp=mod.graph_paths(n,True);new=n-1
        end=sum(p[0]==new or p[-1]==new for p in pp)
        square.append({'order':n,'edges':len(edges),'unoriented_paths':len(pp)//2,'directed_paths_ending_at_new_vertex':end,'directed_paths_with_new_vertex_internal':len(pp)-end,'new_neighbors':[a+1 for a,b in edges if b==new]})
    need({r['order']:r['unoriented_paths'] for r in square if r['unoriented_paths']}=={15:1,16:1,17:1,23:3,25:10},'inherited first25 integer nodes')
    out={'scope':'Exact arithmetic and complete seven-vertex permutation checks. The Paley orientation requires an explicit generator and quadratic-residue gauge; it is not canonically supplied by the unmarked trace. Square-sum counts are inherited replays, not new counts.','unit_polynomial_mod13':'x^2+3*x+1','unit_powers':powers,'Paley7_rows':rows,'Hamiltonian_paths':paths,'automorphisms':len(autos),'anti_automorphisms':len(antis),'old_clock_order_mod13':old_order,'old_trace_classes':len(traces),'old_paired_trace_classes':len(paired),'shifted_divisors12':shifted,'Bernoulli12':str(bern[12]),'square_sum_first25':square}
    save('connection_checks.json.gz',out);print({k:v for k,v in out.items() if k!='square_sum_first25'})
if __name__=='__main__':run()
