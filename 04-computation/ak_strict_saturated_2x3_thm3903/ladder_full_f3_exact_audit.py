#!/usr/bin/env python3
"""Exact algebraic audit of 3+3 traces in the full strict 2x3 ladder.

After the affine slope gauge alpha=0, a three-vertex first wave forces the
two live square-cycle rows to have rank at most one on the three survivors.
Their determinantal varieties have a finite coordinate-zero cover.  On each
cover we impose the full four-row rank-one condition and saturate by the
initial and terminal Gram determinants.  A unit Groebner basis rules out a
rational full-rank initial/terminal pair.
"""

from itertools import combinations, combinations_with_replacement
import hashlib
import json
import sys
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

from ladder_all_cost10_kernel_search import make_topology
from ladder_full_f4_symbolic_audit import generator_matrix
from ladder_full_f3_zero_cover_audit import minimal_zero_covers


def full_cycle_matrix(topology, params):
    return sp.Matrix([
        [sum(sp.Rational(topology[4][i][v][p])*params[p] for p in range(8))
         for v in range(6)]
        for i in range(4)
    ])


def alignment_minors(matrix, rest):
    return tuple(dict.fromkeys(
        sp.expand(matrix[i,u]*matrix[j,v]-matrix[i,v]*matrix[j,u])
        for i,j in combinations(range(4),2)
        for u,v in combinations(rest,2)
        if sp.expand(matrix[i,u]*matrix[j,v]-matrix[i,v]*matrix[j,u]) != 0
    ))


def verify_zero_cover(rest, covers):
    b,c,g,h=sp.symbols("b c g h")
    vals=(b,c,g,h)
    w1=(b,-b,0,-g,g,0)
    w2=(0,c,-c,0,-h,h)
    equations=[sp.expand(w1[i]*w2[j]-w1[j]*w2[i]) for i,j in combinations(rest,2)]
    possible=[]
    z=sp.symbols("cover_z")
    for mask in range(16):
        zero={i for i in range(4) if mask&(1<<i)}
        subs={vals[i]:0 for i in zero}
        eq=[sp.expand(poly.subs(subs)) for poly in equations]
        eq=[poly for poly in eq if poly!=0]
        nonzero=[vals[i] for i in range(4) if i not in zero]
        product=sp.prod(nonzero)
        variables=[z]+nonzero
        ideal=eq+([1-z*product] if product!=1 else [])
        unit=bool(ideal) and sp.groebner(ideal,*variables,order="grevlex")==[1]
        if not unit:
            possible.append(tuple(sorted(zero)))
            if not any(set(cover).issubset(zero) for cover in covers):
                raise RuntimeError(f"zero cover incomplete for rest={rest}, zero={zero}")
    return tuple(possible)


def main():
    cover_certificates={}
    for rest in combinations(range(6),3):
        covers=minimal_zero_covers(rest)
        cover_certificates[rest]=verify_zero_cover(rest,covers)

    branches=generic_rank_fail=unit_saturations=potential=0
    potentials=[]
    for positions in combinations_with_replacement(range(6),3):
        topology=make_topology((0,1,2,3),positions)
        if topology is None:
            continue
        params,A=generator_matrix(topology)
        W=full_cycle_matrix(topology,params)
        for fire in combinations(range(6),3):
            rest=tuple(v for v in range(6) if v not in fire)
            terminal_cols=tuple(2*v+j for v in rest for j in (0,1))
            for cover in minimal_zero_covers(rest):
                branches+=1
                subs={params[0]:0}
                subs.update({params[1+i]:0 for i in cover})
                Ws=W.subs(subs)
                Ts=A[:,terminal_cols].subs(subs)
                if Ws.to_DM().to_field().rank()<4 or Ts.to_DM().to_field().rank()<6:
                    generic_rank_fail+=1
                    continue
                equations=alignment_minors(Ws,rest)
                gram_initial=sp.expand((Ws*Ws.T).det())
                gram_terminal=sp.expand((Ts.T*Ts).det())
                variables=sorted(set().union(
                    *(poly.free_symbols for poly in equations+(gram_initial,gram_terminal))),key=str)
                sat=sp.symbols("sat")
                ideal=list(equations)+[1-sat*gram_initial*gram_terminal]
                gb=sp.groebner(ideal,sat,*variables,order="grevlex")
                if gb==[1]:
                    unit_saturations+=1
                else:
                    potential+=1
                    if len(potentials)<20:
                        potentials.append((positions,fire,cover,len(equations),tuple(str(x) for x in gb.polys[:4])))
    print("LADDER_FULL_F3_EXACT_SATURATION_AUDIT")
    cover_serial = tuple(sorted(cover_certificates.items()))
    cover_hash = hashlib.sha256(json.dumps(cover_serial, separators=(",", ":")).encode()).hexdigest()
    print(f"cover_rest_sets={len(cover_serial)};support_patterns={sum(len(x[1]) for x in cover_serial)};cover_sha256={cover_hash}")
    print(f"branches={branches};generic_rank_fail={generic_rank_fail};unit_saturations={unit_saturations};potential={potential}")
    print(f"potentials={tuple(potentials)}")
    print("status=F3_BRANCH_IMPOSSIBLE" if potential==0 else "status=POTENTIALS_REMAIN")


if __name__=="__main__":main()
