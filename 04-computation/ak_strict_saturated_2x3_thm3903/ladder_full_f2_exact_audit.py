#!/usr/bin/env python3
"""Exact determinantal saturation audit of 2+4 full-ladder traces."""

from itertools import combinations, combinations_with_replacement
import sys
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

from ladder_all_cost10_kernel_search import make_topology
from ladder_full_f4_symbolic_audit import generator_matrix
from ladder_full_f3_exact_audit import full_cycle_matrix


def rank_at_most_two_minors(matrix, rest):
    result=[]
    for rows in combinations(range(4),3):
        for cols in combinations(rest,3):
            poly=sp.factor(matrix[list(rows),list(cols)].det())
            if poly!=0 and poly not in result:
                result.append(poly)
    return tuple(result)


def main():
    cases=generic_rank_fail=unit_saturations=potential=0
    potentials=[]
    for positions in combinations_with_replacement(range(6),3):
        topology=make_topology((0,1,2,3),positions)
        if topology is None: continue
        params,A=generator_matrix(topology)
        W=full_cycle_matrix(topology,params)
        for fire in combinations(range(6),2):
            cases+=1
            rest=tuple(v for v in range(6) if v not in fire)
            cols=tuple(2*v+j for v in rest for j in (0,1))
            subs={params[0]:0}
            Ws=W.subs(subs);Ts=A[:,cols].subs(subs)
            if Ws.to_DM().to_field().rank()<4 or Ts.to_DM().to_field().rank()<8:
                generic_rank_fail+=1;continue
            equations=rank_at_most_two_minors(Ws,rest)
            initial_minors=[]
            for column_set in combinations(range(6),4):
                determinant=sp.factor(Ws[:,list(column_set)].to_DM().det().as_expr())
                if determinant!=0 and determinant not in initial_minors:
                    initial_minors.append(determinant)
            necessary_initial=initial_minors[0]
            for poly in initial_minors[1:]:
                necessary_initial=sp.gcd(necessary_initial,poly)
            terminal_minors=[]
            for row_set in combinations(range(10),8):
                determinant=sp.factor(Ts[list(row_set),:].to_DM().det().as_expr())
                if determinant!=0 and determinant not in terminal_minors:
                    terminal_minors.append(determinant)
            necessary_terminal=terminal_minors[0]
            for poly in terminal_minors[1:]:
                necessary_terminal=sp.gcd(necessary_terminal,poly)
            variables=sorted(set().union(
                *(poly.free_symbols for poly in equations+(necessary_initial,necessary_terminal))),key=str)
            sat=sp.symbols("sat")
            gb=sp.groebner(list(equations)+[1-sat*necessary_initial*necessary_terminal],
                           sat,*variables,order="grevlex")
            if gb==[1]:
                unit_saturations+=1
            else:
                surviving_pair=None
                for initial_minor in initial_minors:
                    if surviving_pair is not None: break
                    for terminal_minor in terminal_minors:
                        pair_vars=sorted(set().union(
                            *(poly.free_symbols for poly in equations+(initial_minor,terminal_minor))),key=str)
                        pair_gb=sp.groebner(list(equations)+[1-sat*initial_minor*terminal_minor],
                                            sat,*pair_vars,order="grevlex")
                        if pair_gb!=[1]:
                            surviving_pair=(str(initial_minor),str(terminal_minor),
                                            tuple(str(x) for x in pair_gb.polys[:4]))
                            break
                if surviving_pair is None:
                    unit_saturations+=1
                else:
                    potential+=1
                    if len(potentials)<20:potentials.append((positions,fire,len(equations),surviving_pair))
    print("LADDER_FULL_F2_EXACT_SATURATION_AUDIT")
    print(f"cases={cases};generic_rank_fail={generic_rank_fail};unit_saturations={unit_saturations};potential={potential}")
    print(f"potentials={tuple(potentials)}")
    print("status=F2_BRANCH_IMPOSSIBLE" if potential==0 else "status=POTENTIALS_REMAIN")


if __name__=="__main__":main()
