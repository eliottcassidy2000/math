#!/usr/bin/env python3
"""Bounded Schur-complement repair of self-complementary weighted modes.

Complete Eulerian and switching-class universes through n=7. A tested
free-pair block is an extra obligation, never inferred from weighted roots.
"""
from collections import deque
import hashlib
import json
import sympy as sp
from sympy.polys.matrices import DomainMatrix
import eulerian_boolean_kernel_overnight_hexagon_sep05 as native
import even_graph_triangle_quotient_spectrum_thm4078_independent_audit as dual


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def exact_rank(matrix):
    if not matrix.rows or not matrix.cols:
        return 0
    return len(DomainMatrix.from_Matrix(matrix).convert_to(sp.QQ).rref()[1])


def apply_linear(mask,columns):
    answer=0
    while mask:
        bit=mask & -mask
        answer ^= columns[bit.bit_length()-1]
        mask ^= bit
    return answer


def switching_orbits(n):
    edges,index=dual.gauge_edges(n)
    m=len(edges)
    maps=[]
    for v in range(n-1):
        perm=list(range(n))
        perm[v],perm[v+1]=perm[v+1],perm[v]
        maps.append([dual.regauge_after_permutation(1<<j,n,tuple(perm),edges,index)
                     for j in range(m)])
    labels=[-1]*(1<<m)
    orbits=[]
    for start in range(1<<m):
        if labels[start]>=0:
            continue
        label=len(orbits)
        labels[start]=label
        queue=deque([start])
        members=[]
        while queue:
            state=queue.popleft()
            members.append(state)
            for columns in maps:
                target=apply_linear(state,columns)
                if labels[target]<0:
                    labels[target]=label
                    queue.append(target)
                else:
                    require(labels[target]==label,"complete dual generator orbit")
        orbits.append(members)
    complements=[labels[members[0]^((1<<m)-1)] for members in orbits]
    require(all(complements[complements[i]]==i for i in range(len(orbits))),"complement orbit involution")
    return orbits,complements,edges


def main():
    report=[]
    for n in range(3,8):
        reps,orbit,sizes,index=native.generator_orbits(n)
        orbits,complement,gauge_edges=switching_orbits(n)
        require(len(orbits)==len(reps),"complete primal/dual count")
        compressed=[sum(((rep >> index[edge])&1)<<j for j,edge in enumerate(gauge_edges)) for rep in reps]
        psi=sp.Matrix([[sum(1 if (h & x).bit_count()%2==0 else -1 for h in members)
                        for members in orbits] for x in compressed])
        require(psi.det(method="domain-ge")!=0,"complete invariant Fourier basis")
        selfdual=[i for i,j in enumerate(complement) if i==j]
        pairs=[(i,j) for i,j in enumerate(complement) if i<j]
        even=[i for i,r in enumerate(reps) if r.bit_count()%2==0]
        odd=[i for i,r in enumerate(reps) if r.bit_count()%2]
        require(len(selfdual)==len(even)-len(odd) and len(pairs)==len(odd),"parity basis dimensions")
        fixed=psi.extract(even,selfdual)
        free=sp.Matrix([[psi[i,a]+psi[i,b] for a,b in pairs] for i in even])
        odd_basis=sp.Matrix([[psi[i,a]-psi[i,b] for a,b in pairs] for i in odd])
        require(fixed.row_join(free).det(method="domain-ge")!=0 and odd_basis.det(method="domain-ge")!=0,"square parity Fourier bases")
        triangles=native.old.simple_cycles(n,3,index)
        M=sp.Matrix(native.old.weighted_operator(reps,orbit,triangles))
        B=sp.Matrix(native.old.boolean_support(M.tolist()))
        C=B.extract(odd,even)
        require(M.extract(odd,even)*fixed==sp.zeros(len(odd),len(selfdual)),"weighted self-complementary zero modes")
        defect=C*fixed
        block=C*free
        rank=exact_rank(block)
        row={"n":n,"q":len(reps),"fixed_modes":len(selfdual),
             "free_pairs":len(pairs),"free_block_rank":rank,
             "free_block_determinant":str(block.det(method="domain-ge")),
             "unrepaired_boolean_defect_rank":exact_rank(defect)}
        if rank==len(pairs):
            coefficients=(DomainMatrix.from_Matrix(block).convert_to(sp.QQ).inv().to_Matrix()*defect
                          if selfdual else sp.zeros(len(odd),0))
            repaired=fixed-free*coefficients
            require(C*repaired==sp.zeros(len(odd),len(selfdual)),"actual Boolean repaired modes")
            require(exact_rank(repaired)==len(selfdual),"repaired fixed modes independent")
            require(exact_rank(C)==len(odd),"free block forces saturation")
            row["repair_coefficient_denominator_lcm"]=str(sp.ilcm(1,1,*[sp.denom(v) for v in coefficients]))
            if n<=6:
                row["repair_coefficients"]=[[str(v) for v in coefficients.row(i)] for i in range(coefficients.rows)]
        else:
            row["free_block_hostile"]=True
        print("EXACT FOURIER REPAIR",json.dumps(row,sort_keys=True),flush=True)
        report.append(row)
    print("SEMANTIC SHA256",hashlib.sha256(json.dumps(report,sort_keys=True).encode()).hexdigest())
    print("FINITE-EXACT n3..7; general free-pair transversality OPEN")


if __name__=="__main__":
    main()
