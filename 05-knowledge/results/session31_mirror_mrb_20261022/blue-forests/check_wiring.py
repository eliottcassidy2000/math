#!/usr/bin/env python3
"""Independent explicit controls for inherited leg/parity laws and a new forest witness."""
from itertools import permutations
from pathlib import Path
import json
from probe import edge_orbits,get_lap,forest_coeff

def energy(L):return sum((2*L[i][i]-len(L)+1)**2 for i in range(len(L)))
def c3(L):
    n=len(L)
    return sum((L[i][j]<0 and L[j][k]<0 and L[k][i]<0) or
               (L[j][i]<0 and L[k][j]<0 and L[i][k]<0)
               for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
def isomorphism(L,M):
    n=len(L)
    for p in permutations(range(n)):
        if any(L[i][i]!=M[p[i]][p[i]] for i in range(n)):continue
        if all(L[i][j]==M[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j):return p
    return None

def controls(n,orbs,expected_blue=False):
    full=(1<<len(orbs))-1;blue=True
    for mask in range(full+1):
        L=get_lap(n,orbs,mask);M=get_lap(n,orbs,mask^full)
        W,V=energy(L),energy(M)
        residue=n if n%2==0 else 0
        assert (W-residue)%8==0
        flux=L[0][0]-L[-1][-1]-1
        assert W-V==8*flux
        C,E=c3(L),c3(M)
        assert E-C==flux
        assert W==n*(n*n-1)//3-8*C
        assert ((W-residue)//8-(V-residue)//8-flux)%2==0
        if expected_blue:
            assert flux==2*L[0][0]-n
            assert (E-C)%2==n%2
        blue &= all(L[n-1-i][n-1-j]==L[j][i] for i in range(n) for j in range(n) if i!=j)
    if expected_blue:assert blue
    return full+1

out={'status':'PASS','inherited_law_checks':[],'zero_flux_loop_checks':[]}
for n in range(3,10):
    edges,orbs=edge_orbits(n)
    nb=controls(n,orbs,expected_blue=True)
    row={'n':n,'blue_presentations':nb}
    if n<=7:row['all_presentations']=controls(n,[[e] for e in edges])
    out['inherited_law_checks'].append(row)
for n in (4,6):
    edges,orbs=edge_orbits(n);full=(1<<len(orbs))-1;candidates=loops=0
    witnesses=[]
    for m in range(full+1):
        L=get_lap(n,orbs,m)
        if L[0][0]!=n//2:continue
        candidates+=1;M=get_lap(n,orbs,m^full);p=isomorphism(L,M)
        loops+=p is not None
        if p is None and not witnesses:
            witnesses.append({'orbit_mask':m,'orbits':orbs,'adjacency':[[int(x<0) for x in row] for row in L],
                              'complement_adjacency':[[int(x<0) for x in row] for row in M],
                              'outdegrees':[L[i][i] for i in range(n)],'complement_outdegrees':[M[i][i] for i in range(n)],
                              'forest_coefficients':forest_coeff(L),'complement_forest_coefficients':forest_coeff(M),
                              'energy':energy(L),'cyclic_triangles':c3(L),'isomorphic':False})
    assert (candidates,loops)==((2,2) if n==4 else (24,4))
    out['zero_flux_loop_checks'].append({'n':n,'candidates':candidates,'loop_endpoints':loops,'first_nonloop':witnesses})
(HERE:=Path(__file__).resolve().parent).joinpath('wiring_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
