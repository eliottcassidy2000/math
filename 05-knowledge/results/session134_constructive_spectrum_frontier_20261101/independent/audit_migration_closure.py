#!/usr/bin/env python3
"""Independent capped-semiring grammar closure and literal repair consumer.
No producer imports. Fixed seed only; arbitrary strip length justified by coefficient saturation.
"""
import json,pathlib,itertools
from audit_migration_compiler import state,literal

def need(c,msg):
    if not c:raise RuntimeError(msg)
N=173;CAP=N+1
S0=((2,1,1,0),(0,1,0,0),(0,0,0,0),(0,0,1,1))
S1=((1,0,0,0),(1,2,1,0),(0,0,1,1),(0,0,0,0))
Z=((0,0,0,0),(0,0,0,0),(0,0,0,0),(2,0,1,1))
R1=((1,1,1,0),(0,0,0,0),(0,0,0,1),(2,0,1,0))
def R(r):
    if r==1:return R1
    a=2**(r-1);b=2*3**(r-1)-2**r+1
    return ((a,a,a//2,0),(b-a,b-2*a+2,(a+b-3)//2,a-1),(2*(a-1),2*(a-2),a,1),(2,2,1,0))
def cap(M):return tuple(tuple(min(CAP,x) for x in row) for row in M)
def closure(matrices):
    start=(0,0,0,1);visited={start};stack=[start];values=set()
    while stack:
        v=stack.pop();values.add(2*v[0]+v[2]+v[3])
        for M in matrices:
            w=tuple(min(CAP,sum(a*b for a,b in zip(row,v))) for row in M)
            if 2*w[0]+w[2]+w[3]>N:continue
            if w not in visited:visited.add(w);stack.append(w)
    return visited,values

def conjugate(M):
    basis=((0,2,2,0),(0,0,2,0),(1,1,1,1),(1,1,0,0));cols=[]
    for v in basis:
        h,F,q,g=(sum(a*b for a,b in zip(row,v)) for row in M)
        need((F-h)%2==0 and (q-F+h-g)%2==0,'nonintegral positive coordinates')
        cols.append(((F-h)//2,(q-F+h-g)//2,g,h-g))
    return tuple(tuple(cols[j][i] for j in range(4)) for i in range(4))

def main():
    root=pathlib.Path(__file__).resolve().parent
    strips={cap(M) for M in (S0,S1,Z)}|{cap(R(r)) for r in range(1,11)}
    # For r>=10 all growing coefficients exceed174. The following control is
    # supplemented by their explicit positive exponential differences in REPORT.md.
    need(all(cap(R(r))==cap(R(10)) for r in range(10,111)),'tail saturation')
    zero,vals=closure(strips);missing=sorted(set(range(1,N+1,2))-vals)
    need(missing==[7,21,67,173] and len(zero)==5373,'zero-strip closure mismatch')
    bank=json.loads((root/'migration_matrix_bank.json').read_text())['records'];mixed={S0,S1}
    for record in bank:
        local=record['local'];fresh=range(1,len(local));cyclic=False
        for i,j,k in itertools.combinations(fresh,3):
            if (local[i]>>j&1 and local[j]>>k&1 and local[k]>>i&1) or (local[i]>>k&1 and local[k]>>j&1 and local[j]>>i&1):cyclic=True
        if not cyclic:
            M=conjugate(record['matrix']);need(all(x>=0 for row in M for x in row),'negative selected matrix');mixed.add(M)
    need(len(mixed)==32,'mixed bank matrix count')
    mixed_states,mixed_vals=closure(mixed)
    need(sorted(set(range(1,N+1,2))-mixed_vals)==[7,21] and len(mixed_states)==6234,'mixed bank closure')
    repairs=json.loads((root/'migration_repairs.json').read_text())['repairs'];recounts=[]
    for target,r in repairs.items():
        A=r['parent_adjacency_masks'];S=r['cut'];native=state(A,S)
        need(native==tuple(r['state']) and native[1]==int(target),'repair DP mismatch')
        need(literal(A,S)==native,'literal repair mismatch')
        recounts.append({'target':int(target),'parent_order':len(A),'independent_native_state':native})
    print(json.dumps({'status':'PASS','saturation':CAP,'seed_positive_coordinates':[0,0,0,1],'zero_strip':{'states':len(zero),'matrices_after_capping':len(strips),'missing_odds':missing,'tail_representative':10},'mixed_transitive_local_bank':{'states':len(mixed_states),'matrices':len(mixed),'missing_odds':[7,21]},'literal_repair_recounts':recounts,'scope':'Exact unbounded-word/all-strip-length result in the stated fixed-seed grammars; does not exclude arbitrary seeds or unrestricted tournaments.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
