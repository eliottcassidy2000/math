"""Literal interval decomposition and transport to the incoming norm chart."""
from bulk import *
from escape_lifts import source4
from return_program import induced
from sympy import Matrix,symbols,expand
import json


def norm(s,u,d):return d**3+2*d*d*s+d*d*u+2*d*s*s-2*d*s*u-d*u*u+2*s**3-2*s*u*u+u**3


def run():
    s,u,d=symbols('s u d')
    need(expand(norm(u,s+u+d,s)-norm(s,u,d))==0,'incoming cubic norm invariant')
    M=Matrix([[1,1,0],[1,2,1],[3,3,2]])
    C=Matrix([[0,1,0],[1,1,1],[1,0,0]])
    L=Matrix([[0,1,0],[0,0,1],[1,1,1]])
    need(M*C==L*M and M.det()==2,'exact index-two guard observation transport')
    rows=[]
    for kind,h in [('straight',(291,485,873)),('revisit',(255,487,889))]:
        v=tuple(M.inv()*Matrix(h));need(all(x.q==1 for x in v),'guard lattice membership')
        v=tuple(map(int,v));chain=[v]
        while v[1]>=v[0]+v[2]:
            s,u,d=v;v=(d,s,u-s-d);chain.append(v)
        rows.append(dict(kind=kind,H_triple=h,guard_coordinates=chain[0],descent=chain,norm=int(norm(*v))))
    checks=0
    def interval(n):return tuple(sum(1<<j for j in range(n) if j==i+1 or j<=i-2) for i in range(n))
    for m in range(1,5):
      for n in range(1,5):
        T,_=source4(m,1,1,n);U=native.native(T,tuple(range(len(T))))
        blocks=(tuple(range(m+1)),(m+1,m+2),tuple(range(m+3,m+n+4)))
        need(induced(U,blocks[0])==interval(m+1) and induced(U,blocks[2])==interval(n+1),'two literal interval factors')
        need(all(U[a]>>b&1 for i,p in enumerate(blocks) for q in blocks[i+1:] for a in p for b in q),'actual ordinal join')
        checks+=1
    result=dict(norm_transport=rows,observation_matrix=[list(map(int,M.row(i))) for i in range(3)],literal_join_checks=checks)
    (HERE/'incoming.json').write_text(json.dumps(result,indent=2)+'\n')
    print('incoming norm/source integration',result,flush=True)
    return result


if __name__=='__main__':run()
