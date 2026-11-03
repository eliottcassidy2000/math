"""Exact joint-observer discovery for actual balanced-source macro returns."""
from pathlib import Path
import sys,json
from fractions import Fraction
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PREV=HERE.parent/'rooted_carry_closure_20261103';sys.path.insert(0,str(PREV))
import phase_recurrence as cap
from observer import old,observer,need,junction

STATES=cap.STATES;DIM=len(STATES)
def dense(command):
    M=s.zeros(DIM)
    for i,row in enumerate(cap.matrix(command)):
        for j,c in row:M[i,j]+=c
    return M

def macro(j):
    if j==0:commands=[(0,1,0),(0,1,2)];swap=False
    else:
        commands=[(0,1,1),(0,2,0)]+[(0,2,2)]*(j-1)+[(0,1,1),(0,2,0)];swap=True
    return compose(commands,swap)

def compose(commands,swap=False):
    M=s.eye(DIM)
    for command in commands:M=M*dense(command)
    if swap:
        P=s.zeros(DIM)
        for i,(sk,m) in enumerate(STATES):
            new=tuple(sorted((1-a if a>=0 else a,1-b if b>=0 else b) for a,b in sk))
            P[i,cap.INDEX[new,m]]=1
        M=M*P
    return M

def named_queries():
    out={name:s.zeros(DIM,1) for name in ('h','s','q','F','E1','S2','P2','D2','E3')}
    for i,(sk,m) in enumerate(STATES):
        k=len(sk)+m
        out['h'][i]=int(k==1)
        out['s'][i]=int(k==1 and any(a==0 for a,b in sk))
        out['q'][i]=int(k==2 and (0,0) in sk)
        out['F'][i]=int(k in (1,2) and any(b==0 for a,b in sk))
        out['E1'][i]=int(k==1 and any(b==0 for a,b in sk))
        out['S2'][i]=int(k in (1,2) and any(a==0 for a,b in sk))
        out['P2'][i]=2 if k in (1,2) else 0
        out['D2'][i]=2 if k in (2,3) and (0,0) in sk else 0
        out['E3'][i]=(1 if k==1 else 2 if k in (2,3) else 0)*int(any(b==0 for a,b in sk))
    return out

def basis(columns):
    if not columns:return []
    Q=s.Matrix.hstack(*columns);return [columns[i] for i in Q.rref()[1]]

def close(queries,operators):
    B=basis(queries);depth=0
    while True:
        nxt=basis(B+[M*q for M in operators for q in B])
        if len(nxt)==len(B):return s.Matrix.hstack(*B),depth
        B=nxt;depth+=1

def data(M):return [[int(x) if x.q==1 else str(x) for x in row] for row in M.tolist()]

def main():
    ops={'3':macro(3),'5':macro(5),'N':macro(0)};q=named_queries()
    Q,depth=close([q[name] for name in ('h','q','E1','S2')],list(ops.values()))
    print('universal observable dimension',Q.cols,'depth',depth,flush=True)
    transforms={name:Q.gauss_jordan_solve(M*Q)[0] for name,M in ops.items()}
    # All macro endpoints satisfy balanced linear equalities, so quotient the
    # source space by the joint image of the universal physical macros.
    seals=[ops['N'],compose([(0,1,1),(0,2,2)],True),compose([(0,1,1),(0,2,0)],True)]
    image=s.Matrix.vstack(*seals);pivots=image.T.rref()[1];B=image[list(pivots),:]
    need(all(s.Matrix.vstack(B,B*M).rank()==B.rows for M in ops.values()),'sealed source span is invariant under every macro')
    need(all(M*(q['s']-q['q'])==s.zeros(DIM,1) and M*(q['F']-q['h'])==s.zeros(DIM,1) for M in ops.values()),'universal balanced outputs')
    effective=B*Q;rank=effective.rank();print('macro image rank',B.rows,'effective observable rank',rank,flush=True)
    natural=s.Matrix.hstack(*[q[name] for name in q]);nat= B*natural
    natural_pivots=nat.rref()[1]
    print('natural rank',len(natural_pivots),'basis',[list(q)[i] for i in natural_pivots],flush=True)
    R=s.Matrix.hstack(*[q[list(q)[i]] for i in natural_pivots])
    if (B*R).rank()==rank and (B*s.Matrix.hstack(R,Q)).rank()==rank:
        reduced={name:(B*R).gauss_jordan_solve(B*M*R)[0] for name,M in ops.items()}
        print('closed natural matrices',{name:data(M) for name,M in reduced.items()},flush=True)
    else:reduced={}
    out={'status':'EXACT_OPERATOR_CERTIFICATE','states':DIM,'observable_dimension':Q.cols,'macro_image_dimension':B.rows,
         'effective_observable_dimension':rank,'natural_queries':list(q),'natural_basis':[list(q)[i] for i in natural_pivots],
         'natural_matrices':{name:data(M) for name,M in reduced.items()},'universal_observer_basis':data(Q),
         'universal_matrices':{name:data(M) for name,M in transforms.items()}}
    (HERE/'effects.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
if __name__=='__main__':main()
