"""Finite cap discovery for shared-port birth/seal returns."""
import json,sys
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
HERE=Path(__file__).resolve().parent;PREV=HERE.parent/'lambda_emergent_20261103';sys.path.insert(0,str(PREV))
import effects as e

def pivots(M):return DomainMatrix.from_Matrix(M).rref()[1]
def rank(M):return len(pivots(M))
def close_labeled(queries,ops):
    columns=[q for name,q in queries];labels=[{'word':[],'query':name} for name,q in queries];depth=0
    while True:
        nxt=columns+[M*q for M in ops.values() for q in columns]
        labs=labels+[{'word':[k]+lab['word'],'query':lab['query']} for k in ops for lab in labels]
        ps=pivots(s.Matrix.hstack(*nxt));new=[nxt[i] for i in ps];newlabels=[labs[i] for i in ps]
        if len(new)==len(columns):return s.Matrix.hstack(*new),newlabels,depth
        columns,labels=new,newlabels;depth+=1

def main():
    base={'3':e.macro(3),'5':e.macro(5),'N':e.macro(0)}
    q=e.named_queries();queries=[(k,q[k]) for k in ('h','s','E1','S2')]
    seals=[base['N'],e.compose([(0,1,1),(0,2,2)],True),e.compose([(0,1,1),(0,2,0)],True)]
    formal=s.Matrix.vstack(*seals);BF=formal[list(pivots(formal.T)),:]
    # The empty-cover coordinate is zero on every nonempty physical source.
    # Keeping it as an arbitrary formal input creates a two-singleton ghost.
    empty=e.STATES.index(((),0))
    image=s.Matrix.vstack(*[M[[i for i in range(M.rows) if i!=empty],:] for M in seals])
    B=image[list(pivots(image.T)),:];e.need(B.rows==31,'nonempty sealed physical source span')
    cases=[]
    for masks in ((3,),(0,),(1,),(2,),(0,1,2,3)):
        ops={**base,**{'J'+str(m):e.dense((0,m,2))*base['N'] for m in masks}}
        Q,labels,depth=close_labeled(queries,ops);effective=B*Q
        columns=list(pivots(effective));R=Q[:,columns];V=B*R
        print('masks',masks,'universal rank',Q.cols,'sealed rank',V.cols,'depth',depth,flush=True)
        transforms={k:V.gauss_jordan_solve(B*M*R)[0] for k,M in ops.items()}
        e.need(all(rank(s.Matrix.vstack(B,B*M))==B.rows for M in ops.values()),'physical sealed source space is invariant')
        cases.append({'masks':masks,'universal_rank':Q.cols,'sealed_rank':V.cols,'depth':depth,
                      'formal_rank_with_empty_ghost':rank(BF*Q),'query_labels':[labels[i] for i in columns],
                      'basis':e.data(R),'matrices':{k:e.data(M) for k,M in transforms.items()}})
    (HERE/'joint_probe.json').write_text(json.dumps({'status':'EXACT_OPERATOR_CERTIFICATE','physical_source_basis':e.data(B),'cases':cases},indent=2)+'\n',encoding='utf8',newline='\n')
if __name__=='__main__':main()
