"""Exact legal-word observation growth and algebraic error annihilators."""
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix
from support import HERE,ROOT,need,save,source,observer
import effects as e
def pivots(M):return DomainMatrix.from_Matrix(M).rref()[1]
def rank(M):return len(pivots(M))
def data(M):return [[int(x) if x.q==1 else str(x) for x in row] for row in M.tolist()]
def modular_rank(M,p=1009):
    rows=[[int(x.p)*pow(int(x.q),-1,p)%p for x in row] for row in M.tolist()];r=0
    for j in range(M.cols):
        pivot=next((i for i in range(r,len(rows)) if rows[i][j]),None)
        if pivot is None:continue
        rows[r],rows[pivot]=rows[pivot],rows[r];inv=pow(rows[r][j],-1,p);rows[r]=[x*inv%p for x in rows[r]]
        for i in range(r+1,len(rows)):
            if rows[i][j]:
                z=rows[i][j];rows[i]=[(a-z*b)%p for a,b in zip(rows[i],rows[r])]
        r+=1
        if r==len(rows):break
    return r
def main():
    d=json.loads((HERE.parent/'return701_algebra_20261103/algebra_exact.json').read_text())
    B=s.Matrix(d['source_basis']);ops={k:s.Matrix(M) for k,M in d['generators'].items()}
    algebra=[s.Matrix(27,27,row) for row in d['basis']];q=e.named_queries();results=[]
    for initial in [('h',),('h','s','E1','S2')]:
        cols=[B*q[k] for k in initial];labels=[{'word':[],'query':k} for k in initial];levels=[];depth=0
        while True:
            Q=s.Matrix.hstack(*cols);ps=pivots(Q);cols=[cols[i] for i in ps];labels=[labels[i] for i in ps];Q=s.Matrix.hstack(*cols)
            constraints=s.Matrix([list(A*Q) for A in algebra]);err=95-rank(constraints)
            need(modular_rank(Q)==Q.cols and modular_rank(constraints)==95-err,'independent modular rank lower certificates')
            levels.append({'depth':depth,'response_rank':Q.cols,'invisible_source_dimension':27-Q.cols,'algebra_error_dimension':err})
            nxt=cols+[M*c for M in ops.values() for c in cols]
            labs=labels+[{'word':[k]+lab['word'],'query':lab['query']} for k in ops for lab in labels]
            keep=pivots(s.Matrix.hstack(*nxt))
            if len(keep)==len(cols):break
            cols=[nxt[i] for i in keep];labels=[labs[i] for i in keep];depth+=1
        need(all(rank(s.Matrix.hstack(Q,M*Q))==Q.cols for M in ops.values()),'all future queries closed')
        entry={'initial_queries':initial,'levels':levels,'basis_labels':labels,'basis':data(Q),
          'scalar_completion_number_for_all_future_linear_responses':Q.cols-levels[0]['response_rank']}
        if Q.cols==27:
            T,r,b,names=source();V=observer(T,(r,b));row=s.Matrix([[V.get(key,0) for key in e.STATES]])
            x=row[:,d['source_pivots']];need(x*B==row,'actual sealed source lies in the declared span')
            reads=x*Q;decoded=reads*Q.inv();need(decoded==x and decoded*B==row,'exact actual-source cap recovery')
            entry['actual6325_decode']={'names':names,'root':r,'spare':b,'context_scalar_responses':data(reads),'recovered_source_coordinates':data(decoded)}
        results.append(entry)
        print('queries',initial,'levels',levels,flush=True)
    save('observer_filtration.json',{'status':'EXACT_RATIONAL_ALL_WORD_CLOSURE','cases':results,
      'scope':'inherited persistent-spare seven-macro language on its 27-dimensional rational source span; word labels are actual macros, arbitrary algebra units are not legal words'})
if __name__=='__main__':main()
