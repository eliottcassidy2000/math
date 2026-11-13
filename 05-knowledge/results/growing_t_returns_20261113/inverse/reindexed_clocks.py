"""Exact coefficient gates for affine reindexing of a single uniform clock."""
from pathlib import Path
import gzip,json,math,hashlib,sys
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'mrb_seven_cutoffs_20261113'/'module'))
import module_certificates as m

def main():
    src=RESULTS/'multicore_metagraph_20261113'/'inflation_descent.json.gz'
    data=json.loads(gzip.decompress(src.read_bytes()));rows=[];independent=0
    for d in (1,2):
        cs=range(0,6) if d==1 else range(1,5)
        for c in cs:
            k=d+c;h=405*(1<<k)+71;target=405*(1<<(2*d+c))+71
            sources=[r for r in data['records'] if r['H']==h]
            m.need(bool(sources),'nonempty inherited complete ready fibre')
            values=[]
            for row in sources:
                T=tuple(row['rows']);m.need(m.count(T)==h,'independent literal source count');independent+=1
                vals=row['clone2_counts'];m.need(target not in vals,'complete inherited clone2 miss')
                # Recount the nearest candidate without sharing its original
                # weighted inflation counter.
                v=min(range(len(vals)),key=lambda v:abs(vals[v]-target))
                U,_=m.inflate(T,v,2)
                m.need(m.count(U)==vals[v],'independent nearest-clone full count');independent+=1
                values.append(dict(parent_order=row['parent_order'],parent_rank=row['parent_rank'],source_rows=T,clone2_counts=vals,nearest_vertex=v,nearest_H=vals[v]))
            rows.append(dict(d=d,c=c,source_exponent=k,source_H=h,required_clone2_H=target,sources=values))
    # Positivity/normalization controls: q^m=sum_j binom(q,j) K_j(TT_m).
    for q in (2,4,8):
        for length in range(1,10):
            m.need(q**length==sum(math.comb(q,j)*m.kj(length,j) for j in range(1,q+1)),'ordered-cover coefficient normalization')
    m.need(math.factorial(7)%7==0 and 405%7!=0,'uniform q>=8 top factorial obstruction')
    m.need(405%6!=0 and 810%6==0,'q4 phase boundary')
    rec=dict(status='PASS',scope='Inherited all-order ready fibres and clone tables; independently counted source and nearest clone in each row; does not regenerate fibres',source=dict(path=str(src.relative_to(RESULTS.parent.parent)).replace('\\','/'),sha256=hashlib.sha256(src.read_bytes()).hexdigest()),independent_full_counts=independent,excluded_charts=rows,remaining_scope='Only d1,c>=6 or d2,c>=5 can survive these necessary gates; this is not an existence claim')
    suffix='' if __debug__ else '.optimized'
    (HERE/('reindexed_clocks_receipt'+suffix+'.json')).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS',len(rows),'clock phases;',independent,'independent full counts')

if __name__=='__main__':main()
