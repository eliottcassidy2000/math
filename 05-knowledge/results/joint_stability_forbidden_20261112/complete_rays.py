"""Literal constructors for two infinite coefficient families with every head.

Run without options to audit. To emit a tournament use, for example,
--family crossing --parameter 8 --exponent 1.
"""
import argparse,json
import numpy as np
from census_nine import need,save,tables,defect,ib
from joint_nine import cones
from construct import tail

FIVE=[2,4,9,19,7]
NINETY_NINE=[42,44,9,16,39,76,31]
EXCEPTION=[320,25,35,389,141,27,62,103,246]
HEAD1=[6,200,50,293,267,83,157,61,231]
HEAD2=[94,156,24,16,32,79,158,61]
Q=[2,100,49,7,75,409,173,31,223]
S1039=[4,17,26,19,673,15,63,111,255,495]

def join(A,B):
    n=len(A);return [r|(((1<<len(B))-1)<<n) for r in A]+[r<<n for r in B]

def construct(family,parameter,exponent):
    need(parameter>=1 and exponent>=1,'positive family parameter and dyadic exponent')
    a=parameter;s=exponent
    if family=='1039':
        c=49+495*2**a
        if s==1:
            rows=EXCEPTION if a==1 else join(NINETY_NINE,tail(FIVE,[0,1],a-1))
            method='separate2079 witness' if a==1 else '99 times coefficient5 ray'
        else:
            rr,oo,_=cones.expand(S1039,[2,3,4],[5],[a]);rows=tail(rr,oo,s-1);method='unit stable-port family'
    elif family=='crossing':
        c=5841*2**a-241
        if s==1:rows=tail(HEAD1,[3,6,7],a);method='defect-481 endpoint interface'
        elif s==2:rows=join(FIVE,tail(HEAD2,[1,6,7],a+1));method='9 times defect-107 endpoint interface'
        else:
            rr,oo,_=cones.expand(Q,[1,4,5],[0,8],[2,a+2]);rows=tail(rr,oo,s-2);method='unit crossing line'
    else:raise ValueError(family)
    return {'family':family,'parameter':a,'exponent':s,'c':c,'H':1+c*2**s,'order':len(rows),'method':method,'rows':rows}

def audit():
    interfaces=[]
    for rows,O,D,U in [(HEAD1,[3,6,7],-481,11682),(HEAD2,[1,6,7],-107,1298),(FIVE,[0,1],1,20)]:
        h,u,w,_=tables(np.array(rows,dtype=np.int64));need(defect(h,w,O)==D and sum(int(u[i]) for i in O)==U,'compiled endpoint interface')
        ih,iu,iw,*_=ib.tc.path_tables(rows);need(ih==int(h) and iu==u.tolist() and iw==w.tolist(),'independent full Python endpoint tables')
        interfaces.append({'rows':rows,'selected':O,'H':int(h),'defect':D,'U':U})
    need(int(tables(np.array(FIVE,dtype=np.int64))[0])==9,'literal factor9')
    controls=[]
    for family in ['1039','crossing']:
        for a in range(1,7):
            for s in range(1,7):
                r=construct(family,a,s)
                if r['order']<=18:
                    h,_,_,_=tables(np.array(r['rows'],dtype=np.int64));need(int(h)==r['H'],'literal entire constructed tournament')
                    controls.append({k:v for k,v in r.items() if k!='rows'})
    result={'scope':'All positive parameter and exponent values for c=49+495*2^a and c=5841*2^a-241. Full lower-head constructions, with actual adjacency constructors. Algebraic all-length tail and SCC-product proofs; literal entire graphs for the declared grid where order<=18. No claim of every odd coefficient.',
            'interfaces':interfaces,'literal_controls':controls,'all_heads_closed_for_both_families':True,
            'crossing_head1':'-481+11682*2^k=1+2*(5841*2^k-241), tail length k>=1',
            'crossing_head2':'9*(-107+1298*2^(k+1))=1+4*(5841*2^k-241), tail length k+1>=2',
            '1039_head':'99*(1+20*2^(a-1)) for a>=2; separate2079 at a=1'}
    save('complete_rays.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='literal_controls'}),flush=True);print('literal whole-graph controls',len(controls),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--family',choices=['1039','crossing']);p.add_argument('--parameter',type=int,default=1);p.add_argument('--exponent',type=int,default=1);a=p.parse_args()
    if a.family:print(json.dumps(construct(a.family,a.parameter,a.exponent),indent=2))
    else:audit()
