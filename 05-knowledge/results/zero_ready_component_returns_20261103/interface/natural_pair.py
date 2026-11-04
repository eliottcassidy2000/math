"""Exact root-child factor in the universal static17 basis and ten-state controls."""
from pathlib import Path
import hashlib,json
from fractions import Fraction
from gap import count,nine_fields,induced,need
from check import birth,rows_code
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent

def rank(matrix):
    a=[list(map(Fraction,row)) for row in matrix];r=0
    for c in range(len(a[0])):
        i=next((i for i in range(r,len(a)) if a[i][c]),None)
        if i is None:continue
        a[r],a[i]=a[i],a[r];v=a[r][c];a[r]=[x/v for x in a[r]]
        for j in range(len(a)):
            if j!=r and a[j][c]:
                v=a[j][c];a[j]=[x-v*y for x,y in zip(a[j],a[r])]
        r+=1
        if r==len(a):break
    return r

def ten(T,fb,w):
    h,s,q,F,e,P2,S2,E3,Q2=nine_fields(T,fb,w)
    vs=[v for v in range(len(T)) if v!=w]
    Q3=count(induced(T,vs),tuple(vs.index(v) for v in fb if v!=w),3)
    return [h,e,S2,F,Q2,q,P2,count(T,fb,3,starts=((0,w),)),E3,Q3]

def main():
    src=RESULTS/'adaptive_state_induction_20261103/two_port_envelope.json';common=json.loads(src.read_text())
    states=[tuple(map(tuple,s)) for s in common['states']];index={s:i for i,s in enumerate(states)}
    words=common['rank']['suffixes'];columns=[]
    for word in words:
        col=common['observer'][:]
        for mask in reversed(word):
            nxt=[0]*25
            for i,j,c in common['matrices'][3*mask+2]:nxt[i]+=c*col[j]
            col=nxt
        columns.append(col)
    C=[list(x) for x in zip(*columns)]
    old=[((-1,-1),),((-1,0),),((0,-1),),((-1,-1),(0,0)),((-1,-1),(-1,0))]
    B=[[0]*25 for _ in old]
    # Command4: incoming slot0; discard old slot1; newborn takes slot1.
    for k,s in enumerate(old):
        for i,j,c in common['matrices'][4]:
            if i==index[s]:B[k][j]+=c
    swap=[index[tuple(sorted((1-a if a>=0 else a,1-b if b>=0 else b) for a,b in s))] for s in states]
    factors={};cancellation={}
    for name,cols in [('mark_root',C),('root_mark',[C[j] for j in swap])]:
        BC=[[sum(row[i]*cols[i][j] for i in range(25)) for j in range(17)] for row in B]
        e=[-BC[0][j]+BC[1][j]-BC[4][j] for j in range(17)]
        need(e==[0]*17,'forgotten ending-root count cancels on all static17 columns')
        F=[BC[0],[BC[2][j]-BC[0][j] for j in range(17)],
           [BC[3][j]-BC[4][j] for j in range(17)],BC[4]]
        need(rank(F)==4 and rank([F[i] for i in (0,2,3)])==3,'four-profile factor and TT-prefix formal rank')
        factors[name]=F;cancellation[name]=e
    checks=[]
    targets=json.loads((HERE/'phase.json').read_text())
    need(words==[[int(c) for c in w] for w in targets['query_words']],'same universal query order')
    for obj in targets['records']:
        h,s,q,F,e,P2,S2,E3,Q2=obj['component_nine'];k=2
        z=[h,0,q,(1<<k)*F]
        predict=[sum(z[i]*factors['root_mark'][i][j] for i in range(4)) for j in range(17)]
        pair=next(x for x in obj['joint_pairs'] if x['ports']==[39,31])
        need(predict==pair['observer17'],'actual root/mark static17 factor')
        root=[e,h,h,h+(1<<k)*F,(1<<k)*P2-S2]
        mark=[e,q,q,(1<<(k+1))*F-q,S2+(1<<k)*Q2]
        need(pair['port_data']==[root,mark] and pair['ordered_W2']==[0,0,(1<<k)*F,0],
             'actual natural-pair endpoint and W2 formula')
        need(pair['P2']==2*(3**k*E3+(1<<k)*P2-S2),'actual natural pair P2')
        checks.append(dict(H=obj['H'],three_coordinates=[h,q,(1<<k)*F],observer17=predict))
    scalar_controls=0
    for code in (0,13,42,63):
        T=rows_code(4,code);v=ten(T,(0,1),0);h,e,S2,F,Q2,q,P2,S3,E3,Q3=v
        expectedO=[h+S2-e,q,2*S2-q,F+Q2,2*Q2,q,2*(P2-F+S3),3*S3-Q2,2*E3+Q3,3*Q3]
        expectedD=[h,0,2*S2,F,2*Q2,q,2*P2,3*S3,2*E3,3*Q3]
        for incoming,expected in ((14,expectedO),(15,expectedD)):
            need(ten(birth(T,incoming),(0,1,4),0)==expected,'ten-state O/D update')
            scalar_controls+=1
    sources=[src,HERE/'phase.json',HERE/'natural_pair.py',HERE/'NATURAL_PAIR.md',HERE/'gap.py',HERE/'check.py',
             RESULTS/'two_port_return_induction_20261103/ROOTED_RESETS.md',
             RESULTS/'two_port_return_induction_20261103/rooted_resets.json']
    print(json.dumps(dict(scope='Exact coefficient identity on five universal one-port signatures; eight tiny ten-state controls and two actual natural pairs',
        query_words=words,factors_in_profile_order_h_s_q_F=factors,annihilated_e_coefficients=cancellation,
        profile_factor_rank=4,TT_prefix_formal_factor_rank=3,actual_target_controls=checks,
        ten_state_controls=scalar_controls,source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}),indent=2,sort_keys=True))
if __name__=='__main__':main()
