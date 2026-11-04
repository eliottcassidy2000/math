"""Independent direct endpoint calculation of the positive root-star potential."""
from itertools import combinations
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'zero_ready_strong_switch_20261103'/'literal'

def need(ok,message):
    if not ok:raise RuntimeError(message)

def dp(rows):
    n=len(rows);full=(1<<n)-1
    end=[[0]*n for _ in range(full+1)];h=[0]*(full+1);h[0]=1
    incoming=[sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)]
    for mask in range(1,full+1):
        todo=mask
        while todo:
            bit=todo&-todo;todo-=bit;v=bit.bit_length()-1;prev=mask^bit
            if not prev:end[mask][v]=1
            else:
                u_set=prev&incoming[v]
                while u_set:
                    ub=u_set&-u_set;u_set-=ub
                    end[mask][v]+=end[prev][ub.bit_length()-1]
        h[mask]=sum(end[mask])
    return h,end

def all_data(rows):
    n=len(rows);full=(1<<n)-1
    h,e=dp(rows);_,s=dp([full^row^(1<<i) for i,row in enumerate(rows)])
    f=[sum(e[mask][i]*h[full^mask] for mask in range(1,full+1)) for i in range(n)]
    j=[[sum(e[mask][i]*s[full^mask][k] for mask in range(1,full+1))
        for k in range(n)] for i in range(n)]
    return h[full],e[full],s[full],f,j

def all_tournaments(n):
    pairs=list(combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        rows=[0]*n
        for k,(a,b) in enumerate(pairs):
            u,v=(a,b) if mask>>k&1 else (b,a)
            rows[u]|=1<<v
        yield rows

def birth(rows,I):
    n=len(rows);return [row|(((I>>v)&1)<<n) for v,row in enumerate(rows)]+[((1<<n)-1)^I]

def main():
    cases=0
    for n in range(5):
        for rows in all_tournaments(n):
            q,e,s,f,j=all_data(rows);values=[]
            for I in range(1<<n):
                graph=birth(rows,I);gh,_,_,gf,_=all_data(graph)
                b=sum(s[i] for i in range(n) if I>>i&1)
                b+=sum(j[i][k]+j[k][i] for i,k in combinations(range(n),2)
                       if I>>i&1 and I>>k&1)
                need(gf[n]==q+sum(f[i] for i in range(n) if I>>i&1),'linear endpoint formula')
                need(gf[n]-gh==b,'quadratic defect formula')
                values.append(b);cases+=1
            for I in range(1<<n):
                available=[v for v in range(n) if not I>>v&1]
                for v in available:
                    need(values[I|1<<v]-values[I]==s[v]+sum(j[i][v]+j[v][i] for i in range(n) if I>>i&1),'first difference')
                for v,w in combinations(available,2):
                    need(values[I|1<<v|1<<w]-values[I|1<<v]-values[I|1<<w]+values[I]==j[v][w]+j[w][v],'second difference')
    seed=json.loads((OLD/'realized_generation2_seed4_depth3_all_prefixes.json').read_text())[1]
    names=seed['suffix_names'];r=names.index(32);vs=[i for i in range(14) if i!=r];dn=[names[i] for i in vs]
    D=[sum(((seed['suffix_rows'][u]>>v)&1)<<i for i,v in enumerate(vs)) for u in vs]
    q,e,s,f,j=all_data(D);a,b=dn.index(22),dn.index(33)
    need((q,s[a],s[b],j[a][b],j[b][a])==(13043,0,0,890,2070),'actual named deleted collar')
    need(s[a]+s[b]+j[a][b]+j[b][a]==2960,'actual pure interaction defect')
    bh,_,_,bf,bj=all_data(seed['suffix_rows'])
    need((bh,bf[r],bj[r][names.index(33)],bj[r][names.index(22)])==(80319,83279,890,2070),'separate full-suffix queries')
    receipt={'status':'FINITE-EXACT','all_labelled_deletion_orders':[0,4],
             'all_incoming_set_cases':cases,'actual_deletion_names':dn,'actual_deletion_rows':D,
             'actual_H_deletion':q,'start_counts':[s[a],s[b]],'directed_joint_counts':[j[a][b],j[b][a]],
             'full_suffix_H_F32':[bh,bf[r]],'full_suffix_J32_to33_to22':[bj[r][names.index(33)],bj[r][names.index(22)]],
             'defect':2960,'source':'previous packet literal generation2_seed4 index1'}
    (HERE/'defect_check.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))

if __name__=='__main__':main()
