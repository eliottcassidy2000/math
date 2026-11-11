"""Exact next-apex supports of the count-neutral two-arc switch pairs."""
import importlib.util
import json
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'session251_creative_completion_20261111/ternary/compiler.py'
spec=importlib.util.spec_from_file_location('inherited_compiler',SOURCE)
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
core=(2,4,9,19,7)

def require(p,msg):
    if not p:raise RuntimeError(msg)

def capacities(rows):
    n=len(rows);full=(1<<n)-1;hs=c.hp_subsets(rows);C=[[0]*n for _ in rows];F=[0]*n
    for start in range(n):
        dp=[[0]*n for _ in range(full+1)];dp[1<<start][start]=1
        for mask in range(1,full+1):
            if not(mask>>start&1):continue
            for v,val in enumerate(dp[mask]):
                if not val:continue
                F[v]+=val*hs[full^mask]
                if not(mask.bit_count()&1):C[start][v]+=val*hs[full^mask]
                av=rows[v]&~mask&full
                while av:
                    b=av&-av;av-=b;dp[mask|b][b.bit_length()-1]+=val
    return hs[full],C,[f-hs[full] for f in F]

def append_apex(rows,O):
    n=len(rows)
    return [r|(0 if O>>v&1 else 1<<n) for v,r in enumerate(rows)]+[O]

def output(h,C,O):
    n=len(C)
    return h+2*sum(C[i][j] for i in range(n) if O>>i&1 for j in range(n) if not(O>>j&1))

def main():
    report=[];direct=0
    for t in range(8):
        graphs=[c.construct(core,O,t+1,t) for O in (1,2)]
        responses=[];hs=[];root_gaps=[]
        for rows in graphs:
            h,C,gaps=capacities(rows);hs.append(h);root_gaps.append(gaps)
            by_value={}
            for O in range(1<<len(rows)):
                v=output(h,C,O);by_value.setdefault(v,[]).append(O)
                if t==0:
                    require(c.hp_subsets(append_apex(rows,O))[-1]==v,'literal apex control')
                    direct+=1
            responses.append(by_value)
        require(hs[0]==hs[1],'neutral switch')
        h=hs[0];p,q=responses
        left=sorted(set(p)-set(q));right=sorted(set(q)-set(p))
        require(left and right,'changed response support')
        row={'t':t,'A':t+1,'B':t,'H':h,'source_rows':graphs,
             'support_sizes':[len(p),len(q)],
             'support_difference_sizes':[len(left),len(right)],
             'first_source_only':{'target':left[0],'cut':p[left[0]][0]},
             'second_source_only':{'target':right[0],'cut':q[right[0]][0]},
             'root_defects_F_minus_H':root_gaps,
             'ready_roots':[[v for v,gap in enumerate(gaps) if gap in (1,2)] for gaps in root_gaps],
             'branch_cut_counts':[[len(z.get(3*h+2,[])),len(z.get(3*h+4,[]))] for z in responses]}
        require(row['branch_cut_counts']==([[1,0],[1,0]] if t==0 else [[0,0],[0,0]]),'one-apex branch result')
        require(row['ready_roots']==([[0],[0]] if t==0 else [[],[]]),'rooted two-birth readiness')
        if t==0:require(all(gaps[0]==1 for gaps in root_gaps),'H77 root is1-ready')
        report.append(row)
    return {'status':'PASS','direct_new_apex_controls':direct,'pairs':report,
            'scope':'All next-apex cuts and every rooted two-bin defect for each of16 named sources; one-apex branch attainment and two-birth readiness are separate tests. No claim over all same-count sources or all t.'}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
