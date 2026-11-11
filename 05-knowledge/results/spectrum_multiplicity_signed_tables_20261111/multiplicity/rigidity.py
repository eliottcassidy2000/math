"""Structural controls only: inherited templates, degree identities, path witnesses.

No tournament enumeration, Hamilton-path counting, or automorphism search.
The all-order arguments are in RIGIDITY.md.
"""
from collections import defaultdict
from pathlib import Path
import hashlib,json,subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
SOURCES=[
    '05-knowledge/results/berggren_target_return_20261111/reset/target_ready_return.json',
    '05-knowledge/results/berggren_target_return_20261111/reset/TARGET_RETURN.md',
    '05-knowledge/results/two_entry_funnels_20261103/READY.md',
    '05-knowledge/results/polarized_path_calculus_20261103/FUNNELS.md',
]


def need(condition,label):
    if not condition:raise RuntimeError(label)


def inflate(T,v,k):
    labels=[(i,j) for i in range(len(T)) for j in range(k if i==v else 1)]
    return tuple(sum(1<<b for b,(x,y) in enumerate(labels)
                     if (i==x and j<y) or (i!=x and (T[i]>>x)&1)) for i,j in labels)


def indegrees(T):return [sum(row>>i&1 for row in T) for i in range(len(T))]


def main():
    raw={p:subprocess.check_output(['git','show','HEAD:'+p],cwd=ROOT) for p in SOURCES}
    inherited=json.loads(raw[SOURCES[0]])
    S=tuple(inherited['controls'][0]['source'])
    P=tuple(inherited['parent_templates'][0])
    need(S==(222,76,16,116,226,199,132,14),'inherited S template')
    need(P==(34,244,249,243,225,192,129,1),'inherited parent template')
    controls=[]
    for k in (1,2,7,16):
        T=inflate(S,0,k);degrees=indegrees(T)
        expected=list(range(1,k+1))+[k+3,k+5,k+2,k+2,2,k+4,k+3]
        need(degrees==expected,'symbolic S degree formula')
        inn=[i for i,row in enumerate(T) if row&1]
        need(inn==[k+4],'unique in-neighbour of first prefix')
        fixed={0,k+4};parts=defaultdict(list)
        for i,d in enumerate(degrees):
            if i not in fixed:parts[d].append(i)
        # Each residual invariant degree class has at most two vertices;
        # its internal oriented edge fixes a pair pointwise.
        need(all(len(v)<=2 for v in parts.values()),'degree-pair rigidity certificate')
        paths=[]
        for tail in ((1,2,4,6,7,3),(1,2,4,7,3,6)):
            path=[k+4]+list(range(1,k))+[k+v-1 for v in tail]
            need(sorted(path)==list(range(1,len(T))),'spanning parent witness')
            need(all(T[a]>>b&1 for a,b in zip(path,path[1:])),'parent path arcs')
            paths.append(path)
        parent=inflate(P,5,k);n=len(parent)
        ready=tuple(row|((1<<n) if i==0 else 0) for i,row in enumerate(parent))+(((1<<n)-1)^1,)
        need(indegrees(parent)[0]==5,'R parent root indegree')
        need([i for i,d in enumerate(indegrees(ready)) if d==1]==[n],'unique R indegree-one root')
        controls.append(dict(k=k,source_order=len(T),ready_order=len(ready),
                             source_indegrees=degrees,two_parent_paths=paths,
                             rigid_degree_parts=list(parts.values())))
    low=json.loads((HERE/'all_order_H_le375.json').read_text())['221']
    need(low['ordinary_polynomial_ascending']==['-262','20','2'],'inherited H221 polynomial')
    need(low['labelled_over_factorial_polynomial_ascending']==['-262','20','2'],'H221 automorphism mass')
    report=dict(status='PASS structural controls; all-order proofs in RIGIDITY.md',
                no_tournament_enumeration=True,no_path_count_or_automorphism_search=True,
                source_sha256={p:hashlib.sha256(data).hexdigest() for p,data in raw.items()},
                controls=controls,H221_core_words=low['core_words'],
                inherited_ready_H221='Exactly one strong ready class, R1 of order9, from TARGET_RETURN section4',
                inherited_ready_H1997='No strong1-ready class at any order, from TARGET_RETURN section4')
    (HERE/'rigidity_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: S/R degree identities, explicit non-ready path witnesses, H221 weighted fibre; no new enumeration')


if __name__=='__main__':main()
