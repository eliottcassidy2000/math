"""Exact local-partner gluing versus direct global-host enumeration."""
from itertools import combinations,permutations
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent;COUNT=0
PAIRS=list(combinations(range(4),2))
def need(b,why):
    global COUNT
    COUNT+=1
    if not b:raise RuntimeError(why)
def matrix(word):
    A=[[0]*4 for _ in range(4)]
    for k,(u,v)in enumerate(PAIRS):A[u][v]=(word>>k)&1;A[v][u]=1-A[u][v]
    return A
def records(A,cell):
    out=[]
    for p in permutations(cell):
        if not all(A[u][v]for u,v in zip(p,p[1:])):continue
        flags={tuple(sorted((u,v)))for u,v in zip(p,p[1:])}
        arcs={(u,v):A[u][v]if(u,v)in flags else 1-A[u][v]for u,v in combinations(sorted(cell),2)}
        out.append((p,flags,arcs))
    return out
def direct_hosts(local):
    out=[]
    for w in range(64):
        A=matrix(w)
        if all(all(A[u][v]==a for(u,v),a in r[2].items())for r in local):out.append(w)
    return out
def main():
    cases=positive=negative=0;cells=list(combinations(range(4),3));sample=None
    for w in range(64):
        A=matrix(w);op=[[A[v][u]for v in range(4)]for u in range(4)]
        for c,d in combinations(cells,2):
            overlap=set(combinations(sorted(set(c)&set(d)),2))
            covered=set(combinations(c,2))|set(combinations(d,2))
            for r in records(A,c):
                for s in records(A,d):
                    flag_ok=all((e in r[1])==(e in s[1])for e in overlap)
                    hosts=direct_hosts([r,s]);need(bool(hosts)==flag_ok,'flag criterion iff actual host existence')
                    need(len(hosts)==(2**(6-len(covered))if flag_ok else 0),'exact free-crossarc extension count')
                    rr=next(z for z in records(op,c)if z[0]==r[0][::-1]);ss=next(z for z in records(op,d)if z[0]==s[0][::-1])
                    need(rr[1]==r[1] and ss[1]==s[1],'literal converse preserves every shared-pair flag')
                    need(set(direct_hosts([rr,ss]))=={63^h for h in hosts},'literal converse of all host extensions')
                    cases+=1;positive+=flag_ok;negative+=not flag_ok
                    if sample is None and not flag_ok:sample=dict(source_word=w,cells=[c,d],paths=[r[0],s[0]],shared_pairs=sorted(overlap))
    # Equal-size transitive hostile: owner order a,d,b,c is indexed0,1,2,3.
    A=matrix(63);r=records(A,(0,2,3))[0];s=records(A,(0,1,2))[0]
    need(r[0]==(0,2,3) and s[0]==(0,1,2) and not direct_hosts([r,s]),'actual fourowner transitive hostile')
    short=records(A,(0,2))[0]
    need(not direct_hosts([short,s]),'minimal unequal-size threeowner hostile')
    # Compatible local pairs need not extend to one global observer.
    star=[records(A,(0,j))[0]for j in (1,2,3)]
    need(len(direct_hosts(star))==8,'three free crossarcs in glued star example')
    allflags=set().union(*(r[1]for r in star))
    observer_extensions=[p for p in permutations(range(4))if allflags<={tuple(sorted((u,v)))for u,v in zip(p,p[1:])}]
    need(not observer_extensions,'three forced neighbours cannot form one global path')
    out=dict(status='FINITE-EXACT controls for elementary PROVED partner-gluing lemma',source_tournaments=64,
        distinct_threecell_pairs=6,legal_observer_cases=cases,compatible_cases=positive,incompatible_cases=negative,
        first_incompatible=sample,equal_size_hostile=dict(source_order=['a','d','b','c'],paths=[r[0],s[0]]),
        compatible_no_global_observer=dict(source_word=63,paths=[z[0]for z in star],partner_hosts=direct_hosts(star)),checks=COUNT)
    out=json.loads(json.dumps(out));path=HERE/'partner_gluing_receipt.json'
    if '--write'in sys.argv:path.write_text(json.dumps(out,indent=2))
    else:need(json.loads(path.read_text())==out,'frozen gluing receipt')
    print('PASS',COUNT,'checks;',cases,'legal observer pairs;',positive,'glue;',negative,'fail')
if __name__=='__main__':main()
