"""Independent literal controls for the canonical prefix-tail classification."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json
import audit_defect as audit

HERE=Path(__file__).resolve().parent


def components(rows,mask):
    reach={v:(rows[v]&mask)|(1<<v) for v in range(len(rows)) if mask>>v&1}
    for u in reach:
        for v in reach:
            if reach[v]>>u&1:reach[v]|=reach[u]
    unseen=mask;parts=[]
    while unseen:
        v=(unseen&-unseen).bit_length()-1
        comp=sum(1<<w for w in reach if reach[v]>>w&1 and reach[w]>>v&1)
        parts.append(comp);unseen&=~comp
    parts.sort(key=lambda S:-reach[(S&-S).bit_length()-1].bit_count())
    return parts


def submasks(mask):
    x=mask
    while True:
        yield x
        if not x:break
        x=(x-1)&mask


def check_graph(rows,roots=None):
    n=len(rows);full=(1<<n)-1
    if len(components(rows,full))!=1:return []
    starts=audit.starts(rows)
    ends=audit.starts([sum(((rows[w]>>v)&1)<<w for w in range(n)) for v in range(n)])
    H=[1]+[sum(x) for x in starts[1:]]
    reports=[]
    for r in range(n) if roots is None else roots:
        I=sum(1<<v for v in range(n) if rows[v]>>r&1);d=I.bit_count()
        if d<2:continue
        R=full^I^(1<<r);parts=components(rows,R)
        P=0
        for S in parts:
            if all(rows[v]&I==I for v in range(n) if S>>v&1):P|=S
            else:break
        F=sum(ends[S][r]*H[full^S] for S in range(1,full+1) if S>>r&1)
        q=H[full^(1<<r)];s=starts[full][r];b=F-H[full]
        audit.check((q==s)==bool(P),'canonical prefix iff neutrality')
        W=full^P^(1<<r);V=W|(1<<r)
        F0=sum(ends[S][r]*H[V^S] for S in submasks(V) if S>>r&1)
        q0=H[W];H0=H[V];s0=starts[V][r]
        C=F0-H0-q0+s0;E=F0-q0
        nonordinary=lateordinary=0
        for S in parts:
            if S&P or S.bit_count()!=1:continue
            v=(S&-S).bit_length()-1
            if rows[v]&I==I:lateordinary+=1
            else:nonordinary+=1
        audit.check(nonordinary<=b//2+(2 if d==2 else 0),'exceptional singleton bound')
        audit.check(lateordinary<=b//2,'late ordinary bound')
        audit.check(nonordinary+lateordinary<=b+(2 if d==2 else 0),'total tail bound')
        if P:
            Q=sum(H[S]*H[P^S] for S in submasks(P))
            audit.check(b==H[P]*C and C%2==0,'constant neutral defect')
            audit.check(H[full]==E*Q+(q0-C)*H[P],'exact prefix Hamilton response')
            audit.check(E>0 and components(rows,W)[-1]&I,'positive growth and strong tail')
            audit.check(b%2==0,'odd defect has no ordinary prefix')
        reports.append({'d':d,'b':b,'prefix_vertices':P.bit_count(),'prefix_H':H[P],
                        'H':H[full],'q0':q0,'C':C,'E':E})
    return reports


def labelled(n):
    pairs=list(combinations(range(n),2))
    for code in range(1<<len(pairs)):
        rows=[0]*n
        for k,(a,b) in enumerate(pairs):
            if code>>k&1:rows[a]|=1<<b
            else:rows[b]|=1<<a
        yield rows


def ordinal(parts):
    size=sum(map(len,parts));rows=[];offset=0
    for part in parts:
        later=((1<<size)-1)^((1<<(offset+len(part)))-1)
        rows.extend((r<<offset)|later for r in part);offset+=len(part)
    return rows


def prefix_insert(P,V):
    p=len(P);v=len(V);root=p+v-1
    rows=[row|(((1<<(v-1))-1)<<p) for row in P]+[row<<p for row in V]
    rows[root]|=(1<<p)-1
    return rows


def main():
    strong=roots=neutral=0;hist=Counter()
    for n in range(1,6):
        for rows in labelled(n):
            if len(components(rows,(1<<n)-1))!=1:continue
            strong+=1
            report=check_graph(rows);roots+=len(report)
            for r in report:hist[r['b']]+=1;neutral+=bool(r['prefix_vertices'])
    extra=[]
    for core in ((2,4,1),(2,4,9,3)):
        for left,right in ((0,0),(1,0),(0,1),(1,1)):
            P=ordinal([tuple(sum(1<<j for j in range(i+1,left)) for i in range(left)),core,
                       tuple(sum(1<<j for j in range(i+1,right)) for i in range(right))])
            for tail in ((),(1,),(2,),(1,2),(2,1)):
                rows=prefix_insert(P,audit.graph(tail))
                report=check_graph(rows,[len(rows)-1])
                audit.check(len(report)==1 and report[0]['prefix_H']==(3 if len(core)==3 else 5),'nontransitive canonical prefix')
                extra.append({'core_H':report[0]['prefix_H'],'left':left,'right':right,'tail':tail,**report[0]})
    # Explicit convention hostile: kernel W+r is not the unrooted W.
    E,B,_=__import__('tail_squares').coefficients((1,))
    audit.check((E,B)==(6,1),'q0 intercept rather than rooted H0')
    report={'status':'PASS independent literal controls for the proved canonical decomposition',
            'all_labelled_strong_graphs_orders1_through5':strong,'all_branching_roots':roots,
            'neutral_roots':neutral,'defect_histogram':dict(hist),'nontransitive_prefix_controls':extra,
            'scope':'All strong labelled tournaments through order5 and every root with indegree>=2;40 additional C3/S4-prefix graphs with specified root. Universal classification follows from the proof, not finite extrapolation.'}
    out=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (HERE/'dyadic_defect.json').write_text(out,encoding='utf-8')
    print(json.dumps({'status':report['status'],'strong_graphs':strong,'branching_roots':roots,'nontransitive_controls':len(extra)}))


if __name__=='__main__':main()
