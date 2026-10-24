"""Same-agent audit: original matrices, whole census, sparse moves and literal rows."""
from common import *
from itertools import permutations,combinations
from collections import defaultdict,Counter
from math import comb
CHECKS=0
def check(v,m):
    global CHECKS
    CHECKS+=1;need(v,m)
def main():
    inv=inventory();ids=[r['source'] for r in inv];nodes,owner,merged,X,Y=data();labelmap={s:label(X,Y,s) for s in ids}
    costs=[];reservoir=set();raw_count=0
    for r in inv:
        s=r['source'];B=matrix(nodes[s]['canon_a']);p=r['before']['path'];q=r['after']['path'];fixed=[i for i in range(8) if p[i]==q[i]]
        start,end,j,u,v=r['pins'];check(p[0]==q[0]==start and p[-1]==q[-1]==end and p[j:j+2]==q[j:j+2]==[u,v],'literal shared owner roles')
        check(labelmap[s]==(sum(b<<j for j,b in enumerate(r['source_label'][:4])),r['source_label'][4]),'original paired potentials agree with inventory')
        D=[];auts=0
        for h in permutations(range(8)):
            if all(B[a][b]==B[h[a]][h[b]] for a,b in PAIRS):auts+=1
            if all(h[i]==p[i] for i in fixed) and all(B[a][b] for a,b in zip(h,h[1:])):D.append(h)
        check(auts==1,'no ignored source automorphism lifts')
        targetpaths=defaultdict(list)
        for h in D:
            w=word(B,h);t=merged[owner[FULL^w]];O=set(orbit(w))
            check(merged[owner[w]]==s and len(O)==4 and not O&reservoir,'source, closed ports, injective actual resource')
            reservoir.update(O);targetpaths[t].append(h)
        raw_count+=len(D)
        check(targetpaths[1965]==[tuple(p)] and targetpaths[2234]==[tuple(q)],'complete margin-eligible pinned paths')
        check(set(targetpaths)==({1559,1965,2234} if s==1923 else {1965,2234}),'exact raw domain, no new path hidden')
        costs.append(sum(a!=b for a,b in zip(p,q)))
    check(raw_count==57 and len(reservoir)==228,'complete inventory normalized capacity')
    check(label(X,Y,1965)==label(X,Y,2234)==label(X,Y,570),'whole target charge equality')
    report=[]
    for budget in (6,8):
        A=[json.loads(s) for s in (HERE/f'fibres{budget}.jsonl').read_text().splitlines()]
        B=[json.loads(s) for s in (HERE/f'audit_fibres{budget}.jsonl').read_text().splitlines()]
        check(len(A)==len(B)==(15 if budget==6 else 2),'complete grades of both enumeration implementations')
        for a,b in zip(A,B):
            check(a['rows']==b['rows'] and a['states']==b['states']==comb(28,a['k']),'independent lexicographic/intersection census agrees in every fibre')
            check(sum(r[1] for r in a['rows'])==a['states'],'every enumerated subset belongs to a fibre')
            check(a['forest_edges']==a['states']-sum(r[2] for r in a['rows']),'forest covers every component')
            bad=sum(r[2]>1 for r in a['rows']);check(bad==({4:197,5:11}.get(a['k'],0) if budget==6 else 0),'exact exceptional grades')
            report.append({k:a[k] for k in ('k','budget','states','fibres','disconnected','max_root_depth')})
    check(2*sum(comb(28,k) for k in range(14))+comb(28,14)==1<<28,'complement covers all binary assignments')
    def charge(m):
        y=x=0
        for i,s in enumerate(ids):
            if m>>i&1:y^=labelmap[s][0];x+=labelmap[s][1]
        return y,x%128
    def work(m):return sum(c for i,c in enumerate(costs) if m>>i&1)
    small=set()
    for k in (2,3,4):
        bank=defaultdict(list)
        for I in combinations(range(28),k):
            m=sum(1<<i for i in I);bank[charge(m)].append(m)
        stepmax=0
        for S in bank.values():
            for a,b in combinations(S,2):
                A=a&~b;B=b&~a
                if A>B:A,B=B,A
                if k<4:small.add((A,B))
                stepmax=max(stepmax,work(A|B))
        if k==3:max6=stepmax
        if k==4:max8=stepmax
    cert=json.loads((HERE/'bridges.json').read_text());check(cert['sources']==ids,'sparse bridge source identities')
    check(len(small)==1944 and len(cert['moves'])==220,'full small alphabet and sharp extra alphabet')
    check(max(work(a|b) for a,b in cert['moves'])==22,'sharp-clock sparse completion')
    # Authentic harmonic generators, accumulated on physical owners, not role copies.
    chars=[X]+[[64*y for y in row] for row in Y[:4]]
    for A,B in list(small)+[tuple(x) for x in cert['moves']]:
        check(not A&B and A.bit_count()==B.bit_count() and charge(A)==charge(B),'trade margins and full charge')
        for chi in chars:
            row=defaultdict(int)
            for i,s in enumerate(ids):
                delta=int(bool(A>>i&1))-int(bool(B>>i&1))
                if not delta:continue
                # -2 delta on sU and +2 delta on sV, including diagonals.
                dU=-2*delta*(chi[s]-chi[1965]);dV=2*delta*(chi[s]-chi[2234])
                row[s]+=dU+dV;row[1965]-=dU;row[2234]-=dV
            check(all(x%256==0 for x in row.values()),'every physical row annihilates actual dyadic generator')
    # Rebuild small boundary graphs without producer's component helper or forest selection.
    max_sparse_depth=0
    for k in (4,5):
        bank=defaultdict(list)
        for I in combinations(range(28),k):
            m=sum(1<<i for i in I);bank[charge(m)].append(m)
        if k==4:
            check(sorted(bank[6,123])==[12864,69272576] and work(12864^69272576)==22,'independent whole-inventory two-state sharp-clock hostile')
        for S in bank.values():
            unused=set(S);root=min(unused);unused.remove(root);todo=[(root,0)]
            for a,depth in todo:
                hits={b for b in unused if (a^b).bit_count()<=6}
                for U,V in cert['moves']:
                    if a&U==U and not a&V:hits.add((a^U)|V)
                    if a&V==V and not a&U:hits.add((a^V)|U)
                hits&=unused;unused.difference_update(hits)
                if hits:max_sparse_depth=max(max_sparse_depth,depth+1)
                todo.extend((b,depth+1) for b in sorted(hits))
            check(not unused,'independent sparse-alphabet reachability in every boundary fibre')
    # Recover the old score-complement law and its all-context target obstruction.
    S={s:sorted(map(sum,matrix(nodes[s]['canon_a']))) for s in (570,1965,2234)}
    check(S[570] not in [S[1965],sorted(7-x for x in S[1965]),S[2234],sorted(7-x for x in S[2234])],'target570 outside fixed-endpoint score sector')
    check(sum(x*x for x in S[570])==112 and all(sum(x*x for x in S[s])==106 for s in (1965,2234)),'even converse-invariant scalar separates targets')
    return {'status':'PASS, same-agent independent implementation and original-native audit','checks':CHECKS,'census':report,
            'complete_binary_assignments':1<<28,'raw_paths':raw_count,'raw_closed_ports':len(reservoir),
            'small_moves':len(small),'sparse_eight_slot_moves':len(cert['moves']),'max_sparse_boundary_root_depth':max_sparse_depth,
            'max_owner_clocks_six_slot_move':max6,'max_owner_clocks_eight_slot_move':max8,
            'max_owner_clocks_sparse_alphabet':22,
            'target_scores':S,'input_revision':REV,'inputs':INPUTS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
