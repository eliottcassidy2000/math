"""Original-path audit plus adversarial all-grade route tests. No asserts."""
from route import *
from itertools import permutations
import random

def main():
    receipt=json.loads((HERE/'menus.json').read_text());rows=receipt['rows']
    nodes,owner,merged,X,Y=data();autsizes={};liftcounts={};checks=0
    def check(v,m):
        nonlocal checks
        checks+=1;need(v,m)
    # Independently permute all 8! orders; producer uses adjacency-pruned DFS.
    permutations8=list(permutations(range(8)));reserved={};literal=0
    for r in rows:
        s=r['source'];B=matrix(nodes[s]['canon_a']);pins=r['pins'];eligible=[];raw=[];auts=[];lifts=Counter()
        for p in permutations8:
            if all(B[i][j]==B[p[i]][p[j]] for i,j in PAIRS):auts.append(p)
            if not all(B[a][b] for a,b in zip(p,p[1:])):continue
            w=word(B,p);lifts[w]+=1
            if not all(p[i]==v for i,v in pins):continue
            t=merged[owner[FULL^w]];raw.append(t)
            if t in receipt['targets']:eligible.append((p,t,w))
        expected=sorted((tuple(o['path']),o['target'],o['word']) for o in r['options'])
        check(sorted(eligible)==expected,'complete pin domain: '+str(s))
        check(dict(Counter(raw))=={int(k):v for k,v in r['raw_targets'].items()},'raw domain before margin filter')
        autsizes[s]=len(auts)
        for o in r['options']:
            p=o['path'];w=o['word'];t=o['target'];O=tuple(o['closed_words']);literal+=1
            liftcounts[str(s)+':'+str(w)]=lifts[w]
            check(lifts[w]==len(auts),'all canonical-source kernel lifts retained as one resource')
            K=[[B[i][j] for j in range(8)] for i in range(8)]
            arcs={frozenset((a,b)) for a,b in zip(p,p[1:])}
            for i,j in PAIRS:
                if frozenset((i,j)) not in arcs:K[i][j]^=1;K[j][i]^=1
            check(word(K,p)==FULL^w,'actual partner word')
            C=[list(x) for x in zip(*B)];rev=p[::-1]
            check(word(C,rev)==mirror(w),'actual converse and reversed integer clock')
            check(merged[owner[w]]==s and merged[owner[word(K,p)]]==t,'source-target lookup')
            check(O==tuple(orbit(w)) and len(O)==4,'full inverse/converse resource')
            check(o['resource'] not in reserved,'no normalized resource reuse');reserved[o['resource']]=s
            check((X[s]-X[t])%8==(X[s]-146)%8,'selected ordinary mod8 X phase')
    check(literal==171 and len(reserved)==171,'85 occupied slots,171 distinct options')
    chars=[X]+[[64*y for y in row] for row in Y[:4]]
    W=[r['source'] for r in rows if any(o['target']==570 for o in r['options'])]
    physical=defaultdict(list)
    for m in range(1<<len(W)):
        currents=tuple(sum(2*(chi[570]-chi[s]) for i,s in enumerate(W) if m>>i&1)%256 for chi in chars)
        physical[(m.bit_count(),currents)].append(m)
    claimed=json.loads((HERE/'bank.json').read_text())['W_fibres']
    check(sorted(physical.values())==sorted(f['masks'] for f in claimed),'entire W partition from authentic physical target currents')
    rng=random.Random(57020261024);totals=[];example=None;event_seen=set()
    for remove in (False,True):
        R=Router(remove);tested=events=0;max_support=max_work=max_length=0
        lowmax=[max(c for c,m in d.values()) for d in R.low]
        for k in range(R.n+1):
            for r in range(max(0,k-len(R.E)),min(20,k)+1):
                start=sum(1<<i for i in rng.sample(R.R,r)+rng.sample(R.E,k-r));q=R.q(start)
                if k<=4:finish=R.low[k][q][1]
                elif R.n-k<=4:finish=R.full^R.low[R.n-k][R.q(R.full^start)][1]
                else:
                    t=max(5,k-len(R.E));outside=sum(1<<i for i in rng.sample(R.E,k-t))
                    finish=R.replace(outside,t,q)
                P=R.path(start,finish);check(P[0]==start and P[-1]==finish,'prescribed full endpoints')
                tested+=1;max_length=max(max_length,len(P)-1)
                for m in P:
                    # Independent unpacked source sums.
                    y=x=0
                    for i,row in enumerate(R.rows):
                        if m>>i&1:y^=row['label'][0];x+=row['label'][1]
                    check(m.bit_count()==k and (y,x%128)==(q&15,q>>4),'tuple arrival charge and integer count')
                for a,b in zip(P,P[1:]):
                    events+=1;support=(a^b).bit_count();work=0
                    max_support=max(max_support,support)
                    for i,row in enumerate(R.rows):
                        if not (a^b)>>i&1:continue
                        opt={o['target']:o for o in row['options']}
                        p=opt[2234 if a>>i&1 else 1965]['path'];h=opt[2234 if b>>i&1 else 1965]['path']
                        work+=sum(u!=v for u,v in zip(p,h))
                        check(all(p[j]==h[j]==v for j,v in row['pins']),'shared fixed clock survives')
                    max_work=max(max_work,work);check(support<=20 and work<=45,'literal event bounds')
                    if (remove,a,b) in event_seen:continue
                    event_seen.add((remove,a,b))
                    for chi in chars:
                        debt=defaultdict(int)
                        for i,s in enumerate(R.ids):
                            if not (a^b)>>i&1:continue
                            u=2234 if a>>i&1 else 1965;v=2234 if b>>i&1 else 1965
                            old=2*(chi[s]-chi[u]);new=2*(chi[s]-chi[v])
                            debt[s]+=new-old;debt[u]+=old;debt[v]-=new
                        check(all(v%256==0 for v in debt.values()),'all physical owners, authentic generators')
                if example is None and k==38 and r==0:example=R.certificate(start,finish)
        totals.append({'binary_sources':R.n,'endpoint_pairs':tested,'events':events,'max_observed_support':max_support,
                       'max_observed_clocks':max_work,'max_observed_events':max_length,'max_low_grade_witness_cost':lowmax})
    # Rejected checksum shortcut: parity carries must not reach the cyclic factor.
    check(add(15,15)==0,'packed-law hostile');check(add(16*127+15,16+1)==14,'independent wrap and xor')
    return {'status':'PASS same-agent independent native enumeration and exact route checks',
            'checks':checks,'aut_size_histogram':dict(Counter(autsizes.values())),
            'kernel_lift_histogram':dict(Counter(liftcounts.values())), 'route_tests':totals,
            'example':example,'inputs':INPUTS,
            'scope':'routes tested across every binary grade and reservoir occupancy; universal coverage follows the graded-section proof, not sampled endpoints'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
