"""Smallest-label-count counterexample to naive charge-coverage padding."""
import json,sys
from itertools import product,combinations,combinations_with_replacement
from pathlib import Path
HERE=Path(__file__).resolve().parent
COUNT=0
def need(x,why):
    global COUNT
    COUNT+=1
    if not x:raise RuntimeError(why)
def coverage(labels,k,moduli):
    return {tuple(sum(labels[i][j]for i in ids)%m for j,m in enumerate(moduli))for ids in combinations(range(len(labels)),k)}
def main():
    groups=[(1,),(2,),(3,),(4,),(2,2)];cases=0;hostiles=[];padding_cases=0
    for moduli in groups:
        G=list(product(*(range(m)for m in moduli)));allG=set(G)
        for h in range(1,5):
            for labels in combinations_with_replacement(G,h):
                for k in range(1,(h+1)//2):
                    cases+=1;S=coverage(labels,k,moduli);T=coverage(labels,k+1,moduli)
                    if S==allG and T!=allG:hostiles.append(dict(moduli=moduli,labels=labels,k=k,missing=sorted(allG-T)))
                    for f in range(h):
                        deleted=labels[:f]+labels[f+1:];small=coverage(deleted,k,moduli)
                        forced={tuple((labels[f][j]+x[j])%m for j,m in enumerate(moduli))for x in small}
                        direct={tuple(sum(labels[i][j]for i in ids)%m for j,m in enumerate(moduli))for ids in combinations(range(h),k+1)if f in ids}
                        need(forced==direct,'forced-label translation identity')
                        need((forced==allG)==(small==allG),'deletion coverage iff forced padding coverage');padding_cases+=1
    need(len(hostiles)==1,'unique minimal-universe multiset hostile')
    z=hostiles[0];need(z['moduli']==(2,2) and z['k']==1 and len(z['labels'])==4 and z['missing']==[(0,0)],'four-label binarygroup hostile')
    matchings={c:[list(e)for e in combinations(range(4),2)if e[0]^e[1]==c]for c in (1,2,3)}
    for c,edges in matchings.items():need(sorted(v for edge in edges for v in edge)==list(range(4)),'K4 charge classes are perfect matchings')
    need(coverage(z['labels'],3,(2,2))==set(z['labels']),'complementary occupancy3 covers all charges')
    result=dict(status='PROVED tiny counterexample and forced-padding equivalence, FINITE-EXACT minimality control',
        group_invariant_factors=groups,label_counts=[1,2,3,4],multiset_candidate_cases=cases,forced_padding_cases=padding_cases,hostiles=hostiles,K4_charge_matchings=matchings,checks=COUNT)
    result=json.loads(json.dumps(result));out=HERE/'padding_probe_receipt.json'
    if '--write' in sys.argv:out.write_text(json.dumps(result,indent=2))
    else:need(json.loads(out.read_text())==result,'frozen receipt')
    print('PASS',COUNT,'checks;',cases,'multiset candidates;',padding_cases,'forced-padding controls; sole hostile F2^2 all four labels')
if __name__=='__main__':main()
