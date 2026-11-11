"""Verify finite all-depth closure certificate WITHOUT running optimizer.
Uses exact matrix transitions from morphic_optimizer; all checks explicitly raise.
"""
from collections import defaultdict
from pathlib import Path
import json
from morphic_optimizer import U,V,pack,step,output,dominates,literal

def check(test,message):
    if not test:raise ValueError(message)

def verify(path=None):
    root=Path(__file__).parent
    path=Path(path) if path else root/'morphic_mod16_certificate.json'
    cert=json.loads(path.read_text());q=cert['modulus']
    check(q==16,'unexpected modulus')
    records=cert['representatives'];bank=defaultdict(list)
    for item in records:
        v=tuple(item['vector']);word=item['reverse_address']
        check(len(v)==32 and all(type(x)==int and x>=0 for x in v),'invalid vector')
        check(all(a in 'ABC' for a in word),'invalid address')
        reconstructed=pack(U,V)
        for a in word:reconstructed=step(reconstructed,a)
        check(reconstructed==v,'word does not reconstruct saved vector')
        bank[tuple(x%q for x in v)].append(v)
    seed=pack(U,V)
    def covered(v):return any(dominates(u,v) for u in bank.get(tuple(x%q for x in v),[]))
    check(covered(seed),'seed not dominated')
    transitions=0
    for item in records:
        v=tuple(item['vector'])
        for a in 'ABC':
            check(covered(step(v,a)),'transition not covered in its residue label');transitions+=1
    # Check genuine antichains too, though not needed for the closure proof.
    for vectors in bank.values():
        for i,v in enumerate(vectors):
            check(not any(dominates(w,v) for j,w in enumerate(vectors) if j!=i),'redundant antichain')
    best={}
    for item in records:
        h=output(tuple(item['vector']));r=h%q;w=item['reverse_address'];candidate=(h,len(w),w)
        if r not in best or candidate<best[r]:best[r]=candidate
    expected={1:17,3:3,5:5,7:408071,9:9,11:107,13:3181,15:1615}
    check({r:x[0] for r,x in best.items()}==expected,'wrong minima')
    witnesses=[]
    for r,(h,k,w) in sorted(best.items()):
        rec=literal(w[::-1]);check(rec['H']==h,'literal witness mismatch')
        rec['residue']=r;witnesses.append(rec)
    stats=cert['statistics']
    check(len(records)==stats['retained']==633,'retained count mismatch')
    check(len(bank)==stats['labels']==377,'label count mismatch')
    check(cert['depth_or_size_cap'] is None,'unexpected cap')
    result={'scope':'FINITE-EXACT COMPLETE SEMIGROUP CERTIFICATE',
            'modulus':q,'reachable_labels':len(bank),'retained':len(records),
            'construction_additions':stats['additions'],
            'checked_transition_closures':transitions,'checked_word_receipts':len(records),
            'depth_or_size_cap':None,'minimum_by_residue':witnesses}
    (root/'morphic_certificate_verified.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='minimum_by_residue'},indent=2))
    return result
if __name__=='__main__':verify()
