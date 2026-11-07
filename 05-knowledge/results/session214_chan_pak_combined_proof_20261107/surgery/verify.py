"""Exact certificate: a two-apex class with <=1 chain exit omits 7477.
No third-party packages. Run normally and under -O. Output is deterministic.
"""
from itertools import permutations, product
import hashlib
import json
from counting import count, coeff, scan, adj, dp

def need(test,message):
    if not test: raise RuntimeError(message)

def direct_permutations(out):
    return sum(all(out[p[i]] >> p[i+1] & 1 for i in range(len(out)-1))
               for p in permutations(range(len(out))))

def run():
    need(all(7477%d for d in range(2,87)),'7477 primality')
    need(3*5**5==9375,'Busch order18 evaluation')
    dpchecks=permchecks=polynomialchecks=0
    for L in range(1,7):
        for r in range(-1,L):
            for ab in (0,1):
                c,lin,quad=coeff(L,r,ab)
                for q in range(1<<L):
                    a=[int(i==r) for i in range(L)];b=[q>>i&1 for i in range(L)]
                    out=adj(a,b,ab);h=count(a,b,ab)
                    need(h==dp(out),('subset DP',L,r,ab,q,h));dpchecks+=1
                    pred=c+sum(lin[i]*b[i] for i in range(L))+sum(v*b[i]*b[j] for (i,j),v in quad.items())
                    need(h==pred,('quadratic',L,r,ab,q));polynomialchecks+=1
                    if L<=4:
                        need(h==direct_permutations(out),('permutations',L,r,ab,q));permchecks+=1
    witnesses={1:(1,0,1,0)};cases=0;levels=[]
    for L in range(1,18):
        for r in range(L):
            for ab in (0,1):
                for h,q in scan(L,r,ab,7477):
                    need(h>0 and h%2==1,('parity',L,r,ab,q,h))
                    witnesses.setdefault(h,(L,r,ab,q))
                cases+=1<<L
        missing=[h for h in range(1,7478,2) if h not in witnesses]
        levels.append({'chain_length':L,'covered_through_7477':len(witnesses),'first_missing_other_than_7_21':next((h for h in missing if h not in (7,21)),None)})
    missing=[h for h in range(1,7478,2) if h not in witnesses]
    need(missing==[7,21,7477],('complete omitted list',missing))
    # Replay every positive witness with transfer, not its interpolated polynomial.
    for h,(L,r,ab,q) in witnesses.items():
        a=[int(i==r) for i in range(L)];b=[q>>i&1 for i in range(L)]
        need(count(a,b,ab)==h,('positive witness replay',h))
    # A two-exit control crosses the obstruction.
    a=[int(i in (0,3)) for i in range(9)];b=[199>>i&1 for i in range(9)]
    out=adj(a,b,1)
    need(out==[510,1020,1016,1520,2016,1984,896,768,1536,1033,199],'control rows')
    need(count(a,b,1)==7477 and dp(out)==7477,'two-exit positive control')
    raw=json.dumps([[h,*witnesses[h]] for h in sorted(witnesses)],separators=(',',':')).encode()
    return {'status':'PASS','scope':'FINITE-EXACT census; all-order exclusion additionally uses the written SCC proof and CITED Busch bound',
       'prime_target':7477,'busch_f18':9375,'universes':{'chain_lengths':[1,17],'apex_chain_outdegree':1,'apex_arc_orientations':2,'other_apex_cuts':'all','parameter_cases':cases},
       'omitted_odd_values_through_7477':missing,'positive_values':len(witnesses),
       'positive_witness_sha256':hashlib.sha256(raw).hexdigest(),
       'positive_witness_replays':len(witnesses),'subset_dp_controls':dpchecks,'permutation_controls':permchecks,'quadratic_controls':polynomialchecks,
       'two_exit_control':{'chain_length':9,'a_out_chain':[0,3],'a_beats_b':True,'b_cut':199,'outgoing_rows':out,'H':7477},
       'selected_witnesses':{str(h):witnesses[h] for h in (13,35,47,7475)},'levels':levels}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
