#!/usr/bin/env python3
"""Exact all-order identity certificate plus finite actual-word/hostile checks."""
import itertools,json,random
from word_recovery import *


def require(ok,message):
    if not ok: raise RuntimeError(message)


def occupancy_rows():
    rows=[[0]*D for _ in range(8)]
    for j,(o,p,q,r) in enumerate(STATES):
        po,qo,ro=int(p!=0),int(q!=0),int(r!=0)
        group=4*po+2*(qo if not o else ro)+(ro if not o else qo)
        rows[group][j]=1 if not o else -1
    return rows


def group_combination(coeff,rows):
    return tuple(sum(c*r[j] for c,r in zip(coeff,rows)) for j in range(D))


def in_occupancy_span(row,rows):
    # Every coordinate belongs to exactly one signed, disjoint occupancy group.
    for r in rows:
        vals={row[j]*r[j] for j in range(D) if r[j]}
        if len(vals)!=1:return False
    return True


def direct_hamilton(word):
    n=len(word)+2;adj=[0]*n;adj[0]|=2
    for i,c in enumerate(word,2):
        for a in range(2):
            if int(c)>>a&1:adj[a]|=1<<i
            else:adj[i]|=1<<a
        for j in range(2,i):adj[j]|=1<<i
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            if dp[mask][i]:
                nxt=adj[i]&~mask
                while nxt:
                    bit=nxt&-nxt;nxt-=bit;j=bit.bit_length()-1
                    dp[mask|bit][j]+=dp[mask][i]
    return sum(dp[-1])


def cyclic_chain_core(word):
    n=len(word)+2;reach=[1<<i for i in range(n)];reach[0]|=2
    for i,c in enumerate(word,2):
        for a in range(2):
            if int(c)>>a&1:reach[a]|=1<<i
            else:reach[i]|=1<<a
        for j in range(2,i):reach[j]|=1<<i
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    keep=[]
    for i in range(2,n):
        if any(j!=i and reach[i]>>j&1 and reach[j]>>i&1 for j in range(n)):
            keep.append(word[i-2])
    return ''.join(keep)


def rejected(thunk):
    try:thunk()
    except ValueError:return True
    return False


def main():
    E=occupancy_rows()
    require(all(dot(r,INITIAL)==0 for r in E),'initial occupancy relation')
    invariant=0
    for r in E:
        for a in range(4):
            require(in_occupancy_span(adjoint(r,a),E),'occupancy invariance')
            invariant+=1
    qx=tuple(int(o==0 and q!=0 and r==0) for o,p,q,r in STATES)
    qy=tuple(int(o==0 and q==0 and r!=0) for o,p,q,r in STATES)
    ell=tuple(int(o==0 and q==0 and r==0) for o,p,q,r in STATES)
    require(dot(ell,INITIAL)==1,'constant initial value')
    require(all(adjoint(ell,a)==ell for a in range(4)),'constant preservation')
    e_coeff=((0,-6,-12,0,0,-6,-12,0),(0,-6,0,0,0,-6,0,0))
    identities=0
    for bit,(numer,den,qrow,eco) in enumerate(zip((X_NUMERATOR,Y_NUMERATOR),
                                  DENOMINATORS,(qx,qy),e_coeff)):
        lhs=group_combination(numer,PROBE_ROWS)
        rhs=tuple(den*qrow[j]+group_combination(eco,E)[j] for j in range(D))
        require(lhs==rhs,'exact 36-coordinate probe identity')
        identities+=D
        require(dot(qrow,INITIAL)==0,'binary code initial value')
        for a in range(4):
            expected=tuple(2*qrow[j]+((a>>bit)&1)*ell[j] for j in range(D))
            require(adjoint(qrow,a)==expected,'binary digital recurrence')
            identities+=D
    total=0;direct=0
    frontier=[('',INITIAL)]
    for length in range(8):
        nxt=[]
        for w,v in frontier:
            p=tuple(dot(r,v) for r in PROBE_ROWS)
            require(extract_codes(p)==direct_codes(w),'code recovery')
            require(recover_word(p)==w.lstrip('0'),'normalized reconstruction')
            if length<=4:
                require(hamilton(w)==direct_hamilton(w),'independent subset DP')
                require(recover_full_profile(full_profile(w))==w.lstrip('0'),'full state decode')
                core=cyclic_chain_core(w)
                require(direct_hamilton(core)==direct_hamilton(w),'singleton SCC deletion')
                require(len(core)<=bounded_spine_length(hamilton(w)),'spine cutoff')
                direct+=1
            if length<7:
                for a in range(4):nxt.append((w+str(a),step(v,a)))
            total+=1
        frontier=nxt
    rng=random.Random(134)
    for length in (8,11,31,100,257,1000):
        for _ in range(4):
            w=''.join(str(rng.randrange(4)) for _ in range(length))
            require(recover_word(probes(w))==w.lstrip('0'),'long reconstruction')
            total+=1
    hostile=0
    base=probes('');full=full_profile('')
    for multiple in (3,7,21):
        candidate=tuple(multiple*x for x in base)
        require(extract_codes(candidate)==(0,0),'hostile integer codes')
        require(rejected(lambda c=candidate:recover_word(c)),'reject positive scaled probes')
        require(rejected(lambda m=multiple:recover_full_profile(tuple(m*x for x in full))),
                'reject positive scaled full state')
        hostile+=2
    for values in (base[:-1],tuple([True]+list(base[1:])),tuple([-1]+list(base[1:]))):
        require(rejected(lambda v=values:recover_word(v)),'reject malformed input')
        hostile+=1
    # The x/y convention and the empty-middle guard are independently visible.
    require(direct_codes('12')==(2,1) and hamilton('12')==5,'letter convention')
    require(hamilton('120')==13 and hamilton('210')==15,'apex arc control')
    for nums in (tuple([X_NUMERATOR[0]+1]+list(X_NUMERATOR[1:])),Y_NUMERATOR):
        lhs=group_combination(nums,PROBE_ROWS)
        rhs=tuple(36*qx[j]+group_combination(e_coeff[0],E)[j] for j in range(D))
        require(lhs!=rhs,'reject mutated extraction identity')
        hostile+=1
    require(hamilton('100311')==25 and bounded_spine_length(25)==6,
            'tight two-apex cutoff control')
    require(hamilton('3333')==1 and hamilton('33330')==47 and hamilton('0')==1,
            'scalar SCC deletion does not preserve continuation profile')
    decisions={}
    for h in (1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31):
        witness=find_count(h)
        require((witness is None)==(h in (7,21)),'bounded scalar decision control')
        if witness is not None:require(direct_hamilton(witness)==h,'decision witness')
        decisions[str(h)]={'bound':bounded_spine_length(h),'word':witness}
    for h in (0,-1,True,1.0):
        require(rejected(lambda t=h:bounded_spine_length(t)),'bound input guard')
        hostile+=1
    print(json.dumps({'status':'PASS','claim':'all-order exact word recovery from twelve continuation probes',
      'probe_suffixes':SUFFIXES,'occupancy_invariance_rows':invariant,
      'exact_coordinate_identity_checks':identities,'actual_word_tests':total,
      'independent_subset_dp_tests':direct,'singleton_scc_cutoff_tests':direct,'hostile_rejections':hostile,
      'finite_max_exhaustive_word_length':7,'finite_random_max_length':1000,
      'ordinary_spectrum_completeness':'OPEN','small_complete_decisions':decisions},indent=2,sort_keys=True))

if __name__=='__main__':main()
