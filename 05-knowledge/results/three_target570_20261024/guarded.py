"""Actual moving176/368 owner equation and complete conditioned charge tables."""
from route import *
from itertools import product
import random

def main():
    rng=random.Random(368176);checks=0;summary=[];example=None
    def check(v,m):
        nonlocal checks
        checks+=1;need(v,m)
    rows={r['source']:r for r in json.loads((HERE/'menus.json').read_text())['rows']}
    rho=[8,9,10,2,6,11,12,13]
    nodes,_,_,_,_=data();M=matrix(nodes[176]['canon_a']);N=matrix(nodes[368]['canon_a'])
    check(M[6][2]==N[4][3]==1,'literal shared pair has coherent source orientation')
    A={o['target']:o['path'] for o in rows[176]['options']};B={o['target']:o['path'] for o in rows[368]['options']}
    table=[]
    for a,b in product((0,1),repeat=2):
        x=A[2234 if a else 1965][2];y=rho[B[2234 if b else 1965][1]]
        check((x==y)==(a!=b),'literal incoming shared-owner equality');table.append([a,b,x,y])
        pa=A[2234 if a else 1965];pb=B[2234 if b else 1965]
        ea={frozenset(z) for z in zip(pa,pa[1:])};eb={frozenset(z) for z in zip(pb,pb[1:])}
        check(frozenset((6,2)) not in ea and frozenset((4,3)) not in eb,'both local partners reverse the shared source pair')
    for remove in (False,True):
        R=Router(remove,True);tested=events=maxcost=maxsupport=0
        for k in range(1,R.n):
            for r in range(max(0,k-len(R.E)),min(20,k)+1):
                feasible=[]
                for bit in (0,1):
                    # Reservoir176=bit, outside368=1-bit.
                    if 0<=r-bit<=19 and 0<=k-r-(1-bit)<=len(R.E)-1:feasible.append(bit)
                if not feasible:continue
                bit=rng.choice(feasible)
                I=rng.sample([i for i in R.R if i!=R.pin],r-bit)
                J=rng.sample([i for i in R.E if i!=R.partner],k-r-(1-bit))
                start=sum(1<<i for i in I+J+[R.pin if bit else R.partner]);q=R.q(start)
                if k<7:finish=R.low[k][q][1]
                elif R.n-k<7:finish=R.full^R.low[R.n-k][R.q(R.full^start)][1]
                else:
                    t=max(7,k-len(R.E));outside=sum(1<<i for i in rng.sample(R.E,k-t));finish=R.replace(outside,t,q)
                P=R.path(start,finish);tested+=1;events+=len(P)-1
                for m in P:
                    a=int(bool(m>>R.pin&1));b=int(bool(m>>R.partner&1))
                    check(A[2234 if a else 1965][2]==rho[B[2234 if b else 1965][1]],'actual shared owner at each arrival')
                    check(m.bit_count()==k and R.q(m)==q,'full margin and charge')
                for a,b in zip(P,P[1:]):
                    maxcost=max(maxcost,R.work(a,b));maxsupport=max(maxsupport,(a^b).bit_count())
                if example is None and len({int(bool(m>>R.pin&1)) for m in P})==2 and len(P)>2:
                    example=R.certificate(start,finish)
        counts=[[len(R.conditional[b].get(k,{})) for k in range(21)] for b in (0,1)]
        # Independently enumerate the 19 unpinned-bank coordinates by Gray tuples.
        RR=[i for i in R.R if i!=R.pin];hist=[set() for _ in range(20)];x=y=k=prev=0
        for j in range(1<<19):
            m=j^(j>>1)
            if j:
                b=(m^prev).bit_length()-1;sgn=1 if m>>b&1 else -1
                ell=R.rows[RR[b]]['label'];y^=ell[0];x=(x+sgn*ell[1])%128;k+=sgn
            hist[k].add((y,x));prev=m
        check(counts[0]==[len(h) for h in hist]+[0],'independent conditional-zero tuple census')
        check(counts[1]==[0]+[len(h) for h in hist],'independent conditional-one translate census')
        summary.append({'binary_sources':R.n,'endpoint_pairs':tested,'events':events,'max_observed_cost':maxcost,
                        'max_observed_support':maxsupport,'conditioned_coverage':counts,
                        'boundary_witness_maxima':[max(c for c,m in R.low[k].values()) for k in range(1,7)]})
    check(example is not None,'nonconstant moving-owner positive control')
    return {'status':'PASS actual moving-owner equality, full conditional tables, and all-grade tests',
            'checks':checks,'literal_equality_table':table,'cases':summary,'example':example}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
