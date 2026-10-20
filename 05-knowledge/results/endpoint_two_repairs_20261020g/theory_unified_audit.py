"""Independent full-partition replay of the endpoint-two theorem at n9,10.

The all-order proofs are distributed across theory.md, last_port.md and
lift.md plus the root open-middle proof. This checker does not prove an
all-order statement by enumeration. Explicit checks survive python -O.
"""
from collections import Counter
from pathlib import Path
import json
from theory_probe import G,endpoint_two
from theory_verify import ledger,replay

HERE=Path(__file__).resolve().parent


def choose(n,u,closed,epsilon):
    a,b,x,c,z=0,1,2,n-2,n-1;y=n-3;D=tuple(range(3,n-2))
    if u==c:
        if epsilon:
            q=(a,x)+D+(c,z,b)
            s=(a,b,x,c)+D[::-1]+(z,)
            v=(b,z)+D+(c,x,a)
            return 'last_port_forward_chord',[q,s[::-1],v]
        return 'last_port_backward_chord',[(a,x)+D+(c,z,b),D+(a,x,c,z,b),D+(a,b,x,c,z)]
    if closed:
        branch,words=ledger(n,u,epsilon)
        return 'closed_'+branch,words
    if u==y:
        return 'open_penultimate_port',[(a,b,c,z,x)+D,(b,c,a,z,x)+D,(b,x)+D+(c,a,z)]
    W=tuple(range(2,n-1))
    return 'open_other_port',[(b,)+W+(z,a),(z,)+W+(b,a),(a,z)+W+(b,)]


def main():
    out=[]
    for n in (9,10):
        counts=Counter();coincident=0
        for middle in G.blue_tilings(n-2):
            closed=bool(middle[-1]&1);epsilon=bool(middle[0]>>(n-4)&1)
            for u in range(2,n-1):
                a=endpoint_two(middle,u)
                G.need(G.is_path(a,tuple(range(n))) and G.blue(a,tuple(range(n))),'marked blue source')
                G.need(a[-1]&1 and a[0].bit_count()==n-3 and a[-1].bit_count()==2,'exact closed endpoint layer')
                branch,words=choose(n,u,closed,epsilon);replay(a,words);counts[branch]+=1
                coincident+=int(2*u==n-1)
        out.append(dict(n=n,total=sum(counts.values()),coincident_port_presentations=coincident,branches=dict(counts)))
    result=dict(status='FINITE-EXACT independent complete layer replay',
                universe='Every standard-blue marked path at orders9,10 with closed endpoint scores(n-3,2); no partner-strongness or mixedness filter',
                results=out)
    (HERE/'theory_unified_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
