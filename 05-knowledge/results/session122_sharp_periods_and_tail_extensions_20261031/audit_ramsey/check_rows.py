#!/usr/bin/env python3
"""Independent replay of exact all-B rows5..8 analytic cutoff certificates.
Imports only this auditor's own geometry/weight routines, never producer code.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import sys
from check import weight, patterns, vertices, solve, require, capacities

BASE={1:(1,(0,)),2:(2,(0,2)),3:(4,(0,3,7,11)),4:(4,(0,4,10,16))}

def ceiling(x):return (x.numerator+x.denominator-1)//x.denominator

def audit(data):
    require(data.get('status')=='PASS','receipt status')
    certs=data['certificates']
    require([c['a'] for c in certs]==[5,6,7,8],'complete ordered row coverage')
    rows=dict(BASE)
    vv=vertices()
    totalfinite=totalpatterns=totalseeds=0
    output=[]
    for c in certs:
        a=c['a'];L=c['period'];hs=tuple(c['residue_values']);s=weight(a)
        require(L==1<<(a-1).bit_length() and len(hs)==L and hs[0]==0,'period/seed format')
        require(Q(c['slope'])==s,'slope')
        rows[a]=(L,hs)
        def row(x,b):
            period,h=rows[x]
            z=weight(x)*period*(b//period)+h[b%period]
            require(z.denominator==1,'integral row')
            return int(z)
        require(all(row(a,b+1)>=row(a,b) for b in range(1,2*L+1)),'row monotonicity')
        D=max(s*r-hs[r] for r in range(L))
        require(Q(c['max_deficit'])==D and D>0,'maximum deficit')
        pp=set(patterns(a,1,True))
        sharp={tuple(t['pattern']):t for t in c['sharp_patterns']}
        other={tuple(t['pattern']):t for t in c['other_patterns']}
        require(len(sharp)==len(c['sharp_patterns']) and len(other)==len(c['other_patterns']),'duplicate pattern')
        require(not(sharp.keys()&other.keys()) and sharp.keys()|other.keys()==pp,'pattern completeness')
        cutoff=1
        for p in sorted(pp):
            co=[weight(x) for x in p]
            half=sum(co)/2
            vertex=max(sum(x*y for x,y in zip(co,v)) for v in vv)
            if half==s:
                require(p in sharp,'sharp pattern label')
                # b-port edges are i,i+2: y_i+y_(i-2)=co_i.
                matrix=[]
                for i in range(5):
                    line=[0]*5;line[i]=line[(i-2)%5]=1;matrix.append(line)
                dual=solve(matrix,co)
                require(dual is not None and min(dual)>0,'strict positive dual')
                require(all((L//2)%rows[x][0]==0 for x in p),'child period compatibility')
                S=ceiling(D/min(dual))-1;K=S+L+1
                item=sharp[p]
                require(item['slack_bound']==S and item['cutoff']==K,'sharp cutoff')
                # Producer lists dual in cycle order; compare multiset as a
                # rotation/permutation invariant and audit our own equations.
                require(sorted(map(Q,item['dual']))==sorted(dual),'dual values')
            else:
                require(p in other and vertex<s,'strict nonsharp gap')
                gap=s-vertex;K=ceiling(D/gap)-1
                require(Q(other[p]['gap'])==gap and other[p]['cutoff']==K,'nonsharp cutoff')
            cutoff=max(cutoff,K);totalpatterns+=1
        joins={z['split']:z for z in c['joins']}
        require(set(joins)==set(range(1,a)) and len(joins)==len(c['joins']),'join coverage')
        for x in range(1,a):
            gap=s-weight(x)-weight(a-x)
            require(gap>0,'strict join gap')
            K=ceiling(D/gap)-1
            require(Q(joins[x]['gap'])==gap and joins[x]['cutoff']==K,'join cutoff')
            cutoff=max(cutoff,K)
        require(c['cutoff']==cutoff,'maximum cutoff')
        require([z[0] for z in c['finite_upper']]==list(range(1,cutoff+1)),'finite coverage')
        for b,claimedc,claimedj,h in c['finite_upper']:
            # Entire maximal port universe, no path-DP algorithm shared.
            cr=max((sum(row(x,y) for x,y in zip(aa,bb)) for aa in pp for bb in patterns(b,2,True)),default=-1)
            j=max(row(x,b)+row(a-x,b) for x in range(1,a))
            require((cr,j,row(a,b))==(claimedc,claimedj,h),'independent finite root values')
            require(cr<=h and j<=h,'root upper bound')
            totalfinite+=1
        require([z['b'] for z in c['seeds']]==list(range(1,L+1)),'seed coverage')
        for se in c['seeds']:
            b=se['b'];pl=se['construction'];aa=pl['a'];bb=pl['b'];kind=pl['kind']
            require(se['a']==a and se['n']==row(a,b),'seed target')
            require(len(aa)==len(bb) and all(x>0 for x in aa+bb),'positive seed children')
            require(all(x<a or x==a and y<b for x,y in zip(aa,bb)),'seed dependency decreases')
            require(sum(row(x,y) for x,y in zip(aa,bb))==se['n'],'seed child order')
            if kind=='c5':
                require(len(aa)==5,'five nonempty children')
                om=max(aa[i]+aa[(i+2)%5] for i in range(5))
                al=max(bb[i]+bb[(i+1)%5] for i in range(5))
            elif kind=='join':
                require(len(aa)==2,'binary join');om=sum(aa);al=max(bb)
            elif kind=='union':
                require(len(aa)==2,'binary union');om=max(aa);al=sum(bb)
            else:raise RuntimeError('unknown seed kind')
            require(om<=a and al<=b,'actual seed budget feasibility');totalseeds+=1
        for r in range(L):
            for t in range(L):
                require(hs[r]+hs[t]<=s*L*((r+t)//L)+hs[(r+t)%L],'all union residues')
        require(c['union_checks']==L*L,'union check count')
        output.append({'a':a,'period':L,'cutoff':cutoff,'maximal_patterns':len(pp),'sharp_patterns':len(sharp)})
    return {'status':'PASS','scope':'independent all-B cutoff-certificate replay for exact rows5..8',
            'rows':output,'finite_root_values':totalfinite,'clique_patterns':totalpatterns,
            'seed_constructions':totalseeds,'union_residue_tests':256}

def main():
    path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent.parent/'ramsey'/'rows.normal.json'
    data=json.loads(path.read_text())
    print(json.dumps(audit(data),sort_keys=True,indent=2))

if __name__=='__main__':main()
