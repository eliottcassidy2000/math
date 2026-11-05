#!/usr/bin/env python3
"""Exact cut-cube intervals and ordered-block sumset certificates. Stdlib only."""
from itertools import combinations, permutations, product
from math import gcd, lcm, prod
from pathlib import Path
import hashlib,json

def need(x,msg):
    if not x: raise RuntimeError(msg)
def bits(x):
    while x:
        b=x&-x;yield b.bit_length()-1;x-=b

def validate(r):
    n=len(r);full=(1<<n)-1
    for i in range(n):
        need(not(r[i]>>i&1) and not(r[i]&~full),'invalid row')
        for j in range(i):need(((r[i]>>j)&1)+((r[j]>>i)&1)==1,'invalid pair')
def rows(n,code):
    r=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        a,b=(i,j)if code>>k&1 else(j,i);r[a]|=1<<b
    return r

def extend(r,S):
    n=len(r);return [a|((not(S>>i&1))<<n)for i,a in enumerate(r)]+[S]
def tables(r):
    n=len(r);N=1<<n
    e=[[0]*n for _ in range(N)];s=[[0]*n for _ in range(N)]
    inc=[sum(1<<j for j in range(n)if r[j]>>i&1)for i in range(n)]
    for i in range(n):s[1<<i][i]=e[1<<i][i]=1
    for S in range(1,N):
        for i in bits(S):
            A=S^(1<<i)
            if A:
                e[S][i]=sum(e[A][j]for j in bits(A&inc[i]))
                s[S][i]=sum(s[A][j]for j in bits(A&r[i]))
    return s,e

def hp(r):return sum(tables(r)[1][-1]) if r else 1

def data(r):
    validate(r);n=len(r);N=1<<n;full=N-1;s,e=tables(r)
    H=[1]+[sum(e[S])for S in range(1,N)]
    B=sum(H[S]*H[full^S]for S in range(N))
    R=[sum(s[S][i]*H[full^S]for S in range(1,N))for i in range(n)]
    L=[sum(e[S][i]*H[full^S]for S in range(1,N))for i in range(n)]
    Q=[[0]*n for _ in r]
    for S in range(1,full):
        for i in bits(S):
            for j in bits(full^S):Q[i][j]+=e[S][i]*s[full^S][j]
    C=[[Q[i][j]+Q[j][i]for j in range(n)]for i in range(n)]
    need(all(C[i][j]>0 and C[i][j]%2==0 for i,j in combinations(range(n),2)),'positive even capacity')
    return {'r':r,'h':H[-1],'B':B,'s':s[-1],'e':e[-1],'R':R,'L':L,'C':C}

def cut(d,S):return sum(d['C'][i][j]for i in bits(S)for j in bits(d['r'][i]&~S))
def strong(r):
    full=(1<<len(r))-1
    for i in range(len(r)):
        S=1<<i
        while True:
            T=S
            for j in bits(S):T|=r[j]
            if T==S:break
            S=T
        if S!=full:return False
    return True

def literal(r):
    n=len(r);H=0;Q=[[0]*n for _ in r]
    for p in permutations(range(n)):
        bad=[i for i in range(n-1)if not(r[p[i]]>>p[i+1]&1)]
        H+=not bad
        for i in range(n-1):
            if not bad or bad==[i]:Q[p[i]][p[i+1]]+=1
    return H,[[Q[i][j]+Q[j][i]for j in range(n)]for i in range(n)]

def join(blocks):
    n=sum(map(len,blocks));r=[];off=0
    for b in blocks:
        later=((1<<n)-1)^((1<<(off+len(b)))-1)
        r.extend((a<<off)|later for a in b);off+=len(b)
    return r

def joined_formula(ds):
    sizes=[len(d['r'])for d in ds];P=prod(d['h']for d in ds)
    offsets=[sum(sizes[:i])for i in range(len(ds))];n=sum(sizes);C=[[0]*n for _ in range(n)]
    for a,d in enumerate(ds):
        for i in range(sizes[a]):
            for j in range(sizes[a]):C[offsets[a]+i][offsets[a]+j]=(P//d['h'])*d['C'][i][j]
    for a,b in combinations(range(len(ds)),2):
        outside=prod(ds[t]['h']for t in range(len(ds))if t<a or t>b)
        middle=prod(ds[t]['B']for t in range(a+1,b))
        for i in range(sizes[a]):
            for j in range(sizes[b]):
                c=outside*(ds[a]['R'][i]*ds[b]['L'][j]*middle+(ds[a]['e'][i]*ds[b]['s'][j]if b==a+1 else 0))
                C[offsets[a]+i][offsets[b]+j]=C[offsets[b]+j][offsets[a]+i]=c
    return P,C

def cube(d,base,pairs):
    C=d['C'];need(len(set(sum(([a,b]for a,b in pairs),[])))==2*len(pairs),'overlapping exchanges')
    for a,b in pairs:need((base>>a&1)and not(base>>b&1),'exchange orientation')
    for(a,b),(c,e)in combinations(pairs,2):need(C[a][c]-C[a][e]-C[b][c]+C[b][e]==0,'nonzero mixed interaction')
    F=d['h']+cut(d,base)
    steps=[cut(d,base^(1<<a)^(1<<b))-cut(d,base)for a,b in pairs]
    vals=[];masks=[]
    for es in product((0,1),repeat=len(pairs)):
        S=base
        for x,(a,b)in zip(es,pairs):
            if x:S^=(1<<a)|(1<<b)
        v=d['h']+cut(d,S);need(v==F+sum(x*b for x,b in zip(es,steps)),'cube formula')
        vals.append(v);masks.append(S)
    coins=sorted(abs(v)//2 for v in steps if v);reach=0
    for c in coins:
        if c>reach+1:break
        reach+=c
    return {'base_mask':base,'pairs':pairs,'base_H':F,'steps':steps,'masks':masks,'support':sorted(set(vals)),'covered_normalized_prefix':reach}

def main():
    report={'status':'FINITE-EXACT checks; all-order claims have separate proofs','counts':{}}
    parent=[50,8,3,5,14,30];d=data(parent);need(strong(parent),'cube parent strong')
    need(literal(parent)==(d['h'],d['C']),'literal cube parent')
    cr=cube(d,37,[(0,4),(2,3)])
    need(cr['support']==[129,131,133,135],'four-odd cube interval')
    for S in cr['masks']:
        child=extend(parent,S);need(strong(child),'strong cube child')
        need(hp(child)==literal(child)[0]==d['h']+cut(d,S),'literal cube child')
    cr.update({'rows':parent,'H':d['h'],'capacities':d['C']});report['interval_cube']=cr
    # H13 remains holey, despite having unit gcd of energy differences.
    h13=data([6,8,26,17,3]);support=sorted({h13['h']+cut(h13,S)for S in range(1,31)})
    need(support==[23,25,29,31,33,37,43,45],'H13 support')
    need(gcd(*( (v-support[0])//2 for v in support))==1,'H13 gcd one')
    report['H13_hostile']={'support':support,'missing_between_extrema':[v for v in range(23,46,2)if v not in support]}
    # Independent universal finite check: every ordered pair of labelled blocks
    # of orders 1..3 (121 pairs), and all eight ordered triples from {1,C3}.
    bank=[data(rows(n,c))for n in range(1,4)for c in range(1<<(n*(n-1)//2))]
    pair_checks=0;child_checks=0
    for a,b in product(bank,repeat=2):
        rr=join([a['r'],b['r']]);actual=data(rr);expected=joined_formula([a,b])
        need((actual['h'],actual['C'])==expected,'ordered-pair capacity law')
        for S in range(1<<len(rr)):
            need(hp(extend(rr,S))==actual['h']+cut(actual,S),'ordered-pair child direct DP');child_checks+=1
        pair_checks+=1
    one=data([0]);c3=data([2,4,1])
    for ds in product((one,c3),repeat=3):
        rr=join([d['r']for d in ds]);actual=data(rr)
        need((actual['h'],actual['C'])==joined_formula(ds),'ordered-triple capacity law')
    report['counts'].update({'ordered_block_pairs':pair_checks,'ordered_block_pair_children':child_checks,'ordered_block_triples':8})
    c5rows=[6,12,24,17,3];c5=data(c5rows)
    need((c5['h'],c5['B'],c5['R'],c5['L'],c5['s'],c5['e'])==(15,120,[21]*5,[21]*5,[3]*5,[3]*5),'C5 interface')
    support5=sorted({cut(c5,S)//2 for S in range(32)if S.bit_count()==2})
    need(support5==[13,14],'C5 size-two support')
    report['C5']={'rows':c5rows,'H':15,'B':120,'R_and_L':21,'start_and_end':3,'size_two_energy_support':support5}
    repeated=[]
    for m in range(1,4):
        rr=join([c5rows]*m);h,C=joined_formula([c5]*m);dd={'r':rr,'h':h,'C':C}
        if m<=2:need((data(rr)['h'],data(rr)['C'])==(h,C),'repeated C5 direct capacity')
        # Exhaust all 10^m local size-two choices, via proved matrix formula.
        choices=[s for s in range(32)if s.bit_count()==2]
        vals=set()
        for ss in product(choices,repeat=m):vals.add(h+cut(dd,sum(S<<(5*j)for j,S in enumerate(ss))))
        Ksum=sum(15**(m-2)*450 if gap==1 else 441*120**(gap-1)*15**(m-gap-1) for gap in range(1,m)for _ in range(m-gap))
        base=15**m+26*m*15**(m-1)+6*Ksum
        expected=[base+2*15**(m-1)*t for t in range(m+1)]
        need(sorted(vals)==expected,'repeated C5 exact progression')
        for t in range(m+1):
            S=sum((5 if j<t else 3)<<(5*j)for j in range(m));child=extend(rr,S)
            need(hp(child)==expected[t]and strong(child),'C5 progression direct child DP')
        g=15**(m-1);need(gcd(*(v-expected[0]for v in expected))==2*g,'cofactor lattice')
        repeated.append({'m':m,'support':expected,'difference_gcd':2*g,'all_local_size_two_masks':10**m})
    report['repeated_C5']=repeated
    padded=[]
    for t in range(9):
        rr=join([c5rows]+[[0]]*t);ds=[c5]+[one]*t;h,C=joined_formula(ds);dd={'r':rr,'h':h,'C':C}
        vals=[hp(extend(rr,S))for S in (3,5)]
        expected=[41,43]if t==0 else[42*(1<<t)+5,42*(1<<t)+7]
        need(vals==expected,'padded C5 consecutive pair')
        need(all(h+cut(dd,S)==expected[i]and strong(extend(rr,S))for i,S in enumerate((3,5))),'padded C5 cut and strongness')
        padded.append({'sink_padding':t,'order':6+t,'H':vals})
    report['padded_consecutive_pairs']=padded
    dyadic_records=[([20,1,2,7,14],10,29,32,4,1),([20,1,2,7,14],17,33,32,2,3),([18,8,3,5,14],6,31,32,6,5),([18,8,3,5,14],3,33,32,6,7)]
    dyadic=[]
    for rr,S,H0,A,E,b in dyadic_records:
        dd=data(rr);need(strong(rr),'dyadic core strong')
        need(literal(rr)==(dd['h'],dd['C']),'dyadic core literal capacities')
        need(literal(extend(rr,S))[0]==H0==dd['h']+cut(dd,S),'dyadic seed child literal')
        need(sum(dd['R'][i]for i in bits(S))==A and sum(dd['e'][i]for i in bits(S))==E,'dyadic boundary signature')
        need(H0-A+E==b,'dyadic intercept')
        vals=[]
        for k in range(6,13):
            t=k-5;parent=join([rr]+[[0]]*t);child=extend(parent,S)
            need(len(child)==k+1 and strong(child),'dyadic strong order')
            v=hp(child);need(v==(1<<k)+b==H0+A*((1<<t)-1)+E,'dyadic exact padding law')
            vals.append({'k':k,'H':v})
        dyadic.append({'rows':rr,'mask':S,'H_child':H0,'R_sum':A,'end_sum':E,'intercept':b,'verified_children':vals})
    report['dyadic_four_odd_window']=dyadic
    # The generic cofactor identity, without invoking real-valued normalization.
    for hs in product((1,3,15,175,189,3267),repeat=4):
        P=prod(hs);need(gcd(*(P//h for h in hs))==P//lcm(*hs),'cofactor gcd identity')
    report['counts']['cofactor_gcd_tuples']=6**4
    report['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
