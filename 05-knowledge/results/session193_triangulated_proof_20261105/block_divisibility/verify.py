#!/usr/bin/env python3
"""Exact module-cover divisibility and an explicit two-apex escape.
Standard-library verifier. The endpoint subset DP is adapted from
session193-cut-digits/verify.py with attribution; all new checks are below.
"""
import hashlib,itertools,json,math,random
from pathlib import Path
BASE=Path(__file__).resolve().parent

def need(test,message):
    if not test: raise ValueError(message)
def bits(mask):
    while mask:
        b=mask&-mask; yield b.bit_length()-1; mask-=b

def valid(rows):
    n=len(rows); full=(1<<n)-1
    for i,row in enumerate(rows):
        need(row>=0 and not(row&~full) and not(row>>i&1),'invalid row')
        for j in range(i):need((row>>j&1)+(rows[j]>>i&1)==1,'not tournament')
def reverse(rows):return [((1<<len(rows))-1)^row^(1<<i) for i,row in enumerate(rows)]
def tables(rows):
    n=len(rows); size=1<<n; ends=[[0]*n for _ in range(size)]; H=[0]*size;H[0]=1
    for i in range(n):ends[1<<i][i]=1
    for mask in range(1,size):
        line=ends[mask]
        for i in bits(mask):
            count=line[i];H[mask]+=count
            if count:
                for j in bits(rows[i]&~mask):ends[mask|1<<j][j]+=count
    return H,ends
def hp(rows):return tables(rows)[0][-1]
def literal(rows):
    return sum(all(rows[u]>>v&1 for u,v in zip(p,p[1:])) for p in itertools.permutations(range(len(rows))))
def strong(rows):
    for rr in (rows,reverse(rows)):
        seen=1
        while True:
            nxt=seen
            for i in bits(seen):nxt|=rr[i]
            if nxt==seen:break
            seen=nxt
        if seen!=(1<<len(rows))-1:return False
    return True

def join(blocks):
    rows=[]; total=sum(map(len,blocks)); offset=0
    for block in blocks:
        later=((1<<total)-1)^((1<<(offset+len(block)))-1)
        rows += [(x<<offset)|later for x in block];offset+=len(block)
    valid(rows);return rows

def apex(rows,mask):
    n=len(rows); out=[a|((1-(mask>>i&1))<<n)for i,a in enumerate(rows)]+[mask]
    valid(out);return out

def covers(rows,qmax):
    H,_=tables(rows);size=len(H);last=[0]*size;last[0]=1;Bs=[]
    for q in range(qmax):
        nxt=[0]*size
        for M in range(size):
            A=M
            while True:
                nxt[M]+=H[A]*last[M^A]
                if A==0:break
                A=(A-1)&M
        Bs.append(nxt[-1]);last=nxt
    return Bs

def nonempty_cover_coefficients(rows):
    B=[0]+covers(rows,len(rows));out=[]
    for k in range(1,len(rows)+1):
        out.append(sum((-1)**(k-j)*math.comb(k,j)*B[j]for j in range(k+1)))
    return out

def profile(rows):
    H,E=tables(rows);_,S=tables(reverse(rows));full=(1<<len(rows))-1
    R=[sum(S[A][i]*H[full^A]for A in range(1,full+1)if A>>i&1)for i in range(len(rows))]
    L=[sum(E[A][i]*H[full^A]for A in range(1,full+1)if A>>i&1)for i in range(len(rows))]
    C=[[0]*len(rows)for _ in rows]
    for i in range(len(rows)):
        for j in range(i):
            c=0
            for A in range(1,full):
                if A>>i&1 and not(A>>j&1):c+=E[A][i]*S[full^A][j]
                if A>>j&1 and not(A>>i&1):c+=E[A][j]*S[full^A][i]
            C[i][j]=C[j][i]=c
    return dict(h=H[-1],B=sum(H[A]*H[full^A]for A in range(full+1)),R=R,L=L,e=E[-1],s=S[-1],C=C)

def joined_capacities(blocks):
    ps=[profile(b)for b in blocks]; offsets=[0]
    for b in blocks:offsets.append(offsets[-1]+len(b))
    C=[[0]*offsets[-1]for _ in range(offsets[-1])];P=math.prod(p['h']for p in ps)
    for a,p in enumerate(ps):
        for i in range(len(blocks[a])):
            for j in range(i):C[offsets[a]+i][offsets[a]+j]=C[offsets[a]+j][offsets[a]+i]=(P//p['h'])*p['C'][i][j]
        for b in range(a+1,len(blocks)):
            O=math.prod(ps[t]['h']for t in range(len(blocks))if t<a or t>b)
            M=math.prod(ps[t]['B']for t in range(a+1,b))
            for i in range(len(blocks[a])):
                for j in range(len(blocks[b])):
                    c=O*(p['R'][i]*ps[b]['L'][j]*M+(p['e'][i]*ps[b]['s'][j]if b==a+1 else 0))
                    C[offsets[a]+i][offsets[b]+j]=C[offsets[b]+j][offsets[a]+i]=c
    return P,C

def cut_count(rows,h,C,mask):
    return h+sum(C[i][j]for i in bits(mask)for j in bits(rows[i]&~mask))

def divisor(gs,r):
    # gcd of all products after exactly min(2r,m) omissions; equivalent to
    # subtracting the largest 2r valuations prime by prime.
    keep=max(0,len(gs)-2*r)
    return math.gcd(*(math.prod(gs[i]for i in I)for I in itertools.combinations(range(len(gs)),keep)))

def uniform_two_apex_count(block,word):
    # Independent physical run-allocation transfer for entire modules.
    # State (order,x,y,z) as in S193 counters, but nonempty subsets of slots
    # are filled simultaneously, weighted by exact nonempty cover counts.
    A=nonempty_cover_coefficients(block); A += [0]*(3-len(A))
    states={(0,0,0,0):1,(1,0,0,0):1}
    for c in word:
        nxt={}
        for (o,x,y,z),v in states.items():
            p=(c>>o)&1;q=(c>>(1-o))&1
            for used in range(1,8):
                if used&2 and not(y or p):continue
                if used&4 and not(z or q):continue
                key=(o,(2 if p else 1)if used&1 else x,(2 if q else 1)if used&2 else y,1 if used&4 else z)
                nxt[key]=nxt.get(key,0)+v*A[used.bit_count()-1]
        states=nxt
    return sum(v for (o,x,y,z),v in states.items()if x!=2 and (y==1 or(y==0 and o==0)))

def uniform_graph(block,word):
    rows=join([block]*len(word));n=len(rows);a=n;b=n+1
    rows += [1<<b,0]
    for k,c in enumerate(word):
        for i in range(k*len(block),(k+1)*len(block)):
            for t in range(2):
                if c>>t&1:rows[n+t]|=1<<i
                else:rows[i]|=1<<(n+t)
    valid(rows);return rows

C3=[2,4,1];C5=[6,12,24,17,3]
V9=[20,1,2,7,14];V13=[18,8,3,5,14];V25=[2,12,9,49,7,23]

def main():
    report={'status':'PASS: all-order proofs are in PROOF.md; enumerations have named finite scope'}
    inventories={}
    for name,rows in [('C3',C3),('C5',C5),('V9',V9),('V13',V13),('V25',V25)]:
        valid(rows);B=covers(rows,max(8,len(rows)));A=nonempty_cover_coefficients(rows)
        need(all(x>0 for x in A),'nonempty cover positivity')
        need(A[-1]==math.factorial(len(rows)),'singleton path covers')
        need(B[1]>2**len(rows),'strict nontransitive entropy margin')
        for j,b in enumerate(B,1):need(b==sum(math.comb(j,k)*a for k,a in enumerate(A,1)if k<=j),'binomial cover expansion')
        inventories[name]={'rows':rows,'H':B[0],'B_1_through_8':B,'A_nonempty':A,'g_infinity':math.gcd(*A),'g_1':math.gcd(*B[:2]),'g_2':math.gcd(*B[:3])}
    need(inventories['C3']['A_nonempty']==[3,6,6],'C3 profile')
    need(inventories['C5']['A_nonempty']==[15,90,210,240,120],'C5 profile')
    need(inventories['V9']['g_infinity']==3,'asymmetric persistent divisor')
    need(inventories['V13']['g_1']==13 and inventories['V13']['g_2']==1,'second apex profile escape')
    need(inventories['V25']['g_1']==1,'gcd-one target survives')
    report['module_profiles']=inventories
    record_counts=[]
    for name,rows in [('C3',C3),('C5',C5),('V13',V13),('V25',V25)]:
        p=profile(rows);records=set()
        for mask in range(1<<len(rows)):
            records.add((cut_count(rows,p['h'],p['C'],mask)-p['h'],)+tuple(sum(p[key][i]for i in bits(mask))for key in ['R','L','e','s']))
        record_counts.append({'core':name,'q':len(records),'B2':p['B']})
    need([z['q']for z in record_counts]==[4,8,32,64],'lossless mixed-library record counts')
    from fractions import Fraction
    total=sum(Fraction(z['q'],z['B2'])for z in record_counts)
    need(total==Fraction(1134,1235),'mixed-library reciprocal ratio')
    need(17**100>304*16**100 and Fraction(17,16)*total<1,'exact exponent 99/100 certificate')
    report['mixed_core_sparsity_certificate']={'records':record_counts,'sum_q_over_B2':str(total),'exponent':'99/100','integer_power_check':True,'ratio_upper_bound':str(Fraction(17,16)*total)}
    # Exhaustive one-apex cut universes, independently checking sampled actual
    # Hamilton counts (all small cases) and every full capacity/cut expression.
    cases=[('C3^1',[C3]),('C3^2',[C3,C3]),('C3^3',[C3,C3,C3]),('C3-singleton-C3-C3',[C3,[0],C3,C3]),('V9^2',[V9,V9]),('C5^2',[C5,C5])]
    cut_reports=[];direct=0; rng=random.Random(193)
    for name,blocks in cases:
        rows=join(blocks);h,C=joined_capacities(blocks)
        D=divisor([math.gcd(p['h'],p['B'])for p in map(profile,blocks)],1)
        values=[]
        checks=set(range(1<<len(rows)))if len(rows)<=6 else {0,(1<<len(rows))-1,*[rng.randrange(1<<len(rows))for _ in range(10)]}
        for mask in range(1<<len(rows)):
            H=cut_count(rows,h,C,mask);need(H%D==0,'all-face divisor');values.append(H)
            if len(blocks)>=2 and all(block==blocks[0]for block in blocks):
                block=blocks[0];width=len(block);last=((1<<width)-1)<<(len(rows)-width)
                if mask&((1<<width)-1) and (mask&last)!=last:
                    pp=profile(block);need(H>=pp['h']**2*pp['B']**(len(blocks)-2),'strong all-cut growth lower bound')
            if mask in checks:
                need(H==hp(apex(rows,mask)),'physical subset-DP cut count');direct+=1
        cut_reports.append({'name':name,'base_order':len(rows),'cuts':len(values),'predicted_divisor':D,'actual_value_gcd':math.gcd(*values),'direct_DP_checks':len(checks)})
    report['one_apex_all_face_universes']=cut_reports
    # Sharpness: choose one vertex in each C3 block, except change the last
    # cardinality to 2 when m=1 mod3. All these actual cuts are strong.
    sharp=[]
    for m in range(2,13):
        blocks=[C3]*m;rows=join(blocks);h,C=joined_capacities(blocks)
        sizes=[1]*m
        if m%3==1:sizes[-1]=2
        mask=sum(((1<<k)-1)<<(3*i)for i,k in enumerate(sizes));H=cut_count(rows,h,C,mask)
        need(H%3**(m-2)==0 and H%3**(m-1)!=0,'sharp exponent')
        need(strong(apex(rows,mask)),'sharp strong cut')
        sharp.append({'m':m,'cardinalities':sizes,'H':H,'v3':m-2})
    need(sharp[0]['H']==41,'two-block prime witness')
    report['sharp_C3_examples']=sharp
    # Arbitrary two-apex incidences, not uniform block patterns. These physical
    # witnesses test the general 2r boundary count beyond the cut identity.
    multi=[]
    for m,repetitions in [(3,2),(4,2),(5,2)]:
        for sample in range(repetitions):
            rows=join([C3]*m);rows=apex(rows,rng.randrange(1<<len(rows)));rows=apex(rows,rng.randrange(1<<len(rows)))
            H=hp(rows);D=divisor([3]*m,2);need(H%D==0,'arbitrary two-apex divisor');direct+=1
            multi.append({'m':m,'rows':rows,'H':H,'divisor':D})
    report['arbitrary_two_apex_physical_checks']=multi
    # Unbounded physical escape from the H13 one-apex divisor. Complete blocks
    # have bits3 then bits0; transfer is independent of full-graph subset DP.
    escape=[]
    for m in range(2,15):
        word=[3]+[0]*(m-1);H=uniform_two_apex_count(V13,word)
        formula=552*471**(m-1)-675*104**(m-1)-7*13**(m-1)
        need(H==formula,'all-clock escape formula finite controls')
        need(H%13==2*pow(3,m,13)%13!=0,'escape residue')
        rows=uniform_graph(V13,word);need(strong(rows),'escape strongness')
        if m<=3:need(hp(rows)==H,'independent physical escape DP');direct+=1
        escape.append({'m':m,'H':H,'residue_mod13':H%13,'order':5*m+2})
    # The same module transfer is checked on every word of lengths one/two for
    # blocks with one/two/three vertices; literal paths provide a third route.
    transfer_checks=0;literal_checks=0
    for block in [[0],[2,0],C3]:
        for m in (1,2):
            for w in itertools.product(range(4),repeat=m):
                rows=uniform_graph(block,w);H=uniform_two_apex_count(block,w)
                need(H==hp(rows),'module transfer physical check');transfer_checks+=1
                if len(rows)<=6:need(H==literal(rows),'literal permutation control');literal_checks+=1
    report['two_apex_escape']=escape
    report['physical_subset_DP_checks']=direct+transfer_checks
    report['module_transfer_checks']=transfer_checks
    report['literal_permutation_checks']=literal_checks
    report['source_sha256']={'verify.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
