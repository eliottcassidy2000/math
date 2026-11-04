#!/usr/bin/env python3
"""Exhaust original root-containing module families on C3[TT5,TT2,TT2]."""
from pathlib import Path
from itertools import combinations
import json,subprocess,hashlib
import bank_audit as B
P=Path(__file__).resolve().parent

def tt(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def op(q):return tuple(((1<<len(q))-1)^r^(1<<i) for i,r in enumerate(q))
def sub(Q,D):
    blocks=[];n=0
    for d in D:blocks.append(list(range(n,n+len(d))));n+=len(d)
    out=[0]*n
    for i,d in enumerate(D):
        for u,row in enumerate(d):
            for v in range(len(d)):
                if row>>v&1:out[blocks[i][u]]|=1<<blocks[i][v]
        for j,b in enumerate(blocks):
            if Q[i]>>j&1:
                for u in blocks[i]:out[u]|=sum(1<<v for v in b)
    return tuple(out),blocks

def H(rows):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for mask in range(1,1<<n):
        for j in range(n):
            if mask>>j&1:
                prev=mask^(1<<j);dp[mask][j]=1 if not prev else sum(dp[prev][i] for i in range(n)if prev>>i&1 and rows[i]>>j&1)
    return sum(dp[-1])
def basis(n):return [0]+[1<<i for i in range(n)]+[s for s in range(1<<n)if s.bit_count()==2]
def module(rows,S):
    names=[i for i in range(len(rows))if S>>i&1]
    return all(len({rows[v]>>i&1 for i in names})==1 for v in range(len(rows))if not S>>v&1)
def toggle(rows,sets):
    out=list(rows)
    for S in sets:
        for i,j in combinations([i for i in range(len(rows))if S>>i&1],2):out[i]^=1<<j;out[j]^=1<<i
    return tuple(out)
def transp(mask,p):return sum(1<<p[i]for i in range(len(p))if mask>>i&1)
def evaluate(bank,a,c):
    m=bank['masks'];co=bank['coefficients'];w=len(m)
    return sum(co[i*w+j]for i,u in enumerate(m)if u&a==u for j,v in enumerate(m)if v&c==v)
def coeff(bank,a,c):
    w=len(bank['masks']);ix={s:i for i,s in enumerate(bank['masks'])};return bank['coefficients'][ix[a]*w+ix[c]]

def main():
    Q=(2,4,1);T,blocks=sub(Q,[tt(5),tt(2),tt(2)]);U,_=sub(op(Q),[tt(5),tt(2),tt(2)]);n=len(T);contexts=[];graphs=[]
    for ambient in [False,True]:
        if ambient:
            src,_=sub(Q,[T,(0,),(0,)]);tar,_=sub(Q,[U,(0,),(0,)]);outside=1<<(n+1)
        else:src,tar,outside=T,U,0
        contexts.append({'ambient':ambient,'source_rows':src,'target_canonical_rows':tar,'outside_one_mask':outside,'graph_start':len(graphs)});graphs.extend([src,tar])
    queries=[];tags={}
    def include(a,c,tag):
        if (a,c) not in tags:queries.append((a,c));tags[a,c]=tag
    include(0,0,'scalar')
    quartics=[]
    for i in range(5):
        for j,k in combinations([j for j in range(5)if j!=i],2):
            a=(1<<i)|(1<<j);c=(1<<i)|(1<<k);quartics.append((a,c))
            for s in B.subs(a):
                for t in B.subs(c):include(s,t,'internal_quartic')
    groups=[sum(1<<i for i in b)for b in blocks[1:]]
    for r in [4,2]:
        for G in groups:
            for a in [0,1,2,3]:
                for c in [0,1<<r,G,G|(1<<r)]:include(a,c,'mixed_quartic')
    compressed=[[i]for i in range(5)]+blocks[1:];d=0;anchor=[]
    for C in compressed[1:]:
        a=sum(1<<v for v in C)
        if T[0]>>C[0]&1:a|=1;d+=1
        anchor.append((a,1))
    if d%2:anchor.append((1,1))
    if (len(compressed)-1+d)%2:anchor.append((0,1))
    anchor.append((1,0))
    for a,c in anchor:include(a,c,'anchor')
    intervals=[sum(1<<i for i in range(l,r+1))for l in range(5)for r in range(l+1,5)];selected_modules=intervals+groups
    B.need(all(module(T,S)for S in selected_modules),'every selected set original module');B.need(len(intervals)==10 and len(set(intervals))==10,'all source intervals')
    # Recount four complete physical named banks.
    input_lines=[];banks=[]
    for g in graphs:
        bs=basis(len(g));start=len(input_lines)
        for a in bs:
            for c in bs:
                aug=B.ports(g,a,c);input_lines.append(' '.join(map(str,[len(aug),*aug])))
        banks.append({'rows':g,'masks':bs,'query_start':start})
    # One scalar physical query for every possible internal tournament, in both contexts.
    families=[];scalar_start=len(input_lines);seen=set()
    for fmask in range(1<<10):
        A=toggle(op(tt(5)),[S for j,S in enumerate(intervals)if fmask>>j&1]);B.need(A not in seen,'interval XOR enumeration is injective');seen.add(A);target,_=sub(op(Q),[A,tt(2),tt(2)]);hc=H(A);order=sorted(range(5),key=lambda i:-A[i].bit_count());istrans=all(A[order[i]]>>order[j]&1 for i in range(5)for j in range(i+1,5));B.need((hc==1)==istrans,'scalar equality/transitive boundary')
        record={'interval_mask':fmask,'target_child_rows':A,'child_H':hc,'transitive_order':order if istrans else None,'scalar_start':len(input_lines)}
        for C in contexts:
            ambient_target=sub(Q,[target,(0,),(0,)])[0] if C['ambient']else target;z=B.ports(ambient_target,C['outside_one_mask'],C['outside_one_mask']);input_lines.append(' '.join(map(str,[len(z),*z])))
        families.append(record)
    subprocess.run(['g++','-O3','-std=c++17',str(P/'count.cpp'),'-o',str(P/'count')],check=True);output=subprocess.check_output([str(P/'count')],input='\n'.join(input_lines)+'\n',text=True);values=list(map(int,output.split()));B.need(len(values)==len(input_lines),'family full-bank/count outputs')
    for bank in banks:
        bs=bank['masks'];w=len(bs);start=bank.pop('query_start');vals=values[start:start+w*w];ix={m:i for i,m in enumerate(bs)};co=[sum((-1)**((a^s).bit_count()+(c^t).bit_count())*vals[ix[s]*w+ix[t]] for s in B.subs(a)for t in B.subs(c))for a in bs for c in bs];bank['values']=vals;bank['coefficients']=co
    statistics={};physical_checks=[]
    for F in families:
        A=F['target_child_rows'];F['contexts']=[]
        for ci,C in enumerate(contexts):
            src=banks[2*ci];tar=banks[2*ci+1];ext=C['outside_one_mask'];s0=evaluate(src,ext,ext);t0=values[F['scalar_start']+ci];B.need(t0>=s0,'scalar monotonicity')
            if F['child_H']!=1:
                B.need(t0>s0,'strict scalar target');F['contexts'].append({'query':[ext,ext],'source':s0,'target':t0,'branch':'scalar'});continue
            B.need(t0==s0,'profile-equal transitive scalar');p=list(range(len(C['source_rows'])))
            for rank,v in enumerate(F['transitive_order']):p[v]=rank
            internal_equal=all(coeff(src,a,c)==coeff(tar,transp(a,p),transp(c,p))for a,c in quartics)
            original=tuple(A)==tt(5);reversed_=tuple(A)==op(tt(5));B.need(internal_equal==(original or reversed_),'exact order reconstruction up to converse')
            if reversed_:
                choices=anchor;branch='conversed_child_anchor'
            elif original:
                choices=[z for z in queries if tags[z]=='mixed_quartic'];branch='unchanged_child_mixed'
            else:
                choices=[z for z in queries if tags[z]=='internal_quartic'];branch='internal_order'
            chosen=None
            for a,c in choices:
                aa=a|ext;cc=c|ext;s=evaluate(src,aa,cc);t=evaluate(tar,transp(aa,p),transp(cc,p))
                if s!=t:
                    chosen={'query':[aa,cc],'transported_target_query':[transp(aa,p),transp(cc,p)],'transport':p,'source':s,'target':t,'branch':branch};break
            B.need(chosen is not None,'fixed branch battery detects');F['contexts'].append(chosen)
            if reversed_:B.need(any((evaluate(src,a|ext,c|ext)-evaluate(tar,transp(a|ext,p),transp(c|ext,p)))%4==2 for a,c in anchor),'anchor residue2')
        # Each of the two pair modules may independently be selected. Query masks are uniform on both.
        F['pair_transports']=[]
        for pairmask in range(4):
            chosen=[(1<<n)-1]+[S for j,S in enumerate(intervals)if F['interval_mask']>>j&1]+[S for j,S in enumerate(groups)if pairmask>>j&1]
            target=toggle(T,chosen)
            for j,G in enumerate(groups):
                b=blocks[j+1];source_arc=target[b[0]]>>b[1]&1;expected=bool(pairmask>>j&1);B.need(source_arc==expected,'actual pair descendant XOR')
            B.need(tuple(x&31 for x in target[:5])==tuple(A),'actual original module XOR target')
            transport=list(range(n))
            if F['child_H']==1:
                for rank,v in enumerate(F['transitive_order']):transport[v]=rank
            for j,b in enumerate(blocks[1:]):
                if not(pairmask>>j&1):transport[b[0]],transport[b[1]]=transport[b[1]],transport[b[0]]
            canonical=sub(op(Q),[tt(5) if F['child_H']==1 else A,tt(2),tt(2)])[0]
            B.need(all(((target[i]>>j)&1)==((canonical[transport[i]]>>transport[j])&1)for i in range(n)for j in range(n)),'complete physical target transport')
            F['pair_transports'].append(transport)
            for rec in F['contexts']:
                for a in rec['query']:B.need(all((a&G) in [0,G]for G in groups),'physical uniform pair query')
        for rec in F['contexts']:statistics[rec['branch']]=statistics.get(rec['branch'],0)+4
    data={'source_core_rows':T,'source_modules':selected_modules,'source_prime_quotient':Q,'child_blocks':blocks,'fixed_query_battery':queries,'anchor_battery':anchor,'internal_quartics':quartics,'contexts':contexts,'families':families,'pair_selection_options_per_record':4,'statistics':statistics}
    (P/'family_controls.json').write_text(json.dumps(data,separators=(',',':'))+'\n');(P/'family_banks.json').write_text(json.dumps(banks,separators=(',',':'))+'\n')
    receipt={'status':'PASS','original_root_families':4*len(families),'ambient_contexts':len(contexts),'family_context_separations':4*len(families)*len(contexts),'fixed_battery_size':len(queries),'full_physical_named_banks':len(banks),'full_named_coefficients':sum(len(b['coefficients'])for b in banks),'literal_graph_queries':len(input_lines),'branch_counts':statistics,'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['family_audit.py','bank_audit.py','count.cpp']},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['family_controls.json','family_banks.json']}}
    (P/'family_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
