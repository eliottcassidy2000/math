"""Finite exact marked block-swap probe on the inherited n<=7 bank."""
from collections import defaultdict
from fractions import Fraction
from itertools import permutations
from math import comb,gcd
import hashlib
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent
BANK='05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/extended_census.json'

def ends(t):
    n=len(t);z=[[0]*n for _ in range(1<<n)]
    for i in range(n):z[1<<i][i]=1
    for mask,row in enumerate(z):
        for i,v in enumerate(row):
            if not v:continue
            q=t[i]&~mask
            while q:
                b=q&-q;q-=b;z[mask|b][b.bit_length()-1]+=v
    return z

def converse(t):return tuple(((1<<len(t))-1)^row^(1<<i) for i,row in enumerate(t))

def ordinal(*parts):
    n=sum(map(len,parts));out=[];offset=0
    for part in parts:
        later=((1<<n)-1)^((1<<(offset+len(part)))-1)
        out.extend((row<<offset)|later for row in part);offset+=len(part)
    return tuple(out)

def transitive(n):return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def delete(t,v):
    ids=[i for i in range(len(t)) if i!=v]
    return tuple(sum(1<<j for j,w in enumerate(ids) if t[u]>>w&1) for u in ids)

def cyclic(t):
    n=len(t)
    return tuple(row|(1<<n) for row in t)+(1<<(n+1),(1<<n)-1)

def root_state(t,v):
    e=ends(t);h=[sum(x) for x in e];h[0]=1;full=len(h)-1
    return (h[-1],ends(converse(t))[-1][v],h[full^(1<<v)],
            sum(e[sub][v]*h[full^sub] for sub in range(full+1)))

def two_slots(t):
    h=[sum(x) for x in ends(t)];h[0]=1;full=len(h)-1
    return h[-1],sum(h[sub]*h[full^sub] for sub in range(full+1))

def data(t):
    e=ends(t);es=ends(converse(t));h=[sum(x) for x in e];h[0]=1
    f2=[]
    for mask in range(len(h)):
        value=0;sub=mask
        while True:
            value+=h[sub]*h[mask^sub]
            if not sub:break
            sub=(sub-1)&mask
        f2.append(value)
    full=len(h)-1
    f3=sum(h[sub]*f2[full^sub] for sub in range(full+1))
    roots=[]
    for v in range(len(t)):
        e2=sum(e[sub][v]*h[full^sub] for sub in range(full+1))
        s2=sum(es[sub][v]*h[full^sub] for sub in range(full+1))
        e3=sum(e[sub][v]*f2[full^sub] for sub in range(full+1))
        roots.append(dict(root=v,E1=e[full][v],S1=es[full][v],q=h[full^(1<<v)],
                          E2=e2,S2=s2,E3=e3,
                          alpha=str(Fraction(2*e2+s2,2*e3))))
    return dict(rows=t,n=len(t),H=h[-1],K=f2[-1],L=f3,roots=roots)

def controls():
    empty=dict(rows=(),H=1,K=1,L=1)
    c3=(2,4,1);sides=[empty,data((0,)),data(c3)];kernels=[data((0,)),data(c3)]
    count=0
    for left in sides:
        for core in kernels:
            for right in sides:
                t=ordinal(left['rows'],core['rows'],right['rows']);d=data(t)
                for r in core['roots']:
                    v=len(left['rows'])+r['root'];z=d['roots'][v]
                    expected=(left['H']*r['E1']*(not right['rows']),
                              left['K']*r['E2']*right['H'],
                              left['L']*r['E3']*right['K'])
                    assert (z['E1'],z['E2'],z['E3'])==expected
                    assert z['S2']==left['H']*r['S2']*right['K']
                    dh,dk=two_slots(delete(t,v))
                    expected_state=(d['H']+d['K'],z['S2'],dh+dk-(len(t)==1),
                                    2*z['E2']+2*z['E3']-z['E1'])
                    assert root_state(cyclic(t),v)==expected_state
                    count+=1
    h17=[root_state(cyclic(transitive(4)),j) for j in range(4)]
    assert h17==[(17,8,9,18),(17,4,9,28),(17,2,9,44),(17,1,9,69)]
    return dict(marked_factor_and_literal_wrapper_root_controls=count,H17_ranks=h17)

def named_modules(a,k,b,reverse=False):
    parts=[a,k,b];offset=[0,len(a),len(a)+len(k)];out=[0]*(len(a)+len(k)+len(b))
    for part,start in zip(parts,offset):
        for i,row in enumerate(part):out[start+i]=row<<start
    order=[2,1,0] if reverse else [0,1,2]
    for i,x in enumerate(order):
        for y in order[i+1:]:
            mask=((1<<len(parts[y]))-1)<<offset[y]
            for v in range(offset[x],offset[x]+len(parts[x])):out[v]|=mask
    return tuple(out)

def literal_profile3(t):
    out=[0,0,0];n=len(t)
    for path in permutations(range(n)):
        backwards=sum(not (t[a]>>b&1) for a,b in zip(path,path[1:]))
        for x in (1,2,3):
            top=x+n-1-backwards
            if top>=n:out[x-1]+=comb(top,n)
    return out

def pair_receipt(hit):
    a,b,k,r=hit['A'],hit['B'],hit['K'],hit['root']
    assert literal_profile3(a['rows'])==[a['H'],a['K'],a['L']]
    assert literal_profile3(b['rows'])==[b['H'],b['K'],b['L']]
    dh,dk=two_slots(delete(k['rows'],r['root']))
    h=a['H']*k['H']*b['H']+a['K']*k['K']*b['K']
    q=a['H']*dh*b['H']+a['K']*dk*b['K'];states=[];returns=[];rows=[]
    for rev,(left,right) in enumerate(((a,b),(b,a))):
        s=left['H']*r['S2']*right['K']
        f=2*left['K']*r['E2']*right['H']+2*left['L']*r['E3']*right['K']
        states.append((h,s,q,f));returns.append((h+f-s,q,q,2*f-q))
        rows.append(cyclic(named_modules(a['rows'],k['rows'],b['rows'],bool(rev))))
    assert returns[0][0]==returns[1][0]
    diff=tuple(x-y for x,y in zip(*states));assert diff[0]==diff[2]==0 and diff[1]==diff[3]!=0
    return dict(A_H=a['H'],B_H=b['H'],order=hit['order'],alpha=hit['alpha'],
                root=a['n']+r['root'],named_source_rows=rows,states=states,R_states=returns,
                difference=diff,root_gcds=[gcd(*z) for z in states],
                R_root_gcds=[gcd(*z) for z in returns],
                direct_large_graph_count=False,
                certificate='proved marked product/wrapper laws; independent permutation/cut triples for A and B')

def main():
    revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    raw=subprocess.check_output(['git','show',revision+':'+BANK])
    records=[r for r in json.loads(raw)['records'] if r['n']<=7 and r['strong']]
    blocks={};rooted={};root_count=0;all_roots=[]
    extra=[(56,117,129,22,100,76,141,59),(56,53,129,22,100,76,143,59),
           (56,117,137,18,100,76,141,59),(56,53,137,18,100,76,143,59),
           (212,137,186,97,170,195,22,72)]
    for row in [(0,)]+[tuple(r['adjacency_outmasks']) for r in records]+extra:
        for t in {tuple(row),converse(row)}:
            d=data(t);key=(d['H'],d['K'],d['L'])
            blocks.setdefault(key,d)
            for root in d['roots']:
                root_count+=1;alpha=Fraction(root['alpha'])
                all_roots.append((d,root))
                old=rooted.get(alpha)
                if old is None or d['n']<old['block']['n']:
                    rooted[alpha]=dict(block=d,marked=root)
    bs=list(blocks.values());slopes=[];hits=[]
    for i,a in enumerate(bs):
        for b in bs[i+1:]:
            den=a['H']*b['K']-b['H']*a['K']
            if not den:continue
            alpha=Fraction(a['L']*b['K']-b['L']*a['K'],den)
            if not 0<alpha<=Fraction(3,2):continue
            entry=dict(A=a,B=b,alpha=str(alpha))
            slopes.append(entry)
    strong_root_hits=sum(Fraction(x['alpha']) in rooted for x in slopes)
    targets={Fraction(x['alpha']) for x in slopes}
    padded_matches=0
    for core,oldroot in all_roots:
        for a in range(7):
            for b in range(7):
                e2=2**a*oldroot['E2'];s2=2**b*oldroot['S2'];e3=3**a*2**b*oldroot['E3']
                alpha=Fraction(2*e2+s2,2*e3)
                if alpha not in targets:continue
                padded_matches+=1;n=core['n']+a+b;old=rooted.get(alpha)
                if old is not None and old['block']['n']<=n:continue
                root=dict(root=a+oldroot['root'],E1=oldroot['E1'] if b==0 else 0,
                          E2=e2,S2=s2,E3=e3,alpha=str(alpha))
                block=dict(rows=ordinal(transitive(a),core['rows'],transitive(b)),n=n,
                           H=core['H'],K=2**(a+b)*core['K'],L=3**(a+b)*core['L'],
                           roots=[root],construction='transitive prefix/core/suffix',
                           prefix=a,suffix=b,core=core)
                rooted[alpha]=dict(block=block,marked=root)
    empty=dict(rows=(),n=0,H=1,K=1,L=1,roots=[])
    constructed=0
    for left in [empty]+bs:
        for right in [empty]+bs:
            e2=left['K']*right['H'];s2=left['H']*right['K'];e3=left['L']*right['K']
            alpha=Fraction(2*e2+s2,2*e3)
            if alpha not in targets:continue
            constructed+=1
            n=left['n']+right['n']+1;old=rooted.get(alpha)
            if old is not None and old['block']['n']<=n:continue
            root=dict(root=left['n'],E1=left['H'] if not right['n'] else 0,
                      E2=e2,S2=s2,E3=e3,alpha=str(alpha))
            block=dict(rows=ordinal(left['rows'],(0,),right['rows']),n=n,
                       H=left['H']*right['H'],K=2*left['K']*right['K'],
                       L=3*left['L']*right['L'],roots=[root],
                       construction='L -> singleton root -> R',left=left,right=right)
            rooted[alpha]=dict(block=block,marked=root)
    for entry in slopes:
        alpha=Fraction(entry['alpha'])
        if alpha not in rooted:continue
        k=rooted[alpha];a=entry['A'];b=entry['B']
        hits.append(dict(**entry,K=k['block'],root=k['marked'],
                         order=a['n']+b['n']+k['block']['n']+2))
    hits.sort(key=lambda x:x['order'])
    checks=controls();pairs=[pair_receipt(hit) for hit in hits]
    out=dict(status='FINITE-EXACT; named strong blocks through7 plus singleton and five inherited eight-vertex repair blocks, both converse root gauges',
             revision=revision,bank=BANK,strong_converse_records=len(records),
             additional_literal_eight_vertex_blocks=extra,
             distinct_unmarked_triples=len(bs),root_gauges=root_count,
             distinct_root_slopes=len(rooted),admissible_positive_secants=len(slopes),
             strong_root_hits=strong_root_hits,
             padded_root_matches=padded_matches,padding_bound_each_flank=6,
             constructed_root_matches=constructed,
             exact_compensating_swaps=len(hits),smallest_hits=hits[:10],
             admissible_secants=slopes,small_literal_controls=checks,
             compensated_literal_pairs=pairs,
             source_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                            for name in ('PROOF.md','mixed_probe.py')})
    (HERE/'mixed_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('strong_converse_records','distinct_unmarked_triples',
          'root_gauges','distinct_root_slopes','admissible_positive_secants','strong_root_hits',
          'padded_root_matches','constructed_root_matches','exact_compensating_swaps')},indent=2))
    if pairs:print(json.dumps({k:pairs[0][k] for k in ('order','root','states','R_states','root_gcds','R_root_gcds')},indent=2))

if __name__=='__main__':main()
