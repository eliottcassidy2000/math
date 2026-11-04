#!/usr/bin/env python3
"""Independent composition of all finite child-type coverage certificates."""
import pathlib,json,hashlib,itertools,sys
HERE=pathlib.Path(__file__).resolve().parent
P=pathlib.Path(sys.argv[1]).resolve()if len(sys.argv)>1 else HERE.parent/'session189-profiles'
E=pathlib.Path(sys.argv[2]).resolve()if len(sys.argv)>2 else HERE.parent/'session189-prime7'
def need(v,msg):
    if not v:raise RuntimeError(msg)
def sub(a):
    return[x for x in range(128)if x&a==x]
def load(p):return json.loads(p.read_text())
def inflate(Q,ns,a,c):
    off=[];n=0
    for s in ns:off.append(n);n+=s
    r=[0]*(n+2);r[n]=1<<(n+1)
    for i,s in enumerate(ns):
        for u in range(s):
            v=off[i]+u
            for j,t in enumerate(ns):
                if Q[i]>>j&1:r[v]|=((1<<t)-1)<<off[j]
            for u2 in range(u+1,s):r[v]|=1<<(off[i]+u2)
            for mask,port in((a,n),(c,n+1)):
                if mask>>i&1:r[v]|=1<<port
                else:r[port]|=1<<v
    return r,off
def relabel(r,pi):
    ans=[0]*len(r)
    for u in range(len(r)):
        ans[pi[u]]=sum(1<<pi[v]for v in range(len(r))if r[u]>>v&1)
    return ans
def main():
    reps=load(P/'quotients.json');anchors=load(P/'anchor_exploration.json');res4=load(P/'residual_mod4.json');aff=load(P/'affine_banks.json');safe=load(P/'ambient_analysis.json');wit=load(P/'type_witnesses.json')
    need(len(reps)==197 and len({tuple(x['rows'])for x in reps})==197,'Representative uniqueness')
    nonflat=[z for z in reps if any(z['p'])];flat=[z for z in reps if not any(z['p'])]
    need(len(nonflat)==159 and {tuple(z['rows'])for z in flat}=={tuple(z['rows'])for z in anchors},'Endpoint/anchor partition')
    anchor_weight=0;blind=[]
    for z in anchors:
        need([t['V']for t in z['patterns']]==list(range(128)),'Complete singleton-position domain')
        for t in z['patterns']:
            if 2 in t['first'] or 4 in t['second']:anchor_weight+=4**t['V'].bit_count()
            else:blind.append((tuple(z['rows']),t['V']))
    need(len(blind)==10,'Anchor blind pattern count')
    singleton_weight=0
    for z in res4:
        need((tuple(z['Q']),z['V'])in blind and z['V']!=127,'Unexpected singleton residual')
        need(any(x%4==2 for x in z['contrast']),'Singleton residue witness')
        singleton_weight+=4**z['V'].bit_count()
    residue_rows={tuple(q)for q,v in blind if v==127}
    need(len(res4)==6 and len({(tuple(z['Q']),z['V'])for z in res4})==6 and len(residue_rows)==4,'Residual partition cardinalities')
    need({tuple(z['Q'])for z in aff}=={tuple(z['Q'])for z in safe}=={tuple(z['Q'])for z in wit}==residue_rows,'Full four-quotient correspondence')
    bysafe={tuple(z['Q']):z for z in safe};bywit={tuple(z['Q']):z for z in wit}
    masks=[a for a in range(128)if a.bit_count()<=2];ix={x:i for i,x in enumerate(masks)}
    basis=[(a,c)for a in masks for c in masks]
    contributors=[[29*ix[x]+ix[y]for x in sub(a)for y in sub(c)]for a,c in basis]
    local_count=odd_count=0;uncovered=[];row_details=[]
    for z in aff:
        q=tuple(z['Q']);S=bysafe[q];W=bywit[q];src,dst=z['full_run_coefficients']
        K0=[a-b for a,b in zip(src['K0'],dst['K0'])]
        K2=[[a-b for a,b in zip(x,y)]for x,y in zip(src['K2'],dst['K2'])]
        K3=[[a-b for a,b in zip(x,y)]for x,y in zip(src['K3'],dst['K3'])]
        K22=[[a-b for a,b in zip(x,y)]for x,y in zip(src['K22'],dst['K22'])]
        coef=[[K0[j]+2*sum(v[j]for v in K2)+4*sum(v[j]for v in K22)for j in range(841)]]+[[4*x for x in v]for v in K2]+[[2*x for x in v]for v in K3]
        val=[[sum(row[k]for k in subix)%8 for subix in contributors]for row in coef]
        need([r[0]for r in W['records']]==list(range(16384)),'Complete 14-bit profile domain')
        n_odd=0;n_ok=0
        for code,odd,j,a,c,gap in W['records']:
            need(odd==((code>>7).bit_count()%2),'Wrong actual total-order parity')
            if j<0:
                need(a is None and c is None and gap==0,'Malformed unresolved record')
                need(all((val[0][k]+sum(val[t+1][k]for t in range(14)if code>>t&1))%8==0 for k in range(841)),'False local blindness')
                uncovered.append((q,code));continue
            need(basis[j]==(a,c),'Query index/mask disagreement')
            actual=(val[0][j]+sum(val[t+1][j]for t in range(14)if code>>t&1))%8
            need(actual==gap==4,'Selected actual physical query has wrong gap')
            if odd:need(S['cycle_parity_contrast'][j]==0,'Unsafe odd-total witness');n_odd+=1
            n_ok+=1
        local_count+=n_ok;odd_count+=n_odd
        row_details.append({'Q':q,'covered_mod8_types':n_ok,'odd_types_with_checked_safe_actual_query':n_odd})
    carousel=(14,28,56,112,97,67,7)
    need(uncovered==[(carousel,127)],'The two packets do not meet at exactly all-E carousel')
    rows=load(E/'rows.json');queries=load(E/'even_queries.json');raw=[list(map(int,line.split()))for line in(E/'even_bank.txt').read_text().splitlines()]
    ri=rows.index(list(carousel));mod=[16]*5+[2]*105;cert=[]
    for query,constant in(([0,3],0),([4,68],8)):
        qi=queries.index(query);at=2*(ri*len(queries)+qi);delta=[(x-y)%m for x,y,m in zip(raw[at],raw[at+1],mod)]
        target=[constant,0,0,0,0]+[int(k==2 and i in(3,5)and h==1)for k in(2,3,5)for i in range(7)for h in range(1,6)]
        need(delta==target,'Incomplete all-cover mod16 identity')
        cert.append({'query':query,'all_110_contrast_channels':delta})
    for e3,e5 in itertools.product((0,1),repeat=2):need(sorted([(8*(e3+e5))%16,(8*(1+e3+e5))%16])==[0,8],'Complementary all-profile branch')
    total=159*5**7+anchor_weight+singleton_weight+local_count+1
    need(total==197*5**7 and local_count==65535 and odd_count==32768,'Incomplete formal five-type cover')
    # Physical transport uses ONE common relabeling of both quotient sides,
    # every complete child set and both masks. Ports stay fixed.
    permutations=[[(i+2)%7 for i in range(7)],list(reversed(range(7))),[3,0,6,2,5,1,4]]
    transports=0
    for rec in reps:
        for pi in permutations:
            ns=[1,2,1,1,2,1,1];newns=[0]*7
            for i in range(7):newns[pi[i]]=ns[i]
            for a,c in((1,2),(85,42)):
                aa=sum(1<<pi[i]for i in range(7)if a>>i&1);cc=sum(1<<pi[i]for i in range(7)if c>>i&1)
                for side in range(2):
                    Q=rec['rows']if not side else[127^(1<<i)^x for i,x in enumerate(rec['rows'])]
                    original,off=inflate(Q,ns,a,c);transportedQ=relabel(Q,pi);target,noff=inflate(transportedQ,newns,aa,cc)
                    mp=[0]*len(original)
                    for i in range(7):
                        for u in range(ns[i]):mp[off[i]+u]=noff[pi[i]]+u
                    mp[-2:]=[len(original)-2,len(original)-1]
                    need(relabel(original,mp)==target,'Physical child/mask transport failed')
                    transports+=1
    inputs=[P/name for name in('PROOF.md','quotients.json','anchor_exploration.json','residual_mod4.json','affine_banks.json','ambient_analysis.json','type_witnesses.json','normal.json','optimized.json')]+[E/name for name in('PROOF.md','rows.json','even_queries.json','even_bank.txt','normal.json','optimized.json')]
    result={'status':'PASS_ALL_CHILD_DOMAIN_AND_PHYSICAL_TRANSPORT','formal_mod8_type_assignments':total,'endpoint_branch_types':159*5**7,'anchor_branch_types':anchor_weight,'singleton_residual_types':singleton_weight,'nontrivial_mod8_types':local_count,'mod16_residual_types':1,'odd_total_witnesses_independently_checked':odd_count,'row_details':row_details,'complete_mod16_carousel_channels':cert,'physical_graph_transport_cases':transports,'source_pins':{str(p.relative_to(p.parent.parent)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs}}
    (HERE/'coverage_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: all 15390625 formal five-type assignments partitioned; 32768 odd-total actual witnesses safe; complete mod16 carousel residual; 2364 common physical label transports.')
if __name__=='__main__':main()
