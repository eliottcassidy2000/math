#!/usr/bin/env python3
"""Full exact producer replay; checks remain active under python -O."""
import collections,csv,hashlib,io,json,pathlib,subprocess,tempfile
P=pathlib.Path(__file__).resolve().parent

def check(test,msg):
    if not test:raise RuntimeError(msg)
def j(name):return json.loads((P/name).read_text())
def rows(name):return list(csv.DictReader((P/name).open(),delimiter='\t'))
def decode(code):
    r=[0]*9;b=35
    for k in range(1,9):
        for i in range(k):
            if code>>b&1:r[i]|=1<<k
            else:r[k]|=1<<i
            b-=1
    return r
def expand(m):return sum(7<<(3*i) for i in range(9) if m>>i&1)
def inflate(q):return [(1<<(3*i+(v+1)%3))|expand(q[i]) for i in range(9) for v in range(3)]
def aug(r,a,c):
    n=len(r);f=(1<<n)-1
    return [x|((a>>i&1)<<n)|((c>>i&1)<<(n+1)) for i,x in enumerate(r)]+[(f^a)|(1<<(n+1)),f^c]
def tournament(r):
    n=len(r)
    return all(0<=r[i]<1<<n and not(r[i]>>i&1) for i in range(n)) and all(((r[i]>>k)&1)+((r[k]>>i)&1)==1 for i in range(n) for k in range(i+1,n))
def sha(data):return hashlib.sha256(data).hexdigest()
def main():
    Q=[(0,1),(1,0),(1,2),(1,3),(3,1),(3,3),(3,5),(0,7),(0,18),(40,256),(4,3),(130,128)]
    battery=j('battery.json')
    check(battery==[dict(index=i,a=a,c=c,physical_a=expand(a),physical_c=expand(c)) for i,(a,c) in enumerate(Q)],'battery masks')
    check(len(set(Q+[(1<<i,0) for i in range(9)]))==20,'union size')
    census=rows('census.tsv');filt=rows('filter.tsv')
    check(len(census)==len(filt)==97144,'census size')
    check([int(x['node']) for x in census]==list(range(97144)),'census IDs')
    check(len({x['code'] for x in census})==97144,'distinct source codes')
    candidates=j('candidates.json');C=candidates['candidates'];derived=[];anchors=[];prime=0
    for a,b in zip(census,filt):
        check(all(a[k]==b[k] for k in ('node','code')),'filter provenance')
        r=decode(int(a['code']));check(tournament(r),'decoded tournament')
        m=int(b['proper_module'])
        if m:
            check(1<m<511 and m.bit_count()>1,'proper module size')
            check(all((r[v]&m) in (0,m) for v in range(9) if not(m>>v&1)),'literal nonprime witness')
            check(all(b[k]=='-' for k in ('p_mask','starts','ends')),'nonprime records')
            continue
        prime+=1;s=list(map(int,b['starts'].split(',')));e=list(map(int,b['ends'].split(',')))
        check(len(s)==len(e)==9 and sum(s)==sum(e),'endpoint sums')
        p=sum(((s[i]+e[i])%2)<<i for i in range(9));check(p==int(b['p_mask']),'endpoint parity')
        if not p:derived.append(dict(node=int(a['node']),code=int(a['code']),rows=r,starts=s,ends=e))
        else:
            i=(p&-p).bit_length()-1
            anchors.append(dict(node=int(a['node']),p_mask=p,anchor=i,coarse_a=1<<i,coarse_c=0,physical_a=expand(1<<i),physical_c=0,contrast_mod4=2))
    check(prime==69445 and len(anchors)==67097 and derived==C and len(C)==2348,'complete prime partition')
    check({k:candidates[k] for k in ['total_merged_nodes','prime_nodes','p_nonzero','p_zero']}==dict(total_merged_nodes=97144,prime_nodes=69445,p_nonzero=67097,p_zero=2348),'candidate metadata')
    check([{k:int(v) for k,v in x.items()} for x in rows('endpoint_witnesses.tsv')]==anchors,'endpoint physical witness retention')
    bank=rows('residues.tsv');W=j('physical_witnesses.json');first=j('search_mod8.json')
    check(len(bank)==len(W)==len(first)==len(C),'response record sizes')
    hist=collections.Counter();nonzero=0
    for q,r,w,old in zip(C,bank,W,first):
        check(q['node']==int(r['node'])==w['node']==old['node'],'response node order')
        values=[(int(r[f'S{i}']),int(r[f'T{i}'])) for i in range(12)]
        check(all(0<=s<8 and 0<=t<8 and s%2==t%2==1 for s,t in values),'physical response residues')
        check(all((s-t)%8 in (0,4) for s,t in values),'response gap shape')
        hits=[i for i,(s,t) in enumerate(values) if s!=t];check(bool(hits),'unseparated prime candidate')
        nonzero+=len(hits);i=hits[0];hist[i]+=1;b=battery[i];s,t=values[i];qop=[511^x^(1<<k) for k,x in enumerate(q['rows'])]
        expected=dict(node=q['node'],code=q['code'],query_index=i,coarse_a=b['a'],coarse_c=b['c'],physical_a=b['physical_a'],physical_c=b['physical_c'],source_residue=s,target_residue=t,augmented_source_rows=aug(inflate(q['rows']),b['physical_a'],b['physical_c']),augmented_target_rows=aug(inflate(qop),b['physical_a'],b['physical_c']))
        check(w==expected,'physical augmented graph or mask')
        check(tournament(w['augmented_source_rows']) and tournament(w['augmented_target_rows']),'physical augmented tournaments')
        check(old==dict(node=q['node'],query_index=i,source_residue=s,target_residue=t),'exploratory first witness')
    histogram=[hist[i] for i in range(12)]
    check(histogram==[1159,630,296,133,68,0,30,15,7,4,4,2],'witness histogram')
    controls=j('known_controls.json');controlqueries=[tuple(z['query']) for z in controls[:4]]
    check(all([tuple(z['query']) for z in controls[i:i+4]]==controlqueries for i in range(0,16,4)),'known controls ordering')
    inputbank=str(len(Q))+'\n'+'\n'.join(f'{a} {c}' for a,c in Q)+'\n'+str(len(C))+'\n'+'\n'.join(str(q['node'])+' '+' '.join(map(str,q['rows'])) for q in C)+'\n'
    check(inputbank==(P/'bank_input.txt').read_text(),'bank input provenance')
    with tempfile.TemporaryDirectory(prefix='c3-frontier-') as temp:
        temp=pathlib.Path(temp)
        for name in ('filter','bank'):
            subprocess.run(['g++','-O3','-std=c++17',str(P/(name+'.cpp')),'-o',str(temp/name)],check=True)
        actual=subprocess.run([str(temp/'filter')],input=(P/'census.tsv').read_bytes(),stdout=subprocess.PIPE,check=True).stdout
        check(actual==(P/'filter.tsv').read_bytes(),'clean complete filter replay')
        actual=subprocess.run([str(temp/'bank')],input=inputbank.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True).stdout
        check(actual==(P/'residues.tsv').read_bytes(),'clean full residue bank replay')
        controlinput='4\n'+'\n'.join(f'{a} {c}' for a,c in controlqueries)+'\n4\n'+'\n'.join(str(i)+' '+' '.join(map(str,controls[4*i]['rows'])) for i in range(4))+'\n'
        actual=subprocess.run([str(temp/'bank')],input=controlinput.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True).stdout.decode()
        rec=list(csv.DictReader(io.StringIO(actual),delimiter='\t'))
        for i,r in enumerate(rec):
            for k in range(4):
                check((int(r[f'S{k}']),int(r[f'T{k}']))==(controls[4*i+k]['source'],controls[4*i+k]['target']),'prior exact-profile control')
    names=['PROOF.md','census.tsv','filter.cpp','filter.tsv','bank.cpp','battery.json','candidates.json','residues.tsv','physical_witnesses.json','endpoint_witnesses.tsv','known_controls.json','prepare.py','verify.py','DEPENDENCIES.json']
    receipt=dict(status='PASS',scope='FINITE: every prime quotient of order nine with nine unchanged physical C3 interiors; source-relative <=27 consequence uses pinned reviewed inputs',merged_census=97144,prime_merged=prime,p_nonzero=len(anchors),p_zero=len(C),fixed_residual_queries=12,endpoint_union_queries=20,state_count=94208,complete_residue_count=2*12*len(C),nonzero_query_pairs=nonzero,first_witness_histogram=histogram,physical_witness_graphs=2*len(W),all_prime_endpoint_coordinates=18*prime,known_profile_control_queries=len(controls),all_named_terms='Actual complete physical queries are counted; no quadratic or cubic truncation. This is not a full named-polynomial bank at order 27.',ambient_limit='No new arbitrary-ambient lift is claimed. At <=27 the chosen core is the whole source.',sha256={name:sha((P/name).read_bytes()) for name in names})
    (P/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='sha256'},indent=2))
if __name__=='__main__':main()
