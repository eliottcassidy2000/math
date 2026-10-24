"""Second path enumeration and arithmetic audit; no native570 producer import."""
from common import *
from itertools import combinations
from collections import Counter
def main():
    C=json.loads((HERE/'native570.normal.json').read_text());N,owners,merge,X,Y=data();old=inventory();checks=0
    def test(v,msg):
        nonlocal checks
        checks+=1;need(v,msg)
    test([r['source'] for r in C['literal_sources']]==[1130,1673,1754,2678],'literal selected owners')
    used=set()
    for r in old:
        for key in ('before','after'):used.update(orbit(r[key]['mask']))
    for r in C['literal_sources']:
        s=r['source'];M=matrix(N[s]['canon_a']);F=dict(r['fixed_clock_owners']);D=[]
        for p in paths(M):
            if all(p[k]==v for k,v in F.items()):D.append(p)
        test(len(D)==r['raw_domain_size'],'DFS independently exhausts raw pinned domain')
        eligible=[]
        for p in D:
            w=word(M,p);t=merge[owners[FULL^w]]
            if t in (570,2234):eligible.append((p,w,t))
        test(sorted(eligible)==sorted((tuple(q['path']),q['word'],q['target']) for q in r['domain']),'no missing eligible lift')
        for p,w,t in eligible:
            O=set(orbit(w));test(len(O)==4 and not used&O,'old/new global resource injection');used.update(O)
            E={tuple(sorted(e)) for e in zip(p,p[1:])}
            K=[[0 if i==j else M[i][j] if tuple(sorted((i,j))) in E else M[j][i] for j in range(8)] for i in range(8)]
            test(K==next(q['target_matrix'] for q in r['domain'] if q['target']==t),'target matrix independently reconstructed')
        test(all(p[r['shared_anchor_clock']]==F[r['shared_anchor_clock']] for p,w,t in eligible),'shared anchor at the same actual local clock')
    ids=[r['source'] for r in C['literal_sources']];accept=[];base=[1754,2678]
    for pair in combinations(ids,2):
        x=(sum(X[s] for s in pair)-sum(X[s] for s in base))%128
        yy=[(sum(Y[j][s] for s in pair)-sum(Y[j][s] for s in base))%2 for j in range(4)]
        if x==0 and not any(yy):accept.append(pair)
    test(accept==[(1130,1673),(1754,2678)],'all six balanced states via independent subset arithmetic')
    test(sum(X[s] for s in ids[:2])==382 and sum(X[s] for s in ids[2:])==254,'integer debt before modular projection')
    for chi in [X]+[[64*y for y in row] for row in Y[:4]]:
        u=sum(chi[s] for s in ids[:2])-sum(chi[s] for s in ids[2:]);v=chi[570]-chi[2234]
        test(2*u%256==2*v%256==0,'rank-two physical-row identity annihilates full original generators')
    test(len(used)==256,'old28 plus fresh4 reserved alternatives')
    # Incoming audited40-bank is a different supply; authenticate disjointness.
    path='05-knowledge/results/compatible_charged_helpers_20261023/supply/selected_supply.json'
    raw=subprocess.check_output(['git','show','1d94c1d1d9:'+path],cwd=ROOT);bank=json.loads(raw)
    b=bank['banks'][bank['main_bank_index']]['source_ids'];helpers={r['source']:r for r in bank['helpers']}
    test(len(b)==40 and set(b).isdisjoint(ids) and set(b).isdisjoint(r['source'] for r in old),'different actual source owners')
    all_helpers=b+sorted(set(helpers)-set(b))
    test(len(all_helpers)==53 and set(all_helpers).isdisjoint(ids),'all53 supplied helpers remain different source slots')
    for j,s in enumerate(all_helpers):
        r=helpers[s];M=matrix(N[s]['canon_a'])
        for key in ('before','after'):
            q=r['selected'][key];p=q['path'];w=word(M,p)
            test(all(M[i][j] for i,j in zip(p,p[1:])) and w==q['mask'],'incoming record literal path')
            test(merge[owners[w]]==s and merge[owners[FULL^w]]==q['target'],'incoming original source and target')
            O=set(orbit(w));test(len(O)==4 and not used&O,'incoming40-bank and new570 return reserve distinct resources');used.update(O)
        if j==39:test(len(used)==576,'144 distinct alternatives across72 slots with main40 bank')
    test(len(used)==680,'170 distinct alternatives across85 slots with full53 bank')
    return {'status':'PASS same-agent independent enumeration, current identity and incoming-bank disjointness',
            'checks':checks,'two_state_V_sets':accept,'integer_source_X_sums':[382,254],
            'incoming_bank_source_overlap':[],'old28_plus_new4_resources':64,'with_incoming40_reserved_resources':144,
            'with_incoming40_arrival_resources':72,'incoming_revision':'1d94c1d1d9',
            'with_incoming53_reserved_resources':170,'with_incoming53_arrival_resources':85,
            'incoming_supply_sha256':hashlib.sha256(raw).hexdigest(),'original_inputs':INPUTS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
