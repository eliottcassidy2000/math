"""Small exact certificates for the large targeted tests; no large search."""
from pathlib import Path
from math import gcd
import json,sys
sys.path.insert(0,str(Path(__file__).parent/'dependencies'))
import core

def main():
    r=Path(__file__).parent
    cert=json.loads((r/'PRIME-CERTIFICATES.json').read_text());proved=set();modular_checks=0
    for p in sorted(map(int,cert)):
        c=cert[str(p)];core.need(c['prime']==p,'prime certificate label')
        if p==2:
            core.need(c['factorization']=={} and c['base'] is None,'prime two seed');proved.add(p);continue
        factors={int(q):e for q,e in c['factorization'].items()};product=1
        for q,e in factors.items():
            core.need(q in proved and e>=1,'certified prime factor')
            product*=q**e
        core.need(product==p-1,'complete p-minus-one factorization')
        a=c['base'];core.need(1<a<p and pow(a,p-1,p)==1,'Lucas Fermat condition');modular_checks+=1
        for q in factors:
            core.need(gcd(pow(a,(p-1)//q,p)-1,p)==1,'Lucas full-order condition');modular_checks+=1
        proved.add(p)
    menu={x['word'] for x in json.loads((r/'dependencies/MENUS.json').read_text())['5']}
    targets=json.loads((r/'TARGETED-WITNESSES.json').read_text());out=[]
    for t in targets:
        n=t['target'];core.need(t['prime_triple']==[n,3*n-2,5*n-4] and all(p in proved for p in t['prime_triple']),
                              'exact hostile prime triple')
        w,u,v=t['word'],t['left_index'],t['right_index']
        core.need(w in menu and u>=2 and v>=2,'actual twelve-template indices')
        core.need(core.K(core.tau(u),core.state(w,core.tau(v)))==n,'exact matrix target')
        for ww,uu,vv in t['normalization_trace']:
            core.need(core.K(core.tau(uu),core.state(ww,core.tau(vv)))==n,'physical trace readout')
        word=core.graph_word(w,u,v);s=core.state(word[:-1]);rows=core.graph(word)
        core.need(word==t['graph_word'] and list(s)==t['source_profile'] and core.F(s)==n,'full graph compiler profile')
        core.need(len(rows)==t['vertices']<=n.bit_length()-1+6,'constructed adjacency order')
        out.append({'target':n,'template':w,'indices':[u,v],'vertices':len(rows),
                    'two_rank_one_templates':'excluded at all heights by the prime triple'})
    return {'status':'PASS','Lucas_certified_primes':len(proved),'modular_checks':modular_checks,
            'target_certificates':out,'literal_large_path_enumeration':False,'universal_occupancy':'OPEN'}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
