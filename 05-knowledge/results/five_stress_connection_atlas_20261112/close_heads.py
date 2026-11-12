"""Construct actual below-ray heads using factors and the inherited compiler."""
from pathlib import Path
import importlib.util,json
from sympy import factorint
from tournament_tools import hp
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('c',HERE.parent/'session251_creative_completion_20261111/ternary/compiler.py')
c=importlib.util.module_from_spec(sp);sp.loader.exec_module(c)

def need(b,msg):
    if not b:raise RuntimeError(msg)
def apex_count(word):
    n=len(word)
    if not n:return 1
    answer=int(word[0]==1)+int(word[-1]==0)
    answer+=sum(word[i]==0 and word[i+1]==1 for i in range(n-1))
    answer+=sum(2**(i-j-1) for i in range(n) for j in range(i) if word[j]==1 and word[i]==0)
    return answer
def join(blocks):
    n=sum(map(len,blocks));out=[];offset=0
    for B in blocks:
        after=((1<<n)-1)^((1<<(offset+len(B)))-1)
        out.extend((r<<offset)|after for r in B);offset+=len(B)
    return out

def run():
    inventory=json.loads((HERE.parent/'session173_modular_clocks_20261103/counts527/inventory_375.json').read_text())
    bank={1:[0]}
    for ns,items in inventory['orders'].items():
        for code,H,*_ in items:
            if H not in bank:bank[H]=c.tournament(int(ns),int(code))
    sigs=json.loads((HERE.parent/'chanpak_five_core_bridges_20261110/cover/signatures.json').read_text())
    values=[(393,4),(1359,1),(44825,3),(236287,0)]
    targets=[1+odd*(1<<s) for odd,j in values for s in range(1,j+1)]
    primes=sorted({int(p) for N in targets for p in factorint(N)})
    witnesses={}
    for p in primes:
        if p in bank:
            B=bank[p];need(hp(B)==p,'inventory component');witnesses[p]={'method':'strong inventory, direct DP replay','rows':B,'H':p};continue
        address=next((2*q for q in range(p) if c.g(2*q)==p),None)
        if address is not None:
            B=c.construct([],0,address,0);word=list(map(int,bin(address)[2:]))
            need(apex_count(word)==p,'independent endpoint formula')
            if len(B)<=15:need(hp(B)==p,'small direct apex DP')
            witnesses[p]={'method':'one apex, independent two-endpoint formula','address':address,'word':word,'rows':B,'H':p};continue
        found=None
        for row in sigs:
            for H,a,b in c.generate(row['signature'],p):
                if H==p:found=(row,a,b);break
            if found:break
        need(found is not None,('target inverse failed',p))
        row,a,b=found;B=c.construct(row['core'],row['cut'],a,b)
        need(hp(B)==p,'literal compiler child DP')
        witnesses[p]={'method':'inherited five-core compiler, direct DP replay','core':row['core'],'cut':row['cut'],'A':a,'B':b,'rows':B,'H':p}
    heads=[]
    for N in targets:
        factors=[int(p) for p,e in factorint(N).items() for _ in range(e)]
        B=join([witnesses[p]['rows'] for p in factors])
        product=1
        for p in factors:product*=p
        need(product==N,'SCC product')
        if len(B)<=15:need(hp(B)==N,'small literal whole graph DP')
        heads.append({'H':N,'factors':factors,'rows':B,'order':len(B)})
    # Broad independent control for the one-apex formula.
    for a in range(1,4096):
        need(apex_count(list(map(int,bin(a)[2:])))==c.g(a),'apex formula control')
    out={'scope':'Eight actual finite heads complete four coefficient rays, using the four original C5-core seeds through order12. No new bounded odd-coverage range is claimed.','rays':[{'c':odd,'tail_from_s':j+1,'head_s':list(range(1,j+1))} for odd,j in values],'heads':heads,'prime_component_witnesses':witnesses,'apex_formula_controls':4095}
    dest=HERE/('heads.optimized.json' if not __debug__ else 'heads.json')
    dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'heads':[(h['H'],h['order'],h['factors']) for h in heads],'component_methods':[(p,w['method'],len(w['rows'])) for p,w in witnesses.items()]}))

if __name__=='__main__':run()
