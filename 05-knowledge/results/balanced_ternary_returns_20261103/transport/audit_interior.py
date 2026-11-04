#!/usr/bin/env python3
"""Independent 24-state audit of the interior law and complete H721 menu."""
from pathlib import Path
from itertools import combinations
import hashlib
import importlib.util
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result

def need(ok,message):
    if not ok:
        raise RuntimeError(message)

WORD_PATH=ROOT/'05-knowledge/results/orbit_spectrum_session_20261103/words/check.py'
W=load('independent_24state',WORD_PATH)
T=load('literal_transport_audit',HERE/'check.py')

def all_markings(rows):
    result=[]
    for p,q in combinations(range(len(rows)),2):
        if rows[q] >> p & 1:
            p,q=q,p
        rest=[v for v in range(len(rows)) if v not in (p,q)]
        scores={u:sum(rows[u] >> v & 1 for v in rest) for u in rest}
        rest.sort(key=lambda u:-scores[u])
        if all(rows[u] >> v & 1 for i,u in enumerate(rest) for v in rest[i+1:]):
            word=''.join(str((rows[p] >> v & 1)+2*(rows[q] >> v & 1)) for v in rest)
            result.append(((p,q),tuple(rest),word))
    return result

def superword_menu(old,extra):
    """Distinct superwords using a greedy subsequence automaton, not insertion addresses."""
    length=len(old)+extra
    counts={719:0,721:0}
    total=0
    witness=None
    def visit(prefix,index,state):
        nonlocal total,witness
        remaining=length-len(prefix)
        if remaining < len(old)-index:
            return
        if remaining==0:
            if index!=len(old) or prefix[-1]=='3':
                return
            total+=1
            value=W.dot(W.OUTPUT,state)
            if value in counts:
                counts[value]+=1
                if witness is None:
                    witness=prefix
            return
        for letter in '0123':
            if not prefix and letter=='0':
                continue
            new_index=index+int(index<len(old) and old[index]==letter)
            if remaining-1 < len(old)-new_index:
                continue
            visit(prefix+letter,new_index,W.step(state,int(letter)))
    visit('',0,W.INITIAL)
    return {'extra':extra,'distinct_strong_superwords':total,'hits':counts,'first_positive':witness}

def graded_bridge(P,root,k):
    # Names 0..n-1 are P; a=n, prefix n+1..n+k, new root n+k+1.
    n=len(P)
    a=n
    out=list(P)+[(1 << n)-1]
    for i in range(k):
        v=len(out)
        out.append((1 << v)-1)
        # Reverse the prefix's internal order so it is ascending in names.
        for j in range(n+1,v):
            out[v] ^= 1 << j
            out[j] |= 1 << v
    R=len(out)
    for v in range(R):
        if v in (a,root):
            out[v] |= 1 << R
    out.append(((1 << R)-1) ^ (1 << a) ^ (1 << root))
    return tuple(out),R

def main():
    # Independent state dimension and recurrence; producer uses36.
    certificates=[]
    for prefix in ('2','202'):
        v=W.state(prefix+'0')
        v1=W.step(v,0)
        v2=W.step(v1,0)
        need(all(c-3*b+2*a==0 for a,b,c in zip(v,v1,v2)),'all-state zero recurrence')
        residual_controls=0
        for k in (1,2):
            v=W.state(prefix+'0'*k+'3')
            v1=W.step(v,3)
            v2=W.step(v1,3)
            residual=tuple(c-3*b+2*a for a,b,c in zip(v,v1,v2))
            terminal=W.suffix_row('2')
            for _ in range(24):
                need(W.dot(terminal,residual)==0,'24 exact Cayley-Hamilton observations')
                residual=W.step(residual,3)
                residual_controls+=1
        certificates.append({'prefix':prefix,'zero_state_identity':True,'residual_zero_outputs':residual_controls})
    formula_controls=0
    for m in range(2,11):
        for r in range(1,9):
            a=W.hamilton('202'+'0'*(m-1)+'3'*r+'2')
            b=W.hamilton('2'+'0'*(m-1)+'3'*r+'2')
            need(a==5*(2**r+1)*2**m-1,'balanced H formula')
            need(b==(2**r+1)*(2**m+2)-1,'base H formula')
            formula_controls+=2
    chain=[]
    for word,expected in [('202032',(59,17,17,59)),
                          ('20203332',(179,53,53,179)),
                          ('202030332',(541,147,147,541))]:
        rows=T.graph(word)
        p=T.profile(rows,2)
        need(p==expected,'independent literal named balanced chain')
        chain.append({'word':word,'rows':rows,'root_profile':p})
    rows=T.graph('20200032')
    need(T.hamilton(rows)==239,'m4 literal source H')
    triangles=((1,2,3),(0,8,9))
    for a,b,c in triangles:
        need(rows[a]>>b&1 and rows[b]>>c&1 and rows[c]>>a&1,'two disjoint directed triangles')
    marked=all_markings(rows)
    need(len(marked)==3,'complete old feedback pair count')
    menus=[]
    for pair,chain_order,word in marked:
        menu=[superword_menu(word,e) for e in range(1,4)]
        need(not any(row['hits'][721] for row in menu),'all strong superwords miss721')
        menus.append({'physical_pair':pair,'old_chain':chain_order,'word':word,'menu':menu})
    need(any(row['hits'][719] for menu in menus for row in menu['menu']),'positive719 control')
    graded=[]
    for P,port in [((2,4,1),0),(T.graph('202032'),2)]:
        b,t,q,F=T.profile(P,port)
        need((q,F)==(t,b),'actual balanced input to graded bridge')
        for k in (0,1,2):
            U,r=graded_bridge(P,port,k)
            h=T.count_dp(U)
            if k==0:
                need(h==2*b-t,'empty-prefix native boundary')
            else:
                p=T.profile(U,r)
                need(p==(3*b*2**k-t,b,b,3*b*2**k-t+b+t),'literal graded bridge profile')
            graded.append({'input_pair':[b,t],'k':k,'rows':U,'root':r,'H':h})
    for k in range(1,7):
        for r in range(1,7):
            u=2**r+1
            need(W.hamilton('3'+'0'*k+'3'*r+'1')==3*u*2**k-1,'independent graded two-run formula')
    sources=(WORD_PATH,HERE/'check.py')
    print(json.dumps({'status':'PASS','state_dimension':24,'recurrence_certificates':certificates,
         'formula_controls':formula_controls,'literal_named_chain':chain,
         'm4_triangles':triangles,'complete_m4_superword_menus':menus,
         'graded_bridge_controls':graded,'graded_word_controls':36,
         'scope':'All strong final two-apex extensions of actual C4 at H721 via proved SCC/order reduction; no FVS>=3 exclusion.',
         'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
