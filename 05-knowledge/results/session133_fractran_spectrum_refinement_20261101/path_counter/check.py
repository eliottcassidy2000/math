"""Portable producer checks. All checks survive Python -O; see independent audit too."""
import hashlib
import itertools
import json
import math
from pathlib import Path
from compiler import Compiler, fractions, run_register, run_fractran_exponents, run_fractran_literal
from hamilton_source import hamilton
from square_shell_source import square_shell
import one_apex as oa
from fibres import f,G,word,reflected_word,inverse,bounded_fibres,summary
ROOT=Path(__file__).resolve().parent

def require(c,msg):
    if not c: raise RuntimeError(msg)

def rows_from_word(w):
    n=len(w)+1;rows=[0]*n
    for i,b in enumerate(w,1):
        if b: rows[0]|=1<<i
        else: rows[i]|=1
        for j in range(i+1,n): rows[i]|=1<<j
    return tuple(rows)

def rows_from_bits(n,b):
    rows=[0]*n
    for j in range(n):
        for i in range(j):
            if b>>(j*(j-1)//2+i)&1: rows[i]|=1<<j
            else: rows[j]|=1<<i
    return tuple(rows)

def check_rows(rows):
    n=len(rows)
    require(all(type(r) is int and 0<=r<2**n and not(r>>i&1) for i,r in enumerate(rows)),'bad row domain')
    require(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for j in range(n) for i in range(j)),'not a tournament')

def path_dp(rows):
    check_rows(rows);n=len(rows)
    if n==0:return 1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if mask>>v&1:
                sub=mask^(1<<v)
                if sub: dp[mask][v]=sum(dp[sub][u] for u in range(n) if sub>>u&1 and rows[u]>>v&1)
    return sum(dp[-1])

def odd_cycles(rows):
    n=len(rows);cycles=[]
    for size in range(3,n+1,2):
        for subset in itertools.combinations(range(n),size):
            first=subset[0]
            for tail in itertools.permutations(subset[1:]):
                cyc=(first,)+tail
                if all(rows[cyc[i]]>>cyc[(i+1)%size]&1 for i in range(size)):
                    cycles.append(sum(1<<v for v in cyc))
    return cycles

def saved_program(stem):
    p=json.loads((ROOT/(stem+'.ir.json')).read_text())
    fresh=Compiler().compile((ROOT/(stem+'_source.py')).read_text())
    for key in fresh: require(p[key]==json.loads(json.dumps(fresh[key])),'IR differs from source '+stem)
    emitted=''.join(f'{a}/{b}\n' for a,b,_,_ in fractions(fresh)['fractions'])
    require(emitted==(ROOT/(stem+'.fractions.00.txt')).read_text(),'fraction list differs '+stem)
    return p

def main():
    hp=saved_program('hamilton');gp=saved_program('square_shell')
    require((ROOT/'one_apex.fractions.txt').read_text()==''.join(f'{a}/{b}\n' for a,b in oa.FRACTIONS),'apex list mismatch')
    native=0;compiled=[]
    for n in range(5):
        for b in range(1<<(n*(n-1)//2)):
            code=(1<<n)*(2*b+1)-1;answer=path_dp(rows_from_bits(n,b))
            require(hamilton(code)==answer,'source path counter mismatch');native+=1
            if n<=3:
                rr=run_register(hp,code);fr=run_fractran_exponents(hp,code)
                require(rr==fr and rr['output']==answer and rr['error']==0,'compiled path mismatch')
                compiled.append({'order':n,'bits':b,'H':answer,'steps':rr['steps']})
    invalid=[]
    for n in range(5):
        for b in (1<<(n*(n-1)//2),(1<<(n*(n-1)//2))+1):
            code=(1<<n)*(2*b+1)-1
            require(hamilton(code)==0,'invalid graph not rejected')
            rr=run_register(hp,code)
            require(rr['output']==0 and rr['error']==0,'compiled invalid graph not rejected');invalid.append(code)
    universal_literal=run_fractran_literal(hp,0)
    require(universal_literal['output']==1,'literal empty graph')
    shell_runs=[]
    for q in range(1001):require(square_shell(q)==(G(q) if q else 0),'square-shell source')
    for q in range(26):
        rr=run_register(gp,q);fr=run_fractran_exponents(gp,q)
        require(rr==fr and rr['output']==(G(q) if q else 0) and rr['error']==0,'compiled square-shell')
        shell_runs.append({'q':q,'G':rr['output'],'steps':rr['steps']})
    shell_literal=[]
    for q in range(5):
        rr=run_fractran_literal(gp,q,limit=1_000_000)
        require(rr['output']==(G(q) if q else 0),'literal square-shell');shell_literal.append({'q':q,**rr})
    apex_n=apex_graph_n=apex_cycles_n=apex_literal_n=0
    max_apex_steps=0
    for m in range(9):
        for w in itertools.product((0,1),repeat=m):
            H,c=oa.expected(w);Q=0
            for b in w:Q=2*Q+b
            r=oa.run_exponents(w);E=oa.encode(w)
            exact_steps=4*E+3*Q-6*Q.bit_count()-(Q%2)+4*m-2
            require((r['H'],r['odd_cycles'])==(H,c),'apex exponent mismatch')
            require(r['steps']==exact_steps,'exact step law')
            require((r['T'],r['last'])==(Q//2,Q%2) and H==f(Q),'faithful final address')
            apex_n+=1;max_apex_steps=max(max_apex_steps,r['steps'])
            if m<=7:
                rows=rows_from_word(w);require(path_dp(rows)==H,'apex graph DP');apex_graph_n+=1
            if m<=6:
                cs=odd_cycles(rows_from_word(w))
                require(len(cs)==c and all(mask&1 for mask in cs),'apex cycle conflict invariant');apex_cycles_n+=1
            if m<=5:
                rr=oa.run_integer(w)
                require((rr['H'],rr['odd_cycles'],rr['steps'])==(H,c,exact_steps),'apex literal');apex_literal_n+=1
            wc=tuple(1-b for b in reversed(w))
            require(oa.expected(wc)==(H,c),'converse word action')
    # Positive controls and explicit quotient losses.
    h13=(22,12,24,17,2)
    require(path_dp(h13)==13 and inverse(13)==[],'one-apex scope hostile H13')
    h5=rows_from_word((1,0,0));anti=(1,0,3,2)
    require(h5==(2,12,9,1) and path_dp(h5)==5,'H5 hostile graph')
    require(all((h5[u]>>v&1)==(h5[anti[v]]>>anti[u]&1) for u in range(4) for v in range(4) if u!=v),'H5 anti-isomorphism')
    require(set(inverse(5))=={(1,0,0),(1,1,0)},'pointed H5 fibre')
    require(all(w!=tuple(1-b for b in reversed(w)) for w in inverse(5)),'pointed H5 falsely fixed')
    # First-eligible order really matters, and the carry coordinate cannot vanish.
    broken=list(oa.ROWS);broken[0],broken[2]=broken[2],broken[0]
    try:
        detected=oa.run_exponents((1,0),rows=broken)['H']!=3
    except RuntimeError:
        detected=True
    require(detected,'priority hostile undetected')
    broken=list(oa.ROWS);pos,neg=broken[14];pos=dict(pos);pos[17]=1;broken[14]=(pos,neg)
    require(oa.run_exponents((1,0,1,0),rows=broken)['H']!=9,'doubling hostile undetected')
    # Fully bounded graph enumeration is independent of the digital inverse.
    small={h:[] for h in range(3,102,2)}
    fixed={h:[] for h in range(3,102,2)}
    for m in range(2,9):
        for middle in itertools.product((0,1),repeat=m-2):
            w=(1,)+middle+(0,);H=path_dp(rows_from_word(w))
            if H<=101:
                small[H].append(w)
                if w==tuple(1-b for b in reversed(w)):fixed[H].append(w)
    for h,ws in small.items():
        require(set(ws)==set(inverse(h)),'all-order inverse finite graph check')
        require(len(fixed[h])<=1 and len(ws)%2==len(fixed[h]),'pointed fibre parity')
    reflection=[]
    for q in range(1,32):
        w=reflected_word(q);H=path_dp(rows_from_word(w))
        require(H==G(q),'direct reflected graph G');reflection.append([q,H])
    for q in range(1,10001):
        g=G(q)
        require(q*q<g<=(q+1)*(q+1),'square shell')
        require((g==(q+1)*(q+1))==(q>=2 and q&(q-1)==0),'shell endpoint iff')
    report=summary(1000000)
    require(report==json.loads((ROOT/'FIBRES.json').read_text()),'frozen fibre receipt differs')
    return {'status':'PASS','native_general_graphs':native,'compiled_general_graphs':compiled,
            'invalid_general_codes':invalid,'universal_literal':universal_literal,
            'square_shell_compiled_runs':shell_runs,'square_shell_literal':shell_literal,
            'one_apex_fraction_words':apex_n,'one_apex_graph_DP_words':apex_graph_n,
            'one_apex_cycle_words':apex_cycles_n,'one_apex_literal_words':apex_literal_n,
            'maximum_apex_steps':max_apex_steps,'direct_reflected_graphs':reflection,
            'shell_endpoint_checks':10000,'inverse_graph_bound_H':101,
            'million_inverse_summary':{k:v for k,v in report.items() if k!='first_odd_H_rows'},
            'hostiles':['priority inversion','lost doubling','unpointed H5 quotient','global H13 grammar boundary'],
            'alternating_H':[oa.expected((1,0)*m)[0] for m in range(12)],
            'balanced_H':[oa.expected((1,)*m+(0,)*m)[0] for m in range(1,12)]}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
