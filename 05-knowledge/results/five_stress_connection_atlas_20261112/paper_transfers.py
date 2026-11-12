"""Independent finite transfer checks of the supplied king-strip paper.

Only the displayed local conditions are inputs. We do not treat the paper's
automaton as a classifier of globally resolving sets. Also check the new
marked-pair seed language and its aggregate counting laws.
"""
from itertools import product
from collections import defaultdict
from pathlib import Path
import json,hashlib
import sympy as s
HERE=Path(__file__).resolve().parent
def need(x,msg):
    if not x:raise RuntimeError(msg)
def bits(t):return ((t>>2)&1,(t>>1)&1,t&1)
def legal(u,v,w):
    p,m,q=bits(v);pu,mu,qu=bits(u);pw,mw,qw=bits(w)
    return (p,pu+pw)!=(q,qu+qw) and (p!=m or qu+q+qw>=1) and (q!=m or pu+p+pw>=1)
def reach(starts,edges):
    seen=set(starts);todo=list(starts)
    while todo:
        v=todo.pop()
        for w in edges.get(v,[]):
            if w not in seen:seen.add(w);todo.append(w)
    return seen
def resolving(word):
    pts=[(x,y) for x,t in enumerate(word) for y,b in enumerate(bits(t)) if b]
    observed=set()
    for x in range(len(word)):
        for y in range(3):
            code=tuple(sorted(max(abs(x-a),abs(y-b)) for a,b in pts))
            if code in observed:return False
            observed.add(code)
    return True
def run():
    states=list(product(range(8),repeat=2));infty=10**9
    transitions=[(u,v,w) for u,v,w in product(range(8),repeat=3) if legal(u,v,w)]
    V={st:(st[1].bit_count() if st[0]==0 else infty) for st in states};vectors={1:V}
    for n in range(2,6):
        nxt={st:infty for st in states}
        for u,v,w in transitions:
            if V[u,v]<infty:nxt[v,w]=min(nxt[v,w],V[u,v]+w.bit_count())
        V=nxt;vectors[n]=V
    V4=vectors[4]
    need(all(vectors[5][q]==(V4[q]+1 if V4[q]<infty else infty) for q in states),'min-plus eigenvector')
    tight=[((u,v),(v,w)) for u,v,w in transitions if V4[u,v]<infty and V4[v,w]<infty and w.bit_count()+V4[u,v]-V4[v,w]==1]
    need(len(tight)==184,'tight transition count')
    prefixes=[]
    for word in product(range(8),repeat=4):
        aug=(0,)+word
        if all(legal(*aug[i:i+3]) for i in range(3)) and sum(v.bit_count() for v in word)==V4[word[-2:]]:prefixes.append(word)
    need(len(prefixes)==320,'minimum prefixes')
    accept={q for q in states if legal(*q,0) and V4[q]==4}
    edges=defaultdict(list);back=defaultdict(list)
    for u,v in tight:edges[u].append(v);back[v].append(u)
    live=reach({w[-2:] for w in prefixes},edges)&reach(accept,back)
    need(len(live)==19 and len(accept)==12,'trimmed automaton')
    need(all(reach({v},edges)&live==live for v in live),'strongly connected core')
    ordered=sorted(live);idx={q:i for i,q in enumerate(ordered)}
    A=s.zeros(19)
    for u,v in tight:
        if u in live and v in live:A[idx[u],idx[v]]+=1
    eig=(A-3*s.eye(19)).nullspace();need(len(eig)==1,'Perron eigenspace')
    ev=eig[0];den=s.ilcm(*(x.q for x in ev));ev=ev*den
    if ev[0]<0:ev=-ev
    need(all(x>0 for x in ev),'positive Perron vector')
    # Reflection of top and bottom rows is an actual commuting involution.
    def flip(t):
        p,m,q=bits(t);return (q<<2)|(m<<1)|p
    J=s.zeros(19)
    for q,i in idx.items():J[idx[tuple(flip(v) for v in q)],i]=1
    need(J*J==s.eye(19) and J*A==A*J,'reflection action')
    even=s.Matrix.hstack(*(J-s.eye(19)).nullspace())
    odd=s.Matrix.hstack(*(J+s.eye(19)).nullspace())
    evenA=even.gauss_jordan_solve(A*even)[0]
    oddA=odd.gauss_jordan_solve(A*odd)[0]
    z=s.Symbol('z')
    ranks=[(A**k).rank() for k in range(1,5)]
    need(ranks==[9,7,5,5],'persistent state ranks')
    start=s.zeros(1,19)
    for w in prefixes:
        if w[-2:] in live:start[idx[w[-2:]]]+=1
    finish=s.Matrix([int(q in accept) for q in ordered])
    need(start*J==start and J*finish==finish,'even input and observer')
    generating=s.factor((start*(s.eye(19)-z*A).inv()*finish)[0])
    expected=4*(15+17*z+30*z*z)/((1-3*z)*(1+z+2*z*z))
    need(s.cancel(generating-expected)==0,'exact accepted-word generating function')
    # Enumeration of every accepted word, followed by the global metric test.
    words=[w for w in prefixes if w[-2:] in live];counts={};global_counts={}
    for n in range(4,11):
        accepted=[w for w in words if w[-2:] in accept]
        counts[n]=len(accepted)
        global_counts[n]=sum(resolving(w) for w in accepted)
        if n<10:words=[w+(v[1],) for w in words for v in edges[w[-2:]] if v in live]
    need([counts[n] for n in range(5,9)]==[188,556,1660,5004],'accepted counts')
    need([global_counts[n] for n in range(5,9)]==[0,32,320,1328],'global resolving counts')
    # Five-state minimal DFA for 0*1(epsilon|0)10*; fifth state is dead.
    dfa=((0,1),(2,3),(4,3),(3,4),(4,4));classes=[{3},{0,1,2,4}]
    while True:
        cl={q:i for i,c in enumerate(classes) for q in c};parts={}
        for q in range(5):parts.setdefault((cl[q],cl[dfa[q][0]],cl[dfa[q][1]]),set()).add(q)
        new=list(parts.values())
        if len(new)==len(classes):break
        classes=new
    need(len(classes)==5,'minimal seed DFA')
    N=z/(1-z)**2+z**3/(1-z)+2*z**5*(1+z)/(1-z)**2
    E=2*z**4+z**5+2*z**6+2*z**8+2*z**9+2*z**11
    M=z/((1-z)*(1-2*z))+6*z**3/(1-2*z)+z**5*((80+208*z)/((1-2*z)*(1-3*z))-4*(1+z)/(1-z))
    ME=27*z**4+37*z**5+210*z**6+1437*z**8+9006*z**9+594887*z**11
    receipt=json.loads((HERE/'tournament_receipt.json').read_text());records=receipt['seed_records']
    if (HERE/'round_extension.json').exists():
        ext=json.loads((HERE/'round_extension.json').read_text());records+=ext['additional_seed_records']
    max_n=max(r['n'] for r in records)
    nc=s.series(N+E,z,0,max_n+1).removeO().expand();mc=s.series(M+ME,z,0,max_n+1).removeO().expand()
    seq=[]
    for n in range(1,max_n+1):
        rr=[r for r in records if r['n']==n]
        need(len(rr)==int(nc.coeff(z,n)),'family-count conjecture finite control')
        need(sum(r['U'] for r in rr)==int(mc.coeff(z,n)),'slope-moment conjecture finite control')
        seq.append({'n':n,'types':len(rr),'sum_U':sum(r['U'] for r in rr)})
    pdf=Path(r'C:\Users\Eliott\Downloads\The_Multiset_Dimension_of_Graphs_Extremal_Values_a.pdf')
    out={'scope':'Independent checks of supplied paper sections14-16 on local necessary conditions; global resolving check on all accepted words through10. No audit of the billion-graph census. Seed-family generating laws checked against complete round census only to stated order.','pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'tight_edges':len(tight),'minimum_prefixes':len(prefixes),'core_states':ordered,'core_minimum_prefixes':sum(w[-2:] in live for w in prefixes),'accepting_states':sorted(accept),'core_characteristic_polynomial':str(A.charpoly(z).as_expr().factor()),'even_characteristic_polynomial':str(evenA.charpoly(z).as_expr().factor()),'odd_characteristic_polynomial':str(oddA.charpoly(z).as_expr().factor()),'ranks_of_A_powers_1_to_4':ranks,'accepted_word_GF_from_n4':str(generating),'positive_eigenvector_for_3':[int(x) for x in ev],'accepted_counts':counts,'global_resolving_counts':global_counts,'seed_dfa':dfa,'minimal_dfa_states':len(classes),'family_count_GF':str(s.factor(N)),'finite_exception_count_polynomial':str(E),'family_sum_U_GF':str(s.factor(M)),'finite_exception_U_polynomial':str(ME),'round_count_and_moment_controls':seq}
    tag='.optimized' if not __debug__ else ''
    (HERE/('paper_transfer_receipt'+tag+'.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)
if __name__=='__main__':run()
