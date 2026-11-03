"""Complete weighted finite birth windows; no sampling and no asserted gates."""
from pathlib import Path
from collections import Counter,deque
from math import isqrt
import json

HERE=Path(__file__).resolve().parent
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def graph(n):
    return [(i,j) for j in range(2,n+1) for i in range(1,j) if isqrt(i+j)**2==i+j]
def bank(n):
    E=graph(n);adj={i:[] for i in range(1,n+1)}
    for a,b in E:adj[a].append(b);adj[b].append(a)
    leaves=[i for i in adj if len(adj[i])==1]
    need(leaves==[18],'declared path-bank source range')
    paths=[]
    def visit(p,mask):
        if len(p)==n:paths.append(p);return
        for v in adj[p[-1]]:
            if not (mask>>(v-1)&1):visit(p+(v,),mask|(1<<(v-1)))
    visit((18,),1<<17)
    paths=sorted(paths+[p[::-1] for p in paths]);idx={e:k for k,e in enumerate(E)}
    pats=[]
    for p in paths:
        m=v=0
        for a,b in zip(p,p[1:]):
            k=idx[tuple(sorted((a,b)))];m|=1<<k
            if a<b:v|=1<<k
        pats.append((m,v))
    active=set()
    for mask,value in pats:
        words=[value]
        for k in range(len(E)):
            if not(mask>>k&1):words += [w|(1<<k) for w in words]
        active.update(words)
    legal={x:tuple(i for i,(m,v) in enumerate(pats) if x&m==v) for x in sorted(active)}
    return {'n':n,'E':E,'P':paths,'patterns':pats,'active':active,'legal':legal,'full':(1<<len(E))-1}
def compress(keys):
    table={}
    return [table.setdefault(k,len(table)) for k in keys]
def refine(obs,succ):
    labels=compress(obs);rounds=[max(labels)+1]
    while True:
        nxt=compress((labels[i],tuple((a,labels[j]) for a,j in succ[i])) for i in range(len(obs)))
        rounds.append(max(nxt)+1)
        if nxt==labels:return labels,rounds
        labels=nxt
def split_refine(obs,succ):
    """Independent block splitting by inverse images of named transitions."""
    init=compress(obs);blocks={k:set() for k in set(init)}
    for i,k in enumerate(init):blocks[k].add(i)
    loc=init[:];incoming=[[] for _ in succ]
    for i,row in enumerate(succ):
        for a,j in row:incoming[j].append((a,i))
    pending=deque(blocks);queued=set(blocks)
    while pending:
        b=pending.popleft();queued.discard(b)
        hits={}
        for j in blocks[b]:
            for a,i in incoming[j]:hits.setdefault(a,set()).add(i)
        for pre in hits.values():
            touched={}
            for i in pre:touched.setdefault(loc[i],set()).add(i)
            for c,part in touched.items():
                if len(part)==len(blocks[c]):continue
                other=blocks[c]-part
                if len(part)>len(other):part,other=other,part
                blocks[c]=other;d=max(blocks)+1;blocks[d]=part
                for i in part:loc[i]=d
                if c not in queued:pending.append(c);queued.add(c)
                pending.append(d);queued.add(d)
    return loc
def same_partition(a,b):return len(set(zip(a,b)))==len(set(a))==len(set(b))
def compile_stage(g,target=None,independent=True):
    m=len(g['E']);full=g['full'];active=g['active']
    relevant=active if target is None else active|{z&full for z in target['relevant']}
    domain=[-1]+sorted(relevant);index={x:i for i,x in enumerate(domain)}
    new_edges=[] if target is None else target['g']['E'][m:]
    if target is not None:need(target['g']['E'][:m]==g['E'],'literal nested edge order')
    targets=[()] * len(domain)
    if target is not None:
        target_dead=target['class_by_source'][-1]
        targets=[tuple(target['class_by_source'].get(x|(a<<m),target_dead) for a in range(1<<len(new_edges)))
                 if x!=-1 else (target_dead,)*(1<<len(new_edges)) for x in domain]
    names=[g['legal'].get(x,()) for x in domain]
    succ=[tuple((a,index[x^(full^g['patterns'][a][0])]) for a in L) for x,L in zip(domain,names)]
    obs=list(zip(names,targets));labels,rounds=refine(obs,succ)
    if independent:
        alt=split_refine(obs,succ);need(same_partition(labels,alt),'independent partition refinement')
    weights=[(1<<m)-len(relevant)]+[1]*len(relevant)
    need(weights[0]>0,'symbolic residual exists')
    groups={}
    for i,c in enumerate(labels):groups.setdefault(c,[]).append(i)
    rows=[];mass=Counter()
    for c,members in groups.items():
        i=members[0];x=domain[i];row={'class':c,'representative':x,'weight':sum(weights[j] for j in members),
             'legal_paths':names[i],'complements':[(a,labels[j]) for a,j in succ[i]],'birth_targets':targets[i],
             'converse':labels[index[x^full if x!=-1 else -1]]}
        for j in members:
            need(names[j]==names[i] and targets[j]==targets[i],'observation closure')
            need([(a,labels[k]) for a,k in succ[j]]==row['complements'],'command closure')
            y=domain[j];need(labels[index[y^full if y!=-1 else -1]]==row['converse'],'converse descends')
        rows.append(row)
        for t in targets[i]:mass[t]+=row['weight']
    w={r['class']:r['weight'] for r in rows}
    need(sum(w.values())==1<<m,'source mass')
    for r in rows:
        need(w[r['converse']]==r['weight'],'converse preserves weight')
        for a,t in r['complements']:need(w[t]==r['weight'],'partial involution preserves weight')
    if target is not None:
        need(mass==Counter({r['class']:r['weight'] for r in target['rows']}),'all birth pushforward weights')
        for i,x in enumerate(domain):
            y=index[x^full if x!=-1 else -1]
            for a,t in enumerate(targets[i]):
                need(targets[y][a^((1<<len(new_edges))-1)]==target['rows'][t]['converse'],'birth converse square')
    result={'g':g,'relevant':relevant,'domain':domain,'labels':labels,'class_by_source':dict(zip(domain,labels)),
            'rows':rows,'rounds':rounds,'new_edges':new_edges,'residual_weight':weights[0],
            'explicit_relevant':len(relevant),'active_prebirth_classes':len({labels[index[x]] for x in active})}
    return result
def export(stage):
    return {'n':stage['g']['n'],'edges':stage['g']['E'],'directed_paths':stage['g']['P'],
       'active_sources':len(stage['g']['active']),'explicit_relevant':stage['explicit_relevant'],
       'symbolic_residual_weight':stage['residual_weight'],'new_edges':stage['new_edges'],
       'refinement':stage['rounds'],'states':len(stage['rows']),
       'active_prebirth_classes':stage['active_prebirth_classes'],'controller':stage['rows']}

def projected_support(old,banks):
    """Construct old support directly, never completing future free edge bits."""
    full=old['full'];result=set()
    for future in banks:
        for mask,value in future['patterns']:
            mask &= full;value &= full;words=[value]
            for k in range(len(old['E'])):
                if not(mask>>k&1):words += [x|(1<<k) for x in words]
            result.update(words)
    return result
def main():
    G={n:bank(n) for n in (25,26,27)}
    need([len(G[n]['P']) for n in G]==[20,24,70],'inherited complete banks')
    fixed={n:compile_stage(G[n]) for n in G}
    need([len(fixed[n]['rows']) for n in G]==[33,71,637],'inherited fixed-order minima')
    one25=compile_stage(G[25],fixed[26]);need(len(one25['rows'])==97,'inherited first birth')
    one26=compile_stage(G[26],fixed[27]);two25=compile_stage(G[25],one26)
    need([len(one26['rows']),len(two25['rows'])]==[761,903],'new finite minima')
    need(projected_support(G[26],[G[26],G[27]])==one26['relevant'],
         'direct old forest support for one birth')
    need(projected_support(G[25],[G[25],G[26],G[27]])==two25['relevant'],
         'direct old forest support for two births')
    # A literal same-old-controller witness separated by the next birth.
    witness=None
    for x in sorted(G[26]['active']):
        for k in range(len(G[26]['E'])):
            y=x^(1<<k)
            if y<=x or y not in G[26]['active']:continue
            if fixed[26]['class_by_source'][x]!=fixed[26]['class_by_source'][y]:continue
            if one26['class_by_source'][x]==one26['class_by_source'][y]:continue
            m=len(G[26]['E']);answers=[]
            for a in range(4):
                xx=x|(a<<m);yy=y|(a<<m)
                Lx=G[27]['legal'].get(xx,());Ly=G[27]['legal'].get(yy,())
                if Lx!=Ly:
                    # Independently inspect directed graph arcs, without the bit pattern masks.
                    def literal(z):
                        arcs={e if z>>b&1 else e[::-1] for b,e in enumerate(G[27]['E'])}
                        return tuple(i for i,p in enumerate(G[27]['P']) if all(e in arcs for e in zip(p,p[1:])))
                    need(literal(xx)==Lx and literal(yy)==Ly,'literal hostile path check')
                    answers.append({'new_bits':a,'left_paths':Lx,'right_paths':Ly})
            if answers:
                witness={'x':x,'y':y,'flipped_edge':G[26]['E'][k],
                         'old_state':fixed[26]['class_by_source'][x],'birth_answers':answers};break
        if witness:break
    need(witness is not None,'one-edge birth hostile')
    incoming=sum((z&G[26]['full']) in G[26]['active'] for z in G[27]['active'])
    old25_runs=Counter()
    for p in G[27]['P']:
        inside=False;runs=0
        for v in p:
            if v>25:inside=False
            elif not inside:runs+=1;inside=True
        old25_runs[runs]+=1
    need(old25_runs==Counter({3:62,2:8}),'literal deletion to old25 forests')
    result={'status':'FINITE-EXACT full guarded orientation universe with symbolic dead mass',
       'fixed':{n:export(s) for n,s in fixed.items()},'birth25':export(one25),
       'birth26':export(one26),'two_birth25':export(two25),'witness26':witness,
       'new27_active_from_old26_active':incoming,'new27_active_from_old26_dead':len(G[27]['active'])-incoming,
       'new27_paths_by_old25_runs':dict(old25_runs)}
    (HERE/'birth.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'fixed':{n:len(s['rows']) for n,s in fixed.items()},
         'birth25':len(one25['rows']),'birth26':len(one26['rows']),'two_birth25':len(two25['rows']),
         'rounds26':one26['rounds'],'rounds25_two':two25['rounds'],
         'relevant26':len(one26['relevant']),'relevant25_two':len(two25['relevant']),
         'new27_active':len(G[27]['active']),'new27_from_old_active':incoming,
         'witness26':witness},sort_keys=True))
if __name__=='__main__':main()
