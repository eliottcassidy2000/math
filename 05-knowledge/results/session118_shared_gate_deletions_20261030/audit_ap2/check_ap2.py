"""Independent literal shared-DAG AP2 sharpness check; no producer imports."""
import json
from itertools import product

def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def bitgate(t,a,b):return (t>>(a+2*b))&1
def make(m):
    ins=['p','w','r','h','z','t','u']+[f'{c}{i}' for i in range(m) for c in 'xy']
    gates=[('B','p','r',2),('D','w','h',2),('G','p','w',8),('Hr','r','z',8),('Hh','h','t',8),('K','G','u',6),('J','Hr','Hh',8),('O','J','K',6)]
    for i in range(m):gates.append((f'A{i}',f'x{i}',f'y{i}',8))
    acc='A0'
    for i in range(1,m):
        nxt=f'P{i}';gates.append((nxt,acc,f'A{i}',6));acc=nxt
    gates.extend([('Q','O',acc,6),('C','Q','B',6),('E','C','D',6)])
    return ins,gates,'E'
def ev(gates,out,assignment):
    values=dict(assignment)
    for name,a,b,t in gates:values[name]=bitgate(t,values[a] if isinstance(a,str) else a,values[b] if isinstance(b,str) else b)
    return values[out] if isinstance(out,str) else out

def simplify(gates,out,fixed):
    # Constants and positive literal rerouting suffice for this constructed trace.
    alias=dict(fixed);kept=[];deleted=[]
    def get(v):
        while isinstance(v,str) and v in alias:v=alias[v]
        return v
    for name,a,b,t in gates:
        a,b=get(a),get(b)
        if name in fixed:
            deleted.append(name);continue
        variables=sorted(set(x for x in [a,b] if isinstance(x,str)))
        if len(variables)<=1:
            if not variables:alias[name]=bitgate(t,a,b)
            else:
                v=variables[0]
                vals=[bitgate(t,z if a==v else a,z if b==v else b) for z in [0,1]]
                if vals==[0,1]:alias[name]=v
                elif vals[0]==vals[1]:alias[name]=vals[0]
                else:raise RuntimeError('unexpected complemented literal')
            deleted.append(name)
        else:kept.append((name,a,b,t))
    return kept,get(out),deleted

def usages(gates):
    result={}
    for name,a,b,t in gates:
        for p in set([a,b]):
            if isinstance(p,str):result.setdefault(p,set()).add(name)
    return result

def run(m):
    ins,g0,o0=make(m)
    g1,o1,d1=simplify(g0,o0,{'B':0})
    g2,o2,d2=simplify(g1,o1,{'D':0})
    g3,o3,d3=simplify(g2,o2,{'r':1,'h':1,'w':0})
    require(d1==['B','C'],d1);require(d2==['D','E'],d2)
    require(set(d3)=={'G','Hr','Hh','K'},d3)
    child=usages(g2)
    require(child['p']==child['w']=={'G'},'paired first parents')
    require(child['r']=={'Hr'} and child['h']=={'Hh'},'paired partner children')
    require(child['G']=={'K'},'shared gate single child')
    require(len(g2)==len(ins)-1,'minimum core count')
    require(len(g3)==2*m+2,'minimum residual count')
    witnesses=set();residual_witnesses=set();checks=0
    for bits in product([0,1],repeat=len(ins)):
        a=dict(zip(ins,bits)); f=ev(g0,o0,a)
        allowed1=not(a['p']==1 and a['r']==0)
        allowed2=allowed1 and not(a['w']==1 and a['h']==0)
        if allowed1:require(f==ev(g1,o1,a),'first GP semantics')
        if allowed2:
            val=ev(g2,o2,a);require(f==val,'second GP semantics')
            for var in ins:
                b=dict(a);b[var]^=1
                if (not(b['p']==1 and b['r']==0)) and (not(b['w']==1 and b['h']==0)) and ev(g2,o2,b)!=val:witnesses.add(var)
        if a['r']==a['h']==1 and a['w']==0:
            val=ev(g3,o3,a);require(f==val,'AP2 semantics')
            for var in ins:
                if var in ['r','h','w']:continue
                b=dict(a);b[var]^=1
                if ev(g3,o3,b)!=val:residual_witnesses.add(var)
        checks+=1
    require(witnesses==set(ins),'all core variables essential on paired domain')
    require(residual_witnesses==set(ins)-{'p','w','r','h'},'residual essential variables')
    d=m+7
    a_values=[len(ins),len(ins)-1,len(ins)-2,len(ins)-3]
    require(all(a>=d+2 for a in a_values[:3]),'loop guards')
    require(a_values[-1]>d,'nonconstant terminal guard')
    return {'m':m,'input_bits':len(ins),'full_truth_table_assignments':checks,'gate_counts':[len(g) for g in [g0,g1,g2,g3]],'deleted':[d1,d2,d3],'a_values':a_values,'disperser_dimension_upper_bound':d,'core_essential_inputs':len(witnesses),'residual_essential_inputs':len(residual_witnesses)}
print(json.dumps({'status':'PASS','scope':'literal reachable minimum-core AP2 sharpness; general lemma proved separately','cases':[run(4)]},indent=2,sort_keys=True))
