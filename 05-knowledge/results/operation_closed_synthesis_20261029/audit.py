"""Independent-path checks, literal operators, and source-sensitive controls.

This is a same-agent audit; it is not an independent referee report.
Run with normal Python and -O. All checks are explicit, never assert.
"""
from pathlib import Path
from fractions import Fraction
import collections,hashlib,importlib.util,itertools,json,sys
import numpy as np
import sympy as sp
from flint import fmpz_mat
import tournament_probe as tp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
checks=0
def check(ok,msg):
    global checks
    checks+=1
    if not ok:raise RuntimeError(msg)
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
alib=load('inherited_prepend',ROOT/'04-computation/experiments/growth_gf_runs_fermat_20261029/audit2/alib.py')
runs=load('inherited_permutation',ROOT/'04-computation/experiments/growth_gf_runs_fermat_20261029/runs.py')
data=json.loads((HERE/'tournament_results.json').read_text())
audit={'scope':'all 4180 grammar classes orders 4..12; specified independent subsets below','independent_matrices':{},'native_determinants':{},'controls':{}}
for ns,layer in data['layers'].items():
    n=int(ns);bykey={r['bits']:r for r in layer['records']}
    independent=0;native=0
    for idx,r in enumerate(layer['records']):
        w=r['word'];A=alib.adj_from_ascii(n,r['bits'])
        check(bool(r['det'])==('be' not in w and 'ee' not in w),'forbidden subword criterion')
        check((abs(r['det_plus'])==abs(r['det_minus'])==1) if r['det'] else (r['det_minus']==0 and abs(r['det_plus'])>1),'sector refinement')
        reverse=[0]+list(range(n-1,0,-1))
        key,_=tp.canonical(A.T,reverse)
        check(key in bykey,'converse closes the generated universe')
        check(bool(bykey[key]['det'])==bool(r['det']),'converse preserves criterion')
        # All small classes; deterministic spread of eight at each larger order.
        selected=(n<=8 or idx in {i*(len(layer['records'])-1)//7 for i in range(8)})
        if selected:
            st,M=tp.exact_stats(A)
            H=alib.path_table(A)
            Q,pr=alib.m2_matrix(A,H)
            check(np.array_equal(M,Q),'append/prepend independent path DP')
            check(alib.ham_cycles(A,H)==1,'independent HC count')
            check(int(fmpz_mat(Q.tolist()).det())==r['det'],'native integer determinant vs rational sectors')
            independent+=1;native+=1
            if n<=6:
                _,brute=runs.run_matrix(A.tolist(),2,pr)
                check(np.array_equal(Q,np.array(brute)),'literal permutation cutting vs both DPs')
    audit['independent_matrices'][ns]=independent
    audit['native_determinants'][ns]=native

# Exact six-state Moore machine: three productive, two legal-but-rejected,
# and one illegal. The measured determinant equivalence is only finite.
states=['A','B0','B1','FA','FB','X']
transition={
 'A':{'a':'A','b':'B1'},'B0':{'c':'A','d':'B0','e':'B1'},
 'B1':{'c':'A','d':'B0','e':'FB'},'FA':{'a':'FA','b':'FB'},
 'FB':{'c':'FA','d':'FB','e':'FB'},'X':{}}
output={s:(1 if s in ['A','B0','B1'] else -1 if s=='X' else 0) for s in states}
def run(s,w):
    for a in w:s=transition[s].get(a,'X')
    return output[s]
for layer in data['layers'].values():
    for r in layer['records']:check(run('B0',r['word'])==bool(r['det']),'six-state experimental detector')
witnesses={}
for s,t in itertools.combinations(states,2):
    found=None
    for k in range(3):
        for letters in itertools.product('abcde',repeat=k):
            w=''.join(letters)
            if run(s,w)!=run(t,w):found=w;break
        if found is not None:break
    check(found is not None,'all six states distinguishable')
    witnesses[s+'/'+t]=found
C=sp.Matrix([[1,0,1],[1,1,1],[1,1,0]])
x=sp.Symbol('x')
check(C.charpoly(x).as_expr()==x**3-2*x*x-x+1,'three productive states yield exact recurrence')
gen=(sp.Matrix([[0,1,0]])*(sp.eye(3)-x*C).inv()*sp.ones(3,1))[0]
check(sp.cancel(gen-(1+x)/(1-2*x-x*x+x**3))==0,'exact rational language GF')
audit['automaton']={'transition':transition,'output':output,'distinguishers':witnesses,'productive_matrix':[[int(v) for v in row] for row in C.tolist()],'generating_function':str(sp.factor(gen))}

# Literal order-five same-current-invariant / different-future witness.
r5={r['word']:r for r in data['layers']['5']['records']}
for word in ['d','e']:
    A=alib.adj_from_ascii(5,r5[word]['bits'])
    graphs=tp.children(A)
    label,key,B=next(item for item in graphs if item[0]=='e')
    st,M=tp.exact_stats(B)
    diag=[abs(int(v)) for v in fmpz_mat(M.tolist()).snf().entries()[::len(M)+1]]
    audit['controls'][word+'+e']={'source':r5[word]['bits'],'target':key,'det':st['det'],'native_smith':dict(collections.Counter(diag))}
check(audit['controls']['d+e']['det']==1 and audit['controls']['e+e']['det']==0,'one-step hostile')

# Exact endpoint measure, with owner multiplicity, for closed loneliness.
def endpoint_measure(V,q):
    points=collections.Counter()
    for v in V:
        for k in range(v):
            for sign in [-1,1]:points[Fraction(q*k+sign,q*v)%1]+=1
    good=[];mass=0
    for t,owners in sorted(points.items()):
        if all(min((v*t)%1,1-(v*t)%1)>=Fraction(1,q) for v in V):
            good.append(str(t));mass+=owners
    return {'speeds':V,'q':q,'owner_tokens':sum(points.values()),'distinct_endpoints':len(points),'good_endpoints':good,'good_owner_tokens':mass}
ap=endpoint_measure(list(range(1,14)),14)
check(ap['good_endpoints']==['1/14','3/14','5/14','9/14','11/14','13/14'],'AP exact boundary support')
check(ap['owner_tokens']==182 and ap['good_owner_tokens']==12,'AP endpoint owner multiplicities')
far=endpoint_measure([26*i for i in range(1,13)]+[339],14)
check(len(far['good_endpoints'])>0,'first-good-period-27 hostile has endpoint witness')
audit['lrc']={'tight_AP13':ap,'period27_hostile':far}

# A full-support arithmetic source measure retains a 5x+1 cycle.
def odd_step(a,n):
    v=a*n+1
    while v%2==0:v//=2
    return v
cyc=[13,33,83]
check([odd_step(5,n) for n in cyc]==[33,83,13],'positive nonroot 5x+1 cycle')
weights=[Fraction(1,2*4**(((n+1)//2).bit_length()-1)) for n in cyc]
check(sum(weights)==Fraction(69,2048),'source atomic lower bound')
audit['collatz']={'map':'(5n+1)/2^v2(5n+1)','cycle':cyc,'source_weights':list(map(str,weights)),'survivor_mass_lower_bound':str(sum(weights))}

old=load('inherited_81_27',HERE.parent/'four_vertex_glue_fermat_20261029/verify.py')
audit['controls']['81_27']=old.hostile_control()
check(audit['controls']['81_27']['automorphisms']==[81,27],'decorated triple hostile retained')
def fixed(word):
    A,C,D=1,0,1
    for a in word:A,C,D=3*A,3*C+D,D*(1<<a)
    return Fraction(C,D-A)
pair=[fixed((1,1,5,1)),fixed((1,3,3,1))]
check(pair==[Fraction(37,35),Fraction(221,175)],'same future different arithmetic response')
audit['controls']['35_175']=[str(x) for x in pair]

audit['checks']=checks
audit['tournament_data_sha256']=hashlib.sha256((HERE/'tournament_results.json').read_bytes()).hexdigest()
(HERE/'audit_results.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'checks':checks,'matrix_comparisons':sum(audit['independent_matrices'].values()),'lrc_good':{k:len(v['good_endpoints']) for k,v in audit['lrc'].items()}},indent=2))
