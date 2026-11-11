"""Independent return-lane audit of the companion Berggren statements."""
from pathlib import Path
from itertools import permutations,product
import importlib.util,json,sys
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('berggren_checks',HERE.parent/'berggren/check.py')
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
sys.path.insert(0,str(HERE.parent.parent/'uniform_return_metagraph_20261111'/'reset'))
from classify_modules import H,need,converse,inherited as a
from ready_rebase import full_state

def det(M):
    out=0
    for p in permutations(range(len(M))):
        sign=(-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
        term=sign
        for i,j in enumerate(p):term*=M[i][j]
        out+=term
    return out

def main():
    prefixes=['','O','OI','OO'];suffixes=['','I','II','OI']
    hankel=[[H(b.graph(p+s)) for s in suffixes] for p in prefixes]
    need(hankel==[[1,1,1,3],[1,3,5,5],[3,5,9,9],[1,5,11,9]] and det(hankel)==16,
         'independent actual-path full-rank certificate')
    controls=0
    for m in range(3,8):
        for middle in product('OI',repeat=m-3):
            w='OO'+''.join(middle)+'I';T=b.graph(w)
            G=a.birth_beaten_only_by(a.birth_beaten_only_by(T,0),len(T))
            h=H(T);P=b.state(w)[1];Q=b.reverse_complement_value(w)
            need(full_state(G,m)==[3*h+2,3,1,3*(P+1),6*Q+9,3,6*(1<<m)],
                 'independent complete old-apex response on every OO-start strong cut through7')
            controls+=1
    # Concatenation leaves only these four coefficients at an arbitrary suffix.
    hp,pp,ep,_=b.state('OOOO');hn,pn,en,_=b.state('OIOOIO')
    residual=[-2,pn-3*pp,en-3*ep,hn-3*hp+2]
    need(residual==[-2,0,-2,14],'all-suffix symbolic eligibility residual')
    dual=[]
    for t in range(1,4):
        word='OIOOIO'+'I'*t+'OII';T=b.graph(word);C=converse(T);root=len(word)-1
        p=a.root_profiles(C)[root]
        need(p['gap']==1 and p['h']==315*(1<<t)-25,'independent converse-ready escape')
        child=a.birth_beaten_only_by(a.birth_beaten_only_by(C,root),len(C))
        need(H(child)==945*(1<<t)-73,'actual reverse-gauge additional+2 exit')
        dual.append(dict(t=t,root=root,H=p['h'],target_H=H(child)))
    report=dict(status='PASS',hankel=hankel,determinant=det(hankel),
                full_old_apex_controls=controls,arbitrary_suffix_coefficients=residual,dual_continuations=dual,
                logical_scope={'rank':'homogeneous linear append-word series only; nonlinear two-register carrier remains possible',
                               'future_order':'common O/Z/J continuations at retained oldapex only; no unproved seven-state context order',
                               'readiness':'forward no-ready claim is orientation-specific; OII has converse-ready exit'})
    (HERE/'audit_berggren.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: independent Berggrenrank/interface/suffix audit; converseescape explicitly retained')

if __name__=='__main__':main()
