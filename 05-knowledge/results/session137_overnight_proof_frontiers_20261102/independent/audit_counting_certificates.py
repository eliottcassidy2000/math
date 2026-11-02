"""Audit a supplied verifier against independent gaps and all one-bit mutants.

Reproduce: python3 audit_counting_certificates.py ../global_referee/counting_gaps.py
The target verifier/encoder is intentionally loaded as the implementation under
test. Expected numerical gaps use this reviewer's separate kernels, not target
literal_expected or target determinant.
"""
import importlib.util,json,sys
from pathlib import Path
from hashlib import sha256
from tournament_kernel import *
from audit_feedback2 import determinant

def expected(a):
    n=len(a)
    if not n:return [0]*4
    h=literal_hp(a);u=determinant([[a[i][j]+(i==j) for j in range(n)] for i in range(n)])
    v=(-1)**n*determinant(a);cycles=proper=0
    for k in range(3,n+1):
        for vs in combinations(range(n),k):
            for tail in permutations(vs[1:]):
                p=(vs[0],)+tail
                if all(a[p[i]][p[(i+1)%k]] for i in range(k)):cycles+=1;proper+=k<n
    delta=h-u
    return list(map(int,(delta,delta-cycles,delta-abs(v),delta-abs(v)-proper)))

def ordinal(a,b):
    n,m=len(a),len(b)
    return tuple(row+(1,)*m for row in a)+tuple((0,)*n+row for row in b)

def run():
    source=Path(sys.argv[1] if len(sys.argv)>1 else '../global_referee/counting_gaps.py').resolve()
    spec=importlib.util.spec_from_file_location('counting_verifier_under_test',source);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    inputs=certificates=mutants=accepted_mutants=0;controls=[]
    def check(a):
        nonlocal inputs,certificates,mutants,accepted_mutants
        S=module.GapSystem(a);values=expected(a)
        for mode in range(4):
            words=S.certificates(mode);bank=set(words)
            require(len(words)==len(bank)==values[mode],'independent certificate count')
            for word in words:
                require(S.verify(mode,word),'genuine certificate rejected');certificates+=1
                for i in range(len(word)):
                    mutation=word[:i]+str(1-int(word[i]))+word[i+1:]
                    accepted=S.verify(mode,mutation);mutants+=1;accepted_mutants+=accepted
                    require(accepted==(mutation in bank),'one-bit parser/counting mismatch')
                require(not S.verify(mode,word+'0') and not S.verify(mode,word[:-1]),'length control')
                require(not S.verify(mode,'x'+word[1:]),'nonbinary control')
        inputs+=1;return values
    for n in range(5):
        for code in range(1<<(n*(n-1)//2)):check(matrix(n,code))
    c3=matrix(3,2);two=ordinal(c3,c3);source_two=ordinal(matrix(1,0),two)
    controls.append({'name':'C3=>C3','gaps':check(two)})
    controls.append({'name':'source=>C3=>C3','gaps':check(source_two)})
    require(controls[0]['gaps']==[5,3,4,2] and controls[1]['gaps']==[5,3,5,3],'positive gamma controls')
    for bad in ([[0,1],[1,0]],[[1]],[[0],[0,0]],[]):
        S=module.GapSystem(bad)
        for mode in range(4):require(not S.verify(mode,'0'*(S.n*S.record_bits)),'malformed/empty zero extension')
    print(json.dumps({'scope':'implementation-under-test audit; independent literal expected gaps for every labelled order0..4 plus two positive-Gamma controls, and every single-bit mutation of every valid certificate',
                      'verifier_source_sha256':sha256(source.read_bytes()).hexdigest(),
                      'inputs':inputs,'valid_certificates':certificates,'single_bit_mutants':mutants,
                      'mutants_that_are_other_valid_certificates':accepted_mutants,'controls':controls,
                      'verdict':'PASS exact expected counts, canonical uniqueness, exhaustive local parser mutations, and zero extensions'},sort_keys=True,indent=2))
if __name__=='__main__':run()
