"""Test constant signed chart-multiplicity formulas on a declared finite range."""
import importlib.util
import json
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
source = ROOT/'05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py'
spec = importlib.util.spec_from_file_location('core_compiler',source)
compiler = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compiler)
rows = json.loads((ROOT/'05-knowledge/results/chanpak_five_core_bridges_20261110/cover/signatures.json').read_text())


def exact_rank(matrix):
    """Primitive-integer elimination, exactly the rational row rank."""
    basis = {}
    for row in matrix:
        v = list(row)
        for p, b in sorted(basis.items()):
            if v[p]:
                q = gcd(v[p],b[p])
                u,t = b[p]//q,v[p]//q
                v = [u*x-t*y for x,y in zip(v,b)]
                z = 0
                for x in v:
                    z = gcd(z,x)
                if z:
                    v = [x//z for x in v]
        p = next((j for j,x in enumerate(v) if x),None)
        if p is not None:
            basis[p] = v
    return len(basis)


def first_inconsistency(matrix,target,values):
    basis={}
    chosen=[]
    columns=len(matrix[0])
    for row,t,value in zip(matrix,target,values):
        v=list(row)+[t]
        for p,b in sorted(basis.items()):
            if v[p]:
                q=gcd(v[p],b[p]);u,t=b[p]//q,v[p]//q
                v=[u*x-t*y for x,y in zip(v,b)]
                z=0
                for x in v:z=gcd(z,x)
                if z:v=[x//z for x in v]
        p=next((j for j,x in enumerate(v) if x),None)
        if p==columns:
            return {'first_inconsistent_output':value,'equation_outputs':chosen+[value]}
        if p is not None:
            basis[p]=v;chosen.append(value)
    return None


def main():
    bound = 1023
    values = list(range(1,bound+1,2))
    matrix = [[0]*(4*len(rows)) for _ in values]
    for ri,row in enumerate(rows):
        for h,a,b in compiler.generate(row['signature'],bound):
            matrix[h//2][4*ri+2*(a&1)+(b&1)] += 1
    result = {'bound':bound,'universe':'44 canonical signature rows; inherited finite inverse; degenerate axes normalized to even addresses'}
    target = [int(h not in (7,21)) for h in values]
    for label,A in [('row_only',[[sum(r[j:j+4]) for j in range(0,len(r),4)] for r in matrix]),('row_and_parity',matrix)]:
        # Rational row-space arithmetic: not a modular-rank shortcut.
        rank = exact_rank(A)
        augmented = exact_rank([r+[t] for r,t in zip(A,target)])
        result[label] = {'variables':len(A[0]),'equations':len(A),'coefficient_rank':rank,'augmented_rank':augmented,'consistent':rank==augmented,'contradiction':first_inconsistency(A,target,values)}
        print(label,result[label],flush=True)
    Path(__file__).with_name('chart_counts.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':
    main()
