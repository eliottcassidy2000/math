"""Small literal controls for the new rooted constructions and loop formula.

No spectrum coverage scan. Inherited Q-polynomials are not re-audited.
Run from any directory: python3 -B /absolute/path/to/check_positive.py
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import json
from find_intersections import g, rev, append_coeff

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
COMPILER = ROOT / '05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py'
spec = importlib.util.spec_from_file_location('compiler', COMPILER)
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)

def state(rows, root):
    n = len(rows)
    full = (1 << n) - 1
    h = c.hp_subsets(rows)
    dp = [[0] * n for _ in range(full + 1)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full + 1):
        for v, x in enumerate(dp[mask]):
            if not x:
                continue
            avail = rows[v] & ~mask & full
            while avail:
                bit = avail & -avail
                avail -= bit
                dp[mask | bit][bit.bit_length() - 1] += x
    E = dp[full][root]
    start = [[0] * n for _ in range(full + 1)]
    start[1 << root][root] = 1
    for mask in range(1, full + 1):
        for v, x in enumerate(start[mask]):
            if not x:
                continue
            avail = rows[v] & ~mask & full
            while avail:
                bit = avail & -avail
                avail -= bit
                start[mask | bit][bit.bit_length() - 1] += x
    U = sum(sum(start[mask]) * h[full ^ mask] for mask in range(full + 1))
    return (h[full], E, U, h[full ^ (1 << root)])

def q(s):
    H, E, U, d = s
    return H, U, U-d, H+U-E

def join(a, b):
    shift = len(a)
    return [x | (((1 << len(b))-1) << shift) for x in a] + [x << shift for x in b]

def source_root(rows):
    return list(rows) + [(1 << len(rows))-1]

def Z(rows, root):
    n = len(rows)
    return [x | (0 if v == root else 1 << n) for v, x in enumerate(rows)] + [1 << root]

def delete(rows, v):
    keep = [i for i in range(len(rows)) if i != v]
    return [sum(((rows[i] >> j) & 1) << jj for jj, j in enumerate(keep)) for i in keep]

def loop_prediction(p, r, length, k):
    M = 1 << length
    Q = M-1
    a, b, z = append_coeff(r, length)
    D = F(b*(r&1)+z)-F(a*r, Q)
    n = F(p) + F(r,Q)
    return F(a,Q)*n*M**k + D*k + g(p)-F(a,Q)*n+b*((p&1)-(r&1))

def main():
    A = (56,9,51,36,10,18)
    Ast = (56,9,51,4,10,26)
    wrapped = [source_root(x) for x in (A,Ast)]
    states = [state(x,6) for x in wrapped]
    assert states == [(37,0,376,37),(37,0,352,37)]
    zgraphs = [Z(x,6) for x in wrapped]
    assert [c.hp_subsets(x)[-1] for x in zgraphs] == [413,389]
    # Positive nonunit control: R=C3, so H(R)=3. No mod3-unit hypothesis.
    contexts = [source_root(join(x,(2,4,1))) for x in (A,Ast)]
    context_states = [state(x,9) for x in contexts]
    assert context_states == [(111,0,4512,111),(111,0,4224,111)]
    assert tuple(x-y for x,y in zip(q(context_states[0]),q(context_states[1]))) == (0,288,288,288)
    carrier_states = [state(join((2,4,1),x),0) for x in (A,Ast)]
    assert carrier_states == [(111,0,1128,37),(111,0,1056,37)]
    oldcore = (18,8,3,5,14)
    newcores = ((50,8,3,5,46,14),(50,40,35,5,46,8))
    samples = []
    for a,b in ((0,1),(1,1),(1,2),(1,4),(2,2),(3,2)):
        old = [c.construct(oldcore,cut,a,b) for cut in (6,3)]
        new = [c.construct(core,cut,a,b) for core,cut in zip(newcores,(38,35))]
        # The new vertex is at index bit_length(a)+5, before the right flank.
        assert [delete(x,a.bit_length()+5) for x in new] == old
        so = [state(x,len(x)-1) for x in old]
        sn = [state(x,len(x)-1) for x in new]
        assert so[0][1:] == so[1][1:]
        assert sn[0][1:] == sn[1][1:]
        assert so[1][0]-so[0][0] == 3*a-(a&1)+2
        assert sn[1][0]-sn[0][0] == -3*a-(a&1)+2*(b&1)
        assert sn[0][1:] == (25,(147+304*a)*(1<<b.bit_length())-25*rev(b,b.bit_length()),25)
        samples.append({'A':a,'B':b,'old_states':so,'new_states':sn})
    loop_tests = 0
    for r in (5,629):
        for p in (0,1,4,79,102):
            n = p
            for k in range(1,5):
                n = (n << 14)+r
                assert loop_prediction(p,r,14,k) == g(n)
                loop_tests += 1
    receipt = {'scope':'Literal positive controls, not a coverage census; arbitrary-R and general-address conclusions proved algebraically in PROOF.md.', 'H37_wrapped_states':states,'H37_Z_rows':zgraphs,'nonunit_C3_context_states':context_states,'rooted_C3_carrier_states':carrier_states,'paired_birth_samples':samples,'exact_loop_formula_controls':loop_tests}
    (HERE/'positive_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__ == '__main__':
    main()
