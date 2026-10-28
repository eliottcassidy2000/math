"""Exact small controls for the information/code and Jacobian comparisons."""
from pathlib import Path
import json
import sympy as s
import reconstruct as r

HERE=Path(__file__).resolve().parent

def main():
    x,y,z=s.symbols('x y z');u=1+x*y
    F=s.Matrix([u**3*z+y*y*u*(4+3*x*y),
                y+3*x*u*u*z+3*x*y*y*(4+3*x*y),2*x-3*x*x*y-x**3*z])
    determinant=s.expand(F.jacobian((x,y,z)).det())
    points=[(0,0,-s.Rational(1,4)),(1,-s.Rational(3,2),s.Rational(13,2)),
            (-1,s.Rational(3,2),s.Rational(13,2))]
    images=[tuple(f.subs(dict(zip((x,y,z),p))) for f in F) for p in points]
    r.need(determinant==-2 and images==[(-s.Rational(1,4),0,0)]*3,'inherited Jacobian certificate')
    masks=[1072,4176]
    initial=[m^r.core.free_mask(r.core.sigma(r.core.presentation(7,m))) for m in masks]
    extended=[r.adjoin(r.core.presentation(7,m),True,False) for m in masks]
    syndrome=[r.core.free_mask(t)^r.core.free_mask(r.core.sigma(t)) for t in extended]
    r.need(all(syndrome) and syndrome[0]!=syndrome[1],'same blue zero-test, different syndrome')
    result=dict(status='FINITE-EXACT',jacobian=dict(determinant=str(determinant),
                points=[[str(v) for v in p] for p in points],images=[[str(v) for v in p] for p in images],
                source='THM-1300-jacobian-counterexample-dixmier-A3-explicit.md; independent verification, not discovery'),
                blue=dict(presentation_masks=masks,context='J_L',initial_syndromes=initial,linear_syndromes=syndrome,
                          zero_tests=[False,False],difference=masks[0]^masks[1]))
    (HERE/'cross-domain-receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result));return result

if __name__=='__main__':main()
