"""Exact minimum common polynomial recurrences for the complete joint observer."""
import json
import sympy as s
from return_probe import HERE,need

def evaluate(poly,M):
    R=s.zeros(M.rows);I=s.eye(M.rows)
    for c in poly.all_coeffs():R=R*M+c*I
    return R

def minimum(M,x):
    p=s.Poly(M.charpoly(x).as_expr(),x)
    for factor,multiplicity in s.factor_list(p)[1]:
        for _ in range(multiplicity):
            q=p.exquo(factor)
            if evaluate(q,M)==s.zeros(M.rows):p=q
            else:break
    need(evaluate(p,M)==s.zeros(M.rows),'exact matrix annihilator')
    return p

def main():
    d=json.loads((HERE/'split_algebra.json').read_text());x=s.Symbol('x');rows=[]
    for name,raw in d['diagonal_matrices'].items():
        M=s.Matrix(raw);blocks=[M[:14,:14],M[14:26,14:26],M[26:,26:]]
        ps=[minimum(B,x) for B in blocks];p=ps[0]
        for q in ps[1:]:p=s.lcm(p,q)
        need(evaluate(p,M)==s.zeros(31),'all-source all-query recurrence')
        rows.append({'operation':name,'block_minimal_polynomials':[str(s.factor(q.as_expr())) for q in ps],
                     'minimal_polynomial':str(s.factor(p.as_expr())),'descending_coefficients':[int(c) for c in p.all_coeffs()]})
    (HERE/'recurrences.json').write_text(json.dumps({'status':'EXACT_ALL_WORD_STATE_RECURRENCES','generators':rows},indent=2)+'\n',encoding='utf8',newline='\n')
    for r in rows:print(r['operation'],r['minimal_polynomial'],flush=True)
if __name__=='__main__':main()
