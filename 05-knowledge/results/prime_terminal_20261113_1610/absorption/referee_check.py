"""Read-only independent checks of the repairs lane's induced ternary family."""
from itertools import combinations, product
from pathlib import Path
import json
import check as base
from cone_check import rows_of_code, hub_of

HERE = Path(__file__).resolve().parent
S = ((44,17,2,6,13,30),2)
T2 = ((110,92,24,112,97,70,4),6)
T4 = ((2,8,11,97,15,87,151,63),2)

def G(rows, root, subset_h=None):
    h = base.all_h(rows) if subset_h is None else subset_h
    full = (1 << len(rows))-1
    def visit(v, used):
        result = h[full ^ used]
        todo = rows[v] & ~used
        while todo:
            bit = todo & -todo
            todo -= bit
            result += visit(bit.bit_length()-1, used | bit)
        return result
    return visit(root, 1 << root)

def extension(rows, root):
    n = len(rows)
    out = [row | 1 << n | (0 if v == root else 1 << (n+1)) for v,row in enumerate(rows)]
    return tuple(out + [1 << (n+1), 1 << root])

def main():
    base.need(all(base.valid(rows) and base.strong(rows) and hub_of(rows, port) is not None
                  for rows,port in (S,T2,T4)), "valid strong two-run frames")
    embeddings = [(T2, [4,2,6,5,3,0]), (T4,[6,3,2,4,5,7])]
    for (target,port), mapping in embeddings:
        base.need(mapping[S[1]] == port, "module ports align")
        base.need(all(((S[0][i]>>j)&1) == ((target[mapping[i]]>>mapping[j])&1)
                      for i,j in combinations(range(6),2)), "induced old frame embedding")
    modules = [('TT1',base.tt(1)), ('TT2',base.tt(2)), ('TT3',base.tt(3)),
               ('C3',(2,4,1)), ('Q4',base.Q)]
    controls = []
    for name,module in modules:
        mh = base.all_h(module)
        h = mh[-1]
        p = sum(mh[x]*mh[(len(mh)-1)^x] for x in range(len(mh)))
        inflations = [base.inflate(frame,port,module) for frame,port in (S,T2,T4)]
        counts = [base.all_h(t)[-1] for t in inflations]
        base.need(counts == [12*p+h,36*p+5*h,36*p+7*h], "literal full ternary responses")
        # Inflation label 1 remains source vertex1 because the substituted port is2.
        g = G(inflations[0],1)
        base.need(g == 12*p+3*h, "direct path-enumeration marked response")
        base.need(base.all_h(extension(inflations[0],1))[-1] == counts[2], "two births give plus4 target")
        controls.append({"module": name, "h": h, "p": p, "full_counts": counts, "G_source_root1": g})
    generic = 0
    for n in range(1,5):
        for code in range(1 << (n*(n-1)//2)):
            rows = rows_of_code(n,code)
            h = base.all_h(rows)
            for root in range(n):
                base.need(base.all_h(extension(rows,root))[-1] == h[-1]+2*G(rows,root,h),
                          "generic two-birth lemma on all small labelled roots")
                generic += 1
    failure = []
    for name,(rows,port) in [('T2',T2),('T4',T4)]:
        h = base.all_h(rows)
        diffs = [G(rows,u,h)-h[-1] for u in range(len(rows))]
        base.need(2 not in diffs, "chosen child has no automatically rearmed root")
        failure.append({"name": name,"G_minus_H": diffs})
    alternate = []
    for rows, root in [((22,124,88,113,96,69,1),6), ((64,17,67,7,13,31,58),0)]:
        h = base.all_h(rows)[-1]
        g = G(rows,root)
        base.need(g-h == 2, "alternative finite head rearm")
        alternate.append({"rows": rows,"root": root,"H": h,"G": g})
    report = {"status": "ACCEPTED independent finite controls; proof review REFEREE.md",
              "generic_small_root_cases": generic,"literal_module_controls": controls,
              "failure_of_automatic_rearm": failure,"alternate_heads": alternate}
    # Audit every universal sign/congruence conclusion of the saved exact
    # two-port coefficient certificate, without using its category checker.
    incoming = json.loads((HERE.parent/'repairs/receipt.json').read_text())
    matrices = {tuple(r['ports']): {tuple(t['runs']):t['coefficient'] for t in r['delta']}
                for r in incoming['two_port_obstruction']['pairs']}
    negative = {(0,1),(0,2),(0,3),(0,5),(1,2),(1,3),(1,5),(2,3),(2,5),(3,5)}
    positive = {(4,6)}
    lower_y = {(2,4),(2,6),(3,4),(3,6),(5,6)}
    lower_x = {(4,5)}
    modular = {(0,4),(0,6),(1,4),(1,6)}
    base.need(set(matrices) == negative|positive|lower_y|lower_x|modular, "all21 pairs partitioned")
    for pair in negative|positive:
        sign = -1 if pair in negative else 1
        poly = matrices[pair]
        base.need(all(sign*c > 0 for c in poly.values()), "strict uniform coefficient sign")
        base.need(sign*poly[(1,2)] > 0 and sign*poly[(2,1)] > 0, "nonsingleton detector terms")
    for pair in lower_y|lower_x:
        transformed = {}
        axis = 1 if pair in lower_y else 0
        for exponents, coefficient in matrices[pair].items():
            if exponents[axis] == 2:
                base_exp = list(exponents); base_exp[axis] = 1; base_exp=tuple(base_exp)
                slack = list(exponents); slack[axis] = 0; slack=tuple(slack)
                transformed[base_exp] = transformed.get(base_exp,0)+2*coefficient
                transformed[slack] = transformed.get(slack,0)+coefficient
            else:
                transformed[exponents] = transformed.get(exponents,0)+coefficient
        base.need(all(c >= 0 for c in transformed.values()) and transformed[(1,1)] > 0,
                  "physical K2 lower bound turns mixed response positive")
        single = {e:c for e,c in matrices[pair].items() if e[axis] == 1}
        expected = {(2,1):-2} if axis == 1 else {(1,2):-2}
        base.need(single == expected, "singleton branch exact negative term")
    for pair in modular:
        for x1,x3,y1,y2,y3 in product((1,3),(0,2),(1,3),(0,2),(0,2)):
            x,y = {1:x1,2:2,3:x3},{1:y1,2:y2,3:y3}
            residue = sum(c*x.get(i,0)*y.get(j,0) for (i,j),c in matrices[pair].items()) % 4
            base.need(residue == 2, "physical factorial/2-adic obstruction")
        single = {j:c for (i,j),c in matrices[pair].items() if i == 1}
        base.need(single in ({2:22,3:12},{2:10}), "remaining singleton-X branch")
    report['two_port_universal_certificate'] = {"pairs":len(matrices),"sign_cases":11,"lower_bound_cases":6,"mod4_cases":4}
    bridge_a=(44,17,2,6,13,30,63)
    bridge_b=(110,16,2,6,13,30,62)
    old_core=(44,17,2,6,13,30)
    three_port=[]
    for length,module in [(1,base.tt(1)),(2,old_core),(3,base.Q)]:
        mh=base.all_h(module)
        h=mh[-1]
        p=sum(mh[x]*mh[(len(mh)-1)^x] for x in range(len(mh)))
        literal=[]
        for bridge in (bridge_a,bridge_b):
            first=base.inflate(bridge,1,base.tt(length))
            graph=base.inflate(first,6+length-1,module)
            literal.append(base.all_h(graph)[-1])
        expected=[(14*3**length-17)*h,2*(4*2**length-3)*p+(3*2**length-1)*h]
        base.need(literal==expected,"independent three-port chart counts")
        three_port.append({"b":length,"h":h,"p":p,"counts":literal})
    base.need(three_port[1]['counts']==[2725,8179] and all(8179%d for d in range(2,91)),
              "literal ternary prime8179 instance")
    report['three_port_controls']=three_port
    q=(6,12,8,1)
    bezout_frames=[((2,4,25,19,3),1),((2,12,73,81,103,79,3),5),((2,12,9,17,103,79,15),5)]
    bezout=[base.all_h(base.inflate(f,r,q))[-1] for f,r in bezout_frames]
    base.need(bezout==[199,599,601],"root frame Bezout instances")
    shifted_source=((16,5,9,35,46,7),0)
    shifted_target=((258,4,265,51,39,135,191,287,122),1)
    shifted=[]
    for k in range(3,6):
        a=base.all_h(base.inflate(*shifted_source,base.tt(k)))[-1]
        b=base.all_h(base.inflate(*shifted_target,base.tt(k-2)))[-1]
        base.need(a==27*2**k-17 and b==3*a+2,"literal shifted all-size chart controls")
        shifted.append({"k":k,"source":a,"target":b})
    report['root_frame_controls']={'bezout':bezout,'shifted':shifted}
    (HERE/'referee_receipt.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({"status":"ACCEPTED","generic_cases":generic,"controls":controls,"child_defects":failure}))

if __name__ == '__main__':
    main()
