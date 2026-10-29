"""Exact faithful unit-distance theta(3,3,3)+2K1 and geometric hostiles."""
from pathlib import Path
from itertools import combinations
import json
import sympy as s
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(b,msg):
    global CHECKS
    CHECKS+=1
    if not b:raise RuntimeError(msg)

def distance2(a,b):
    return s.simplify(sum((x-y)**2 for x,y in zip(a,b)))

def main():
    points=[(s.Integer(0),s.Integer(0)),(s.Integer(2),s.Integer(0))];edges=set()
    for t in (s.Rational(1,4),s.Rational(1,3),s.Rational(1,2)):
        a=((1-t*t)/(1+t*t),2*t/(1+t*t));w=points[1]
        dx,dy=w[0]-a[0],w[1]-a[1];d2=dx*dx+dy*dy
        h=s.sqrt((4-d2)/(4*d2))
        b=(s.simplify((a[0]+w[0])/2-h*dy),s.simplify((a[1]+w[1])/2+h*dx))
        ia=len(points);points.extend([a,b])
        edges.update({tuple(sorted(e)) for e in ((0,ia),(ia,ia+1),(ia+1,1))})
    points.extend([(s.Integer(10),s.Integer(10)),(s.Integer(20),s.Integer(20))])
    separation=[]
    for i,j in combinations(range(len(points)),2):
        d2=distance2(points[i],points[j]);eq=(s.simplify(d2-1)==0)
        need(d2.is_zero is False,'distinct exact algebraic points')
        need(eq==((i,j) in edges),'faithful unit-distance relation')
        if not eq:
            need(s.simplify(d2-1).is_zero is False,'certified nonedge distance is not one')
            separation.append(abs(float(d2-1)))
    # An exact rational degree-seven star, refuting an unrestricted degree-six bound.
    circle=[]
    for t in map(s.Rational,range(7)):
        p=((1-t*t)/(1+t*t),2*t/(1+t*t));circle.append(p)
        need(distance2((0,0),p)==1,'rational unit neighbor')
    need(len(set(circle))==7,'seven distinct rational neighbors')
    # A rigid oriented excess-minus-one carrier containing K_(2,3).
    # Its eight distinct (outdegree,indegree) pairs prove trivial automorphism group.
    arcs=[]
    for z,pattern in zip((2,3,4),((0,0),(1,0),(1,1))):
        for centre,bit in enumerate(pattern):arcs.append((centre,z) if bit else (z,centre))
    arcs.append((5,6))
    degrees=[(sum(a==v for a,b in arcs),sum(b==v for a,b in arcs)) for v in range(8)]
    need(len(arcs)==7 and len(set(degrees))==8,'rigid E=V-1 geometric obstruction')
    report={'status':'FINITE-EXACT','geometry':'faithful Euclidean unit distances in R2',
       'theta_paths':3,'theta_path_length':3,'isolates':2,'vertices':len(points),'edges':sorted(edges),
       'coordinates':[[str(x),str(y)] for x,y in points],
       'nonedge_squared_distance_gap_float_diagnostic':min(separation),
       'degree_seven_control':[[str(x),str(y)] for x,y in circle],
       'rigid_nonunit_carrier':{'vertices':8,'arcs':arcs,'degree_pairs':degrees,
           'reason_no_unit_realization':'vertices 0 and 1 have three distinct common neighbors 2,3,4'},
       'checks':CHECKS,'sympy_version':s.__version__}
    need(len(edges)==len(points)-1,'sparse carrier Euler characteristic one')
    report['checks']=CHECKS
    (HERE/'unit_distance.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report))
if __name__=='__main__':main()
