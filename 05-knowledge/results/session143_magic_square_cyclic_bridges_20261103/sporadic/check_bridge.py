#!/usr/bin/env python3
"""Finite exact checks for the scoped C5 group-algebra bridge (Python stdlib).
No Baby Monster matrices are constructed; its character inputs are CITED.
Usage: python3 check_bridge.py [--sources-dir PATH] [--mutate NAME]
Mutation modes deliberately corrupt inputs and MUST fail, including under -O.
"""
import argparse, itertools, json, math, pathlib, re
from fractions import Fraction
from html.parser import HTMLParser

checks_executed=0
def require(condition, message):
    """Cannot be disabled by Python optimization flags."""
    global checks_executed
    checks_executed+=1
    if not condition:
        raise SystemExit("FAIL: "+message)

parser=argparse.ArgumentParser()
parser.add_argument('--sources-dir',type=pathlib.Path)
parser.add_argument('--mutate',choices=['bad-eigenvalue','bad-character','bad-transcription'])
args=parser.parse_args()
if args.mutate=='bad-transcription' and args.sources_dir is None:
    parser.error('bad-transcription requires --sources-dir')

def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a

def add(a,b):
    return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)

def rem(a,b,p=None):
    a=trim(a)
    while len(a)>=len(b) and a!=[0]:
        c=a[-1] if p else Fraction(a[-1],b[-1])
        if p:c=c*pow(b[-1],-1,p)%p
        for i,x in enumerate(b):a[i+len(a)-len(b)]-=c*x
        if p:a=[x%p for x in a]
        a=trim(a)
    return a

def cp(A):
    n=len(A);out=[0]
    for pi in itertools.permutations(range(n)):
        sign=(-1)**sum(pi[i]>pi[j] for i in range(n) for j in range(i+1,n))
        t=[sign]
        for i,j in enumerate(pi):t=mul(t,[-A[i][j],1] if i==j else [-A[i][j]])
        out=add(out,t)
    return out

def totient(n):return sum(math.gcd(k,n)==1 for k in range(1,n+1))

q=[1,3,4,2,1] # ascending coefficient convention
A=[[int((j-i)%5 in (1,2)) for j in range(5)] for i in range(5)]
require(all(A[i][i]==0 for i in range(5)), 'all(A[i][i]==0 for i in range(5))')
require(all(A[i][j]+A[j][i]==1 for i in range(5) for j in range(i)), 'all(A[i][j]+A[j][i]==1 for i in range(5) for j in range(i))')
require(all(sum(row)==2 for row in A), 'all(sum(row)==2 for row in A)')
require(cp(A)==mul([-2,1],q)==[-2,-5,-5,0,0,1], 'cp(A)==mul([-2,1],q)==[-2,-5,-5,0,0,1]')
phi5=[1]*5
alpha=[0,1] if args.mutate=='bad-eigenvalue' else [0,1,1]
qalpha=[0]
for c in reversed(q):qalpha=add(mul(qalpha,alpha),[c])
require(rem(qalpha,phi5)==[0], 'rem(qalpha,phi5)==[0]')
qmod=[x%2 for x in q]
# A reducible quartic over F2 has a degree1 or degree2 irreducible factor.
require(rem(qmod,[0,1],2)!=[0], 'rem(qmod,[0,1],2)!=[0]')
require(rem(qmod,[1,1],2)!=[0], 'rem(qmod,[1,1],2)!=[0]')
require(rem(qmod,[1,1,1],2)!=[0], 'rem(qmod,[1,1,1],2)!=[0]')
# Every nontrivial power alpha_r is another root, checked in Q[z]/Phi5.
for r in range(1,5):
    ar=[0]*(2*r+1);ar[r]+=1;ar[2*r]+=1
    value=[0]
    for c in reversed(q):value=rem(add(mul(value,ar),[c]),phi5)
    require(value==[0], 'value==[0]')
# Exact DFT inversion in the same quotient, with chi=[d,c,c,c,c].
rows=[]
for name,d,c,expected in [('C5 regular',5,0,[1,1,1,1,1]),('B 5A',4371,21,[891,870,870,870,870]),('B 5B',4371,-4,[871,875,875,875,875])]:
    if args.mutate=='bad-character' and name=='B 5A':c+=5
    multiplicities=[]
    for r in range(5):
        v=[d]
        for j in range(1,5):
            term=[0]*((-r*j)%5+1);term[-1]=c;v=add(v,term)
        v=rem(v,phi5)
        require(len(v)==1 and v[0]%5==0, 'len(v)==1 and v[0]%5==0')
        multiplicities.append(int(v[0]/5))
    require(multiplicities==expected, 'multiplicities==expected')
    require(sum(multiplicities)==d, 'sum(multiplicities)==d')
    require(2*multiplicities[0]-2*multiplicities[1]==2*c, '2*multiplicities[0]-2*multiplicities[1]==2*c')
    rows.append({'case':name,'degree':d,'nonidentity_character':c,'multiplicities':multiplicities,'operator_trace':2*c})
require(totient(5)==4 and totient(17)==16 and totient(85)==64, 'totient(5)==4 and totient(17)==16 and totient(85)==64')
# Genus formula for Gamma0(13): index14, two elliptic points of each kind, two cusps.
require(1+Fraction(14,12)-Fraction(2,4)-Fraction(2,3)-Fraction(2,2)==0, '1+Fraction(14,12)-Fraction(2,4)-Fraction(2,3)-Fraction(2,2)==0')
# Hostile control: nonrational one-dimensional character chi(g^j)=z^j.
# Its exact DFT is supported only at r=1, not equal on nontrivial roots.
hostile=[]
for r in range(5):
    v=[0]
    for j in range(5):
        term=[0]*(((1-r)*j)%5+1);term[-1]=1;v=add(v,term)
    v=rem(v,phi5)
    require(len(v)==1 and v[0]%5==0, 'len(v)==1 and v[0]%5==0')
    hostile.append(int(v[0]/5))
require(hostile==[0,1,0,0,0], 'hostile==[0,1,0,0,0]')
# Characteristic-five hostile: all five roots coalesce, invalidating separability.
require([v%5 for v in cp(A)]==[v%5 for v in [-32,80,-80,40,-10,1]], '[v%5 for v in cp(A)]==[v%5 for v in [-32,80,-80,40,-10,1]]')

source_check='not requested: imported Baby Monster data remain CITED'
if args.sources_dir:
    class Text(HTMLParser):
        def __init__(self):super().__init__();self.parts=[]
        def handle_data(self,d):self.parts.append(d)
    text=Text();text.feed((args.sources_dir/'ctblbm.html').read_text())
    s=''.join(text.parts)
    marker='gap> chi:= chi{ bfuslabels };'; tail=s[s.index(marker)+len(marker):]
    vector=re.search(r'\[([^]]+)\]',tail).group(1)
    chi=[None if x.strip()=='fail' else int(x) for x in vector.split(',')]
    if args.mutate=='bad-transcription':chi[17]+=5
    atlas=(args.sources_dir/'atlas_b.html').read_bytes().decode('latin1')
    classes=re.findall(r'<TD CLASS="ccl">(\d+[A-Z]+)<TD',atlas)
    require(len(chi)==len(classes)==184, 'len(chi)==len(classes)==184')
    require(classes[17:19]==['5A','5B'], "classes[17:19]==['5A','5B']")
    require(chi[17:19]==[21,-4], 'chi[17:19]==[21,-4]')
    require('Set( irr_atlas ) = Set( irr );' in s, "'Set( irr_atlas ) = Set( irr );' in s")
    # Both classes are absent from the complete exceptional Galois-class list.
    gal=s[s.index('gap> galoisinfo:= rec('):s.index('gap> galoisinfo.irrats:=')]
    require('"5A"' not in gal and '"5B"' not in gal, '\'"5A"\' not in gal and \'"5B"\' not in gal')
    require("# The 'i'-th class is rational." in s, '"# The \'i\'-th class is rational." in s')
    source_check='PASS: 184 aligned class positions, 5A/5B values, rational-class code'
print(json.dumps({'status':'PASS','checks_executed':checks_executed,'coefficient_order':'ascending','q':q,'regular_charpoly':cp(A),'cases':rows,'cyclotomic_degrees':[4,16,64],'source_transcription_check':source_check,'boundary':'Algebra checked; cited Baby Monster character table not independently reconstructed.'},indent=2))
