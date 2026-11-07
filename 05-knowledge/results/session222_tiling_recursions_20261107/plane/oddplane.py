"""Exact, standard-library checks. Run python3 oddplane.py (also under -O)."""
from math import isqrt
import json
from pathlib import Path

def require(p, msg):
    if not p: raise RuntimeError(msg)

def address(n):
    if n < 1: raise ValueError('positive integers only')
    k=(n & -n).bit_length()-1
    return ((n >> k)-1)//2,k

def label(j,k): return (2*j+1)<<k

def shell(t):
    r=isqrt(t); s=t-r*r
    x,y=(r,s) if s<=r else (2*r-s,r)
    return (x,y) if r%2==0 else (y,x)

def rank(j,k):
    r=max(j,k)
    x,y=(j,k) if r%2==0 else (k,j)
    s=y if x==r else 2*r-x
    return r*r+s

def plane(t):
    r=isqrt(t); s=t-r*r; m=r//2
    if r%2==0:
        return (-m,m-s) if s<=r else (-m+s-r,-m)
    return (m+1,-m+s) if s<=r else (m+1-s+r,m+1)

def plane_rank(x,y):
    # Smallest r whose nested square [-floor(r/2),ceil(r/2)]^2 contains point.
    r=max(2*max(-x,-y,0),2*max(x,y,0)-1)
    m=r//2
    if r%2==0: s=m-y if x==-m else r+x+m
    else: s=y+m if x==m+1 else r+m+1-x
    return r*r+s

def main():
    checks=0
    for n in range(1,100001):
        j,k=address(n)
        require(label(j,k)==n,'dyadic inverse'); checks+=1
        require(address(2*n)==(j,k+1),'doubling'); checks+=1
        require(address(n+(1<<(k+1)))==(j+1,k),'scaled odd-core successor'); checks+=1
        require(shell(rank(j,k))==(j,k),'address shell inverse'); checks+=1
    for t in range(40000):
        j,k=shell(t)
        require(rank(j,k)==t,'shell rank'); checks+=1
        require(plane_rank(*plane(t))==t,'plane inverse'); checks+=1
        if t:
            require(sum(abs(a-b) for a,b in zip(shell(t),shell(t-1)))==1,'quadrant adjacency'); checks+=1
            require(sum(abs(a-b) for a,b in zip(plane(t),plane(t-1)))==1,'plane adjacency'); checks+=1
    for r in range(100):
        require({shell(t) for t in range((r+1)**2)}=={(j,k) for j in range(r+1) for k in range(r+1)},'quadrant coverage'); checks+=1
        a=-(r//2); b=(r+1)//2
        require({plane(t) for t in range((r+1)**2)}=={(x,y) for x in range(a,b+1) for y in range(a,b+1)},'plane coverage'); checks+=1
    hostile=[1,3,6,4,2,1]
    require(all(abs(a-b)==2 or a==2*b or b==2*a for a,b in zip(hostile,hostile[1:])), '5-cycle hostile')
    require(address(2)==(0,1) and address(3)==(1,0),'natural successor nonadjacency')
    out={'status':'PASS','checks':checks+2,'integer_universe':'1..100000','traversal_universe':'t=0..39999','complete_square_universe':'r=0..99','first_25_dyadic_path_labels':[label(*shell(t)) for t in range(25)],'literal_plus2_times2_odd_cycle':hostile}
    Path(__file__).with_name('checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
    # A labelled chart of 49 cells; path is the first seven complete shells.
    size=7; cell=80; margin=45
    svg=['<svg xmlns="http://www.w3.org/2000/svg" width="650" height="670" viewBox="0 0 650 670">','<rect width="650" height="670" fill="white"/>','<text x="35" y="28" font-family="sans-serif" font-size="18">Odd-core × dyadic-depth grid: nearest-neighbor shell path</text>']
    for j in range(size):
        for k in range(size):
            x=margin+j*cell;y=margin+(size-1-k)*cell
            svg.append(f'<rect x="{x}" y="{y}" width="80" height="80" fill="hsl({max(j,k)*40},70%,92%)" stroke="#ddd"/>')
    pts=' '.join(f'{margin+j*cell+40},{margin+(size-1-k)*cell+40}' for j,k in map(shell,range(size*size)))
    svg.append(f'<polyline points="{pts}" fill="none" stroke="#157080" stroke-width="3"/>')
    for j in range(size):
        for k in range(size):
            x=margin+j*cell+40;y=margin+(size-1-k)*cell+40
            svg.append(f'<circle cx="{x}" cy="{y}" r="22" fill="white"/><text x="{x}" y="{y+5}" text-anchor="middle" font-family="sans-serif" font-size="15">{label(j,k)}</text>')
    svg+=['<text x="45" y="635" font-family="sans-serif" font-size="15">Right: odd core +2 (label +2^(k+1)); up: label ×2.</text>','<text x="45" y="658" font-family="sans-serif" font-size="15">Numbers are labels, not visit times. Shell sizes: 1, 3, 5, 7, ...</text>','</svg>']
    Path(__file__).with_name('oddplane.svg').write_text('\n'.join(svg))
if __name__=='__main__': main()
