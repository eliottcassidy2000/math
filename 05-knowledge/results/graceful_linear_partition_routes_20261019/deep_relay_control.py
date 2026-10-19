"""Stress an iterative cycle tree deeper than Python's default recursion limit."""
from verify_partition_routes import route
import json

def main():
    r=1100
    classes=[tuple(range(r))]
    for j in range(1,r-1):
        excluded={0,r-1}|set(range(1,j-1))
        classes.append(tuple(x for x in range(r) if x not in excluded))
    classes[1]=tuple(range(1,r));classes.append((0,r-1))
    target=tuple(classes[i].index(i) for i in range(r));source=list(target)
    source[r-3]=classes[r-3].index(r-2);source[r-2]=classes[r-2].index(r-3)
    path=route(classes,tuple(source),target)
    if len(path)-1!=2*r-3:raise RuntimeError('long relay route count')
    print(json.dumps({'status':'PASS','rank':r,'route_length':len(path)-1,'arborescence_depth':r-2,'note':'Buffer-only routing stress; the full width-two distance here is one, so this is not a diameter lower bound.'},indent=2))

if __name__=='__main__':main()
