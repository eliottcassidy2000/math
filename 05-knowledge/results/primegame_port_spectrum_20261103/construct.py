"""Return an actual saved recursive-tree tournament for a requested H."""
import argparse,json
from pathlib import Path
from port_scout import add,h
def main():
    p=argparse.ArgumentParser();p.add_argument('target',type=int);args=p.parse_args()
    data=json.loads(Path(__file__).with_name('tree_scan11.json').read_text(encoding='utf-8-sig'))
    word=data['witnesses'].get(str(args.target))
    if word is None:raise SystemExit('No witness in this finite grammar bank; this is not a global nonexistence certificate.')
    a=[0]
    for x in word:a=add(a,x)
    actual=h(a)
    if actual!=args.target:raise RuntimeError('saved construction failed exact count')
    print(json.dumps({'target':args.target,'verified_H':actual,'order':len(a),'birth_ports':word,'outgoing_rows':a},indent=2))
if __name__=='__main__':main()
