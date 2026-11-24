"""Replay the closed-cycle antipode counts and the converse-paired diamonds."""
from itertools import permutations, product
from math import gcd
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent


def totient(n):
    return sum(gcd(i,n)==1 for i in range(1,n+1))


def fixed_merged(n):
    if n % 2:
        return 0
    m = n//2
    numerator = sum(totient(d)*(1 << ((n//d)*(m-1)))
                    for d in range(1,n+1) if n%d==0 and d%4==2)
    assert numerator % (2*n) == 0
    return numerator//(2*n)


def heights(rows):
    n = len(rows)
    return sum(all(rows[a] >> b & 1 for a,b in zip(word,word[1:]))
               for word in permutations(range(n)))


def iso_counts(rows):
    n = len(rows)
    result = []
    for reverse in (0, 1):
        # Degrees restrict all possible relabelings. Count each exactly once.
        groups = [[i for i,r in enumerate(rows) if r.bit_count()==s] for s in range(n)]
        image = groups[::-1] if reverse else groups
        if list(map(len,groups)) != list(map(len,image)):
            result.append(0)
            continue
        count = 0
        source = sum(groups,[])
        for blocks in product(*(permutations(g) for g in image)):
            target = sum(map(list,blocks),[])
            if all((rows[source[i]] >> source[j] & 1) ==
                   ((rows[target[i]] >> target[j] & 1) ^ reverse)
                   for i in range(n) for j in range(i+1,n)):
                count += 1
        result.append(count)
    return tuple(result)


def flip(rows, i):
    rows = rows.copy()
    j = (i+1)%len(rows)
    rows[i] ^= 1 << j
    rows[j] ^= 1 << i
    return rows


def main():
    assert [fixed_merged(n) for n in (5,6,7,8,10)] == [0,6,0,256,52432]
    examples = [
        ([248,125,241,212,96,8,160,50],3,(173,175,177)),
        ([154,224,227,246,70,209,129,16],0,(217,219,221)),
        ([154,100,233,210,70,153,33,82],3,(541,539,537)),
    ]
    for rows, i, (a,b,c) in examples:
        j = (i+4)%8
        graphs = [rows,flip(rows,i),flip(rows,j),flip(flip(rows,i),j)]
        assert [heights(g) for g in graphs] == [a,b,b,c]
        assert [iso_counts(g) for g in graphs] == [(1,1),(1,0),(1,0),(1,1)]
        print("diamond",a,b,c,"middle merged mass",2*b)
    expected = {
        4: ([3],[1]),
        6: ([23,29,41,43],[1,29,41]),
        8: ([117,197,245,257,369,409,415,417,497,507,519,521,531,545,641,649,651,657,659],
            [1,29,41,173,217,537]),
    }
    with tempfile.TemporaryDirectory() as temp:
        exe = Path(temp)/"rim_cube"
        subprocess.run(["c++","-O3","-std=c++17",str(HERE/"rim_cube.cpp"),"-o",str(exe)],check=True)
        for n, rows in expected.items():
            output = subprocess.check_output([str(exe),str(n)],text=True).splitlines()
            actual = tuple(list(map(int,line.split())) for line in output)
            assert actual == rows, (n,actual)
            print("order",n,"unit lines",actual[0],"balanced diamonds",actual[1])


if __name__ == "__main__":
    main()
