"""Wave 58 §1 checks.  (1) Hub-arc flip: H(u2->u1) - H(u1->u2) = 2(K(x,y) - K(y,x)) on every word pair, spine <= 8.
(2) Neutrality: a shared position i (x_i = y_i = g) contributes g*A_i + (1-g)*B_i to K(x,y) - K(y,x), where
    A_i = sum_{live j >= i+2} (x_j - y_j) 2^(j-i-2) and B_i = sum_{live j <= i-2} (x_j - y_j) 2^(i-j-2); so it is a free
    neutral choice iff A_i = B_i.  Also: the largest number of free positions in one maximal run of shared positions."""
import itertools
exec(open('w58_frames.py').read().split('def frames(')[0])      # count_order, H
def K(x, y):
    m = len(x)
    return sum(1 for i in range(m) if x[i] == 1 and y[i] == 0) + sum(2 ** (j - i - 2) for i in range(m) for j in range(i + 2, m) if x[i] == 1 and y[j] == 0)
bad = tot = 0
for m in range(1, 9):
    for x in itertools.product((0, 1), repeat=m):
        for y in itertools.product((0, 1), repeat=m):
            tot += 1
            if H(list(x), list(y), 0) - H(list(x), list(y), 1) != 2 * (K(x, y) - K(y, x)): bad += 1
print(f"(1) gap formula on {tot} word pairs (spine <= 8): mismatches {bad}")
viol = cases = 0; maxfree = 0; witness = None
from collections import Counter
runstat = Counter()
for m in range(1, 9):
    for x in itertools.product((0, 1), repeat=m):
        for y in itertools.product((0, 1), repeat=m):
            live = [j for j in range(m) if x[j] != y[j]]
            freeset = []
            for i in range(m):
                if x[i] != y[i]: continue
                A = sum((x[j] - y[j]) * 2 ** (j - i - 2) for j in live if j >= i + 2)
                B = sum((x[j] - y[j]) * 2 ** (i - j - 2) for j in live if j <= i - 2)
                x2 = list(x); y2 = list(y); x2[i] = y2[i] = 1 - x[i]
                same = (K(x, y) - K(y, x)) == (K(x2, y2) - K(y2, x2))
                cases += 1
                if same != (A == B): viol += 1
                if same: freeset.append(i)
            if live:
                # maximal runs of shared positions
                run = []; runs = []
                for i in range(m):
                    if x[i] == y[i]: run.append(i)
                    else:
                        if run: runs.append(run)
                        run = []
                if run: runs.append(run)
                for r in runs:
                    k = sum(1 for i in r if i in freeset); runstat[k] += 1
                    if k > maxfree: maxfree = k; witness = (x, y, r)
print(f"(2) shared-bit flips: {cases} cases; criterion 'free iff A_i = B_i' violated {viol} times")
print(f"    free positions per maximal shared run (live set nonempty): {dict(sorted(runstat.items()))}; max {maxfree}, e.g. {witness}")
