#!/usr/bin/env python3
# mac-mini-2026-11-08-firstrealizers: refinements of the first-creator fibre (least order n(h)) and whether the owner's
# dichotomy 'one self-converse class or exactly one converse pair' then holds. Input: directory with cls_n.txt (n <= 9).
# Output: refine_dichotomy.out
import sys, os
D = sys.argv[1] if len(sys.argv) > 1 else '.'
import collections
fibers = {}           # h -> list of (n, code, ccode, SC, aut, c3, strong, scores?)
seen_values = set()
for n in range(1, 10):
    cur = collections.defaultdict(list)
    with open(os.path.join(D, f"cls_{n}.txt")) as f:
        for line in f:
            head = line.split("|")[0].split()
            code, H, c3, strong, SC, aut, cc = int(head[0]), int(head[1]), int(head[2]), int(head[3]), int(head[4]), int(head[5]), int(head[6])
            if H not in seen_values:
                cur[H].append(dict(n=n, code=code, cc=cc, SC=SC, aut=aut, c3=c3, strong=strong))
    for H, lst in cur.items():
        fibers[H] = lst
    seen_values |= set(cur)
def dichotomy(lst):
    # one SC class, or exactly two classes that are converses of each other
    if len(lst) == 1: return lst[0]["SC"] == 1
    if len(lst) == 2: return lst[0]["SC"] == 0 and lst[0]["cc"] == lst[1]["code"]
    return False
vals = sorted(h for h in fibers if h >= 3)
print("values with birth order <= 9:", len(vals))
def test(name, key):
    ok = 0
    for h in vals:
        lst = fibers[h]
        best = key(lst)
        sub = [c for c in lst if c in best]
        ok += dichotomy(sub)
    print(f"{name:42s}: dichotomy holds for {ok}/{len(vals)} values")
test("no refinement", lambda L: L)
test("max |Aut|", lambda L: [c for c in L if c["aut"] == max(x["aut"] for x in L)])
test("min |Aut|", lambda L: [c for c in L if c["aut"] == min(x["aut"] for x in L)])
test("min c3", lambda L: [c for c in L if c["c3"] == min(x["c3"] for x in L)])
test("max c3", lambda L: [c for c in L if c["c3"] == max(x["c3"] for x in L)])
test("max |Aut| then max c3", lambda L: (lambda A: [c for c in A if c["c3"] == max(x["c3"] for x in A)])([c for c in L if c["aut"] == max(x["aut"] for x in L)]))
test("max |Aut| then min c3", lambda L: (lambda A: [c for c in A if c["c3"] == min(x["c3"] for x in A)])([c for c in L if c["aut"] == max(x["aut"] for x in L)]))
# size distribution of the max-|Aut| refined fiber
dist = collections.Counter()
for h in vals:
    L = fibers[h]; A = [c for c in L if c["aut"] == max(x["aut"] for x in L)]
    dist[(sum(c['SC'] for c in A), (len(A) - sum(c['SC'] for c in A)) // 2)] += 1
print("max|Aut| refined (s,p) distribution (top 10):", dist.most_common(10))
