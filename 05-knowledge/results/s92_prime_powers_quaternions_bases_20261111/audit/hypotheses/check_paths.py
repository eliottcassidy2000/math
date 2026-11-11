#!/usr/bin/env python3
# Auditor "hypotheses": existence of every path named in the frontmatter (source / scripts / related) of the four
# canon files, and of the paths named in the session note. Read-only.
import os, re, glob, sys

ROOT = "C:/tmp/wt-s92"
FILES = [
    "01-canon/theorems/THM-5134-f1-local-densities-at-every-odd-modulus-prime-power-factor-and-the-four-square-singular-series.md",
    "01-canon/theorems/THM-5135-class-laws-of-bilinear-seams-table-powers-four-squares-e8-kernels-and-the-eight-extensions-of-f2-squared.md",
    "05-knowledge/hypotheses/HYP-9665-f1-densities-beyond-odd-moduli-even-moduli-single-values-primes-of-the-base.md",
    "05-knowledge/hypotheses/HYP-9666-seams-of-rank-four-and-the-octonion-law-equidistributed-families-rank-three-kernels-degenerate-primes.md",
    "05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/README.md",
]
PAT = re.compile(r"[A-Za-z0-9_./*\-]+\.(?:md|py|out|json|txt|c|npz)|(?:01-canon|05-knowledge|04-computation)/[A-Za-z0-9_./*\-]+")
EXPDIR = "04-computation/experiments/s92_prime_powers_20261111"


def resolve(tok, base):
    cands = []
    t = tok.strip("`.,;:()")
    if "*" in t:
        hits = glob.glob(os.path.join(ROOT, t + ("" if t.endswith(".md") else "*")))
        if not hits and "/" not in t:
            hits = glob.glob(os.path.join(ROOT, "**", t + "*"), recursive=True)
        return ("glob", hits)
    for prefix in ("", EXPDIR + "/", os.path.dirname(base) + "/", "05-knowledge/results/", "01-canon/theorems/", "05-knowledge/hypotheses/"):
        p = os.path.join(ROOT, prefix + t)
        if os.path.exists(p):
            return ("ok", [prefix + t])
    return ("missing", [])


def main():
    print("# check_paths.py (auditor hypotheses)")
    for f in FILES:
        text = open(os.path.join(ROOT, f), encoding="utf-8").read()
        if f.endswith("README.md"):
            head = text
            part = "whole note"
        else:
            # frontmatter only
            m = re.match(r"---\n(.*?)\n---\n", text, re.S)
            head = m.group(1)
            # skip title and status blocks (prose): keep from 'session:' on
            i = head.find("\nsession:")
            head = head[i:]
            part = "frontmatter (session .. related)"
        toks = sorted(set(PAT.findall(head)))
        print("== %s  [%s]: %d path-like tokens" % (os.path.basename(f)[:40], part, len(toks)))
        for t in toks:
            kind, hits = resolve(t, f)
            if kind == "missing":
                # try a repo-wide search by basename
                bn = os.path.basename(t.strip("`.,;:()"))
                found = glob.glob(os.path.join(ROOT, "**", bn), recursive=True) if bn else []
                found = [x.replace("\\", "/").replace(ROOT + "/", "") for x in found][:3]
                print("   MISSING as written: %s   %s" % (t, ("(basename found at: %s)" % found) if found else "(basename not found anywhere)"))
            elif kind == "glob":
                if not hits:
                    print("   GLOB WITHOUT MATCH: %s" % t)
                else:
                    print("   ok (glob, %d match) %s" % (len(hits), t))
            else:
                print("   ok %s" % t)


if __name__ == "__main__":
    main()
