#!/usr/bin/env python3
"""Audit helper (audit_ssj): list id / published date / primary category / title of the entries
of an (untrusted) arXiv Atom API response.  Parsed as XML data only.
Usage: python3 -I atom_titles.py <file.xml> [keyword ...]   (keywords searched in title+abstract)
"""
import sys
import xml.etree.ElementTree as ET

NS = {"a": "http://www.w3.org/2005/Atom", "x": "http://arxiv.org/schemas/atom"}


def main():
    tree = ET.parse(sys.argv[1])
    kws = [k.lower() for k in sys.argv[2:]]
    entries = tree.getroot().findall("a:entry", NS)
    print(f"{len(entries)} entries")
    for e in entries:
        eid = (e.findtext("a:id", default="", namespaces=NS) or "").rsplit("/", 1)[-1]
        title = " ".join((e.findtext("a:title", default="", namespaces=NS) or "").split())
        pub = (e.findtext("a:published", default="", namespaces=NS) or "")[:10]
        prim = e.find("x:primary_category", NS)
        prim = prim.get("term") if prim is not None else "?"
        summ = " ".join((e.findtext("a:summary", default="", namespaces=NS) or "").split()).lower()
        hits = [k for k in kws if k in summ or k in title.lower()]
        print(f"  {eid:14s} {pub} {prim:10s} {title[:110]}" + (f"   [hits: {hits}]" if hits else ""))


if __name__ == "__main__":
    main()
