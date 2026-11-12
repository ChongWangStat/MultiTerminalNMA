#!/usr/bin/env python3
"""Word count for the Research in Brief main text (Sections 1 to Discussion),
excluding figure environments (captions), display equations, citations and
references. Inline math counts as one word. Also reports caption words."""
import re, sys
s = open(sys.argv[1] if len(sys.argv) > 1 else "rib_body.tex", encoding="utf-8").read()
s = re.sub(r"(?m)^%.*$", "", s)
figs = re.findall(r"\\begin\{figure\}.*?\\end\{figure\}", s, flags=re.S)
def words(t):
    t = re.sub(r"\\begin\{(equation|align\*?)\}.*?\\end\{\1\}", " ", t, flags=re.S)
    t = re.sub(r"\\cite\{[^}]*\}", "", t)
    t = re.sub(r"\$[^$]*\$", " X ", t)
    t = re.sub(r"\\(label|ref|eqref)\{[^}]*\}", "1", t)
    t = re.sub(r"\\(section|paragraph)\*?\{([^}]*)\}", r" \2 ", t)
    t = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", t)
    t = re.sub(r"[{}~]", " ", t)
    return len(re.findall(r"[A-Za-z0-9][A-Za-z0-9'\u2019\-:.]*", t))
body = re.sub(r"\\begin\{figure\}.*?\\end\{figure\}", " ", s, flags=re.S)
heads = len(re.findall(r"\\(section|paragraph)\*?\{", body))
secs = re.split(r"\\section\{", body)
for sec in secs[1:]:
    print("%5d  %s" % (words("\\section{" + sec), sec[:sec.find("}")]))
caps = sum(words(c) for f in figs for c in re.findall(r"\\caption\{(.*)\}\\label|\\caption\{(.*)\}\s*\n\\label", f, flags=re.S) for c in c if c)
total = words(body)
print("MAIN TEXT (Sections 1-6, incl. headings, excl. captions/equations/references): %d words" % total)
print("figure captions: %d words; main text + captions: %d" % (caps, total + caps))
