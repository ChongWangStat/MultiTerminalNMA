#!/usr/bin/env python3
"""Merge the two updated LaTeX PDFs preserving intra-document citation/equation jumps.

Standalone TeX links can be named destinations; pdfunite silently drops their
resolved page targets. Recreate links as page-specific GoTo annotations.
"""
import sys
from pathlib import Path
import fitz
if len(sys.argv)!=4: raise SystemExit('usage: merge_preserve_pdf_links.py main.pdf supp.pdf combined.pdf')
inputs=[fitz.open(sys.argv[1]),fitz.open(sys.argv[2])];out=fitz.open(); offset=0
for doc in inputs:
    out.insert_pdf(doc,links=False)
    for i,src in enumerate(doc):
        dst=out[offset+i]
        for link in src.get_links():
            kind=link['kind']
            if kind in (fitz.LINK_GOTO, fitz.LINK_NAMED) and 'page' in link and 0<=link['page']<len(doc):
                dst.insert_link({'kind':fitz.LINK_GOTO,'from':link['from'], 'page':offset+link['page'], 'to':link.get('to',fitz.Point(0,0))})
            elif kind==fitz.LINK_URI:
                dst.insert_link({'kind':fitz.LINK_URI,'from':link['from'],'uri':link['uri']})
            elif kind==fitz.LINK_LAUNCH:
                raise RuntimeError('Unresolved file-launch link from '+str(i))
            elif kind==fitz.LINK_NAMED and 'page' not in link:
                raise RuntimeError('Unresolved named destination on page '+str(i))
    offset+=len(doc)
out.set_toc([[1,'Main manuscript',1],[1,'Supplementary Material',len(inputs[0])+1]])
out.set_metadata({'title':'Multi-terminal study resistances in network meta-analysis - Main and Supplement','author':'Chong Wang, Zhezhen Jin, Annette M. O\'Connor','subject':'Research Synthesis Methods - Research in Brief'})
out.save(sys.argv[3],garbage=4,deflate=True)
check=fitz.open(sys.argv[3]); links=[x for pg in check for x in pg.get_links()]
assert any(x['kind']==fitz.LINK_GOTO for x in links),'Internal links missing after merge'
assert not any(x['kind']==fitz.LINK_NAMED for x in links),'Unresolved named destinations'
assert not any(x['kind']==fitz.LINK_GOTO and not 0<=x['page']<len(check) for x in links),'Invalid goto'
print('Combined',len(check),'pages;',sum(x['kind']==fitz.LINK_GOTO for x in links),'internal links;',sum(x['kind']==fitz.LINK_URI for x in links),'URL links')
