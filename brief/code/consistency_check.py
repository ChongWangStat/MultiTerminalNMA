#!/usr/bin/env python3
"""Pre-submission check of the Research in Brief package: RSM format limits,
numbers quoted versus recomputed outputs, cross-references between main text
and supplement, and terminology/spelling consistency. Run from rib/ after
build_all.sh."""
import json, re, csv, subprocess, sys, os
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
main = open('main_rib.tex').read(); body = open('rib_body.tex').read(); supp = open('supp_rib_body.tex').read()
fails = []
def must(cond, msg):
    print(("OK  " if cond else "FAIL") + " " + msg)
    if not cond: fails.append(msg)
def nwords(t):
    t = re.sub(r"\$[^$]*\$", " X ", t); t = re.sub(r"\\cite\{[^}]*\}", "", t)
    t = re.sub(r"\\[a-zA-Z]+\*?", " ", t); t = re.sub(r"[{}~]", " ", t)
    return len(re.findall(r"[A-Za-z0-9][A-Za-z0-9'\-:.]*", t.replace("--", " ")))
# ---- RSM Research in Brief format
wc = subprocess.run([sys.executable, "code/wordcount.py", "rib_body.tex"], capture_output=True, text=True).stdout
n = int(re.search(r"MAIN TEXT.*?: (\d+) words", wc).group(1))
must(n <= 1500, "main text <= 1,500 words excluding references (%d)" % n)
ab = re.search(r"\\abstract\{(.*?)\n\}", main, re.S).group(1)
must(nwords(ab) <= 250, "abstract <= 250 words (%d)" % nwords(ab))
must("\\cite" not in ab and not re.search(r"\(\d{4}\)", ab), "abstract cites no works")
kw = re.search(r"\\keywords\{(.*?)\}", main, re.S).group(1).split(",")
must(4 <= len(kw) <= 6, "4-6 keywords (%d)" % len(kw))
for h in ["What is already known?", "What is new?", "Potential impact for RSM readers"]:
    must(h in main, "Highlights heading '%s'" % h)
must(main.find("\\section*{Highlights}") > main.find("\\end{Frontmatter}") and main.find("\\section*{Highlights}") < main.find("\\input{rib_body.tex}"), "Highlights after abstract, before main text")
for h in ["Author contributions", "Competing interest statement", "Data availability statement", "Funding statement", "Ethical standards", "Supplementary material"]:
    must("\\paragraph{%s}" % h in main, "back matter: %s" % h)
must("The authors declare none." in main, "competing interests wording")
must("\\articletype{Research-in-Brief}" in main and "journal=rsm" in main, "RSM template, article type Research-in-Brief")
bbl = open('main_rib.bbl').read(); nref = len(re.findall(r"\\entry\{", bbl))
must(10 <= nref <= 25, "reference count within precedent range 10-22 (+3) (%d)" % nref)
must("Kansas State University" in main and "Michigan State University},\n\\orgaddress" not in main, "O'Connor affiliation Kansas State")
# ---- numbers vs outputs
d = json.load(open('data/depression_diagnostic_output.json'))
for key, txt in [("fitted_NMA_common_effect_1v3_published_orientation", "-0.119725"), ("theta_sub_shortest_path", "-0.035099"), ("theta_sub_random_walk", "-0.032287")]:
    must(abs(d[key] - float(txt)) < 5e-7 and txt in body and txt in supp, "depression %s in main and supplement" % txt)
must(abs(d["discrepancy_shortest"] - 0.084626) < 5e-7 and abs(d["discrepancy_random_walk"] - 0.087438) < 5e-7 and "0.084626" in supp, "depression discrepancies")
must(d["max_abs_diff_hat_row_vs_published_flow"] < 1.45e-8 and "1.44\\times10^{-8}" in supp, "hat row vs published flow")
out = subprocess.run([sys.executable, "-I", "code/covid_currents.py", "figures/data/covid_AE_study_polygon_currents.csv"], capture_output=True, text=True).stdout
must(float(re.search(r"= ([0-9.e-]+)", out).group(1)) < 1e-15, "A:E displayed currents reproduce from arm counts")
be = re.search(r"\('B', 'E'\) (.*)", out).group(1)
for v in ["0.578807", "0.057231", "0.031858", "0.133519", "0.534377"]:
    must(v in be and v in supp, "B:E flow %s" % v)
must("-0.046109" in out and "-0.046109" in supp and "0.004829" in supp, "A:E fitted value and variance")
for fn, key in [('audit_output/covid_netmeta_audit_results.csv', 'S:tab:netmeta}'), ('audit_output/covid_netmeta_audit_results_labelorientation.csv', 'S:tab:netmeta_u}')]:
    for r in csv.DictReader(open(fn)):
        f = float(r["Fitted"]); line = "%s & %s & " % (r["Target"], ("$%.6f$" if f < 0 else "%.6f") % f)
        vals = " & ".join(("$%.6f$" if float(r[m]) < 0 else "%.6f") % float(r[m]) for m in ["shortestpath", "randomwalk", "cccp", "pseudoinverse"])
        must(line + vals in supp, "netmeta row %s (%s)" % (r["Target"], key[:-1]))
rows = [l for l in open('data/depression_percomparison_rows.tex').read().strip().split("\n")]
def fmt_row(l):
    parts = [p.strip() for p in l.replace("\\\\", "").split("&")]
    return " & ".join([parts[0]] + [("$%s$" % p if p.startswith("-") else p) for p in parts[1:]]) + " \\\\"
must(all(fmt_row(l) in supp for l in rows), "depression per-comparison table equals script output")
# counts quoted in main text equal supplement
# v18: main-text counts are for the terminal-current graph and the aggregate hat-matrix flow graph
ra = list(csv.DictReader(open('audit_output/route_audit_terminal_covid.csv')))
n_study = sum(r["study_unique"] == "True" for r in ra); n_trt = sum(r["trt_unique"] == "True" for r in ra)
must(len(ra) == 10 and n_study == 0 and "all ten targets have non-unique study-level routes" in body and "No target has unique study-level routes" in supp, "terminal-current audit: 0 of 10 unique at study level")
must(n_trt == 2 and "eight do so after treatment-level reduction" in body and "two of ten have unique treatment-level paths" in supp, "aggregate-flow audit: 8 of 10 non-unique at treatment level")
for r in ra:
    row = "%s & %s & %s & %s & %s & %s & %s\\\\" % (r["target"], r["study_routes"], r["study_dim"], "yes" if r["study_unique"] == "True" else "no", r["trt_routes"], r["trt_dim"], "yes" if r["trt_unique"] == "True" else "no")
    must(row in supp, "terminal audit table row %s" % r["target"])
# New independent B:E route-bound audit, separate from the CSP pairing convention.
routes = list(csv.DictReader(open('audit_output/route_bounds_terminal_BE.csv')))
keyed = {r['route']: r for r in routes}
markov_sum = sum(float(r['markov']) for r in routes)
for r in routes:
    must(float(r['lower']) - 1e-10 <= float(r['markov']) <= float(r['upper']) + 1e-10,
         'B:E proportional route falls in sharp interval: ' + r['route'])
must(len(routes) == 8 and abs(markov_sum - 1) < 1e-9 and
     sum(abs(float(r['upper']) - float(r['lower'])) > 1e-9 for r in routes) == 6,
     'B:E terminal route audit: 8 routes, 6 nonunique, Markov mass 1')
rr = keyed.get('B-Rec-A-Act-E', {})
must(rr and abs(float(rr['markov']) - 0.115709) < 1e-6 and
     abs(float(rr['lower']) - 0.044430) < 1e-6 and
     abs(float(rr['upper']) - 0.133519) < 1e-6 and
     '0.1157' in body and '0.0444' in body and '0.1335' in body,
     'B:E illustrative route interval and proportional selection checked')
must('\\label{S:tab:terminalroutes}' in supp and
     'audit_output/route_bounds_terminal_BE_rows.tex' in supp,
     'New B:E bounds table referenced and generated from code')
# CSP edge-graph audits stay in the supplement only
for b in ["four targets have unique study provenance", "nine of ten appear unique", "8/45 & 8/45", "33/45 & 33/45"]:
    must(b in supp, "CSP edge-graph count in supplement: %s" % b)
co = open('audit_output/correlated_arm_examples.txt').read()
must("1.5866" in co and "1.5866" in supp and "1.59 units" in body, "correlated-arm example: 1.59 units at terminal B (B:D)")
must("'A:B': 0.3644" in co and "'B:D': 14.0186" in co and "$14.019$ (B:D)" in supp, "correlated-arm example reproduces published network variances")
must("1.0042" in co and "$+1.0042$" in supp and "0.250434" in co, "real correlated-arm example values")
# ---- cross-references
ma = open('main_rib.aux').read(); sa = open('supplement_rib.aux').read()
meq = dict(re.findall(r"newlabel\{(eq:[^}]*)\}\{\{([0-9]+)\}", ma))
expect = {"eq:bridge": "1", "eq:twostep": "2", "eq:element": "3", "eq:studycurrent": "4", "eq:poly": "5", "eq:family": "6"}
must(meq == expect, "main equation numbers %s" % meq)
for name, num in [("current", "1"), ("twostep", "3"), ("poly", "4"), ("family", "5"), ("element", "6"), ("studycurrent", "7")]:
    pass
cited = set(re.findall(r"main Equation~\((\d)\)", supp))
must(cited <= set(expect.values()), "supplement cites existing main equations %s" % sorted(cited))
ctx = {"1": "GBL_B^+(e_a-e_b)]_{uv}", "3": "Q_{jl}=w_{jl}/d_j^B", "4": "A\\pi=f\\}$ (main Equation~(4))", "5": "main Equation~(5)", "6": "\\one\\T i_k=0$ (main Equation~(6))", "7": "L^+(e_a-e_b)$ (main Equation~(7))"}
must("for $j\\ne l$ (main Equation~(2))" in supp and "A\\pi=f\\}$ (main Equation~(5))" in supp and supp.count("main Equation~(6)") == 2 and "\\one\\T i_k=0$ (main Equation~(3))" in supp and "L^+(e_a-e_b)$ (main Equation~(4))" in supp, "supplement equation citations point to the right main equations")
sl = dict(re.findall(r"newlabel\{(S:[^}]*)\}\{\{([^}]*)\}", sa))
for sec, num in [("S:bridge", "S1"), ("S:routes", "S2"), ("S:element", "S3"), ("S:audits", "S4"), ("S:recon", "S5"), ("S:prov", "S6"), ("S:fig:covid", "S1")]:
    must(sl.get(sec) == num, "supplement %s = %s" % (sec, num))
mains = re.findall(r"Supplement(?:, | )Sections?~(S\d)(?: and (S\d))?|Supplement, Figure~(S\d)", body)
must(all(any(x in ("S1", "S2", "S3", "S4", "S5") for x in m if x) for m in mains) and len(mains) >= 5, "main-text supplement citations resolve (%s)" % mains)
must("Supplement Section~S1." in body, "S1 cited for the proofs of the conjectures")
must("Supplement Sections~S2 and S4" in body and "Supplement, Section~S3)" in body and "Supplement Section~S3)" in body and body.count("Supplement, Section~S5)") == 1, "S2/S3/S4/S5 citations")
order = []
for m in re.finditer(r"\\ref\{(fig:[a-z]+)\}", body):
    if m.group(1) not in order: order.append(m.group(1))
must(order == ["fig:circuit", "fig:flows", "fig:covid"], "figures first cited in order 1, 2, 3 (%s)" % order)
labs = set(re.findall(r"\\label\{([^}]*)\}", body)); refs = set(re.findall(r"\\ref\{([^}]*)\}", body))
must(refs <= labs, "all main \\ref targets defined")
labs = set(re.findall(r"\\label\{([^}]*)\}", supp)); refs = set(re.findall(r"\\ref\{([^}]*)\}", supp))
must(refs <= labs, "all supplement \\ref targets defined")
# ---- terminology, spelling, acronyms
alltxt = main + body + supp
uk = re.findall(r"\b(labelled|labelling|behaviour|centre|modelling|analysed|neighbour|colour|realisation|normalise)\b", alltxt)
must(not uk, "US spelling throughout %s" % set(uk))
b = re.sub(r"\\cite\{[^}]*\}", "", body)
for ac, full in [("NMA", "Network meta-analysis (NMA)"), ("CSP", "contrast-space projection (CSP)"), ("GLS", "generalized least-squares (GLS)"), ("COVID-19", "coronavirus disease 2019 (COVID-19)")]:
    i = b.find(ac); j = b.find(full)
    must(j >= 0 and j <= i, "%s defined at first use in main text" % ac)
must("PASI-75" not in body, "PASI-75 not used in main text (defined in supplement)")
must("L_B" not in body and "L_B=B\\T GB" in supp and "L=\\sum_k T_k\\T K_kT_k" in body, "graph Laplacian L_B (supplement only) distinct from NMA information matrix L")

# ---- figures are code outputs
import subprocess as _sp
tikz = open('figures/covid_AE_study_polygons.tex').read()
raw = dict(re.findall(r"\\newcommand\{\\(\w+?)Raw\}\{([0-9.]+)\}", tikz))
ae = {}
for r in csv.DictReader(open('figures/data/covid_AE_study_polygon_currents.csv')):
    ae[(r['study'], r['t'])] = abs(float(r['signed_current']))
pairs = {"RecoveryA": ("RECOVERY", "A"), "RecoveryB": ("RECOVERY", "B"), "RecoveryC": ("RECOVERY", "C"), "RecoveryD": ("RECOVERY", "D"), "SolidarityA": ("SOLIDARITY", "A"), "SolidarityE": ("SOLIDARITY", "E"), "ActtA": ("ACTT-1", "A"), "CodexA": ("CoDEX", "A"), "RemapA": ("REMAP-CAP", "A")}
must(all(abs(float(raw[k]) - ae[v]) < 1e-9 for k, v in pairs.items()), "Figure S1 TikZ current values equal the A:E data file")
for f in ["COVID19_Circuit_Structure_NoFlow.pdf", "flow_polygon_overview.pdf", "COVID19_MultiTerminal_MultiComparison.pdf", "covid_AE_study_polygons.pdf"]:
    info = _sp.run(["pdfinfo", "figures/" + f], capture_output=True, text=True).stdout
    imgs = _sp.run(["pdfimages", "-list", "figures/" + f], capture_output=True, text=True).stdout.strip().split("\n")[2:]
    must(("Matplotlib" in info or "pdfTeX" in info) and not [x for x in imgs if x.strip()], "%s is a vector code output (matplotlib/pdfTeX, no embedded raster)" % f)
for f in re.findall(r"includegraphics(?:\[[^\]]*\])?\{figures/([^}]+)\}", body + supp):
    must(os.path.exists("figures/" + f), "figure file exists: %s" % f)
# Revision-specific conceptual and provenance safeguards.
must("total signed study contributions" in body and "direct--indirect splits require a pairing rule" in body, "direct/indirect caveat in main article")
must("direct--indirect split is not determined by terminal currents alone" in supp and "conditional on the CSP pairing rule" in supp, "direct/indirect caveat in supplement")
must("should not be interpreted as clinical treatment effects" in body and "not a validated clinical NMA" in supp and "contemporaneous randomized" in main, "CSP-derived COVID inputs explicitly distinguished from clinical NMA")
must("ChatGPT GPT-6" in main and "Claude Opus 5.5" in main and "November 2026" in main, "CUP generative AI disclosure updated")
print("\n%d failed" % len(fails))
if fails: raise SystemExit(1)
