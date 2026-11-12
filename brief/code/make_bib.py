#!/usr/bin/env python3
r"""Build an AMA-style thebibliography ordered by first citation.
Usage: make_bib.py OUT.tex FILE1.tex [FILE2.tex ...]
Keys are taken from \cite{...} / \citep{...} in order of first appearance.
"""
import re, sys

AMA = {
"Rucker2012": r"R\"ucker G. Network meta-analysis, electrical networks and graph theory. \emph{Res Synth Methods}. 2012;3(4):312--324. doi:10.1002/jrsm.1058",
"Konig2013": r"K\"onig J, Krahn U, Binder H. Visualizing the flow of evidence in network meta-analysis and characterizing mixed treatment comparisons. \emph{Stat Med}. 2013;32(30):5414--5429. doi:10.1002/sim.6001",
"Papakonstantinou2018": r"Papakonstantinou T, Nikolakopoulou A, R\"ucker G, et al. Estimating the contribution of studies in network meta-analysis: paths, flows and streams. \emph{F1000Research}. 2018;7:610. doi:10.12688/f1000research.14770.3",
"Davies2022": r"Davies AL, Papakonstantinou T, Nikolakopoulou A, R\"ucker G, Galla T. Network meta-analysis and random walks. \emph{Stat Med}. 2022;41(12):2091--2114. doi:10.1002/sim.9346",
"Rucker2024": r"R\"ucker G, Papakonstantinou T, Nikolakopoulou A, Schwarzer G, Galla T, Davies AL. Shortest path or random walks? A framework for path weights in network meta-analysis. \emph{Stat Med}. 2024;43(22):4287--4304. doi:10.1002/sim.10177",
"Davies2026": r"Davies AL. The bipartite structure of treatment-trial networks reveals the flow of information in network meta-analysis. \emph{J R Stat Soc Ser A Stat Soc}. Published online February 18, 2026. doi:10.1093/jrsssa/qnag012",
"RuckerDavies2026": r"R\"ucker G, Davies AL, Schwarzer G. Network meta-analysis and diffusion. \emph{Res Synth Methods}. Published online October 15, 2026. doi:10.1017/rsm.2026.10123",
"DoyleSnell1984": r"Doyle PG, Snell JL. \emph{Random Walks and Electric Networks}. Carus Mathematical Monographs 22. Mathematical Association of America; 1984.",
"RuckerSchwarzer2014": r"R\"ucker G, Schwarzer G. Reduce dimension or reduce weights? Comparing two approaches to multi-arm studies in network meta-analysis. \emph{Stat Med}. 2014;33(25):4353--4369. doi:10.1002/sim.6236",
"RuckerSchwarzer2025": r"R\"ucker G, Schwarzer G. Trials and triangles: network meta-analysis of multi-arm trials with correlated arms. \emph{Res Synth Methods}. 2025;16(6):961--974. doi:10.1017/rsm.2025.10026",
"Balduzzi2023": r"Balduzzi S, R\"ucker G, Nikolakopoulou A, et al. netmeta: an R package for network meta-analysis using frequentist methods. \emph{J Stat Softw}. 2023;106(2):1--40. doi:10.18637/jss.v106.i02",
"WangCSP2026": r"Wang C, Zhang Y, Jin Z, O'Connor AM. Contrast-space projection for network meta-analysis: an exact and invariant study-based decomposition of direct and indirect contributions. \emph{arXiv}. Preprint posted online April 23, 2026 (version 2, August 30, 2026). arXiv:2604.21994. doi:10.48550/arXiv.2604.21994",
"ZhangWangOConnor2026": r"Zhang Y, Wang C, O'Connor AM. NMAforest: an R package for generating forest plots in network meta-analysis with decomposition of direct, indirect, and study-level evidence contributions. \emph{R J}. 2026;18(2):175--186. doi:10.32614/RJ-2026-030",
"DorflerBullo2013": r"D\"orfler F, Bullo F. Kron reduction of graphs with applications to electrical networks. \emph{IEEE Trans Circuits Syst I Regul Pap}. 2013;60(1):150--163. doi:10.1109/TCSI.2012.2215780",
"Millani2020": r"Garlet Millani M, Molter H, Niedermeier R, Sorge M. Efficient algorithms for measuring the funnel-likeness of DAGs. \emph{J Comb Optim}. 2020;39(1):216--245. doi:10.1007/s10878-019-00464-4",
"Khan2022": r"Khan S, Kortelainen M, C\'aceres M, Williams L, Tomescu AI. Improving RNA assembly via safety and completeness in flow decompositions. \emph{J Comput Biol}. 2022;29(12):1270--1287. doi:10.1089/cmb.2022.0261",
"DesoerKuh1969": r"Desoer CA, Kuh ES. \emph{Basic Circuit Theory}. McGraw-Hill; 1969.",
"vanDerSchaft2010": r"van der Schaft A. Characterization and partial synthesis of the behavior of resistive circuits at their terminals. \emph{Syst Control Lett}. 2010;59(7):423--428. doi:10.1016/j.sysconle.2010.05.005",
"Siemieniuk2020": r"Siemieniuk RAC, Bartoszko JJ, Ge L, et al. Drug treatments for covid-19: living systematic review and network meta-analysis. \emph{BMJ}. 2020;370:m2980. doi:10.1136/bmj.m2980",
"Sbidian2021": r"Sbidian E, Chaimani A, Garcia-Doval I, et al. Systemic pharmacological treatments for chronic plaque psoriasis: a network meta-analysis. \emph{Cochrane Database Syst Rev}. 2021;(4):CD011535. doi:10.1002/14651858.CD011535.pub4",
"Barrett2000": r"Barrett C, Jacob R, Marathe M. Formal-language-constrained path problems. \emph{SIAM J Comput}. 2000;30(3):809--837. doi:10.1137/S0097539798337716",
}

out, files = sys.argv[1], sys.argv[2:]
order = []
for fn in files:
    s = open(fn, encoding="utf-8").read()
    s = re.sub(r"(?m)^%.*$", "", s)
    for m in re.finditer(r"\\cite[pt]?\*?(?:\[[^\]]*\])?\{([^}]*)\}", s):
        for k in m.group(1).split(","):
            k = k.strip()
            if k and k not in order:
                order.append(k)
missing = [k for k in order if k not in AMA]
if missing:
    sys.exit("missing AMA entries: " + ", ".join(missing))
with open(out, "w", encoding="utf-8") as f:
    f.write("\\begin{thebibliography}{99}\\setlength{\\itemsep}{1pt}\\small\n")
    for k in order:
        f.write("\\bibitem{%s} %s\n" % (k, AMA[k]))
    f.write("\\end{thebibliography}\n")
print("wrote", out, "with", len(order), "entries")
