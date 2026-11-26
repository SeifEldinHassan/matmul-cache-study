import csv, statistics as st
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, AQUA = "#257de0", "#f47230", "#00b894"
INK, MUTED, GRID = "#191919", "#666666", "#e6e6e6"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#999999",
    "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False,
})


def load(f):
    rows = []
    for l in open(f):
        l = l.strip()
        if not l or l.startswith("#") or l.startswith("version"):
            continue
        v, n, b, r, s = l.split(",")
        rows.append((v, int(n), int(b), int(r), float(s)))
    return rows


sizes = load("results_sizes.csv")
blocks = load("results_blocks.csv")
d = defaultdict(list)
for v, n, b, r, s in sizes:
    d[(v, n)].append(s)
Ns = sorted({k[1] for k in d})
style = {"ijk": (BLUE, "o", "-", "Plain (i-j-k)"),
         "ikj": (ORANGE, "s", "--", "Reordered (i-k-j)"),
         "blocked": (AQUA, "^", ":", "Blocked (block = 32)")}

# ---- Figure 1: runtime vs N ----
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
for v, (c, m, ls, lab) in style.items():
    med = [st.median(d[(v, n)]) for n in Ns]
    lo = [min(d[(v, n)]) for n in Ns]
    hi = [max(d[(v, n)]) for n in Ns]
    ax.fill_between(Ns, lo, hi, color=c, alpha=0.15, lw=0)
    ax.plot(Ns, med, color=c, marker=m, ls=ls, lw=2, ms=6, mec="white", mew=1.0, label=lab)
ax.set_xscale("log", base=2); ax.set_yscale("log")
ax.set_xticks(Ns); ax.set_xticklabels([str(n) for n in Ns])
ax.set_xlabel("Matrix size N (N x N, double precision)")
ax.set_ylabel("Runtime (seconds, log scale)")
ax.set_title("Runtime of matrix multiplication by version", loc="left", fontsize=11,
             fontweight="bold", color=INK)
ax.legend(frameon=False, loc="upper left")
ax.annotate("%.2f s" % st.median(d[("ijk", 1024)]), (1024, st.median(d[("ijk", 1024)])), xytext=(-8, 6),
            textcoords="offset points", ha="right", color=INK, fontsize=9)
ax.annotate("%.2f s" % st.median(d[("ikj", 1024)]), (1024, st.median(d[("ikj", 1024)])), xytext=(-8, 8),
            textcoords="offset points", ha="right", color=INK, fontsize=9)
ax.annotate("%.2f s" % st.median(d[("blocked", 1024)]), (1024, st.median(d[("blocked", 1024)])), xytext=(-8, -14),
            textcoords="offset points", ha="right", color=INK, fontsize=9)
fig.text(0.01, 0.005, "Lines: median of 10 runs. Shaded band: minimum to maximum.",
         fontsize=8, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig("fig1_runtime_vs_size.png", dpi=220, facecolor="white")
plt.close(fig)

# ---- Figure 2: speedup vs plain ----
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.grid(True, axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
w = 0.38
xs = range(len(Ns))
for off, v in ((-w / 2 - 0.01, "ikj"), (w / 2 + 0.01, "blocked")):
    c, m, ls, lab = style[v]
    vals = [st.median(d[("ijk", n)]) / st.median(d[(v, n)]) for n in Ns]
    bars = ax.bar([x + off for x in xs], vals, width=w, color=c, label=lab,
                  edgecolor="white", linewidth=1.5,
                  hatch="//" if v == "blocked" else None)
    for x, val in zip(xs, vals):
        ax.text(x + off, val + 0.1, "%.1f" % val, ha="center", va="bottom",
                fontsize=8.5, color=INK)
ax.axhline(1, color="#999999", lw=1)
ax.set_xticks(list(xs)); ax.set_xticklabels([str(n) for n in Ns])
ax.set_xlabel("Matrix size N")
ax.set_ylabel("Speedup compared with plain version (x)")
ax.set_title("How much faster than the plain loop", loc="left", fontsize=11,
             fontweight="bold", color=INK)
ax.legend(frameon=False, loc="upper left")
ax.set_ylim(0, max(st.median(d[("ijk", n)]) / st.median(d[("ikj", n)]) for n in Ns) * 1.18)
fig.text(0.01, 0.005, "Speedup = median runtime of plain version / median runtime of the other version.",
         fontsize=8, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig("fig2_speedup.png", dpi=220, facecolor="white")
plt.close(fig)

# ---- Figure 3: block size sweep ----
e = defaultdict(list)
for v, n, b, r, s in blocks:
    e[b].append(s)
Bs = sorted(e)
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.grid(True, axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
med = [st.median(e[b]) for b in Bs]
lo = [m - min(e[b]) for m, b in zip(med, Bs)]
hi = [max(e[b]) - m for m, b in zip(med, Bs)]
ax.bar(range(len(Bs)), med, width=0.55, color=AQUA, edgecolor="white", linewidth=1.5,
       yerr=[lo, hi], error_kw=dict(ecolor=INK, lw=1, capsize=4))
ikj = st.median(d[("ikj", 1024)])
ax.axhline(ikj, color=ORANGE, lw=2, ls="--")
ax.text(-0.4, max(med) * 1.3, "Dashed line: reordered version without blocking (%.2f s)" % ikj,
        ha="left", va="bottom", color=INK, fontsize=9)
for x, m in zip(range(len(Bs)), med):
    ax.text(x, max(med) * 0.03, "%.2f" % m, ha="center", va="bottom", color="white" if False else INK,
            fontsize=9, fontweight="bold")
ax.set_xticks(range(len(Bs))); ax.set_xticklabels([str(b) for b in Bs])
ax.set_xlabel("Block size (elements per side)")
ax.set_ylabel("Runtime at N = 1024 (seconds)")
ax.set_ylim(0, max(med) * 1.45)
ax.set_title("Effect of block size (N = 1024)", loc="left", fontsize=11,
             fontweight="bold", color=INK)
fig.text(0.01, 0.005, "Bars: median of 10 runs. Whiskers: minimum to maximum.",
         fontsize=8, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig("fig3_block_size.png", dpi=220, facecolor="white")
plt.close(fig)

# ---- summary table CSV ----
with open("summary_table.csv", "w", newline="") as f:
    w_ = csv.writer(f)
    w_.writerow(["N", "version", "median_s", "mean_s", "sd_s", "min_s", "max_s", "speedup_vs_plain"])
    for n in Ns:
        base = st.median(d[("ijk", n)])
        for v in ("ijk", "ikj", "blocked"):
            x = d[(v, n)]
            w_.writerow([n, v, "%.6f" % st.median(x), "%.6f" % st.mean(x), "%.6f" % st.stdev(x),
                         "%.6f" % min(x), "%.6f" % max(x), "%.2f" % (base / st.median(x))])
print("done")
