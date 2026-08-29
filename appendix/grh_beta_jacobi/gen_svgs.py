# -*- coding: utf-8 -*-
"""Generate the two zoomable SVGs for the GRH page from beta_J50_results.json:
   1) structure SVG: 50x50 tridiagonal (blue diag, gold off-diag), white in-cell values
   2) spectrum SVG: 50 diagonal cells (green locked, gold converging), 1/sqrt(lambda)
"""
import json, math

d = json.load(open("/tmp/grh_beta/beta_J50_results.json"))
N = d["N"]
a = [float(x) for x in d["diag_a"]]
b = [float(x) for x in d["offdiag_b"]]
ev = d["eigenvalues"]
inv = d["inv_sqrt"]
zeros = d["beta_zeros_independent"]
rows = {r["k"]: r for r in d["comparison"]}

def fmt(v):
    s = f"{v:.2e}"
    return s.replace("e-0", "e-").replace("e+0", "e+").replace("e-","e-")

VIEW = 574; M0 = 11; CELL = 11.04  # same geometry as RH page
def xy(i): return 12.6 + i*CELL

def text(x, y, t, size=2.7, col="#fff", weight="bold"):
    return (f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" fill="{col}" '
            f'text-anchor="middle" dominant-baseline="central" font-weight="{weight}" '
            f'font-family="Consolas,Monaco,monospace" pointer-events="none">{t}</text>')

def rect(x, y, col, title, label):
    return (f'<rect x="{x:.1f}" y="{y:.1f}" width="9.8" height="9.8" fill="{col}">'
            f'<title>{title}</title>{text(x+4.9, y+4.9, label)}</rect>')

# ---------- structure SVG ----------
parts = [f'<svg xmlns="http://www.w3.org/2000/svg" class="zoomable" viewBox="0 0 {VIEW} {VIEW}" '
         f'style="width:100%;min-width:760px;max-width:760px;height:auto;background:#fff;'
         f'border:1px solid #e0e0e0;border-radius:6px;display:block;margin:.6em auto">',
         f'<rect x="11" y="11" width="552" height="552" fill="none" stroke="#888" stroke-width="1"/>']
BLUE = "#1a4d8f"; GOLD = "#b8860b"
for i in range(N):
    # diagonal
    x = y = xy(i)
    parts.append(rect(x, y, BLUE, f"a[{i+1},{i+1}] = {a[i]:.6e}", fmt(a[i])))
    if i < N-1:
        bo = b[i]
        parts.append(rect(xy(i+1), xy(i), GOLD, f"b[{i+1},{i+2}] = {bo:.6e}", fmt(bo)))
        parts.append(rect(xy(i), xy(i+1), GOLD, f"b[{i+2},{i+1}] = {bo:.6e}", fmt(bo)))
# axis labels
for lab, idx in [(1,0),(10,9),(20,19),(30,29),(40,39),(50,49)]:
    cx = xy(idx)+4.9; cy = xy(idx)+4.9
    parts.append(text(cx, 8, str(lab), size=8, col="#666", weight="normal"))
    parts.append(text(8, cy, str(lab), size=8, col="#666", weight="normal"))
parts.append("</svg>")
open("/tmp/grh_web/svg_structure.html","w").write("".join(parts))

# ---------- spectrum SVG: 50 diagonal cells stacked ----------
# compare 1/sqrt(lambda_k) with actual beta zeros; locked if rel err < 1e-4
parts = [f'<svg xmlns="http://www.w3.org/2000/svg" class="zoomable" viewBox="0 0 {VIEW} {VIEW}" '
         f'style="width:100%;min-width:760px;max-width:760px;height:auto;background:#fff;'
         f'border:1px solid #e0e0e0;border-radius:6px;display:block;margin:.6em auto">',
         f'<rect x="11" y="11" width="552" height="552" fill="none" stroke="#888" stroke-width="1"/>']
GREEN = "#0a5c23"
nlocked = 0
for k in range(N):
    x = y = xy(k)
    r = rows.get(k+1)
    lam = ev[k]
    isl = inv[k]
    if r and r["rel_err"] < 1e-4:
        col = GREEN; nlocked += 1
        t = (f"lambda_{k+1} = {lam:.6e} -&gt; 1/sqrt(lambda) = {isl:.10f}  "
             f"LOCKED to beta gamma_{k+1} = {zeros[k]:.10f}")
    else:
        col = GOLD
        g = zeros[k] if k < len(zeros) else float("nan")
        t = (f"lambda_{k+1} = {lam:.6e} -&gt; 1/sqrt(lambda) = {isl:.6f}  "
             f"(beta gamma_{k+1} = {g:.6f}, converging with order)")
    parts.append(rect(x, y, col, t, fmt(lam)))
for lab, idx in [(1,0),(10,9),(20,19),(30,29),(40,39),(50,49)]:
    cx = xy(idx)+4.9; cy = xy(idx)+4.9
    parts.append(text(cx, 8, str(lab), size=8, col="#666", weight="normal"))
    parts.append(text(8, cy, str(lab), size=8, col="#666", weight="normal"))
parts.append("</svg>")
open("/tmp/grh_web/svg_spectrum.html","w").write("".join(parts))
print("locked:", nlocked, "/", min(N, len(zeros)))
print("gamma1 =", zeros[0], " inv_sqrt1 =", inv[0])
