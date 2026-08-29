# -*- coding: utf-8 -*-
"""双函数零点交替锁定视图：生成 HTML 片段并注入 grh_zh.html / grh.html，主篇加入口。"""
import json, io

ZETA_JSON = '/tmp/hpgit/data/J50_results.json'
BETA_JSON = '/tmp/grh_beta/beta_J50_results.json'
GAMMA_MAX = 100.0

z = json.load(open(ZETA_JSON))
z50 = [c for c in z if c['N'] == 50][0]
beta = json.load(open(BETA_JSON))
bcmp = beta['comparison']
bref = [float(x) for x in beta['beta_zeros_independent']]

dots = []
for row in z50['zeros']:
    g = float(row['known'])
    if g < GAMMA_MAX:
        dots.append((g, 'zeta', row['idx'], bool(row['locked']),
                     float(row['gamma_est']), g, abs(float(row['err']) / g)))
for row in bcmp:
    g = row['gamma']; k = row['k']
    if g < GAMMA_MAX:
        dots.append((g, 'beta', k, k <= beta['n_locked'],
                     row['inv_sqrt'], bref[k - 1], row['rel_err']))
dots.sort(key=lambda d: (d[0], 0 if d[1] == 'beta' else 1))
locked = [d for d in dots if d[3]]
print(f'dots={len(dots)} locked={len(locked)} (beta {sum(1 for d in locked if d[1]=="beta")}, zeta {sum(1 for d in locked if d[1]=="zeta")})')

W, H = 940, 150
X0, X1 = 30, 910
AXIS_Y = 75
ZETA_Y, BETA_Y = 38, 112
def xof(g): return X0 + (g / GAMMA_MAX) * (X1 - X0)
GOLD, TEAL = '#b8860b', '#0a5c23'

parts = [f'<svg id="dual-rail" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
         'style="width:100%;min-width:760px;max-width:940px;height:auto;background:#fff;'
         'border:1px solid #e0e0e0;border-radius:6px;display:block;margin:.6em auto;overflow:visible">']
parts.append(f'<line x1="{X0}" y1="{AXIS_Y}" x2="{X1}" y2="{AXIS_Y}" stroke="#999" stroke-width="1"/>')
for t in range(0, 101, 10):
    x = xof(t)
    parts.append(f'<line x1="{x:.1f}" y1="{AXIS_Y-4}" x2="{x:.1f}" y2="{AXIS_Y+4}" stroke="#999" stroke-width="1"/>')
    parts.append(f'<text x="{x:.1f}" y="{AXIS_Y+18}" font-size="11" fill="#666" text-anchor="middle" font-family="Consolas,Monaco,monospace">{t}</text>')
parts.append(f'<text x="{X0}" y="{ZETA_Y-14}" font-size="12" fill="{GOLD}" font-weight="bold">&#950;(s) 黎曼零点（J50 锁定前 25 个）</text>')
parts.append(f'<text x="{X0}" y="{BETA_Y+22}" font-size="12" fill="{TEAL}" font-weight="bold">&#946;(s)=L(s,&#967;&#8324;) 零点（J50 锁定前 27 个）</text>')

n = 0
for g, fn, idx, lk, est, ref, rel in dots:
    x = xof(g)
    y = ZETA_Y if fn == 'zeta' else BETA_Y
    color = GOLD if fn == 'zeta' else TEAL
    if lk:
        n += 1
        tip = f'第 {n} 个被锁定：{"&#950;" if fn=="zeta" else "&#946;"} &#947;_{idx} = {ref:.10f}'
        parts.append(f'<circle class="dual-dot locked" data-seq="{n}" data-lane="{fn}" cx="{x:.2f}" cy="{y}" r="4.5" '
                     f'fill="{color}" fill-opacity="0" stroke="{color}" stroke-width="2" style="cursor:pointer">'
                     f'<title>{tip}</title></circle>')
    else:
        nm = '&#950;' if fn == 'zeta' else '&#946;'
        tip = f'{nm} &#947;_{idx} &#8776; {ref:.6f}（未锁定，随阶数收敛中）'
        parts.append(f'<circle class="dual-dot" cx="{x:.2f}" cy="{y}" r="4.0" fill="#fff" '
                     f'stroke="{color}" stroke-width="1.2" stroke-opacity="0.45"><title>{tip}</title></circle>')
parts.append('</svg>')
svg = ''.join(parts)

def row_html(d):
    g, fn, idx, lk, est, ref, rel = d
    name = '&#950;' if fn == 'zeta' else '&#946;'
    fname = '黎曼 &#950;(s)' if fn == 'zeta' else 'L 函数 &#946;(s)'
    cls = 'tr-zeta' if fn == 'zeta' else 'tr-beta'
    return (f'<tr class="{cls}" style="cursor:pointer" data-g="{g:.6f}"><td>{name}</td>'
            f'<td>{fname}</td><td>{idx}</td><td>{g:.6f}</td><td>{rel:.2e}</td></tr>')
rows = '\n'.join(row_html(d) for d in locked)

CSS = """
<style>
.dual-dot.locked{transition:fill-opacity .25s}
.dual-dot.locked.flash{filter:drop-shadow(0 0 6px currentColor)}
.dual-bar{display:flex;gap:12px;flex-wrap:wrap;align-items:center;margin:.6em 0;font-size:.92em}
.dual-bar button{padding:6px 16px;border-radius:6px;border:1px solid #1a4d8f;background:#1a4d8f;color:#fff;cursor:pointer;font-size:.95em}
.dual-bar .lg{display:inline-flex;align-items:center;gap:6px;color:var(--subtle,#666)}
.dual-bar .sw{width:12px;height:12px;border-radius:50%;display:inline-block}
#dual-readout{min-height:2.4em;margin:.5em 0;padding:.5em .9em;background:rgba(26,77,143,.06);
  border-left:3px solid #1a4d8f;border-radius:4px;font-size:.95em;font-family:Consolas,Monaco,monospace}
.dual-table tr.tr-beta td{background:rgba(10,92,35,.05)}
.dual-table tr.tr-zeta td{background:rgba(184,134,11,.06)}
</style>
"""

JS = """
<script>
(function(){
  var svg=document.getElementById('dual-rail'); if(!svg) return;
  var locked=[].slice.call(svg.querySelectorAll('.dual-dot.locked'));
  var ro=document.getElementById('dual-readout');
  var timer=null, running=false;
  function msg(t){ if(ro) ro.textContent=t; }
  function reset(){
    clearInterval(timer); running=false;
    locked.forEach(function(d){d.setAttribute('fill-opacity','0'); d.classList.remove('flash');});
    var btn=document.getElementById('dual-play'); if(btn) btn.textContent=BTN_PLAY;
    msg(HINT_IDLE);
  }
  function play(){
    if(running){reset(); return;}
    reset(); running=true;
    var btn=document.getElementById('dual-play'); if(btn) btn.textContent=BTN_REPLAY;
    var i=0;
    timer=setInterval(function(){
      if(i>=locked.length){clearInterval(timer); running=false; msg(MSG_DONE); return;}
      var d=locked[i];
      d.setAttribute('fill-opacity','1'); d.classList.add('flash');
      var ti=d.querySelector('title');
      msg(''+(i+1)+' / '+locked.length+'  '+ (ti?ti.textContent:''));
      (function(dot){setTimeout(function(){dot.classList.remove('flash');},350);})(d);
      i++;
    },420);
  }
  var btn=document.getElementById('dual-play');
  if(btn) btn.addEventListener('click',play);
  svg.querySelectorAll('.dual-dot').forEach(function(c){
    c.addEventListener('click',function(){ var ti=c.querySelector('title'); if(ti) msg(ti.textContent); });
  });
  document.querySelectorAll('.dual-table tbody tr').forEach(function(tr){
    tr.addEventListener('click',function(){
      var g=+tr.getAttribute('data-g');
      var tx=__X0__+g/__GM__*(__X1__-__X0__);
      var best=null,bd=1e9;
      svg.querySelectorAll('.dual-dot').forEach(function(c){
        var dd=Math.abs((+c.getAttribute('cx'))-tx);
        if(dd<bd){bd=dd;best=c;}
      });
      if(best){var ti=best.querySelector('title'); if(ti) msg(ti.textContent);}
    });
  });
  reset();
})();
</script>""".replace('__X0__', str(X0)).replace('__X1__', str(X1)).replace('__GM__', str(GAMMA_MAX))

ZH_JS = JS.replace('BTN_PLAY', "'▶ 播放锁定过程'").replace('BTN_REPLAY', "'⟲ 重播'") \
          .replace('HINT_IDLE', "'点击播放：零点按高度 γ 从小到大，逐个被两个 J50 矩阵锁定（金色 ζ，绿色 β）。点圆点或表格行可读全精度。'") \
          .replace('MSG_DONE', "'锁定完成：52 个零点（ζ 25 个 + β 27 个），两个 50×50 矩阵的谱与独立零点全部对上。'")
EN_JS = JS.replace('BTN_PLAY', "'▶ Play the locking sequence'").replace('BTN_REPLAY', "'⟲ Replay'") \
          .replace('HINT_IDLE', "'Press play: zeros are locked one by one in ascending γ by the two J50 matrices (gold ζ, green β). Click a dot or a table row for full precision.'") \
          .replace('MSG_DONE', "'Locking complete: 52 zeros (25 ζ + 27 β) — spectra of two 50×50 matrices match independent zero locations.'")

ZH = """
<h2 id="dual" style="border-left:4px solid #b8860b;padding-left:10px">★ 双函数零点交替锁定：ζ（RH）与 β（GRH）同台</h2>
<p>两个 50×50 Jacobi 矩阵，一台用 \\(\\zeta(s)\\) 的矩序列构造（主篇），一台用 \\(\\beta(s)=L(s,\\chi_4)\\) 的矩序列构造（本页）。把两边锁定的零点<strong>按虚部高度 \\(\\gamma\\) 排在同一条轴上</strong>——金色是黎曼零点（RH），绿色是 L 函数零点（GRH）：它们沿着临界线<strong>交替出现、交替被锁定</strong>。上排金点 25 个锁定，下排绿点 27 个锁定，空圈是矩阵阶数尚在收敛中的更高零点。</p>
__CSS__
<div class="dual-bar">
  <button id="dual-play">▶ 播放锁定过程</button>
  <span class="lg"><span class="sw" style="background:#b8860b"></span>ζ 黎曼零点（RH）·锁定 25</span>
  <span class="lg"><span class="sw" style="background:#0a5c23"></span>β = L(s,χ₄) 零点（GRH）·锁定 27</span>
  <span class="lg"><span class="sw" style="background:#fff;border:1.5px solid #999"></span>收敛中</span>
</div>
<div class="table-scroll">__SVG__</div>
<div id="dual-readout"></div>
<p class="figcap">图：γ ∈ [0,100] 区间内，ζ 有 29 个非平凡零点、β 有 50 个；两个矩阵共锁定 52 个（ζ 25 + β 27），其余空圈随矩阵阶数收敛。第一个被锁定的是 β 的 γ₁=6.02095，第一个 ζ 零点 γ₁=14.13473 排第四位。</p>
<p>锁定顺序全表（按 γ 高度，两函数交替）：</p>
<div class="table-scroll">
<table class="dual-table"><thead><tr><th>函数</th><th>名称</th><th>第 n 个零点</th><th>γ（虚部）</th><th>相对误差</th></tr></thead>
<tbody>
__ROWS__
</tbody></table>
</div>
<p style="font-size:.9em;color:var(--subtle);text-align:center">两个函数、两台矩阵、同一种锁定。这不是 GRH 的解析证明——谱与零点重合是<strong>可独立复算的数值事实</strong>；角单调性解析论证见 Zenodo 存档论文，本页只摆能复算的东西。</p>
__JS__
""".replace('__CSS__', CSS).replace('__SVG__', svg).replace('__ROWS__', rows).replace('__JS__', ZH_JS)

EN = """
<h2 id="dual" style="border-left:4px solid #b8860b;padding-left:10px">★ Interlaced zeros of two functions: ζ (RH) and β (GRH) on one rail</h2>
<p>Two 50×50 Jacobi matrices — one built from the moments of \\(\\zeta(s)\\) (main paper), the other from those of \\(\\beta(s)=L(s,\\chi_4)\\) (this page). Place the zeros of <strong>both</strong> on one rail, sorted by height \\(\\gamma\\): gold dots are Riemann zeros (RH), green dots are L-function zeros (GRH) — <strong>interlaced along the critical line and locked in alternating order</strong>. The gold lane locks 25, the green lane 27; open rings are higher zeros still converging with matrix order.</p>
__CSS__
<div class="dual-bar">
  <button id="dual-play">▶ Play the locking sequence</button>
  <span class="lg"><span class="sw" style="background:#b8860b"></span>ζ Riemann zeros (RH) · 25 locked</span>
  <span class="lg"><span class="sw" style="background:#0a5c23"></span>β = L(s,χ₄) zeros (GRH) · 27 locked</span>
  <span class="lg"><span class="sw" style="background:#fff;border:1.5px solid #999"></span>converging</span>
</div>
<div class="table-scroll">__SVG__</div>
<div id="dual-readout"></div>
<p class="figcap">Fig: in γ ∈ [0,100], ζ has 29 nontrivial zeros and β has 50; the two matrices lock 52 in total (25 ζ + 27 β). The first locked zero is β γ₁=6.02095; the first ζ zero γ₁=14.13473 is locked fourth.</p>
<p>Full locking sequence (sorted by γ, alternating between the two functions):</p>
<div class="table-scroll">
<table class="dual-table"><thead><tr><th>Fn</th><th>Name</th><th>zero #</th><th>γ (imag. part)</th><th>rel. error</th></tr></thead>
<tbody>
__ROWS__
</tbody></table>
</div>
<p style="font-size:.9em;color:var(--subtle);text-align:center">Two functions, two matrices, one locking pattern. This is not an analytic proof of GRH — the spectrum-zero match is an <strong>independently reproducible numerical fact</strong>; the analytic angular-monotonicity argument is in the Zenodo-archived paper. This page shows only what can be re-computed.</p>
__JS__
""".replace('__CSS__', CSS).replace('__SVG__', svg).replace('__ROWS__',
      rows.replace('黎曼 &#950;(s)', 'Riemann &#950;(s)').replace('L 函数 &#946;(s)', 'L-function &#946;(s)')).replace('__JS__', EN_JS)

ANCHOR_ZH = '<h2>五、46 个本原特征：角单调性网格</h2>'
ANCHOR_EN = '<h2>5. The 46 primitive characters: angular-monotonicity grids</h2>'
for path, anchor, frag in [('/tmp/hpgit/grh_zh.html', ANCHOR_ZH, ZH),
                           ('/tmp/hpgit/grh.html', ANCHOR_EN, EN)]:
    s = io.open(path, encoding='utf-8').read()
    assert anchor in s, path
    if 'dual-rail' in s:
        print('skip (already injected) ->', path); continue
    s = s.replace(anchor, frag + '\n' + anchor, 1)
    io.open(path, 'w', encoding='utf-8').write(s)
    print('injected ->', path)

ZH_LINK = ' <a href="grh_zh.html#dual" style="white-space:nowrap">★ 双函数零点交替锁定（ζ/β 动画）→</a>'
EN_LINK = ' <a href="grh.html#dual" style="white-space:nowrap">★ Interlaced ζ/β zero locking (animation) →</a>'
p = '/tmp/hpgit/zh.html'
s = io.open(p, encoding='utf-8').read()
if 'grh_zh.html' in s and '交替锁定' not in s:
    i = s.find('grh_zh.html'); j = s.find('</a>', i) + len('</a>')
    s = s[:j] + ZH_LINK + s[j:]
    io.open(p, 'w', encoding='utf-8').write(s); print('main zh link added')
p = '/tmp/hpgit/index.html'
s = io.open(p, encoding='utf-8').read()
if 'grh.html' in s and 'Interlaced' not in s:
    i = s.find('grh.html'); j = s.find('</a>', i) + len('</a>')
    s = s[:j] + EN_LINK + s[j:]
    io.open(p, 'w', encoding='utf-8').write(s); print('main en link added')
