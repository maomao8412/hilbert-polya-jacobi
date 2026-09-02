#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build goldbach_complete_proof.html"""

OUTPUT = "/Coze/Drive/黎曼猜想论文审核/所有对话/主对话/goldbach_complete_proof.html"

# ── helper: make bilingual block ──
def bi_block(zh_html, en_html):
    """Wrap zh and en content for the bi (bilingual) view."""
    return f'<div class="bi-zh">{zh_html}</div>\n<div class="bi-en">{en_html}</div>\n'

# ── CSS ──
CSS = r"""
:root{--bg:#fdfdfb;--text:#1a1a1a;--heading:#1a4a8a;--subtle:#555;--rule:#c9c9c9;--accent:#8a2b1a;--code-bg:#f4f3ef;--green:#0a6b2b}
*{box-sizing:border-box}
body{font-family:"Segoe UI","Helvetica Neue",Arial,"Microsoft YaHei","PingFang SC","Hiragino Sans GB",sans-serif;background:var(--bg);color:var(--text);max-width:none;margin:0;padding:34px 22px 80px;line-height:1.85;font-size:17px}
h1,h2,h3{font-family:"Segoe UI","Helvetica Neue",Arial,"Microsoft YaHei","PingFang SC","Hiragino Sans GB",sans-serif;color:var(--heading);line-height:1.4}
h1{font-size:1.85em;text-align:center;margin:0 0 .3em;border-bottom:2px solid var(--heading);padding-bottom:.4em}
.subtitle{text-align:center;color:var(--subtle);font-style:italic;margin-bottom:2em}
h2{font-size:1.35em;margin-top:2.2em;border-left:4px solid var(--heading);padding-left:.5em}
h3{font-size:1.1em;margin-top:1.6em}
p{margin:.7em 0;text-align:justify}
.formula{text-align:center;margin:1.2em 0;overflow-x:auto;padding:.3em 0}
mjx-container[display="false"]{vertical-align:-0.06em}
mjx-container[display="true"]{overflow-x:auto;overflow-y:hidden;max-width:100%}
code{font-family:"SFMono-Regular",Consolas,monospace;background:var(--code-bg);padding:.1em .35em;border-radius:3px;font-size:.92em}
pre{background:var(--code-bg);padding:.9em 1.1em;border-radius:5px;overflow-x:auto;font-size:.9em;line-height:1.5}
pre code{background:none;padding:0}
table{border-collapse:collapse;margin:1.2em auto;font-size:.92em}
th,td{border:1px solid var(--rule);padding:.4em .7em;text-align:center}
th{background:#eef2f8;color:var(--heading)}
tr:nth-child(even) td{background:#fafaf7}
.fig{text-align:center;margin:1.6em 0}
.fig img{max-width:100%;border:1px solid var(--rule);border-radius:4px}
.figcap{color:var(--subtle);font-size:.88em;margin-top:.4em;font-style:italic}
.links{background:#f0f4fa;border:1px solid #c9d6e8;border-radius:5px;padding:1em 1.4em;margin:1.5em 0}
.links p{margin:.35em 0}
.boundary{background:#f7f5e8;border:1px solid #d8d2a8;border-radius:5px;padding:.8em 1.3em;margin:1.2em 0}
.boundary ul{margin:.4em 0;padding-left:1.5em}
.boundary li{margin:.3em 0}
a{color:var(--heading);text-decoration:none;border-bottom:1px dotted var(--heading)}
a:hover{border-bottom-style:solid}
hr{border:none;border-top:1px solid var(--rule);margin:2em 0}
.hero{background:linear-gradient(135deg,#f8f6f0 0%,#f0f4fa 100%);border:1px solid #d8d2a8;border-radius:8px;padding:1.5em 2em;margin:1.5em 0;text-align:center}
.hero .chain{font-size:1.02em;line-height:2.3;color:var(--heading)}
.hero .chain .arrow{color:var(--accent);font-weight:bold;margin:0 .25em}
.proof-box{background:#eef7f0;border:1px solid #b8d8be;border-radius:6px;padding:.8em 1.4em;margin:1.2em auto;max-width:1000px;text-align:center}
.proof-box .label{font-weight:bold;color:var(--green);font-size:1.05em;display:inline}
.consequence{background:#fdf4e8;border:1px solid #e0c896;border-radius:6px;padding:.6em 1.2em;margin:1em 0}
body{max-width:none}
h1,h2,h3,p,blockquote,.boundary,.links,.hero,.chain-box,.fig,.formula,pre,.subtitle,.zh-link,ol,ul{max-width:1000px;margin-left:auto;margin-right:auto}
.table-scroll{overflow-x:auto;-webkit-overflow-scrolling:touch;touch-action:pan-x pan-y;margin:.8em 0;border:1px solid #eee;border-radius:6px;padding:.4em .5em;background:#fff}
.table-scroll table{min-width:560px;border-collapse:collapse;margin:0 auto}
.table-scroll td,.table-scroll th{white-space:nowrap}
img{display:block;margin:1.2em auto;max-width:100%;width:min(100%,1000px);height:auto;border:1px solid var(--rule);border-radius:4px}
.lang-switcher{display:flex;gap:8px;margin:16px auto;justify-content:center}
.lang-switcher button{padding:6px 16px;border:1px solid var(--rule);background:var(--code-bg);cursor:pointer;border-radius:4px;font-size:14px;font-family:inherit;color:var(--text)}
.lang-switcher button.active{background:var(--heading);color:#fff;border-color:var(--heading)}
.lang-switcher button:hover:not(.active){background:#e8e7e3}
.lang-zh,.lang-en,.lang-bi{display:none}
.lang-zh.visible,.lang-en.visible,.lang-bi.visible{display:block}
.lang-bi .bi-zh,.lang-bi .bi-en{display:block}
.lang-bi .bi-zh{margin-bottom:1.8em;padding-bottom:1.2em;border-bottom:1px dashed var(--rule)}
.lang-bi .bi-en{opacity:.88}
.meta-block{text-align:center;font-size:.92em;color:var(--subtle);margin:1em auto 2em;line-height:2}
.meta-block a{color:var(--heading)}
@media (max-width:768px){
  body{padding:16px 10px 60px;font-size:16px;line-height:1.8}
  h1{font-size:1.45em}h2{font-size:1.18em}h3{font-size:1.05em}
  .hero{padding:1em .7em}.hero .chain{font-size:.9em;line-height:2.05}
  .table-scroll{margin-left:-6px;margin-right:-6px;padding:.3em .2em}
  .formula{font-size:.9em}
  .boundary,.links,.proof-box,.consequence{padding:.7em .9em}
}
"""

NAV = '''<!-- top-nav -->
<nav style="text-align:center;font-size:.85em;margin:0 auto 1.5em;padding:.7em .6em;border-bottom:1px solid #e6e6e6;color:#999;line-height:2.2;max-width:960px;font-family:system-ui,sans-serif"><a href="index.html" style="color:#1a5fb4;text-decoration:none;margin:0 .15em">RH证明(EN)</a> &nbsp;<span style="color:#ccc">·</span>&nbsp; <a href="zh_v4.html" style="color:#1a5fb4;text-decoration:none;margin:0 .15em">RH证明(中文)</a> &nbsp;<span style="color:#ccc">·</span>&nbsp; <a href="grh.html" style="color:#1a5fb4;text-decoration:none;margin:0 .15em">GRH(EN)</a> &nbsp;<span style="color:#ccc">·</span>&nbsp; <a href="grh_zh.html" style="color:#1a5fb4;text-decoration:none;margin:0 .15em">GRH(中文)</a> &nbsp;<span style="color:#ccc">·</span>&nbsp; <a href="physics_web/index.html" style="color:#1a5fb4;text-decoration:none;margin:0 .15em">量子混沌物理论文</a> &nbsp;<span style="color:#ccc">·</span>&nbsp; <a href="qcao_zh.html" style="color:#1a5fb4;text-decoration:none;margin:0 .15em">物理论文(中文)</a> &nbsp;<span style="color:#ccc">·</span>&nbsp; <a href="three_term_viz/index.html" style="color:#1a5fb4;text-decoration:none;margin:0 .15em">三通道可视化</a> &nbsp;<span style="color:#ccc">·</span>&nbsp; <a href="constants.html" style="color:#1a5fb4;text-decoration:none;margin:0 .15em">数论常数</a> &nbsp;<span style="color:#ccc">·</span>&nbsp; <strong style="color:#333;font-weight:700">哥德巴赫猜想</strong> &nbsp;<span style="color:#ccc">·</span>&nbsp; <a href="paper/q_unifies_four_conjectures.html" style="color:#1a5fb4;text-decoration:none;margin:0 .15em">八字曲线之钥</a> &nbsp;<span style="color:#ccc">·</span>&nbsp; <a href="rebuttal.html" style="color:#1a5fb4;text-decoration:none;margin:0 .15em">审稿回应</a> &nbsp;<span style="color:#ccc">·</span>&nbsp; <a href="https://github.com/maomao8412/hilbert-polya-jacobi" style="color:#1a5fb4;text-decoration:none;margin:0 .15em" target="_blank" rel="noopener">GitHub代码仓库</a></nav>'''

LANG_SWITCHER = '''<!-- Language Switcher -->
<div class="lang-switcher">
  <button onclick="switchLang('zh')" class="active" id="btn-zh">中文</button>
  <button onclick="switchLang('en')" id="btn-en">English</button>
  <button onclick="switchLang('bi')" id="btn-bi">中英对照</button>
</div>'''

JS = '''<script>
function switchLang(lang) {
  document.querySelectorAll('.lang-zh,.lang-en,.lang-bi').forEach(el => el.classList.remove('visible'));
  document.querySelectorAll('.lang-switcher button').forEach(b => b.classList.remove('active'));
  document.querySelector('.lang-' + lang).classList.add('visible');
  document.getElementById('btn-' + lang).classList.add('active');
  if (window.MathJax && MathJax.typeset) { try { MathJax.typeset(); } catch(e){} }
}
document.addEventListener('DOMContentLoaded', function(){ switchLang('zh'); });
</script>'''

MATHJAX = '''<script>
window.MathJax = {
  tex: { inlineMath: [['\\\\(', '\\\\)']], displayMath: [['$$','$$'], ['\\\\[','\\\\]']], tags: 'none' },
  svg: { fontCache: 'global' },
  options: { skipHtmlTags: ['script','noscript','style','textarea','pre','code'] }
};
</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js" async></script>'''

# We'll assemble the file from parts
parts = []

# Part 1: head
parts.append(f'''<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>哥德巴赫猜想的分圆塔证明 | Cyclotomic Tower Proof of Goldbach's Conjecture</title>
{MATHJAX}
<style>{CSS}</style>
</head>
<body>
{NAV}
{LANG_SWITCHER}
''')

# Now I need to define all content and wrap in lang-zh, lang-en, lang-bi
# Due to file size, I'll write sections to a file and assemble

# Write the complete HTML via Python
with open(OUTPUT, 'w', encoding='utf-8') as f:
    f.write(parts[0])
    # Will continue writing sections below

print("Head written OK")
