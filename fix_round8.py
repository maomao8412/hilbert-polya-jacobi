# -*- coding: utf-8 -*-
import io
FILES = ["zh.html", "index.html"]

def patch(fn):
    with io.open(fn, encoding="utf-8") as f:
        s = f.read()
    orig = s
    def rep(old, new, n=1):
        nonlocal s
        c = s.count(old)
        assert c == n, (fn, "anchor count %d != %d for: %s" % (c, n, old[:70]))
        s = s.replace(old, new)

    # 1. copyright block: center, full gold border
    rep('<div style="margin:1.4em 0;padding:.85em 1.15em;border-left:4px solid #b8860b;background:#fbf8ef;border-radius:4px;font-size:.92em;line-height:1.75">',
        '<div style="margin:1.4em auto;padding:.85em 1.15em;border:1px solid #b8860b;background:#fbf8ef;border-radius:6px;font-size:.92em;line-height:1.75;text-align:center;max-width:920px">')

    # 1b. no Chinese name, no word "proof"
    if fn == "zh.html":
        rep('本页所述证明与构造由<strong>陈倬（Zhuo Chen）</strong>完成，全文已于 2026 年 8 月在 Zenodo 公开存档（带时间戳 DOI，可独立核验优先权）：RH 证明 <a href="https://zenodo.org/records/22113226">',
            '本页所述工作由 <strong>Zhuo Chen</strong> 完成，全文已于 2026 年 8 月在 Zenodo 公开存档（带时间戳 DOI，可独立核验优先权）：RH 篇 <a href="https://zenodo.org/records/22113226">')
        assert "陈倬" not in s
    else:
        rep(' The proof and construction presented here are by <strong>Zhuo Chen (陈倬)</strong>, archived on Zenodo in August 2026 with timestamped DOIs establishing priority: RH proof <a href="https://zenodo.org/records/22113226">',
            ' The work presented here is by <strong>Zhuo Chen</strong>, archived on Zenodo in August 2026 with timestamped DOIs establishing priority: RH paper <a href="https://zenodo.org/records/22113226">')
        assert "陈倬" not in s

    # 2. font stack: serif -> system sans
    if fn == "zh.html":
        rep('body{font-family:Georgia,"Times New Roman","Songti SC","SimSun",serif;',
            'body{font-family:"Segoe UI","Helvetica Neue",Arial,"Microsoft YaHei","PingFang SC","Hiragino Sans GB",sans-serif;')
        rep('h1,h2,h3{font-family:Georgia,"Songti SC",serif;',
            'h1,h2,h3{font-family:"Segoe UI","Helvetica Neue",Arial,"Microsoft YaHei","PingFang SC","Hiragino Sans GB",sans-serif;')
    else:
        rep('body{font-family:Georgia,"Times New Roman",serif;',
            'body{font-family:"Segoe UI","Helvetica Neue",Arial,"Microsoft YaHei","PingFang SC",sans-serif;')
        rep('h1,h2,h3{font-family:Georgia,serif;',
            'h1,h2,h3{font-family:"Segoe UI","Helvetica Neue",Arial,"Microsoft YaHei","PingFang SC",sans-serif;')

    # 2b. inline MathJax baseline nudge
    rep('  .formula{text-align:center;margin:1.2em 0;overflow-x:auto;padding:.3em 0}',
        '  .formula{text-align:center;margin:1.2em 0;overflow-x:auto;padding:.3em 0}\n  mjx-container[display="false"]{vertical-align:-0.06em}\n  mjx-container[display="true"]{overflow-x:auto;overflow-y:hidden;max-width:100%}')

    # 3. legend: four colours
    if fn == "zh.html":
        rep('；灰色 = 仍在收敛；&mdash; = 尚未出现。',
            '；<span style="background:#fff8e1;color:#856404;font-weight:700;padding:.05em .35em;border-radius:3px">金色</span> = 临近锁定；<span style="color:#999;font-weight:700;padding:.05em .2em;border-radius:3px">灰色</span> = 仍在收敛；<span style="background:#f0f0f0;color:#aaa;padding:.05em .35em;border-radius:3px">&mdash;</span> = 尚未出现。')
    else:
        rep('; grey = still converging; &mdash; = not yet present.',
            '; <span style="background:#fff8e1;color:#856404;font-weight:700;padding:.05em .35em;border-radius:3px">gold</span> = near-locked; <span style="color:#999;font-weight:700;padding:.05em .2em;border-radius:3px">grey</span> = still converging; <span style="background:#f0f0f0;color:#aaa;padding:.05em .35em;border-radius:3px">&mdash;</span> = not yet present.')

    # 4a. touch-action on horizontal scrollers
    rep('.table-scroll{overflow-x:auto;-webkit-overflow-scrolling:touch;',
        '.table-scroll{overflow-x:auto;-webkit-overflow-scrolling:touch;touch-action:pan-x pan-y;')
    rep('.img-scroll{overflow-x:auto;-webkit-overflow-scrolling:touch;',
        '.img-scroll{overflow-x:auto;-webkit-overflow-scrolling:touch;touch-action:pan-x pan-y;')

    # 4b. matrix SVGs: wider min-width for meaningful sideways scroll on phones
    rep('style="width:100%;min-width:520px;max-width:760px;height:auto;',
        'style="width:100%;min-width:760px;max-width:760px;height:auto;', n=2)

    # 4c. tap-a-cell tooltip JS
    JS = (
        '<script>\n'
        '(function(){\n'
        '  var tip=document.createElement("div");\n'
        '  tip.style.cssText="display:none;position:fixed;z-index:9999;max-width:72vw;padding:.5em .8em;background:#222;color:#fff;font-size:.85em;line-height:1.4;border-radius:6px;box-shadow:0 2px 10px rgba(0,0,0,.35);word-break:break-word;pointer-events:none;font-family:Consolas,Menlo,monospace;";\n'
        '  document.body.appendChild(tip);\n'
        '  function hide(){tip.style.display="none";}\n'
        '  document.addEventListener("click",function(e){\n'
        '    var r=e.target&&e.target.closest?e.target.closest("rect"):null;\n'
        '    if(r){var t=r.querySelector("title");if(t&&t.textContent){tip.textContent=t.textContent;tip.style.display="block";var w=window.innerWidth||document.documentElement.clientWidth;var x=e.clientX,y=e.clientY;tip.style.left=Math.max(8,Math.min(x+12,w-tip.offsetWidth-12))+"px";tip.style.top=Math.max(8,y+14)+"px";return;}}\n'
        '    hide();\n'
        '  });\n'
        '  window.addEventListener("scroll",hide,{passive:true});\n'
        '})();\n'
        '</script>\n'
    )
    rep('</body>', JS + '</body>')

    with io.open(fn, "w", encoding="utf-8") as f:
        f.write(s)
    print(fn, "patched, delta:", len(s) - len(orig))

for fn in FILES:
    patch(fn)
print("ALL PATCHED")
