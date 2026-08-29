# -*- coding: utf-8 -*-
"""Round-4 fixes:
1. All subtle notes / figcaptions centered.
2. Section 9 PNGs wrapped in .img-scroll (phone horizontal scroll; tap opens full-size PNG).
3. Real 50x50 J50 tridiagonal data table (148 non-zero entries) inserted after tri-SVG caption.
"""
import csv, re

def center_notes(html):
    a = '<p style="font-size:.9em;color:var(--subtle)">'
    b = '<p style="font-size:.9em;color:var(--subtle);text-align:center">'
    n1 = html.count(a); html = html.replace(a, b)
    c = '<figcaption style="font-size:.88em;color:var(--subtle);margin-top:.4em">'
    d = '<figcaption style="font-size:.88em;color:var(--subtle);margin-top:.4em;text-align:center">'
    n2 = html.count(c); html = html.replace(c, d)
    assert n1 == 7 and n2 == 5, (n1, n2)
    return html

IMG_RE = re.compile(r'<img src="(data/[^"]+)" alt="([^"]*)" style="width:100%;max-width:1000px;border:1px solid #e0e0e0;border-radius:6px"\s*/>')

def wrap_images(html, hint):
    count = {'n': 0}
    def rep(m):
        count['n'] += 1
        src, alt = m.group(1), m.group(2)
        return ('<div class="img-scroll"><a href="{src}" target="_blank" rel="noopener" '
                'title="{hint}"><img src="{src}" alt="{alt}" '
                'style="border:1px solid #e0e0e0;border-radius:6px;cursor:zoom-in"/></a></div>'
               ).format(src=src, alt=alt, hint=hint)
    html = IMG_RE.sub(rep, html)
    assert count['n'] == 5, count
    return html

def add_caption_hint(html, hint):
    parts = html.split('</figcaption>')
    assert len(parts) == 6, len(parts)
    return (' ' + hint + '</figcaption>').join(parts)

CSS_ADD = """
.img-scroll{overflow-x:auto;-webkit-overflow-scrolling:touch;margin:1.2em auto;max-width:1000px;border:1px solid #eee;border-radius:6px;padding:.4em;background:#fff;}
.img-scroll img{display:block;width:1000px;max-width:none;height:auto;margin:0 auto;}
table.j50-data{min-width:3400px;font-size:.62em;border-collapse:collapse;font-variant-numeric:tabular-nums;margin:0 auto;}
table.j50-data th,table.j50-data td{border:1px solid var(--rule);padding:.18em .3em;text-align:center;white-space:nowrap;}
table.j50-data th{background:#eef2f8;color:var(--heading);}
table.j50-data td.diag{background:#f0f4fa;color:#1a4d8f;font-weight:bold;}
table.j50-data td.off{background:#fbf6e6;color:#8a6d1a;}
table.j50-data td.rl,table.j50-data th.rl{background:#dde6f2;font-weight:bold;position:sticky;left:0;z-index:2}
.table-scroll table.j50-data{min-width:3400px;}
"""

def add_css(html):
    anchor = '.table-scroll td,.table-scroll th{white-space:nowrap;}'
    assert html.count(anchor) == 1
    return html.replace(anchor, anchor + CSS_ADD)

def build_j50_table(corner):
    rows = []
    with open('data/J50_matrix.csv', newline='') as f:
        for row in csv.reader(f):
            rows.append([float(x) for x in row])
    assert len(rows) == 50 and all(len(x) == 50 for x in rows)
    nz = sum(1 for i in range(50) for j in range(50) if rows[i][j] != 0.0)
    assert nz == 148, nz
    def fmt(v):
        s = '%.3e' % v
        return s.replace('e+0', 'e+').replace('e-0', 'e-')
    out = ['<div class="table-scroll"><table class="j50-data">']
    out.append('<tr><th class="rl">' + corner + '</th>'
               + ''.join('<th>%d</th>' % (j + 1) for j in range(50)) + '</tr>')
    for i in range(50):
        tds = []
        for j in range(50):
            v = rows[i][j]
            if v == 0.0:
                tds.append('<td></td>')
            elif i == j:
                tds.append('<td class="diag">%s</td>' % fmt(v))
            else:
                tds.append('<td class="off">%s</td>' % fmt(v))
        out.append('<tr><td class="rl">%d</td>%s</tr>' % (i + 1, ''.join(tds)))
    out.append('</table></div>')
    return ''.join(out)

def insert_j50(html, lead, anchor, corner):
    assert html.count(anchor) == 1, anchor
    block = ('<p style="font-size:.9em;color:var(--subtle);text-align:center">' + lead + '</p>'
             + build_j50_table(corner))
    return html.replace(anchor, anchor + block)

# ---------------- zh ----------------
zh = open('zh.html', encoding='utf-8').read()
zh = center_notes(zh)
zh = wrap_images(zh, '点击在新标签页打开全尺寸图片，可双指放大查看细节')
zh = add_caption_hint(zh, '（点按图片可打开全尺寸大图、双指放大查看细节）')
zh = add_css(zh)
zh = insert_j50(
    zh,
    ('下表给出 J<sub>50</sub> 的全部 2500 个矩阵元的正算数值：只有主对角 50 个 '
     '&alpha;<sub>n</sub>（蓝）与两条次对角 98 个 b<sub>n</sub>（金）非零，共 148 个；'
     '其余 2352 格是严格的零，留空。手机上可左右滑动查看整表、双指放大：'),
    '它的谱就是下一节的黎曼零点。</p>',
    '行\\列')
open('zh.html', 'w', encoding='utf-8').write(zh)
print('zh.html done, size', len(zh))

# ---------------- en ----------------
en = open('index.html', encoding='utf-8').read()
en = center_notes(en)
en = wrap_images(en, 'Tap to open the full-size image in a new tab; pinch to zoom for details')
en = add_caption_hint(en, '(Tap the image to open the full-size picture and pinch-zoom for detail.)')
en = add_css(en)
en = insert_j50(
    en,
    ('The table below gives the forward-computed values of all 2500 entries of '
     'J<sub>50</sub>: only the 50 main-diagonal &alpha;<sub>n</sub> (blue) and the 98 '
     'off-diagonal b<sub>n</sub> (gold) are non-zero &mdash; 148 entries in all; the other '
     '2352 cells are exact zeros, left blank. Scroll sideways on a phone and pinch to zoom:'),
    'its spectrum is the Riemann zeros of the next section.</p>',
    'r\\c')
open('index.html', 'w', encoding='utf-8').write(en)
print('index.html done, size', len(en))
