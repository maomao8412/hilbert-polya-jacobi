# -*- coding: utf-8 -*-
"""
Round-2 webpage fixes for zh.html / index.html (Hilbert-Polya-Jacobi Pages).
A2: logic chain rewrite — RH is proved in §4 directly (angular monotonicity ->
    Herglotz -> pole contradiction); D_n>0 is a post-RH theorem serving the
    Jacobi construction; spectrum identification comes later still.
B : Chinese translation fixes (sign row, high-precision table header).
C : wide tables -> horizontal scroll containers, no digit squashing.
D : inline-SVG full 50x50 tridiagonal picture of J50 inserted in §8'.
E : hero chain: Jacobi matrix -> Riemann zeros -> primes.
F : §1 write out the Mobius map explicitly; correct the T-conjugation claim.
All replacements are anchored and asserted count==1 (zh+en independently).
"""
import csv, html, sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
GIT = "/tmp/hpgit"

# ---------------------------------------------------------------- J50 SVG
def build_svg(lang):
    n = 50
    cellsize = 11.0
    pad = 12
    W = pad * 2 + n * cellsize
    with open(os.path.join(GIT, "data", "J50_matrix.csv"), newline="") as f:
        rows = [[float(x) for x in r] for r in csv.reader(f)]
    assert len(rows) == n and all(len(r) == n for r in rows), "J50 csv shape"
    parts = []
    parts.append(
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %.0f %.0f" '
        'style="width:100%%;min-width:520px;max-width:760px;height:auto;'
        'background:#fff;border:1px solid #e0e0e0;border-radius:6px;'
        'display:block;margin:.6em auto">' % (W, W))
    # outer frame
    parts.append('<rect x="%g" y="%g" width="%g" height="%g" fill="none" '
                 'stroke="#888" stroke-width="1"/>'
                 % (pad - 1, pad - 1, n * cellsize + 2, n * cellsize + 2))
    blue = "#1a4d8f"
    gold = "#b8860b"
    for i in range(n):
        for j in range(n):
            v = rows[i][j]
            x = pad + j * cellsize
            y = pad + i * cellsize
            if abs(v) < 1e-300:
                continue
            diag = (i == j)
            fill = blue if diag else gold
            tip = "J[%d,%d] = %.4e" % (i + 1, j + 1, v)
            parts.append(
                '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" '
                'fill="%s"><title>%s</title></rect>'
                % (x + 0.6, y + 0.6, cellsize - 1.2, cellsize - 1.2,
                   fill, html.escape(tip)))
    # index labels 1,10,...,50 along top & left
    for k in [1, 10, 20, 30, 40, 50]:
        cx = pad + (k - 1) * cellsize + cellsize / 2
        parts.append('<text x="%.1f" y="%g" font-size="8" fill="#666" '
                     'text-anchor="middle">%d</text>' % (cx, pad - 4, k))
        cy = pad + (k - 1) * cellsize + cellsize / 2 + 3
        parts.append('<text x="%g" y="%.1f" font-size="8" fill="#666" '
                     'text-anchor="end">%d</text>' % (pad - 4, cy, k))
    parts.append("</svg>")
    svg = "".join(parts)

    if lang == "zh":
        lead = ('<p><strong>⓪ J<sub>50</sub> 的全貌：50×50 三对角结构</strong>'
                '（蓝=主对角 α<sub>n</sub>，金=两条次对角 b<sub>n</sub>；'
                '其余 2352 个位置全是严格的零，图中留白）。手机上可左右滑动、'
                '双指放大；点按（或悬停）任一彩色格子可读出该矩阵元的精确数值：</p>')
        cap = ('<p style="font-size:.9em;color:var(--subtle)">一条主对角线加两条'
               '次对角线——50×50 共 2500 个位置中只有 148 个非零，其余严格为零，'
               '这就是三对角 Jacobi 矩阵的形状；它与 c<sub>n</sub> 二阶递推'
               '（同样只有“前两项”）同构。矩阵长这样，它的谱就是下一节的黎曼零点。</p>')
    else:
        lead = ('<p><strong>(0) The full J<sub>50</sub>: its 50&times;50 '
                'tridiagonal shape</strong> (blue = main diagonal '
                '&alpha;<sub>n</sub>, gold = the two off-diagonals '
                'b<sub>n</sub>; the other 2352 entries are exact zeros, left '
                'blank). Scroll sideways and pinch-zoom on a phone; tap (or '
                'hover) any coloured cell to read the exact value of that '
                'entry:</p>')
        cap = ('<p style="font-size:.9em;color:var(--subtle)">One main diagonal '
               'and two off-diagonals &mdash; only 148 of the 2500 entries are '
               'non-zero, all others exactly zero: the shape of a tridiagonal '
               'Jacobi matrix, isomorphic to the second-order c<sub>n</sub> '
               'recurrence (which likewise keeps only the two preceding '
               'terms). This is what the matrix looks like &mdash; and its '
               'spectrum is the Riemann zeros of the next section.</p>')
    return ('<div class="table-scroll">' + lead + svg + "</div>" + cap)

# ---------------------------------------------------------------- helpers
def rep(text, old, new, tag):
    c = text.count(old)
    assert c == 1, "[%s] anchor count=%d (need 1)" % (tag, c)
    return text.replace(old, new)

def wrap_tables(htmltext):
    """Wrap every bare <table>...</table> (not already in a scroll div) in
    <div class="table-scroll"> ... </div>, stack-style pairing."""
    out = []
    i = 0
    wrapped = 0
    while True:
        t = htmltext.find("<table", i)
        if t == -1:
            out.append(htmltext[i:])
            break
        ct = htmltext.find(">", t)
        te = htmltext.find("</table>", ct)
        assert te != -1, "unclosed <table>"
        te_end = te + len("</table>")
        # already in an overflow-x scroll container?
        before = htmltext[max(0, t - 220):t]
        if "overflow-x:auto" in before or 'class="table-scroll"' in before:
            out.append(htmltext[i:te_end])
        else:
            out.append(htmltext[i:t])
            out.append('<div class="table-scroll">')
            out.append(htmltext[t:te_end])
            out.append("</div>")
            wrapped += 1
        i = te_end
    return "".join(out), wrapped

CSS = (
".table-scroll{overflow-x:auto;-webkit-overflow-scrolling:touch;"
"margin:.8em 0;border:1px solid #eee;border-radius:6px;padding:.4em .5em;"
"background:#fff;}\n"
".table-scroll table{min-width:560px;border-collapse:collapse;margin:0 auto;}\n"
".table-scroll table.conv-table{min-width:900px;}\n"
".table-scroll td,.table-scroll th{white-space:nowrap;}\n"
"</style>"
)

# =============================================================== ZH
zhp = os.path.join(GIT, "zh.html")
zh = open(zhp, encoding="utf-8").read()
z0 = len(zh)

# ---- C: CSS
zh = rep(zh, "</style>", CSS, "zh-css")

# ---- E: hero chain add Riemann zeros
zh = rep(
    zh,
    '<span class="arrow">→</span> Jacobi 矩阵\n<span class="arrow">→</span> 素数',
    '<span class="arrow">→</span> Jacobi 矩阵\n'
    '<span class="arrow">→</span> 黎曼零点\n'
    '<span class="arrow">→</span> 素数',
    "zh-hero")

# ---- F: §1 Mobius explicit + corrected T-conjugation
zh_old_p = ('很少有人注意到，这正是关于 \\(\\tan B\\) 的一个 <strong>Möbius'
            '（分式线性）变换</strong>。取 \\(A=\\pi/8\\)，令 '
            '\\(q=\\tan(\\pi/8)=\\sqrt2-1\\)（银比共轭，方程 '
            '\\(q^2+2q-1=0\\) 的正根，源自 Euler 1752 年的双纽线加法定理）。'
            'T-共轭 \\(\\sqrt2\\mapsto-\\sqrt2\\) 生成的 Möbius 对合产生一对'
            '代数共轭数：</p>')
zh_new_p = ('很少有人注意到，把 \\(\\tan B\\) 记作变量 \\(t\\)，这正是一个'
            '<strong>Möbius（分式线性）变换</strong>。固定 \\(A=\\pi/8\\)，'
            '令 \\(q=\\tan(\\pi/8)=\\sqrt2-1\\)（方程 \\(q^2+2q-1=0\\) 的正根，'
            '源自 Euler 1752 年的双纽线加法定理），加法公式写成显式形式</p>')
zh = rep(zh, zh_old_p, zh_new_p, "zh-mob-p")

zh_old_formula = ('<div class="formula">\n$$z_1=2-\\sqrt2\\approx0.5858,'
                  '\\qquad z_2=\\sqrt2-1\\approx0.4142,\\qquad z_1+z_2=1.$$\n'
                  '</div>')
zh_new_formula = ('<div class="formula">\n$$M(t)=\\tan(B+\\pi/8)'
                  '=\\frac{t+q}{1-qt},\\qquad q=\\sqrt2-1.$$\n</div>\n'
                  '<p>\\(M\\) 是双纽线 Möbius 算术给出的分式线性变换；'
                  '它的银比共轭为 \\(q^{-1}=1+\\sqrt2\\)。把 \\(q\\) 与它的'
                  '互补值 \\(1-q\\) 并置，得到一对<strong>互补对</strong>'
                  '（由对合 \\(z\\mapsto1-z\\) 配对，二者是二次方程 '
                  '\\(x^2-x+q(1-q)=0\\) 的两根）：</p>\n'
                  '<div class="formula">\n$$z_1=1-q=2-\\sqrt2\\approx0.5858,'
                  '\\qquad z_2=q=\\sqrt2-1\\approx0.4142,\\qquad z_1+z_2=1.$$\n'
                  '</div>')
zh = rep(zh, zh_old_formula, zh_new_formula, "zh-mob-formula")

zh_old_back = ('这个互补关系是后续一切的代数主干。Möbius 结构不是任意选择'
               '——它是正切加法公式强制给出的。')
zh_new_back = ('这个互补关系是后续一切的代数主干（论文中变形族的完整 Möbius '
               '对称性，正是由银比共轭 \\(q=\\sqrt2-1\\) 与其互补 \\(1-q\\) '
               '生成）。Möbius 结构不是任意选择——它是正切加法公式强制给出的。')
zh = rep(zh, zh_old_back, zh_new_back, "zh-mob-back")

# ---- A2-1: §4 closing paragraph
zh = rep(
    zh,
    '论证从角单调性直接走到 Herglotz 解析性再到极点矛盾——不经过任何测度或矩阵。'
    '下文的 Jacobi 构造是另一条<strong>独立路线</strong>：它的拱心石 '
    '\\(D_n&gt;0\\) 由同一三项和代数<strong>无条件</strong>证明'
    '（第五、六节），RH 在谱的一端作为结论重新出现，而非作为前提登场。',
    'RH 就在本节证毕：论证从角单调性直接走到 Herglotz 解析性再到极点矛盾'
    '——不经过任何测度或矩阵。下文的 Jacobi 构造是 RH <strong>之后</strong>'
    '的独立构造：它的拱心石 \\(D_n&gt;0\\) 由同一三项和代数'
    '<strong>无条件</strong>证明（第五、六节），是 RH 证毕之后交出的定理，'
    '支撑 Stieltjes 测度与 Jacobi 矩阵的构造；谱被识别为 '
    '\\(\\{1/\\gamma_n^2\\}\\) 更在其后。',
    "zh-a2-132")

# ---- A2-2: "两个方向" box
zh = rep(
    zh,
    '<strong>两个方向。</strong>构造篇无条件证明的是正向——三项和代数 '
    '\\(\\to\\) 角单调性 \\(\\to\\) Herglotz \\(\\to\\) \\(D_n&gt;0\\) '
    '对所有 \\(n\\) 成立——矩序列正定<em>即</em> RH——Hankel 全正定等价于'
    '零点全在临界线上；此后才把谱读作 \\(\\{1/\\gamma_n^2\\}\\)，那是 RH '
    '之后的识别。反方向是经典结论（Hadamard＋Gram）：<em>若</em>全部零点在'
    '临界线上，则 \\(S_k=\\sum_n 1/\\gamma_n^{2k}\\) 使 Hankel 矩阵',
    '<strong>两个方向。</strong>RH 的证明在第四节完成——三项和代数 '
    '\\(\\to\\) 角单调性 \\(\\to\\) Herglotz 解析性 \\(\\to\\) 极点矛盾，'
    '直接推出全部零点在临界线上。RH 之后，同一三项和代数<strong>无条件'
    '</strong>证明 \\(D_n&gt;0\\) 对所有 \\(n\\) 成立（Hankel 行列式 = '
    'Hankel 矩阵的顺序主子式，全部为正），作为 Stieltjes–Jacobi 构造的地基；'
    '谱读作 \\(\\{1/\\gamma_n^2\\}\\) 是更后面的识别。反方向是经典结论'
    '（Hadamard＋Gram）：<em>若</em>全部零点在临界线上，则 '
    '\\(S_k=\\sum_n 1/\\gamma_n^{2k}\\) 使 Hankel 矩阵',
    "zh-a2-146")
zh = rep(
    zh,
    '成为向量 \\((\\gamma_n^{-i-j-1})\\) 的 Gram 矩阵，故每个顺序主子式 '
    '\\(D_n=\\det(H_n)\\) 严格为正。前人从未拿到的是正向；下面的 Gram '
    '图景解释的是反向，不是证书本身。',
    '成为向量 \\((\\gamma_n^{-i-j-1})\\) 的 Gram 矩阵，故每个顺序主子式'
    '（Hankel 行列式）\\(D_n=\\det(H_n)\\) 严格为正。前人从未拿到的是'
    'RH 的直接证明与无条件的 \\(D_n&gt;0\\)；下面的 Gram 图景解释的是'
    '反向，不是 RH 的证明本身。',
    "zh-a2-150")

# ---- B: D_n table caption + sign row
zh = rep(
    zh,
    '<table style="font-size:.78em">\n<tr><th style="text-align:left">\\(n\\)</th>',
    '<p style="font-size:.92em;margin-bottom:.3em">下表：Hankel 矩阵的'
    '顺序主子式 \\(D_n=\\det H_n\\)（又称 Hankel 行列式）——其对数值与符号：</p>\n'
    '<table style="font-size:.78em">\n<tr><th style="text-align:left">\\(n\\)</th>',
    "zh-dn-caption")
zh = rep(zh, '<td style="text-align:left">sign</td>',
         '<td style="text-align:left">符号</td>', "zh-sign")

# ---- B: high precision table header
zh = rep(
    zh,
    '<tr><th>n</th><th>1/√λ<sub>n</sub> from J<sub>50</sub></th>'
    '<th>known γ<sub>n</sub></th><th>rel. error</th></tr>',
    '<tr><th>n</th><th>由 J<sub>50</sub> 算出的 1/√λ<sub>n</sub>'
    '（=γ<sub>n</sub> 估计）</th><th>已知 γ<sub>n</sub></th>'
    '<th>相对误差</th></tr>',
    "zh-hp-header")

# ---- A2-3: §6 closing chain paragraph
zh = rep(
    zh,
    '我们的链条方向相反：<strong>三项和恒等式 → 角单调性 → Herglotz → '
    '所有 \\(D_n&gt;0\\)</strong>（无条件证明）——矩序列正定<em>本身就是 '
    'RH</em>：Hankel 矩阵全部正定，等价于全部零点落在临界线上。RH 在这里'
    '落地，不在开头，也不在矩阵的谱上。此后矩阵才作为<em>推论</em>出场：'
    '正定性交出 Stieltjes 测度，Gram–Schmidt 造出 Jacobi 算子'
    '（Carleman 自伴条件随 \\(q=\\sqrt2-1\\) 白送），其谱再被<em>识别'
    '</em>为 \\(\\{1/\\gamma_n^2\\}\\)——这是 RH 之后的确认，不是证明 RH '
    '的箭头（谱识别用的 Hadamard 乘积本就以零点全在实轴为前提）。'
    '矩阵正是从三项和恒等式里长出来的——这就是为什么没有它，再多的经典'
    '机器也造不出这个矩阵。',
    '我们的链条方向相反：<strong>三项和恒等式 → 角单调性 → Herglotz 解析性'
    ' → 极点矛盾</strong>，RH 在第四节由此直接证毕（全部零点落在临界线上）；'
    '<strong>此后</strong>才轮到矩量正定性——同一三项和代数<em>无条件'
    '</em>证明所有 \\(D_n&gt;0\\)，作为 Jacobi 构造的地基：正定性交出 '
    'Stieltjes 测度，Gram–Schmidt 造出 Jacobi 算子（Carleman 自伴条件随 '
    '\\(q=\\sqrt2-1\\) 白送），其谱再被<em>识别</em>为 '
    '\\(\\{1/\\gamma_n^2\\}\\)（谱识别用的 Hadamard 乘积本就以零点全在'
    '实轴为前提，所以它只能在 RH 之后）。矩阵正是从三项和恒等式里长出来的'
    '——这就是为什么没有它，再多的经典机器也造不出这个矩阵。',
    "zh-a2-192")

# ---- D: J50 full SVG before the 7x7 literal block
zh = rep(
    zh,
    '<p><strong>① 前 7 阶的字面矩阵</strong>',
    build_svg("zh") + '\n<p><strong>① 前 7 阶的字面矩阵</strong>',
    "zh-svg")

# ---- C: wrap all bare tables
zh, nw = wrap_tables(zh)

open(zhp, "w", encoding="utf-8").write(zh)
print("zh.html: %d -> %d bytes, %d tables wrapped" % (z0, len(zh), nw))

# =============================================================== EN
enp = os.path.join(GIT, "index.html")
en = open(enp, encoding="utf-8").read()
e0 = len(en)

en = rep(en, "</style>", CSS, "en-css")

en = rep(
    en,
    '<span class="arrow">→</span> Jacobi matrix\n'
    '<span class="arrow">→</span> primes',
    '<span class="arrow">→</span> Jacobi matrix\n'
    '<span class="arrow">→</span> Riemann zeros\n'
    '<span class="arrow">→</span> primes',
    "en-hero")

en_old_p = ('What is less commonly noticed is that this is a <strong>Möbius'
            ' (fractional-linear) transformation</strong> in \\(\\tan B\\). '
            'Fix \\(A=\\pi/8\\) and set \\(q=\\tan(\\pi/8)=\\sqrt2-1\\), the '
            'silver-ratio conjugate (positive root of \\(q^2+2q-1=0\\), '
            "arising from Euler's lemniscate addition theorem of 1752). "
            'The Möbius involution generated by the T-conjugation '
            '\\(\\sqrt2\\mapsto-\\sqrt2\\) produces an algebraic conjugate '
            'pair</p>')
en_new_p = ('What is less commonly noticed is that, writing \\(t=\\tan B\\), '
            'this is a <strong>Möbius (fractional-linear) transformation'
            '</strong>. Fix \\(A=\\pi/8\\) and set '
            '\\(q=\\tan(\\pi/8)=\\sqrt2-1\\) (positive root of '
            '\\(q^2+2q-1=0\\), arising from Euler’s lemniscate addition '
            'theorem of 1752); the addition formula is the explicit map</p>')
en = rep(en, en_old_p, en_new_p, "en-mob-p")

en_old_formula = ('<div class="formula">\n$$z_1=2-\\sqrt2\\approx0.5858,'
                  '\\qquad z_2=\\sqrt2-1\\approx0.4142,\\qquad z_1+z_2=1.$$\n'
                  '</div>')
en_new_formula = ('<div class="formula">\n$$M(t)=\\tan(B+\\pi/8)'
                  '=\\frac{t+q}{1-qt},\\qquad q=\\sqrt2-1.$$\n</div>\n'
                  '<p>\\(M\\) is the fractional-linear map of the lemniscatic '
                  'Möbius arithmetic; its silver-ratio conjugate is '
                  '\\(q^{-1}=1+\\sqrt2\\). Placing \\(q\\) beside its '
                  'complement \\(1-q\\) gives a <strong>complementary pair'
                  '</strong> (paired by the involution \\(z\\mapsto1-z\\); '
                  'the two are the roots of \\(x^2-x+q(1-q)=0\\)):</p>\n'
                  '<div class="formula">\n$$z_1=1-q=2-\\sqrt2\\approx0.5858,'
                  '\\qquad z_2=q=\\sqrt2-1\\approx0.4142,\\qquad z_1+z_2=1.$$\n'
                  '</div>')
en = rep(en, en_old_formula, en_new_formula, "en-mob-formula")

en_old_back = ('This complement relation is the algebraic backbone of '
               'everything that follows. The Möbius structure is not an '
               'arbitrary choice — it is forced by the tangent addition '
               'formula.')
en_new_back = ('This complement relation is the algebraic backbone of '
               'everything that follows (the full Möbius symmetry of the '
               'deformation family is generated precisely by the '
               'silver-ratio conjugate \\(q=\\sqrt2-1\\) and its complement '
               '\\(1-q\\)). The Möbius structure is not an arbitrary choice '
               '— it is forced by the tangent addition formula.')
en = rep(en, en_old_back, en_new_back, "en-mob-back")

en = rep(
    en,
    'The argument routes directly from angular monotonicity to Herglotz '
    'analyticity to the pole contradiction — it does not pass through any '
    'measure or matrix. The Jacobi construction below is a separate, '
    'independent route: its keystone \\(D_n&gt;0\\) is proved '
    '<strong>unconditionally</strong> from the same three-term algebra '
    '(Sections 5–6), and RH reappears as its spectral conclusion rather '
    'than standing as its premise.',
    'RH is proved in this very section: the argument routes directly from '
    'angular monotonicity to Herglotz analyticity to the pole contradiction '
    '— it does not pass through any measure or matrix. The Jacobi '
    'construction below is an independent construction <strong>after '
    'RH</strong>: its keystone \\(D_n&gt;0\\) is proved '
    '<strong>unconditionally</strong> from the same three-term algebra '
    '(Sections 5–6) as a theorem delivered once RH is established; it '
    'supports the Stieltjes measure and the Jacobi matrix, and the '
    'identification of the spectrum as \\(\\{1/\\gamma_n^2\\}\\) comes later '
    'still.',
    "en-a2-135")

en = rep(
    en,
    '<strong>The two directions.</strong> The construction paper proves the '
    'forward implication <strong>unconditionally</strong> &mdash; three-term '
    'algebra \\(\\to\\) angular monotonicity \\(\\to\\) Herglotz '
    '\\(\\to\\) \\(D_n&gt;0\\) for all \\(n\\) ; positivity of the moment '
    'sequence <em>is</em> RH &mdash; Hankel positivity is equivalent to '
    'every zero lying on the line &mdash; and only afterwards is the '
    'spectrum read off as \\(\\{1/\\gamma_n^2\\}\\), an identification that '
    'RH itself licenses. The converse is classical (Hadamard + Gram): '
    '<em>if</em> all zeros lie on the line, then '
    '\\(S_k=\\sum_n 1/\\gamma_n^{2k}\\) makes the Hankel matrices',
    '<strong>The two directions.</strong> RH is proved in Section 4 — '
    'three-term algebra \\(\\to\\) angular monotonicity \\(\\to\\) Herglotz '
    'analyticity \\(\\to\\) the pole contradiction, which places every zero '
    'on the critical line. After RH, the same three-term algebra proves '
    '\\(D_n&gt;0\\) for all \\(n\\) <strong>unconditionally</strong> '
    '(the Hankel determinants — the leading principal minors — all '
    'positive), as the foundation of the Stieltjes–Jacobi construction; '
    'reading the spectrum as \\(\\{1/\\gamma_n^2\\}\\) is a still later '
    'identification. The converse is classical (Hadamard + Gram): '
    '<em>if</em> all zeros lie on the line, then '
    '\\(S_k=\\sum_n 1/\\gamma_n^{2k}\\) makes the Hankel matrices',
    "en-a2-149")
en = rep(
    en,
    'Gram matrices of the vectors \\((\\gamma_n^{-i-j-1})\\), so every '
    'leading principal minor \\(D_n=\\det(H_n)\\) is strictly positive. The '
    'forward direction is what was never available before; the Gram picture '
    'explains the converse, it is not the certificate.',
    'Gram matrices of the vectors \\((\\gamma_n^{-i-j-1})\\), so every '
    'leading principal minor (Hankel determinant) \\(D_n=\\det(H_n)\\) is '
    'strictly positive. What was never available before is the direct proof '
    'of RH and the unconditional \\(D_n&gt;0\\); the Gram picture explains '
    'the converse — it is not the proof of RH itself.',
    "en-a2-150")

en = rep(
    en,
    'Our chain runs the other way: <strong>three-term identity &rarr; '
    'angular monotonicity &rarr; Herglotz &rarr; every \\(D_n&gt;0\\)'
    '</strong>, proved unconditionally &mdash; a positive moment sequence '
    '<em>is</em> RH: positivity of all Hankel matrices is equivalent to '
    'every zero lying on the critical line. RH lands here &mdash; not at '
    'the beginning of the construction, and not at the matrix spectrum. The '
    'matrix comes only afterwards, as a <em>corollary</em>: positivity '
    'hands over the Stieltjes measure, Gram–Schmidt produces the Jacobi '
    'operator (self-adjoint; Carleman’s condition comes for free with '
    '\\(q=\\sqrt2-1\\)), and its spectrum is then <em>identified</em> as '
    '\\(\\{1/\\gamma_n^2\\}\\) &mdash; a confirmation after RH, not an '
    'arrow that proves it (the Hadamard product used in the identification '
    'presupposes real zeros). That is why the matrix grew out of the '
    'three-term identity &mdash; and why, without it, no amount of '
    'classical machinery could produce it.',
    'Our chain runs the other way: <strong>three-term identity &rarr; '
    'angular monotonicity &rarr; Herglotz analyticity &rarr; the pole '
    'contradiction</strong>, by which RH is proved directly in Section 4 '
    '(every zero lies on the critical line); <strong>only then</strong> '
    'does moment positivity take its place — the same three-term algebra '
    'proves every \\(D_n&gt;0\\) <em>unconditionally</em> as the foundation '
    'of the Jacobi construction: positivity hands over the Stieltjes '
    'measure, Gram–Schmidt produces the Jacobi operator (self-adjoint; '
    'Carleman’s condition comes for free with \\(q=\\sqrt2-1\\)), and its '
    'spectrum is then <em>identified</em> as \\(\\{1/\\gamma_n^2\\}\\) (the '
    'Hadamard product used in that identification presupposes real zeros, '
    'so it can only come after RH). That is why the matrix grew out of the '
    'three-term identity &mdash; and why, without it, no amount of '
    'classical machinery could produce it.',
    "en-a2-195")

en = rep(
    en,
    '<p><strong>(1) The literal first seven orders</strong>',
    build_svg("en") + '\n<p><strong>(1) The literal first seven orders'
    '</strong>',
    "en-svg")

en, nw2 = wrap_tables(en)

open(enp, "w", encoding="utf-8").write(en)
print("index.html: %d -> %d bytes, %d tables wrapped" % (e0, len(en), nw2))
print("ROUND2 OK")
