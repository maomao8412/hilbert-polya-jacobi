# -*- coding: utf-8 -*-
"""Round 9: (1) 首屏副标题去掉"证明/proof"，黎曼保留（黎曼零点）；(2) 版权块署名加回中文名 陈倬。
正文（§2 起）所有"证明/proof"表述不动。"""
import io

def patch(fn, pairs):
    with io.open(fn, encoding='utf-8') as f:
        s = f.read()
    for old, new, expect in pairs:
        c = s.count(old)
        assert c == expect, f"{fn}: anchor count {c} != {expect} for: {old[:60]}"
        s = s.replace(old, new)
    with io.open(fn, 'w', encoding='utf-8') as f:
        f.write(s)
    print(f"{fn}: OK, {len(pairs)} patches applied")

# ---- 中文 ----
patch('zh.html', [
    # 1) 副标题：去"证明"，改为 长出矩阵->解出黎曼零点->解码素数
    ('<p class="subtitle">一个高中生都认识的三角恒等式，如何长出黎曼猜想的证明和一台解码素数的矩阵机器</p>',
     '<p class="subtitle">一个高中生都认识的三角恒等式，如何长出矩阵、解出黎曼零点，并解码素数</p>', 1),
    # 2) 版权块署名加回中文名
    ('本页所述工作由 <strong>Zhuo Chen</strong> 完成',
     '本页所述工作由 <strong>陈倬（Zhuo Chen）</strong> 完成', 1),
])

# ---- 英文 ----
patch('index.html', [
    # 1) 副标题：去 proof
    ('<p class="subtitle">How a high-school trigonometric identity grows into a proof of the Riemann hypothesis and an explicit Jacobi matrix whose eigenvalues decode the primes</p>',
     '<p class="subtitle">How a high-school trigonometric identity grows into a matrix that yields the Riemann zeros &mdash; and decodes the primes</p>', 1),
    # 2) 版权块署名加中文名
    ('The work presented here is by <strong>Zhuo Chen</strong>',
     'The work presented here is by <strong>Zhuo Chen (陈倬)</strong>', 1),
])
print("ROUND 9 DONE")
