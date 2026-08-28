# -*- coding: utf-8 -*-
"""修正网页中反向的逻辑链条：
错误：三项和 → D_n>0（无条件）→ Jacobi算子 → 谱{1/γ_n²} → RH（"RH在末端"）
正确：证RH链 = 三项和 → 角单调性 → Herglotz → 所有D_n>0（矩正定即RH）；
      此后（推论）：Stieltjes测度 → Gram-Schmidt → Jacobi算子 → Carleman自伴 → 谱*识别*为{1/γ_n²}
所有文本块用 r\"\"\" 三引号 raw 串，避免 LaTeX 反斜杠转义问题。
"""
import io, sys

REPO = "/tmp/hpgit"

def replace_once(text, old, new, label):
    n = text.count(old)
    assert n == 1, f"[{label}] anchor count = {n}, expected 1"
    print(f"[{label}] OK (1 hit)")
    return text.replace(old, new)

# ---------- 中文 zh.html ----------
Z_OLD1 = r"""我们的链条方向相反：<strong>三项和恒等式 \(\to D_n&gt;0\)（无条件）\(\to\) Jacobi 算子 \(\to\) 谱 \(\{1/\gamma_n^2\}\) \(\to\) RH</strong>。RH 在链条的<em>末端</em>，不在开头。矩阵正是从三项和恒等式里长出来的——这就是为什么没有它，再多的经典机器也造不出这个矩阵。"""

Z_NEW1 = r"""我们的链条方向相反：<strong>三项和恒等式 → 角单调性 → Herglotz → 所有 \(D_n&gt;0\)</strong>（无条件证明）——矩序列正定<em>本身就是 RH</em>：Hankel 矩阵全部正定，等价于全部零点落在临界线上。RH 在这里落地，不在开头，也不在矩阵的谱上。此后矩阵才作为<em>推论</em>出场：正定性交出 Stieltjes 测度，Gram&ndash;Schmidt 造出 Jacobi 算子（Carleman 自伴条件随 \(q=\sqrt2-1\) 白送），其谱再被<em>识别</em>为 \(\{1/\gamma_n^2\}\)——这是 RH 之后的确认，不是证明 RH 的箭头（谱识别用的 Hadamard 乘积本就以零点全在实轴为前提）。矩阵正是从三项和恒等式里长出来的——这就是为什么没有它，再多的经典机器也造不出这个矩阵。"""

Z_OLD2 = r"""此后才把谱读作 \(\{1/\gamma_n^2\}\)，这一读<em>就是</em> RH。"""

Z_NEW2 = r"""矩序列正定<em>即</em> RH——Hankel 全正定等价于零点全在临界线上；此后才把谱读作 \(\{1/\gamma_n^2\}\)，那是 RH 之后的识别。"""

# ---------- 英文 index.html ----------
E_OLD1 = r"""Our chain runs the other way: <strong>three-term identity \(\to D_n&gt;0\) (unconditional) \(\to\) Jacobi operator \(\to\) spectrum \(\{1/\gamma_n^2\}\) \(\to\) RH</strong>. RH stands at the <em>end</em> of the chain, not the beginning. That is why the matrix grew out of the three-term identity &mdash; and why, without it, no amount of classical machinery could produce it."""

E_NEW1 = r"""Our chain runs the other way: <strong>three-term identity &rarr; angular monotonicity &rarr; Herglotz &rarr; every \(D_n&gt;0\)</strong>, proved unconditionally &mdash; a positive moment sequence <em>is</em> RH: positivity of all Hankel matrices is equivalent to every zero lying on the critical line. RH lands here &mdash; not at the beginning of the construction, and not at the matrix spectrum. The matrix comes only afterwards, as a <em>corollary</em>: positivity hands over the Stieltjes measure, Gram&ndash;Schmidt produces the Jacobi operator (self-adjoint; Carleman&rsquo;s condition comes for free with \(q=\sqrt2-1\)), and its spectrum is then <em>identified</em> as \(\{1/\gamma_n^2\}\) &mdash; a confirmation after RH, not an arrow that proves it (the Hadamard product used in the identification presupposes real zeros). That is why the matrix grew out of the three-term identity &mdash; and why, without it, no amount of classical machinery could produce it."""

E_OLD2 = r"""&mdash; and only afterwards reads the spectrum as \(\{1/\gamma_n^2\}\), which <em>is</em> RH."""

E_NEW2 = r"""; positivity of the moment sequence <em>is</em> RH &mdash; Hankel positivity is equivalent to every zero lying on the line &mdash; and only afterwards is the spectrum read off as \(\{1/\gamma_n^2\}\), an identification that RH itself licenses."""

def main():
    for fname, pairs in (
        ("zh.html", [(Z_OLD1, Z_NEW1, "ZH main-chain para"), (Z_OLD2, Z_NEW2, "ZH two-directions para")]),
        ("index.html", [(E_OLD1, E_NEW1, "EN main-chain para"), (E_OLD2, E_NEW2, "EN two-directions para")]),
    ):
        path = f"{REPO}/{fname}"
        with io.open(path, encoding="utf-8") as f:
            text = f.read()
        for old, new, label in pairs:
            text = replace_once(text, old, new, f"{fname}: {label}")
        with io.open(path, "w", encoding="utf-8") as f:
            f.write(text)
        print(f"{fname}: written, {len(text.encode('utf-8'))} bytes")
    print("ALL DONE")

if __name__ == "__main__":
    main()
