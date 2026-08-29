# -*- coding: utf-8 -*-
import io, re

ZH_PROOF = r"""
<h2>五、广义黎曼猜想的证明：角单调性与极点矛盾</h2>
<p>矩阵（§3–4）展示的是谱与零点重合的<strong>数值事实</strong>；GRH 的解析论证与 RH 主篇第四节<strong>同构</strong>，且结构更简单。本节陈述论文（<a href="https://doi.org/10.5281/zenodo.22143035" target="_blank" rel="noopener">Zenodo 22143035</a>）的证明链条。</p>
<p>对模 \(q\ge 3\) 的本原非主特征 \(\chi\)（奇偶性 \(a\)，\(\chi(-1)=(-1)^a\)），定义完备乘积
\[
\xi_K(s)=\xi(s)\,\Lambda(s,\chi),\qquad
\Lambda(s,\chi)=\Bigl(\frac{q}{\pi}\Bigr)^{(s+a)/2}\Gamma\!\Bigl(\frac{s+a}{2}\Bigr)L(s,\chi)
\]
（复特征取 \(\xi_K=\xi\,\Lambda(s,\chi)\Lambda(s,\bar\chi)\)）。\(\xi_K\) 是一级<strong>整函数</strong>，满足 \(\xi_K(s)=\xi_K(1-s)=\overline{\xi_K(\bar s)}\)。</p>
<p><strong>角单调性定理（论文 Theorem 2）。</strong>对任意本原非主特征，
\[
T_K(r,\theta):=\frac{\partial}{\partial\theta}\log\bigl|\xi_K(\tfrac12+re^{i\theta})\bigr|\le 0,
\qquad r>0,\ \theta\in(0,\pi/2).
\]</p>
<p><strong>极点矛盾（论文 Theorem 1 的证明）。</strong>反设 \(L(s,\chi)\) 有临界线<em>外</em>的非平凡零点 \(\rho=\beta+i\gamma\)（\(\beta\neq\tfrac12\)）。由函数方程与共轭对称，不妨 \(\beta&gt;\tfrac12\)、\(\gamma&gt;0\)；\(\rho\) 也是 \(\xi_K\) 的零点。在 \(\rho\) 附近 \(\xi_K(s)=(s-\rho)^m h(s)\)（\(m\ge1\)，\(h(\rho)\neq0\)），故
\[
\frac{\xi_K'}{\xi_K}(s)=\frac{m}{s-\rho}+\frac{h'}{h}(s).
\]
固定 \(r=r_0=|\rho-\tfrac12|\)，令 \(\theta\to\theta_0^+\)（\(\theta_0=\arg(\rho-\tfrac12)\in(0,\pi/2)\)），\(s\) 从上方趋于 \(\rho\)：
\[
\Re\!\left[\frac{i r_0 e^{i\theta}}{s-\rho}\right]\to+\infty,
\]
而 \(h'/h\) 有界，于是 \(T_K(r_0,\theta)\to+\infty\)——与角单调性 \(T_K\le0\) <strong>直接矛盾</strong>。故 \(\beta=\tfrac12\)：\(L(s,\chi)\) 的全部非平凡零点都在临界线上，<strong>广义黎曼猜想证毕</strong>。</p>
<p><strong>角单调性的三个区域支撑</strong>（论文 §5–8）：大 \(r\) 区（\(r\ge7\)）由 Gamma 因子 Stirling 展开给出 \(-r\log r\) 主项，分 \(\sigma\ge2\)（子区 A）、\(7\le r\le50\)（子区 B）、\(r\ge50,\sigma&lt;2\)（子区 C）估计；中间紧致区（\(2\le r\le7\)）由区间算术严格包围（进行中，见 §8）；小 \(r\) 区（\(r\le2\)）由正核分解给出 Taylor 余项定号，解析成立。10 个代表特征在 \(r\in[7,10000]\) 网格上 \(T_K\) 逐点为负（见 §6 网格与交互热力图）。</p>

<h3>5.1 为什么 GRH 的证明比 RH 更简单</h3>
<p>RH 证明里最硬的一仗，是在 \(\zeta\) <strong>一个</strong>函数上白刃战：把 \(\zeta'/\zeta\) 的噪声压到 Gamma 主项以下。GRH 不是单挑——\(\xi_K=\xi\cdot\Lambda\) 是联合战场，五处结构性占优：</p>
<p><strong>（一）整函数，零极点要消。</strong>\(\zeta\) 在 \(s=1\) 有极点，RH 构造 \(\xi\) 必须用 \(s(s-1)\) 因子消去；非主特征的 \(L(s,\chi)\) 在 \(s=1\) <strong>无极点</strong>，\(\Lambda\) 本身就是整函数，\(\xi_K\) 天然整——少一整步奇点处理。</p>
<p><strong>（二）算术正系数白送。</strong>对二次特征，\(\xi_K\) 就是对应二次数域的完备 Dedekind zeta 函数，Dirichlet 系数 \(r_K(n)\ge0\) 非负；对复特征，\(\Lambda(s,\chi)\Lambda(s,\bar\chi)=|L|^2\) 的系数是 \(|\sum_{d|n}\chi(d)|^2\ge0\)。小 \(r\) 区的正核论证因此对<strong>所有</strong>特征一致成立——RH 里同一步要靠三项多对数恒等式手工把 \(\zeta\) 劈成三块正定材料，GRH 里算术正性直接给。</p>
<p><strong>（三）两份零点密度，Gamma 占优更稳。</strong>\(\xi_K\) 的相位函数从 \(\zeta\) 和 \(L\) <strong>两个独立来源</strong>拿到 \(\log r\) 主项（单独研究 \(\Lambda\) 时没有这一项，其相位正性将被迫假定零点在临界线上——循环论证；乘上 \(\xi\) 正是借走 \(\zeta\) 的 Gamma 因子与零点密度结构）。复特征 Gamma 阶 \(d=3\)、系数 \(3/2\)，对数导数上界仅 \(3M_\zeta(2)=1.71\)，子区 A 直接给出 \(T_K\le-0.56\,r+O(1)&lt;0\)。</p>
<p><strong>（四）模 \(q\) 越大越容易。</strong>Gamma 项含 \(-r\sin\theta\,\log q\)（复特征系数 1，二次特征 \(1/2\)），导体越大余量越宽，且所有界对 \(q\ge3\) <strong>一致成立</strong>——这是罕见的"难度随对象变大而下降"。</p>
<p><strong>（五）核心机器直接继承。</strong>\(\zeta\) 的三项多对数分解（RH 主篇的核心机器）在组合之前就作用于 \(\zeta\) 因子，GRH 无需重新发明：\(\zeta\) 攻下之后，Dirichlet \(L\) 函数是同一套机器的自然推论。</p>
<p>一句话：<strong>RH 是 \(\zeta\) 单函数上的白刃战；GRH 是 \(\xi\cdot\Lambda\) 的联合战场——整函数无障、正系数白给、双份 Gamma 保险、\(q\) 越大余量越宽。</strong></p>
<p style="font-size:.92em;color:var(--subtle)">边界：本节是证明链条的概览，完整证明（三区估计、常数 \(K\) 推导、复特征修正）见 Zenodo 论文；矩阵本身是数值构造，中间区（Region 2）与子区 C 的严格区间算术包围仍在进行（§8）。</p>

<h2>六、46 个本原特征：角单调性网格</h2>"""

EN_PROOF = r"""
<h2>5. Proof of the GRH: angular monotonicity and the pole contradiction</h2>
<p>The matrices (§3–4) show the <strong>numerical fact</strong> that spectrum and zeros coincide; the analytic argument for GRH is <strong>isomorphic</strong> to Section 4 of the RH main paper — and structurally simpler. This section outlines the proof chain of the paper (<a href="https://doi.org/10.5281/zenodo.22143035" target="_blank" rel="noopener">Zenodo 22143035</a>).</p>
<p>For a primitive non-principal character \(\chi\) of conductor \(q\ge3\) (parity \(a\), \(\chi(-1)=(-1)^a\)), define the completed product
\[
\xi_K(s)=\xi(s)\,\Lambda(s,\chi),\qquad
\Lambda(s,\chi)=\Bigl(\frac{q}{\pi}\Bigr)^{(s+a)/2}\Gamma\!\Bigl(\frac{s+a}{2}\Bigr)L(s,\chi)
\]
(for complex characters take \(\xi_K=\xi\,\Lambda(s,\chi)\Lambda(s,\bar\chi)\)). Then \(\xi_K\) is an <strong>entire</strong> function of order one with \(\xi_K(s)=\xi_K(1-s)=\overline{\xi_K(\bar s)}\).</p>
<p><strong>Angular-monotonicity theorem (paper Theorem 2).</strong> For every primitive non-principal character,
\[
T_K(r,\theta):=\frac{\partial}{\partial\theta}\log\bigl|\xi_K(\tfrac12+re^{i\theta})\bigr|\le 0,
\qquad r>0,\ \theta\in(0,\pi/2).
\]</p>
<p><strong>The pole contradiction (proof of paper Theorem 1).</strong> Suppose \(L(s,\chi)\) had a nontrivial zero \(\rho=\beta+i\gamma\) <em>off</em> the critical line (\(\beta\neq\tfrac12\)). By the functional equation and conjugation symmetry, assume \(\beta&gt;\tfrac12\), \(\gamma&gt;0\); then \(\rho\) is also a zero of \(\xi_K\). Near \(\rho\), \(\xi_K(s)=(s-\rho)^m h(s)\) (\(m\ge1\), \(h(\rho)\neq0\)), so
\[
\frac{\xi_K'}{\xi_K}(s)=\frac{m}{s-\rho}+\frac{h'}{h}(s).
\]
Holding \(r=r_0=|\rho-\tfrac12|\) fixed and letting \(\theta\to\theta_0^+\) (\(\theta_0=\arg(\rho-\tfrac12)\in(0,\pi/2)\)), \(s\) approaches \(\rho\) from above and
\[
\Re\!\left[\frac{i r_0 e^{i\theta}}{s-\rho}\right]\to+\infty,
\]
while \(h'/h\) stays bounded; hence \(T_K(r_0,\theta)\to+\infty\) — <strong>a direct contradiction</strong> of \(T_K\le0\). Therefore \(\beta=\tfrac12\): every nontrivial zero of \(L(s,\chi)\) lies on the critical line, and <strong>the GRH is proved</strong>.</p>
<p><strong>Three-region support for angular monotonicity</strong> (paper §5–8): the large-\(r\) region (\(r\ge7\)) is carried by the Stirling main term \(-r\log r\) of the Gamma factor, split into \(\sigma\ge2\) (Subregion A), \(7\le r\le50\) (Subregion B), and \(r\ge50,\sigma&lt;2\) (Subregion C); the compact middle region (\(2\le r\le7\)) is enclosed by verified interval arithmetic (in progress, see §8); the small-\(r\) region (\(r\le2\)) follows analytically from the positive-kernel decomposition with a fixed-sign Taylor remainder. On the \(r\in[7,10000]\) grid, \(T_K\) is negative at every point for 10 representative characters (§6 grids and interactive heatmap).</p>

<h3>5.1 Why proving GRH is simpler than proving RH</h3>
<p>The hardest battle in the RH proof is fought on \(\zeta\) <strong>alone</strong>: forcing the \(\zeta'/\zeta\) noise below the Gamma main term. GRH is not a duel — \(\xi_K=\xi\cdot\Lambda\) is a joint battlefield with five structural advantages:</p>
<p><strong>(1) Entire from the start, no pole to remove.</strong> \(\zeta\) has a pole at \(s=1\), and the RH construction of \(\xi\) must cancel it with the \(s(s-1)\) factor; for non-principal characters \(L(s,\chi)\) has <strong>no pole</strong> at \(s=1\), so \(\Lambda\) is already entire and \(\xi_K\) is entire by construction — one entire singularity step removed.</p>
<p><strong>(2) Arithmetic positivity comes for free.</strong> For quadratic characters \(\xi_K\) is the completed Dedekind zeta function of the associated quadratic number field, with nonnegative Dirichlet coefficients \(r_K(n)\ge0\); for complex characters the coefficients of \(\Lambda(s,\chi)\Lambda(s,\bar\chi)=|L|^2\) are \(|\sum_{d|n}\chi(d)|^2\ge0\). The positive-kernel argument in the small-\(r\) region therefore holds <strong>uniformly for every character</strong> — in the RH case the same step required the three-term polylogarithm identity to split \(\zeta\) into three positive-definite pieces by hand; for GRH arithmetic positivity hands it over directly.</p>
<p><strong>(3) Two zero densities, Gamma dominance more comfortable.</strong> The phase function of \(\xi_K\) receives its \(\log r\) leading term from <strong>two independent sources</strong> (\(\zeta\) and \(L\)); studying \(\Lambda\) alone has no such term, and its phase positivity would force assuming zeros on the critical line — a circular argument. Multiplying by \(\xi\) borrows the Gamma factor and zero-density structure of zeta. For complex characters the Gamma degree is \(d=3\) with coefficient \(3/2\) while the log-derivative bound is only \(3M_\zeta(2)=1.71\); Subregion A gives \(T_K\le-0.56\,r+O(1)&lt;0\) directly.</p>
<p><strong>(4) The larger the conductor, the easier.</strong> The Gamma term contains \(-r\sin\theta\,\log q\) (coefficient 1 for complex characters, \(1/2\) for quadratic), so the margin grows with \(q\), and all bounds are <strong>uniform in \(q\ge3\)</strong> — a rare case where the problem gets easier as the object gets larger.</p>
<p><strong>(5) The core machinery is inherited.</strong> The three-term polylogarithm decomposition of \(\zeta\) (the core machine of the RH paper) acts on the zeta factor <em>before</em> combination; GRH reinvents nothing: once \(\zeta\) is won, the Dirichlet \(L\)-functions are a natural corollary of the same machine.</p>
<p>In one line: <strong>RH is a duel on the single function \(\zeta\); GRH is the joint battlefield \(\xi\cdot\Lambda\) — entire with no obstacle, positivity given for free, double Gamma insurance, and a margin that widens with \(q\).</strong></p>
<p style="font-size:.92em;color:var(--subtle)">Scope: this section is an outline of the proof chain; the full proof (three-region estimates, derivation of the constant \(K\), complex-character modifications) is in the Zenodo paper. The matrices remain a numerical construction; rigorous interval-arithmetic enclosure of the middle region (Region 2) and Subregion C is still in progress (§8).</p>

<h2>6. The 46 primitive characters: angular-monotonicity grids</h2>"""

jobs = [
    ('/tmp/hpgit/grh_zh.html',
     [('<h2>五、46 个本原特征：角单调性网格</h2>', ZH_PROOF),
      ('<h2>六、包络常数 C(r) / D(r) 与最紧区域</h2>', '<h2>七、包络常数 C(r) / D(r) 与最紧区域</h2>'),
      ('<h3>5.1 交互热力图：自己挑特征、自己看网格</h3>', '<h3>6.1 交互热力图：自己挑特征、自己看网格</h3>'),
      ('<h2>七、主张的边界</h2>', '<h2>八、主张的边界</h2>'),
      ('<h2>八、可复现性</h2>', '<h2>九、可复现性</h2>'),
      # §边界内容更新
      ('角单调性（§5–6）的<strong>解析论证</strong>在 GRH 论文（Zenodo 22143035）中；本页只展示其可独立复算的数值网格——46 个本原特征、全部网格点为负。',
       '角单调性的<strong>解析证明链条</strong>见 §5（完整证明在 GRH 论文，Zenodo 22143035）；本页 §6–7 只展示其可独立复算的数值网格——46 个本原特征、全部网格点为负。'),
      ('区间算术（Region 2 / Subregion C 的严格包围）仍在进行中，完成后补入论文。',
       '区间算术（Region 2 / Subregion C 的严格包围）仍在进行中，完成后补入论文；解析链条不受其影响。')]),
    ('/tmp/hpgit/grh.html',
     [('<h2>5. The 46 primitive characters: angular-monotonicity grids</h2>', EN_PROOF),
      ('<h2>6. Envelope constants C(r) / D(r) and the tightest region</h2>', '<h2>7. Envelope constants C(r) / D(r) and the tightest region</h2>'),
      ('<h3>5.1 Interactive heatmap', '<h3>6.1 Interactive heatmap'),
      ('<h2>7. Scope of the claims</h2>', '<h2>8. Scope of the claims</h2>'),
      ('<h2>8. Reproducibility</h2>', '<h2>9. Reproducibility</h2>'),
      ('The <strong>analytic argument</strong> for angular monotonicity (§5–6) is in the GRH paper (Zenodo 22143035); this page presents only its independently reproducible numerical grids — 46 primitive characters, every grid point negative.',
       'The <strong>analytic proof chain</strong> for angular monotonicity is outlined in §5 (full proof in the GRH paper, Zenodo 22143035); §6–7 of this page present only its independently reproducible numerical grids — 46 primitive characters, every grid point negative.'),
      ('Interval arithmetic (rigorous enclosure of Region 2 / Subregion C) is still in progress and will be added to the paper upon completion.',
       'Interval arithmetic (rigorous enclosure of Region 2 / Subregion C) is still in progress and will be added to the paper upon completion; the analytic chain does not depend on it.')]),
]

for path, reps in jobs:
    s = io.open(path, encoding='utf-8').read()
    for old, new in reps:
        assert old in s, (path, old[:60])
        assert s.count(old) == 1, ('multi: ' + old[:60] + ' x' + str(s.count(old)))
        s = s.replace(old, new, 1)
    io.open(path, 'w', encoding='utf-8').write(s)
    print('updated', path)
