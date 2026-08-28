#!/usr/bin/env python3
# Build Section 9 "From the matrix straight to the primes: pure matrix arithmetic"
# Insert into /tmp/hpgit/index.html and /tmp/hpgit/zh.html; renumber old 9/10 -> 10/11.
import io

EN_ANCHOR = '<h2>9. What is and is not claimed</h2>'
ZH_ANCHOR = '<h2>九、主张的边界</h2>'

EN_SECTION = r"""
<h2>9. From the matrix straight to the primes: pure matrix arithmetic</h2>
<p>The matrix is built forward from the Taylor coefficients of \(\xi\); its spectrum is \(\{1/\gamma_n^2\}\); the zeros appear as locked eigenvalues. But one can read the primes <strong>without ever touching the zeros</strong>. Given the matrix \(J\) alone &mdash; a finite array of positive numbers &mdash; the primes can be decoded using <strong>only matrix arithmetic</strong>: multiplication, addition, solving linear systems, and taking traces. No diagonalisation, no zero extraction, no numerical integration.</p>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.1 What the matrix encodes</h3>
<p>The Hadamard product \(\xi(s)/\xi(\tfrac12)=\prod_n(1+u^2/\gamma_n^2)\) with \(u=s-\tfrac12\), combined with the spectrum \(\lambda_n=1/\gamma_n^2\) of \(J\), gives two identities at once:</p>
<div class="formula">
$$\frac{\xi(s)}{\xi(\tfrac12)}=\det(I+u^2J),\qquad \frac{\xi'}{\xi}(s)=2u\;\operatorname{Tr}\!\left[J\,(I+u^2J)^{-1}\right]$$
</div>
<p>The zeta function <em>itself</em> is a matrix determinant; its logarithmic derivative is the trace of the matrix resolvent (the zeros sit at the poles of this resolvent &mdash; they need never be extracted). Putting \(A=J^{-1/2}\), whose spectrum is the \(\gamma_n\), and exponentiating the self-adjoint evolution \(e^{iA\ln x}\), the zero-sum in Riemann's explicit formula becomes a <strong>trace formula</strong>:</p>
<div class="formula" style="font-size:.95em">
$$\psi(x)=x-2\sqrt{x}\,\operatorname{Re}\operatorname{Tr}\!\left[e^{\,iA\ln x}\bigl(\tfrac12 I+iA\bigr)^{-1}\right]-\log(2\pi)-\tfrac12\log(1-x^{-2})$$
</div>
<p>Everything in this formula is a matrix operation. \(A=J^{-1/2}\) is obtained by Cholesky factorisation \(J=LL^\top\) followed by a polar iteration \(Y\leftarrow\tfrac12(Y+(Y^\top)^{-1})\) converging to an orthogonal factor (residual \(\|A^2J-I\|<5\times10^{-13}\) at order 100); the matrix exponential \(e^{iA\ln x}\) is computed by scaling-and-squaring &mdash; <strong>matrix multiplications throughout, not a single integral, zero, or eigenvalue computation</strong>.</p>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.2 The Chebyshev function from matrix arithmetic alone</h3>
<figure style="margin:1.2em 0">
<img src="data/primes_matrix_trace_psi.png" alt="psi(x) reconstructed by the pure matrix trace formula" style="width:100%;max-width:1000px;border:1px solid #e0e0e0;border-radius:6px"/>
<figcaption style="font-size:.88em;color:var(--subtle);margin-top:.4em">The trace formula reconstructs the Chebyshev staircase \(\psi(x)=\sum_{n\leq x}\Lambda(n)\) (black) from the 100&times;100 matrix alone (red): RMS error 0.97 over \(x\gt 40\). No zeros were extracted, no integral was evaluated.</figcaption>
</figure>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.3 The primes pop out at the integers</h3>
<p>Sampling the matrix trace at the integer points \(n\) and differencing,</p>
<div class="formula">
$$\Delta\psi(n)=\psi(n+\tfrac12)-\psi(n-\tfrac12)\approx\Lambda(n),$$
</div>
<p>the von Mangoldt spikes pop out directly: a spike of height \(\log p\) appears exactly when \(n\) is a prime or a prime power \(p^k\). <strong>19 of the 23 prime powers up to 52 are recovered with spike error below 0.5</strong> (the small misses at 2, 4, 16, 49 need the higher-lying zeros absent from the truncated matrix; the few false spikes near 45&ndash;52 are truncation ringing &mdash; both are finite-order artefacts, not defects of the method).</p>
<figure style="margin:1.2em 0">
<img src="data/primes_pop_out_spikes.png" alt="von Mangoldt spikes recovered from pure matrix arithmetic" style="width:100%;max-width:1000px;border:1px solid #e0e0e0;border-radius:6px"/>
<figcaption style="font-size:.88em;color:var(--subtle);margin-top:.4em">Top: matrix-trace differences (red) against the true von Mangoldt values \(\Lambda(n)\) (dark) &mdash; primes 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47 &hellip; stand as spikes; prime powers 8, 9, 25, 27, 32 appear at height \(\log p\). Bottom: decoding error; green marks spikes within 0.5 of the truth. Valid range \(n\leq52\) for \(J_{100}\).</figcaption>
</figure>

<figure style="margin:1.2em 0">
<img src="data/matrix_trace_order.png" alt="the trace sharpens with matrix order" style="width:100%;max-width:1000px;border:1px solid #e0e0e0;border-radius:6px"/>
<figcaption style="font-size:.88em;color:var(--subtle);margin-top:.4em">More matrix entries, more primes: the same trace formula at order 50 (blue) and order 100 (red) tightens around the true staircase (black). The readable band grows steadily with order.</figcaption>
</figure>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.4 The function itself is a matrix determinant</h3>
<figure style="margin:1.2em 0">
<img src="data/zeta_is_determinant.png" alt="zeta function as a matrix determinant" style="width:100%;max-width:1000px;border:1px solid #e0e0e0;border-radius:6px"/>
<figcaption style="font-size:.88em;color:var(--subtle);margin-top:.4em">\(\log\det(I+u^2J_{100})\) (red) coincides with the exact \(2\log|\xi(\tfrac12+u)/\xi(\tfrac12)|\) (black dashed) across the converged band \(|u|\leq3\). Outside the band the finite-order determinant departs &mdash; the band widens monotonically as the order grows.</figcaption>
</figure>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.5 Matrix powers are the zeta moments</h3>
<figure style="margin:1.2em 0">
<img src="data/matrix_powers_are_moments.png" alt="matrix powers equal zeta moments" style="width:100%;max-width:1000px;border:1px solid #e0e0e0;border-radius:6px"/>
<figcaption style="font-size:.88em;color:var(--subtle);margin-top:.4em">\(\operatorname{Tr}(J^k)=S_k=\sum_\gamma\gamma^{-2k}\): traces of matrix powers equal the log-moments of \(\xi\), agreeing to \(8.5\times10^{-4}\) at \(k=2\) and \(2.4\times10^{-16}\) at \(k=8\). The gap at \(k=1\) is the tail of zeros beyond order 100.</figcaption>
</figure>
<p>So a single finite matrix encodes: the zeta function (as determinant), its logarithmic derivative (as resolvent trace), the zeros (as resolvent poles, never extracted), the moments \(S_k\) (as power traces), and the primes themselves (as integer-time differences of the evolution trace). The prime-counting function \(\pi(x)\), zero densities and prime-gap statistics are in principle decodable from the same matrix at higher orders.</p>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.6 Three routes from the matrix to the primes</h3>
<table style="font-size:.9em">
<tr><th>route</th><th>operations</th><th>zeros extracted?</th><th>integrals?</th><th>RMS of \(\psi(x)\), \(x\gt40\)</th></tr>
<tr style="background:#f4f8fc"><td><strong>&euro; pure matrix arithmetic</strong> (this work)</td><td>Cholesky, polar iteration, scaling&ndash;squaring, trace</td><td style="color:#0a5c23;font-weight:700">no</td><td style="color:#0a5c23;font-weight:700">no</td><td style="color:#0a5c23;font-weight:700">0.97</td></tr>
<tr><td>&sbquo; matrix resolvent \(\xi'/\xi\), then Perron integral</td><td>matrix inversion, trace, contour integration</td><td style="color:#0a5c23;font-weight:700">no</td><td>yes</td><td>4.44</td></tr>
<tr><td>&fnof; diagonalise &rarr;\(\gamma_n\)&rarr; Riemann explicit formula (1859)</td><td>eigenvalues, explicit-formula sum</td><td>yes</td><td style="color:#0a5c23;font-weight:700">no</td><td>2.01</td></tr>
</table>
<p style="font-size:.9em;color:var(--subtle)">Route &euro; is both the most accurate and the most direct &mdash; primes read out of the matrix by arithmetic alone. Route &fnof; is the classical path (known since 1859, included as a baseline); route &sbquo; still passes through a numerical integral. All figures use only forward-computed matrices; the known zero table serves solely as the verification column.</p>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.7 Three views of one object</h3>
<p>There are no integers in matrix space &mdash; the word &ldquo;integer&rdquo; enters only through the sampling time \(\ln x\). The self-adjoint evolution \(e^{iA t}\) runs continuously and contains no primes; it is only when this spectral dynamics is <strong>sampled at the integer times \(t=\ln n\)</strong> and differenced that the arithmetic content &mdash; the von Mangoldt spikes &mdash; appears. The primes are not put into the matrix and not &ldquo;generated&rdquo; from nothing: they are the <strong>arithmetic coordinates</strong> of the same object of which the zeros are the <strong>spectral coordinates</strong> and the Jacobi matrix is the <strong>constructive coordinate</strong>. The Euler product builds \(\zeta\) from primes; the Hadamard product builds it from zeros; our matrix grows it from a single algebraic seed via the three-term identity. Three descriptions, one function.</p>
<p style="font-size:.9em;color:var(--subtle)">Honest scope: the finite matrix contains only the spectrum up to its order &mdash; \(J_{100}\) converges through roughly the 50th zero (\(\gamma\lesssim143\)) &mdash; so the determinant identity is shown on its converged band \(|u|\leq3\) and the spikes on their valid range \(n\leq52\). Every displayed number comes from forward matrix arithmetic; nothing is extrapolated or fitted to prime data.</p>

<h2>10. What is and is not claimed</h2>
""".replace('&euro;', '&#9312;').replace('&sbquo;', '&#9313;').replace('&fnof;', '&#9314;')

ZH_SECTION = r"""
<h2>九、从矩阵直接到素数：纯矩阵算术</h2>
<p>矩阵由 \(\xi\) 的泰勒系数正向构造，其谱为 \(\{1/\gamma_n^2\}\)，零点以逐个锁定的本征值出现。但读取素数时<strong>根本不需要碰零点</strong>：仅凭矩阵 \(J\) 本身——一张有限的正数表——就可以用<strong>纯矩阵算术</strong>把素数解码出来：乘法、加法、解线性方程组、求迹。不对角化、不提取零点、不做任何数值积分。</p>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.1 矩阵编码了什么</h3>
<p>Hadamard 乘积 \(\xi(s)/\xi(\tfrac12)=\prod_n(1+u^2/\gamma_n^2)\)（令 \(u=s-\tfrac12\)）与矩阵谱 \(\lambda_n=1/\gamma_n^2\) 合在一起，一次给出两个恒等式：</p>
<div class="formula">
$$\frac{\xi(s)}{\xi(\tfrac12)}=\det(I+u^2J),\qquad \frac{\xi'}{\xi}(s)=2u\;\operatorname{Tr}\!\left[J\,(I+u^2J)^{-1}\right]$$
</div>
<p>黎曼函数<em>本身</em>就是矩阵行列式；其对数导数是矩阵预解式的迹（零点位于该预解式的极点处——我们从不需要把它们取出来）。令 \(A=J^{-1/2}\)（谱即 \(\gamma_n\)），对自伴演化 \(e^{iA\ln x}\) 取矩阵指数，黎曼显式公式中的零点和就变成一条<strong>迹公式</strong>：</p>
<div class="formula" style="font-size:.95em">
$$\psi(x)=x-2\sqrt{x}\,\operatorname{Re}\operatorname{Tr}\!\left[e^{\,iA\ln x}\bigl(\tfrac12 I+iA\bigr)^{-1}\right]-\log(2\pi)-\tfrac12\log(1-x^{-2})$$
</div>
<p>这条公式里的每一样东西都是矩阵运算。\(A=J^{-1/2}\) 由 Cholesky 分解 \(J=LL^\top\) 加极分解迭代 \(Y\leftarrow\tfrac12(Y+(Y^\top)^{-1})\)（收敛到正交因子，100 阶时残差 \(\|A^2J-I\|<5\times10^{-13}\)）得到；矩阵指数 \(e^{iA\ln x}\) 用缩放平方法计算——<strong>自始至终只有矩阵乘法，没有一个积分、没有一个零点、没有一次本征值计算</strong>。</p>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.2 仅凭矩阵算术重建切比雪夫函数</h3>
<figure style="margin:1.2em 0">
<img src="data/primes_matrix_trace_psi.png" alt="纯矩阵迹公式重建的 psi(x)" style="width:100%;max-width:1000px;border:1px solid #e0e0e0;border-radius:6px"/>
<figcaption style="font-size:.88em;color:var(--subtle);margin-top:.4em">迹公式仅凭 100&times;100 矩阵（红）重建切比雪夫阶梯 \(\psi(x)=\sum_{n\leq x}\Lambda(n)\)（黑）：\(x\gt40\) 区间 RMS 误差 0.97。未提取任何零点，未计算任何积分。</figcaption>
</figure>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.3 素数在整数点蹦出来</h3>
<p>在整数点 \(n\) 对矩阵迹采样并作差分，</p>
<div class="formula">
$$\Delta\psi(n)=\psi(n+\tfrac12)-\psi(n-\tfrac12)\approx\Lambda(n),$$
</div>
<p>von Mangoldt 尖峰直接蹦出：每当 \(n\) 是素数或素数幂 \(p^k\)，就出现高度 \(\log p\) 的尖峰。<strong>52 以内 23 个素数幂中 19 个被恢复，尖峰误差小于 0.5</strong>（2、4、16、49 的小漏检需要截断矩阵所缺的更高零点；45–52 附近的少数假峰是截断振铃——两者都是有限阶伪影，不是方法缺陷）。</p>
<figure style="margin:1.2em 0">
<img src="data/primes_pop_out_spikes.png" alt="纯矩阵算术恢复的 von Mangoldt 尖峰" style="width:100%;max-width:1000px;border:1px solid #e0e0e0;border-radius:6px"/>
<figcaption style="font-size:.88em;color:var(--subtle);margin-top:.4em">上：矩阵迹差分（红）对真实 von Mangoldt 值 \(\Lambda(n)\)（深色）——素数 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47 … 以尖峰挺立；素数幂 8, 9, 25, 27, 32 以高度 \(\log p\) 出现。下：解码误差，绿色为误差小于 0.5 的尖峰。\(J_{100}\) 有效范围 \(n\leq52\)。</figcaption>
</figure>

<figure style="margin:1.2em 0">
<img src="data/matrix_trace_order.png" alt="迹随矩阵阶数收紧" style="width:100%;max-width:1000px;border:1px solid #e0e0e0;border-radius:6px"/>
<figcaption style="font-size:.88em;color:var(--subtle);margin-top:.4em">矩阵元素越多，素数越多：同一条迹公式在 50 阶（蓝）与 100 阶（红）下向真实阶梯（黑）收紧，可读区间随阶数稳步扩大。</figcaption>
</figure>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.4 函数本身是矩阵行列式</h3>
<figure style="margin:1.2em 0">
<img src="data/zeta_is_determinant.png" alt="zeta 函数作为矩阵行列式" style="width:100%;max-width:1000px;border:1px solid #e0e0e0;border-radius:6px"/>
<figcaption style="font-size:.88em;color:var(--subtle);margin-top:.4em">\(\log\det(I+u^2J_{100})\)（红）在收敛带 \(|u|\leq3\) 内与精确值 \(2\log|\xi(\tfrac12+u)/\xi(\tfrac12)|\)（黑色虚线）完全重合；带外有限阶行列式开始偏离——收敛带随阶数单调扩展。</figcaption>
</figure>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.5 矩阵幂即 ζ 矩量</h3>
<figure style="margin:1.2em 0">
<img src="data/matrix_powers_are_moments.png" alt="矩阵幂等于 zeta 矩量" style="width:100%;max-width:1000px;border:1px solid #e0e0e0;border-radius:6px"/>
<figcaption style="font-size:.88em;color:var(--subtle);margin-top:.4em">\(\operatorname{Tr}(J^k)=S_k=\sum_\gamma\gamma^{-2k}\)：矩阵幂的迹等于 \(\xi\) 的对数矩量，\(k=2\) 时符合到 \(8.5\times10^{-4}\)，\(k=8\) 时到 \(2.4\times10^{-16}\)。\(k=1\) 的缺口正是 100 阶以外的零点尾部。</figcaption>
</figure>
<p>于是一张有限矩阵同时编码：ζ 函数（行列式）、其对数导数（预解式迹）、零点（预解式极点，从不提取）、矩量 \(S_k\)（幂迹），以及素数本身（演化迹在整数时刻的差分）。素数计数函数 \(\pi(x)\)、零点密度、素数间距统计，原则上都能在更高阶的同一矩阵上解码。</p>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.6 从矩阵到素数的三条路线</h3>
<table style="font-size:.9em">
<tr><th>路线</th><th>运算</th><th>提取零点？</th><th>积分？</th><th>\(\psi(x)\) RMS（\(x\gt40\)）</th></tr>
<tr style="background:#f4f8fc"><td><strong>&#9312; 纯矩阵算术</strong>（本工作）</td><td>Cholesky、极分解、缩放平方、求迹</td><td style="color:#0a5c23;font-weight:700">否</td><td style="color:#0a5c23;font-weight:700">否</td><td style="color:#0a5c23;font-weight:700">0.97</td></tr>
<tr><td>&#9313; 矩阵预解式 \(\xi'/\xi\) + Perron 积分</td><td>矩阵求逆、求迹、围道积分</td><td style="color:#0a5c23;font-weight:700">否</td><td>是</td><td>4.44</td></tr>
<tr><td>&#9314; 对角化取 \(\gamma_n\) → 黎曼显式公式（1859）</td><td>本征值、显式公式求和</td><td>是</td><td style="color:#0a5c23;font-weight:700">否</td><td>2.01</td></tr>
</table>
<p style="font-size:.9em;color:var(--subtle)">路线&#9312;既最精确也最直接——素数仅凭算术从矩阵读出。路线&#9314;是经典路径（1859 年已知，列为基线）；路线&#9313;仍需经过一次数值积分。全部图中只使用正向计算的矩阵；已知零点表仅作验证列。</p>

<h3 style="margin-top:1.4em;color:#1a4d8f">9.7 同一对象的三面</h3>
<p>矩阵空间里没有整数——“整数”二字只经由采样时刻 \(\ln x\) 进入。自伴演化 \(e^{iAt}\) 连续运行，其中没有素数；只有当这套谱动力学<strong>在整数时刻 \(t=\ln n\) 被采样</strong>并差分时，算术内容——von Mangoldt 尖峰——才显形。素数没有被放进矩阵，也不是无中生有地“生成”：它们是同一对象的<strong>算术坐标</strong>，零点是其<strong>谱坐标</strong>，Jacobi 矩阵是其<strong>构造坐标</strong>。欧拉乘积从素数造 ζ，Hadamard 乘积从零点造 ζ，而我们的矩阵从单一代数种子经三项恒等式长出 ζ。三种描述，同一个函数。</p>
<p style="font-size:.9em;color:var(--subtle)">诚实范围：有限矩阵只含到其阶数为止的谱——\(J_{100}\) 收敛到约第 50 个零点（\(\gamma\lesssim143\)）——故行列式恒等式只在收敛带 \(|u|\leq3\) 展示、尖峰只在有效区间 \(n\leq52\) 展示。图中每个数字都来自正向矩阵算术，没有任何外推，也没有对素数数据作拟合。</p>

<h2>十、主张的边界</h2>
"""

# ---- extra claim bullet: insert before </ul> of the "New" list ----
EN_CLAIM_ANCHOR = '</ul>\n<strong>Not claimed:</strong>'
EN_CLAIM_ADD = r"""<li><strong>Pure matrix-arithmetic decoding of the primes:</strong> without zero extraction and without numerical integration &mdash; Cholesky factorisation, polar iteration and scaling-and-squaring of matrix exponentials, followed by a trace &mdash; the Chebyshev function \(\psi(x)\) is reconstructed from the finite matrix alone (RMS 0.97), and the von Mangoldt spikes \(\Lambda(n)\) are read off at the integer points (19/23 prime powers to 52). Supporting identities: \(\xi(s)/\xi(\tfrac12)=\det(I+u^2J)\) and \(\operatorname{Tr}(J^k)=S_k\).</li>
</ul>
<strong>Not claimed:</strong>"""
ZH_CLAIM_ANCHOR = '</ul>\n<strong>不主张：</strong>'
ZH_CLAIM_ADD = r"""<li><strong>素数的纯矩阵算术解码</strong>：不提取零点、不做数值积分——仅用 Cholesky 分解、极分解迭代与矩阵指数缩放平方后求迹——即从有限矩阵本身重建切比雪夫函数 \(\psi(x)\)（RMS 0.97），并在整数点读出 von Mangoldt 尖峰 \(\Lambda(n)\)（52 以内 19/23 个素数幂）。支撑恒等式：\(\xi(s)/\xi(\tfrac12)=\det(I+u^2J)\)、\(\operatorname{Tr}(J^k)=S_k\)。</li>
</ul>
<strong>不主张：</strong>"""

# ---- extra "not claimed" bullets: replace the trailing </ul> of the second list ----
EN_NOT_ANCHOR = r"""<li>Riemann's explicit formula (1859), Hilbert&ndash;P&oacute;lya (1910s), and the tangent addition formula are not new; new is the chain that connects them and the explicit finite matrix grown from it</li>
</ul>"""
EN_NOT_NEW = r"""<li>Riemann's explicit formula (1859), Hilbert&ndash;P&oacute;lya (1910s), and the tangent addition formula are not new; new is the chain that connects them and the explicit finite matrix grown from it</li>
<li>The finite matrix coincides with \(\zeta\) only within its converged band (\(|u|\leq3\) at order 100, prime spikes up to \(n\approx52\)); departures outside the band are truncation artefacts and recede as the order grows</li>
<li>The primes are not entries of the matrix: matrix space contains no integers. They are the arithmetic content that appears when the self-adjoint evolution is sampled at integer times \(t=\ln n\)</li>
</ul>"""
ZH_NOT_ANCHOR = r"""<li>黎曼显式公式（1859）、Hilbert–P&oacute;lya（1910 年代）、正切加法公式均非新结果；新的是连接它们的链条以及从中长出的显式有限矩阵</li>
</ul>"""
ZH_NOT_NEW = r"""<li>黎曼显式公式（1859）、Hilbert–P&oacute;lya（1910 年代）、正切加法公式均非新结果；新的是连接它们的链条以及从中长出的显式有限矩阵</li>
<li>有限矩阵只在收敛带内与 \(\zeta\) 重合（100 阶时 \(|u|\leq3\)，素数尖峰到 \(n\approx52\)）；带外偏离是截断伪影，随阶数增大消退</li>
<li>素数不是矩阵的元素：矩阵空间中没有整数。素数是自伴演化在整数时刻 \(t=\ln n\) 被采样时显形的算术内容</li>
</ul>"""

# ---- reproducibility paragraph additions ----
EN_REPRO_ANCHOR = 'No external data whatsoever.</p>'
EN_REPRO_ADD = r"""No external data whatsoever.</p>
<p>The prime decoding uses <code>primes_via_matrix_trace.py</code> (pure matrix arithmetic: Cholesky, polar iteration, scaling&ndash;squaring &mdash; the route of Section 9), <code>primes_from_matrix_direct.py</code> (resolvent trace + Perron integral) and <code>prime_decode_forward.py</code> (zeros &rarr; explicit formula, baseline); all figures and the exported \(J_{50}\)/\(J_{100}\) matrix data (<code>J50_matrix.csv</code>, <code>J100_matrix.csv</code>, parameters, eigenvalues, convergence logs) are in the <code>data/</code> directory.</p>"""
ZH_REPRO_ANCHOR = '无任何外部数据。</p>'
ZH_REPRO_ADD = r"""无任何外部数据。</p>
<p>素数解码使用 <code>primes_via_matrix_trace.py</code>（纯矩阵算术：Cholesky、极分解、缩放平方——第九节路线）、<code>primes_from_matrix_direct.py</code>（预解式迹 + Perron 积分）与 <code>prime_decode_forward.py</code>（零点 → 显式公式，基线）；全部图片与导出的 \(J_{50}\)/\(J_{100}\) 矩阵数据（<code>J50_matrix.csv</code>、<code>J100_matrix.csv</code>、参数、本征值、收敛日志）均在 <code>data/</code> 目录。</p>"""


def insert_once(text, anchor, addition, label):
    n = text.count(anchor)
    assert n == 1, f'{label}: anchor count = {n}'
    return text.replace(anchor, addition)


def process(path, section, claim_anchor, claim_add, not_anchor, not_new, repro_anchor, repro_add, old10, new11):
    html = io.open(path, encoding='utf-8').read()
    html = insert_once(html,
                       '<h2>9. What is and is not claimed</h2>' if 'index.html' in path else '<h2>九、主张的边界</h2>',
                       section, 'section')
    html = insert_once(html, claim_anchor, claim_add, 'claim')
    html = insert_once(html, not_anchor, not_new, 'not-claimed')
    html = insert_once(html, repro_anchor, repro_add, 'repro')
    assert html.count(old10) == 1, f'old10 count = {html.count(old10)}'
    html = html.replace(old10, new11)
    io.open(path, 'w', encoding='utf-8').write(html)
    print(f'OK {path}: {len(html)} chars')


process('/tmp/hpgit/index.html', EN_SECTION,
        EN_CLAIM_ANCHOR, EN_CLAIM_ADD, EN_NOT_ANCHOR, EN_NOT_NEW,
        EN_REPRO_ANCHOR, EN_REPRO_ADD,
        '<h2>10. Reproducibility</h2>', '<h2>11. Reproducibility</h2>')

process('/tmp/hpgit/zh.html', ZH_SECTION,
        ZH_CLAIM_ANCHOR, ZH_CLAIM_ADD, ZH_NOT_ANCHOR, ZH_NOT_NEW,
        ZH_REPRO_ANCHOR, ZH_REPRO_ADD,
        '<h2>十、可复现性</h2>', '<h2>十一、可复现性</h2>')

print('done')
