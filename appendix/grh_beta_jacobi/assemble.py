# -*- coding: utf-8 -*-
"""Assemble GRH interactive page(s): inject SVGs, tables, heatmap, lightbox JS."""
import json, re, html

R = json.load(open("/tmp/grh_beta/beta_J50_results.json"))
W = json.load(open("/tmp/grh_beta/grh_web_data.json"))
svg_struct = open("/tmp/grh_web/svg_structure.html").read()
svg_spec   = open("/tmp/grh_web/svg_spectrum.html").read()

# extract lightbox JS from main zh page
zh = open("/tmp/hpgit/zh.html").read()
m = re.search(r"(<script>\s*\(function\(\)\{.*?</script>)", zh, re.S)
LB_ZH = m.group(1)
# english hint version
lb_en = LB_ZH.replace("滚轮/双指缩放 · 拖拽平移 · 双击放大 · 点 × 或按 Esc 关闭 · 点彩色格子读完整数值",
                      "Scroll / pinch to zoom · drag to pan · double-click to zoom in · click × or press Esc to close · click a colored cell for full precision")

comp = R["comparison"]
zeros = R["beta_zeros_independent"]
inv = R["inv_sqrt"]
NLOCKED = sum(1 for r_ in comp if r_["rel_err"] < 1e-4)

# ---------- lock table ----------
def lock_table():
    rows = ["<tr><th>n</th><th>矩阵谱 1/√λₙ</th><th>β 零点 γₙ（独立扫描）</th><th>相对误差</th><th>状态</th></tr>"]
    for r_ in comp[:12]:
        ok = r_["rel_err"] < 1e-4
        rows.append(f'<tr><td>{r_["k"]}</td><td>{r_["inv_sqrt"]:.10f}</td><td>{r_["gamma"]:.10f}</td>'
                    f'<td>{r_["rel_err"]:.2e}</td>'
                    f'<td class="{"pass" if ok else ""}">{"✓ 锁定" if ok else "收敛中"}</td></tr>')
    rows.append(f'<tr><td>…</td><td colspan="4" style="color:#666">共 {NLOCKED} 个零点锁定（相对误差 &lt; 10⁻⁴），其余随阶数收敛；ζ 主篇 J₅₀ 锁定 25 个</td></tr>')
    return "<table>" + "".join(rows) + "</table>"

# ---------- grid table (10 chars) ----------
NAME_ZH = {"chi3":"χ₃ (q=3)","chi4":"χ₄=β (q=4)","chi5e":"χ₅ 偶 (q=5)","chi5c":"χ₅ 复 (q=5)",
           "chi7":"χ₇ (q=7)","chi8e":"χ₈ 偶 (q=8)","chi8o":"χ₈ 奇 (q=8)","chi11":"χ₁₁ (q=11)",
           "chi13":"χ₁₃ (q=13)","chi101":"χ₁₀₁ (q=101)"}
NAME_EN = {"chi3":"χ₃ (q=3)","chi4":"χ₄=β (q=4)","chi5e":"χ₅ even (q=5)","chi5c":"χ₅ complex (q=5)",
           "chi7":"χ₇ (q=7)","chi8e":"χ₈ even (q=8)","chi8o":"χ₈ odd (q=8)","chi11":"χ₁₁ (q=11)",
           "chi13":"χ₁₃ (q=13)","chi101":"χ₁₀₁ (q=101)"}
order = ["chi3","chi4","chi5e","chi5c","chi7","chi8e","chi8o","chi11","chi13","chi101"]

def grid_table(names):
    rows = ["<tr><th>特征</th><th>网格点数</th><th>max T(r,θ)</th><th>最危险 r</th><th>最危险 θ</th><th>结果</th></tr>"]
    for k in order:
        g = W["grid"][k]
        rows.append(f'<tr><td>{names[k]}</td><td>{g["n_pts"]}</td><td>{g["max_T"]:.4e}</td>'
                    f'<td>{g["worst_r"]:g}</td><td>{g["worst_th"]:.4f}</td>'
                    f'<td class="pass">{g["status"]}</td></tr>')
    return "<table>" + "".join(rows) + "</table>"

# ---------- C/D table ----------
def cd_table(names):
    rows = ["<tr><th>特征</th><th>C(0.5)</th><th>D(0.5)</th><th>C(2.0)</th><th>D(2.0)</th></tr>"]
    for k in order:
        if k not in W["CD"]: continue
        c = W["CD"][k]
        rows.append(f'<tr><td>{names[k]}</td><td>{c["C05"]:.4f}</td><td>{c["D05"]:.4f}</td>'
                    f'<td>{c["C20"]:.4f}</td><td>{c["D20"]:.4f}</td></tr>')
    return "<table>" + "".join(rows) + "</table>"

# ---------- heatmap ----------
def heatmap_js():
    data = {}
    for ch in W["chars"]:
        data[ch["name"]] = {"detail": ch["detail"], "TK_max": ch["TK_max"], "q": ch["q"],
                            "complex": ch["complex"], "n": ch["n_points"]}
    payload = json.dumps({"theta": W["theta"], "r": W["r_full"], "chars": data, "order": order},
                         ensure_ascii=False, separators=(",", ":"))
    return """
<div id="heatmap-wrap"></div>
<script>
(function(){
  var DATA = __PAYLOAD__;
  var order = DATA.order, theta = DATA.theta, rs = DATA.r;
  var NR = rs.length, NT = theta.length;
  // global min/max for color scale (all negative)
  var gmin = 1e9, gmax = -1e9;
  order.forEach(function(nm){
    DATA.chars[nm].detail.forEach(function(p){
      if(p[2]<gmin)gmin=p[2]; if(p[2]>gmax)gmax=p[2];
    });
  });
  function col(v){
    // all v<0: near 0 -> light blue-green; very negative -> deep teal/blue
    var t = (v - gmin)/(gmax - gmin); // 0=most negative (deep), 1=near zero (light)
    t = Math.max(0, Math.min(1, t));
    // deep blue #0d3b66 -> teal #1b998b -> light #cfe8e0
    function mix(a,b,f){return Math.round(a+(b-a)*f);}
    var c1=[13,59,102], c2=[27,153,139], c3=[222,240,235];
    var r,g,bl;
    if(t<0.55){var f=t/0.55;r=mix(c1[0],c2[0],f);g=mix(c1[1],c2[1],f);bl=mix(c1[2],c2[2],f);}
    else{var f=(t-0.55)/0.45;r=mix(c2[0],c3[0],f);g=mix(c2[1],c3[1],f);bl=mix(c2[2],c3[2],f);}
    return 'rgb('+r+','+g+','+bl+')';
  }
  function build(nm){
    var ch = DATA.chars[nm];
    var det = ch.detail;
    // map [r,theta,T] -> cell
    var cell = {};
    det.forEach(function(p){ cell[p[0].toFixed(3)+'_'+p[1].toFixed(3)] = p[2]; });
    var CW = 720, CH = 300, L = 60, Rm = CW-20, Tp = 14, Bm = CH-34;
    var w = (Rm-L)/NT, h = (Bm-Tp)/NR;
    var s = '<svg viewBox="0 0 '+CW+' '+CH+'" style="width:100%;max-width:760px;height:auto;background:#fff;border:1px solid #e0e0e0;border-radius:6px;display:block;margin:.4em auto">';
    for(var i=0;i<NR;i++){
      for(var j=0;j<NT;j++){
        var key = rs[i].toFixed(3)+'_'+theta[j].toFixed(3);
        var v = cell[key];
        if(v===undefined){ s+='<rect x="'+(L+j*w)+'" y="'+(Tp+i*h)+'" width="'+(w-0.6)+'" height="'+(h-0.6)+'" fill="#f2f2f2"/>'; continue; }
        var ttl = 'r='+rs[i]+', θ='+theta[j]+', T='+v.toExponential(4);
        s+='<rect x="'+(L+j*w)+'" y="'+(Tp+i*h)+'" width="'+(w-0.6)+'" height="'+(h-0.6)+'" fill="'+col(v)+'"><title>'+ttl+'</title></rect>';
      }
    }
    // axis labels
    s+='<text x="'+((L+Rm)/2)+'" y="'+(CH-8)+'" text-anchor="middle" font-size="12" fill="#555">θ →</text>';
    s+='<text x="16" y="'+((Tp+Bm)/2)+'" text-anchor="middle" font-size="12" fill="#555" transform="rotate(-90 16 '+((Tp+Bm)/2)+'))">r →</text>';
    for(var j2=0;j2<NT;j2+=2){ s+='<text x="'+(L+j2*w+w/2)+'" y="'+(Bm+14)+'" text-anchor="middle" font-size="9" fill="#666">'+theta[j2]+'</text>'; }
    for(var i2=0;i2<NR;i2+=2){ s+='<text x="'+(L-6)+'" y="'+(Tp+i2*h+h/2+3)+'" text-anchor="end" font-size="9" fill="#666">'+rs[i2]+'</text>'; }
    s+='</svg>';
    document.getElementById('heatmap-wrap').innerHTML = s;
    var cap = document.getElementById('heat-cap');
    if(cap) cap.textContent = nm+'：'+det.length+' 个网格点，max T = '+ch.TK_max.toExponential(4)+'（全部为负）。深色=更负，浅色=接近 0；点格子读数值。';
  }
  var btns = document.querySelectorAll('.ctrl button[data-char]');
  btns.forEach(function(b){
    b.addEventListener('click', function(){
      btns.forEach(function(x){x.classList.remove('on');});
      b.classList.add('on');
      build(b.getAttribute('data-char'));
    });
  });
  build(order[0]);
})();
</script>
""".replace("__PAYLOAD__", payload)

def char_buttons(names):
    out = []
    for i, k in enumerate(order):
        out.append(f'<button data-char="{k}"{" class=\"on\"" if i==0 else ""}>{names[k]}</button>')
    return "\n".join(out)

TIGHT = f'{W["tight_q3"]["max_T"]:.4f}（{W["tight_q3"]["n_pts"]} 个点，r={W["tight_q3"]["worst_r"]:g}）'

def build(tpl_path, out_path, names, lb, lang):
    t = open(tpl_path).read()
    t = t.replace("__SVG_STRUCTURE__", svg_struct)
    t = t.replace("__SVG_SPECTRUM__", svg_spec)
    t = t.replace("__NLOCKED__", str(NLOCKED))
    t = t.replace("__TABLE_LOCK__", lock_table())
    t = t.replace("__TABLE_GRID__", grid_table(names))
    t = t.replace("__TABLE_CD__", cd_table(names))
    t = t.replace("__TIGHT_MAX__", TIGHT)
    t = t.replace("__CHAR_BUTTONS__", char_buttons(names))
    t = t.replace("__HEATMAP__", heatmap_js())
    t = t.replace("__SCRIPT_END__", lb)
    open(out_path, "w").write(t)
    print("wrote", out_path, len(t), "bytes; nlocked =", NLOCKED)

build("/tmp/grh_web/page_template_zh.html", "/tmp/grh_web/grh_zh.html", NAME_ZH, LB_ZH, "zh")

build("/tmp/grh_web/page_template_en.html", "/tmp/grh_web/grh.html", NAME_EN, lb_en, "en")
