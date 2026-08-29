# -*- coding: utf-8 -*-
"""Round 10: 矩阵 SVG 真正可放大 + 放大后格子内可读数值。
1) 每个有色 rect 内插入白色数值 text（值来自该 rect 的 <title>）；
2) 两张 SVG 加 class=zoomable，点击开灯箱（滚轮/双指缩放、拖拽平移、双击）；
3) 灯箱内点彩色格子弹完整精度浮层；文章页悬停 title 保留；
4) 更新两版引导句。"""
import re, io

LIGHTBOX_JS = r'''<script>
(function(){
  var box=document.createElement('div');
  box.style.cssText='display:none;position:fixed;inset:0;z-index:10000;background:rgba(0,0,0,.9);overflow:hidden;touch-action:none;';
  var hint=document.createElement('div');
  hint.textContent='__HINT__';
  hint.style.cssText='position:absolute;top:10px;left:50%;transform:translateX(-50%);color:#e8e8e8;font-size:13px;background:rgba(0,0,0,.55);padding:5px 14px;border-radius:14px;pointer-events:none;white-space:nowrap;z-index:3;max-width:92vw;overflow:hidden;text-overflow:ellipsis;';
  var close=document.createElement('div');
  close.textContent='\u00d7';
  close.style.cssText='position:absolute;top:6px;right:16px;color:#fff;font-size:34px;cursor:pointer;z-index:4;line-height:1;padding:4px 10px;';
  var view=document.createElement('div');
  view.style.cssText='position:absolute;inset:0;display:flex;align-items:center;justify-content:center;overflow:hidden;';
  var inner=document.createElement('div');
  inner.style.cssText='position:relative;transform-origin:center center;will-change:transform;';
  view.appendChild(inner); box.appendChild(view); box.appendChild(hint); box.appendChild(close);
  document.body.appendChild(box);
  var tip=document.createElement('div');
  tip.style.cssText='display:none;position:fixed;z-index:10002;max-width:72vw;padding:.5em .8em;background:rgba(20,20,20,.94);color:#fff;font-family:Consolas,Monaco,monospace;font-size:13px;line-height:1.5;border-radius:6px;pointer-events:none;white-space:pre-wrap;';
  document.body.appendChild(tip);
  var scale=1,tx=0,ty=0;
  function apply(){inner.style.transform='translate('+tx+'px,'+ty+'px) scale('+scale+')';box.style.cursor='grab';}
  function openBox(svg){
    var c=svg.cloneNode(true);
    c.removeAttribute('style');c.removeAttribute('class');
    c.setAttribute('width','760');c.setAttribute('height','760');
    c.style.cssText='width:760px;height:760px;display:block;background:#fff;border-radius:4px;';
    inner.innerHTML='';inner.appendChild(c);
    box.style.display='block';
    scale=Math.min(1,(window.innerWidth-40)/760,(window.innerHeight-90)/760);
    if(scale<0.15)scale=0.15; tx=0;ty=0; apply();
  }
  function closeBox(){box.style.display='none';tip.style.display='none';}
  close.addEventListener('click',function(e){e.stopPropagation();closeBox();});
  box.addEventListener('click',function(e){if(e.target===box||e.target===view)closeBox();});
  document.addEventListener('keydown',function(e){if(e.key==='Escape'&&box.style.display==='block')closeBox();});
  document.addEventListener('click',function(e){
    if(box.style.display!=='block')return;
    var r=e.target.closest('rect');
    if(r&&box.contains(r)){
      var t=r.querySelector('title');
      if(t){tip.textContent=t.textContent;
        tip.style.left=Math.min(e.clientX+14,window.innerWidth-280)+'px';
        tip.style.top=Math.min(e.clientY+14,window.innerHeight-70)+'px';
        tip.style.display='block';}
    } else if(!e.target.closest('svg')){tip.style.display='none';}
  });
  box.addEventListener('wheel',function(e){
    e.preventDefault();
    var f=e.deltaY<0?1.18:1/1.18,ns=Math.min(10,Math.max(0.15,scale*f));
    var r=view.getBoundingClientRect();
    var cx=e.clientX-r.left-r.width/2,cy=e.clientY-r.top-r.height/2;
    tx=cx-(cx-tx)*(ns/scale);ty=cy-(cy-ty)*(ns/scale);scale=ns;apply();
  },{passive:false});
  var dragging=false,sx=0,sy=0,px=0,py=0,moved=false;
  var actives={},pinch=false,pd=0;
  box.addEventListener('pointerdown',function(e){
    if(e.target.closest('rect')||e.target===close)return;
    dragging=true;moved=false;sx=e.clientX;sy=e.clientY;px=tx;py=ty;
    try{box.setPointerCapture(e.pointerId);}catch(_){}
    actives[e.pointerId]={x:e.clientX,y:e.clientY};
    if(Object.keys(actives).length===2){pinch=true;var p=Object.values(actives);pd=Math.hypot(p[0].x-p[1].x,p[0].y-p[1].y);}
  });
  box.addEventListener('pointermove',function(e){
    if(actives[e.pointerId])actives[e.pointerId]={x:e.clientX,y:e.clientY};
    if(pinch&&Object.keys(actives).length===2){
      var p=Object.values(actives),d=Math.hypot(p[0].x-p[1].x,p[0].y-p[1].y);
      if(pd>0){var ns=Math.min(10,Math.max(0.15,scale*d/pd));
        var cx2=(p[0].x+p[1].x)/2,cy2=(p[0].y+p[1].y)/2,r=view.getBoundingClientRect();
        var ax=cx2-r.left-r.width/2,ay=cy2-r.top-r.height/2;
        tx=ax-(ax-tx)*(ns/scale);ty=ay-(ay-ty)*(ns/scale);scale=ns;apply();}
      pd=d;return;
    }
    if(dragging){if(Math.abs(e.clientX-sx)+Math.abs(e.clientY-sy)>4)moved=true;tx=px+(e.clientX-sx);ty=py+(e.clientY-sy);apply();}
  });
  function endP(e){delete actives[e.pointerId];if(Object.keys(actives).length<2){pinch=false;pd=0;}dragging=false;}
  box.addEventListener('pointerup',endP);box.addEventListener('pointercancel',endP);
  inner.addEventListener('dblclick',function(e){
    e.preventDefault();
    var ns=Math.min(10,scale*1.8),r=view.getBoundingClientRect();
    var cx=e.clientX-r.left-r.width/2,cy=e.clientY-r.top-r.height/2;
    tx=cx-(cx-tx)*(ns/scale);ty=cy-(cy-ty)*(ns/scale);scale=ns;apply();
  });
  document.querySelectorAll('svg.zoomable').forEach(function(s){
    s.style.cursor='zoom-in';
    s.addEventListener('click',function(e){if(!moved)openBox(s);});
  });
})();
</script>'''

HINT_ZH = '滚轮/双指缩放 \u00b7 拖拽平移 \u00b7 双击放大 \u00b7 点 \u00d7 或按 Esc 关闭 \u00b7 点彩色格子读完整数值'
HINT_EN = 'Wheel/pinch to zoom \u00b7 drag to pan \u00b7 double-click to zoom in \u00b7 \u00d7 or Esc to close \u00b7 click a coloured cell for its full value'

def add_cell_text(svg):
    """为每个带 <title> 的有色 rect 插入白色数值 text。"""
    n=0
    def rep(m):
        nonlocal n
        attrs, inner = m.group(1), m.group(2)
        title = m.group(3)
        xv=re.search(r'\bx="([\d.]+)"',attrs); yv=re.search(r'\by="([\d.]+)"',attrs)
        wv=re.search(r'\bwidth="([\d.]+)"',attrs); hv=re.search(r'\bheight="([\d.]+)"',attrs)
        if not(xv and yv and wv and hv): return m.group(0)
        num=re.search(r'(-?\d+\.\d+[eE][+-]?\d+)',title)
        if not num: return m.group(0)
        v=float(num.group(1))
        lab='{:.2e}'.format(v).replace('e+0','e+').replace('e-0','e-')
        cx=float(xv.group(1))+float(wv.group(1))/2
        cy=float(yv.group(1))+float(hv.group(1))/2
        txt=('<text x="%.2f" y="%.2f" font-size="2.7" fill="#ffffff" '
             'text-anchor="middle" dominant-baseline="central" font-weight="bold" '
             'font-family="Consolas,Monaco,monospace" pointer-events="none">%s</text>'%(cx,cy,lab))
        n+=1
        return '<rect'+attrs+'>'+inner+txt+'</rect>'
    out=re.sub(r'<rect\b([^>]*)>(<title>([^<]*)</title>)</rect>', rep, svg)
    return out,n

def patch(fn, hint, guide_old, guide_new, guide2_old, guide2_new):
    s=io.open(fn,encoding='utf-8').read()
    # 1) 格子数值：逐 SVG 处理
    def svg_rep(m):
        svg=m.group(0)
        svg2,n=add_cell_text(svg)
        svg_rep.cnt+=n
        return svg2
    svg_rep.cnt=0
    s=re.sub(r'<svg\b[^>]*>.*?</svg>', svg_rep, s, flags=re.S)
    print(f"{fn}: cell texts inserted = {svg_rep.cnt}")
    assert svg_rep.cnt>=190, f"{fn}: 插入文本数异常 {svg_rep.cnt}"
    # 2) svg 加 class=zoomable
    c=s.count('viewBox="0 0 574 574"')
    assert c==2, f"{fn}: viewBox count {c}"
    s=s.replace('viewBox="0 0 574 574"','class="zoomable" viewBox="0 0 574 574"')
    # 3) 替换旧点按 JS 块
    m=re.search(r'<script>(?:(?!</script>).)*?closest\("rect"\).*?</script>', s, flags=re.S)
    assert m, f"{fn}: 旧点按JS块未找到"
    s=s[:m.start()]+LIGHTBOX_JS.replace('__HINT__',hint)+s[m.end():]
    # 4) 引导句
    assert s.count(guide_old)==1, f"{fn}: guide1 count {s.count(guide_old)}"
    s=s.replace(guide_old,guide_new)
    assert s.count(guide2_old)==1, f"{fn}: guide2 count {s.count(guide2_old)}"
    s=s.replace(guide2_old,guide2_new)
    io.open(fn,'w',encoding='utf-8').write(s)
    print(f"{fn}: OK")

patch('zh.html', HINT_ZH,
      '手机上可左右滑动、双指放大；点按（或悬停）任一彩色格子可读出该矩阵元的精确数值：',
      '点击图片即可放大：电脑上滚轮缩放、拖拽平移，手机上双指缩放；放大后每个彩色格子内直接标注数值，点按（或悬停）格子还可读出完整精确值：',
      '同样可滑动、放大、点按读数：',
      '同样点击图片放大：格子内直接标注数值，点按格子读完整精确值：')

patch('index.html', HINT_EN,
      'Scroll sideways and pinch-zoom on a phone; tap (or hover) any coloured cell to read the exact value of that entry:',
      'Click the figure to zoom: mouse-wheel to zoom and drag to pan on a computer (pinch on phones); values are printed inside every coloured cell, and tapping (or hovering) a cell shows its full exact value:',
      'Scroll, pinch-zoom and tap for values:',
      'Click to zoom here as well: values are printed inside each cell, and tapping a cell shows its full value:')
print("ROUND 10 DONE")
