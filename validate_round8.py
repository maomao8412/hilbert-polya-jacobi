# -*- coding: utf-8 -*-
from html.parser import HTMLParser
import io

VOID = {"meta","link","img","br","hr","input","area","base","col","embed","source","track","wbr","rect","path","circle","line","polyline","polygon","ellipse","use"}

class P(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack=[]; self.errors=[]
    def handle_starttag(self, tag, attrs):
        if tag not in VOID:
            self.stack.append((tag, self.getpos()))
    def handle_startendtag(self, tag, attrs):
        pass
    def handle_endtag(self, tag):
        if tag in VOID: return
        # tolerate implicit closes of p/li/td/th/option
        IMPLICIT = {"p","li","td","th","tr","option","thead","tbody","dt","dd"}
        if self.stack and self.stack[-1][0]==tag:
            self.stack.pop()
        elif tag in IMPLICIT:
            # pop until match, allowing implicit closures
            popped=[]
            while self.stack and self.stack[-1][0]!=tag:
                t,pos=self.stack.pop(); popped.append(t)
                if t not in IMPLICIT:
                    self.errors.append(("mismatch-close",tag,"while open",t,pos)); return
            if self.stack: self.stack.pop()
        else:
            self.errors.append(("stray-close",tag,self.getpos()))

for fn in ["zh.html","index.html"]:
    s=io.open(fn,encoding="utf-8").read()
    p=P(); p.feed(s)
    unclosed=[t for t,_ in p.stack if t not in ("html","body")]
    print(fn, "errors:", p.errors[:5], "| unclosed(non html/body):", unclosed[:5])
    # key assertions
    checks = {
      "center block": 'text-align:center;max-width:920px' in s,
      "no Chinese name 陈倬": "陈倬" not in s,
      "sans font": "Segoe UI" in s and "Georgia" not in s,
      "mathjax nudge": 'mjx-container[display="false"]{vertical-align:-0.06em}' in s,
      "gold legend": ('金色' if fn=="zh.html" else 'gold') + '"' or True,
      "touch-action": s.count("touch-action:pan-x pan-y")==2,
      "svg min-width 760": s.count("min-width:760px")==2,
      "no 520px": "min-width:520px" not in s,
      "tap JS": 'closest("rect")' in s and "querySelector(\"title\")" in s,
      "RH label": ("RH 篇" if fn=="zh.html" else "RH paper") in s,
      "no word proof zh/en": ("证明" not in s.split("Zenodo")[0] if fn=="zh.html" else True),
    }
    if fn=="zh.html":
        checks["gold legend"]=('<span style="background:#fff8e1;color:#856404' in s and '金色</span> = 临近锁定' in s)
        checks["grey legend"]=('灰色</span> = 仍在收敛' in s)
        checks["dash legend"]=('= 尚未出现。' in s and 'background:#f0f0f0' in s)
        checks["copyright text"]="本页所述工作由 <strong>Zhuo Chen</strong>" in s
        # ensure 证明 not in copyright block region
        import re
        m=re.search(r"版权与引用.*?</p>", s, re.S)
        checks["no 证明 in block"]= m is not None and "证明" not in m.group(0)
    else:
        checks["gold legend"]=('<span style="background:#fff8e1;color:#856404' in s and 'gold</span> = near-locked' in s)
        checks["grey legend"]=('grey</span> = still converging' in s)
        checks["dash legend"]=('= not yet present.' in s and 'background:#f0f0f0' in s)
        checks["copyright text"]="The work presented here is by <strong>Zhuo Chen</strong>" in s
        import re
        m=re.search(r"Copyright &amp; citation.*?</p>", s, re.S)
        checks["no proof in block"]= m is not None and "proof" not in m.group(0).lower()
    bad=[k for k,v in checks.items() if not v]
    print("   checks:", "ALL OK" if not bad else ("FAIL: "+str(bad)))
