from pathlib import Path
import re
p=Path(__file__).resolve().parent
for name in ["main","supplement"]:
    source=p/(name+".aux")
    if source.exists():
        lines=[s for s in source.read_text().splitlines() if s.startswith(r"\newlabel{")]
        (p/("xrefs_"+name+".aux")).write_text("\n".join(lines)+"\n")
