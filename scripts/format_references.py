from pathlib import Path
import re,json
M=Path(__file__).resolve().parents[1]/'manuscript'
s=(M/'references.bib').read_text();entries=[]
for m in re.finditer(r'@(\w+)\{([^,]+),',s):
 i=m.end();depth=1;j=i
 while j<len(s) and depth:
  if s[j]=='{' and s[j-1]!='\\':depth+=1
  if s[j]=='}' and s[j-1]!='\\':depth-=1
  j+=1
 body=s[i:j-1];fields={};pos=0
 while True:
  fm=re.search(r'(\w+)\s*=\s*\{',body[pos:])
  if not fm:break
  a=pos+fm.end();d=1;b=a
  while b<len(body) and d:
   if body[b]=='{' and body[b-1]!='\\':d+=1
   if body[b]=='}' and body[b-1]!='\\':d-=1
   b+=1
  fields[fm.group(1).lower()]=body[a:b-1];pos=b
 entries.append((m.group(1),m.group(2),fields))
def initials(given):
 parts=given.split();result=[]
 for v in parts:
  if '.' in v:result.append(v)
  else:
   ch=re.search('[A-Za-z]',v).group();result.append(ch+'.')
 return '~'.join(result)
def name(a):
 family,given=a.split(',',1);return family.strip()+', '+initials(given.strip())
def family(a):return a.split(',')[0].strip()
def joined(xs):
 if len(xs)==1:return xs[0]
 if len(xs)==2:return xs[0]+r', \& '+xs[1]
 return ', '.join(xs[:-1])+r', \& '+xs[-1]
out=['% Publisher name-year / alphabetical style. Generated from references.bib.','\\begin{thebibliography}{99}']
entries.sort(key=lambda e:re.sub(r'[^a-z]','',e[2]['author'].split(' and ')[0].lower()))
for kind,key,f in entries:
 authors=f['author'].split(' and ');fam=[family(a) for a in authors]
 short=(fam[0] if len(fam)==1 else ' and '.join(fam) if len(fam)==2 else fam[0]+' et~al.')
 title=f['title'];title=title.replace('A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems','A fast iterative shrinkage-thresholding algorithm for linear inverse problems').replace('Proximal Algorithms','Proximal algorithms')
 text=joined([name(a) for a in authors])+' ('+f['year']+'). '
 if kind=='article':text+=title+'. '+r'\textit{'+f['journal']+'}, '+r'\textit{'+f['volume']+'}'+('('+f['number']+')' if 'number'in f else '')+', '+f['pages']+'.'
 elif kind=='book':text+=r'\textit{'+title+'}. '+f['publisher']+'.'
 else:
  eds=f['editor'].split(' and ');ednames=[initials(a.split(',')[1].strip())+' '+family(a) for a in eds]
  text+=title+'. In '+joined(ednames)+r' (Eds.), \textit{'+f['booktitle']+'} (pp. '+f['pages']+'). '+f['publisher']+'.'
 text+=' '+r'\url{https://doi.org/'+f['doi']+'}'
 out+=['',r'\bibitem['+short+'('+f['year']+')]'+r'{'+key+'}',text]
out+=['',r'\end{thebibliography}']
(M/'references_frozen.bbl').write_text('\n'.join(out)+'\n')
print('references',len(entries))

# Online Resource 1 contains only references cited in its own text.
cited=set()
for p in (M/'appendices').glob('*.tex'):
 for v in re.findall(r'\\cite\w*\*?(?:\[[^]]*\])*\{([^}]+)\}',p.read_text()): cited.update(v.split(','))
bbl=(M/'references_frozen.bbl').read_text(); parts=[]
for part in re.split(r'(?=\\bibitem\[)',bbl)[1:]:
 key=re.search(r'\\bibitem\[[^]]*\]\{([^}]+)\}',part).group(1)
 if key in cited: parts.append(part.split(r'\end{thebibliography}')[0].rstrip())
(M/'references_supplement.bbl').write_text('\\begin{thebibliography}{99}\n\n'+'\n\n'.join(parts)+'\n\n\\end{thebibliography}\n')
