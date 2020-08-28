r2=3
r1=6
for i in range(3,8):
  r2n=r1
  r1n=r2*2+r1*2
  r1=r1n
  r2=r2n
  print(r1+r2,i)
print(r1+r2)

a=['R','G','B']
s=[]
for a1 in a:
    for a2 in a:
        for a3 in a:
            for a4 in a:
                for a5 in a:
                    for a6 in a:
                        for a7 in a:
                            s.append([a1,a2,a3,a4,a5,a6,a7])
bad=0
for s1 in s:
    ibad=0
    for k in range(7-2):
        if s1[k]==s1[k+1] and s1[k]==s1[k+2] and s1[k+1]==s1[k+2]:
            ibad=1
    bad+=ibad
print(3**7-bad)
