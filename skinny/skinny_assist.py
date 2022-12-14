from cao_pylib import sbox
from cao_pylib import sbox_model
from cao_pylib import ine
from gurobipy import *
from cao_pylib import operation

def int2tri(width,n):
    #将10进制数转为3进制数
    assert 3**width-1>=n,"int2tri:number too big"
    res=[0 for i in range(width)]

    for i in range(width-1,-1,-1):
        res[i]=n%3
        n=int(n/3)
        if n<3:
            res[i-1]=n
            break
    return res

def trilist2binlist(x):
    res=[]
    for i in range(len(x)):
        if x[i]==0:
            res+=[0,0]
        if x[i]==1:
            res+=[0,1]
        if x[i]==2:
            res+=[1,1]
    return res


def gen_ine_from_cnf(cnf):
    """
    从cnf生成不等式，并存入best_ine文件,文字的命名：0,1，2，3...
    """
    vars=set()
    all_clause=cnf[1:len(cnf)-1].split(')(')
    vars={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
    ine_all=[]
    for clause in all_clause:
        liter=clause.split('+')
        const=-1
        tmp=[0 for i in range(len(vars)+2)]
        for index in vars:
            if str(index) in liter:
                tmp[index]=1
            if str(index)+"'" in liter:
                const+=1
                tmp[index]=-1
        tmp[len(vars)]=const
        ine_all.append(tmp)
    return ine_all

def gen_ine_cnf(cnf,best_ine):
    ine0=gen_ine_from_cnf(cnf[0])
    ine1=gen_ine_from_cnf(cnf[1])

    for i in range(len(ine0)):
        print(ine0[i])
        ine0[i]=ine0[i][:len(ine0[i])-2]+[17]+ine0[i][len(ine0[i])-2:]
        print(ine0[i])
    for i in range(len(ine1)):
        ine1[i]=ine1[i][:len(ine1[i])-2]+[-17]+[ine1[i][len(ine1[i]) - 2 ]+17]+[ine1[i][len(ine1[i])-1]]
    ine_all=ine0+ine1

    fp=open(best_ine,'w+')
    res=''
    for line in ine_all:
        s=ine.list2ine(line)
        s+='\n'
        res+=s
    fp.write(res[:len(res)-1])
    fp.close()

#####################################
'''
sbox_bit=[]
fp=open('skinny//skinny_sbox_8.csv')
s=fp.read()
fp.close()
s=s.split('\n')[1:]

for i in range(len(s)):
    s[i]=s[i].split(',')[:32]
    for j in range(32):
        s[i][j]=eval(s[i][j])
for i in range(16):
    tmp=[]
    for j in range(len(s)):
        tmp.append(s[j][:16]+[s[j][i+16]])
    sbox_bit.append(tmp)
#print(sbox_bit[0])
'''
######################################

cnf1=["(1+3+4+7+9+11+15)(12'+13)(4'+5)(14'+15)(10'+11)(8'+9)(0'+1)(2'+3)(6'+7)",
     "(1+3+5+7+9+11+15)(0'+1)(2'+3)(4'+5)(6'+7)(8'+9)(10'+11)(14'+15)(12'+13)",
     "(1+3+6)(0'+1)(2'+3)(4'+5)(8'+9)(10'+11)(12'+13)(14'+15)(6'+7)",
     "(1+3+7)(14'+15)(0'+1)(2'+3)(6'+7)(12'+13)(4'+5)(8'+9)(10'+11)",
     "(9+11+14)(8'+9)(10'+11)(0'+1)(2'+3)(4'+5)(6'+7)(12'+13)(14'+15)",
     "(9+11+15)(12'+13)(8'+9)(6'+7)(0'+1)(2'+3)(4'+5)(10'+11)(14'+15)",
     "(1+3+5+7+9+11+15)(1+3+4+5'+7+9+11+14+15')(12'+13)(4'+5)(14'+15)(10'+11)(8'+9)(0'+1)(2'+3)(6'+7)",
     "(1+3+5+7+9+11+15)(1+3+4+5'+7+9+11+14+15')(12'+13)(4'+5)(14'+15)(10'+11)(8'+9)(0'+1)(2'+3)(6'+7)",
     "(9+11+12+15)(9+10+11'+12+14+15')(8'+9)(0'+1)(2'+3)(4'+5)(6'+7)(10'+11)(12'+13)(14'+15)",
     "(9+11+13+15)(9+10+11'+13+14+15')(8'+9)(12'+13)(0'+1)(2'+3)(4'+5)(6'+7)(10'+11)(14'+15)",
     "(2+11+13)(10'+11)(12'+13)(0'+1)(4'+5)(6'+7)(8'+9)(14'+15)(2'+3)",
     "(3+11+13)(14'+15)(2'+3)(8'+9)(0'+1)(4'+5)(6'+7)(10'+11)(12'+13)",
     "(1+3+5+7+9+11+13+15)(0'+1)(2'+3)(4'+5)(6'+7)(8'+9)(10'+11)(12'+13)(14'+15)",
     "(1+3+5+7+9+11+13+15)(0'+1)(2'+3)(4'+5)(6'+7)(8'+9)(10'+11)(12'+13)(14'+15)",
     "(1+3+5+7+9+11+13+15)(0'+1)(2'+3)(4'+5)(6'+7)(8'+9)(10'+11)(12'+13)(14'+15)",
     "(1+3+5+7+9+11+13+15)(0'+1)(2'+3)(4'+5)(6'+7)(8'+9)(10'+11)(12'+13)(14'+15)"]

cnf0=["(15')(14')(11')(10')(9')(8')(7')(6')(4')(3')(2')(1')(0')(12'+13)",
      "(15')(14')(11')(10')(9')(8')(7')(6')(5')(4')(3')(2')(1')(0')(12'+13)",
      "(6')(3')(2')(1')(0')(4'+5)(8'+9)(10'+11)(12'+13)(14'+15)",
      "(7')(6')(3')(2')(1')(0')(4'+5)(8'+9)(10'+11)(12'+13)(14'+15)",
      "(14')(11')(10')(9')(8')(0'+1)(2'+3)(4'+5)(6'+7)(12'+13)",
      "(15')(14')(11')(10')(9')(8')(0'+1)(2'+3)(4'+5)(6'+7)(12'+13)",
      "(14')(11')(10')(9')(8')(7')(6')(4')(3')(2')(1')(0')(12'+13)(5+15')(5'+15)",
      "(14')(11')(10')(9')(8')(7')(6')(4')(3')(2')(1')(0')(12'+13)(5+15')(5'+15)",
      "(14')(12')(10')(9')(8')(0'+1)(2'+3)(4'+5)(6'+7)(11+15')(11'+15)",
      "(14')(13')(12')(10')(9')(8')(0'+1)(2'+3)(4'+5)(6'+7)(11+15')(11'+15)",
      "(13')(12')(11')(10')(2')(0'+1)(4'+5)(6'+7)(8'+9)(14'+15)",
      "(13')(12')(11')(10')(3')(2')(0'+1)(4'+5)(6'+7)(8'+9)(14'+15)",
      "(15')(14')(13')(12')(11')(10')(9')(8')(7')(6')(5')(4')(3')(2')(1')(0')",
      "(15')(14')(13')(12')(11')(10')(9')(8')(7')(6')(5')(4')(3')(2')(1')(0')",
      "(15')(14')(13')(12')(11')(10')(9')(8')(7')(6')(5')(4')(3')(2')(1')(0')",
      "(15')(14')(13')(12')(11')(10')(9')(8')(7')(6')(5')(4')(3')(2')(1')(0')"]

active_cnf0=["(0'+2'+4'+6'+8'+10'+12'+14')(1+2'+4'+6'+8'+10'+12'+14')(0'+3+4'+6'+8'+10'+12'+14')(0'+2'+5+6'+8'+10'+12'+14')(0'+2'+4'+7+8'+10'+12'+14')(0'+2'+4'+6'+9+10'+12'+14')(0'+2'+4'+6'+8'+11+12'+14')(0'+2'+4'+6'+8'+10'+13+14')(0'+2'+4'+6'+8'+10'+12'+15)(1+3+4'+6'+8'+10'+12'+14')(1+2'+5+6'+8'+10'+12'+14')(0'+3+5+6'+8'+10'+12'+14')(1+2'+4'+7+8'+10'+12'+14')(0'+3+4'+7+8'+10'+12'+14')(0'+2'+5+7+8'+10'+12'+14')(1+2'+4'+6'+9+10'+12'+14')(0'+3+4'+6'+9+10'+12'+14')(0'+2'+5+6'+9+10'+12'+14')(0'+2'+4'+7+9+10'+12'+14')(1+2'+4'+6'+8'+11+12'+14')(0'+3+4'+6'+8'+11+12'+14')(0'+2'+5+6'+8'+11+12'+14')(0'+2'+4'+7+8'+11+12'+14')(0'+2'+4'+6'+9+11+12'+14')(1+2'+4'+6'+8'+10'+13+14')(0'+3+4'+6'+8'+10'+13+14')(0'+2'+5+6'+8'+10'+13+14')(0'+2'+4'+7+8'+10'+13+14')(0'+2'+4'+6'+9+10'+13+14')(0'+2'+4'+6'+8'+11+13+14')(1+2'+4'+6'+8'+10'+12'+15)(0'+3+4'+6'+8'+10'+12'+15)(0'+2'+5+6'+8'+10'+12'+15)(0'+2'+4'+7+8'+10'+12'+15)(0'+2'+4'+6'+9+10'+12'+15)(0'+2'+4'+6'+8'+11+12'+15)(0'+2'+4'+6'+8'+10'+13+15)(1+3+5+6'+8'+10'+12'+14')(1+3+4'+7+8'+10'+12'+14')(1+2'+5+7+8'+10'+12'+14')(0'+3+5+7+8'+10'+12'+14')(1+3+4'+6'+9+10'+12'+14')(1+2'+5+6'+9+10'+12'+14')(0'+3+5+6'+9+10'+12'+14')(1+2'+4'+7+9+10'+12'+14')(0'+3+4'+7+9+10'+12'+14')(0'+2'+5+7+9+10'+12'+14')(1+3+4'+6'+8'+11+12'+14')(1+2'+5+6'+8'+11+12'+14')(0'+3+5+6'+8'+11+12'+14')(1+2'+4'+7+8'+11+12'+14')(0'+3+4'+7+8'+11+12'+14')(0'+2'+5+7+8'+11+12'+14')(1+2'+4'+6'+9+11+12'+14')(0'+3+4'+6'+9+11+12'+14')(0'+2'+5+6'+9+11+12'+14')(0'+2'+4'+7+9+11+12'+14')(1+3+4'+6'+8'+10'+13+14')(1+2'+5+6'+8'+10'+13+14')(0'+3+5+6'+8'+10'+13+14')(1+2'+4'+7+8'+10'+13+14')(0'+3+4'+7+8'+10'+13+14')(0'+2'+5+7+8'+10'+13+14')(1+2'+4'+6'+9+10'+13+14')(0'+3+4'+6'+9+10'+13+14')(0'+2'+5+6'+9+10'+13+14')(0'+2'+4'+7+9+10'+13+14')(1+2'+4'+6'+8'+11+13+14')(0'+3+4'+6'+8'+11+13+14')(0'+2'+5+6'+8'+11+13+14')(0'+2'+4'+7+8'+11+13+14')(0'+2'+4'+6'+9+11+13+14')(1+3+4'+6'+8'+10'+12'+15)(1+2'+5+6'+8'+10'+12'+15)(0'+3+5+6'+8'+10'+12'+15)(1+2'+4'+7+8'+10'+12'+15)(0'+3+4'+7+8'+10'+12'+15)(0'+2'+5+7+8'+10'+12'+15)(1+2'+4'+6'+9+10'+12'+15)(0'+3+4'+6'+9+10'+12'+15)(0'+2'+5+6'+9+10'+12'+15)(0'+2'+4'+7+9+10'+12'+15)(1+2'+4'+6'+8'+11+12'+15)(0'+3+4'+6'+8'+11+12'+15)(0'+2'+5+6'+8'+11+12'+15)(0'+2'+4'+7+8'+11+12'+15)(0'+2'+4'+6'+9+11+12'+15)(1+2'+4'+6'+8'+10'+13+15)(0'+3+4'+6'+8'+10'+13+15)(0'+2'+5+6'+8'+10'+13+15)(0'+2'+4'+7+8'+10'+13+15)(0'+2'+4'+6'+9+10'+13+15)(0'+2'+4'+6'+8'+11+13+15)(1+3+5+7+8'+10'+12'+14')(1+3+5+6'+9+10'+12'+14')(1+3+4'+7+9+10'+12'+14')(1+2'+5+7+9+10'+12'+14')(0'+3+5+7+9+10'+12'+14')(1+3+5+6'+8'+11+12'+14')(1+3+4'+7+8'+11+12'+14')(1+2'+5+7+8'+11+12'+14')(0'+3+5+7+8'+11+12'+14')(1+3+4'+6'+9+11+12'+14')(1+2'+5+6'+9+11+12'+14')(0'+3+5+6'+9+11+12'+14')(1+2'+4'+7+9+11+12'+14')(0'+3+4'+7+9+11+12'+14')(0'+2'+5+7+9+11+12'+14')(1+3+5+6'+8'+10'+13+14')(1+3+4'+7+8'+10'+13+14')(1+2'+5+7+8'+10'+13+14')(0'+3+5+7+8'+10'+13+14')(1+3+4'+6'+9+10'+13+14')(1+2'+5+6'+9+10'+13+14')(0'+3+5+6'+9+10'+13+14')(1+2'+4'+7+9+10'+13+14')(0'+3+4'+7+9+10'+13+14')(0'+2'+5+7+9+10'+13+14')(1+3+4'+6'+8'+11+13+14')(1+2'+5+6'+8'+11+13+14')(0'+3+5+6'+8'+11+13+14')(1+2'+4'+7+8'+11+13+14')(0'+3+4'+7+8'+11+13+14')(0'+2'+5+7+8'+11+13+14')(1+2'+4'+6'+9+11+13+14')(0'+3+4'+6'+9+11+13+14')(0'+2'+5+6'+9+11+13+14')(0'+2'+4'+7+9+11+13+14')(1+3+5+6'+8'+10'+12'+15)(1+3+4'+7+8'+10'+12'+15)(1+2'+5+7+8'+10'+12'+15)(0'+3+5+7+8'+10'+12'+15)(1+3+4'+6'+9+10'+12'+15)(1+2'+5+6'+9+10'+12'+15)(0'+3+5+6'+9+10'+12'+15)(1+2'+4'+7+9+10'+12'+15)(0'+3+4'+7+9+10'+12'+15)(0'+2'+5+7+9+10'+12'+15)(1+3+4'+6'+8'+11+12'+15)(1+2'+5+6'+8'+11+12'+15)(0'+3+5+6'+8'+11+12'+15)(1+2'+4'+7+8'+11+12'+15)(0'+3+4'+7+8'+11+12'+15)(0'+2'+5+7+8'+11+12'+15)(1+2'+4'+6'+9+11+12'+15)(0'+3+4'+6'+9+11+12'+15)(0'+2'+5+6'+9+11+12'+15)(0'+2'+4'+7+9+11+12'+15)(1+3+4'+6'+8'+10'+13+15)(1+2'+5+6'+8'+10'+13+15)(0'+3+5+6'+8'+10'+13+15)(1+2'+4'+7+8'+10'+13+15)(0'+3+4'+7+8'+10'+13+15)(0'+2'+5+7+8'+10'+13+15)(1+2'+4'+6'+9+10'+13+15)(0'+3+4'+6'+9+10'+13+15)(0'+2'+5+6'+9+10'+13+15)(0'+2'+4'+7+9+10'+13+15)(1+2'+4'+6'+8'+11+13+15)(0'+3+4'+6'+8'+11+13+15)(0'+2'+5+6'+8'+11+13+15)(0'+2'+4'+7+8'+11+13+15)(0'+2'+4'+6'+9+11+13+15)(1+3+5+7+9+10'+12'+14')(1+3+5+7+8'+11+12'+14')(1+3+5+6'+9+11+12'+14')(1+3+4'+7+9+11+12'+14')(1+2'+5+7+9+11+12'+14')(0'+3+5+7+9+11+12'+14')(1+3+5+7+8'+10'+13+14')(1+3+5+6'+9+10'+13+14')(1+3+4'+7+9+10'+13+14')(1+2'+5+7+9+10'+13+14')(0'+3+5+7+9+10'+13+14')(1+3+5+6'+8'+11+13+14')(1+3+4'+7+8'+11+13+14')(1+2'+5+7+8'+11+13+14')(0'+3+5+7+8'+11+13+14')(1+3+4'+6'+9+11+13+14')(1+2'+5+6'+9+11+13+14')(0'+3+5+6'+9+11+13+14')(1+2'+4'+7+9+11+13+14')(0'+3+4'+7+9+11+13+14')(0'+2'+5+7+9+11+13+14')(1+3+5+7+8'+10'+12'+15)(1+3+5+6'+9+10'+12'+15)(1+3+4'+7+9+10'+12'+15)(1+2'+5+7+9+10'+12'+15)(0'+3+5+7+9+10'+12'+15)(1+3+5+6'+8'+11+12'+15)(1+3+4'+7+8'+11+12'+15)(1+2'+5+7+8'+11+12'+15)(0'+3+5+7+8'+11+12'+15)(1+3+4'+6'+9+11+12'+15)(1+2'+5+6'+9+11+12'+15)(0'+3+5+6'+9+11+12'+15)(1+2'+4'+7+9+11+12'+15)(0'+3+4'+7+9+11+12'+15)(0'+2'+5+7+9+11+12'+15)(1+3+5+6'+8'+10'+13+15)(1+3+4'+7+8'+10'+13+15)(1+2'+5+7+8'+10'+13+15)(0'+3+5+7+8'+10'+13+15)(1+3+4'+6'+9+10'+13+15)(1+2'+5+6'+9+10'+13+15)(0'+3+5+6'+9+10'+13+15)(1+2'+4'+7+9+10'+13+15)(0'+3+4'+7+9+10'+13+15)(0'+2'+5+7+9+10'+13+15)(1+3+4'+6'+8'+11+13+15)(1+2'+5+6'+8'+11+13+15)(0'+3+5+6'+8'+11+13+15)(1+2'+4'+7+8'+11+13+15)(0'+3+4'+7+8'+11+13+15)(0'+2'+5+7+8'+11+13+15)(1+2'+4'+6'+9+11+13+15)(0'+3+4'+6'+9+11+13+15)(0'+2'+5+6'+9+11+13+15)(0'+2'+4'+7+9+11+13+15)(1+3+5+7+9+11+12'+14')(1+3+5+7+9+10'+13+14')(1+3+5+7+8'+11+13+14')(1+3+5+6'+9+11+13+14')(1+3+4'+7+9+11+13+14')(1+2'+5+7+9+11+13+14')(0'+3+5+7+9+11+13+14')(1+3+5+7+9+10'+12'+15)(1+3+5+7+8'+11+12'+15)(1+3+5+6'+9+11+12'+15)(1+3+4'+7+9+11+12'+15)(1+2'+5+7+9+11+12'+15)(0'+3+5+7+9+11+12'+15)(1+3+5+7+8'+10'+13+15)(1+3+5+6'+9+10'+13+15)(1+3+4'+7+9+10'+13+15)(1+2'+5+7+9+10'+13+15)(0'+3+5+7+9+10'+13+15)(1+3+5+6'+8'+11+13+15)(1+3+4'+7+8'+11+13+15)(1+2'+5+7+8'+11+13+15)(0'+3+5+7+8'+11+13+15)(1+3+4'+6'+9+11+13+15)(1+2'+5+6'+9+11+13+15)(0'+3+5+6'+9+11+13+15)(1+2'+4'+7+9+11+13+15)(0'+3+4'+7+9+11+13+15)(0'+2'+5+7+9+11+13+15)(1+3+5+7+9+11+13+14')(1+3+5+7+9+11+12'+15)(1+3+5+7+9+10'+13+15)(1+3+5+7+8'+11+13+15)(1+3+5+6'+9+11+13+15)(1+3+4'+7+9+11+13+15)(1+2'+5+7+9+11+13+15)(0'+3+5+7+9+11+13+15)(14'+15)(12'+13)(0'+1)(2'+3)(4'+5)(6'+7)(8'+9)(10'+11)","(15')(14')(13')(12')(11')(10')(9')(8')(7')(6')(5')(4')(3')(2')(1')(0')"]
active_cnf1=["(0+1')(0'+1)(2+3')(2'+3)(4+5')(4'+5)(6+7')(6'+7)(8+9')(8'+9)(10+11')(10'+11)(12+13')(12'+13)(14+15')(14'+15)(1+3+5+7+9+11+13+15)","(1+3+5+7+9+11+13+15)(0'+1)(2'+3)(4'+5)(6'+7)(8'+9)(10'+11)(12'+13)(14'+15)"]

'''
for i in range(16):
    gen_ine_cnf([cnf0[i],cnf1[i]],'skinny//skinny'+str(i)+'.txt')
#test=gen_ine_from_cnf(cnf0[0],'skinny0.txt',1)
#print(test)
'''
res=[]
for i in range(3**8):
    itri=int2tri(8,i)
    if i==0:
        res.append(itri+[0])
    if 1 in itri:
        res.append(itri+[1])
    else:
        res.append(itri+[2])
for i in res:
    print(i)
