from gurobipy import *
from cao_pylib import operation
from cao_pylib import sbox_model
import copy 
import sys

def print_table(x):
    for i in x:
        print(i)

def _state(x):
    res=[]
    for row in range(5):
        tmpx=[]
        for col in range(64):
            tmpx.append([x[(row*64+col)*2],x[(row*64+col)*2+1]])
        res.append(tmpx)
    return res


def print_state(x):
    for i in x:
        s=''
        for j in i:
            s+=str(j)
        print(s)

def _get_value(ascon):
    xnum=len(ascon.x)
    snum=len(ascon.sout)
    rounds=int(snum/640)
    x=[0 for i in range(xnum)]
    sout=[0 for i in range(snum)]
    ix=[0 for i in range(640)]

    for v in ascon.m.getVars():
        if v.varName[0]=='x':
            x[eval(v.varName[2:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='s':
            sout[eval(v.varName[5:len(v.varName)-1])]=int(v.x)
        if v.varName[0]=='i':
            ix[eval(v.varName[3:len(v.varName)-1])]=int(v.x)
    
    iix=[]
    for row in range(5):
        tmpiix=[]
        for col in range(64):
            tmp=[ix[(row*64+col)*2],ix[(row*64+col)*2+1]]
            if tmp==[0,0]:
                tmpiix.append(0)
            if tmp==[0,1]:
                tmpiix.append(1)
            if tmp==[1,1]:
                tmpiix.append('?')
        iix.append(tmpiix)
    print_state(iix)
    print('--------------------------------------------------------')
    
    xx=[]
    for r in range(int(xnum/640)):
        tmpx=[]
        for row in range(5):
            tmpxx=[]
            for col in range(64):
                tmp=[x[r*640+(row*64+col)*2],x[r*640+(row*64+col)*2+1]]
                if tmp==[0,0]:
                    tmpxx.append(0)
                if tmp==[0,1]:
                    tmpxx.append(1)
                if tmp==[1,1]:
                    tmpxx.append('?')
            tmpx.append(tmpxx)
        xx.append(tmpx)
    
    ss=[]
    for r in range(int(snum/640)):
        tmps=[]
        for row in range(5):
            tmpss=[]
            for col in range(64):
                tmp=[sout[r*640+(row*64+col)*2],sout[r*640+(row*64+col)*2+1]]
                if tmp==[0,0]:
                    tmpss.append(0)
                if tmp==[0,1]:
                    tmpss.append(1)
                if tmp==[1,1]:
                    tmpss.append('?')
            tmps.append(tmpss)
        ss.append(tmps)
    
    for r in range(rounds):
        print_state(xx[r])
        print()
        print_state(ss[r])
        print()
        print()
    if xnum>snum:
        print_state(xx[rounds])

def inv_xor_constr(m,x,y,t):
    remain_x=(len(x)-1)%3
    ine_num=int((len(x)-1)/3)

    tx=[]
    ty=[]
    for pos in range(4):
        tmpx=[x[pos]]
        if pos==0:
            for i in range(ine_num-1):
                tmpx.append(t[i])
        else:
            for i in range(1,ine_num):
                tmpx.append(x[i*3+pos])
        tx.append(tmpx)
    
    for i in range(ine_num-1):
        ty.append(t[i])
    if remain_x==0:
        ty.append(y)
    if remain_x!=0:
        ty.append(t[ine_num-1])
    
    for i in range(len(tx[0])):
        #print(tx[0][i]+tx[1][i]+tx[2][i]+tx[3][i]+ty[i])
        sbox_model.gen_model_from_ine(m,tx[0][i]+tx[1][i]+tx[2][i]+tx[3][i]+ty[i],'ine/_4xor_ine.txt')
    
    xx=[]
    yy=[]
    if remain_x==2:
        xx.append(t[ine_num-1])
        xx.append(x[len(x)-2])
        xx.append(x[len(x)-1])
        yy.append(y)
        sbox_model.gen_model_from_ine(m,xx[0]+xx[1]+xx[2]+yy[0],'ine/3xor_ine.txt')
        #print(xx[0]+xx[1]+xx[2]+yy[0])
    if remain_x==1:
        xx.append(t[ine_num-1])
        xx.append(x[len(x)-1])
        yy.append(y)
        sbox_model.gen_model_from_ine(m,xx[0]+xx[1]+yy[0],'ine/xor_ine.txt')
        #print(xx[0]+xx[1]+yy[0])

class Ascon():
    def __init__(self,round,block_size,index,value):
        self.r=round
        self.block_size=block_size
        #self.flag=flag  #flag=1,线性层多于非线性层(L->| SL |->| SL |->| SL |(3轮))，flag=0,线性层等于非线性层(L->| SL |->| SL |->| S |(3轮))
        self.flag=1
        self.index=index
        if value==1:
            self.value=[0,1]
        if value==0:
            self.value=[0,0]
        self.rotate_right=[[0,19,28],[0,61,39],[0,1,6],[0,10,17],[0,7,41]]
        
        self.m=Model()
        if self.flag==0:
            self.x=self.m.addVars([i for i in  range(640*(self.r+1))],vtype=GRB.BINARY,name='x')
            self.sout=self.m.addVars([i for i in  range(640*(self.r))],vtype=GRB.BINARY,name='sout')
        if self.flag==1:
            self.x=self.m.addVars([i for i in  range(640*(self.r))],vtype=GRB.BINARY,name='x')
            self.sout=self.m.addVars([i for i in  range(640*self.r)],vtype=GRB.BINARY,name='sout')
        self.q=self.m.addVars([i for i in range(320)],vtype=GRB.BINARY,name='q')  #输入中undisturbed bit的表示
        self.ix=self.m.addVars([i for  i in range(640)],vtype=GRB.BINARY,name='ix') #第一轮线性层的输入
        
        constr_ix=LinExpr()
        for i in range(64,320):
            constr_ix+=self.q[i]
        self.m.setObjective(constr_ix,GRB.MAXIMIZE)
        
        #输入和undisturbed 数量之间的约束
        for i in range(320):
            self.m.addConstr(self.ix[i*2]>=self.q[i])
            self.m.addConstr(self.ix[i*2+1]>=self.q[i])
            self.m.addConstr(self.ix[i*2]+self.x[i*2+1]-self.q[i]<=1)

    def gen_constr(self):
        #第一个线性层
        for row in range(5):
            inputx=[]
            outputx=[]
            for col in range(64):
                inputx.append([self.ix[(row*64+col)*2],self.ix[(row*64+col)*2+1]])
                outputx.append([self.x[(row*64+col)*2],self.x[(row*64+col)*2+1]])
            tmpx1=operation.roate_right(inputx,self.rotate_right[row][1])
            tmpx2=operation.roate_right(inputx,self.rotate_right[row][2])
            for col in range(64):
                #self.m.update()
                #print(inputx[col]+tmpx1[col]+tmpx2[col]+outputx[col])
                sbox_model.gen_model_from_ine(self.m,inputx[col]+tmpx1[col]+tmpx2[col]+outputx[col],
                        'ine/3xor_ine.txt')

        if self.flag==0:
            for r in range(self.r):
                indexStart=640*r
                xstate=[]
                sout_state=[]
                xorout_state=[]
                for row in range(5):
                    tmpx=[]
                    tmpxor=[]
                    tmps=[]
                    for col in range(64):
                        tmpx.append([self.x[(row*64+col)*2+indexStart],self.x[(row*64+col)*2+indexStart+1]])
                        tmps.append([self.sout[(row*64+col)*2+indexStart],self.sout[(row*64+col)*2+indexStart+1]])
                        tmpxor.append([self.x[(row*64+col)*2+indexStart+640],self.x[(row*64+col)*2+indexStart+641]])
                    xstate.append(tmpx)
                    xorout_state.append(tmpxor)
                    sout_state.append(tmps)

                #非线性层
                for col in range(64):
                    tmpx=[]
                    tmps=[]
                    for row in range(5):
                        tmpx+=xstate[row][col]
                        tmps.append(sout_state[row][col])
                    for i in range(5):
                        tmp=tmpx+tmps[i]
                        sbox_model.gen_model_from_ine(self.m,tmp,'ascon/ascon_sbox_'+str(i)+'.ine')
                
                #线性层
                for row in range(5):
                    tmp1=operation.roate_right(sout_state[row],self.rotate_right[row][1])
                    tmp2=operation.roate_right(sout_state[row],self.rotate_right[row][2])
                    for col in range(64):
                        sbox_model.gen_model_from_ine(self.m,sout_state[row][col]+tmp1[col]+tmp2[col]+xorout_state[row][col],
                        'ine/3xor_ine.txt')
                
            index=640*self.r+self.index*2
            self.m.addConstr(self.x[index]==self.value[0],name='output0')
            
            self.m.addConstr(self.x[index+1]==self.value[1],name='output1')
            
        if self.flag==1:
            for r in range(self.r-1):
                indexStart=640*r
                xstate=[]
                sout_state=[]
                xorout_state=[]
                for row in range(5):
                    tmpx=[]
                    tmpxor=[]
                    tmps=[]
                    for col in range(64):
                        tmpx.append([self.x[(row*64+col)*2+indexStart],self.x[(row*64+col)*2+indexStart+1]])
                        tmps.append([self.sout[(row*64+col)*2+indexStart],self.sout[(row*64+col)*2+indexStart+1]])
                        tmpxor.append([self.x[(row*64+col)*2+indexStart+640],self.x[(row*64+col)*2+indexStart+641]])
                    xstate.append(tmpx)
                    xorout_state.append(tmpxor)
                    sout_state.append(tmps)

                #非线性层
                for col in range(64):
                    tmpx=[]
                    tmps=[]
                    for row in range(5):
                        tmpx+=xstate[row][col]
                        tmps.append(sout_state[row][col])
                    for i in range(5):
                        tmp=tmpx+tmps[i]
                        sbox_model.gen_model_from_ine(self.m,tmp,'ascon/ascon_sbox_'+str(i)+'.ine')
                
                #线性层
                for row in range(5):
                    tmp1=operation.roate_right(sout_state[row],self.rotate_right[row][1])
                    tmp2=operation.roate_right(sout_state[row],self.rotate_right[row][2])
                    for col in range(64):
                        sbox_model.gen_model_from_ine(self.m,sout_state[row][col]+tmp1[col]+tmp2[col]+xorout_state[row][col],
                        'ine/3xor_ine.txt')
            
            #额外的非线性层
            xstate=[]
            sout_state=[]
            for row in range(5):
                tmpx=[]
                tmps=[]
                for col in range(64):
                    tmpx.append([self.x[(self.r-1)*640+(row*64+col)*2],self.x[(self.r-1)*640+(row*64+col)*2+1]])
                    tmps.append([self.sout[(self.r-1)*640+(row*64+col)*2],self.sout[(self.r-1)*640+(row*64+col)*2+1]])
                xstate.append(tmpx)
                sout_state.append(tmps)

            for col in range(64):
                tmpx=[]
                tmps=[]
                for row in range(5):
                    tmpx+=xstate[row][col]
                    tmps.append(sout_state[row][col])
                for i in range(5):
                    tmp=tmpx+tmps[i]
                    #sbox_model.gen_model_from_ine(self.m,tmp,'ascon/ascon_ine_'+str(i)+'.txt')
                    sbox_model.gen_model_from_ine(self.m,tmp,'ascon/ascon_sbox_'+str(i)+'.ine')
            index=640*(self.r-1)+self.index*2
            self.m.addConstr(self.sout[index]==self.value[0],name='output0')
            self.m.addConstr(self.sout[index+1]==self.value[1],name='output1')


if __name__=='__main__':
    for i in range(5):
        t1=Ascon(3,320 , i*64, 0)
        t1.gen_constr()
        t1.m.optimize()
        _get_value(t1)
    
