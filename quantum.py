import numpy as np

quantum_simulator_version="1.0"
resizeable = False

def set_resizeable(value):
    global resizeable
    resizeable = value

def mod(a,b):
    c=a/b
    c=(c.real-np.floor(c.real)) +1j* (c.imag-np.floor(c.imag))
    return c*b

def gcd(a,b,n=0):
    if (n>100):
        return 1
    if np.abs(b.real)<0.000001:
        b=1j*b.imag
    if np.abs(b.imag)<0.000001:
        b=b.real
    if np.abs(a.real)<0.000001:
        a=1j*a.imag
    if np.abs(a.imag)<0.000001:
        a=a.real

    if np.abs(b.real)<0.000001:
        return a
    return gcd(b,mod(a,b),n+1)

def lcm(a,b):
    return a*b/gcd(a,b)

def showTable(records):
    keys=records[0].keys()
    srecords=[]
    ml=0
    for i in records:
        d={}
        for k,v in i.items():
            s=str(v).replace('\n','')
            
            d[k]=s
            if len(s)>ml:
                ml=len(s)
        srecords.append(d)
    for i in keys:
        if len(i)>ml:
            ml=len(i)
    for i in keys:
        print(i,end='\t'*(1+(ml-len(i))//8))
    print()
    for i in srecords:
        for j in keys:
            print(i[j],end='\t'*(1+(ml-len(i[j]))//8))
        print()

class QMatrix:
    def __init__(self,m=-1,n=-1,s=None):
        if type(s) is type(None):
            self.mat=np.array([[0 for j in range(n)]for i in range(m)])
        else:
            #print(s)
            self.mat=np.array(s)
            m=len(self.mat)
            n=len(self.mat[0])
        self.m=m
        self.n=n

    def get(self,i,j):
        return self.mat[i][j]

    def set(self,i,j,v):
        print(self.mat[i][j])
        self.mat[i][j]=v
        print(self.mat[i][j],v)

    def clone(self):
        return QMatrix(s=self.mat.copy())

    def __mul__(self,other): #Matrix Product
        if type(other) is QMatrix:
            return QMatrix(s=np.matrix(self.mat)*np.matrix(other.mat))
        return QMatrix(s=self.mat*other)

    def __pow__(self,other): #Tensor Product
        s=[[0 for j in range(self.n*other.n)] for i in range(self.m*other.m)]
        for i in range(self.m):
            for k in range(other.m):
                for j in range(self.n):
                    for l in range(other.n):
                        s[i*other.m+k][j*other.n+l]=self.get(i,j)*other.get(k,l)
        return QMatrix(s=s)

    def isValid(self):
        l=np.log2(self.m)
        if l<1 or l!=l//1:
            return False
        for j in range(self.n):
            p=-1
            for i in range(self.m):
                p+=np.abs(self.get(i,j))**2
            if p>0.01 or p<-0.01:
                return False
        return True

    def factorize(self):
        f0=self.get(0,0)
        for i in range(1,self.m//2):
            f0=gcd(f0,self.get(i,0))
        f1=self.get(self.m//2,0)
        for i in range(self.m//2,self.m):
            f1=gcd(f1,self.get(i,0))
        c=(f0**2 + f1**2)**0.5
        if c==0:
            return None
        f0/=c
        f1/=c
        m1=QMatrix(s=[[f0],[f1]])
        s=[[0] for i in range(self.m//2)]
        for i in range(self.m//2):
            i1=self.get(i,0)
            i2=self.get(i+self.m//2,0)
            if f0!=0 and f1!=0 and i1/f0 != i2/f1:
                return None
            if f0!=0:
                s[i][0]=i1/f0
            else:
                s[i][0]=i2/f1
            if np.abs(s[i][0].real)<0.000001:
                s[i][0]=1j*s[i][0].imag
            if np.abs(s[i][0].imag)<0.000001:
                s[i][0]=s[i][0].real
        m2=QMatrix(s=s)
        if m1.isValid() and m2.isValid():
            return [m1,m2]
        return None
        
    

    def __str__(self):
        s=''
        for i in range(self.m):
            s+='|'
            for j in range(self.n):
                if j!=0:
                    s+='\t'
                s+=str(self.mat[i][j])
            s+='|'
            if i!=self.m-1:
                s+='\n'
        return s
    
                
class QState:
    def __init__(self,q,parent,s=[[1],[0]]):
        self.index=[q.index]
        self.mat=QMatrix(s=s)
        self.parent=parent
        
    def entangle(self,q):
        #print('e')
        if q.state!=self:
                q.entangled=True
                il=len(self.index)
                index2 = []
                for i in range(il):
                    index2.append(self.index[i])
                for i in range(len(q.state.index)):
                    index2.append(q.state.index[i])
                self.mat=self.mat**q.state.mat
                f=False
                for i in range(len(q.state.index)):
                    #print(i,len(q.state.index),q.index)
                    if q.index!=q.state.index[i]:
                       q.parent.getQBit(q.state.index[i]).state=self
                    else:
                        f=True
                if f:
                    q.parent.getQBit(q.index).state=self
                self.index=index2

    def swap(self,i,j):
        t=self.index[i]
        self.index[i]=self.index[j]
        self.index[j]=t
        #
        m=[[0] for i in range(self.mat.m)]
        sl=len(self.index)-1
        for a in range(self.mat.m):
            bi=(a>>(sl-i))%2
            bj=(a>>(sl-j))%2
            b=a&~( (1<<(sl-i)) | (1<<(sl-j)) )
            c=b|(bi<<(sl-j))|(bj<<(sl-i))
            m[c][0]=self.mat.get(a,0)
        self.mat=QMatrix(s=m)

    def factorize(self):
        il=len(self.index)
        i=0
        while i<il:
            self.swap(0,i)
            s=self.mat.factorize()
            if s!=None:
                    q=self.parent.getQBit(self.index[0])
                    qs=QState(q,self.parent,s=s[0].mat)
                    q.state=qs
                    q.entangled=False
                    ind=[self.index[j] for j in range(1,len(self.index))]
                    self.index=ind
                    self.mat=s[1]
                    il-=1
                    i-=1
                    if il<=1 and il>0:
                        il=0
                        q2=self.parent.getQBit(self.index[0])
                        q2.entangled=False
            i+=1
        for i in range(il-1,0,-1):
            self.swap(0,i)

    def measure(self):
        r=np.random.random()
        i=0
        while i<self.mat.m:
            p=np.abs(self.mat.get(i,0))**2
            if r<p:
                break
            r-=p
            i+=1
        il=len(self.index)
        for j in range(0,il):
            b=(i>>(il-1-j))%2
            q=self.parent.getQBit(self.index[j])
            if b==1:
                q.state=QState(q,self.parent,s=[[0],[1]])
                q.entangled=False
            else:
                q.state=QState(q,self.parent,s=[[1],[0]])
                q.entangled=False

    def __str__(self):
        return 'QState:\n'+str(self.mat)
        
               
class QBit:
    def __init__(self,index,parent):
        self.index=index
        self.name=str(index)
        self.state=QState(self,parent)
        self.entangled=False
        self.parent=parent

    def setName(self,name):
        self.name=name

    def getName(self):
        return self.name

    def isEntangled(self):
        return self.entangled

    def measure(self):
        self.state.measure()
        if self.state.mat.get(1,0)>0.5:
            return True
        return False
        


class QGate:
    def __init__(self,operands):
        self.operand=operand
        self.parent=operand[0].parent
        self.operation=QMatrix(s=[
            [1,0],
            [0,1]
            ])

    def run(self):
        e=self.operand[0].state
        for i in range(len(self.operand)):
            e.entangle(self.operand[i])
            for j in range(len(e.index)):
                if e.index[j]==self.operand[i].index:
                    e.swap(i,j)
                    break
        s=self.operation.clone()
        for i in range(len(e.index)-len(self.operand)):
            s=s**QMatrix(s=[[1,0],[0,1]])
        e.mat=s*e.mat
        e.factorize()
    

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                s+='?'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

    def __str__(self):
        return self.__repr__()

class Measure(QGate):
    def __init__(self,qbit,name):
        self.operand=[qbit]
        self.name=name
        self.parent=qbit.parent

    def run(self):
        if self.name!=None:
            self.parent.measurement[self.name]=self.operand[0].measure()

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                n = str(self.name)
                if self.name==None:
                    n='Measure'
                if not resizeable and len(n)>7:
                    n=n[:7]
                s+=n
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class Numeric(QGate):
    def __init__(self,qbits,name):
        self.operand=qbits
        self.name=name
        self.parent=qbits[0].parent

    def run(self):
        n=0
        for q in self.operand:
            n *= 2
            if q.measure():
                n += 1
        self.parent.measurement[self.name]=n

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                n = str(self.name)
                if not resizeable and len(n)>7:
                    n=n[:7]
                s+=n
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class ViewMatrix(QGate):
    def __init__(self,qbit,name):
        self.operand=[qbit]
        self.name=name
        self.parent=qbit.parent

    def run(self):
        self.parent.measurement[self.name]=self.operand[0].state.mat
        #print(self.name,":",self.operand[0].state.mat)

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                n = str(self.name)
                if not resizeable and len(n)>7:
                    n=n[:7]
                s+=n
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class ChangeMatrix(QGate):
    def __init__(self,qbit,mat):
        self.operand=[qbit]
        self.mat=mat
        self.parent=qbit.parent

    def run(self):
        self.operand[0].state.mat = self.mat

    def __str__(self):
        return str('QMatrix')

class QCircuit:
    def __init__(self,qbits,name=None):
        self.qbit=[]
        self.defaultstate=[]
        self.parent=None
        self.name=name
        if type(qbits) is int:
            for i in range(qbits):
                q=QBit(i,self)
                self.qbit.append(q)
                self.defaultstate.append(QMatrix(s=[[1],[0]]))
        else:
            for qbit in qbits:
                self.qbit.append(qbit)
            self.parent=self.qbit[0].parent
        self.gate=[]
        self.measurement={}
        
            
    def getQBit(self,ind):
        #print(ind,self.parent,len(self.qbit))
        #if self.parent!=None:
            #return self.parent.getQBit(ind)
        return self.qbit[ind]

    def setName(self,name):
        self.name=name

    def getName(self):
        return self.name

    def add(self,component,index=-1):
        component.parent=self
        if self.parent!=None:
            component.parent=self.parent
            #self.parent.add(component)
        if index == -1:
            self.gate.append(component)
        else:
            self.gate.insert(index,component)

    def remove(self,component):
        component.parent=None
        self.gate.remove(component)

    def run(self):
        if self.parent == None:
            self.measurement = {}
            for i in range(len(self.qbit)):
                #print(self.defaultstate[i])
                self.qbit[i].state.mat=self.defaultstate[i].clone()
        for i in range(len(self.gate)):
            #print(type(self.gate[i]))
            self.gate[i].run()
        if self.parent != None:
            self.parent.measurement.update(self.measurement)
        else:
            for q in self.qbit:
                q.measure()

    def showMeasurements(self):
        for key in self.measurement:
            extra=''
            if type(self.measurement[key]) is QMatrix:
                extra='\n'
            print(key,':',extra+str(self.measurement[key]))

    def record(self,n):
        m=[]
        for i in range(n):
            self.run()
            m.append(self.measurement.copy())
        return m
    
    def __repr__(self):
        s = ''
        qb=self.qbit
        if self.parent!=None:
            qb=self.parent.qbit
        else:
            for q in range(len(qb)):
                name=qb[q].getName()
                if not resizeable and len(name)>7:
                    name=name[:7]
                s+=name+'\t'
                if q!=len(qb)-1:
                    s+='\t'
            s+='\n'
            for q in range(len(qb)):
                s+='|\t'
                if q!=len(qb)-1:
                    s+='\t'
            s+='\n'
        for gate in self.gate:
            s+=gate.__repr__()
            if gate!=self.gate[-1]:
                s+='\n'
                for q in qb:
                    s+='|\t'
                    if q!=qb[-1]:
                        s+='\t'
                s+='\n'
        return s

    def __str__(self):
        return self.__repr__()


class X(QGate):
    def __init__(self,qbit):
        self.operand=[qbit]
        self.operation=QMatrix(s=[
            [0,1],
            [1,0]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                s+='X'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class Y(QGate):
    def __init__(self,qbit):
        self.operand=[qbit]
        self.operation=QMatrix(s=[
            [0 ,-1j],
            [1j, 0 ]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                s+='Y'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class Z(QGate):
    def __init__(self,qbit):
        self.operand=[qbit]
        self.operation=QMatrix(s=[
            [1, 0],
            [0,-1]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                s+='Z'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class S(QGate):
    def __init__(self,qbit):
        self.operand=[qbit]
        self.operation=QMatrix(s=[
            [1,0 ],
            [0,1j]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                s+='S'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class T(QGate):
    def __init__(self,qbit):
        self.operand=[qbit]
        self.operation=QMatrix(s=[
            [1,0],
            [0,np.e**(1j*np.pi/4)]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                s+='T'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class R(QGate):
    def __init__(self,qbit,phi):
        self.operand=[qbit]
        phi = float(phi)
        self.operation=QMatrix(s=[
            [1,0],
            [0,np.e**(2j*np.pi*phi)]
            ])
        self.phi=phi

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                p=str(self.phi)
                if not resizeable and len(p)>6:
                    p=p[:6]
                s+='R'+p
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class H(QGate):
    def __init__(self,qbit):
        self.operand=[qbit]
        d=1/np.sqrt(2)
        self.operation=QMatrix(s=[
            [d,d],
            [d,-d]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                s+='H'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class SQRTX(QGate):
    def __init__(self,qbit):
        self.operand=[qbit]
        d=(1/2)**1/2
        a=d+1j*d
        b=d-1j*d
        self.operation=QMatrix(s=[
            [a,b],
            [b,a]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                s+='X(1/2)'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class CX(QGate):
    def __init__(self,control,qbit):
        self.operand=[control,qbit]
        self.operation=QMatrix(s=[
            [1,0,0,0],
            [0,1,0,0],
            [0,0,0,1],
            [0,0,1,0]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i == self.operand[0]:
                s+='O'
            elif i == self.operand[1]:
                s+='X'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class CY(QGate):
    def __init__(self,control,qbit):
        self.operand=[control,qbit]
        self.operation=QMatrix(s=[
            [1,0,0 , 0],
            [0,1,0 , 0],
            [0,0,0 ,-1j],
            [0,0,1j, 0]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i == self.operand[0]:
                s+='O'
            elif i == self.operand[1]:
                s+='Y'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class CZ(QGate):
    def __init__(self,control,qbit):
        self.operand=[control,qbit]
        self.operation=QMatrix(s=[
            [1,0,0, 0],
            [0,1,0, 0],
            [0,0,1, 0],
            [0,0,0,-1]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i == self.operand[0]:
                s+='O'
            elif i == self.operand[1]:
                s+='Z'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class CR(QGate):
    def __init__(self,control,qbit,phi):
        #print(phi)
        self.operand=[control,qbit]
        phi = float(phi)
        self.operation=QMatrix(s=[
            [1,0,0,0],
            [0,1,0,0],
            [0,0,1,0],
            [0,0,0,np.e**(1j*phi)]
            ])
        self.phi=phi

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i == self.operand[0]:
                s+='O'
            elif i in self.operand:
                p=str(self.phi)
                if not resizeable and len(p)>6:
                    p=p[:6]
                s+='R'+p
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class SWAP(QGate):
    def __init__(self,qbit0,qbit1):
        self.operand=[qbit0,qbit1]
        self.operation=QMatrix(s=[
            [1,0,0,0],
            [0,0,1,0],
            [0,1,0,0],
            [0,0,0,1]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i == self.operand[0]:
                s+='X'
            elif i == self.operand[1]:
                s+='X'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class SQRTSWAP(QGate):
    def __init__(self,qbit0,qbit1):
        self.operand=[qbit0,qbit1]
        d=(1/2)**1/2
        a=d+1j*d
        b=d-1j*d
        self.operation=QMatrix(s=[
            [1,0,0,0],
            [0,a,b,0],
            [0,b,a,0],
            [0,0,0,1]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i == self.operand[0]:
                s+='X(1/2)'
            elif i == self.operand[1]:
                s+='X(1/2)'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class CCX(QGate):
    def __init__(self,control0,control1,qbit):
        self.operand=[control0,control1,qbit]
        self.operation=QMatrix(s=[
            [1,0,0,0,0,0,0,0],
            [0,1,0,0,0,0,0,0],
            [0,0,1,0,0,0,0,0],
            [0,0,0,1,0,0,0,0],
            [0,0,0,0,1,0,0,0],
            [0,0,0,0,0,1,0,0],
            [0,0,0,0,0,0,0,1],
            [0,0,0,0,0,0,1,0]
            ])

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i == self.operand[0] or i == self.operand[1]:
                s+='O'
            elif i == self.operand[2]:
                s+='X'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s

class QFT(QGate):
    def __init__(self,qbits):
        self.operand=qbits
        l=len(qbits)
        n=2**l
        f=n**0.5
        w=np.e**(2j*np.pi/n)
        #print(l,n,f,w)
        m=[[(w**-(i*j))/f for j in range(n)] for i in range(n)]
        self.operation=QMatrix(s=m)
        #print(self.operation)

    def run(self):
        e=self.operand[0].state
        for i in range(len(self.operand)):
            e.entangle(self.operand[i])
            for j in range(len(e.index)):
                if e.index[j]==self.operand[i].index:
                    e.swap(i,j)
                    break
        s=self.operation.clone()
        for i in range(len(e.index)-len(self.operand)):
            s=s**QMatrix(s=[[1,0],[0,1]])
        e.mat=s*e.mat
        e.factorize()

    def __repr__(self):
        s=''
        for i in self.parent.qbit:
            if i in self.operand:
                s+='QFT'
            else:
                s+='-'
            s+='\t'
            if i != self.parent.qbit[-1]:
                s+='-\t'
        return s


class BasicQCircuit:

    def SET0(qbit_in,qbit_out):
        c=QCircuit([qbit_in,qbit_out],"Set0")
        return c
    
    def SET1(qbit_in,qbit_out):
        c=QCircuit([qbit_in,qbit_out],"Set1")
        c.add(X(qbit_out))
        return c

    def IDENTITY(qbit_in,qbit_out):
        c=QCircuit([qbit_in,qbit_out],"Identity")
        c.add(CX(qbit_in,qbit_out))
        return c
    
    def NOT(qbit_in,qbit_out):
        c=QCircuit([qbit_in,qbit_out],"NOT")
        c.add(X(qbit_out))
        c.add(CX(qbit_in,qbit_out))
        return c
    
    def AND(qbit_in_0,qbit_in_1,qbit_out):
        c=QCircuit([qbit_in_0,qbit_in_1,qbit_out],"AND")
        c.add(CCX(qbit_in_0,qbit_in_1,qbit_out))
        return c

    def AND_3(qbit_in_0,qbit_in_1,qbit_in_2,qbit_extra,qbit_out):
        c=QCircuit([qbit_in_0,qbit_in_1,qbit_in_2,qbit_extra,qbit_out],"AND 3")
        c.add(BasicQCircuit.AND(qbit_in_0,qbit_in_1,qbit_extra))
        c.add(BasicQCircuit.AND(qbit_extra,qbit_in_2,qbit_out))
        return c
    

    def NAND(qbit_in_0,qbit_in_1,qbit_out):
        c=QCircuit([qbit_in_0,qbit_in_1,qbit_out],"NAND")
        c.add(CCX(qbit_in_0,qbit_in_1,qbit_out))
        c.add(X(qbit_out))
        return c

    def XOR(qbit_in_0,qbit_in_1,qbit_out):
        c=QCircuit([qbit_in_0,qbit_in_1,qbit_out],"XOR")
        c.add(CX(qbit_in_0,qbit_out))
        c.add(CX(qbit_in_1,qbit_out))
        return c

    def XOR_3(qbit_in_0,qbit_in_1,qbit_in_2,qbit_extra,qbit_out):
        c=QCircuit([qbit_in_0,qbit_in_1,qbit_in_2,qbit_extra,qbit_out],"XOR 3")
        c.add(BasicQCircuit.XOR(qbit_in_0,qbit_in_1,qbit_extra))
        c.add(BasicQCircuit.XOR(qbit_extra,qbit_in_2,qbit_out))
        return c

    def XNOR(qbit_in_0,qbit_in_1,qbit_out):
        c=QCircuit([qbit_in_0,qbit_in_1,qbit_out],"XNOR")
        c.add(CX(qbit_in_0,qbit_out))
        c.add(CX(qbit_in_1,qbit_out))
        c.add(X(qbit_out))
        return c

    def OR(qbit_in_0,qbit_in_1,qbit_out):
        c=QCircuit([qbit_in_0,qbit_in_1,qbit_out],"OR")
        c.add(CX(qbit_in_0,qbit_out))
        c.add(CX(qbit_in_1,qbit_out))
        c.add(CCX(qbit_in_0,qbit_in_1,qbit_out))
        return c

    def OR_3(qbit_in_0,qbit_in_1,qbit_in_2,qbit_extra,qbit_out):
        c=QCircuit([qbit_in_0,qbit_in_1,qbit_in_2,qbit_extra,qbit_out],"OR 3")
        c.add(BasicQCircuit.OR(qbit_in_0,qbit_in_1,qbit_extra))
        c.add(BasicQCircuit.OR(qbit_extra,qbit_in_2,qbit_out))
        return c

    def NOR(qbit_in_0,qbit_in_1,qbit_out):
        c=QCircuit([qbit_in_0,qbit_in_1,qbit_out],"NOR")
        c.add(CX(qbit_in_0,qbit_out))
        c.add(CX(qbit_in_1,qbit_out))
        c.add(CCX(qbit_in_0,qbit_in_1,qbit_out))
        c.add(X(qbit_out))
        return c

    def HALF_ADDER(A,B,Carry,Sum):
        c=QCircuit([A,B,Carry,Sum],"HALF ADDER")
        c.add(BasicQCircuit.XOR(A,B,Sum))
        c.add(BasicQCircuit.AND(A,B,Carry))
        return c

    def FULL_ADDER(A,B,C,qbit_extra_0,qbit_extra_1,qbit_extra_2,Carry,Sum):
        c=QCircuit([A,B,C,qbit_extra_0,qbit_extra_1,qbit_extra_2,Carry,Sum],"FULL ADDER")
        c.add(BasicQCircuit.XOR(A,B,qbit_extra_0))
        c.add(BasicQCircuit.XOR(C,qbit_extra_0,Sum))
        c.add(BasicQCircuit.AND(A,B,qbit_extra_1))
        c.add(BasicQCircuit.AND(C,qbit_extra_0,qbit_extra_2))
        c.add(BasicQCircuit.OR(qbit_extra_1,qbit_extra_2,Carry))
        return c

    def QUANTUM_TELEPORTATION(qbit_in,qbit_extra,qbit_out):
        c=QCircuit([qbit_in,qbit_extra,qbit_out],"QUANTUM TELEPORTATION")
        c.add(H(qbit_out))
        c.add(CX(qbit_out,qbit_extra))
        c.add(CX(qbit_in,qbit_extra))
        c.add(H(qbit_in))
        c.add(CX(qbit_extra,qbit_out))
        c.add(CZ(qbit_in,qbit_out))
        return c

    def DEUTSCH(qbit_in,qbit_out,f):
        if type(qbit_in) is QBit:
            c = QCircuit([qbit_in,qbit_out],"DEUTCH")
            c.add(X(qbit_out))
            c.add(H(qbit_in))
            c.add(H(qbit_out))
            c.add(f)
            c.add(H(qbit_in))
            return c
    def PHASE_ESTIMATOR(qbits_in,qbit_out,f,para=()):
        c=QCircuit([*qbits_in,qbit_out],"PHASE ESTIMATOR")
        c.add(X(qbit_out))
        for q in qbits_in:
            c.add(H(q))
        n=1
        for q in qbits_in[::-1]:
            for i in range(n):
                c.add(f(q,qbit_out,*para))
            n*=2
        c.add(QFT(qbits_in))
        return c
        

    def encode(qbits,num):
        n=int(num)
        c=QCircuit(qbits,"ENCODE")
        for i in range(len(qbits)):
            if n%2==1:
                c.add(X(qbits[-i-1]))
            n /= 2
        return c
        

def test():
    c=QCircuit(3)
    q=[c.getQBit(i) for i in range(3)]
    c.add(H(q[0]))
    c.add(H(q[1]))
    c.add(CCX(q[0],q[1],q[2]))
    c.add(Measure(q[0],'M0'))
    c.add(Measure(q[1],'M1'))
    c.add(Measure(q[2],'M2'))
    print(c)
    showTable(c.record(10))

def test1():
    c=QCircuit(4)
    q=[c.getQBit(i) for i in range(4)]
    q[0].setName("A")
    q[1].setName("B")
    q[2].setName("CARRY")
    q[3].setName("SUM")
    c.add(H(q[0]))
    c.add(H(q[1]))
    c.add(BasicQCircuit.HALF_ADDER(q[0],q[1],q[2],q[3]))
    c.add(Measure(q[0],'A'))
    c.add(Measure(q[1],'B'))
    c.add(Measure(q[2],'CARRY'))
    c.add(Measure(q[3],'SUM'))
    print(c)
    showTable(c.record(10))

def test2():
    c=QCircuit(8)
    q=[c.getQBit(i) for i in range(8)]
    q[0].setName("A")
    q[1].setName("B")
    q[2].setName("C")
    q[3].setName("CARRY")
    q[4].setName("SUM")
    q[5].setName("EXTRA 1")
    q[6].setName("EXTRA 2")
    q[7].setName("EXTRA 3")
    c.add(H(q[0]))
    c.add(H(q[1]))
    c.add(H(q[2]))
    c.add(BasicQCircuit.FULL_ADDER(q[0],q[1],q[2],q[5],q[6],q[7],q[3],q[4]))
    c.add(Measure(q[0],'A'))
    c.add(Measure(q[1],'B'))
    c.add(Measure(q[2],'C'))
    #c.add(BasicQCircuit.FULL_ADDER(q[0],q[1],q[2],q[5],q[6],q[7],q[3],q[4]))
    c.add(Measure(q[3],'CARRY'))
    c.add(Measure(q[4],'SUM'))
    print(c)
    showTable(c.record(10))

def test3():
    c=QCircuit(3)
    q=[c.getQBit(i) for i in range(3)]
    q[0].setName("IN")
    q[1].setName("EXTRA")
    q[2].setName("OUT")
    c.add(H(q[0]))
    c.add(Measure(q[0],"IN"))
    c.add(BasicQCircuit.QUANTUM_TELEPORTATION(q[0],q[1],q[2]))
    c.add(Measure(q[2],"OUT"))
    print(c)
    showTable(c.record(10))

def test4():
    c=QCircuit(3)
    q=[c.getQBit(i) for i in range(3)]
    q[0].setName("IN")
    q[1].setName("EXTRA")
    q[2].setName("OUT")
    f=BasicQCircuit.NOT(q[0],q[2])
    c.add(BasicQCircuit.DEUTSCH(q[0],q[2],f))
    c.add(Measure(q[0],"IN"))
    c.add(Measure(q[2],"OUT"))
    print(c)
    showTable(c.record(10))

def test5():
    c=QCircuit(3)
    q=[c.getQBit(i) for i in range(3)]
    q[0].setName("A")
    q[1].setName("B")
    q[2].setName("C")
    c.add(X(q[0]))
    c.add(X(q[1]))
    c.add(QFT(q[:-1]))
    #c.add(Z(q[0]))
    #c.add(Z(q[1]))
    c.add(QFT(q[:-1]))
    c.add(ViewMatrix(q[0],"MA"))
    c.add(ViewMatrix(q[1],"MB"))
    c.add(Measure(q[0],"A"))
    c.add(Measure(q[1],"B"))
    #c.add(Measure(q[2],"C"))
    c.add(Numeric([q[0],q[1]],'D'))
    print(c)
    #print(QFT(q[:-1]).operation*QFT(q[:-1]).operation)
    showTable(c.record(10))

def test6():
    tau = 2*np.pi
    c=QCircuit(4)
    q=[c.getQBit(i) for i in range(4)]
    q[0].setName("A")
    q[1].setName("B")
    q[2].setName("C")
    q[3].setName("D")
    f=CR
    c.add(BasicQCircuit.PHASE_ESTIMATOR(q[:-1],q[-1],f,(5*tau/8,)))
    c.add(Measure(q[0],"A"))
    c.add(Measure(q[1],"B"))
    c.add(Measure(q[2],"C"))
    c.add(Measure(q[3],"D"))
    c.add(Numeric([q[0],q[1],q[2]],'E'))
    print(c)
    #print(QFT(q[:-1]).operation*QFT(q[:-1]).operation)
    showTable(c.record(10))

def test7():
    c=QCircuit(6)
    q=[c.getQBit(i) for i in range(6)]
    q[0].setName("A")
    q[1].setName("B")
    q[2].setName("C")
    q[3].setName("D")
    q[4].setName("E")
    q[5].setName("F")
    f=CR
    c.add(BasicQCircuit.PHASE_ESTIMATOR(q[:-1],q[-1],f,(1,)))
    c.add(Measure(q[0],"A"))
    c.add(Measure(q[1],"B"))
    c.add(Measure(q[2],"C"))
    c.add(Measure(q[3],"D"))
    c.add(Measure(q[4],"E"))
    c.add(Measure(q[5],"F"))
    c.add(Numeric(q[:-1],'M'))
    print(c)
    #print(QFT(q[:-1]).operation*QFT(q[:-1]).operation)
    s=0
    n=100
    for i in range(n):
        c.run()
        #print(c.measurement['M'])
        s+=c.measurement['M']
    s /= n*32
    print(1/(2*s))

#test()
