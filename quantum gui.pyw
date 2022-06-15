from quantum import *
import pandas as pd
import tkinter as tk
import tkinter.messagebox as tk_messagebox
import tkinter.filedialog as filedialog
import tkinter.simpledialog as simpledialog
import tkinter.ttk as ttk
import inspect
from functools import partial

about_str ="""Quantum Circuit Simulator


GUI Version: 1.0
Simulator Version: """+quantum_simulator_version+"""


Created by Aditya Gaikwad
"""
set_resizeable(True)

components = [H,X,Y,Z,S,T,R,CX,CY,CZ,CR,SWAP,SQRTSWAP,CCX,QFT,Measure]
inbuilt_circuit =      [BasicQCircuit.SET0,BasicQCircuit.SET1,BasicQCircuit.IDENTITY,BasicQCircuit.NOT,BasicQCircuit.OR,BasicQCircuit.AND,BasicQCircuit.XOR,BasicQCircuit.NOR,BasicQCircuit.NAND,BasicQCircuit.XNOR,BasicQCircuit.OR_3,BasicQCircuit.AND_3,BasicQCircuit.HALF_ADDER,BasicQCircuit.FULL_ADDER,BasicQCircuit.QUANTUM_TELEPORTATION]
root = tk.Tk(className='QCS')
root.title("Quantum Circuit Simulator")

class GUIQCircuit:
    def __init__(self,qbits):
        self.qbits=qbits
        self.qcircuit=QCircuit(qbits)
        self.qbit=[self.qcircuit.getQBit(i) for i in range(qbits)]
        self.component = self.qcircuit.gate
        self.canvas=None
        self.scrollx=0
        self.max_srollx=0
        self.scrollbar_x=None
        self.scrolly=0
        self.max_srolly=0
        self.scrollbar_y=None
        self.sw = 0
        self.sh = 0
        self.scroll_update_X = False
        self.scroll_update_Y = False
        self.component_size=[]
        self.para = []
        self.para_name = []

    def add(self,component,index):
        self.qcircuit.add(component,index)
        self.render()

    def pop(self,index):
        self.qcircuit.remove(self.qcircuit.gate[index])
        self.render()

    def set_canvas(self,canvas):
        self.canvas=canvas

    def draw_component(self,index,line,y):
        #print(type(component))
        component = circuit.component[index]
        para = circuit.para[index]
        para_name = circuit.para_name[index]
        canvas =self.canvas
        qbits = self.qbit
        qi = [q.index for q in qbits]
        if type(component) is not QCircuit:
            text = (component.__repr__()).split('\t-\t')
            #print(text)
            gp = sum([1 for i in text if i.strip()!='-'])
            if gp>1:
                fqi = 0
                for i in range(len(text)):
                    if text[i].strip()!='-':
                        fqi = i
                        break
                lqi = 0
                for i in range(len(text)-1,-1,-1):
                    if text[i].strip()!='-':
                        fqi = i
                        break
                canvas.create_rectangle(line[fqi]-50,y,line[lqi]+50,y+100,fill='yellow')
                canvas.create_line(line[fqi],y+50,line[lqi],y+50)
            for i in range(len(text)):
                c=text[i].strip()
                if c!='-':
                    canvas.create_rectangle(line[i]-50,y,line[i]+50,y+100,fill='yellow')
                    size = np.ceil(75/len(c))
                    canvas.create_text(line[i],y+50,text=c,font=("Purisa",int(size)))
        else:
            name = component.name.replace('_',' ')
            cq = component.qbit
            fqi = 0
            for i in range(len(qbits)):
                if qbits[i] in cq:
                    fqi=i
                    break
            lqi = 0
            for i in range(len(qbits)-1,-1,-1):
                if qbits[i] in cq:
                    lqi=i
                    break
            mp = (line[fqi]+line[lqi])//2
            fs = int(50/np.ceil(len(name)/((lqi-fqi+1)*3)))
            canvas.create_rectangle(line[fqi]-50,y,line[lqi]+50,y+100,fill='red')
            canvas.create_text(mp,y+50,text=name,font=("Purisa",fs))
            for i in range(len(para)):
                if type(para[i]) is QBit:
                    qi = para[i].index
                    canvas.create_text(line[qi],y+10,text=para_name[i].replace('_',' '))
                    canvas.create_text(line[qi],y+90,text=para_name[i].replace('_',' '))
        self.component_size.append([y,y+100])
        

    def render(self):
        cw=int(self.canvas['width'])
        ch=int(self.canvas['height'])
        #print(cw,ch)
        self.max_scrollx=cw
        self.max_scrolly=ch
        self.sw = cw
        self.sh = ch
        scrollx=self.scrollx
        scrolly=self.scrolly
        self.canvas.create_rectangle(0,0,cw+10,ch+10,fill='white')
        #print(self.canvas.bbox("all"))
        self.canvas.config(scrollregion=(0,0,cw,ch))
        line=[((i+1)/(4))*cw-scrollx for i in range(self.qbits)]
        self.component_size=[]
        for i in range(self.qbits):
            self.canvas.create_line(line[i],20-scrolly,line[i],ch)
        cy = 50-scrolly
        for i in range(len(self.component)):
            self.canvas.create_line(0,cy-10,cw+10,cy-10,fill='grey')
            component=self.component[i]
            self.draw_component(i,line,cy)
            self.canvas.create_line(0,cy+110,cw+10,cy+110,fill='grey')
            cy += 150
        self.canvas.create_rectangle(0,0,cw+10,20,fill='pink')
        for i in range(self.qbits):
            name=self.qbit[i].getName()
            self.canvas.create_text(line[i],10,text=name)
        self.max_scrolly = cy+ch
        self.max_scrollx = (cw*(self.qbits+1))//4
        self.scroll_update_X = True
        self.scroll_update_Y = True
        mx = (self.scrollx+self.sw)/self.max_scrollx
        if mx>1:
            mx=1
        my = (self.scrolly+self.sh)/self.max_scrolly
        if my>1:
            my=1
        self.scrollbar_x.set(self.scrollx/self.max_scrollx,mx)
        self.scrollbar_y.set(self.scrolly/self.max_scrolly,my)

    def config(self,x):
        #print('State:',x.state)
        #print('Type:',x.type)
        #print('X:',x.x)
        #print('Y:',x.y)
        #print('Width:',x.width)
        #print('Height:',x.height)
        #print('Delta:',x.delta)
        #print(x.width,x.height)
        cw=int(self.canvas['width'])
        ch=int(self.canvas['height'])
        if np.absolute(x.width-cw)>10 or np.absolute(x.height-ch)>10:
            self.canvas.config(width=x.width,height=x.height)
        self.render()
            
    def scroll_x(self,x,y,z=None):
        #print(x,y,z)
        #print('X',self.scroll_update_X)
        if self.scroll_update_X:
            self.scroll_update_X=False
            return
        if x=='scroll':
            y2=int(y)
            if (y2<0 and self.scrollx>0) or (y2>0 and self.scrollx<self.max_scrollx):
                self.scrollx += 10*y2
        elif x=='moveto':
            self.scrollx = float(y)*self.max_scrollx
        self.render()
        #return "moveto",self.scrollx/self.max_scrollx
    
    def scroll_y(self,x,y,z=None):
        #print(x,y,z)
        #print('Y',self.scroll_update_Y)
        if self.scroll_update_Y:
            self.scroll_update_Y = False
            return
        if x=='scroll':
            y2=int(y)
            if (y2<0 and self.scrolly>0) or (y2>0 and self.scrolly<self.max_scrolly):
                self.scrolly += 10*y2
        elif x=='moveto':
            self.scrolly = float(y)*self.max_scrolly
        self.render()
        #return "moveto",self.scrolly/self.max_scrolly
        

circuit = GUIQCircuit(3)


tframe=tk.Frame(root)
tframe.pack(side="top",fill="both",expand=True)
canvas = tk.Canvas(tframe,bg="white",height=500,width=500)
canvas.pack(side="left",fill="both",expand=True)
circuit.set_canvas(canvas)
vscrollbar=tk.Scrollbar(tframe)
vscrollbar.pack(side="right",fill="y",expand=True)
vscrollbar.config(command=circuit.scroll_y,orient='vertical')
circuit.scrollbar_y = vscrollbar
hscrollbar=tk.Scrollbar(root)
hscrollbar.pack(side="bottom",fill="x",expand=True)
hscrollbar.config(command=circuit.scroll_x,orient='horizontal')
circuit.scrollbar_x = hscrollbar
canvas.bind('<Configure>',circuit.config)

def choose_qbit(entry,ind,parent):
    top = tk.Toplevel(parent)
    top.title('Select QBit')
    frame = tk.Frame(top)
    l = tk.Label(frame,text='Select QBit:')
    l.pack()
    qbits = circuit.qbit
    qv = tk.IntVar(value=entry[ind])
    nr = tk.Radiobutton(frame,text='None',variable=qv,value=-1)
    nr.pack(anchor='w')
    for i in range(len(qbits)):
        qbit = qbits[i]
        R = tk.Radiobutton(frame,text=qbit.name,variable=qv,value=i)
        R.pack(anchor='w')
    def update_choice():
        i = qv.get()
        #print(ind,i)
        entry[ind]=i
        top.destroy()
    b = tk.Button(frame,text='Choose',command=update_choice)
    b.pack(anchor='s')
    frame.pack(anchor='w')

def add_component_with_index(ci,index,allow_destroy=False):
    c=None
    args = None
    args_args = None
    if ci<len(components):
        c=components[ci]
        #print(c.__name__,index)
        args = inspect.getargspec(c.__init__)
        args_args = args.args
    elif ci<len(inbuilt_circuit)+len(components):
        c=inbuilt_circuit[ci-len(components)]
        args = inspect.getargspec(c)
        args_args = args.args
        args_args.insert(0,"")
    #print(args.args)
    top = tk.Toplevel(root)
    top.title("Parameters")
    frame = tk.Frame(top)
    l = tk.Label(frame,text='Parameters for Gate:')
    l.pack()
    entries = []
    qbit_entries = [-1 for i in range(1,len(args.args))]
    for i in range(1,len(args_args)):
        a = args_args[i]
        f = tk.Frame(frame)
        l = tk.Label(f,text=a.replace('_',' ')+':')
        l.pack(side='left')
        b = tk.Button(f,text='Select QBit',command=partial(choose_qbit,qbit_entries,i-1,top))
        b.pack(side='right')
        te = tk.Entry(f)
        te.pack(side='right')
        entries.append(te)
        f.pack(anchor='w')
    if allow_destroy:
        para = circuit.para[index].copy()
        for i in range(len(para)):
            p = para[i]
            if type(p) is QBit:
                qbit_entries[i] = p.index
            else:
                entries[i].insert(0,p)
        def destroy_component():
            confirm = tk_messagebox.askyesno('Delete','Are you sure you want to delete the component')
            if confirm:
                circuit.para.pop(index)
                circuit.para_name.pop(index)
                circuit.pop(index)
                top.destroy()
        db = tk.Button(frame,text='Delete',command=destroy_component)
        db.pack(anchor='s')
    def create_component():
        #print(entries)
        #print(qbit_entries)
        qbits = circuit.qbit
        para = [None for i in range(1,len(args_args))]
        for i in range(len(qbit_entries)):
            if qbit_entries[i]!=-1:
                #print(qbit_entries[i])
                para[i] = qbits[qbit_entries[i]]
            else:
                para[i] = entries[i].get()
        #print(para)
        #f=True
        try:
            if allow_destroy:
                circuit.pop(index)
                #circuit.para.pop(index)
                circuit.para_name.pop(index)
            comp = c(*para)
            circuit.para.insert(index,para)
            circuit.para_name.insert(index,args_args[1:])
            circuit.add(comp,index)
            #print(comp)
            top.destroy()
        except:
            tk_messagebox.showerror('Parameter Error',"Invalid Parameter")
            
    b = tk.Button(frame,text='Done',command=create_component)
    b.pack(anchor='s')
    frame.pack(anchor='w')
    


def show_component_choose(index):
    top = tk.Toplevel(root)
    top.title('Add new component')
    frame = tk.Frame(top)
    label = tk.Label(top)
    label['text']="Add circuit component:"
    label.pack()
    fg = tk.Frame(frame)
    fg.pack(side='left')
    ci = tk.IntVar()
    for i in range(len(components)):
        #print(components[i].__name__)
        R = tk.Radiobutton(fg,text=components[i].__name__,variable=ci,value=i)
        R.pack(anchor='w')
    fic = tk.Frame(frame)
    fic.pack(side='right',anchor='n')
    for i in range(len(inbuilt_circuit)):
        R = tk.Radiobutton(fic,text=inbuilt_circuit[i].__name__.replace('_',' '),variable=ci,value=i+len(components))
        R.pack(anchor='w')
    def component_selected():
        i=ci.get()
        if i<len(components):
            add_component_with_index(i,index)
        elif i<len(components)+len(inbuilt_circuit):
            add_component_with_index(i,index)
        top.destroy()
    b = tk.Button(top,text="Next",command=component_selected)
    b.pack(anchor='s',side='bottom')
    frame.pack(anchor='w')

def click(event):
    widget=event.widget
    x=widget.canvasx(event.x)
    y=widget.canvasx(event.y)
    #print(x,y)
    add_new = True
    index = 0
    for i in range(len(circuit.component_size)):
        if circuit.component_size[i][0]<y<circuit.component_size[i][1]:
            #print(i)
            add_new = False
            index = i
            break
        elif i>0:
            if circuit.component_size[i-1][1]<y<circuit.component_size[i][0]:
                index = i
                break
    if len(circuit.component_size)>0 and circuit.component_size[-1][1]<y:
        index = len(circuit.component_size)

    if add_new:
        #print("New",index)
        show_component_choose(index)
    else:
        #print("Edit",index)
        comp = circuit.component[index]
        if type(comp) is not QCircuit:
            for i in range(len(components)):
                if type(comp) is components[i]:
                    add_component_with_index(i,index,True)
                    break
        else:
            cn = comp.name.upper().replace(' ','_')
            for i in range(len(inbuilt_circuit)):
                if cn == inbuilt_circuit[i].__name__:
                    add_component_with_index(i+len(components),index,True)
                    break
                

def nothing():
    print("Nothing")

def about():
    tk_messagebox.showinfo(title="About",message=about_str)

def circuit_settings():
    print("Circuit Settings")
    #top = tk.Toplevel(root)
    #frame = tk.Frame(top)
    #frame.pack()

def change_qbit_count():
    global circuit
    reset = tk_messagebox.askyesno('Change qbits',"Change number of qbits in the circuit?\nThis will reset the circuit.")
    if reset:
        n = simpledialog.askinteger('Number of qbits',"Enter the new number of qbits:",parent=root,minvalue=1)
        #print(n)
        if n!=None or n!=circuit.qbits:
            circuit = GUIQCircuit(n)
            circuit.set_canvas(canvas)
            vscrollbar.config(command=circuit.scroll_y)
            circuit.scrollbar_y = vscrollbar
            hscrollbar.config(command=circuit.scroll_x)
            circuit.scrollbar_x = hscrollbar
            canvas.bind('<Configure>',circuit.config)
            circuit.render()
            circuit.scroll_x(0,0)
            circuit.scroll_y(0,0)
            

def cir_sim():
    top = tk.Toplevel(root)
    top.title('Measurements')
    frame = tk.Frame(top)
    data = []
    tree = ttk.Treeview(frame)
    tree.pack(side='left',expand=True,fill='both')
    vsb = tk.Scrollbar(top,orient='vertical')
    vsb.pack(side='right',fill='y')
    tree['yscrollcommand']=vsb.set
    vsb['command']=tree.yview
    hsb = tk.Scrollbar(top,orient='horizontal')
    hsb.pack(side='bottom',fill='x')
    tree['xscrollcommand']=hsb.set
    hsb['command']=tree.xview
    columns = []
    def run_cir():
        #print("Run")
        circuit.qcircuit.run()
        data.append(circuit.qcircuit.measurement.copy())
        for k in data[-1].keys():
            if k not in columns:
                columns.append(k)
        tree['columns']=columns
        for c in columns:
            tree.heading(c,text=c)
        item = []
        for c in columns:
            if c in data[-1].keys():
                item.append(data[-1][c])
            else:
                item.append("")
        tree.insert("",'end',text=str(len(data)-1),values=item)
        #print(data)
    fn = [""]
    def save_data():
        if fn[0] == "":
            saveas_data()
            return
        #print(fn[0])
        save_data_with_name(fn[0])
    def saveas_data():
        prev = fn[0]
        fn[0]=filedialog.asksaveasfilename()
        if fn[0]=="":
            fn[0]=prev
            #print("Cancel")
            return
        #print(fn[0])
        top.title('Measurements -'+fn[0])
        save_data_with_name(fn[0])
    def save_data_with_name(f):
        df = pd.DataFrame()
        for c in columns:
            cd = []
            for i in range(len(data)):
                if c in data[i].keys():
                    cd.append(data[i][c])
                else:
                    cd.append(None)
            cd = pd.Series(cd.copy())
            df[c]=cd
        df.to_csv(f,index=False)
    menubar = tk.Menu(top)
    file_menu = tk.Menu(menubar,tearoff=0)
    #file_menu.add_command(label="New",command=nothing)
    #file_menu.add_command(label="Open",command=nothing)
    file_menu.add_command(label="Save",command=save_data)
    file_menu.add_command(label="Save As",command=saveas_data)
    file_menu.add_separator()
    file_menu.add_command(label="Close",command=top.destroy)
    menubar.add_cascade(label="File",menu=file_menu)
    menubar.add_command(label="Run",command=run_cir)
    top.config(menu=menubar)
    frame.pack(side='left',expand=True,fill='both')

canvas.bind("<Button-1>",click)

menubar=tk.Menu(root)
file_menu=tk.Menu(menubar,tearoff=0)
file_menu.add_command(label="New",command=nothing)
file_menu.add_command(label="Open",command=nothing)
file_menu.add_command(label="Save",command=nothing)
file_menu.add_command(label="Save As...",command=nothing)
file_menu.add_separator()
file_menu.add_command(label="Quit",command=root.destroy)
menubar.add_cascade(label="File",menu=file_menu)
root.config(menu=menubar)
options_menu=tk.Menu(menubar,tearoff=0)
options_menu.add_command(label="Circuit Settings",command=circuit_settings)
options_menu.add_command(label="!Change QBit count",command=change_qbit_count)
menubar.add_cascade(label="Options",menu=options_menu)
window_menu=tk.Menu(menubar,tearoff=0)
window_menu.add_command(label="Measurements",command=cir_sim)
menubar.add_cascade(label="Window",menu=window_menu)
help_menu=tk.Menu(menubar,tearoff=0)
help_menu.add_command(label="About",command=about)
menubar.add_cascade(label="Help",menu=help_menu)
circuit.scroll_x(0,0)
circuit.scroll_y(0,0)

circuit.render()
root.config(menu=menubar)
root.mainloop()
