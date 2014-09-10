import Tkinter as tk, tkFileDialog, tkMessageBox
    
root = tk.Tk()

directory = tk.StringVar()
x = tk.StringVar()
y = tk.StringVar()
z = tk.StringVar()
skip = tk.StringVar()
proc = tk.StringVar()
step = tk.StringVar()
plot = tk.BooleanVar()
test = tk.BooleanVar()

tk.Entry(root,textvariable=directory).grid(row=0,column=1)
tk.Entry(root,textvariable=x).grid(row=1,column=1)
tk.Entry(root,textvariable=y).grid(row=2,column=1)
tk.Entry(root,textvariable=z).grid(row=3,column=1)
tk.Entry(root,textvariable=skip).grid(row=4,column=1)
tk.Entry(root,textvariable=proc).grid(row=5,column=1)
tk.Entry(root,textvariable=step).grid(row=6,column=1)

try:
    params = open('config.div','r').readline().split()    
    directory.set(params[0])
    x.set(params[1])
    y.set(params[2])
    z.set(params[3])
    skip.set(params[4])
    proc.set(params[5])
    step.set(params[6])
    plot.set((params[7] == 'True'))
    test.set((params[8] == 'True'))
    print 'El archivin existe :D'
except:
    print 'El archivin no existe D:'

tk.Label(root, text='Path').grid(row=0)
tk.Label(root, text='X Column').grid(row=1)
tk.Label(root, text='Y Column').grid(row=2)
tk.Label(root, text='Z Column').grid(row=3)
tk.Label(root, text='Lines to Skip').grid(row=4)
tk.Label(root, text='Proceses').grid(row=5)
tk.Label(root, text='MCMC Steps').grid(row=6)
tk.Label(root, text='Plotting Mode').grid(row=7)
tk.Label(root, text='Test Mode').grid(row=8)

def change_path():
    directory.set(tkFileDialog.askdirectory())
    
def save():
    export = open('config.div','w')
    line =  directory.get()+' '+x.get()+' '+y.get()+' '+z.get()+' '+skip.get()+' '+proc.get()+' '+step.get()+' '+('%r' % (plot.get()))+' '+('%r' % (test.get()))+'\n'
    export.write(line)
    export.close()
    tkMessageBox.showinfo('Save', 'Configuration file was saved successfully')


tk.Button(root,text='Save',command=save).grid(row=9)
tk.Button(root,text='Change Path',command=change_path).grid(row=9,column=1)

tk.Checkbutton(root,text='Yes',variable=plot).grid(row=7, column=1)
tk.Checkbutton(root,text='Yes',variable=test).grid(row=8,column=1)
    
root.mainloop()

#root.withdraw()


