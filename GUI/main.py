

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import os, subprocess, sys, json
from datetime import datetime
from load_file import LoadFileHandler
from save_file import SaveFileHandler
from new_file import NewFileHandler


class QuDSimGUI:
    APP_TITLE = "Welcome to QuDSim: 3D FEM Tool"
    VERSION = "1.0.0"
    CONFIG_FILE = "qudsim_config.json"
    MODULES = ["Schrodinger-Poisson","Tunneling Current","Dirac Equation","Poisson Only","Open Boundary"]
    BG_PRIMARY="#f0f4f8"; BG_SECONDARY="#ffffff"; ACCENT="#2e6da4"; ACCENT_DARK="#1a4d7e"
    TEXT_PRIMARY="#2c3e50"; TEXT_SECONDARY="#5d6d7e"; BORDER="#d5dfe8"

    def __init__(self, root):
        self.root = root
        self.root.title(self.APP_TITLE)
        self.root.geometry("920x720")
        self.root.minsize(800,600)
        self.root.configure(bg=self.BG_PRIMARY)
        self.project_name=tk.StringVar(value="")
        self.modules_required=tk.StringVar(value="Schrodinger-Poisson")
        self.project_version=tk.StringVar(value="1.0")
        self.maintainer_email=tk.StringVar(value="")
        self.current_geo_file=None; self.current_msh_file=None; self.current_vtu_file=None
        self.mpi_processes=tk.IntVar(value=4); self.input_param_file=None
        self.load_handler=LoadFileHandler(self); self.save_handler=SaveFileHandler(self); self.new_handler=NewFileHandler(self)
        self._create_menu_bar(); self._create_main_layout(); self._create_status_bar(); self._load_config()

    def _create_menu_bar(self):
        mb=tk.Menu(self.root)
        fm=tk.Menu(mb,tearoff=0)
        fm.add_command(label="New Project",command=self.new_handler.new_project,accelerator="Ctrl+N")
        fm.add_command(label="Open .geo File",command=lambda:self.load_handler.load_file("geo"),accelerator="Ctrl+O")
        fm.add_command(label="Open .msh File",command=lambda:self.load_handler.load_file("msh"))
        fm.add_command(label="Open .vtu File",command=lambda:self.load_handler.load_file("vtu"))
        fm.add_separator(); fm.add_command(label="Save",command=self.save_handler.save_current,accelerator="Ctrl+S")
        fm.add_command(label="Save As...",command=self.save_handler.save_as)
        fm.add_separator(); fm.add_command(label="Exit",command=self._on_close)
        mb.add_cascade(label="File",menu=fm)
        sm=tk.Menu(mb,tearoff=0)
        sm.add_command(label="Generate Mesh (GMSH)",command=self._run_gmsh)
        sm.add_command(label="Run Solver",command=self._run_solver)
        sm.add_separator(); sm.add_command(label="View Results (ParaView)",command=self._open_paraview)
        mb.add_cascade(label="Simulation",menu=sm)
        tm=tk.Menu(mb,tearoff=0)
        tm.add_command(label="Configure MPI",command=self._configure_mpi)
        tm.add_command(label="Input Parameters",command=self._edit_input_params)
        tm.add_command(label="Material Database",command=self._open_material_db)
        mb.add_cascade(label="Tools",menu=tm)
        hm=tk.Menu(mb,tearoff=0)
        hm.add_command(label="Documentation",command=self._show_docs)
        hm.add_command(label="About QuDSim",command=self._show_about)
        mb.add_cascade(label="Help",menu=hm)
        self.root.config(menu=mb)
        self.root.bind("<Control-n>",lambda e:self.new_handler.new_project())
        self.root.bind("<Control-o>",lambda e:self.load_handler.load_file("geo"))
        self.root.bind("<Control-s>",lambda e:self.save_handler.save_current())

    def _create_main_layout(self):
        outer=tk.Frame(self.root,bg=self.BG_PRIMARY); outer.pack(fill=tk.BOTH,expand=True,padx=20,pady=10)
        tf=tk.Frame(outer,bg=self.ACCENT,height=50); tf.pack(fill=tk.X,pady=(0,15)); tf.pack_propagate(False)
        tk.Label(tf,text=self.APP_TITLE,font=("Helvetica",16,"bold"),fg="white",bg=self.ACCENT).pack(expand=True)
        body=tk.Frame(outer,bg=self.BG_PRIMARY); body.pack(fill=tk.BOTH,expand=True)
        left=tk.LabelFrame(body,text="  Project Information  ",font=("Helvetica",11,"bold"),fg=self.ACCENT_DARK,bg=self.BG_SECONDARY,bd=1,relief=tk.GROOVE,padx=15,pady=12)
        left.pack(side=tk.LEFT,fill=tk.BOTH,expand=True,padx=(0,8))
        self._add_form_field(left,"Name of the project:",self.project_name,0)
        self._add_form_field(left,"Modules required:",self.modules_required,1,combo_values=self.MODULES)
        self._add_form_field(left,"Project version:",self.project_version,2)
        self._add_form_field(left,"Maintainer email address:",self.maintainer_email,3)
        right=tk.LabelFrame(body,text="  Computational Modules  ",font=("Helvetica",11,"bold"),fg=self.ACCENT_DARK,bg=self.BG_SECONDARY,bd=1,relief=tk.GROOVE,padx=15,pady=12)
        right.pack(side=tk.RIGHT,fill=tk.BOTH,expand=True,padx=(8,0))
        mf=tk.Frame(right,bg=self.BG_SECONDARY); mf.pack(fill=tk.X,pady=(5,10))
        bs=dict(font=("Helvetica",10),width=14,height=2,relief=tk.RAISED,bd=1,cursor="hand2")
        for text,r,c,col in [("DUNE",0,0,"#3498db"),("GMSH",0,1,"#2ecc71"),("PETSc",1,0,"#e74c3c"),("SLEPc",1,1,"#9b59b6"),("LAPACK",2,0,"#f39c12"),("SuperLU",2,1,"#1abc9c")]:
            tk.Button(mf,text=text,bg=col,fg="white",activebackground=col,activeforeground="white",command=lambda m=text:self._check_module(m),**bs).grid(row=r,column=c,padx=6,pady=4,sticky="ew")
        mf.columnconfigure(0,weight=1); mf.columnconfigure(1,weight=1)
        tk.Button(right,text="Input Parameters",font=("Helvetica",10,"bold"),bg="#34495e",fg="white",height=2,activebackground="#2c3e50",activeforeground="white",cursor="hand2",command=self._edit_input_params).pack(fill=tk.X,pady=(5,0))
        wf=tk.LabelFrame(outer,text="  Simulation Workflow  ",font=("Helvetica",11,"bold"),fg=self.ACCENT_DARK,bg=self.BG_SECONDARY,bd=1,relief=tk.GROOVE,padx=15,pady=10)
        wf.pack(fill=tk.X,pady=(15,0))
        wi=tk.Frame(wf,bg=self.BG_SECONDARY); wi.pack(fill=tk.X)
        for text,cmd,col in [("GMSH\n(Generate Mesh)",self._run_gmsh,"#2ecc71"),("Solve\n(Run Simulation)",self._run_solver,"#e74c3c"),("ParaView\n(Visualize Results)",self._open_paraview,"#3498db")]:
            tk.Button(wi,text=text,font=("Helvetica",10,"bold"),bg=col,fg="white",activebackground=col,activeforeground="white",width=22,height=3,cursor="hand2",command=cmd).pack(side=tk.LEFT,padx=8,pady=5,expand=True,fill=tk.X)

    def _create_status_bar(self):
        sf=tk.Frame(self.root,bg=self.ACCENT,height=28); sf.pack(side=tk.BOTTOM,fill=tk.X); sf.pack_propagate(False)
        self.status_label=tk.Label(sf,text="  QuDSim v{} | Ready".format(self.VERSION),font=("Helvetica",9),fg="white",bg=self.ACCENT,anchor=tk.W)
        self.status_label.pack(side=tk.LEFT,padx=10)

    def update_status(self,msg):
        t=datetime.now().strftime("%H:%M:%S")
        self.status_label.config(text="  QuDSim v{} | {} [{}]".format(self.VERSION,msg,t)); self.root.update_idletasks()

    def _add_form_field(self,parent,label_text,var,row,combo_values=None):
        tk.Label(parent,text=label_text,font=("Helvetica",10),fg=self.TEXT_PRIMARY,bg=self.BG_SECONDARY,anchor=tk.W).grid(row=row,column=0,sticky=tk.W,pady=6)
        if combo_values:
            w=ttk.Combobox(parent,textvariable=var,values=combo_values,font=("Helvetica",10),width=30,state="readonly")
        else:
            w=tk.Entry(parent,textvariable=var,font=("Helvetica",10),width=32,relief=tk.SOLID,bd=1)
        w.grid(row=row,column=1,sticky=tk.EW,pady=6,padx=(10,0)); parent.columnconfigure(1,weight=1)

    def _check_module(self,name):
        self.update_status("Checking {} ...".format(name))
        checks={"DUNE":"dunecontrol","GMSH":"gmsh"}; cmd=checks.get(name)
        found=False
        if cmd:
            try: found=subprocess.run(["which",cmd],capture_output=True).returncode==0
            except: pass
        else:
            env_map={"PETSc":"PETSC_DIR","SLEPc":"SLEPC_DIR","LAPACK":"LAPACK_DIR","SuperLU":"SUPERLU_DIR"}
            found=bool(os.environ.get(env_map.get(name,""),""))
        if found: self.update_status("{} available".format(name)); messagebox.showinfo("Module Check","{} is installed.".format(name))
        else: self.update_status("{} not found".format(name)); messagebox.showwarning("Module Check","{} not detected.\nPlease install it and set environment variables.".format(name))

    def _run_gmsh(self):
        if not self.current_geo_file:
            f=filedialog.askopenfilename(title="Select .geo File",filetypes=[("GMSH Geometry","*.geo"),("All","*.*")])
            if not f: return
            self.current_geo_file=f
        self.update_status("Generating mesh..."); out=self.current_geo_file.replace(".geo",".msh")
        try:
            r=subprocess.run(["gmsh",self.current_geo_file,"-3","-o",out,"-format","msh2"],capture_output=True,text=True,timeout=300)
            if r.returncode==0: self.current_msh_file=out; self.update_status("Mesh generated"); messagebox.showinfo("GMSH","Mesh generated:\n{}".format(out))
            else: messagebox.showerror("GMSH Error",r.stderr[:500])
        except FileNotFoundError: messagebox.showerror("Error","GMSH not found. Install from https://gmsh.info")
        except subprocess.TimeoutExpired: messagebox.showerror("Timeout","Mesh generation timed out.")

    def _run_solver(self):
        if not self.project_name.get(): messagebox.showwarning("Info","Enter a project name first."); return
        if not self.input_param_file:
            self.input_param_file=filedialog.askopenfilename(title="Select Input File",filetypes=[("Text","*.txt *.inp *.dat"),("All","*.*")])
            if not self.input_param_file: return
        np=self.mpi_processes.get(); mod=self.modules_required.get()
        self.update_status("Running {} with {} MPI procs...".format(mod,np))
        solver_map={"Schrodinger-Poisson":"./qudsim_sp","Tunneling Current":"./qudsim_tunnel","Dirac Equation":"./qudsim_dirac","Poisson Only":"./qudsim_poisson","Open Boundary":"./qudsim_open"}
        cmd=["mpirun","-np",str(np),solver_map.get(mod,"./qudsim_sp"),self.input_param_file]
        lw=tk.Toplevel(self.root); lw.title("Solver Output"); lw.geometry("700x450")
        tw=tk.Text(lw,font=("Courier",10),bg="#1e1e1e",fg="#00ff00",wrap=tk.WORD); tw.pack(fill=tk.BOTH,expand=True,padx=5,pady=5)
        tw.insert(tk.END,"$ {}\n\n".format(" ".join(cmd)))
        try:
            p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
            def rd():
                l=p.stdout.readline()
                if l: tw.insert(tk.END,l); tw.see(tk.END); lw.after(50,rd)
                else:
                    rc=p.poll()
                    if rc is not None: tw.insert(tk.END,"\n--- Exit code {} ---\n".format(rc)); self.update_status("Solver done")
                    else: lw.after(100,rd)
            lw.after(100,rd)
        except FileNotFoundError: tw.insert(tk.END,"ERROR: mpirun or solver not found.\n")

    def _open_paraview(self):
        if not self.current_vtu_file:
            f=filedialog.askopenfilename(title="Select .vtu",filetypes=[("VTK","*.vtu *.vtk"),("All","*.*")])
            if not f: return
            self.current_vtu_file=f
        try: subprocess.Popen(["paraview",self.current_vtu_file]); self.update_status("ParaView launched")
        except FileNotFoundError: messagebox.showerror("Error","ParaView not found. Install from https://www.paraview.org")

    def _configure_mpi(self):
        d=tk.Toplevel(self.root); d.title("MPI Config"); d.geometry("350x200"); d.resizable(False,False); d.transient(self.root); d.grab_set()
        f=tk.Frame(d,bg=self.BG_SECONDARY,padx=20,pady=20); f.pack(fill=tk.BOTH,expand=True)
        tk.Label(f,text="Number of MPI Processes:",font=("Helvetica",11),bg=self.BG_SECONDARY).pack(anchor=tk.W,pady=(0,5))
        tk.Scale(f,from_=1,to=64,orient=tk.HORIZONTAL,variable=self.mpi_processes,font=("Helvetica",10),length=280,tickinterval=16).pack(fill=tk.X,pady=5)
        tk.Button(f,text="OK",font=("Helvetica",10,"bold"),bg=self.ACCENT,fg="white",width=10,command=d.destroy).pack(pady=10)

    def _edit_input_params(self):
        if not self.input_param_file or not os.path.exists(self.input_param_file or ""):
            self.input_param_file=filedialog.askopenfilename(title="Select Input File",filetypes=[("Text","*.txt *.inp *.dat"),("All","*.*")])
            if not self.input_param_file: return
        ed=tk.Toplevel(self.root); ed.title("Params - {}".format(os.path.basename(self.input_param_file))); ed.geometry("650x500")
        tb=tk.Frame(ed,bg=self.BG_PRIMARY); tb.pack(fill=tk.X)
        ta=tk.Text(ed,font=("Courier",11),wrap=tk.NONE,undo=True); ta.pack(fill=tk.BOTH,expand=True)
        tk.Button(tb,text="Save",command=lambda:(open(self.input_param_file,"w").write(ta.get("1.0",tk.END)),self.update_status("Saved"))).pack(side=tk.LEFT,padx=5,pady=3)
        if os.path.exists(self.input_param_file):
            with open(self.input_param_file) as f: ta.insert(tk.END,f.read())

    def _open_material_db(self):
        w=tk.Toplevel(self.root); w.title("Material Database"); w.geometry("600x400")
        mats=[("Si","11.7","1.12","0.916/0.190","6"),("SiO2","3.9","9.0","0.50","-"),("HfO2","25.0","5.8","0.18","-"),("Al2O3","9.0","6.5","0.35","-"),("InAs","15.15","0.354","0.023","1"),("In0.53Ga0.47As","13.9","0.74","0.041","1"),("GaN","8.9","3.39","0.20","1"),("Graphene","-","0.0","0 (Dirac)","-")]
        tr=ttk.Treeview(w,columns=("name","eps","bg","mass","v"),show="headings",height=15)
        for c,h,wd in [("name","Material",130),("eps","Dielectric",90),("bg","Bandgap(eV)",90),("mass","m*/m0",130),("v","Valleys",70)]:
            tr.heading(c,text=h); tr.column(c,width=wd,anchor=tk.CENTER)
        for m in mats: tr.insert("",tk.END,values=m)
        tr.pack(fill=tk.BOTH,expand=True,padx=10,pady=10)

    def _show_docs(self):
        import webbrowser; webbrowser.open("https://github.com/qudsimcnrnano/QuDSimprojects")

    def _show_about(self):
        messagebox.showinfo("About","QuDSim: 3D FEM Tool v{}\n\nOpen-source simulation tool for nanoscale devices.\nIntegrates: DUNE, GMSH, PETSc, SLEPc, LAPACK, SuperLU\n\ngithub.com/qudsimcnrnano/QuDSimprojects\n\nSupported by SERB Grant CRG/2022/009394".format(self.VERSION))

    def _save_config(self):
        try: json.dump({"project_name":self.project_name.get(),"modules_required":self.modules_required.get(),"project_version":self.project_version.get(),"maintainer_email":self.maintainer_email.get(),"mpi_processes":self.mpi_processes.get(),"geo":self.current_geo_file,"msh":self.current_msh_file,"vtu":self.current_vtu_file,"params":self.input_param_file},open(self.CONFIG_FILE,"w"),indent=2)
        except: pass

    def _load_config(self):
        if not os.path.exists(self.CONFIG_FILE): return
        try:
            c=json.load(open(self.CONFIG_FILE))
            self.project_name.set(c.get("project_name","")); self.modules_required.set(c.get("modules_required","Schrodinger-Poisson"))
            self.project_version.set(c.get("project_version","1.0")); self.maintainer_email.set(c.get("maintainer_email",""))
            self.mpi_processes.set(c.get("mpi_processes",4)); self.current_geo_file=c.get("geo"); self.current_msh_file=c.get("msh")
            self.current_vtu_file=c.get("vtu"); self.input_param_file=c.get("params"); self.update_status("Session restored")
        except: pass

    def _on_close(self): self._save_config(); self.root.destroy()


def main():
    root=tk.Tk(); app=QuDSimGUI(root); root.protocol("WM_DELETE_WINDOW",app._on_close); root.mainloop()

if __name__=="__main__": main()
