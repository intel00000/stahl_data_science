# Adopt these values as required 

# functional/basis set choice, inclusion of D3 and/or other method-specific options

nbo = " pop=nbo7" #other choice: "pop=nbo"
nmr = " nmr=giao" 
polar = " polar"
prop = " prop=efg"
volume = " volume"
solvent = " SCRF=(Solvent=dichloromethane,SMD)"
tddft = " td=(50-50,nstates=12)"

#!new stuff
nbo7read = " pop=NBO7read"
nrt_string = "\n$nbo nrt $end\n"

method = {
"opt": "# B3LYP/6-31G(d,p) symmetry=loose empiricaldispersion=GD3BJ int=(grid=ultrafine)",
"sp1":  "# M062X/def2tzvp int=(grid=ultrafine)" + nbo + polar + prop + volume,
"sp2":  "# M062X/def2tzvp int=(grid=ultrafine)" + nmr,
"sp3": "# M062X/def2tzvp int=(grid=ultrafine)" + nbo7read + polar + prop + volume,
"splig": "# M06/def2TZVP int=(grid=ultrafine)" + solvent + nmr + nbo + prop
}

#ECP to be used in optimization
ecp = "LANL2DZ" #other common choice: SDD
#ECP to be used in single points
specp = "SDD"

chkfolder = "%chk="
link0 = "%nprocs=16\n%mem=16GB\n"

keywords = {
0:" opt freq=noraman",
1:" guess=read geom=check",
2:" guess=read geom=check",
}

#  
# 
#-----------------------------------------------------------------------------------------------------
# no changes should be required below here

import os,sys,copy

#for basis/ECP selection:
heavy = ["K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr"] #,"Br"
        
def get_coms(dir):
    coms = [file for file in os.listdir(dir) if file.endswith((".gjf")) or file.endswith((".com"))]
    return(coms)

def writefile(filecontent,file,dir):    
    f = open(dir+file[:-4] + ".com",'w', newline='\n')
    for line in filecontent:
        f.write("%s" %line)
    f.close()
    return

class Input:    
    def __init__(self,file,dir,route,jobs):
        self.file = dir+file
        self.name = file[:-4].split("/")[-1]
        if self.get_coords():
            self.sroute = copy.copy(route)
            self.basis = {}
            for i in range(jobs):    
                print(f"Job {i+1}: {self.sroute[i].strip()}")
                self.basis[i] = self.sroute[i].split("/")[1].split(" ")[0]
            self.metaldetector(jobs)
        else:
            self.name = False
        
    def get_coords(self):
        filecont = open(self.file,'r').readlines()
        if ".mae" in filecont[0]:
            return(False)
        self.coords,self.elements = "",[]
        start = len(filecont)
        for l in range(len(filecont)):
            if len(filecont[l].split()) == 2 and len(filecont[l+1].split()) == 4:
                start = l+1
                self.chsp = filecont[l].split()
                self.title = filecont[l-2]
            if len(filecont[l+1].split()) < 2 and l > start:
                end = l+1
                break
            if l == len(filecont)-2:
                return(False)
        for i in range (start,end,1):
            self.elements.append(filecont[i].split()[0].split("-")[0])
            self.coords += "    ".join(filecont[i].split()[:4])+"\n"
        return(True)
        
    def metaldetector(self,jobs):
        self.heavyelems = " ".join(set([x for x in self.elements if x in heavy]))
        self.otherelems = " ".join(set([x for x in self.elements if x not in heavy]))
        if len(self.heavyelems) != 0:
            for i in range(jobs):
                if "6-31G" in self.basis[i] or "6-31+G" in self.basis[i] or "6-31++G" in self.basis[i]:
                    self.sroute[i] = self.sroute[i].replace("/"+self.basis[i]," gen 6D pseudo=read")
                else:
                    self.sroute[i] = self.sroute[i].replace("/"+self.basis[i]," gen pseudo=read")
            self.heavy = 1
        else:
            self.heavy = 0
        return

def convert_gaussin(choice,dir):
    setup = {
    "1":["   Opt+freq at " + method["opt"].split(" ")[1] + "\n\n",["opt"]],
    "2":["   For Amide Project:\n    Opt+freq at " + method["opt"].split(" ")[1] + "\n    then Single point at " + method["sp1"].replace("# ","") + "\n    then Single point at " + method["sp2"].replace("# ","") + "\n\n",["opt","sp1","sp2"]],
    "3":["   For the complexes:\n    Single point at " + method["sp1"].replace("# ","") + "\n    then Single point at " + method["sp2"].replace("# ","") + "\n\n",["sp1","sp2"]],
    "4":["   For the ligands:\n    Single point at " + method["splig"].replace("# ","") + "\n\n",["splig"]],
    "5":["   For Chlorination Project:\n    Opt+freq at " + method["opt"].split(" ")[1] + "\n    then Single point at " + method["sp3"].replace("# ","") + " (WITH NRT)" + "\n    then Single point at " + method["sp2"].replace("# ","") + "\n\n",["opt","sp3","sp2"]],
    }
    
    text = "by default, ECP is selected for fourth and higher row elements\n" 
    for opt in sorted(list(setup.keys()))[::-1]:
        text = opt + setup[opt][0] + text
    
    while choice not in setup.keys():
        choice = input(text)
    jobs = len(setup[choice][1]) # number of subjobs
    # level = int(level)
    route = {}    
        
    if choice == "3":
        route[0] = method[setup[choice][1][0]] + "\n"*2
        route[1] = method[setup[choice][1][1]] + keywords[1] + "\n"*2
    elif choice == "4":
        route[0] = method[setup[choice][1][0]] + "\n"*2
    else:
        for i in range(jobs):
            route[i] = method[setup[choice][1][i]] + keywords[i] + "\n"*2

    coms = get_coms(dir)
    
    for com in coms:
        print(f"---------------------------------")
        print(f"Processing file: {com}")
        data = Input(com,dir,route,jobs)
        if not data.name:  # is not a Gaussian input file
            print(f"-------!File {com} is not a Gaussian input file!-------")
            continue
            
        filecontent = ""
        for i in range(jobs):
            if i == 0:
                putcoords = 1
            else:
                putcoords = 0
            if jobs-i-1 == 0:
                putlink = 0
            else:
                putlink = 1
            if "opt" in route[i]:
                pecp = ecp
            else:
                pecp = specp
                
            # !new stuff
            if nbo7read in route[i]: # check if nbo7read string is in the route
                put_nrt = 1
            else:
                put_nrt = 0
            
            filecontent += link0 + chkfolder + data.name + ".chk\n" + data.sroute[i] + data.title + "\n" + " ".join(data.chsp) + "\n" + data.coords*putcoords + nrt_string * put_nrt + "\n" + data.heavy*(data.otherelems + " 0\n" + data.basis[i] + "\n****\n" + data.heavyelems + " 0\n" + pecp + "\n****\n\n" + data.heavyelems + " 0\n" + pecp + "\n\n") + "--Link1--\n"*putlink
        with open(dir+com[:-4] + ".com",'w', newline='\n') as f:
            f.write(filecontent)

if __name__ == '__main__': 
    dir = "./"    
    choice = sys.argv[-1]
    convert_gaussin(choice,dir)
            
