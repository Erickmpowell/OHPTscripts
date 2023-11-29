import numpy as np


def parsedat(path):
    with open(path, "r") as f:
        lines = f.readlines()

    if "variable" in lines[0].lower():
        beginningline = 0
    else:
        beginningline = 1 
    
    vars = []
    headerlen = 0
    for line_i in lines[beginningline:]:

        if "ZONE" == line_i[:4]:
            break

        headerlen += 1 
        varsplit = line_i.split('"')[-2]
        vars.append(varsplit.lower())

    for line_i in lines[headerlen:]:
        headerlen += 1
        if "DT=(" in line_i:
            break
    return vars, headerlen

def parseout(path):
    with open(path, "r") as f:
        lines = f.readlines()

    vars_line = lines[4].lower()
    vars = vars_line.split()

    beginningline = 1 
    headerlen = 5
    return vars, headerlen

def getSWMFdata(Path, Configuration="BATSRUS"):

    if ".dat" in Path[-4:]:
        varlist, headerlen = parsedat(Path)
    if ".out" in Path[-4:]:
        varlist, headerlen = parseout(Path)
        
    print(varlist, headerlen)

    data = np.loadtxt(Path, unpack=True, skiprows=headerlen)
    print(data)

    if Configuration == "OHPT":
        data_class = OHPTdata(data,varlist)
        
    if Configuration == "BATSRUS":
        data_class = BATSRUSdata_SI(data,varlist)
        
    return data_class


class PlasmaFluid_OH:
    def __init__(self, data, varlist):

        den_i = varlist.index("rho amu/cm3")
        vx_i  = varlist.index("u_x km/s")
        vy_i  = varlist.index("u_y km/s")
        vz_i  = varlist.index("u_z km/s")
        p_i   = varlist.index("p dyne/cm^2")
        bx_i  = varlist.index("b_x nt")
        by_i  = varlist.index("b_y nt")
        bz_i  = varlist.index("b_z nt")
        HP_i  = varlist.index( "hplim amu/cm3")
        
        self.den = data[den_i]
        self.vx = data[vx_i]
        self.vy = data[vy_i]
        self.vz = data[vz_i]
        self.p =  data[p_i]
        self.bx = data[bx_i]
        self.by = data[by_i]
        self.bz = data[bz_i]
        self.HP = data[HP_i]/data[den_i]

    def B(self):
        return np.sqrt((self.bx)**2 + (self.by)**2 + (self.bz)**2)


    def vel(self):
        return np.sqrt((self.vx)**2 + (self.vy)**2 + (self.vz)**2)

    def temp(self):
        return (7.2464E15 * self.p/self.den)

    def Cs(self):
        return 1E-5 * np.sqrt( 5/3 * self.p / (self.den * 1.673E-24) )

    def Ca(self):
        return (self.B()*1E-10 / (np.sqrt(12.56 * self.den * 1.673E-24) ))
    
    def M_sonic(self):
        return self.vel()/self.Cs()
    
    def M_alfen(self):
        return self.vel()/self.Ca()
    
    def M_ms(self):
        return self.vel()/ np.sqrt(self.Ca()**2 + self.Cs()**2)

    
class NeutralFluid_OH:
    def __init__(self, data, varlist,pop_i):

        if pop_i == 1:
            pop_i = "u"
        else:
            pop_i = str(pop_i)
            
        den_i = varlist.index("rho^ne" + pop_i + " amu/cm3")
        vx_i  =varlist.index("u_x^ne" + pop_i + " km/s")
        vy_i =varlist.index("u_y^ne" + pop_i + " km/s")
        vz_i =varlist.index("u_z^ne" + pop_i + " km/s")
        p_i = varlist.index("p^ne" + pop_i + " dyne/cm^2")
        
        
        self.den = data[den_i]
        self.vx = data[vx_i]
        self.vy = data[vy_i]
        self.vz = data[vz_i]
        self.p =  data[p_i]


    def vel(self):
        return np.sqrt((self.vx)**2 + (self.vx)**2 + (self.vx)**2)

    def temp(self):
        return (self.p/self.den)


class NeutralFluid_FLEKS:
    def __init__(self, data, varlist,pop_i):

        pop_i = str(pop_i)
            
        den_i = varlist.index("rhopop" + pop_i )
        vx_i  = varlist.index("uxpop" + pop_i )
        vy_i = varlist.index("uypop" + pop_i )
        vz_i = varlist.index("uzpop" + pop_i)
        p_xx_i = varlist.index("pxxpop" + pop_i ) 
        p_yy_i = varlist.index("pyypop" + pop_i )
        p_zz_i = varlist.index("pzzpop" + pop_i )

        try:
            ppc_i = varlist.index("ppcpop"+pop_i)
        except:
            ppc_i = varlist.index("numpop"+pop_i)
    
        self.pxx = data[p_xx_i]
        self.pyy = data[p_yy_i]
        self.pzz = data[p_zz_i]
        self.popi = pop_i
        self.den = data[den_i]
        self.vx = data[vx_i]
        self.vy = data[vy_i]
        self.vz = data[vz_i]

        try:
            p_i = varlist.index("ppop" + pop_i )
            
            self.p = data[p_i]

        except:
            self.p =  1/3*(data[p_xx_i] + data[p_yy_i] + data[p_zz_i])

        

    def vel(self):
        return np.sqrt((self.vx)**2 + (self.vy)**2 + (self.vz)**2)

    def temp(self,pop2factor=1):
        if self.popi == "2":
            temp =np.nan_to_num( pop2factor * (self.p /self.den) * 1E8/ 1.38)
        else:
            temp =np.nan_to_num( (self.p /self.den) * 1E8/ 1.38)
        return temp
    

class BATSRUSdata_SI:
    def __init__(self, data,varlist):
        
        self.x = data[varlist.index("x au")]
        self.y = data[varlist.index("y au")]
        self.z = data[varlist.index("z au")]
        self.plasma = PlasmaFluid_OH(data,varlist)
        self.n1 = NeutralFluid_OH(data,varlist,pop_i = 1)
        self.n2 = NeutralFluid_OH(data,varlist,pop_i = 2)
        self.n3 = NeutralFluid_OH(data,varlist,pop_i = 3)
        self.n4 = NeutralFluid_OH(data,varlist,pop_i = 4)

    
class OHPTdata:
    def __init__(self, data,varlist):
        try:
            self.x = data[varlist.index("x")]
        except:
            pass
        try:
            self.y = data[varlist.index("y")]
        except:
            pass
        try:
            self.z = data[varlist.index("z")]
        except:
            pass
        self.n1 = NeutralFluid_FLEKS(data,varlist,pop_i = 1)
        self.n2 = NeutralFluid_FLEKS(data,varlist,pop_i = 2)
        self.n3 = NeutralFluid_FLEKS(data,varlist,pop_i = 3)
        self.n4 = NeutralFluid_FLEKS(data,varlist,pop_i = 4)

    def p_total(self):
        return self.n1.p + self.n2.p + self.n3.p + self.n4.p
    
    def den_total(self):
        return self.n1.den + self.n2.den + self.n3.den + self.n4.den
    
    def vel_total(self):
        velocity = (self.n1.den * self.n1.vel() + self.n2.den* self.n2.vel() + self.n3.den * self.n3.vel() + self.n4.den * self.n4.vel())/self.den_total()
        return velocity 

    def vx_total(self):
        velocity = (self.n1.den * self.n1.vx + self.n2.den* self.n2.vx + self.n3.den * self.n3.vx + self.n4.den * self.n4.vx)/self.den_total()
        return velocity

    
    def vy_total(self):
        velocity = (self.n1.den * self.n1.vy + self.n2.den* self.n2.vy + self.n3.den * self.n3.vy + self.n4.den * self.n4.vy)/self.den_total()
        return velocity 

    
    def vz_total(self):
        velocity = (self.n1.den * self.n1.vz + self.n2.den* self.n2.vz + self.n3.den * self.n3.vz + self.n4.den * self.n4.vz)/self.den_total()
        return velocity 

    
    def true_Pxx(self):
        Full_pxx = (self.n1.pxx + self.n2.pxx + self.n3.pxx + self.n4.pxx)*1E-9   +\
        (self.n1.den * (self.n1.vx - self.vx_total())**2 +\
        self.n2.den * (self.n2.vx - self.vx_total())**2 +\
        self.n3.den * (self.n3.vx - self.vx_total())**2 +\
        self.n4.den * (self.n4.vx - self.vx_total())**2)* 1.67e-27* (1000**2)* (100**3)

        return Full_pxx

    def true_Pyy(self):
        Full_pyy = (self.n1.pyy + self.n2.pyy + self.n3.pyy + self.n4.pyy)*1E-9  +\
        (self.n1.den * (self.n1.vy - self.vy_total())**2 +\
        self.n2.den * (self.n2.vy - self.vy_total())**2 +\
        self.n3.den * (self.n3.vy - self.vy_total())**2 +\
        self.n4.den * (self.n4.vy - self.vy_total())**2)* 1.67e-27* (1000**2)* (100**3)

        return Full_pyy
    
    def true_Pzz(self):
        Full_pzz = (self.n1.pzz + self.n2.pzz + self.n3.pzz + self.n4.pzz)*1E-9  +\
        (self.n1.den * (self.n1.vz - self.vz_total())**2 +\
        self.n2.den * (self.n2.vz - self.vz_total())**2 +\
        self.n3.den * (self.n3.vz - self.vz_total())**2 +\
        self.n4.den * (self.n4.vz - self.vz_total())**2)* 1.67e-27* (1000**2)* (100**3)

        return Full_pzz

    def true_temp(self):
        total_pressure = 1/3*(self.true_Pxx() + self.true_Pyy() + self.true_Pzz()) 
        temp = total_pressure/(self.den_total()* (100**3)* 1.38E-23) 
        return temp

    def temp_total(self,pop2factor=1):
        temp = (self.n1.den * self.n1.temp() + self.n2.den* self.n2.temp(pop2factor) + self.n3.den * self.n3.temp() + self.n4.den * self.n4.temp())/self.den_total()
        return temp

try:
    #For debugging and so I dont break it when I push the code lol
    #BATS = getSWMFdata("BATSRUS.dat")
    #FLEKS1 = getSWMFdata("FLEKS_line_1.dat","OHPT")
    #FLEKS2 = getSWMFdata("FLEKS_Kin_line.dat","OHPT")
    #FLEKS_out = getSWMFdata("cut.out","OHPT")
    #print(FLEKS1.n1.den[:10],FLEKS2.n1.den[:10])
    pass
except:
    pass
