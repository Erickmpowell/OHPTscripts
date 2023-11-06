import numpy as np


def parsedat(path):
    with open(path, "r") as f:
        lines = f.readlines()

    firstline = lines[0].split('"')
    
    vars = []
    headerlen = 0
    for line_i in lines:

        if "ZONE" == line_i[:4]:
            break

        headerlen += 1 
        varsplit = line_i.split('"')[-2]
        vars.append(varsplit.lower())

    for line_i in lines[headerlen:]:
        headerlen += 1
        if "DT=(" in line_i:
            break
    return vars, headerlen+2


def getSWMFdata(Path, Configuration="BATSRUS"):

    varlist, headerlen = parsedat(Path)
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
        try:
            ppc_i = varlist.index("ppcpop"+pop_i)
        except:
            ppc_i = varlist.index("numpop"+pop_i)
        p_xx_i = varlist.index("pxxpop" + pop_i ) 
        p_yy_i = varlist.index("pyypop" + pop_i )
        p_zz_i = varlist.index("pzzpop" + pop_i )

        
        self.den = data[den_i]
        self.vx = data[vx_i]
        self.vy = data[vy_i]
        self.vz = data[vz_i]
        self.p =  1/3*(data[p_xx_i] + data[p_xx_i] + data[p_xx_i])


    def vel(self):
        return np.sqrt((self.vx)**2 + (self.vx)**2 + (self.vx)**2)

    def temp(self):
        return (self.p * 1E-9/(self.den*100**3 * 1.38E-23))    


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


    def den_total(self):
        return self.n1.den + self.n2.den + self.n3.den + self.n4.den
    
    def vel_total(self):
        velocity = (self.n1.den * self.n1.vel() + self.n2.den* self.n2.vel() + self.n3.den * self.n3.vel() + self.n4.den * self.n4.vel())/self.den_total()
        return velocity 
    
    def temp_total(self):
        temp = (self.n1.den * self.n1.temp() + self.n2.den* self.n2.temp() + self.n3.den * self.n3.temp() + self.n4.den * self.n4.temp())/self.den_total()
        return temp
        
#BATS = getSWMFdata("BATSRUS.dat")
#FLEKS = getSWMFdata("FLEKS_Smooth_10_0_5","OHPT")
