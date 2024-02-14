import numpy as np
import time
import tecplot as tp

"""
A python script to spherically sample populations and their den vel temp for Lyman Alpha glow

If you are running this script you need to make sure your enviroment variables are correcly set
For me, i was able to get this to run by using:

"tec360-env -- python spherical_grid.py"

"""

all_pop = True
mydata = "BU_fastslow_m_1.plt"


#in this coordinate system we have
#R as radius from sun
#Theta as the angle between X and R
#Phi as the rotation around X (YZ plane)
                
R = np.loadtxt("R.txt",delimiter=",")[:81]
#R = np.loadtxt("R.txt",delimiter=",")[:31]
Theta = np.linspace(0,180,19)*np.pi/180
#Theta = np.array([0,np.pi])
Phi = np.linspace(4.5,351.5,40)*np.pi/180
#Phi = np.array([0])

R_upper_sample = 36
R_lower_sample = 2
R_traceback = 70 * 1.496e11

print("1")

tp.session.acquire_license()
print("2")
data3d = tp.data.load_tecplot(mydata)
print("3")
#print(data3d)






def indexs(pop):
    #this function returns the indexes of the neutral populations in a list
    pop_dict = {1:3,2:2,3:1,4:4}
    pop_i = pop_dict[pop]
    base = 28 + (pop_i-1)*5
    
    den_i = base
    vx_i, vy_i, vz_i = base+1,base+2,base+3
    temp_i = base+4

    return [den_i,vx_i,vy_i,vz_i,temp_i]


def all_pops():
    with open("sph_grid.txt","w") as f, open("pop1.txt","w") as p1, open("pop2.txt","w") as p2, open("pop3.txt","w") as p3, open("pop4.txt","w") as p4:
        #opening files to write to:
        #f for the grid, p# for each population

        f.write("X,Y,Z,R,Theta,Phi\n")
        p1.write("R,Theta,Phi,Den,Vx,Vy,Vz,P\n")
        p2.write("R,Theta,Phi,Den,Vx,Vy,Vz,P\n")
        p3.write("R,Theta,Phi,Den,Vx,Vy,Vz,P\n")
        p4.write("R,Theta,Phi,Den,Vx,Vy,Vz,P\n")

        #iterating through each position variables

        for i,Ri in enumerate(R):
            #checking progress
            print("R = ", + str(R),"AU:  "str(i/len(R) * 100)[:4],"%",sep="")
            Ri = Ri
            for Ti in Theta:
                for Pi in Phi:
                    
                    X = Ri* np.cos(Ti) * 1.496e11
                    Y = Ri* np.sin(Ti) * np.sin(Pi) * 1.496e11
                    Z = Ri* np.sin(Ti) * np.cos(Pi) * 1.496e11

                    #print("\n",X/1.496e11,Y,Z)
                    textgrid = str(X)+","+str(Y)+","+str(Z)+","+str(Ri)+","+str(Ti)+","+str(Pi)+ "\n"
                    f.write(textgrid)

                    allpop_string = []
                    pop_string_base = str(Ri)+","+str(Ti*180/np.pi)+","+str(Pi*180/np.pi)+ ","
                    
                    sampled_data = np.array(tp.data.query.probe_at_position(X,Y,Z,dataset=data3d)[0])

                    for pop_i in range(1,5):

                        #take the indexes for the varibles
                        #make array of data
                        pop_is = indexs(pop_i)
                        pop_sampled_data = sampled_data[pop_is]
                        
                        if (Ri < R_upper_sample) and (pop_i == 1 or pop_i == 4):
                            ###SIMPLE SOLUTION HERE
                            simple = evolve([X,Y,Z],R_traceback,pop_i)
                            fraction = (np.log(Ri)-np.log(R_lower_sample))/np.log(R_upper_sample/R_lower_sample)
                            if fraction>1:
                                fraction =1
                            if fraction<0:
                                fraction = 0
                            #print("sampled density -->\t",pop_sampled_data[0]/(100**3))
                            #print("Simple algo -->\t", simple/(100**3))
                            #print("Fraction of sampled to algo -->\t",fraction)
                            #print(pop_sampled_data[0],simple,fraction)
                            #fraction = 0
                            pop_sampled_data[0] = pop_sampled_data[0] * fraction + simple* (1-fraction)
                            #print("final density -->\t",pop_sampled_data[0]/(100**3))
                        
                        popstringdata = pop_string_base
                        for data_i in pop_sampled_data:
                            popstringdata+=(str(data_i)+",")
                        #remove final "," and replace with "endline (\n)"
                        popstringdata = popstringdata[:-1]+"\n"
                        allpop_string.append(popstringdata)

                    

                        
                    p1.write(allpop_string[0])
                    p2.write(allpop_string[1])
                    p3.write(allpop_string[2])
                    p4.write(allpop_string[3])



    print("Done")

def traceback(pos,Rlim,pop):

    Xi,Yi,Zi = pos
    sample = np.array(tp.data.query.probe_at_position(Xi,Yi,Zi,dataset=data3d)[0])

    pop_is = indexs(pop)
    Den,Vx,Vy,Vz,T = sample[pop_is]

    #    Vx_Id,Vy_Id,Vz_Id = id("vel")
    #    Vx,Vy,Vz = sample[Vx_Id],sample[Vy_Id],sample[Vz_Id]
    
    V_total = np.sqrt(Vx**2 + Vy**2 + Vz**2)
    '''
    print(sample)
    print(X/1.496e11,Y/1.496e11,Z/1.496e11)
    print(Vx/1000,Vy/1000,Vz/1000)
    '''

    
    R = np.sqrt(Xi**2 + Yi**2 + Zi**2)
    R_list = [R]
    Beta_list = [beta(R/1.496e11)]
    V_list = [V_total]
    
    
    dt = 1E6 #seconds
    while R < Rlim:
        #Propogate position backwards (minus V*dt)
        Xi,Yi,Zi = Xi - Vx*dt,Yi - Vy*dt,Zi - Vz*dt
        R = np.sqrt(Xi**2 + Yi**2 + Zi**2)
        R_list.append(R)
        V_list.append(V_total)
        Beta_list.append(beta(R/1.496e11))
        
    last_sample = np.array(tp.data.query.probe_at_position(Xi,Yi,Zi,dataset=data3d)[0])
    Den = last_sample[pop_is[0]]
    #print(len(sample_list))
    #print(len(R_list))
    return R_list,V_list,Beta_list,Den

def evolve(Pos_intial,R_lim,pop):

    R_L, V_L, B_L,Den = traceback(Pos_intial,R_lim,pop) 
        
    #subtract first position from second position
    Delta_r = np.abs(np.array(R_L[:-1]) -np.array( R_L[1:]))
    Delta_t = Delta_r/V_L[:-1]
    #print(Delta_t)
    
    #print(V_L[0],R_L[0]/1.496e11,sum(Delta_t))
    Beta_r =B_L[:-1]

    Beta_total = np.sum(Delta_t*Beta_r)
    extinction = np.exp(-Beta_total)
    Den_final = extinction* Den
    '''
    for R_i,V_i,B_i in zip(R_L,V_L,B_L):
        print(R_i/1.496e11,V_i/1000,B_i)
    '''
    #print(extinction)
    return Den_final
    
def beta(r):
    #Need to impliment CX fr fr
    B_E = 6.2E-7
    B_r = B_E* (1/(r+.01)**2)
    return B_r

    

    
all_pops()
