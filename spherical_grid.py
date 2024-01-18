import numpy as np
import time
import tecplot as tp

"""
A python script to spherically sample populations and their den vel temp for Lyman Alpha glow

If you are running this script you need to make sure your enviroment variables are correcly set
For me, i was able to get this to run by using:

"tec360-env -- python spherical_grid.py"

"""

#in this coordinate system we have
#R as radius from sun
#Theta as the angle between X and R
#Phi as the rotation around X (YZ plane)
                
R = np.loadtxt("R.txt",delimiter=",")
Theta = np.linspace(0,180,19) 
Phi = np.linspace(4.5,351.5,40)

tp.session.acquire_license()    
data3d = tp.data.load_tecplot("3d.plt")
print(data3d)



def indexs(pop_i):
    #this function returns the indexes of the neutral populations in a list
    base = 29 + (pop_i-1)*5
    
    den_i = base
    vx_i, vy_i, vz_i = base+1,base+2,base+3
    temp_i = base+3

    return [den_i,vx_i,vy_i,vz_i,temp_i]
            

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
        print(str(i/len(R) * 100)[:4],"%",sep="")
        
        for Ti in Theta:
            for Pi in Phi:

                X = Ri* np.cos(Ti)
                Y = Ri* np.sin(Ti) * np.sin(Pi)
                Z = Ri* np.sin(Ti) * np.cos(Pi)

                
                textgrid = str(X)+","+str(Y)+","+str(Z)+","+str(Ri)+","+str(Ti)+","+str(Pi)+ "\n"
                f.write(textgrid)
                
                sampled_data = np.array(tp.data.query.probe_at_position(X,Y,Z,dataset=data3d)[0])
                allpop_string = []
                pop_string_base = str(Ri)+","+str(Ti)+","+str(Pi)+ ","
                for pop_i in range(1,5):

                    #take the indexes for the varibles
                    #make array of data
                    pop_is = indexs(pop_i)
                    pop_data = sampled_data[pop_is]
                    
                    popstringdata = pop_string_base
                    for data_i in pop_data:
                        popstringdata+=(str(data_i)+",")
                    #remove final "," and replace with "endline (\n)"
                    popstringdata = popstringdata[:-1]+"\n"
                    allpop_string.append(popstringdata)
                    
                p1.write(allpop_string[0])
                p2.write(allpop_string[1])
                p3.write(allpop_string[2])
                p4.write(allpop_string[3])
                    
                
                
print("Done")
