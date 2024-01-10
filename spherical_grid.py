import numpy as np
import time


R = np.loadtxt("R.txt",delimiter=",")
Theta = np.linspace(0,180,19) 
Phi = np.linspace(4.5,351.5,40)

with open("sph_grid.txt","w") as f:
    f.write("X,Y,Z,R,Theta,Phi\n")
    for Ri in R:
        for Ti in Theta:
            for Pi in Phi:
                X = Ri* np.cos(Ti)
                Y = Ri* np.sin(Ti) * np.sin(Pi)
                Z = Ri* np.sin(Ti) * np.cos(Pi)

                text = str(X)+","+str(Y)+","+str(Z)+","+str(Ri)+","+str(Ti)+","+str(Pi)+ "\n"
                f.write(text)
                #print(text)
            
                #print(X,Y,Z)
                #time.sleep(.2)
print("Done")
