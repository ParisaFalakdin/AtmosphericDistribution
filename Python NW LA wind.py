#This code has been initially written in VB6 environment as part of GriddedSoilPlusVeg model (Falakdin et al., 2022)
#This is a module that calculates atmospheric mole distribution among cells
#Several variables such as gridRow and gridColumn are global variables that are defined elsewhere, they are commented within the code for more context


import numpy as np
from math import cos, sin, pi

MOls2AddLAx = np.zeros((gridRows + 2, gridColumns + 2, 8760)) #array for the amount of mole added to the adjacent cell along x
MOls2AddLAy = np.zeros((gridRows + 2, gridColumns + 2, 8760)) #array for the amount of mole added to the adjacent cell along y
Motx_LA = np.zeros((gridRows + 2, gridColumns + 2)) 
Moty_LA = np.zeros((gridRows + 2, gridColumns + 2))
Ct_LA = np.zeros((gridRows + 1, gridColumns + 1)) #Atmospheric chemical concentration of the cell (temporary) without considering other sources

Gx_LA_prim = 0 #Temporary volumetric flow_rate along x direction 
Gy_LA_prim = 0 #Temporary volumetric flow_rate along y direction 

rtlax = 0 #Variable calculating the fraction volumetric flow rate along x over cell lower air volume 
rtlay = 0 #Variable calculating the fraction volumetric flow rate along y over cell lower air volume

nvala = 0 
length = 0

length = SimArea ** 0.5  #length of the cells
angle = alfa(tp)         #input time dependent variable: the angle in which wind enters the study area
wind = WindDirection(tp) #input time dependent variable: the direction in which wind enters the study area(e.g., NW, SW, etc.)

#Gx_LA Volumetric flow_rate along x direction
#Gy_LA Volumetric flow_rate along y direction
#gridColumn is the total number of columns
#gridRow is the total number of rows
#ii is the row number
#jj is the column number
#tp is the current hour

#**************************************************Lower Air***********************************************

if wind == "NW" or wind == "NE" or wind == "SE" or wind == "SW":

    rtlax = Gx_LA / VolLA(tp) # Shows how many times Gx_LA is larger than the cell volume i.e = 1/retention-time
    rtlay = Gy_LA / VolLA(tp) # Shows how many times Gy_LA is larger than the cell volume i.e = 1/retention-time

    MOls2AddLAx[ii, jj, tp] = Mols2AddLA[ii, jj, tp] * (cos(angle * pi / 180) ** 2) # mole going out of the cell lower air along x for the current hour
    MOls2AddLAy[ii, jj, tp] = Mols2AddLA[ii, jj, tp] * (sin(angle * pi / 180) ** 2) # mole going out of the cell lower air along y for the current hour

    # Take the integer part of the number add one to calculate the number of cells getting affected by wind along each direction
    if rtlax > 1:
        alax = int(rtlax) + 1
    else:
        alax = 1
    if rtlay > 1:
        alay = int(rtlay) + 1
    else:
        alay = 1

    Gx_LA_prim = Gx_LA / alax
    Gy_LA_prim = Gy_LA / alay

#******************************************NW***********************************************

    if wind == "NW":

        if jj != gridColumns:
            Mt_LA[ii, jj + 1] += MOls2AddLAx[ii, jj, tp]  # The first cell that chemical enters along x
        Motx_LA[ii, jj + 1] = MOls2AddLAx[ii, jj, tp]

        if ii != gridRows:
            Mt_LA[ii + 1, jj] += MOls2AddLAy[ii, jj, tp]  # The first cell that chemical enters along y
        Moty_LA[ii + 1, jj] = MOls2AddLAy[ii, jj, tp]
        
		#if the number of cells getting affected by wind along both x and y directions is more than one
		
        if alax > 1 or alay > 1:

            pp = min(alax, alay)
            qq = abs(alax - alay)
			
			#first distribute in both directions equally (also in case wind is entering with angle 45 degree)
			
            for p in range(1, pp):
				for i in range(1, gridRows + 1):
					for j in range(1, gridColumns + 1):
						Mt_LA[i][j] = Mt_LA[i][j] + Motx_LA[i][j] + Moty_LA[i][j]
						Ct_LA[i][j] = Mt_LA[i][j] / VolLA(tp)
						Motx_LA[i][j + 1] = Ct_LA[i][j] * Gx_LA_prim
						Moty_LA[i + 1][j] = Ct_LA[i][j] * Gy_LA_prim
						nvala = Motx_LA[i][j + 1] + Moty_LA[i + 1][j]
				
            if nvala > Mt_LA[i][j]:
                Motx_LA[i][j + 1] = Mt_LA[i][j] * (cos(angle * dblPi / 180) ** 2)
                Moty_LA[i + 1][j] = Mt_LA[i][j] * (sin(angle * dblPi / 180) ** 2)
                Mt_LA[i][j] = Mt_LA[i][j] - Motx_LA[i][j + 1] - Moty_LA[i + 1][j]
            else:
                Mt_LA[i][j] = Mt_LA[i][j] - Motx_LA[i][j + 1] - Moty_LA[i + 1][j]
			
			
			#second if along one direction wind is more powerful distribute along that direction (angles != 45)
			
			if qq > 0:

				if alax < alay:

					for q in range(1, qq):
						for i in range(1, gridRows + 1):
							for j in range(1, gridColumns + 1):	
								if i != 1:
									Mt_LA[i][j] = Mt_LA[i][j] + Moty_LA[i][j]  #mole existing in the cell + mole coming in along y axis
								Ct_LA[i][j] = Mt_LA[i][j] / VolLA(tp)
								Moty_LA[i + 1][j] = Ct_LA[i][j] * Gy_LA_prim   # mole going out along y axis
								nvala = Moty_LA[i + 1][j]
								if nvala > Mt_LA[i][j]:
									Moty_LA[i + 1][j] = Mt_LA[i][j]
									Mt_LA[i][j] = Mt_LA[i][j] - Moty_LA[i + 1][j]
								else:
									Mt_LA[i][j] = Mt_LA[i][j] - Moty_LA[i + 1][j]

				else:  

					for q in range(1, qq):
						for i in range(1, gridRows + 1):
							for j in range(1, gridColumns + 1):
								if j != 1:
									Mt_LA[i][j] = Mt_LA[i][j] + Motx_LA[i][j]
								Ct_LA[i][j] = Mt_LA[i][j] / VolLA(tp)
								Motx_LA[i][j + 1] = Ct_LA[i][j] * Gx_LA_prim
								nvala = Motx_LA[i][j + 1]
								if nvala > Mt_LA[i][j]:
									Motx_LA[i][j + 1] = Mt_LA[i][j]
									Mt_LA[i][j] = Mt_LA[i][j] - Motx_LA[i][j + 1]
								else:
									Mt_LA[i][j] = Mt_LA[i][j] - Motx_LA[i][j + 1]

		else            
		
		#Calculate the remaining in each cell
		
			for i in range(1, gridRows + 1):
				for j in range(1, gridColumns + 1):
					nvala = Motx_LA[i][j + 1] + Moty_LA[i + 1][j]
					if nvala > Mt_LA[i][j]:
						Motx_LA[i][j + 1] = Mt_LA[i][j] * (math.cos(angle * math.pi / 180) ** 2)
						Moty_LA[i + 1][j] = Mt_LA[i][j] * (math.sin(angle * math.pi / 180) ** 2)
						Mt_LA[i][j] = Mt_LA[i][j] - Motx_LA[i][j + 1] - Moty_LA[i + 1][j]
					else:
						Mt_LA[i][j] = Mt_LA[i][j] - Motx_LA[i][j + 1] - Moty_LA[i + 1][j]



# This block of code is written for one direction and only Lower Air compartment
#At the end of this module we will have is Mt_LA() for all cells at time tp!


