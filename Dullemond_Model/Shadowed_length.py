import matplotlib.pyplot as plt

# Create a figure and an axes
fig, ax = plt.subplots()


 
# single line
plt.vlines(x = 0.5036194069000989,ymin = -0.06850714268787533, ymax = 0.06850714268787533, 
           colors = 'purple',
           #label = 'Blackbody'
           )


plt.vlines(x = 0.38343241097708336, ymin = -0.046270833860451235, ymax = 0.046270833860451235,
           colors = 'Red',
           #label = '[Me/H] =-0.3'
           )
 
plt.vlines(x = 0.29863670625168576, ymin = -0.0321643566662279, ymax = 0.0321643566662279,
           colors = 'teal',
           #label = '[Me/H] =0'
           )
 
plt.vlines(x = 0.17377332077474755, ymin = -0.014617283481862628, ymax = 0.014617283481862628,
           colors = 'green',
           #label = '[Me/H] =0.3'
           )
 



plt.hlines(y = -0.06850714268787533,xmin = 0.5036194069000989, xmax =0.5036194069000989+1.2661840790799868,
           colors = 'purple')



plt.hlines(y = -0.06850714268787533,xmin = 0.5036194069000989, xmax =0.5036194069000989+1.2661840790799868,
           colors = 'purple',
           label = 'No absorption')

plt.hlines( y = -0.046270833860451235, xmin = 0.38343241097708336, xmax = 0.38343241097708336+0.6198597764820348,
           colors = 'Red',
           label = '[Me/H] =-0.3')

plt.hlines(y = -0.0321643566662279,xmin = 0.29863670625168576, xmax = 0.29206866153546635 + 0.29863670625168576 , 
           colors = 'teal',
           label = '[Me/H] =0')
 
plt.hlines( y = -0.014617283481862628,xmin =0.17377332077474755, xmax =  0.01934359390471044+0.17377332077474755,
           colors = 'green',
           label = '[Me/H] =0.3')

ax.set_xlabel('R')
ax.set_ylabel('$H_s$')

# Add a legend in the upper right corner of the plot
plt.legend(loc='upper right')

 #bbox_to_anchor = (1.0, 1)
plt.show()

#print(Metallicity, R_rim/AU,H_rim/AU, shadowed_length/AU, gray_const)
#0 0.29863670625168576 0.0321643566662279 0.29206866153546635 0.3612614275488969
##BB  0.5036194069000989 0.06850714268787533 1.2661840790799868
#0.3 0.17377332077474755  0.01934359390471044 0.12498405485216646

#-0.3 0.38343241097708336 0.046270833860451235 0.6198597764820348 0.5886475189298194