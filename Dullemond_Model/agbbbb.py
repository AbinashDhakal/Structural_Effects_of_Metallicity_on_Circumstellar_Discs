import matplotlib.pyplot as plt

# set the x-coordinate for the line
x = 5
x2 =10

# set the y-coordinates for the line
y1 = 0
y2 = 10

# create the line
plt.plot([x, x], [y1, y2], 'k--', lw=2)

# set the label text
label = 'Line Label'
label2 = 'Line Label2'

# set the font size for the label
fontsize = 12

# create the label and place it on top of the line
plt.text(x, y2 + 0.5, label, ha='center', va='bottom', fontsize=fontsize)
plt.text(x2, y2 + 0.7, label2, ha='center', va='bottom', fontsize=fontsize)

# show the plot
plt.show()