import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Wedge

# initial map
map = Basemap(projection='robin', resolution='l', lon_0=0)

plt.figure(figsize=(12, 6))

map.drawcoastlines(linewidth=0.1)
map.drawcountries(linewidth=0.1) 
map.drawmapboundary(fill_color='none')
map.fillcontinents(color='#DCDDDD')

parallels = np.arange(-90., 91., 30.)
meridians = np.arange(-180., 181., 60.)
map.drawparallels(parallels, labels=[1, 0, 0, 0],dashes= [4,1], linewidth=0.8) 
map.drawmeridians(meridians, labels=[0, 0, 0, 1],dashes= [4,1],linewidth=0.8) 

# data
with open('D:/onedrive/Lake_virome/VT/piemap.txt', 'r', encoding='utf-8') as file:
    for line in file:

        data = line.strip().split('\t')
        print(data)
        if len(data) < 7:
            continue

        latitude, longitude = float(data[0]), float(data[1])
        pie_data = [float(value) for value in data[2:]]

        x, y = map(longitude, latitude)

        colors = ['red', 'green', 'blue', 'yellow',"purple"]

        radius = 800000
        start_angle = 0

        for i, value in enumerate(pie_data):
            end_angle = start_angle + (360 * value / 100)
            wedge = Wedge((x, y), radius, start_angle, end_angle, fc=colors[i], zorder=2)
            plt.gca().add_patch(wedge)
            start_angle = end_angle

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
plt.show()