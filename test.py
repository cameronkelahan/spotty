import numpy as np
import sys

surface_map = np.load("./spotty/ProxCen/flatMap.npy")
print("Loaded")
print(surface_map)
print("Type of surface map = ", type(surface_map))
print("Type of surface_map[0] = ", type(surface_map[0]))
print("Type of surface_map[0][0] = ", type(surface_map[0][0]))
print("Size of surfacemap[0][0] = ", sys.getsizeof(surface_map[0][0]))
print("Change")
rowcount = 0
surface_map = surface_map.astype(np.int8)
# for row in surface_map:
#     #print(rowcount)
#     surface_map[rowcount] = surface_map[rowcount].astype(np.int8)
#     # count=0
#     # for value in row:
#     #     surface_map[rowcount][count] = surface_map[rowcount][count].astype(np.int8)
#     #     count+= 1
#     rowcount += 1
print(surface_map)
print("Type of surface_map = ", type(surface_map))
print("Type of surface_map[0] = ", type(surface_map[0]))
print("Type of surface_map[0][0] = ", type(surface_map[0][0]))
print("Done")
np.save('./spotty/ProxCen/flatMap.npy', surface_map)