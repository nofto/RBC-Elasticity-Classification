import random
import espressomd
#import geometry
#from geometry import *

import object_in_fluid as oif
from espressomd import lb
from espressomd import lbboundaries
from espressomd import shapes
from espressomd.shapes import Rhomboid
from espressomd.shapes import Cylinder
from espressomd.interactions import OifLocalForces
from espressomd.interactions import OifGlobalForces
from espressomd.interactions import OifOutDirection
from espressomd.interactions import SoftSphereInteraction
from espressomd.interactions import MembraneCollisionInteraction

import numpy as np
import os, sys, time, os

# from oif_classes import OifCellType
# from oif_classes import OifCell
# from oif_classes import Edge
# from oif_classes import Triangle
# from oif_classes import PartPoint
# from oif_classes import FixedPoint
# from oif_classes import Mesh

import shutil
# from oif_utils import *
import time
import datetime


#simulation of rbcs in direct channel without obstacles with two different set of elastic coefficients
#####################################

#PROCEDURES AND FUNCTIONS
def FillBoundaries():
    #walls
    bottom = shapes.Rhomboid(corner=[0, 0, 0], a=[boxX, 0, 0], b=[0, boxY, 0], c=[0, 0, 1], direction = 1)
    oif.output_vtk_rhomboid(rhom_shape=bottom, out_file=vtk_directory+"/bottom.vtk")

    top = shapes.Rhomboid(corner=[0, 0, boxZ-1], a=[boxX, 0, 0], b=[0, boxY, 0], c=[0, 0, 1], direction = 1)
    oif.output_vtk_rhomboid(rhom_shape=top, out_file=vtk_directory + "/top.vtk")

    # left = shapes.Rhomboid(corner=[0, boxY-1, 1], a=[boxX, 0, 0], b=[0, 1, 0], c=[0, 0, boxZ-2], direction = 1)
    # oif.output_vtk_rhomboid(rhom_shape=left, out_file=vtk_directory + "/left.vtk")
    #
    # right = shapes.Rhomboid(corner=[0, 0, 1], a=[boxX, 0, 0], b=[0, 1, 0], c=[0, 0, boxZ-2], direction = 1)
    # oif.output_vtk_rhomboid(rhom_shape=right, out_file=vtk_directory + "/right.vtk")

    # obstacles = shapes.Rhomboid(corner=[rhom_x0, 0, 0], a=[rhom_x, 0, 0], b=[0, rhom_y, 0], c=[0, 0, rhom_z], direction = 1)
    # oif.output_vtk_rhomboid(rhom_shape=obstacles, out_file=vtk_directory + "/obstacles.vtk")

    # obstacles20 = shapes.Rhomboid(corner=[40, 0, 0], a=[20, 0, 0], b=[0, 20, 0], c=[0, 0, 20], direction = 1)
    # oif.output_vtk_rhomboid(rhom_shape=obstacles20, out_file=vtk_directory + "/obstacles20.vtk")

    # obstacles10 = shapes.Rhomboid(corner=[40, 0, 0], a=[20, 0, 0], b=[0, 20, 0], c=[0, 0, 10], direction = 1)
    # oif.output_vtk_rhomboid(rhom_shape=obstacles10, out_file=vtk_directory + "/obstacles10.vtk")

    # cylinder1 = shapes.Cylinder(center=[0, 0, boxZ/2], axis=[0.0, 0.0, 1.0], length=boxZ, radius = obst_radius, direction=1)
    # oif.output_vtk_cylinder(cyl_shape=cylinder1, n = 50, out_file=vtk_directory + "/cylinder1.vtk")

    # cylinders
    for id, pos in enumerate(obst_centers):
        obstacles = shapes.Cylinder(center=pos, axis=[0.0, 0.0, 1.0], length=boxZ, radius=obst_radius, direction=1)
        oif.output_vtk_cylinder(cyl_shape=obstacles, n=50, out_file=vtk_directory + "/cylinder" + str(id) + ".vtk")

        # stara verzia:
        # boundaries.append(shapes.Cylinder(center=pos, axis=[0.0, 0.0, 1.0], length=boxZ, radius=obst_radius, direction=1))
        # output_vtk_cylinder(center=pos, axis=[0.0, 0.0, 1.0], length=boxZ/2, radius=obst_radius, n=50, out_file=vtk_directory + "/cylinder" + str(id) + ".vtk")

        boundaries.append(obstacles)

    boundaries.append(bottom)
    boundaries.append(top)
    # boundaries.append(left)
    # boundaries.append(right)
    # boundaries.append(cylinder1)
    # boundaries.append(cylinder2)
    # boundaries.append(cylinder3)
    # boundaries.append(cylinder4)
    # boundaries.append(cylinder5)
    # boundaries.append(obstacles)

def DistanceBetween(x1,x2,x3,y1,y2,y3):
    return np.sqrt((x1 - y1)*(x1 - y1) + (x2 - y2)*(x2 - y2) + (x3 - y3)*(x3 - y3))

# check whether point is outside all obstacles
def IsOutsideObst (x,y):
    for i in obst_centers:
        ox = i[0]
        oy = i[1]
        dist = DistanceBetween(ox, oy, 0, x, y, 0)
        if dist < obst_radius+radius_rbc:
            return 0
        else:
            return 1

# find distance to the closest obstacle
def DistToClosestObst (x,y):
    min_dist = 100000000.0
    for i in obst_centers:
        ox = i[0]
        oy = i[1]
        dist = DistanceBetween(ox, oy, 0, x, y, 0)
        if dist<min_dist:
            min_dist = dist
    if min_dist < obst_radius:
        min_dist = 0
    else:
        min_dist = min_dist - obst_radius
    return min_dist

#find the smallest distance between the cells and the closest obstacle
def MinDistanceCellObstacle():
    dist_cell_obst = list()
    min_dist_cell_obst = list()
    for id, rbc in enumerate(rbcs):
        cell_obst_0 = 1000
        for j in obst_centers:
            x_cell = rbcs[id].get_origin()
            y_cell = rbcs[id].get_origin()
            x_cell = x_cell.tolist()
            y_cell = y_cell.tolist()
            x_cell = x_cell[0]
            y_cell = y_cell[1]
            x_obst = j[0]
            y_obst = j[1]
            cell_obst = DistanceBetween(x_cell,y_cell,0,x_obst,y_obst,0)
            if cell_obst_0 > cell_obst:
                cell_obst_0 = cell_obst
        min_dist_cell_obst.append(cell_obst_0)
    return min_dist_cell_obst

#INPUT PARAMETERS FROM COMMAND LINE
if len(sys.argv)!= 12:
    print ("Error: Incorrect number of input arguments.")
    print (" ")
else:
    sim_name = sys.argv[1]
    type_1_cell = int(sys.argv[2])
    type_2_cell = int(sys.argv[3])
    type_3_cell = int(sys.argv[4])
    type_4_cell = int(sys.argv[5])
    type_5_cell = int(sys.argv[6])
    type_6_cell = int(sys.argv[7])
    type_7_cell = int(sys.argv[8])
    type_8_cell = int(sys.argv[9])
    type_9_cell = int(sys.argv[10])
    type_10_cell = int(sys.argv[11])
    # num_cells = int(sys.argv[2])
    # diffCells = int(sys.argv[3])
    # ks_rbc = sys.argv[3]
    # kb_rbc = sys.argv[4]
    # kal_rbc = sys.argv[5]
    # kag_rbc = sys.argv[6]
    # kv_rbc = sys.argv[7]
    # fluid_force = sys.argv[8]
    # allowed_overlap = sys.argv[9]

# print("input: " + str(sys.argv))

num_cells = 0
for n in range(2,len(sys.argv)):
    num_cells += int(sys.argv[n])
print("pocet buniek: " + str(num_cells))

cell_types = list()
j=0
for i in range(2,len(sys.argv)):
    for k in range(int(sys.argv[i])):
        cell_types.append(i-1)
        j = j + 1
# print("cell_types: " + str(cell_types))


#OUTPUT DIRECTORY
output_dir = "output_NxKSin1_Aa"
# output_dir = "output_1xKSin1_Aa"
# output_dir = "output_2xKSin1_D"
# output_dir = "output_2xKSin1_test"

if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
    # print("bol vytvoreny adresar \"output\" ")

directory = output_dir + "/sim" + str(sim_name)
if not os.path.isdir(directory):
    os.mkdir(directory)
# else:
    # print("simulacia " + str(sim_name) + " uz existuje!" )

seeding_dir = directory + "/seeding"
if not os.path.isdir(seeding_dir):
    os.mkdir(seeding_dir)


vtk_directory = directory+"/vtk"
os.mkdir(vtk_directory)

os.mkdir(directory+"/vtk_failed")

shutil.copy(__file__, directory + "/script_sim"+sim_name+".py") #copy actual simulation script

#GEOMETRY

# boxX = 60.0
# boxY = 40.0
# boxZ = 40.0

# boxX = 260.0
# boxY = 150.0
# boxZ = 100.0

boxX = 104.0
boxY = 60.0
boxZ = 40.0

# rhom_x0 = 200
# rhom_x = 30
# rhom_y = boxY
# rhom_z = 20

obst_radius = 20
obst_centers = list()
obst_centers.append([0, 0, boxZ/2])
obst_centers.append([0, boxY, boxZ/2])
obst_centers.append([boxX/2, boxY/2, boxZ/2])
obst_centers.append([boxX, 0, boxZ/2])
obst_centers.append([boxX, boxY, boxZ/2])

#CELL PARAMETERS
# nnode_rbc = 141

# ks_rbc = 0.007
# kb_rbc = 0.00025
# kal_rbc = 0.001
# kag_rbc = 0.9
# kv_rbc = 0.5

# #CELL PARAMETERS
nnode_rbc = 374
# diffCells = 5

# ks_rbc = float(ks_rbc)
# kb_rbc = float(kb_rbc)
# kal_rbc = float(kal_rbc)
# kag_rbc = float(kag_rbc)
# kv_rbc = float(kv_rbc)

kag_rbc = 0.7
kv_rbc = 0.9

ks_rbc_1 = 0.005
kb_rbc_1 = 0.003
kal_rbc_1 = 0.02

ks_rbc_2 = 0.0133
kb_rbc_2 = 0.01
kal_rbc_2 = 0.01

ks_rbc_3 = 0.0216
kb_rbc_3 = 0.01
kal_rbc_3 = 0.01

ks_rbc_4 = 0.03
kb_rbc_4 = 0.01
kal_rbc_4 = 0.01

ks_rbc_5 = 1.0
kb_rbc_5 = 1.0
kal_rbc_5 = 1.0

ks_rbc_6 = 1.0
kb_rbc_6 = 1.0
kal_rbc_6 = 1.0

ks_rbc_7 = 1.0
kb_rbc_7 = 1.0
kal_rbc_7 = 1.0

ks_rbc_8 = 1.0
kb_rbc_8 = 1.0
kal_rbc_8 = 1.0

ks_rbc_9 = 1.0
kb_rbc_9 = 1.0
kal_rbc_9 = 1.0

# ks_rbc_5 = 0.05
# kb_rbc_5 = 0.01
# kal_rbc_5 = 0.01
#
# ks_rbc_6 = 0.1
# kb_rbc_6 = 0.01
# kal_rbc_6 = 0.01
#
# ks_rbc_7 = 0.15
# kb_rbc_7 = 0.01
# kal_rbc_7 = 0.01
#
# ks_rbc_8 = 0.225
# kb_rbc_8 = 0.01
# kal_rbc_8 = 0.01
#
# ks_rbc_9 = 0.3
# kb_rbc_9 = 0.01
# kal_rbc_9 = 0.01

ks_rbc_10 = 1.0
kb_rbc_10 = 1.0
kal_rbc_10 = 1.0

kvisc_rbc = 0.00
radius_rbc = 3.91
resize_rbc = float(radius_rbc)
mass_rbc = 500.0

#INPUT FILES
rbc_nodes = "input/rbc" + str(nnode_rbc) + "nodes.dat"
rbc_triangles = "input/rbc" + str(nnode_rbc) + "triangles.dat"

#LBFLUID PARAMETERS

density = 1.025
viscosity = 1.3
# fluid_force = float(fluid_force)
fluid_force = 0.001
time_step = 0.05
# time_step = 0.1
lb_time_step = time_step
lb_grid = 1

surface = 4*np.pi*4**2
surface_rbc = (np.sqrt(134.3)/3.91*radius_rbc)**2
friction_rbc =(393.0/nnode_rbc)*np.sqrt(surface_rbc/surface)*((5.6-1.82)/(6.0/1.025-1.5)*(viscosity-1.5)+(10-1.82)/(6-1.025)*(density-1.025)+1.82)
friction = friction_rbc

print ("fluid_force:" + str(fluid_force))
print ("friction:" + str(friction_rbc))
print ("mass:" + str(mass_rbc/nnode_rbc))

#SOFT SPHERE PARAMETERS
soft_a = 0.001
soft_n = 1.2
soft_cut = 0.1
soft_offset = 0.0
particle_type = num_cells + 1

print("soft_a:" + str(soft_a))
print("soft_n:" + str(soft_n))
print("soft_cut:" + str(soft_cut))

#MEMBRANE COLLISION PARAMETERS
membrane_a = 0.01
membrane_n = 1.0
membrane_cut = 0.4
membrane_offset = 0.0

print("membrane_a:" + str(membrane_a))
print("membrane_n:" + str(membrane_n))
print("membrane_cut:" + str(membrane_cut))

# SELF-CELL PARAMETERS
self_a = 0.001
self_n = 1.2
self_cut = 0.5

#ITERATION PARAMETERS
cycle = 0
n_steps = n_steps_vtk = n_steps_rbc_info = n_steps_time_section = 1000
# n_steps = 500
# n_steps_vtk = 500 #how ofted output vtk
# n_steps_rbc_info = 500 #how ofted output rbc info
# n_steps_time_section = 500 #how ofted output info about all cells
sum_cycle = 5000000000

# sim_info.dat
sim_info = open(directory + "/sim_info.dat", "w+")
sim_info.write("start: " + str(datetime.datetime.now()) + "\n")
sim_info.write("simulacia: " + str(sim_name) + "\n")
sim_info.write("\n")
sim_info.write("pocet vsetkych buniek: " + str(num_cells) + "\n")
sim_info.write("pocet buniek typu 1: " + str(sys.argv[2]) + "\n")
sim_info.write("pocet buniek typu 2: " + str(sys.argv[3]) + "\n")
sim_info.write("pocet buniek typu 3: " + str(sys.argv[4]) + "\n")
sim_info.write("pocet buniek typu 4: " + str(sys.argv[5]) + "\n")
sim_info.write("pocet buniek typu 5: " + str(sys.argv[6]) + "\n")
sim_info.write("pocet buniek typu 6: " + str(sys.argv[7]) + "\n")
sim_info.write("pocet buniek typu 7: " + str(sys.argv[8]) + "\n")
sim_info.write("pocet buniek typu 8: " + str(sys.argv[9]) + "\n")
sim_info.write("pocet buniek typu 9: " + str(sys.argv[10]) + "\n")
sim_info.write("pocet buniek typu 10: " + str(sys.argv[11]) + "\n")
sim_info.write("\n")
sim_info.write("geometria: x = " + str(boxX) + ", y = " + str(boxY) + ", z = " + str(boxZ) + "\n")
# sim_info.write("prekazka: x0 = " + str(rhom_x0) + ", x = " + str(rhom_x) + ", y = " + str(rhom_y) + ", z = " + str(rhom_z) + "\n")
sim_info.write("\n")
sim_info.write("fluid_force: " + str(fluid_force) + "\n")
sim_info.write("friction: " + str(friction) + "\n")
sim_info.write("time_step: " + str(time_step) + "\n")
sim_info.write("lb_grid: " + str(lb_grid) + "\n")
sim_info.write("n_steps: " + str(n_steps) + "\n")
sim_info.write("density: " + str(density) + "\n")
sim_info.write("viscosity: " + str(viscosity) + "\n")
sim_info.write("mass_rbc: " + str(mass_rbc) + "\n")
sim_info.write("\n")
sim_info.write("MEMBRANE COLLISION PARAMETERS: a = " + str(membrane_a) + ", n = " + str(membrane_n) + ", cut = " + str(membrane_cut) + "\n")
sim_info.write("SOFT SPHERE PARAMETERS: a = " + str(soft_a) + ", n = " + str(soft_n) + ", cut = " + str(soft_cut) + "\n")
sim_info.write("SELF-CELL INTERACTIONST (SOFT SPHERE PARAMETERS): a = " + str(self_a) + ", n = " + str(self_n) + ", cut = " + str(self_cut) + "\n")
sim_info.write("\n")
sim_info.write("nnode_rbc: " + str(nnode_rbc) + "\n")
sim_info.write("radius_rbc: " + str(radius_rbc) + "\n")
# sim_info.write("elasticke koef.: " + "ks_rbc: " + str(ks_rbc) + ", kb_rbc: " + str(kb_rbc) + ", kal_rbc: " + str(kal_rbc) + ", kag_rbc: " + str(kag_rbc) + ", kv_rbc: " + str(kv_rbc) + "\n")
sim_info.write("bunky typu 1: " + "ks_rbc: " + str(ks_rbc_1) + ", kb_rbc: " + str(kb_rbc_1) + ", kal_rbc: " + str(kal_rbc_1) + ", kag_rbc: " + str(kag_rbc) + ", kv_rbc: " + str(kv_rbc) + "\n")
sim_info.write("bunky typu 2: " + "ks_rbc: " + str(ks_rbc_2) + ", kb_rbc: " + str(kb_rbc_2) + ", kal_rbc: " + str(kal_rbc_2) + ", kag_rbc: " + str(kag_rbc) + ", kv_rbc: " + str(kv_rbc) + "\n")
sim_info.write("bunky typu 3: " + "ks_rbc: " + str(ks_rbc_3) + ", kb_rbc: " + str(kb_rbc_3) + ", kal_rbc: " + str(kal_rbc_3) + ", kag_rbc: " + str(kag_rbc) + ", kv_rbc: " + str(kv_rbc) + "\n")
sim_info.write("bunky typu 4: " + "ks_rbc: " + str(ks_rbc_4) + ", kb_rbc: " + str(kb_rbc_4) + ", kal_rbc: " + str(kal_rbc_4) + ", kag_rbc: " + str(kag_rbc) + ", kv_rbc: " + str(kv_rbc) + "\n")
sim_info.write("bunky typu 5: " + "ks_rbc: " + str(ks_rbc_5) + ", kb_rbc: " + str(kb_rbc_5) + ", kal_rbc: " + str(kal_rbc_5) + ", kag_rbc: " + str(kag_rbc) + ", kv_rbc: " + str(kv_rbc) + "\n")
sim_info.write("bunky typu 6: " + "ks_rbc: " + str(ks_rbc_6) + ", kb_rbc: " + str(kb_rbc_6) + ", kal_rbc: " + str(kal_rbc_6) + ", kag_rbc: " + str(kag_rbc) + ", kv_rbc: " + str(kv_rbc) + "\n")
sim_info.write("bunky typu 7: " + "ks_rbc: " + str(ks_rbc_7) + ", kb_rbc: " + str(kb_rbc_7) + ", kal_rbc: " + str(kal_rbc_7) + ", kag_rbc: " + str(kag_rbc) + ", kv_rbc: " + str(kv_rbc) + "\n")
sim_info.write("bunky typu 8: " + "ks_rbc: " + str(ks_rbc_8) + ", kb_rbc: " + str(kb_rbc_8) + ", kal_rbc: " + str(kal_rbc_8) + ", kag_rbc: " + str(kag_rbc) + ", kv_rbc: " + str(kv_rbc) + "\n")
sim_info.write("bunky typu 9: " + "ks_rbc: " + str(ks_rbc_9) + ", kb_rbc: " + str(kb_rbc_9) + ", kal_rbc: " + str(kal_rbc_9) + ", kag_rbc: " + str(kag_rbc) + ", kv_rbc: " + str(kv_rbc) + "\n")
sim_info.write("bunky typu 10: " + "ks_rbc: " + str(ks_rbc_10) + ", kb_rbc: " + str(kb_rbc_10) + ", kal_rbc: " + str(kal_rbc_10) + ", kag_rbc: " + str(kag_rbc) + ", kv_rbc: " + str(kv_rbc) + "\n")
# sim_info.write("allowed_overlap: " + str(allowed_overlap) + "\n")
sim_info.write("\n")
sim_info.close()

#####################################
#SIMULATION INITIALIZATION
#####################################
system = espressomd.System(box_l = [boxX, boxY, boxZ])      # novsia verzia espressa - spustanie na clustry
# system = espressomd.System()                                # starsia verzia espressa - spustanie na pc
# system.box_l = [boxX, boxY, boxZ]                          # starsia verzia espressa - spustanie na pc
system.time_step = time_step
system.cell_system.skin = 0.3
# system.cell_system.skin = 0.2

# CLUSTER PARAMETERS parall
system.cell_system.node_grid = [2,6,4]
system.min_global_cut = 5.


# ######### start NO-RANDOM SEEDING
# cell_positions = list()
# cell_positions = [[20.0,30.0,10.0],[50.0,25.0,10.0],[80.0,27.0,10.0]]
# #########

######### start NO-RANDOM SEEDING - pre seeding z nodes.dat (zapise sa tolko centier so suradnicami [0.0, 0.0, 0.0] kolko je zadanych buniek)
# cell_positions = list()
# for k in range(0,num_cells,1):
#     cell_positions.append([0.0, 0.0, 0.0])
######### end

######### start NO-RANDOM SEEDING - seeding retrieved from seeding.dat and rotate.dat files
# cell_positions = list()
# in_file = open("input_seeding/seeding.dat", "r").read().split("\n")
# for line in in_file:  # extracts coordinates from the string line
#     line = list([float(x) for x in line.split()])
#     cell_positions.append(line)
# print ("cell_positions: " + str(cell_positions))
#
# cell_rotate = list()
# in_file = open("input_seeding/rotate.dat", "r").read().split("\n")
# for line in in_file:  # extracts coordinates from the string line
#     line = list([float(x) for x in line.split()])
#     cell_rotate.append(line)
# print ("cell_rotate: " + str(cell_rotate))
#
# for k in range(0,num_cells,1):
#     sim_info = open(directory + "/sim_info.dat", "a")
#     ox = cell_positions[k][0]
#     oy = cell_positions[k][1]
#     oz = cell_positions[k][2]
#     sim_info.write("seeding of cell: "+str(k)+" with origin: "+str(ox)+", "+str(oy)+", "+str(oz) + " time: " + str(datetime.datetime.now()) + " " + "\n")
# sim_info.close()
######## end

########## start NO-RANDOM SEEDING - seeding retrieved from txt file
# cell_positions = open("output_sim2000_adjustment/seeding2503.txt").readlines()
#
# def f(l):
#     l=l.replace("[", "").replace("]", "").replace("'", "").replace("\n", "").split(",")
#     return (float(l[0]),float(l[1]),float(l[2]))
#
# cell_positions = list(map(f,cell_positions))
# # print (cell_positions)
#
# for k in range(0,num_cells,1):
#     sim_info = open(directory + "/sim_info.dat", "a")
#     ox = cell_positions[k][0]
#     oy = cell_positions[k][1]
#     oz = cell_positions[k][2]
#     sim_info.write("seeding of cell: "+str(k)+" with origin: "+str(ox)+", "+str(oy)+", "+str(oz) + " time: " + str(datetime.datetime.now()) + " " + "\n")
# sim_info.close()
######### end


# ######### START - RANDOM SEEDING OF CELLS
cell_positions = list()
vacant_layer_x = radius_rbc     # dlzka na zaciatku kanala bez naseedovanych buniek
# seedBoxX = 100
# seedBoxX = boxX
seedBoxX = boxX - radius_rbc    # naseedovanie len do urcitej casti kanala (vzhladom na os x)
# allowed_overlap = float(allowed_overlap)
allowed_overlap = 0.0

for k in range(0,num_cells,1):
    origin_ok = 0
    sim_info = open(directory + "/sim_info.dat", "a")
    while origin_ok !=1:
        #generate random position in channel;
        # ox = random.random() * seedBoxX + radius_rbc - 1      # hodnoty: 2,9 - 102,9
        # ox = random.random() * seedBoxX                         # hodnoty: 0 - 100
        ox = random.random() * (seedBoxX - vacant_layer_x) + vacant_layer_x      # hodnoty: 0 az dlzka boxu minus vacant_layer na zaciatku kanala
        # oy = random.random() * boxY + radius_rbc
        # oy = random.random() * (boxY-4-2*radius_rbc) + 2 + radius_rbc
        # oz = random.random() * (boxZ-4-2*radius_rbc) + 2 + radius_rbc
        oy = random.random() * (boxY-2-2*radius_rbc) + 1 + radius_rbc
        oz = random.random() * (boxZ-2-2*radius_rbc) + 1 + radius_rbc
        origin_ok = 1

        #check whether are outside of cylinders obstacles
        if DistToClosestObst(ox, oy) < radius_rbc + 1:
            origin_ok = 0

        #check whether are outside of Rhomboid obstacles
        # if origin_ok == 1:
        #     for i in range(0, num_cells,1):
        #         if ((ox > (rhom_x0 - radius_rbc)) and (ox < (rhom_x0 + rhom_x + radius_rbc))):
        #             if oz < (rhom_z + radius_rbc):
        #                 origin_ok = 0

        # check that it does not collide with other rbc
        if origin_ok == 1:
            for i in range(0,k,1):
                dist = DistanceBetween(ox, oy, oz, cell_positions[i][0], cell_positions[i][1], cell_positions[i][2])
                if dist < (1.0 - allowed_overlap)*2*radius_rbc + 1.0:
                    origin_ok = 0

        # if everything was ok, remember origin
        if origin_ok == 1:
            print ("seeding of cell: "+str(k)+" with origin: "+str(ox)+", "+str(oy)+", "+str(oz))
            sim_info.write("seeding of cell: "+str(k)+" with origin: "+str(ox)+", "+str(oy)+", "+str(oz) + " time: " + str(datetime.datetime.now()) + " " + "\n")
            cell_positions.append([ox,oy,oz])
    sim_info.close()
########## END - RANDOM SEEDING OF CELLS

# CREATING TEMPLATES FOR CELLS
rbc_type_1 = oif.OifCellType(nodes_file=rbc_nodes, triangles_file=rbc_triangles, system = system, ks=ks_rbc_1, kb=kb_rbc_1, kal=kal_rbc_1, kag=kag_rbc, kv=kv_rbc, kvisc=kvisc_rbc\
, resize =[resize_rbc, resize_rbc, resize_rbc], check_orientation=False, normal= True)

rbc_type_2 = oif.OifCellType(nodes_file=rbc_nodes, triangles_file=rbc_triangles, system = system, ks=ks_rbc_2, kb=kb_rbc_2, kal=kal_rbc_2, kag=kag_rbc, kv=kv_rbc, kvisc=kvisc_rbc\
, resize =[resize_rbc, resize_rbc, resize_rbc], check_orientation=False, normal= True)

rbc_type_3 = oif.OifCellType(nodes_file=rbc_nodes, triangles_file=rbc_triangles, system = system, ks=ks_rbc_3, kb=kb_rbc_3, kal=kal_rbc_3, kag=kag_rbc, kv=kv_rbc, kvisc=kvisc_rbc\
, resize =[resize_rbc, resize_rbc, resize_rbc], check_orientation=False, normal= True)

rbc_type_4 = oif.OifCellType(nodes_file=rbc_nodes, triangles_file=rbc_triangles, system = system, ks=ks_rbc_4, kb=kb_rbc_4, kal=kal_rbc_4, kag=kag_rbc, kv=kv_rbc, kvisc=kvisc_rbc\
, resize =[resize_rbc, resize_rbc, resize_rbc], check_orientation=False, normal= True)

rbc_type_5 = oif.OifCellType(nodes_file=rbc_nodes, triangles_file=rbc_triangles, system = system, ks=ks_rbc_5, kb=kb_rbc_5, kal=kal_rbc_5, kag=kag_rbc, kv=kv_rbc, kvisc=kvisc_rbc\
, resize =[resize_rbc, resize_rbc, resize_rbc], check_orientation=False, normal= True)

rbc_type_6 = oif.OifCellType(nodes_file=rbc_nodes, triangles_file=rbc_triangles, system = system, ks=ks_rbc_6, kb=kb_rbc_6, kal=kal_rbc_6, kag=kag_rbc, kv=kv_rbc, kvisc=kvisc_rbc\
, resize =[resize_rbc, resize_rbc, resize_rbc], check_orientation=False, normal= True)

rbc_type_7 = oif.OifCellType(nodes_file=rbc_nodes, triangles_file=rbc_triangles, system = system, ks=ks_rbc_7, kb=kb_rbc_7, kal=kal_rbc_7, kag=kag_rbc, kv=kv_rbc, kvisc=kvisc_rbc\
, resize =[resize_rbc, resize_rbc, resize_rbc], check_orientation=False, normal= True)

rbc_type_8 = oif.OifCellType(nodes_file=rbc_nodes, triangles_file=rbc_triangles, system = system, ks=ks_rbc_8, kb=kb_rbc_8, kal=kal_rbc_8, kag=kag_rbc, kv=kv_rbc, kvisc=kvisc_rbc\
, resize =[resize_rbc, resize_rbc, resize_rbc], check_orientation=False, normal= True)

rbc_type_9 = oif.OifCellType(nodes_file=rbc_nodes, triangles_file=rbc_triangles, system = system, ks=ks_rbc_9, kb=kb_rbc_9, kal=kal_rbc_9, kag=kag_rbc, kv=kv_rbc, kvisc=kvisc_rbc\
, resize =[resize_rbc, resize_rbc, resize_rbc], check_orientation=False, normal= True)

rbc_type_10 = oif.OifCellType(nodes_file=rbc_nodes, triangles_file=rbc_triangles, system = system, ks=ks_rbc_10, kb=kb_rbc_10, kal=kal_rbc_10, kag=kag_rbc, kv=kv_rbc, kvisc=kvisc_rbc\
, resize =[resize_rbc, resize_rbc, resize_rbc], check_orientation=False, normal= True)

rbc_types = [rbc_type_1, rbc_type_2, rbc_type_3, rbc_type_4, rbc_type_5, rbc_type_6, rbc_type_7, rbc_type_8, rbc_type_9, rbc_type_10]

cell_types2 = list()
for i in range(2,len(sys.argv)):
    for k in range(int(sys.argv[i])):
        cell_types2.append(rbc_types[i-2])

# CREATING CELLS
rbcs = list()

##Pre nahodny seeding
for id, pos in enumerate(cell_positions):
    type = id + 1
    print("id cell = " + str(id) + " -> type_" + str(cell_types[id]))
    sim_info = open(directory + "/sim_info.dat", "a")
    sim_info.write("cell " + str(id) + " -> type_" + str(cell_types[id]) + " | ")
    sim_info.close()
    cell_rbc = oif.OifCell(cell_type=cell_types2[id], particle_type=id, origin=pos, rotate=[random.random()*2*np.pi, random.random()*2*np.pi, random.random()*2*np.pi], particle_mass=mass_rbc/nnode_rbc)
    # cell_rbc = OifCell(cell_type=rbc_type, part_type=id, origin=pos, rotate=[0.0, 0.0, 0.0], part_mass=mass_rbc/nnode_rbc)
    rbcs.append(cell_rbc)


# ## Pre seeding z nodes.dat
# for id, pos in enumerate(cell_positions):
#     type = id + 1
#     print("id cell = " + str(id) + " -> type_" + str(cell_types[id]))
#     sim_info = open(directory + "/sim_info.dat", "a")
#     sim_info.write("cell " + str(id) + " -> type_" + str(cell_types[id]) + " | ")
#     sim_info.close()
#     cell_rbc = oif.OifCell(cell_type=cell_types2[id], particle_type=id, origin=pos, particle_mass=mass_rbc/nnode_rbc)
#     cell_rbc.set_mesh_points("seeding/seeding_hct10_02/nodes" + str(id) + ".dat")
#     rbcs.append(cell_rbc)

sim_info = open(directory + "/sim_info.dat", "a")
sim_info.write("\n\n")
sim_info.write("all cells were created: " + str(datetime.datetime.now()) + "\n")
sim_info.close()

################################################################# start - save seeding of cells to dat file
#  save seeding of centers to dat file
out_file_seeding = open(seeding_dir + "/seeding.dat", "a")
for id in range(0,num_cells,1):
    rbc_position = rbcs[id].get_origin()
    out_file_seeding.write(str(rbc_position[0]) + " " + str(rbc_position[1]) + " " + str(rbc_position[2]) + "\n")
out_file_seeding.close()

# save seeding of each  of particles to dat file
particles_cell = list()
for id, rbc in enumerate(rbcs):
    for particle in range(0,nnode_rbc,1):
        particles_cell = rbcs[id].mesh.points[particle].get_pos()
        out_file_nodes = open(seeding_dir + "/nodes" + str(id) + ".dat", "a")
        out_file_nodes.write(str(particles_cell[0]) + " " + str(particles_cell[1]) + " " + str(particles_cell[2]) + "\n")
    out_file_nodes.close()
################################################################# end - save seeding of cells to dat file

#volume = rbcs[1].mesh.volume()
#print("objem bunky: " + str(volume))

# vypocet friction zo skriptu
# for i in range(0, num_cells):
#     aa = rbcs[i].suggest_LBgamma(viscosity, density)
#     print ("RBC friction - cell " + str(i) + ": " + str(aa))

# aa=rbcs[0].suggest_LBgamma(viscosity, density)
# print ("rbc0 friction:" + str(aa))
#
# aaa=rbcs[1].suggest_LBgamma(viscosity, density)
# print ("rbc1 friction:" + str(aaa))

#CREATING BOUNDARIES
boundaries = list()
FillBoundaries()
for boundary in boundaries:
    system.constraints.add(shape = boundary, particle_type = particle_type)
print ("particle boundaries were added:")
sim_info = open(directory + "/sim_info.dat", "a")
sim_info.write("particle boundaries were added: " + str(datetime.datetime.now()) + "\n")
sim_info.close()

#VTK OF INITIAL CELLS POSITIONS
for id, rbc in enumerate(rbcs):
    rbc.output_vtk_pos_folded(vtk_directory+"/rbc"+str(id)+".vtk")

#CELL-WALL INTERACTIONS
for i in range(num_cells):
    system.non_bonded_inter[i,particle_type].soft_sphere.set_params(a = soft_a, n = soft_n, cutoff = soft_cut, offset = soft_offset)

print ("cell-wall interactions were added")
sim_info = open(directory + "/sim_info.dat", "a")
sim_info.write("cell-wall interactions were added: " + str(datetime.datetime.now()) + "\n")
sim_info.write("cell-cell interactions ... " + "\n")
sim_info.close()

#CELL-CELL INTERACTIONS
ncells = len(rbcs)
i = 0
print("cell-cell interactions ...")
while (i < ncells):
    j = i+1
    while (j < ncells):
        # system.non_bonded_inter[i,j].membrane_collision.set_params(membrane_a = membrane_a, membrane_n = membrane_n, membrane_cut = membrane_cut, membrane_offset = membrane_offset)    # PC - starsia verzia
        system.non_bonded_inter[i,j].membrane_collision.set_params(a = membrane_a, n = membrane_n, cutoff = membrane_cut, offset = membrane_offset)                                       # CL - novsia verzia
        j = j+1
    # sim_info = open(directory + "/sim_info.dat", "a")
    # sim_info.write(str(i) + ",x cell-cell interactions were added " + str(datetime.datetime.now()) + "\n")
    # sim_info.close()
    i = i+1

print ("cell-cell interactions were added")
sim_info = open(directory + "/sim_info.dat", "a")
sim_info.write("all cell-cell interactions were added: " + str(datetime.datetime.now()) + "\n")
sim_info.close()

# SELF-CELL INTERACTIONS
for i in range(0, len(rbcs)):
    system.non_bonded_inter[i,i].soft_sphere.set_params(a = self_a, n = self_n, cutoff = self_cut, offset = 0.0)
    # print ("self-cell:" + str(i) + "-" + str(i))

print ("self-cell interactions were added:")
sim_info = open(directory + "/sim_info.dat", "a")
sim_info.write("all self-cell interactions were added: " + str(datetime.datetime.now()) + "\n")
sim_info.close()

#LOADING FLUID
# lbf = espressomd.lb.LBFluid(agrid=lb_grid, dens=density, visc=viscosity, tau=lb_time_step, fric=friction, ext_force_density=[fluid_force,0.0,0.0])
# system.actors.add(lbf)

lb_params = {'agrid': lb_grid, 'dens': density, 'visc': viscosity, 'tau': time_step, 'ext_force_density': [fluid_force, 0.0, 0.0]}
lbf = espressomd.lb.LBFluid(**lb_params)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=friction)

print ("fluid was added")
sim_info = open(directory + "/sim_info.dat", "a")
sim_info.write("fluid was added: " + str(datetime.datetime.now()) + "\n")
sim_info.close()

#LOAD FLUID BOUDARIES
for boundary in boundaries:
    system.lbboundaries.add(lbboundaries.LBBoundary(shape = boundary))
print ("fluid boundaries were added:")
sim_info = open(directory + "/sim_info.dat", "a")
sim_info.write("fluid boundaries were added: " + str(datetime.datetime.now()) + "\n")
sim_info.close()

#VTK OF INITIAL CELLS POSITIONS
for id, rbc in enumerate(rbcs):
    rbc.output_vtk_pos_folded(vtk_directory+"/rbc"+str(id)+".vtk")
print ("initial vtk files were created:")

#OUTPUT FILES FOR EACH CELL CREATED
# for id, rbc in enumerate(rbcs):
#     out_file_rbc = open(directory + "/rbc" + str(id) + "_sim" + str(sim_name) + ".dat", "w+")
#     out_file_rbc.write("")
# out_file_rbc.close()

out_file_rbc = open(directory + "/rbc_description_sim" + str(sim_name) + ".dat", "w+")
out_file_rbc.write("1-cycle 2-3-4-rbc_center_position[x y z] 5-6-7-rbc_velocity[x y z]  rbc_cuboid[8-9-10-x_min[x y z] 11-12-13-x_max[x y z] 14-15-16-y_min[x y z] \
17-18-19-y_max[x y z] 20-21-22-z_min[x y z] 23-24-25-z_max[x y z]] 26-27-28-x_min_vel[x y z] 29-30-31-x_max_vel[x y z] 32-33-34-y_min_vel[x y z] 35-36-37-y_max_vel[x y z] \
38-39-40-z_min_vel[x y z] 41-42-43-z_max_vel[x y z] 44-volume 45-surface \n")
out_file_rbc.close()

# out_file_section = open(directory + "/cycle_description_sim" + str(sim_name) + ".dat", "w+")
# out_file_section.write("cell_id cell_center_pos[x y z] x_min x_max y_min y_max z_min z_max \n ")
# out_file_section.close()


###################################################################
                       #MAIN LOOP
###################################################################
while cycle < sum_cycle:
    print (str(cycle))
    sim_info = open(directory + "/sim_info.dat", "a")
    sim_info.write(str(cycle) + ": " + str(datetime.datetime.now()) + "\n")

    # START - output vtk, calculate and output volumetric flow rate and velocity (vtkfluid.vtk)
    if cycle % n_steps_vtk == 0:
        # calculate and output fluid velocity
        lbf.print_vtk_velocity(str(vtk_directory + "fluid.vtk"))
        # START - calculate and output volumetric flow rate
        out_file_fluid = open(directory + "/volFlowRate" + str(sim_name) + ".dat", "a")
        volRate = 0.0
        for i in range(0, int(boxY), 1):
            for j in range(0, int(boxZ), 1):
                volRate += lbf[int(boxX / 2), i, j].velocity[0]
        out_file_fluid.write(str(cycle) + ' ' + str(volRate) + '\n')
        # END - calculate and output fluid velocity
        for id, rbc in enumerate(rbcs):
            #output vtk for rbc
            rbc.output_vtk_pos_folded(vtk_directory+"/rbc"+str(id)+"_"+str(cycle)+".vtk")
    # END - calculate and output volumetric flow rate

    # START H-data, output info about each cell
    if cycle % n_steps_rbc_info == 0:
        for id, rbc in enumerate(rbcs):
            out_file_rbc = open(directory + "/rbc" + str(id) + "_sim" + str(sim_name) + ".dat", "a")
            rbc_position = rbcs[id].get_origin()
            rbc_velocity = rbcs[id].get_velocity()
            bound_points = rbcs[id].point_bound()
            x_min_vel = bound_points[0].get_vel()
            x_max_vel = bound_points[1].get_vel()
            y_min_vel = bound_points[2].get_vel()
            y_max_vel = bound_points[3].get_vel()
            z_min_vel = bound_points[4].get_vel()
            z_max_vel = bound_points[5].get_vel()
            x_min = bound_points[0].get_pos()
            x_max = bound_points[1].get_pos()
            y_min = bound_points[2].get_pos()
            y_max = bound_points[3].get_pos()
            z_min = bound_points[4].get_pos()
            z_max = bound_points[5].get_pos()
            rbc_volume = rbcs[id].volume()
            rbc_surface = rbcs[id].surface()

            out_file_rbc.write(str(cycle) + " " + str(rbc_position[0]) + " " + str(rbc_position[1]) + " " + str(rbc_position[2]) + " " \
            + str(rbc_velocity[0]) + " " + str(rbc_velocity[1]) + " " + str(rbc_velocity[2]) + " " \
            + str(x_min[0]) + " " + str(x_min[1]) + " " + str(x_min[2]) + " " + str(x_max[0]) + " " + str(x_max[1]) + " " + str(x_max[2]) + " " \
            + str(y_min[0]) + " " + str(y_min[1]) + " " + str(y_min[2]) + " " + str(y_max[0]) + " " + str(y_max[1]) + " " + str(y_max[2]) + " " \
            + str(z_min[0]) + " " + str(z_min[1]) + " " + str(z_min[2]) + " " + str(z_max[0]) + " " + str(z_max[1]) + " " + str(z_max[2]) + " " \
            + str(x_min_vel[0]) + " " + str(x_min_vel[1]) + " " + str(x_min_vel[2]) + " " + str(x_max_vel[0]) + " " + str(x_max_vel[1]) + " " + str(x_max_vel[2]) + " " \
            + str(y_min_vel[0]) + " " + str(y_min_vel[1]) + " " + str(y_min_vel[2]) + " " + str(y_max_vel[0]) + " " + str(y_max_vel[1]) + " " + str(y_max_vel[2]) + " " \
            + str(z_min_vel[0]) + " " + str(z_min_vel[1]) + " " + str(z_min_vel[2]) + " " + str(z_max_vel[0]) + " " + str(z_max_vel[1]) + " " + str(z_max_vel[2]) + " " \
            + str(rbc_volume) + " " + str(rbc_surface) + " " + "\n")
            out_file_rbc.close()
    # END H-data

    # # START H-data, output info in time section
    # if cycle % n_steps_time_section == 0:
    #     out_file_section = open(directory + "/cycle" + str(cycle) + "_data_sim" + str(sim_name) + ".dat", "a")
    #     # out_file_section.write("cell_id cell_center_pos[x y z] x_min x_max y_min y_max z_min z_max \n ")
    #     for id, rbc in enumerate(rbcs):
    #         bound_points = rbcs[id].point_bound()
    #         x_min = bound_points[0].get_pos()
    #         x_max = bound_points[1].get_pos()
    #         y_min = bound_points[2].get_pos()
    #         y_max = bound_points[3].get_pos()
    #         z_min = bound_points[4].get_pos()
    #         z_max = bound_points[5].get_pos()
    #         rbc_approx_pos = rbcs[id].get_approx_origin()
    #         out_file_section.write(str(id) + " " + str(rbc_approx_pos[0]) + " " +str(rbc_approx_pos[1]) + " " + str(rbc_approx_pos[2]) + " " + str(x_min[0]) + " " + str(x_max[0]) + " " + str(y_min[1]) + " " + str(y_max[1]) + " " + str(z_min[2]) + " " + str(z_max[2])+ "\n")
    # # END H-data

    # system.integrator.run(n_steps)

    try:
        system.integrator.run(n_steps)
    except Exception:
        for id, rbc in enumerate(rbcs):
            #output vtk for rbc
            rbc.output_vtk_pos_folded(directory+"/vtk_failed/failed_rbc"+str(id)+"_"+str(cycle)+".vtk")
        raise

    cycle = cycle+n_steps
    sim_info.close()
exit()
