#!/usr/bin/env python

##################################################################################################################################################
#                                                                                                                                                #
#                This script gets all opt structues and energies and compares them with each other to erase dublicates!                          #
#                               The Gaussian output files after optimization should be _SMD_opt.log.                                             #
#               The qh(Thrular) free energies are taken from "output.dat" with struture names ending with "_gv.log"!                             #
#               The Dataframe "all_xyz" contains all geometries of all .log files in the folder for geometric analysis!                          #
#               The dataframe "data_opt" contains a list of all E_tot_B3LYP_SMD and qh-G energies sorted by qh-G!.                               #
#                                                                                                                                                #
#                                                       Version: 1.06 (05.2022)                                                                  #
#                                               Done: Distance! String to Float mistake corrected!                                               #
#                                               To Do:                                                                                           #
#                                                                                                                                                #
#                                                                                                                                                #
##################################################################################################################################################

import pandas as pd
import os
import glob
import sys
import re
import fileinput
import shutil           # for copying files
import numpy as np

pd.options.display.float_format = '{:10.6f}'.format

##############################Some Varibles#####################################################################################################

criteria = float(0.026255)               # kJ/mol
criteria_hartree = float(0.00001)        # energy criteria in Hartree 0.0000099999992
hartconv = float(2625.498)               # conversion factor of [Hartree] to [kJ/mol]
R = float(0.008314511)                   # gas constant as [kJ/K*mol]
temp = float(289.15)                     # standard condition temperature
kcalvonv = float(4184)                   # conversion factor [kcal] to [kJ/mol]


elemDict = {"1" : "H", "2" : "He", "3" : "Li", "4" : "Be", "5" : "B", \
"6"  : "C", "7"  : "N", "8"  : "O",  "9" : "F", "10" : "Ne", \
"11" : "Na" , "12" : "Mg" , "13" : "Al" , "14" : "Si" , "15" : "P", \
"16" : "S"  , "17" : "Cl" , "18" : "Ar" , "19" : "K"  , "20" : "Ca", \
"21" : "Sc" , "22" : "Ti" , "23" : "V"  , "24" : "Cr" , "25" : "Mn", \
"26" : "Fe" , "27" : "Co" , "28" : "Ni" , "29" : "Cu" , "30" : "Zn", \
"31" : "Ga" , "32" : "Ge" , "33" : "As" , "34" : "Se" , "35" : "Br", \
"36" : "Kr" , "37" : "Rb" , "38" : "Sr" , "39" : "Y"  , "40" : "Zr", \
"41" : "Nb" , "42" : "Mo" , "43" : "Tc" , "44" : "Ru" , "45" : "Rh", \
"46" : "Pd" , "47" : "Ag" , "48" : "Cd" , "49" : "In" , "50" : "Sn", \
"51" : "Sb" , "52" : "Te" , "53" : "I"  , "54" : "Xe" , "55" : "Cs", \
"56" : "Ba" , "57" : "La" , "58" : "Ce" , "59" : "Pr" , "60" : "Nd", \
"61" : "Pm" , "62" : "Sm" , "63" : "Eu" , "64" : "Gd" , "65" : "Tb", \
"66" : "Dy" , "67" : "Ho" , "68" : "Er" , "69" : "Tm" , "70" : "Yb", \
"71" : "Lu" , "72" : "Hf" , "73" : "Ta" , "74" : "W"  , "75" : "Re", \
"76" : "Os" , "77" : "Ir" , "78" : "Pt" , "79" : "Au" , "80" : "Hg", \
"81" : "Tl" , "82" : "Pb" , "83" : "Bi" , "84" : "Po" , "85" : "At", \
"86" : "Rn" , "87" : "Fr" , "88" : "Ra" , "89" : "Ac" , "90" : "Th", \
"91" : "Pa" , "92" : "U"  , "93" : "Np" , "94" : "Pu" , "95" : "Am", \
"96" : "Cm" , "97" : "Bk" , "98" : "Cf" , "99" : "Es" ,"100" : "Fm", \
"101": "Md" ,"102" : "No" ,"103" : "Lr" ,"104" : "Rf" ,"105" : "Db", \
"106": "Sg" ,"107" : "Bh" ,"108" : "Hs" ,"109" : "Mt" ,"110" : "Ds", \
"111": "Rg" ,"112" : "Uub","113" : "Uut","114" : "Uuq","115" : "Uup", \
"116": "Uuh","117" : "Uus","118" : "Uuo"}


def getCoordinates(dataList) :


        listofstart = []				# a list of indices of 'Standard orientation' in dataList
        listofstop = []					# a list of indices of '--------------------' in dataList
        for index, value in enumerate(dataList):
            if 'Standard orientation' in value:
                #print(index, value)
                index = int(index)
                listofstart.append(index)
            if '----------' in value:
                index = int(index)
                listofstop.append(index)
        start = max(listofstart) + 5			# last occurence of 'Standard orientation' and  xyz starts 5 lines after 'Standard orientation'
        stoptop = start - 1				# get position of  "------" at top
        indexStoptop = listofstop.index(stoptop)	# get the index in dataList of "------" at top
        indexStopbottom = indexStoptop + 1		# get the index in dataList of "------" at bottom
        end = listofstop[indexStopbottom]		# get index of position of "------" at bottom
        xyzlist = dataList[start:end]			# extract last coordinates as list from dataList
        list = pd.Series(xyzlist)				# convert to df
        xyz = list.str.split(expand=True)		# split strings by whitespace and expand as new columns
        xyz.columns = ['Center Number', 'Atomic Number', 'Atomic Type', 'X', 'Y', 'Z']
        atomlist = []
        for i in range(0, xyz.shape[0]):
              value = xyz.iloc[i]['Atomic Number']
              atom = str(elemDict.get(str(value)))	# convert 'Atomic Number' to 'Element'
              atomlist.append(atom)			# apend elements to list
        xyz['Element'] = atomlist			# new column on xyz
        xyz = xyz[['Center Number', 'Element', 'X', 'Y', 'Z']]		# reshape dataframe
        xyz = xyz.astype({"X": float, "Y": float, "Z": float})		# convert  X Y Z to float for later math. operations
        
        return xyz


def getElectronicEnergies(file) :		# input argument file in a loop


        f = open(file)					# Open file on read mode
        dataList = f.read().split("\n")			# Create a list containing all lines
        f.close()
        list_of_SCF_energies = []			# datalist of all indices of lines containing 'SCF Done'
        for index, value in enumerate(dataList):	# fill datalist
            if 'SCF Done:' in value:
                #print(index, value)
                index = int(index)
                list_of_SCF_energies.append(index)
                #print(list_of_SCF_energies)
        last_energy_index = max(list_of_SCF_energies)	# get last 'SCF Done' as highest number in datalist
        #print(last_energy_index)
        line_to_split = dataList[last_energy_index]	# get said line
        #print(line_to_split)
        E_tot  = line_to_split.split()[4]	# get energy value by splitting line by whitespace and returning 4th value
        #print(E_tot_B3LYP_SMD)



        return E_tot


############################Save Methode code as string#############################

directory_in_str = str(os.getcwd())
directory = os.fsencode(directory_in_str)

for file in os.listdir(directory):
    filename = os.fsdecode(file)

    if filename.endswith("_SMD_opt.log"):		#-----------------------------------------
        method_code = "_SMD_opt"			#
        file_end = "_SMD_opt.log"			#
    elif filename.endswith("_SMD_opt_freq.log"):	#
        method_code = "_SMD_opt_freq"			#
        file_end = "_SMD_opt_freq.log"			# my common filnames for opt and opt + freq
    elif filename.endswith("_smd_SP.log"):		#
        method_code = "_smd_SP"				#
        file_end = "_smd_SP.log"			#
    elif filename.endswith("_SP_freq.log"):		#
        method_code = "_SP_freq"			#
        file_end = "_SP_freq.log"			#-----------------------------------------

#    elif filename.endswith("_TZ.out"):			#-----------------------------------------
#        method_code_TZ = "_TZ"				#
#        file_end_TZ = "_TZ.log"				# my common filnames for DPLNO QZ and TZ
#    elif filename.endswith("_QZ.out"):			#
#        method_code_QZ = "_QZ"				#
#        file_end_QZ = "_QZ.out"			#-----------------------------------------#

#    elif filename.endswith("_SMD_freq.log"):		#-----------------------------------------
#        method_code_freq = "_SMD_freq"			# frequency calc ending (no opt)
#        file_end_freq = "_SMD_freq.log"			#-----------------------------------------

#    elif filename.endswith("_SP_gasph.log"):            #-----------------------------------------
#        method_code_SP = "_SP_gasph"                  # gasphase SP calc
#        file_end_SP = "_SP_gasph.log"                 #-----------------------------------------


    else:
        print("No specific file extension found in:", filename)

############################Getting XYZ for all .log files##########################


#directory_in_str = str(os.getcwd())
#directory = os.fsencode(directory_in_str)

all_xyz = pd.DataFrame(columns=['Structure', 'Center Number', 'Element','X','Y','Z'])
for file in os.listdir(directory):
     filename = os.fsdecode(file)
     if filename.endswith(file_end):
         basename  = str(os.path.splitext(filename)[0])
         #basename_xyz_df = basename + "_xyz_df"
         #newfile = str(filename) + ".xyz"          # Gives the name to a possible new file
         f = open(file)                              # Open file on read mode
         dataList = f.read().split("\n")             # Create a list containing all lines
         f.close()
         atom = getCoordinates(dataList)
         atom['Structure'] = basename
         atom = atom[['Structure', 'Center Number', 'Element','X','Y','Z']]
         all_xyz = all_xyz.append(atom)
         continue
     else:
         print("Out:",filename)
         continue

all_xyz = all_xyz.replace(method_code, '',regex=True)


print(all_xyz)

print("---------------------------------------------------------------")


###################################Acessing Structurefile###################################################
#for Structure, single_xyz in all_xyz.groupby('Structure'):
#    if Structure == 'cy2ne1w_kicked07':
#        print(single_xyz)
#        print(single_xyz.info())
#    else:
#        continue
##################################Extracting Total Energy from .log########################################

data_opt = pd.DataFrame(columns=["Structure","E_tot_B3LYP_SMD"])

for file in os.listdir(directory):
     filename = os.fsdecode(file)
     if any(filename.endswith(file_end) for file in os.listdir('.')):
         print("Optimized structures found!")
         #data_opt = pd.DataFrame(columns=["Structure","E_tot_B3LYP_SMD"])
         filename = os.fsdecode(file)
         basename  = str(os.path.splitext(filename)[0])
         E_tot = getElectronicEnergies(file)
         print("Out:", E_tot)
         data_opt = data_opt.append({'Structure': basename, 'E_tot_B3LYP_SMD': E_tot}, ignore_index=True)
         #data_opt = data_opt.append({'E_tot_B3LYP_SMD': E_tot}, ignore_index=True)


data_opt = data_opt.replace(method_code, '',regex=True)

#print(data_opt)



##############################################Add qh(Thrular)-G from output dat######################

data_goodvibes = pd.read_csv("Goodvibes_output.dat",skiprows=11 , delim_whitespace=True, names=["None", "Structure", "E", "ZPE", "H", "S", "T.qh-S", "G", "qh-G"],)


# delim_whitespace=True creates new column on any instance of blank space (tab or space) and is more flexible than sep="\t"
# skips the first 11 rows --> header of .dat file

data_goodvibes[["None", "Structure"]] = data_goodvibes[["None", "Structure"]].astype(str)             #define columns as string
data_goodvibes[["E", "ZPE", "H", "S", "T.qh-S", "G", "qh-G"]] = data_goodvibes[["E", "ZPE", "H", "S", "T.qh-S", "G", "qh-G"]].astype(float)   #define columns as float

data_goodvibes = data_goodvibes[["Structure", "qh-G", "ZPE", "H", "S", "T.qh-S", "G", "E"]]    # cut of first column " o o o o "
data_goodvibes = data_goodvibes[["Structure", "qh-G"]]
data_goodvibes = data_goodvibes.replace("_gv", '',regex=True)
print(data_goodvibes)


data_opt = pd.merge(data_opt, data_goodvibes, on="Structure")
data_opt = data_opt.round(decimals=7)


###################################Acessing Structurefile#############################

print("Before sorting E:", data_opt)

data_opt["E_tot_B3LYP_SMD"] = data_opt["E_tot_B3LYP_SMD"].astype(float)					# converting type: object to float, otherwiese sorting will fail

data_opt = data_opt.sort_values(by="E_tot_B3LYP_SMD")							# sort E
print("sorted by E:", data_opt)
#print(all_xyz)

all_centroids = pd.DataFrame(columns=['Structure', 'Center Number', 'Element','X','Y','Z'])

for i in range(0, data_opt.shape[0]):
    value = data_opt.iloc[i]['Structure']				# value structure name in data_opt
    #print(value)
    for Structure, single_xyz in all_xyz.groupby('Structure'):		# Structure: structure name in all_xyz
        #print(Structure)
        if Structure == value:
            print(single_xyz)
            print("------------------------------ Geometric Centroid----------------------------------------------")
            sum_X = float(single_xyz['X'].sum())			# Sum over whole "X" columns as float
            sum_Y = float(single_xyz['Y'].sum())
            sum_Z = float(single_xyz['Z'].sum())
            num_of_atoms = len(single_xyz.index)			# Atom count as lenght of dataframe "single_xyz"
            centroid = pd.DataFrame()
            centroid['X'] = abs((sum_X / num_of_atoms) - single_xyz['X'])	# Calculating geometric centroid
            centroid['Y'] = abs((sum_Y / num_of_atoms) - single_xyz['Y'])
            centroid['Z'] = abs((sum_Z / num_of_atoms) - single_xyz['Z'])
            centroid['Structure'] = single_xyz['Structure']			# Add "Structure" column to new dataframe "centroid"
            centroid['Center Number'] = single_xyz['Center Number']
            centroid['Element'] = single_xyz['Element']
            centroid = centroid[['Structure', 'Center Number', 'Element','X','Y','Z']]
            print(centroid)
            print("--------------------------------next Structure-------------------------------------------------")
            all_centroids = all_centroids.append(centroid)		#create dataframe with all centroid information by appending "centroid" to each other
        else:
            continue
print("------------------------------------------------------------------------------------------------")
print("------------------------------------------EOF---------------------------------------------------")
print("------------------------------------------------------------------------------------------------")
print(all_centroids)

############################The following block iterates troght the sorted data_opt list and also extracts the following structure

list_to_exclude = []		# list of structures to exlude from data_opt

data_opt_DeltaXYZ = pd.DataFrame(columns=["Structure","E_tot_B3LYP_SMD", "DeltaXYZ"])

for i in range(0, (data_opt.shape[0] - 1)):				# iterating over data_opt 
    #print(data_opt)
    value = data_opt.iloc[i]['Structure']                               # value structure name in data_opt
    #print(value)
    for Structure, single_centroid_1 in all_centroids.groupby('Structure'):          # Structure: structure name in all_centroids
        #print(Structure)
        if Structure == value:						# just extracting centroid information with same name, stored as single_centroid variable
            #print(single_centroid_1)
            break
    value2 = data_opt.iloc[(i + 1)]['Structure']			# get next centroid
    print("Structure:", value, "was compared to:", value2)
    for Structure, single_centroid_2 in all_centroids.groupby('Structure'):          # Structure: structure name in all_centroids
        #print(Structure)
        if Structure == value2:
            #print(single_centroid_2)
            break
    diff_eval = pd.DataFrame()
    diff_eval['X'] = abs(single_centroid_1['X'] - single_centroid_2['X'])	# absolute difference of X coordinate between two consecutive structures
    diff_eval['Y'] = abs(single_centroid_1['Y'] - single_centroid_2['Y'])
    diff_eval['Z'] = abs(single_centroid_1['Z'] - single_centroid_2['Z'])
    diff_eval['Structure'] = single_centroid_1['Structure']
    diff_eval = diff_eval[['Structure','X','Y','Z']]
    sum_diff_eval_X = float(diff_eval['X'].sum())				# sum of all differences in X axis
    sum_diff_eval_Y = float(diff_eval['Y'].sum())
    sum_diff_eval_Z = float(diff_eval['Z'].sum())
    print("Delta_X:", sum_diff_eval_X, "Delta_Y:", sum_diff_eval_Y, "Delta_Z:", sum_diff_eval_Z)
    Delta_XYZ = sum_diff_eval_X + sum_diff_eval_Y + sum_diff_eval_Z

    ##############Delta_XYZ######################
    Delta_XYZ = sum_diff_eval_X + sum_diff_eval_Y + sum_diff_eval_Z
    print("The sum of all differences is:", Delta_XYZ)
    energy = data_opt.iloc[(i + 1)]['E_tot_B3LYP_SMD']
    #print(value2, energy)
    single_Delta_XYZ = pd.DataFrame(columns=["Structure","E_tot_B3LYP_SMD", "DeltaXYZ"])
    single_Delta_XYZ = single_Delta_XYZ.append({'Structure': value2, 'E_tot_B3LYP_SMD': energy, 'DeltaXYZ': Delta_XYZ}, ignore_index=True)
    data_opt_DeltaXYZ = data_opt_DeltaXYZ.append(single_Delta_XYZ, ignore_index=True)



    if sum_diff_eval_X < 0.1 and sum_diff_eval_Y < 0.1 and sum_diff_eval_Z < 0.1:		# Setting criteria of 0.1 Angstroem to define as "Different Structure"
       print("******SAME Structure******")
       value2 = str(value2)
       list_to_exclude += [value2]
    else:
       print("******Different Structure******")
print(data_opt)
print("---------------------------------------------------------------------------")
data_opt['Structure'] = data_opt['Structure'].astype('str')
print(data_opt.info())          # convert  X Y Z to float for later math. operations
data_opt_geom_filt = pd.DataFrame(columns=['Structure', 'E_tot_B3LYP_SMD'])
data_opt_geom_filt  = data_opt[~data_opt['Structure'].isin(list_to_exclude)]
print(data_opt_geom_filt)

########################Filtering Data by Energy Criteria####################
#######################Not so elegant, but works!!###########################

count_column = data_opt_geom_filt.shape[0]   # count items in colums
count_column -= 1              # item count starts with zero

data_opt_geom_filt = data_opt_geom_filt.reset_index(drop=True)
data_opt_geom_filt["RowNumber"] = data_opt_geom_filt.reset_index().index    # crate new column with numbering from 0 --> Rownumber


print(data_opt_geom_filt)

list_to_add = [0]           #generating list of rows to add from old to new data
row = 0                     # start at 0
last_to_compare_defined = False
print("\n")
print("Energy criteria:", criteria_hartree, "Hartree")
print("\n"*2)
while row in data_opt_geom_filt["RowNumber"]:         #iterate over "RowNumber" column
   data_one = abs(data_opt_geom_filt.iloc[row]["E_tot_B3LYP_SMD"])   #store date of E
   row += 1
   data_two = abs(data_opt_geom_filt.iloc[row]["E_tot_B3LYP_SMD"])   #store data of next E
   if ((data_one - data_two) >= criteria_hartree):
        print("Criteria met! For:", row)
        #print(data_one - data_two)
        last_to_compare = row
        last_to_compare_defined = True          #this boolean is important to define the last_to_compare varibale even when the first two structures are no match
   if last_to_compare_defined == True:
       data_last_to_compare = abs(data_opt_geom_filt.iloc[last_to_compare]["E_tot_B3LYP_SMD"])
   else:
       data_last_to_compare = abs(data_opt_geom_filt.iloc[0]["E_tot_B3LYP_SMD"])
   if ((data_last_to_compare - data_two) < criteria_hartree):
       #print("Criteria not met (delete from list)! For:", row)
       #print(last_to_compare, data_last_to_compare)
       if last_to_compare_defined == False:
           list_to_add.append(0)
       else:
           list_to_add.append(last_to_compare)
   if (row == count_column):
       break

print("------------------------Code run smoothly!!------------------------")


list_to_add = list(dict.fromkeys(list_to_add))							# removing dublicates from rows to add to new column
data_opt_geom_filt_criteria = pd.DataFrame(columns=['Structure', 'E_tot_B3LYP_SMD'])
print(list_to_add)
data_opt_geom_filt_criteria = data_opt_geom_filt.iloc[list_to_add]              		# create new dataframe with rows that should be added
data_opt_geom_filt_criteria = data_opt_geom_filt_criteria[["Structure", "E_tot_B3LYP_SMD", "qh-G"]] 	# delete "RowNumber" helper column
data_opt_geom_filt_criteria = data_opt_geom_filt_criteria.reset_index(drop=True)		# resetting index values
data_opt_geom_filt_criteria = data_opt_geom_filt_criteria.sort_values(by=["E_tot_B3LYP_SMD"])

#######################Getting distance between water and acidic center###########

print("\n---------------------------------------------------------------------------")
print("Do you want to determine distance between two atoms in .log files? (y/n)")
dist_statement = input()

if dist_statement == "y" or dist_statement == "Y":
    print("Choose Atomic Centers to measure distances!")
    print("Center Number 1 (Acidic side):")
    cent_num_1 = int(input())
    cent_num_1 -= 1			# to get index
    print("\n")
    print("Center Number 2 (Water oxygen):")
    cent_num_2 = int(input())
    cent_num_2 -= 1			# to get index
    print("Measuring distance between", cent_num_1 + 1, "and", cent_num_2 + 1)

    data_dist = pd.DataFrame(columns=['Structure', 'Distance'])
    columns = list(data_dist)
    data = []

    for i in range(0, data_opt.shape[0]):
        value = data_opt.iloc[i]['Structure']                               # value structure name in data_opt
        #print(value)
        for Structure, single_xyz in all_xyz.groupby('Structure'):          # Structure: structure name in all_xyz
            #print(Structure)
            if Structure == value:
                print(single_xyz)
                #single_xyz["Center Number"] = single_xyz["Center Number"].astype(int)
                X_1 = float(single_xyz.at[cent_num_1, 'X'])
                Y_1 = float(single_xyz.at[cent_num_1, 'Y'])
                Z_1 = float(single_xyz.at[cent_num_1, 'Z'])
                X_2 = float(single_xyz.at[cent_num_2, 'X'])
                Y_2 = float(single_xyz.at[cent_num_2, 'Y'])
                Z_2 = float(single_xyz.at[cent_num_2, 'Z'])
                print("Coordinates of first:",X_1, Y_1, Z_1)
                print("Coordinates of second:",X_2, Y_2, Z_2)
                distance = np.sqrt(((X_2 - X_1)**2) + ((Y_2 - Y_1)**2) + ((Z_2 - Z_1)**2))
                print(distance, Structure)
                values = [Structure, distance]
                zipped = zip(columns, values)
                a_dictionary = dict(zipped)
                print(a_dictionary)
                data.append(a_dictionary)
                print(data)
                print("--------------------------------next Structure-------------------------------------------------")
        else:
            continue
    data_dist = data_dist.append(data, True)
    print(data_dist)

    data_opt_DeltaXYZ = pd.merge(data_opt_DeltaXYZ, data_dist, on="Structure")
    data_opt_geom_filt_criteria = pd.merge(data_opt_geom_filt_criteria, data_dist, on="Structure")
    print(data_opt_DeltaXYZ)
    #delta_xyz = data_opt_DeltaXYZ['Structure','DeltaXYZ']
    #delta_xyz = data_opt_DeltaXYZ[['Structure','DeltaXYZ']].copy()
    #data_opt_geom_filt_criteria = pd.merge(data_opt_geom_filt_criteria, delta_xyz, on="Structure")
    data_opt_geom_filt_criteria.to_csv(r'./getDistance.txt', sep='\t', header='true', index=False, index_label=False)
else:
    print("Proceed without determining distances!!......")


#################################################################################

print("\n")
print("\n")
print("------------------Geometric Centroid and Energy Criteria applied to dataset---------------------")
print("Before refinement:")
print(data_opt)
print("\n")
print("\n")
print("Result sorted by E_tot_B3LYP_SMD:")
print("\n")
print("\n")
print(data_opt_geom_filt_criteria)


######################Some Information###################################
print("\n")
print("\n")
print("------------------Geometric Centroid and Energy and Delta_XYZ---------------------")
print("----------------------------------------------------------------------------------")
print(data_opt_DeltaXYZ)
print("\n")
print("\n")
print("\n")
print("Datapoints before refinement (Num. of .log-files):", len(data_opt))
print("\n")
print("Datapoints after refinement:", len(data_opt_geom_filt_criteria))
print("\n")
print("Structures discarded by energy criteria of", criteria_hartree, "Hartree is:", (len(data_opt) - len(data_opt_geom_filt_criteria)))
print("\n")


####################Copying best 10 .logfiles to subfolder "best10_E"##############################################

print("---------------------------------------------------------------------------")
print("Do you want to copy the best 50 free energies (E) into new folder? (y/n)")
copy_statement = input()

if copy_statement == "y" or copy_statement == "Y":
    print("Copying files..........")
    os.mkdir('best50_E')
    try:
        for i in range(0, 50):                                                  # do  iterate 50 times
            name_of_file = data_opt_geom_filt_criteria.iloc[i]['Structure']     #
            name_of_file = name_of_file + file_end                        # reattach file ending to copy in new folder
            print("copy:",name_of_file, "to best50_E/")
            shutil.copy(name_of_file, 'best50_E')
    except IndexError:
        pass
else:
    print("No files copied!")
