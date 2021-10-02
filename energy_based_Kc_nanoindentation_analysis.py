# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 16:50:15 2021

@author: Jacob Hempel
"""
import numpy as np
import math
import os
from scipy.integrate import trapz
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
#Load_unload_clipping function description:
# A function that takes in the load-displacement data set
# and outputs four data arrays: load_x, load_y, unload_x, and unload_y
# for a specified column. 
# IMPORTANT: assumes that all the columns have the same number of 
# rows. Also the function deletes the last 27 data points of the unloading curve.
#This number is arbitrary, but it's meant to get rid of the thermal drift segment.
def load_unload_clipping(data, sens,column): 
    
    load_x = []
    load_y = []
    unload_x =[]
    unload_y = []
    
    for i in range(len(data)-sens):
        if data[i+sens, column]-data[i,column] >= 0 and i < len(data)-50:
            load_x.append(data[i,column])
            load_y.append(data[i,column+1])
        if data[i+sens,column+1]-data[i,column+1] < 0 and data[i+sens, column]-data[i,column] <= 0 and i<len(data)-15: #Here is the -27 data points.
            unload_x.append(data[i,column])
            unload_y.append(data[i,column+1])
            
    return load_x, load_y, unload_x, unload_y

#A function designed to delete the first N number
#of rows from a data file. This is used to get rid of
#the three rows of header information "Specimen #, load, unload... etc#
def load_curve_fit_function(load_x, load_y, min_y_percentage, max_y_percentage):
    flag_l = False
    clipped_load_x = []
    clipped_load_y = []
    fit_load_x = []
    fit_load_y = []
    max_y = max(load_y)*max_y_percentage
    min_y = max(load_y)*min_y_percentage
    i = 0
    while load_y[i] < max_y:
        if load_y[i]> min_y:
            clipped_load_x.append(load_x[i])
            fit_load_x.append(load_x[i])
            clipped_load_y.append(load_y[i])
        i += 1
        if i >= len(load_y):
            break
    try:
        load_parameters, load_covariance = curve_fit(load_curve,clipped_load_x, clipped_load_y)
        k = load_parameters[0]
        l = load_parameters[1]
    except RuntimeError:
        print("Failed to fit to the loading curve.")
        flag_l = True
        k = 0
        l=0
    for i in range(len(fit_load_x)):
        fit_load_y.append(k*fit_load_x[i]**l)
    
    return fit_load_x, fit_load_y, k, l, flag_l

def unload_curve_fit_function(unload_x, unload_y, min_y_percentage, max_y_percentage):
    flag_u = False
    clipped_unload_x = []
    clipped_unload_y = []
    fit_unload_x = []
    fit_unload_y = []
    max_y = max(unload_y)*max_y_percentage
    min_y = max(unload_y)*min_y_percentage
    i = 0
    while unload_y[i] > min_y:
        if unload_y[i]< max_y:
            clipped_unload_x.append(unload_x[i])
            fit_unload_x.append(unload_x[i])
            clipped_unload_y.append(unload_y[i])
        i += 1
        if i >= len(unload_y):
            break
    try:
        unload_parameters, unload_covariance = curve_fit(unload_curve, clipped_unload_x, clipped_unload_y)
        hp = unload_parameters[0]
        g = unload_parameters[1]
        m = unload_parameters[2]
    except RuntimeError:
        #print("Failed to fit to the unloading curve.")
        flag_u = True
        hp=0
        g=0
        m=0
    for i in range(len(fit_unload_x)):
        fit_unload_y.append(g*(fit_unload_x[i]-hp)**m)
    return fit_unload_x, fit_unload_y, hp, g,m, flag_u

def del_top_row(data, newfilename, N):
    data_file = open(data, "r")
    py_data = data_file.readlines()
    for i in range(3):
        py_data = py_data[:-1]
    data_file.close()
    output = open(newfilename, "w+")
    for i in range(N):
        del py_data[0]

    for line in py_data:
        output.write(line)
    output.close()

    return None

    return None
#A fit function for the unload curve from ref. [1]
def unload_curve2(h,hp,hm,g,m):
    F = g*(((h-hp)/(hm-hp))**m) 
    return F
def unload_curve(h,hp,g,m):
    F = g*(np.abs((h-hp))**m) #The reason for the absolute value is 
    return F

def area_function(m0, m1,m2,m3, h):
    A = m0*h**2+m1*h+m2*h**(0.5)+m3*h**(0.25)
    return A
def load_curve(h,k,l):
    F = k*(h**l)
    return F
###########
#PLEASE be sure to read the readme before using this program and/or watch the youtube tutorial video. 
#This program assumes you have a .txt file with a specific format. There are no guarantees that this
#program will work if your data file is different!!!
#Assumed filename: material_load(units of mN)_measurement#_data.txt
    #Example: LiF_100mN_01_data.txt
material_list = ["fsilica_"] #Insert the material names here in quotes. Include an underscore. Example: "LiF_"
E_measured_list = [70434] #Insert the measured modulus values here IN UNITS OF MPa
H_measured_list = [8703] #Values of hardness IN UNITS OF MPa
loads = [100]

nu = [0.25] #List of poisson ratio for the materials in material_list.
measurement_num = ["01"] #Measurement number for this data set listed here. The first listed value is the measurmeent number
#for the first load data set in the loads list.
ignore_specimens_lists = [[]] #Suppose you find an outlier in the data that is caused by poor fitting. To ignore that
#specific data point, you need a list of lists, where each ignore_specimens[i] represents the specimen number in
#the load[i] that you want to ignore. Example: you have data for loads [100,200,300]. Your ignore_specimens
#variable needs to be re-written as [[],[],[]] since there are 3 load values. Let's say you want to ignore 
#specimen # 4 in the 200 load data set. Then you would have ignore_speciments[[],[4],[]].

m0 = 24.5 #These are your measured area function calibration constants. 
m1 = 0 #The default value is the ideal value for a Berkovich tip.
m2 = 0
m3 = 0
#show_plots will make it so that the python script will plot the load-unloading curves for all of your
#data, + their fits + m, l, and the calculated value of Kc for each curve. This will substantially increase
#the time it takes for the program to finish. Additionally, this will produce a PDF file in the same folder
#that this python script is in showing all these results, named material_Kc_analysis.pdf.
show_plots = False 
#write_kc displays kc in the console. I recommend starting your analysis by making show_plots false and 
#then write_kc as true so that you can first quickly identify any outliers, then see what the fits look like.
write_kc = True
#set ideal_case = true if you want to assume m=2 and l = 1.36 for your loading-unloading curves. Not recommended.
ideal_case = False
#If you would like to explicitly print out the irreversible work to see what the value is, set print_wirrev_fit too be true.
print_wirrev_fit = False
#print_load helps you keep trak of seeing which Kc values were calculated from which load data sets, if you have several load
#data sets.
print_load = True
#Select which material in the list material_list that you want to analyze. material_num = 0 is the first listed material.
material_num = 0

#Below, you define that range of the data you would like to try to fit to. This value should be a fraction of the min and max
#values of the loading or unloading curve. For example, max_fit_load = 0.1 and max_fit_load 0.9 means you'll try fitting the
#data from 10% to 90% of the loading data. This DOES include any holding period data as well.
min_fit_load= 0.1
max_fit_load= 0.9
min_fit_unload= 0.1
max_fit_unload = 0.9

if show_plots == True:
    pdf_filename = ''
    pdf_filename += material_list[material_num] 
    pdf_filename += "Kc_analysis"
    if ideal_case == True:
        pdf_filename+= "_idea.pdf"
    elif ideal_case == False:
        pdf_filename += ".pdf"
    pdf_output = PdfPages(pdf_filename)    
material = material_list[material_num]
global_hm = [] #global variables used to take averages over several load data sets.
global_hm_error = []
global_Am = []
global_Am_error=[]
global_Fm = []
global_Ucrack = []
global_Kc = []
global_Kc_error = []
global_Ucrack_error = []
l_vs_h_error = []
nu_material = nu[material_num]
for l in range(len(loads)):
    ignore_specimens = ignore_specimens_lists[l-1]
    load = loads[l]
    nu_measured = nu[material_num]
    filename = material_list[material_num]+str(load)+"mN_"+measurement_num[l]+"_data"
    if print_load == True:
        print(load)
    data_file_name = filename+".txt" 
    edited_file_name = filename+"_edited.txt"
    del_top_row(data_file_name, edited_file_name, 3) #This deletes the top 3 rows of data. If this is not needed, set 3 = 0 and it will basically do nothing.
    data = np.loadtxt(edited_file_name)
    os.remove(edited_file_name)
    E_measured= E_measured_list[material_num]
    H_measured = H_measured_list[material_num]
    ####Data lists####
    W_total_data = []
    W_rev_data = []
    W_irrev_data = []
    W_irrev_ratio_data = []
    Ucrack_estimate_data = []
    Kc_data = [] 
    hm_data = []
    hp_data = []
    Fm_data = []
    Ac_data = []
    k_data = []
    g_data = []
    m_data = []
    l_data = []
    f_column = len(data[1])
    
    for column in range(0,f_column,2): #This loop will go through all the columns of the excel sheet, plot their P vs h curves + try to fit to unloading curve, and calculate work ratio.
        specimen = int(column/2)
        if specimen in ignore_specimens:
            continue
        load_x, load_y, unload_x, unload_y = load_unload_clipping(data, 3,column) #Data is split into load portion, unload portion.
        fit_load_x, fit_load_y, k, l, flag_l = load_curve_fit_function(load_x, load_y, min_fit_load, max_fit_load)
        fit_unload_x, fit_unload_y, hp, g, m, flag_u = unload_curve_fit_function(unload_x, unload_y, min_fit_unload,max_fit_unload)
        
        if flag_l == False and flag_u == False:
            W_total = trapz(load_y, load_x) #Trapezoidal method of integrating load and unload curves.
            W_rev = trapz(unload_y, unload_x)
            W_irrev = W_total-W_rev
            hm = max(load_x)
            
            #print("m = ",m)
            #print("l = ", l)
            W_total_data.append(W_total)
            W_rev_data.append(W_rev)
            W_irrev_data.append(W_irrev)
            W_irrev_ratio_data.append(W_irrev/W_total)
            Ac_data.append(area_function(m0, m1,m2,m3, hm)*1e-18)
            hm_data.append(hm)
            Fm_data.append(max(load_y))
            hp_data.append(hp)
            k_data.append(k)
            l_data.append(l)
            m_data.append(m)
            g_data.append(g)
            unload_x_fit = np.arange(hp, hm, 0.1)
            unload_y_fit = unload_curve(unload_x_fit, hp, g, m)
            load_x_fit = np.arange(0, hm/1.06, 0.1)
            load_y_fit = load_curve(load_x_fit, k, l)
            W_rev_fit = trapz(unload_y_fit, unload_x_fit)
            W_irrev_fit = W_total - W_rev_fit
            if print_wirrev_fit == True:
                print(W_irrev_fit)
            if ideal_case == False:
                Upp_ratio_theory = (hp/hm)*((l+1)/(m+1))-((l-m)/(m+1)) #See ref. [2]
            if ideal_case == True:
                Upp_ratio_theory = (1-((1-3*(hp/hm)**2+2*(hp/hm)**3)/(1-(hp/hm)**2))) #See ref. [2]
            Ucrack = W_total*(W_irrev_fit/W_total-Upp_ratio_theory)*(1e-12)
            Ucrack_estimate_data.append(Ucrack)
            Gc = (Ucrack/(area_function(m0, m1,m2,m3, hm)*(1e-18)))*(1e-6)
            Kc = np.sqrt(Gc*E_measured/(1-nu_measured**2)) #See ref. [1] to know where these equations come from.
            Kc_data.append(Kc)
            #normalized_Ucrack_estimate_data.append((E_measured*(1e6)/(1-nu**2))*Ucrack/((ref_Kc[x]*1e6)**2))
            #phi.append((E_measured/(1-nu**2))*Ucrack/(ref_Kc[x]*1e6)**2)
#             Plotting the data and displaying work on plot
            if show_plots == True:
                plt.scatter(unload_x, unload_y, color='b', label='Unload data')
                plt.scatter(load_x, load_y, color='r', label='Load data')
                plt.xlim([0, max(load_x)+max(load_x)/5])
                plt.ylim([0, max(load_y)+max(load_y)/5])
                plt.ylabel("Load on sample (mN)")
                plt.xlabel("Displacement into surface (nm)")
                plt.legend()
                plt.text(max(load_x)/2, max(load_y),"Specimen # "+str(int(column/2+1)))
                l_text = "l = "
                m_text = "m = "
                kc_text = "Kc = "
                plt.text(max(load_x)/3.2,max(load_y)/1.3,kc_text+str(round(Kc, 3)))
                plt.text(max(load_x)/3.2,max(load_y)/1.5,m_text+str(round(m, 3)))
                plt.text(max(load_x)/3.2,max(load_y)/2,l_text+str(round(l, 3)))
                plt.plot(unload_x_fit, unload_y_fit, label='Unload fit', color='0.8')
                plt.plot(load_x_fit, load_y_fit, label='Load fit', color='y')
                plt.savefig(pdf_output, format='pdf')
                plt.show()
                plt.close()
            if write_kc == True:
                print("Specimen number "+str(specimen)+" Kc is:"+str((Kc)))
    global_hm.append(np.average(hm_data))
    global_hm_error.append(np.std(hm_data))
    global_Ucrack.append(np.average(Ucrack_estimate_data))
    global_Ucrack_error.append(np.std(Ucrack_estimate_data))
    global_Fm.append(np.average(Fm_data))
    global_Am.append(np.average(Ac_data))
    global_Am_error.append(np.std(Ac_data))
    hm_hf_ratio = []
    Kc_data_no_nan = []
    if len(Kc_data) == 0:
        continue
    for i in range(len(Kc_data)):
        if math.isnan(Kc_data[i]) == False:
            Kc_data_no_nan.append(Kc_data[i])
    global_Kc.append(np.average(Kc_data_no_nan))
    global_Kc_error.append(np.std(Kc_data_no_nan))
print("Fracture toughness for ",material_list[material_num]," is K_c = ", np.average(global_Kc), "+/-", np.average(global_Kc_error))
if show_plots == True:
    pdf_output.close()


#REFERENCES:
#[1] Taha, M. Reda, et al. "Fracture toughness of hydrated cement paste using nanoindentation." Fracture Mechanics of Concrete and Concrete Structures-Recent Advances in Fracture Mechanics of Concrete (2010). 
#[2]Cheng, Yang-Tse, Zhiyong Li, and Che-Min Cheng. "Scaling relationships for indentation measurements." Philosophical Magazine A 82.10 (2002): 1821-1829.
#[3] Liu, Ming, Dongyang Hou, and Chenghui Gao. "Berkovich nanoindentation of Zr55Cu30Al10Ni5 bulk metallic glass 
#at a constant loading rate." Journal of Non-Crystalline Solids 561 (2021): 120750.