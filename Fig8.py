import matplotlib.pyplot as plt
import random
import numpy as np
import pandas as pd
import math

toPrt = True

def to_numpy(file_name, sheet_name):
    """converts CDI data excel file to numpy array""" 
    df = pd.read_excel(file_name, sheet_name)
    df.dropna(inplace = True)
    numparr = df.to_numpy()
    if toPrt: print("Numpy array for ",file_name,"\n=======\n",numparr)
    return numparr

def substrate_blanks(numparr, *blank_rows):
    blank_mean = numparr[[*blank_rows],:].mean()
    if toPrt: print("Blank mean\n=======\n",blank_mean)
    numparr = numparr - blank_mean
    if toPrt: print("Numpy array minus blanks\n=======\n",numparr)
    return numparr

def std_curve(numparr, curve_start_col):
    """CDI standard curve coeff and intercept on y axis"""
    X = np.array([0.25, 0.5, 0.75, 1])
    y = np.array([numparr[4:7, curve_start_col+1],
                  numparr[4:7, curve_start_col],
                  numparr[1:4, curve_start_col+1],
                  numparr[1:4, curve_start_col]],
                 dtype = float)
    if toPrt: print("y vals for standard curve\n=======\n",y)
    result = np.polyfit(X,y,1)
    coeff = result[0].mean()
    intercept = result[1].mean()
    if toPrt: print("CDI value per ug PrP = ",str(coeff))
    if toPrt: print("Intercept on y axis =", str(intercept))
    return coeff, intercept

def relate_to_PrP(numparr, coeff, intercept, ul_analysed = 10):
    """convert CDI values to PrP per g brain, using std curve"""
    numparr = (numparr[1:7,0:10]-intercept)/coeff
    if toPrt: print("Values as ug PrP per 2.5mg brain sample\n========\n",numparr)
    numparr = numparr/(ul_analysed/10000)
    if toPrt: print("Values as ug PrP per gram brain \n========\n",numparr)
    return numparr

def get_vals(file_name, sheet_name, curve_col_start):
    """gets data and converts CDI vals to concentrations of PrP
    Returns as numpy array 
    """
    numparr = to_numpy(file_name, sheet_name)
    numparr = substrate_blanks(numparr, 0,7)
    coeff, intercept = std_curve(numparr, curve_col_start)
    numparr = relate_to_PrP(numparr, coeff, intercept)
    return numparr

def means_SDs(vals):
    vals = np.array(vals, dtype=np.float64)
    means = np.mean(vals, axis=0)
    SDs = np.std(vals, axis=0)
    return means, SDs

def plotBar(title, mean_vals, DNstd, PK_concs = ('0', '1', '2.5', '10', '50')):
    if toPrt: print("for plot entitled ",title)
    if toPrt: print("mean vals\n======\n",mean_vals)
    if toPrt: print("PK concs\n=======\n",PK_concs)
    N = len(mean_vals)
    ind = np.arange(N)    # the x locations for the groups
    width = 0.35       # the width of the bars: can also be len(x) sequence
    plt.bar(ind, mean_vals, width, yerr = DNstd)
    plt.ylabel(chr(956)+"g PrPSc /gram Brain")
    plt.xlabel("PK "+chr(956)+"g/ml")
    plt.title(title)
    plt.xticks(ind, PK_concs)      
    y_scale = 41
    if (max(mean_vals) + max(DNstd)) > 41:
        y_scale = max(mean_vals)+ max(DNstd)    
    plt.yticks(np.arange(0, y_scale, 10))

#get data into numpy arrays
#numparrPMCA = get_vals("CDI13 030B.xls", "Plate", curve_col_start = 10)
##numparrFFIoneP1 = get_vals("CDI13 004.xls", "Plate", curve_col_start = 9)
##numparrFFIoneP2 = get_vals("CDI13 004b.xls", "Plate", curve_col_start = 9)
numparrvCJDFC= get_vals("CDI 12 007 PLATE B.xls", "Plate", curve_col_start = 10)
##numparrvCJDthal= get_vals("CDI 12 016 PLATE B.xls", "Plate", curve_col_start = 9)
##numparrsCJDthal= get_vals("CDI13 004b.xls", "Plate", curve_col_start = 6)
##numparrsCJDMM2TFC= get_vals("CDI 12 008 PLATE B.xls", "Plate", curve_col_start = 10)
##numparrsCJDMM2CFC= get_vals("CDI 12 009 PLATE A.xls","Plate", curve_col_start = 10)

#calculate D-N vals
##DN_vals_PMCA = numparrPMCA[3:,0:5] - numparrPMCA[0:3,0:5]
##DN_vals_PMCA2 = numparrPMCA[3:,5:10] - numparrPMCA[0:3,5:10]
##DN_vals_case1 = numparrFFIoneP1[3:,6:9] - numparrFFIoneP1[0:3,6:9]
##DN_vals_case1 = np.concatenate((DN_vals_case1,(numparrFFIoneP2[3:,0:3] - numparrFFIoneP2[0:3,0:3])), axis = 1) #had to glue data from separate plates
##DN_vals_case2 = numparrFFIoneP2[3:,3:9] - numparrFFIoneP2[0:3,3:9]
##DN_vals_vCJDthal = numparrvCJDthal[3:,3:9] - numparrvCJDthal[0:3,3:9]
DN_vals_vCJDFC = numparrvCJDFC[3:,:10] - numparrvCJDFC[0:3,:10]###have to change vol analysed to 10ul!####
##DN_vals_sCJDMM1FC = numparrFFIoneP1[3:,:6] - numparrFFIoneP1[0:3,:6]
##DN_vals_sCJDMM1Th = numparrsCJDthal[3:,:6] - numparrsCJDthal[0:3,:6]
##DN_vals_sCJDMM2TFC = numparrsCJDMM2TFC[3:,:10] - numparrsCJDMM2TFC[0:3,:10]
##DN_vals_sCJDMM2CFC = numparrsCJDMM2CFC[3:,:10] - numparrsCJDMM2CFC[0:3,:10]

def getResSen(val_array, label, pos2_5 = 2, pos50 = 5): 
    print("D-N values for ", label, val_array)
    totMean, totSD = means_SDs(val_array[:,pos2_5])
    print("Mean tot val and std ",label,":", totMean, totSD)
    resMean, resSD = means_SDs(val_array[:,pos50])
    print("Mean res val and std ",label,":", resMean, resSD)
    senMean, senSD = means_SDs(val_array[:,pos2_5]-val_array[:,pos50])
    print("Mean sen val and std ",label,":", senMean, senSD)
getResSen(DN_vals_vCJDFC, "vCJD FC", pos2_5 = 5, pos50 = 8)

#get values for plotting
##mean_valsPMCA, DNstdPMCA = means_SDs(DN_vals_PMCA)
##print("Mean vals PMCA\n=======\n", mean_valsPMCA)
##mean_valsPMCA2, DNstdPMCA2 = means_SDs(DN_vals_PMCA2)
##print("Mean vals PMCA2\n=======\n", mean_valsPMCA2)
##mean_valsBr, DNstdBr = means_SDs(DN_vals_case1)
##print("Mean vals Br\n=======\n", mean_valsBr)
##mean_valsBr2, DNstdBr2 = means_SDs(DN_vals_case2)
##print("Mean vals Br2\n=======\n", mean_valsBr2)

#put above values in list
##mean_val_list = [
##    ["FFI case 1 Brain", mean_valsBr, DNstdBr, ('0', '1', '2.5','5', '10', '50')],
##    ["FFI case 1 PMCA", mean_valsPMCA, DNstdPMCA],
##    ["FFI case 2 Brain", mean_valsBr2, DNstdBr2, ('0', '1', '2.5','5', '10', '50')],
##    ["FFI case 2 PMCA", mean_valsPMCA2, DNstdPMCA2]]

#from list, plot
def assembleBar(compnts):
    for i in range(len(compnts)):
        plt.subplot((math.ceil(len(compnts)/2)), 2, i+1)
        plotBar(*compnts[i])

#assembleBar(mean_val_list)

#plt.subplots_adjust(wspace = 0.4)
#plt.subplots_adjust(hspace = 0.6)
#plt.show()


