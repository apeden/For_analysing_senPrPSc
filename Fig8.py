import matplotlib.pyplot as plt
import random
import numpy as np
import pandas as pd
import math

toPrt = False

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

def std_curve(numparr, curve_start_col, curve_start, curve_end):
    """CDI standard curve coeff and intercept on y axis"""
    X = np.array([0.25, 0.5, 0.75, 1][curve_start:curve_end+1])
    curve_points = [numparr[4:7, curve_start_col+1],
                  numparr[4:7, curve_start_col],
                  numparr[1:4, curve_start_col+1],
                  numparr[1:4, curve_start_col]][curve_start:curve_end+1]
    y = np.array(curve_points,dtype = float)
    if toPrt: print("y vals for standard curve\n=======\n",y)
    result = np.polyfit(X,y,1)
    coeff = result[0].mean()
    intercept = result[1].mean()
    if toPrt: print("CDI value per ug PrP = ",str(coeff))
    if toPrt: print("Intercept on y axis =", str(intercept))
    return coeff, intercept

def relate_to_PrP(numparr, coeff, intercept, ul_analysed = 25):
    """convert CDI values to PrP per g brain, using std curve"""
    numparr = (numparr[1:7,0:10]-intercept)/coeff
    if toPrt: print("Values as ug PrP per 2.5mg brain sample\n========\n",numparr)
    numparr = numparr/(ul_analysed/10000)
    if toPrt: print("Values as ug PrP per gram brain \n========\n",numparr)
    return numparr

def get_vals(file_name, sheet_name, curve_col_start, curve_start = 0, curve_end= 3):
    """gets data and converts CDI vals to concentrations of PrP
    Returns as numpy array 
    """
    numparr = to_numpy(file_name, sheet_name)
    numparr = substrate_blanks(numparr, 0,7)
    coeff, intercept = std_curve(numparr, curve_col_start, curve_start, curve_end)
    numparr = relate_to_PrP(numparr, coeff, intercept)
    return numparr
         
def means_SDs(vals):
    """takes array of values and returns mean and sd"""
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
    plt.ylabel(chr(956)+"g PrP /gram Brain")
    plt.xlabel("PK "+chr(956)+"g/ml")
    plt.title(title)
    plt.xticks(ind, PK_concs)
    heighest_val = max(mean_vals) + max(DNstd)
    gradation = 10
    y_scale = 100 + gradation
    if heighest_val > 100:
        y_scale = (round(heighest_val,
                         -int(math.log(heighest_val,10))))
        math.log(y_scale, 10)
        gradation = 10**(int(math.log(y_scale, 10)))     
        y_scale += gradation
    plt.yticks(np.arange(0, y_scale, gradation))

def getD_N_vals(arr, first_col,last_col ):
    return arr[3:,first_col:last_col+1]-arr[0:3,first_col:last_col+1] 

def printResSen(DN_array, label, pos2_5 = 2, pos50 = 5): 
    """for printing means and std for:
    total PrPSc (that remaining after 2.5ug/ml PK)
    res PrPSc   (that remaining after 50ug/ml  PK)
    sen PrPSc   (total - res PrPSc)

    DN_array is an array of calibrated D-N values.
    pos2_5 = 2, pos50 = 5 are the positions of the\
    columns containing samples treated with 2.5ug/ml PK\
    and 50ug/ml PK, respectively
    """
    print("D-N values for ", label, DN_array)
    totMean, totSD = means_SDs(DN_array[:,pos2_5])
    print("Mean tot val and std ",label,":", totMean, totSD)
    resMean, resSD = means_SDs(DN_array[:,pos50])
    print("Mean res val and std ",label,":", resMean, resSD)
    senMean, senSD = means_SDs(DN_array[:,pos2_5]-DN_array[:,pos50])
    print("Mean sen val and std ",label,":", senMean, senSD)

def _PrPC_PrPSc(val_array, col_dict):
    """calculates mean/std N value, D value, D-N value\
    for each column in a dictionary of columns (dict keys)\
    with associated str labels (dict values)
    """
    val_array = np.array(val_array, dtype=np.float64)
    data_dict  = {"metric":["N","Nstd","D","Dstd","D-N","D-Nstd"]}
    for column in col_dict:
        data_list = []
        data_list.append(np.mean(val_array[0:3,column]))
        data_list.append(np.std(val_array[0:3,column]))
        data_list.append(np.mean(val_array[3:,column]))
        data_list.append(np.std(val_array[3:,column]))
        data_list.append(np.mean(val_array[3:,column] - val_array[0:3,column]))
        data_list.append(np.std(val_array[3:,column] - val_array[0:3,column]))
        data_dict[str(column)+": "+col_dict[column]] = data_list
        del data_list
    return pd.DataFrame(data_dict)

##############  get data into numpy arrays                 ######################
##numparrPMCA =         get_vals("CDI13 030B.xls",          "Plate", curve_col_start = 10)
numparrFFIoneP1 =     get_vals("CDI13 004.xls",           "Plate", curve_col_start = 9)
numparrFFIoneP2 =     get_vals("CDI13 004b.xls",          "Plate", curve_col_start = 9)
##numparrvCJDFC=        get_vals("CDI 12 007 PLATE B.xls",  "Plate", curve_col_start = 10)
numparrFFIoneThal=    get_vals("CDI 12 016 PLATE A.xls",  "Plate", curve_col_start = 9)
##numparrvCJDthal=      get_vals("CDI 12 016 PLATE B.xls",  "Plate", curve_col_start = 9, curve_start = 1, curve_end =3)
##numparrsCJDthal=      get_vals("CDI13 004b.xls",          "Plate", curve_col_start = 6)
numparrsCJDMM2TFC=    get_vals("CDI 12 008 PLATE B.xls",  "Plate", curve_col_start = 10)
##numparrsCJDMM2CFC=    get_vals("CDI 12 009 PLATE A.xls",  "Plate", curve_col_start = 10)
##numparrNonCJDone =    get_vals("CDI13 007.xls",           "Plate", curve_col_start = 10)
##numparrNonCJDtwo =    get_vals("CDI13 007b.xls",          "Plate", curve_col_start = 10)
numparrsFFIone2011 =  get_vals("CDI 11 013 A FFI vs vCJD vs ALZ.xls", "Plate", curve_col_start = 10)
numparrsFFItwo2011 =  get_vals("CDI 11 014 FFI and vCJD A.xls", "Plate", curve_col_start = 10, curve_end = 2)
numparrsFFIone2011b =  get_vals("CDI 11 015 FFI and sCJD A.xls", "Plate", curve_col_start = 10)
numparrsFFItwo2011b =  get_vals("CDI 11 016 FFI and Alz stability A.xls", "Plate", curve_col_start = 10, curve_end = 2)
numparrsFFIone2011c =  get_vals("CDI 11 017 vCJD and FFI meltcurve A.xls", "Plate", curve_col_start = 10)
numparrsFFItwo2011c =  get_vals("CDI 11 018 FFI vs sCJD plate A.xls", "Plate", curve_col_start = 10)
numparrsFFIone2011d =  get_vals("CDI 11 019 FFI and LBD A.xls", "Plate", curve_col_start = 10)
numparrsFFIone2011e =  get_vals("CDI 11 020 FFI and vCJD A.xls", "Plate", curve_col_start = 10)
numparrsFFItwo2011d =  get_vals("CDI 11 020 FFI and vCJD B.xls", "Plate", curve_col_start = 10)
numparrsFFIone2011f =  get_vals("CDI 11 023 melt curves.xls", "Plate", curve_col_start = 10)
numparrFFIoneThal2 = get_vals("CDI 12 001 PLATE B.xls", "Plate", curve_col_start = 10)

##############  numpy arrays of calibrate D - N CDI vals   ######################
##DN_vals_PMCA        = getD_N_vals(numparrPMCA,          0, 4)
##DN_vals_PMCA        = getD_N_vals(numparrPMCA,          5, 9)
DN_vals_case1       = getD_N_vals(numparrFFIoneP1,      8, 8)
DN_vals_case1       = np.concatenate((DN_vals_case1,getD_N_vals(numparrFFIoneP2, 0, 2)), axis = 1) #had to glue data from separate plates
##DN_vals_case2       = getD_N_vals(numparrFFIoneP2,      5, 8)
##DN_vals_FFIoneThal  = getD_N_vals(numparrFFIoneThal,    0, 5)
##DN_vals_FFItwoThal  = getD_N_vals(numparrFFIoneThal,    6, 8)
##DN_vals_FFItwoThal  = np.concatenate((DN_vals_FFItwoThal, getD_N_vals(numparrvCJDthal, 0, 2)), axis = 1) #had to glue data from separate plates
##DN_vals_vCJDthal    = getD_N_vals(numparrvCJDthal,      3, 8)
##DN_vals_vCJDFC      = getD_N_vals(numparrvCJDFC,        0, 9)###have to change vol analysed to 10ul!####
##DN_vals_sCJDMM1FC   = getD_N_vals(numparrFFIoneP1,      2, 5)
##DN_vals_sCJDMM1Th   = getD_N_vals(numparrsCJDthal,      0, 5)
##DN_vals_sCJDMM2TFC  = getD_N_vals(numparrsCJDMM2TFC,    0, 9)
##DN_vals_sCJDMM2CFC  = getD_N_vals(numparrsCJDMM2CFC,    0, 9)
##DN_vals_58          = getD_N_vals(numparrNonCJDone,     2, 4)
##DN_vals_39          = getD_N_vals(numparrNonCJDone,     7, 9)
##DN_vals_54          = getD_N_vals(numparrNonCJDtwo,     2, 4)
##DN_vals_65          = getD_N_vals(numparrNonCJDtwo,     7, 9)
##DN_nonCJD   = DN_vals_58, DN_vals_39, DN_vals_54, DN_vals_65


########## PrPres and senPrPSc #########################
##printResSen(DN_vals_case1, "FFI case 1")
##printResSen(DN_vals_FFItwoThal, "FFI case 2 Thal")

################ N, D and D-N  #########################    


results = [#a list of dataframes
(_PrPC_PrPSc(numparrFFIoneP1, {0:"FFI cases 1 FC",
                                    6:"FFI cases 2 FC"})),
(_PrPC_PrPSc(numparrsFFIone2011, {3:"FFI cases 1 FC 2011"})),
(_PrPC_PrPSc(numparrsFFItwo2011, {3:"FFI cases 2 FC 2011"})),
(_PrPC_PrPSc(numparrsFFIone2011b, {3:"FFI cases 1 FC 2011b"})),
(_PrPC_PrPSc(numparrsFFItwo2011b, {9:"FFI cases 2 FC 2011b"})),
(_PrPC_PrPSc(numparrsFFIone2011c, {9:"FFI cases 1 FC 2011c"})),
(_PrPC_PrPSc(numparrsFFItwo2011c, {9:"FFI cases 2 FC 2011c"})),
(_PrPC_PrPSc(numparrsFFIone2011d, {9:"FFI cases 1 FC 2011d"})),
(_PrPC_PrPSc(numparrsFFIone2011e, {9:"FFI cases 1 FC 2011e"})),
(_PrPC_PrPSc(numparrsFFIone2011f, {9:"FFI cases 1 FC 2011f"})),
(_PrPC_PrPSc(numparrsFFItwo2011d, {9:"FFI cases 2 FC 2011d"})),
(_PrPC_PrPSc(numparrFFIoneThal, {0:"FFI cases 1 Thal",
                                    6:"FFI cases 2 Thal"})),
(_PrPC_PrPSc(numparrsCJDMM2TFC, {3:"sFI FC"})),
(_PrPC_PrPSc(numparrFFIoneThal2, {3:"FFI cases 2 Thal 2"}))
]

print("An example of one data frame\n",results[0])
master_df = results[0].T
for df in results[1:]:
    master_df = master_df.append(df.loc[:, df.columns != 'metric'].T)
print("All the data transposed\n",master_df)
master_df = master_df.reset_index()# reset index as 0, 1... etc. What was the index is now a column called 'index'
print("All the data, reindexed\n",master_df)

def getMeanN_D_DN(df, subset):
    """for a dataframe and a string describing a subset of
    samples, return of tuple of (mean N, std N, mean D, std D,
    mean D-N, std D-N
    """
    subset_df = df[df['index'].str.contains(subset)]
    print("All data with ",subset)
    data_dict  = {"metric":["N","Nstd","D","Dstd","D-N","D-Nstd","n"]}
    mean_std = []
    for i in range (0,5,2): 
        series = pd.to_numeric(subset_df[i], errors='ignore')
        mean_std += [series.mean(), series.std()]
    mean_std.append(int(len(subset_df)))
    data_dict[subset] = mean_std
    return pd.DataFrame(data_dict)

print(getMeanN_D_DN(master_df, "FFI cases 2"))



#print(_PrPC_PrPSc(numparrFFIoneP2, 0, 8))
###########  get values for plotting ###################
##mean_valsPMCA, DNstdPMCA = means_SDs(DN_vals_PMCA)
##print("Mean vals PMCA\n=======\n", mean_valsPMCA)
##mean_valsPMCA2, DNstdPMCA2 = means_SDs(DN_vals_PMCA2)
##print("Mean vals PMCA2\n=======\n", mean_valsPMCA2)
##mean_valsBr, DNstdBr = means_SDs(DN_vals_case1)
##print("Mean vals Br\n=======\n", mean_valsBr)
##mean_valsBr2, DNstdBr2 = means_SDs(DN_vals_case2)
##print("Mean vals Br2\n=======\n", mean_valsBr2)
##mean_valsSCJDFC, DNstdSCJDFC = means_SDs(DN_vals_sCJDMM1FC)
##print("Mean vals sCJD FC\n=======\n", mean_valsSCJDFC)

############# calculate mean of means and std ###################
##mean_list = []
##for arr in DN_nonCJD:
##    mean_vals, _ = means_SDs(arr)
##    print("Mean vals ",str(arr),"\n=======\n", mean_vals)
##    mean_list.append(mean_vals)
##print("All the means\n=======\n", mean_list)
##mean_of_mean_nonCJD = np.vstack(mean_list)
##mean_nonCJD, DNstd_nonCJD = means_SDs(mean_of_mean_nonCJD)
##print("Mean vals nonCJD FC\n=======\n", mean_nonCJD)  
####put above values in list
##mean_val_list = [
##    ["FFI case 1 FC", mean_valsBr, DNstdBr, ('2.5','5', '10', '50')],
##    ["FFI case 2 FC", mean_valsBr2, DNstdBr2, ('2.5','5', '10', '50')],
##    ["sCJD MM1 FC", mean_valsSCJDFC, DNstdSCJDFC, ('2.5','5', '10', '50')],
##    ["Non CJD", mean_nonCJD, DNstd_nonCJD, ('2.5', '10', '50')]]

#from list, plot
def assembleBar(compnts):
    for i in range(len(compnts)):
        plt.subplot((math.ceil(len(compnts)/2)), 2, i+1)
        plotBar(*compnts[i])

##assembleBar(mean_val_list)
##plt.subplots_adjust(wspace = 0.4)
##plt.subplots_adjust(hspace = 0.6)
##plt.show()


