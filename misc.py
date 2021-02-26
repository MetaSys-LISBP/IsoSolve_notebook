import math
import numpy as np
from itertools import chain
import itertools
import string
import re
import numpy as np
import pandas as pd
from operator import itemgetter
from IPython.display import Markdown, display


####################################################
# CREATE ISOTOPOMER-MEASUREMENTS MAPPING DATAFRAME #
####################################################

def repeat(x):
    _size = len(x)
    repeated = []
    for i in range(_size):
        k = i + 1
        for j in range(k, _size):
            if x[i] == x[j] and x[i] not in repeated:
                repeated.append(x[i])
    return repeated

def build_mapping_dataframe(datasets):
   
    # determine length of the carbon skeleton (based on the first dataset)
    ltmp = next(iter(datasets.items()))
    if type(ltmp[1]) is dict:
        n = len(list(ltmp[1].keys())[0])
    elif type(ltmp[1]) is str:
        n = len(ltmp[1])
    elif type(ltmp[1]) is tuple:
        n = len(ltmp[1][0])
    else:
        raise TypeError("Unknown type of measurements for dataset: '{}', '{}' (expected dict, str or tuple, got '{}').".format(ltmp[0], ltmp[1], type(ltmp[1])))        
    
    # check format of each measurement of each dataset
    pattern_nmr, pattern_ms = '^[01x]+$', '^[xE]+$'
    for k,v in datasets.items():
        if type(v) is dict:
            if len(set(v.values())) < 2:
                raise ValueError("Wrong format for dataset: '{}' (at least two (sets of) isotopic species should be measured).".format(k))
            for i in v.keys():
                if not re.match(pattern_nmr, i):
                    raise ValueError("Wrong format for dataset: '{}', '{}'='{}' (isotopomer/cumomer should contain only '0', '1' and 'x').".format(k, i, v[i]))
                if len(i) != n:
                    raise ValueError("Wrong length of carbon skeleton for dataset: '{}', '{}'='{}' (expected '{}', got '{}').".format(k, i, v[i], n, len(i)))
        elif type(v) is str:
            if not re.match(pattern_ms, v):
                raise ValueError("Wrong format for dataset: '{}', '{}' (EMUs should contain only 'x' and 'E').".format(k, v))
            if len(v) != n:
                raise ValueError("Wrong length of carbon skeleton for dataset: '{}', '{}' (expected '{}', got '{}').".format(k, v, n, len(v)))
        elif type(v) is tuple:
            if not (re.match(pattern_ms, v[0]) and re.match(pattern_ms, v[1])):
                raise ValueError("Wrong format for dataset: '{}', '{}' (EMUs should contain only 'x' and 'E').".format(k, v))
            if len(v[0]) != n:
                raise ValueError("Wrong length of carbon skeleton in parent ion for dataset: '{}', '{}' (expected '{}', got '{}').".format(k, v, n, len(v[0])))
            if len(v[1]) != n:
                raise ValueError("Wrong length of carbon skeleton in daughter ion for dataset: '{}', '{}' (expected '{}', got '{}').".format(k, v, n, len(v[1])))
            if (v[0].count("E") == 0 or v[1].count("E") == 0):
                raise ValueError("Parent and daughter ion must contains at least one carbon atom for dataset: '{}', '{}'.".format(k, v))
            if v[0].count("E") <= v[1].count("E"):
                raise ValueError("Daughter ion must contains less carbon atom than parent ion for dataset: '{}', '{}'.".format(k, v))
            lmap = [l for l,v in enumerate([*v[1]]) if v == "E"]
            if itemgetter(*lmap)([*v[0]]) != itemgetter(*lmap)([*v[1]]):
                raise ValueError("All atoms contained in daughter ion must be in parent ion for dataset: '{}', '{}'.".format(k, v))
        else:
            raise TypeError("Unknown type of measurements for dataset: '{}', '{}' (expected dict, str or tuple, got '{}').".format(k, v, type(v)))
            
    # check no measurements are present in different datasets
    l_meas = [j for v in datasets.values() if type(v) is dict for j in set(v.values())]
    rep = repeat(l_meas)
    if len(rep) != 0:
        raise ValueError("The following symbol(s) is (are) present in different datasets: {}".format(rep))
   
    # list all isotopomers
    isotopomers = []
    for r in range(0, n + 1):
        combinations_object = itertools.combinations(range(0, n), r)
        for i in combinations_object:
            isotopomers.append("".join("1" if k in i else "0" for k in range(0, n)))

    # generate two-letters symbols to map MS measurements (isotopologues) as isotopomers
    letters = [''.join(i) for i in itertools.product(string.ascii_lowercase, repeat = 2)]
    new_meas_sym = np.setdiff1d(letters, l_meas).tolist()
    
    # map measurements for each dataset
    res = {}
    for name, meas in datasets.items():
        # initialize an empty set (i.e. all isotopomers are NaN)
        res[name] = dict(zip(isotopomers, [np.nan]*2**n))
        # map measurements provided as isotopomers and/or cumomers (typically NMR datasets)
        if type(meas) is dict:
            for k in meas.keys():
                lmap = [j for j in range(len(k)) if k[j] != "x"]
                for i in range(len(isotopomers)):
                    if itemgetter(*lmap)([*isotopomers[i]]) == itemgetter(*lmap)([*k]):
                        if type(res[name][isotopomers[i]]) is not type(np.nan):
                            raise ValueError("Isotopomer '{}' should be contained in only one measurement but appears twice ('{}', '{}') in dataset '{}'.".format(isotopomers[i], res[name][isotopomers[i]], meas[k], name))
                        res[name][isotopomers[i]] = meas[k]
        # map measurements provided as elementary metabolite units (typically MS datasets)
        elif type(meas) is str:
            # get index of carbon atoms
            lmap = [l for l,v in enumerate([*meas]) if v == "E"]
            # generate symbol
            symbol = new_meas_sym.pop(0)
            get_symlist = lambda s, lmp: [s + "_m" + str(i) for i in range(len(lmp) + 1)]
            lsymbols = get_symlist(symbol, lmap)
            while any(item in lsymbols for item in l_meas):
                symbol = new_meas_sym.pop(0)
                lsymbols = get_symlist(symbol, lmap)
            # map measurements
            for i in range(len(isotopomers)):
                assert type(res[name][isotopomers[i]]) is type(np.nan), "Isotopomer '{}' should be contained in only one measurement but appears twice ('{}', '{}') in dataset '{}'.".format(isotopomers[i], res[name][isotopomers[i]], sym_meas, name)
                mfrag = sum(int(v) for v in itemgetter(*lmap)([*isotopomers[i]]))
                res[name][isotopomers[i]] = symbol + "_m" + str(mfrag)
        # map measurements provided as elementary metabolite units with established parent-daughter ions relationships (typically MS/MS datasets)
        elif type(meas) is tuple:
            # get index of carbon atoms
            lmap_parent = [l for l,v in enumerate([*meas[0]]) if v == "E"]
            lmap_daughter = [l for l,v in enumerate([*meas[1]]) if v == "E"]
            # generate symbol
            symbol = new_meas_sym.pop(0)
            get_symlist = lambda s, lmp, lmd: [s + "_m" + str(i) + "_" + str(j) for i in range(len(lmp) + 1) for j in range(min(i, len(lmd)) + 1)]
            lsymbols = get_symlist(symbol, lmap_parent, lmap_daughter)
            while any(item in lsymbols for item in l_meas):
                symbol = new_meas_sym.pop(0)
                lsymbols = get_symlist(symbol, lmap_parent, lmap_daughter)
            # map measurements
            for i in range(len(isotopomers)):
                assert type(res[name][isotopomers[i]]) is type(np.nan), "Isotopomer '{}' should be contained in only one measurement but appears twice ('{}', '{}') in dataset '{}'.".format(isotopomers[i], res[name][isotopomers[i]], sym_meas, name)
                mfrag_parent = sum(int(v) for v in itemgetter(*lmap_parent)([*isotopomers[i]]))
                mfrag_daughter = sum(int(v) for v in itemgetter(*lmap_daughter)([*isotopomers[i]]))
                res[name][isotopomers[i]] = symbol + "_m" + str(mfrag_parent) + "_" + str(mfrag_daughter)
                

    return pd.DataFrame.from_dict(res)


#####################
# CALCULATE METRICS #
#####################

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

def get_theoretical_iso(res, ab_13C):
    return [(1-ab_13C)**i.count("0")*ab_13C**i.count("1") for i in res['xbc'].index]

def get_theoretical_iso_sd(res, ab_13C, err_13C):
    return [(1-ab_13C)**i.count("0")*ab_13C**i.count("1")*np.sqrt((err_13C/(1-ab_13C)*i.count("1"))**2 + (err_13C/ab_13C*i.count("0"))**2) for i in res['xbc'].index]

def get_theoretical_emu(nC, ab_13C):
    return [((1-ab_13C)**(nC-i))*(ab_13C**i)*nCr(nC,i) for i in range(nC+1)]

def get_theoretical_emu_full(res, ab_13C):
    emu_th = {}
    for k,v in res['ls']['emu'][0].items():
        nC = k.count("E")
        emu_th[k] = [((1-ab_13C)**(nC-i))*(ab_13C**i)*nCr(nC,i) for i in range(nC+1)]
    return emu_th

def get_theoretical_emu_full_sd(res, ab_13C, err_13C):
    emu_th_sd = {}
    for k,v in res['ls']['emu'][0].items():
        nC = k.count("E")
        emu_th_sd[k] = [((1-ab_13C)**(nC-i))*(ab_13C**i)*nCr(nC,i) * np.sqrt((err_13C/(1-ab_13C)*nCr(nC,i))**2 + ((err_13C/ab_13C*nCr(nC,i))**2)) for i in range(nC+1)]
    return emu_th_sd

def get_theoretical_cumo(res, ab_13C):
    return [ab_13C**i.count("1") for i in res['cusol'].keys()]

def get_theoretical_cumo_sd(res, ab_13C, err_13C):
    return [i.count("1")*(err_13C/ab_13C)*(ab_13C**i.count("1")) for i in res['cusol'].keys()]

def calculate_mean_error_emu(res, ab_13C):
    l_err = []
    theoretical_EMUs = {}
    for k,v in res['ls']['emu'][0].items():
        nC = k.count("E")
        theoretical_EMUs[k] = get_theoretical_emu(nC, ab_13C)
        l_err += list(np.abs([theoretical_EMUs[k][i] - v[i] for i in range(len(theoretical_EMUs[k]))]))
    return np.mean(l_err)

def calculate_mean_sd_emu(res):
    l_sd = list(itertools.chain.from_iterable([i for i in res['ls']['emu'][1].values()]))
    return np.mean(l_sd)

def calculate_mean_error_iso(res, ab_13C):
    theoretical = get_theoretical_iso(res, ab_13C)
    experimental = res["ls"]["iso"][0].flatten().tolist()[0]
    return np.mean(np.abs([theoretical[i] - experimental[i] for i in range(len(theoretical))]))

def calculate_mean_sd_iso(res):
    return np.mean(res["ls"]["iso"][1].flatten().tolist())

def calculate_mean_error_cumo(res, ab_13C):
    theoretical = get_theoretical_cumo(res, ab_13C)
    experimental = res["ls"]["cumo"][0]
    return np.mean(np.abs([theoretical[i] - experimental[i] for i in range(len(res['cusol']))]))

def calculate_mean_sd_cumo(res):
    return np.mean(res["ls"]["cumo"][1])

def calculate_metrics(res, ab_13C, err_13C):
    return {'mean_error_iso':calculate_mean_error_iso(res, ab_13C),
            'mean_sd_iso':calculate_mean_sd_iso(res),
            'mean_error_emu':calculate_mean_error_emu(res, ab_13C),
            'mean_sd_emu':calculate_mean_sd_emu(res),
            'mean_error_cumo':calculate_mean_error_cumo(res, ab_13C),
            'mean_sd_cumo':calculate_mean_sd_cumo(res)}

def get_coverage_relative(res, lk):
    p_ip = {"amino_acid":lk,
            "isotopomers":[100*res[i]["identifiable"]["isotopomers"]/res[i]["total"]["isotopomers"] for i in lk],
            "cumomers":[100*res[i]["identifiable"]['cumomers']/res[i]["total"]['cumomers'] for i in lk],
            "emus":[100*res[i]["identifiable"]['emus']/res[i]["total"]['emus'] for i in lk]}
    return pd.DataFrame(p_ip, columns=['amino_acid', 'isotopomers', 'cumomers', 'emus'])

def get_coverage(res, lk):
    p_ip = {"amino_acid":lk,
            "isotopomers":[res[i]["identifiable"]["isotopomers"] for i in lk],
            "cumomers":[res[i]["identifiable"]['cumomers'] for i in lk],
            "emus":[res[i]["identifiable"]['emus'] for i in lk]}
    return pd.DataFrame(p_ip, columns=['amino_acid', 'isotopomers', 'cumomers', 'emus'])


###################
# PRINT FUNCTIONS #
###################

def print_summary(res):
    # equations for isotopic species that can be quantified
    display(Markdown("**Equations**"))
    df = pd.DataFrame({'isotopomer':[res['sbc'][i] for i in res['xbc'].index], 'equation':res['sols']}).astype('str')
    df = df.query("isotopomer != equation")
    df = df.replace(r"i0", r"0", regex=True).replace(r"i1", r"1", regex=True)
    df.reset_index(inplace=True, drop=True)
    display(df)
    # results summary
    display(Markdown("**Summary**"))
    #   non quantifiable species
    undefined_species = [res['xbc'].index[i] for i in range(len(res['xbc'].index)) if i not in res['metrics']["idef"]]
    display(Markdown("Number of undefined isotopic species: {}".format(len(undefined_species))))
    if len(undefined_species):
        display(Markdown("&nbsp; &nbsp; *Undefined isotopic species: {}*".format(", ".join(i for i in undefined_species))))
    #   redundant measurements
    nb_rdd_meas = res['rdn_meas']
    display(Markdown("Number of redundant measurements: {}".format(len(nb_rdd_meas))))
    if len(nb_rdd_meas):
        display(Markdown("&nbsp; &nbsp; *Redundant measurements: {}*".format(", ".join(i for i in nb_rdd_meas))))
