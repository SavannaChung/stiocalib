# rspc plot function
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import rspc_read_file as rf

import rspc_analyse as ra


def r2(yexp, ymea):
    """ find the r-squared value of yexp and ymea
    ---------------------------------- input ----------------------------------
    yexp, ymea = np.array
    ---------------------------------- output ----------------------------------
    r2 = np.float
    """
    yexp = np.array(yexp)
    ymea = np.array(ymea)
    mean_y = ymea.mean()
    ssreg = np.sum((yexp - ymea)**2)
    sstot = np.sum((ymea - mean_y)**2)
    r2 = 1 - ssreg/ sstot

    return r2

def plot_lsm_pho_en(dict, insert_name):
    """ plots the least square value as a function of photon energy.
        ----------------------------------------------------------------------
        input:
        dict= a dictionary {key = photon energy: value = lsm}
             could be output from calc_mat_ctno_from_lac() in rspc_calc.py

        insert_name = a string ( gammex, old_cirs, new_cirs)
        ----------------------------------------------------------------------
        output: figure
        ----------------------------------------------------------------------
    """
    x = list(dict.keys())
    y = list(dict.values())

    eff_en = min(dict, key = dict.get)
    # max_y = max(dict, key = dict.get(np.arange(0.06, 0.09, 0.001)))
    # print(dict[np.arange(0.06, 0.09, 0.001)])

    plt.rcParams["figure.figsize"] = (5,4)
    plt.plot(x , y, "k+")
    plt.xlabel("photon energy (MeV)")
    plt.ylabel("the sum of the squared difference (a.u.)")
    plt.xticks(rotation = 45)
    plt.title("Insert = %s || Effective energy = %5.4f MeV" % (insert_name, eff_en), fontsize = 10)
    # plt.annotate("effective energy = 0.074 MeV", xy=(eff_en, dt_ctno_lsm[eff_en]), xytext = (eff_en, dt_ctno_lsm[eff_en] +30000),  arrowprops=dict(arrowstyle="->"))
    plt.ylim(0,(dict.get(eff_en))*150)
    plt.xlim(0.06, 0.090)
    plt.tight_layout()
    plt.show()

    return

def plot_sorted_calc_ctno(calc):
    """ plot the calc ctno as a function of inserts. To decide the lung, soft tissue and bone ct number region.
    --------------------------------------------------------------------------
    input: dictionaries
    calc_ctno = {material: calc_ctno, ...}

    --------------------------------------------------------------------------
    output: figure

    """
    calc_x, calc_y = zip(*sorted(calc.items(), key = lambda x: x[1]))

    # plt.rcParams["figure.figsize"] = (5,4)
    plt.plot(calc_x , calc_y, "k+", fillstyle = "none", label = "calc xcom" )
    plt.xlabel("Inserts (a.u.)")
    plt.ylabel("HU (a.u.)")
    plt.xticks(rotation = 45)
    plt.legend()

    plt.show()
    return

def plot_xcalc_meas_ctno(calc, meas, insert_name, eff_en):
    """ plot the calc and meas ctno of different materials
    --------------------------------------------------------------------------
    input: dictionaries
    calc_hu = {material: calc_ctno, ...}
    meas_hu = {material: means_ctno, ...}
    insert_name = string(gammex / old_cirs / new_cirs)
    eff_en = float
    --------------------------------------------------------------------------
    output: figure

    """
    calc_x, calc_y = zip(*sorted(calc.items(), key = lambda x: x[1]))
    meas_x, meas_y = zip(*sorted(meas.items(),  key = lambda x: x[1]))

    print(f"calc_x :{calc_x}")
    print(f"calc_y :{calc_y}")
    print(f"meas_x :{meas_x}")
    print(f"meas_y :{meas_y}")

    # calculate the percentage differences
    perc_diff = {keys:0 for keys in calc_x}
    for mat in perc_diff.keys():
        diff = 100*(meas[mat]-calc[mat])/calc[mat]
        perc_diff[mat] = diff

    # -------------------------- plot -------------------------- #

    fig, axes = plt.subplots(nrows= 2, ncols =1, figsize = (4,6))

    # plt.rcParams["figure.figsize"] = (5,4)
    axes[0].plot(calc_x , calc_y, "k+", fillstyle = "none", label = "calc xcom" )
    axes[0].plot(meas_x, meas_y, "ro", fillstyle = "none", label = "meas")
    axes[0].set_xlabel("Inserts (a.u.)")
    axes[0].set_ylabel("CT no (a.u.)")
    axes[0].tick_params(axis= "x", labelrotation = 90)
    axes[0].set_title("insert= %s || effective energy = %5.4f MeV" % \
                        (insert_name, eff_en), fontsize = 8)
    axes[0].legend()

    perd_x, perd_y = zip(*perc_diff.items())

    axes[1].plot(perd_x, perd_y, "k^", fillstyle = "none", label = "$\%$ diff." )
    axes[1].axhline(y = 0, xmin = 0, xmax =2, color = "gray")
    axes[1].set_xlabel("Inserts (a.u.)")
    axes[1].set_ylabel("percent. diff. ($\%$)")
    axes[1].tick_params(axis= "x", labelrotation = 90)
    axes[1].legend()
    plt.tight_layout()
    plt.show()

    return

def plot_kcalc_meas_ctno(calc, meas, insert_name, eff_en, Kph, Kcoh, KN):
    """ plot the meas_hu as a function of calc_hu
    --------------------------------------------------------------------------
    input: dictionaries
    calc_ctno = {material: calc_ctno, ...}
    meas_ctno = {material: means_ctno, ...}
    insert_name = string(gammex / old_cirs / new_cirs)
    eff_en = float
     Kph, Kcoh, KN = str
    --------------------------------------------------------------------------
    output: figure

    """
    calc_x, calc_y = zip(*sorted(calc.items(), key = lambda x: x[1]))
    meas_x, meas_y = zip(*sorted(meas.items(),  key = lambda x: x[1]))

    # calculate the percentage differences
    perc_diff = {keys:0 for keys in calc_x}
    for mat in perc_diff.keys():
        diff = 100*(meas[mat]-calc[mat])/(calc[mat] + 1000)
        perc_diff[mat] = diff
        # print(f"mat: {mat} || meas[mat]: {meas[mat]} || calc[mat]: {calc[mat]} || perc_diff[mat]: {perc_diff[mat]}")

    # -------------------------- plot -------------------------- #

    fig, axes = plt.subplots(nrows= 2, ncols =1, figsize = (4,6))

    # plt.rcParams["figure.figsize"] = (5,4)
    axes[0].plot(calc_x , calc_y, "k+", fillstyle = "none", label = "calc k")
    axes[0].plot(meas_x, meas_y, "ro", fillstyle = "none", label = "meas")
    axes[0].set_xlabel("Inserts (a.u.)")
    axes[0].set_ylabel("Ctno (a.u.)")
    axes[0].tick_params(axis= "x", labelrotation = 90)
    axes[0].set_title("insert= %s || effective energy = %5.4f MeV \n Kph = %s, Kcoh = %s, KN = %s" % \
                        (insert_name, eff_en, Kph, Kcoh, KN), fontsize = 8)
    axes[0].legend()

    perd_x, perd_y = zip(*perc_diff.items())

    axes[1].plot(perd_x, perd_y, "k^", fillstyle = "none", label = "$\%$ diff." )
    axes[1].axhline(y = 0, xmin = 0, xmax =2, color = "gray")
    axes[1].set_title("based on ctno not HU" , fontsize = 8)
    axes[1].set_xlabel("Inserts (a.u.)")
    axes[1].set_ylabel("percent. diff. ($\%$)")
    axes[1].tick_params(axis= "x", labelrotation = 90)
    axes[1].legend()
    plt.tight_layout()
    plt.show()

    return

def plot_calc_meas_prsp(calc, meas, insert_name):
    """ plot the calc and meas prsp of different materials
    --------------------------------------------------------------------------
    input: dictionaries
    calc_ctno = {material: calc_prsp, ...}
    meas_ctno = {material: means_prsp, ...}
    insert_name = string(gammex / old_cirs / new_cirs)
    eff_en = float
    --------------------------------------------------------------------------
    output: figure

    """
    calc_x, calc_y = zip(*sorted(calc.items(), key = lambda x: x[1]))
    meas_x, meas_y = zip(*sorted(meas.items(),  key = lambda x: x[1]))

    # calculate the percentage differences
    perc_diff = {keys:0 for keys in calc_x}
    for mat in perc_diff.keys():
        diff = 100*(meas[mat]-calc[mat])/calc[mat]
        perc_diff[mat] = diff
        # print(f"mat: {mat} || calc[mat]: {calc[mat]} || mea[mat]: {meas[mat]} || diff:{diff}")

    # -------------------------- plot -------------------------- #

    fig, axes = plt.subplots(nrows= 2, ncols =1, figsize = (4,6))

    # plt.rcParams["figure.figsize"] = (5,4)
    axes[0].plot(calc_x , calc_y, "k+", fillstyle = "none", label = "calc" )
    axes[0].plot(meas_x, meas_y, "ro", fillstyle = "none", label = "meas")
    axes[0].set_xlabel("Inserts (a.u.)")
    axes[0].set_ylabel("prsp (a.u.)")
    axes[0].tick_params(axis= "x", labelrotation = 90)
    axes[0].set_title("insert= %s " % \
                        (insert_name), fontsize = 8)
    axes[0].legend()

    perd_x, perd_y = zip(*perc_diff.items())

    axes[1].plot(perd_x, perd_y, "k^", fillstyle = "none", label = "$\%$ diff." )
    axes[1].axhline(y = 0, xmin = 0, xmax =2, color = "gray")
    axes[1].set_xlabel("Inserts (a.u.)")
    axes[1].set_ylabel("percent. diff. ($\%$)")
    axes[1].tick_params(axis= "x", labelrotation = 90)
    axes[1].legend()
    plt.tight_layout()
    plt.show()

    return

# linear fit for three regions. Lung (0 < ctno < 850) || soft tissue: (1023 < H < 1060) || bone: (H > 1060)
def plot_prsp_ctno(ct, prsp, insert_name):
    """ plot the calc_prsp as a function of calc_ctno
    --------------------------------------------------------------------------
    input: dictionaries
    calc_prsp = {material: calc_prsp, ...}
    calc_prsp = {material: calc_prsp, ...}
    insert_name = string(gammex / old_cirs / new_cirs)
    --------------------------------------------------------------------------
    output: figure

    """
    # store the lists of materials in different ctno range in a dictionary
    dt = {("ls_mat_" + str(i)) :[] for i in np.arange(0,3,1)}

    px = [np.arange(-1000, -300,10), np.arange(-200, 150,1), np.arange(200, max(list(ct.values()))+ 1000,10)] # HU
    linestyles = ["-.", ":", "--"]
    linecolor = [[0, 0, 0] , [0.7, 0.7, 0.7] , [0.3, 0.3, 0.3]]

    for key in ct.keys():
        if ct[key] <= -300: # lung
            dt["ls_mat_0"].append(key)
        elif (ct[key] >= -200) & (ct[key] < 150): # soft Adipose_tissue ( muscle = 1106 ct no) 1150
            dt["ls_mat_1"].append(key)
        elif ct[key] >= 200: # bone
            dt["ls_mat_2"].append(key)

    fig, ax = plt.subplots(dpi = 100)
    pt = plt.scatter(list(ct.values()), list(prsp.values()), s=5, c="k",marker ="^")
    # pt = plt.plot(list(ct.values()), list(prsp.values()), "k+")

    slope = []
    const = []
    for i in np.arange(0,3,1):
        mat = dt["ls_mat_" + str(i)]
        # print(f"mat: {mat}")
        x = [] # ct no
        y = [] # prsp
        for m in mat:
            x.append(ct[m])
            y.append(prsp[m])

        x1 = [] # ct no
        y1 = [] # prsp

        model = list(np.polyfit(x, y, 1))

        slope.append(model[0])
        const.append(model[1])

        predict = np.poly1d(model)

        x1 = np.array(px[i])
        y1 = np.array(predict(x1))

        rsv = r2(predict(x), y)

        ax.plot(x1, y1,  linestyles[i], linewidth= 1, \
            label = "%5.9f x + % 5.9f (r2: %5.3f)" % (slope[i], const[i], rsv))

    annot= ax.annotate("", xy = (0,0), xytext = (-20,20), textcoords = "offset points", \
                       bbox = dict(boxstyle = "sawtooth", fc = "w"), \
                       arrowprops = dict(arrowstyle = "->"))
    annot.set_visible(False)

    def update_annot(ind):
        # Get the x,y coordinates of the selected data point
        pos = pt.get_offsets()[ind["ind"][0]]
        # print(f"pt.get_offsets(): {pt.get_offsets()}") # get the x, y coordinates of all datapoints
        # print(ind["ind"][0], type(ind["ind"][0]) ) # ind["ind"] is an array which contains the index number of the selected datapoint
        annot.xy = pos

        point_x = list(ct.values())[ind["ind"][0]]
        point_y = list(prsp.values())[ind["ind"][0]]
        key_index =  list(ct.keys())[list(ct.values()).index(point_x)]
        text = "{} \n HU:{:.4f} \n prsp:{:.4f})".format(key_index, point_x, point_y)
        # print(f"text: {text}")

        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.6)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = pt.contains(event) # cont = True/ False (is it a data point?) || ind = datapoint index number

            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("button_press_event", hover)
    ax.set_xlabel("HU (a.u.)")
    ax.set_ylabel("prsp (a.u.)")
    ax.set_title("insert= %s || tissue composition: %s" % (insert_name, ra.tc_source), fontsize =8)
    ax.legend()

    plt.tight_layout()
    plt.show()

    return slope, const

# plot calc. ctno and prsp together with mea. HU of different tissue (this function modified "plot_prsp_ctno()")
def plot_rsp_w_mea_hu(ct, prsp, insert_name ):

    """ plot the calc_prsp as a function of calc_ctno
    --------------------------------------------------------------------------
    input: dictionaries
    calc_prsp = {material: calc_prsp, ...}
    calc_prsp = {material: calc_prsp, ...}
    insert_name = string(gammex / old_cirs / new_cirs)
    mean_mea_hu = {tissue: mean HU of a tissue}
    std_mea_hu = {tissue: std of a tissue}
    structure = body structure for the HU validation
    --------------------------------------------------------------------------
    output: figure

    """
    # read the siemens hu dataa
    mean_mea_hu, std_mea_hu, structure = rf.read_siemens_hu_data()
    # store the lists of materials in different ctno range in a dictionary
    dt = {("ls_mat_" + str(i)) :[] for i in np.arange(0,3,1)}

    # px = [np.arange(0,850,10), np.arange(1023,1060,1), np.arange(1060, max(list(ct.values())),10)] # ct no (Schneider 1996)
    # px = [np.arange(-1000, -150,10), np.arange(23,60,1), np.arange(60, max(list(ct.values()))+ 1000,10)] # HU (Schneider 1996)
    px = [np.arange(-1000, -300,10), np.arange(-200, 150,1), np.arange(200, max(list(ct.values()))+ 1000,10)] # HU
    linestyles = ["-.", ":", "--"]
    linecolor = [[0, 0, 0] , [0.7, 0.7, 0.7] , [0.3, 0.3, 0.3]]

    for key in ct.keys():
        if ct[key] <= -300: # lung
            dt["ls_mat_0"].append(key)
        elif (ct[key] >= -200) & (ct[key] < 150): # soft Adipose_tissue ( muscle = 1106 ct no) 1150
            dt["ls_mat_1"].append(key)
        elif ct[key] >= 200: # bone
            dt["ls_mat_2"].append(key)

    fig, ax = plt.subplots()
    pt = plt.scatter(list(ct.values()), list(prsp.values()), s=5, c="k",marker ="^")
    # pt = plt.plot(list(ct.values()), list(prsp.values()), "k+")

    slope = []
    const = []
    for i in np.arange(0,3,1):
        mat = dt["ls_mat_" + str(i)]
        # print(f"mat: {mat}")
        x = [] # ct no
        y = [] # prsp
        for m in mat:
            x.append(ct[m])
            y.append(prsp[m])

        x1 = [] # ct no
        y1 = [] # prsp

        model = np.polyfit(x, y, 1)
        slope.append(model[0])
        const.append(model[1])

        predict = np.poly1d(model)

        x1 = np.array(px[i])
        y1 = np.array(predict(x1))

        rsv = r2(predict(x), y)

        ax.plot(x1, y1,  linestyle = linestyles[i], linewidth= 1, color = [0.5, 0.5, 0.5], \
            label = "%5.9f x + % 5.9f (r2: %5.3f)" % (slope[i], const[i], rsv))

    annot= ax.annotate("", xy = (0,0), xytext = (-20,20), textcoords = "offset points", \
                       bbox = dict(boxstyle = "sawtooth", fc = "w"), \
                       arrowprops = dict(arrowstyle = "->"))
    annot.set_visible(False)

    # loop through each tissue in mean_mea_hu
    vad_hu_t = [t for t in mean_mea_hu.keys()]

    # list of markers and colors used in the plot
    m = ['o', 'p', '*', 'x', '+', 'v', '<', '>', 's', 'd']
    c = ['r','g','b','c','m', 'y', 'k']

    # form arrays of markers and colours according to the number of elements in the vad_hu_t
    ml = [m[i%(len(m))] for i in np.arange(0, len(vad_hu_t), 1) ]
    cl = [c[i%(len(c))] for i in np.arange(0, len(vad_hu_t), 1) ]

    for i, t in enumerate(vad_hu_t):
        ax.plot(mean_mea_hu[t], prsp[t], marker = ml[i], fillstyle = None, color = cl[i], label = t)
        ax.errorbar(mean_mea_hu[t], prsp[t], xerr = std_mea_hu[t],  color = cl[i])

    # print(f"vad_hu_t : {vad_hu_t}")

    def update_annot(ind):
        # Get the x,y coordinates of the selected data point
        pos = pt.get_offsets()[ind["ind"][0]]
        # print(f"pt.get_offsets(): {pt.get_offsets()}") # get the x, y coordinates of all datapoints
        # print(ind["ind"][0], type(ind["ind"][0]) ) # ind["ind"] is an array which contains the index number of the selected datapoint
        annot.xy = pos

        point_x = list(ct.values())[ind["ind"][0]]
        point_y = list(prsp.values())[ind["ind"][0]]
        key_index =  list(ct.keys())[list(ct.values()).index(point_x)]
        text = "{} \n HU:{:.4f} \n prsp:{:.4f})".format(key_index, point_x, point_y)
        # print(f"text: {text}")

        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.6)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = pt.contains(event) # cont = True/ False (is it a data point?) || ind = datapoint index number

            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("button_press_event", hover)
    ax.set_xlabel("HU (a.u.)")
    ax.set_ylabel("prsp (a.u.)")
    ax.set_title("insert= %s || structure = %s  || tissue composition: %s" % (insert_name, structure, ra.tc_source), fontsize =8)
    ax.legend()

    plt.tight_layout()
    plt.show()

    return slope, const
