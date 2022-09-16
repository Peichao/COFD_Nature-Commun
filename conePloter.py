import numpy as np
import pandas as pd
import scipy as sp
from scipy import signal
import scipy.optimize as opt
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import scipy.io as sio
import glob
import os
import csv
import pickle
import math
import hdf5storage
import h5py
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch


# Add details about imaging experiment (unit, expt, plane)
disk = ''
animal = ''  
units = []  

expt_hartley = ['L-cone','M-cone','S-cone','Achromatic']
planeId = ['000','001'] #'001'   #'000'
thres = 0.25
hueaucThres = 0.8

numberofColor = 12
hueList = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
colorSpace = 'HSL'

plot_dominance = False
plot_onoff = True
singleplane = True
plot_diamond = False
white_BK = True    # True for white/gray background, False for aligned 2P image background
BK_alpha = 0.5     # background alpha 
font_size = 15
circle_size = 3
ylimt_hist = 80

## Define the color map for ploting
colormappath = ''
ori_lut = ListedColormap(sio.loadmat(os.path.join(colormappath, 'ori_lut_alpha0.mat'), squeeze_me=True, struct_as_record=False)['lut'])
# cmap_patch = plt.cm.get_cmap('hsv')
cmap_patch_sf = plt.cm.get_cmap('jet')
cmap_patch_cpi = plt.cm.get_cmap('jet')
cmap_patch_osi = plt.cm.get_cmap('jet')
cmap_patch_uniquehue = ListedColormap(hdf5storage.loadmat(os.path.join(colormappath,'uniquehue.mat'), squeeze_me=True, struct_as_record=False)['lut'])
cmap_patch_hsl = ListedColormap(hdf5storage.loadmat(os.path.join(colormappath,'cm_hsl_mshue_l0.4.mat'), squeeze_me=True, struct_as_record=False)['colors'])
cmap_patch_dkl = ListedColormap(hdf5storage.loadmat(os.path.join(colormappath,'cm_dkl_mcchue_l0.mat'), squeeze_me=True, struct_as_record=False)['colors'])
cmap_patch_lidkl = ListedColormap(hdf5storage.loadmat(os.path.join(colormappath,'cm_lidkl_mcchue_l0.mat'), squeeze_me=True, struct_as_record=False)['colors'])
if colorSpace  == 'DKL':
    color_hex_keys = ListedColormap(['#FF8080', '#FF80BF', '#FF80FF', '#BF80FF', '#8080FF', '#80BFFF', '#80FFFF',
            '#80FFBF', '#80FF80', '#BFFF80', '#FFFF80', '#FFBF80', '#808080'])   # DKL hues
elif colorSpace == 'HSL':
    color_hex_keys = ListedColormap(['#af1600', '#8a4600', '#5a5d01', '#2a6600', '#006a00', '#006931', '#006464',
            '#0058b6', '#002DFF', '#6a2ade', '#97209b', '#aa1c50', '#808080'])   # HSL hues
gray = '#B4B4B4'    # '#808080'


main_path = os.path.join(disk, animal, '2P_analysis')
featureData_path = os.path.join(main_path, 'Summary', 'DataExport')  # has ROIs,Bkgs, and data organized in feature dimenstion
result_path_allunit0 = os.path.join(main_path, 'Summary', 'Plots')
result_path_allunit = os.path.join(result_path_allunit0,'STA')

# Make folders for saving data from all units
if not os.path.exists(result_path_allunit0):
    os.mkdir(result_path_allunit0)
if not os.path.exists(result_path_allunit):
    os.mkdir(result_path_allunit)

## Plot cone dominance map
if plot_dominance:
    for j, unit in enumerate(units):
        print('Processing Unit %s.' % unit)
        
        unit_path = os.path.join(main_path, 'U%s' % unit)
        planeData_path = os.path.join(unit_path, '_Summary', 'DataExport')  # has data organized in plane dimension
        result_path_plane0 = os.path.join(unit_path, '_Summary', 'Multiplane')
        result_path_multiplane = os.path.join(result_path_plane0, '0. Original maps')

        # Make folders for saving data from single unit
        if not os.path.exists(result_path_plane0):
            os.mkdir(result_path_plane0)
        if not os.path.exists(result_path_multiplane):
            os.mkdir(result_path_multiplane)

        ## Load cone weights data
        sta_file = os.path.join(planeData_path, '%s_%s_thres%s_sta_dataset.csv' % (animal, unit, thres))
        sta = pd.read_csv(sta_file, index_col=0)
        cellId = sta.cellId.values
        plId = sta.planeId.values
        lmg = sta.lmg.values
        mmg = sta.mmg.values
        smg = sta.smg.values
        ldom = lmg / (lmg + mmg + smg)
        mdom = mmg / (lmg + mmg + smg)
        sdom = smg / (lmg + mmg + smg)
        numCell = np.size(cellId)

        for pl in range(np.size(planeId)):

            plane = planeId[pl]

            # load hue data
            data_file = os.path.join(planeData_path, '%s_%s_%s_sum.csv' % (animal, unit, plane))
            plda = pd.read_csv(data_file, index_col=0)
            xloc = plda.xloc - 1
            yloc = plda.yloc - 1
            print('Processing plane: ' + plane)
            result_path_singleplane = os.path.join(unit_path, '_Summary', 'plane_%s' % (plane), '0. Original maps')
            file_name_singleplane = os.path.join('%s_u%s_plane%s_thres%s_singleplane' % (animal, unit, plane, thres))
            ## load background images and ROIs
            # I use images from orisf_achromatic data as background image. All maps are plotted on this image.
            roibkg_file = os.path.join(planeData_path, '%s_%s_%s_roibkg.jld2' % (animal, unit, plane))
            roibkg = h5py.File(roibkg_file, 'r+')
            bkgimg = np.transpose(roibkg['bkg'][()])
            whitBK = np.ones(bkgimg.shape) * 0.8

            # load roi segements/contours
            rois = roibkg["roi"][()]

            # Transform roi contours to patches
            allROIs = []  # patches of all ROIs
            for i in range(np.size(rois)):
                allROIs.append(0)  # append 0 if there is ROI
                # Find the segment/contour for this cell
                roibkg[rois[i]][()] = roibkg[rois[i]][()].astype(int)
                allROIs[i] = Polygon([tuple(l) for l in list(np.transpose(roibkg[rois[i]][()] - 1))])

            ## Initialize plotting
            fig, ax = plt.subplots()
            if white_BK:
                ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
            else:
                ax.imshow(bkgimg, cmap='gray')

            for vi in np.arange(numCell):

                if plId[vi] == int(plane):
                    ax.add_patch(PolygonPatch(allROIs[cellId[vi] - 1], alpha=1, color=(ldom[vi],mdom[vi],sdom[vi])))
                    # ax.text(xloc[cellId[vi]-1],yloc[cellId[vi]-1], str(cellId[vi]),horizontalalignment='center', fontsize=3)   # add cellId

            ax.set_rasterized(True)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.minorticks_off()
            ax.set_frame_on(False)
            plt.grid(True)
            plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.3)

            # plt.axis('off')
            plt.savefig('%s/%s_hartley_ConeDominance.svg' % (result_path_singleplane, file_name_singleplane), dpi=300,
                        bbox_inches='tight', pad_inches=0, format='svg')
            plt.savefig('%s/%s_hartley_ConeDominance.png' % (result_path_singleplane, file_name_singleplane), dpi=300,
                        bbox_inches='tight', pad_inches=0, format='png')
            # plt.show()




## plotting ON/OFF map
if plot_onoff:
    for j, unit in enumerate(units):
        print('Processing Unit %s.' % unit)

        unit_path = os.path.join(main_path, 'U%s' % unit)
        planeData_path = os.path.join(unit_path,'_Summary', 'DataExport') # has data organized in plane dimension
        result_path_plane0 = os.path.join(unit_path,'_Summary', 'Multiplane')
        result_path_multiplane = os.path.join(result_path_plane0,'0. Original maps')

        # Make folders for saving data from single unit
        if not os.path.exists(result_path_plane0):
            os.mkdir(result_path_plane0)
        if not os.path.exists(result_path_multiplane):
            os.mkdir(result_path_multiplane)

        ## Load cone weights data
        sta_file = os.path.join(planeData_path,'%s_%s_thres%s_sta_dataset.csv' %(animal, unit, thres))
        sta = pd.read_csv(sta_file, index_col=0)

        cellId = sta.cellId.values
        plId = sta.planeId.values
        # lcw = sta.lcw.values
        # mcw = sta.mcw.values
        # scw = sta.scw.values
        # achResp = sta.achroResp.values

        lsign = sta.lsign.values
        msign = sta.msign.values
        ssign = sta.ssign.values
        asign = sta.asign.values

        numCell = np.size(cellId)

        for i, stim in enumerate(expt_hartley):
            print('Processing %s stimulus.' % stim)
            file_name_cone = os.path.join('%s_u%s_%s_thres%s_allplane' % (animal, unit, stim, thres))
            if stim == 'L-cone':
                cs = lsign
            elif stim == 'M-cone':
                cs = msign
            elif stim == 'S-cone':
                cs = ssign
            elif stim == 'Achromatic':
                cs = asign

            
            if singleplane:

                for pl in range(np.size(planeId)):
                    
                    plane = planeId[pl]

                    # load hue data
                    data_file = os.path.join(planeData_path,'%s_%s_%s_sum.csv' %(animal, unit, plane))
                    # data_file = os.path.join(planeData_path, '%s_%s_%s_oriData.csv' % (animal, unit, plane))
                    plda = pd.read_csv(data_file, index_col=0)
                    visResp = plda.visResp
                    # xloc=plda.xloc-1
                    # yloc=plda.yloc-1
                    print('Processing plane: '+plane)
                    result_path_singleplane = os.path.join(unit_path,'_Summary', 'plane_%s' %(plane), '0. Original maps')
                    file_name_singleplane = os.path.join('%s_u%s_plane%s_%s_thres%s_singleplane' % (animal, unit, plane, stim, thres))
                    ## load background images and ROIs
                    # I use images from orisf_achromatic data as background image. All maps are plotted on this image.
                    roibkg_file = os.path.join(planeData_path,'%s_%s_%s_roibkg.jld2' %(animal, unit, plane))
                    # roibkg_file = os.path.join(planeData_path, 'AE6_002_011_roibkg.jld2')   # Sometimes, there was z-shift between Hartley and ori/hue, so the ROIs are different.
                    roibkg = h5py.File(roibkg_file, 'r+')
                    bkgimg = np.transpose(roibkg['bkg'][()])
                    whitBK = np.ones(bkgimg.shape)*0.8

                    # load roi segements/contours
                    rois = roibkg["roi"][()]

                    # Transform roi contours to patches
                    allROIs = []   # patches of all ROIs
                    for i in range(np.size(rois)):
                        allROIs.append(0)    # append 0 if there is ROI
                        #Find the segment/contour for this cell
                        roibkg[rois[i]][()] = roibkg[rois[i]][()].astype(int)
                        allROIs[i] = Polygon([tuple(l) for l in list(np.transpose(roibkg[rois[i]][()]-1))])
                    
                    ## Initialize plotting
                    fig, ax = plt.subplots()
                    if white_BK:
                        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
                    else:
                        ax.imshow(bkgimg, cmap='gray')

                    for vi in np.arange(numCell):

                        if plId[vi]==int(plane):
                            if visResp[cellId[vi]-1]>0:
                                if cs[vi]>0:  # On
                                    ax.add_patch(PolygonPatch(allROIs[cellId[vi]-1], alpha=1, color= 'r'))
                                    # ax.text(xloc[cellId[vi]-1],yloc[cellId[vi]-1], str(cellId[vi]),horizontalalignment='center', fontsize=3)   # add cellId
                                if cs[vi]<0: # Off
                                    ax.add_patch(PolygonPatch(allROIs[cellId[vi]-1], alpha=1, color= 'b'))
                                    # ax.text(xloc[cellId[vi]-1],yloc[cellId[vi]-1], str(cellId[vi]),horizontalalignment='center', fontsize=3)
                    ax.set_rasterized(True)
                    ax.set_xticklabels([])
                    ax.set_yticklabels([])
                    ax.minorticks_off()
                    ax.set_frame_on(False)
                    plt.grid(True)
                    plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.3)

                    # plt.axis('off')
                    plt.savefig('%s/%s_thres%s_hartley_OnOff.svg' % (result_path_singleplane, file_name_singleplane,thres), dpi=300, bbox_inches='tight', pad_inches=0, format='svg')
                    plt.savefig('%s/%s_thres%s_hartley_OnOff.png' % (result_path_singleplane, file_name_singleplane,thres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
                    # plt.show()

                    # fig, ax = plt.subplots()
                    # ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
                    # legend_elements = [Patch(facecolor='b', edgecolor='b', label=stim+' OFF'),
                    #                 Patch(facecolor='r', edgecolor='r', label=stim+' ON')]
                    # ax.legend(handles=legend_elements, loc=4)
                    # plt.axis('off')
                    # plt.savefig('%s/%s_hartley_OnOff_legend.svg' % (result_path, file_name), dpi=300, format='svg')


            ## Initialize plot for multiple planes plotting
            print('Processing Multiple planes')
            fig, ax = plt.subplots()
            
            for pl in range(np.size(planeId)):
            
                plane = planeId[pl]
                # load hue data
                data_file = os.path.join(planeData_path,'%s_%s_%s_sum.csv' %(animal, unit, plane))
                plda = pd.read_csv(data_file, index_col=0)
                visResp = plda.visResp
                ## load background images and ROIs
                # I use images from orisf_achromatic data as background image. All maps are plotted on this image.
                roibkg_file = os.path.join(planeData_path,'%s_%s_%s_roibkg.jld2' %(animal, unit, plane))
                roibkg = h5py.File(roibkg_file, 'r+')
                bkgimg = roibkg['bkg'][()]
                whitBK = np.ones(np.transpose(bkgimg).shape)*0.8
            
                # load roi segements/contours
                rois = roibkg["roi"][()]
            
                # Transform roi contours to patches
                allROIs = []   # patches of all ROIs
                for i in range(np.size(rois)):
                    allROIs.append(0)    # append 0 if there is ROI
                    #Find the segment/contour for this cell
                    roibkg[rois[i]][()] = roibkg[rois[i]][()].astype(int)
                    allROIs[i] = Polygon([tuple(l) for l in list(np.transpose(roibkg[rois[i]][()]-1))])
            
                if white_BK:
                    ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
                else:
                    ax.imshow(bkgimg, cmap='gray')
            
                for vi in np.arange(numCell):
                    if plId[vi]==int(plane):
                        if visResp[cellId[vi]-1]>0:
                            if cs[vi]>0:  # On
                                ax.add_patch(PolygonPatch(allROIs[cellId[vi]-1], alpha=1, color= 'r'))
                            if cs[vi]<0:  # Off
                                ax.add_patch(PolygonPatch(allROIs[cellId[vi]-1], alpha=1, color= 'b'))
            
            ax.set_rasterized(True)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.minorticks_off()
            ax.set_frame_on(False)
            plt.grid(True)
            plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.3)
            
            # plt.axis('off')
            plt.savefig('%s/%s_thres%s_hartley_OnOff.svg' % (result_path_multiplane, file_name_cone,thres), dpi=300, bbox_inches='tight', pad_inches=0, format='svg')
            plt.savefig('%s/%s_thres%s_hartley_OnOff.png' % (result_path_multiplane, file_name_cone,thres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
            # plt.show()

            # fig, ax = plt.subplots()
            # ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
            # legend_elements = [Patch(facecolor='b', edgecolor='b', label=stim+' OFF'),
            #                 Patch(facecolor='r', edgecolor='r', label=stim+' ON')]
            # ax.legend(handles=legend_elements, loc=4)
            # plt.axis('off')
            # plt.savefig('%s/%s_hartley_OnOff_legend.svg' % (result_path, file_name), dpi=300, format='svg')


# #########
if plot_diamond:

    print('Plot Cone Weights of cells have LMS STA......')
    file_name_cw = os.path.join('%s_u%s_thres%s_allplane_cone weight' % (animal, ''.join(units), thres))
    fig, ax = plt.subplots()
    STA_DKLcell = pd.DataFrame(columns=['animal', 'unit', 'plane', 'cellId', 'HueSelec', 'lcw', 'mcw', 'maxhue'])
    cellNum = 0
    Countcell = 0
    for j, unit in enumerate(units):
        print('Processing Unit %s.' % unit)

        unit_path = os.path.join(main_path, 'U%s' % unit)
        planeData_path = os.path.join(unit_path,'_Summary', 'DataExport') # has data organized in plane dimension

        ## Load cone weights data
        sta_file = os.path.join(planeData_path,'%s_%s_thres%s_sta_dataset.csv' %(animal, unit, thres))
        sta = pd.read_csv(sta_file, index_col=0)

        cellId = sta.cellId.values
        plId = sta.planeId.values
        lcw = sta.lcwmn.values
        mcw = sta.mcwmn.values
        scw = sta.scwmn.values
        isl = sta.isl.values
        ism = sta.ism.values
        iss = sta.iss.values

        # lcw = sta.lcwmg.values
        # mcw = sta.mcwmg.values
        # scw = sta.scwmg.values

        # lcw = sta.lcwmgall.values
        # mcw = sta.mcwmgall.values
        # scw = sta.scwmgall.values

        for pl in range(np.size(planeId)):

            plane = planeId[pl]
            # load hue data
            data_file = os.path.join(planeData_path,'%s_%s_%s_sum.csv' %(animal, unit, plane))
            plda = pd.read_csv(data_file, index_col=0)
            max_hue = plda.maxhue
            visResp = plda.visResp
            hue_ax_auc = plda.hueaxauc
            hue_di_auc = plda.huediauc

            for i in range(len(cellId)):
                if (plId[i]==int(plane)) & (~(math.isnan(lcw[i])) | ~(math.isnan(mcw[i]))):
                    if visResp[cellId[i]-1]>0:
                        if ((hue_ax_auc[cellId[i]-1] >hueaucThres) | (hue_di_auc[cellId[i]-1]>hueaucThres)): #& (isl[i]) & (ism[i]):
                        # if  (isl[i]) & (ism[i]):
                            STA_DKLcell.loc[cellNum] = [animal, unit, plane, cellId[i], True, lcw[i], mcw[i], max_hue[cellId[i]-1]]
                            Countcell +=1
                            if colorSpace == 'DKL':
                                plt.scatter(lcw[i], mcw[i], s=circle_size, alpha=1, color=cmap_patch_lidkl(max_hue[cellId[i]-1]/360))
                            elif colorSpace == 'HSL':
                                plt.scatter(lcw[i], mcw[i], s=circle_size, alpha=1, color=color_hex_keys(hueList.index(max_hue[cellId[i] - 1])))
                        # elif ((hue_ax_auc[cellId[i]-1] >hueaucThres) | (hue_di_auc[cellId[i]-1]>hueaucThres)) & (~ism[i] | ~isl[i]):
                        #     if colorSpace == 'DKL':
                        #         plt.scatter(lcw[i], mcw[i], s=circle_size, alpha=0.2, color=cmap_patch_lidkl(max_hue[cellId[i]-1]/360))
                        #     elif colorSpace == 'HSL':
                        #         plt.scatter(lcw[i], mcw[i], s=circle_size, alpha=0.2, color=color_hex_keys(hueList.index(max_hue[cellId[i] - 1])))
                        else:
                            STA_DKLcell.loc[cellNum] = [animal, unit, plane, cellId[i], False, lcw[i], mcw[i], 999]
                            plt.scatter(lcw[i], mcw[i], s=3, alpha=0.4, color='gray')

                        cellNum += 1

    plt.plot([0, 0], [1, 0], 'k-', lw=0.5, linestyle='--')
    plt.plot([0, 0], [-1, 0], 'k-', lw=0.5, linestyle='--')
    plt.plot([0, -1], [0, 0], 'k-', lw=0.5, linestyle='--')
    plt.plot([0, 1], [0, 0], 'k-', lw=0.5, linestyle='--')
    plt.plot([-1, 0], [0, 1], 'k-', lw=0.5, linestyle='-')
    plt.plot([-1, 0], [0, -1], 'k-', lw=0.5, linestyle='-')
    plt.plot([1, 0], [0, 1], 'k-', lw=0.5, linestyle='-')
    plt.plot([1, 0], [0, -1], 'k-', lw=0.5, linestyle='-')

    plt.ylabel('Normalized M-cone Strength', fontsize=font_size)
    plt.xlabel('Normalized L-cone Strength', fontsize=font_size)
    plt.ylim(-1, 1)
    plt.xlim(-1, 1)
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    ax.set_aspect('equal', 'box')

    plt.savefig('%s/%s_CellNum%s_mn_aucThres%s_L&M&S%s_LMScell.svg' % (result_path_allunit, file_name_cw,cellNum,hueaucThres,Countcell), dpi=300, format='svg')
    plt.savefig('%s/%s_CellNum%s_mn_aucThres%s_L&M&S%s_LMScell.png' % (result_path_allunit, file_name_cw,cellNum,hueaucThres,Countcell), dpi=300, bbox_inches='tight', pad_inches=0,format='png')
    # plt.show()
    STA_DKLcell.to_csv('%s/%s_totalcellNum%s__huecellNum%s_aucThres%s_STADKL_List.csv' % (result_path_allunit, file_name_cw,cellNum, Countcell,hueaucThres))


    #### PLot cone weight of cells have strong S-on inputs
    print('Plot Cone Weights of cells if have S cone STA......')
    fig, ax = plt.subplots()
    cellNum = 0
    for j, unit in enumerate(units):
        print('Processing Unit %s.' % unit)

        unit_path = os.path.join(main_path, 'U%s' % unit)
        planeData_path = os.path.join(unit_path,'_Summary', 'DataExport') # has data organized in plane dimension

        ## Load cone weights data
        sta_file = os.path.join(planeData_path,'%s_%s_thres%s_sta_dataset.csv' %(animal, unit, thres))
        sta = pd.read_csv(sta_file, index_col=0)

        cellId = sta.cellId.values
        plId = sta.planeId.values
        lcw = sta.lcwmn.values
        mcw = sta.mcwmn.values
        scw = sta.scwmn.values
        iss = sta.iss.values

        # lcw = sta.lcwmg.values
        # mcw = sta.mcwmg.values
        # scw = sta.scwmg.values

        # lcw = sta.lcwmgall.values
        # mcw = sta.mcwmgall.values
        # scw = sta.scwmgall.values

        for pl in range(np.size(planeId)):

            plane = planeId[pl]

            data_file = os.path.join(planeData_path,'%s_%s_%s_sum.csv' %(animal, unit, plane))
            plda = pd.read_csv(data_file, index_col=0)
            max_hue = plda.maxhue
            visResp = plda.visResp
            hue_ax_auc = plda.hueaxauc
            hue_di_auc = plda.huediauc

            for i in range(len(cellId)):
                # if (plId[i]==int(plane)) & (np.abs(lcw[i])+np.abs(mcw[i])<0.5) & (scw[i] >= 0):
                if (plId[i] == int(plane)) & (iss[i]) & (scw[i] >= 0):
                    if visResp[cellId[i]-1]>0:
                        cellNum += 1
                        if ((hue_ax_auc[cellId[i]-1] >hueaucThres) | (hue_di_auc[cellId[i]-1]>hueaucThres)):
                            if colorSpace == 'DKL':
                                plt.scatter(lcw[i], mcw[i], s=circle_size, alpha=1, color=cmap_patch_lidkl(max_hue[cellId[i]-1]/360))
                            elif colorSpace == 'HSL':
                                plt.scatter(lcw[i], mcw[i], s=circle_size, alpha=1, color=color_hex_keys(hueList.index(max_hue[cellId[i]-1])))
                        else:
                            plt.scatter(lcw[i], mcw[i], s=3, alpha=0.4, color='gray')

    plt.plot([0, 0], [1, 0], 'k-', lw=0.5, linestyle='--')
    plt.plot([0, 0], [-1, 0], 'k-', lw=0.5, linestyle='--')
    plt.plot([0, -1], [0, 0], 'k-', lw=0.5, linestyle='--')
    plt.plot([0, 1], [0, 0], 'k-', lw=0.5, linestyle='--')
    plt.plot([-1, 0], [0, 1], 'k-', lw=0.5, linestyle='-')
    plt.plot([-1, 0], [0, -1], 'k-', lw=0.5, linestyle='-')
    plt.plot([1, 0], [0, 1], 'k-', lw=0.5, linestyle='-')
    plt.plot([1, 0], [0, -1], 'k-', lw=0.5, linestyle='-')

    plt.ylabel('Normalized M-cone Strength', fontsize=font_size)
    plt.xlabel('Normalized L-cone Strength', fontsize=font_size)
    plt.ylim(-1, 1)
    plt.xlim(-1, 1)
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    ax.set_aspect('equal', 'box')

    plt.savefig('%s/%s_cellNum%s_mn_aucThres%s_Son.svg' % (result_path_allunit, file_name_cw, cellNum,hueaucThres), dpi=300, format='svg')
    plt.savefig('%s/%s_cellNum%s_mn_aucThres%s_Son.png' % (result_path_allunit, file_name_cw, cellNum,hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0,format='png')
    # plt.show()

    #### PLot cone weight of cells have strong S-off inputs
    print('Plot Cone Weights of cells if have S cone STA......')
    fig, ax = plt.subplots()
    cellNum = 0
    for j, unit in enumerate(units):
        print('Processing Unit %s.' % unit)

        unit_path = os.path.join(main_path, 'U%s' % unit)
        planeData_path = os.path.join(unit_path, '_Summary', 'DataExport')  # has data organized in plane dimension

        ## Load cone weights data
        sta_file = os.path.join(planeData_path, '%s_%s_thres%s_sta_dataset.csv' % (animal, unit, thres))
        sta = pd.read_csv(sta_file, index_col=0)

        cellId = sta.cellId.values
        plId = sta.planeId.values
        lcw = sta.lcwmn.values
        mcw = sta.mcwmn.values
        scw = sta.scwmn.values
        iss = sta.iss.values

        # lcw = sta.lcwmg.values
        # mcw = sta.mcwmg.values
        # scw = sta.scwmg.values

        # lcw = sta.lcwmgall.values
        # mcw = sta.mcwmgall.values
        # scw = sta.scwmgall.values

        for pl in range(np.size(planeId)):

            plane = planeId[pl]

            data_file = os.path.join(planeData_path, '%s_%s_%s_sum.csv' % (animal, unit, plane))
            plda = pd.read_csv(data_file, index_col=0)
            max_hue = plda.maxhue
            visResp = plda.visResp
            hue_ax_auc = plda.hueaxauc
            hue_di_auc = plda.huediauc

            for i in range(len(cellId)):
                # if (plId[i]==int(plane)) & (np.abs(lcw[i])+np.abs(mcw[i])<0.5) & (scw[i] < 0):
                if (plId[i] == int(plane)) & (iss[i]) & (scw[i] < 0):
                    if visResp[cellId[i]-1]>0:
                        cellNum += 1
                        if ((hue_ax_auc[cellId[i]-1] >hueaucThres) | (hue_di_auc[cellId[i]-1]>hueaucThres)):
                            if colorSpace == 'DKL':
                                plt.scatter(lcw[i], mcw[i], s=circle_size, alpha=1, facecolors='none', edgecolors=cmap_patch_lidkl(max_hue[cellId[i]-1]/360))
                            elif colorSpace == 'HSL':
                                plt.scatter(lcw[i], mcw[i], s=circle_size, alpha=1, facecolors='none', edgecolors=color_hex_keys(hueList.index(max_hue[cellId[i]-1])))
                        else:
                            plt.scatter(lcw[i], mcw[i], s=3, alpha=0.4, facecolors='none', edgecolors='gray')

    plt.plot([0, 0], [1, 0], 'k-', lw=0.5, linestyle='--')
    plt.plot([0, 0], [-1, 0], 'k-', lw=0.5, linestyle='--')
    plt.plot([0, -1], [0, 0], 'k-', lw=0.5, linestyle='--')
    plt.plot([0, 1], [0, 0], 'k-', lw=0.5, linestyle='--')
    plt.plot([-1, 0], [0, 1], 'k-', lw=0.5, linestyle='-')
    plt.plot([-1, 0], [0, -1], 'k-', lw=0.5, linestyle='-')
    plt.plot([1, 0], [0, 1], 'k-', lw=0.5, linestyle='-')
    plt.plot([1, 0], [0, -1], 'k-', lw=0.5, linestyle='-')

    plt.ylabel('Normalized M-cone Strength', fontsize=font_size)
    plt.xlabel('Normalized L-cone Strength', fontsize=font_size)
    plt.ylim(-1, 1)
    plt.xlim(-1, 1)
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    ax.set_aspect('equal', 'box')

    plt.savefig('%s/%s_cellNum%s_mn_aucThres%s_Soff.svg' % (result_path_allunit, file_name_cw, cellNum,hueaucThres), dpi=300, format='svg')
    plt.savefig('%s/%s_cellNum%s_mn_aucThres%s_Soff.png' % (result_path_allunit, file_name_cw, cellNum,hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    # plt.show()


    #### PLot Hue hist of cells have strong S inputs
    print('Plot Hue hists of cells if have S cone STA......')

    sOnCellNum = 0
    sOffCellNum = 0
    scount = 0
    sOnCell = []
    sOffCell = []
    Scell = pd.DataFrame(columns=['animal', 'unit', 'plane', 'cellId', 'sign', 'maxhue'])
    for j, unit in enumerate(units):
        print('Processing Unit %s.' % unit)

        unit_path = os.path.join(main_path, 'U%s' % unit)
        planeData_path = os.path.join(unit_path,'_Summary', 'DataExport') # has data organized in plane dimension

        ## Load cone weights data
        sta_file = os.path.join(planeData_path,'%s_%s_thres%s_sta_dataset.csv' %(animal, unit, thres))
        sta = pd.read_csv(sta_file, index_col=0)

        cellId = sta.cellId.values
        plId = sta.planeId.values
        lcw = sta.lcwmn.values
        mcw = sta.mcwmn.values
        scw = sta.scwmn.values
        iss = sta.iss.values

        # lcw = sta.lcwmg.values
        # mcw = sta.mcwmg.values
        # scw = sta.scwmg.values

        # lcw = sta.lcwmgall.values
        # mcw = sta.mcwmgall.values
        # scw = sta.scwmgall.values

        for pl in range(np.size(planeId)):

            plane = planeId[pl]

            data_file = os.path.join(planeData_path,'%s_%s_%s_sum.csv' %(animal, unit, plane))
            plda = pd.read_csv(data_file, index_col=0)
            max_hue = plda.maxhue
            visResp = plda.visResp
            hue_ax_auc = plda.hueaxauc
            hue_di_auc = plda.huediauc

            for i in range(len(cellId)):
                # if (plId[i] == int(plane)) & (np.abs(lcw[i]) + np.abs(mcw[i]) < 0.5):
                if (plId[i] == int(plane)) & (iss[i]):
                    if ((hue_ax_auc[cellId[i]-1] > hueaucThres) | (hue_di_auc[cellId[i]-1] > hueaucThres)) & (visResp[cellId[i]-1] > 0):
                        if scw[i] >= 0:  # S-on
                            sOnCellNum += 1
                            sOnCell.append(max_hue[cellId[i] - 1])
                        elif scw[i] < 0:  # S-off
                            sOffCellNum += 1
                            sOffCell.append(max_hue[cellId[i] - 1])
                        Scell.loc[scount] = [animal, unit, plane, cellId[i], np.sign(scw[i]), max_hue[cellId[i] - 1]]
                        scount += 1
    Scell.to_csv('%s/%s_cellNum%s_aucThres%s_Sonoff_List.csv' % (result_path_allunit, file_name_cw,scount,hueaucThres))
    # sOffCell.to_csv('%s/%s_cellNum%s_aucThres%s_Soff_List.csv' % (result_path_allunit, file_name_cw, sOffCellNum, hueaucThres))
    ## Son Hist
    fig, ax = plt.subplots()
    histPatch = pd.Series(sOnCell).plot(kind='hist', bins=12, rwidth=0.8)
    # Set colors on hist patches
    for patches in range(0, 12):
        if colorSpace == 'DKL':
            histPatch.patches[patches].set_color(cmap_patch_lidkl(patches/12))
        elif colorSpace == 'HSL':
            histPatch.patches[patches].set_color(color_hex_keys(patches))
        # histPatch.patches[patches].set_color(cmap_patch_hue1(patches/36))


    # ax.axes.get_xaxis().set_visible(False)
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        labelsize='large',
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off
    plt.ylabel('Cell Number', fontsize=font_size)
    plt.xlabel('Preferred Hue', fontsize=font_size)
    plt.ylim(0, ylimt_hist)  # 380
    # plt.grid(b=True, which='major', axis='y', alpha = 0.5, linewidth=1)
    ax.set_axisbelow(True)

    plt.savefig('%s/%s_cellNum%s_aucThres%s_Son_hist.svg' % (result_path_allunit, file_name_cw, sOnCellNum,hueaucThres),dpi=300, format='svg')
    plt.savefig('%s/%s_cellNum%s_aucThres%s_Son_hist.png' % (result_path_allunit, file_name_cw, sOnCellNum,hueaucThres),dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    # plt.show()


    ## Soff hist
    fig, ax = plt.subplots()
    histPatch = pd.Series(sOffCell).plot(kind='hist', bins=12, rwidth=0.8)
    # Set colors on hist patches
    for patches in range(0, 12):
        if colorSpace == 'DKL':
            histPatch.patches[patches].set_color(cmap_patch_lidkl(patches/12))
        elif colorSpace == 'HSL':
            histPatch.patches[patches].set_color(color_hex_keys(patches))

    # ax.axes.get_xaxis().set_visible(False)
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        labelsize='large',
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off
    plt.ylabel('Cell Number', fontsize=font_size)
    plt.xlabel('Preferred Hue', fontsize=font_size)
    plt.ylim(0, ylimt_hist)  # 380
    # plt.grid(b=True, which='major', axis='y', alpha = 0.5, linewidth=1)
    ax.set_axisbelow(True)

    plt.savefig('%s/%s_cellNum%s_aucThres%s_Soff_hist.svg' % (result_path_allunit, file_name_cw, sOffCellNum,hueaucThres),dpi=300, format='svg')
    plt.savefig('%s/%s_cellNum%s_aucThres%s_Soff_hist.png' % (result_path_allunit, file_name_cw, sOffCellNum,hueaucThres), dpi=300,bbox_inches='tight', pad_inches=0, format='png')
    # plt.show()
