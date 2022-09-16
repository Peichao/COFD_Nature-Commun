import numpy as np
import scipy.io as sio
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import pickle
import os
import glob
import warnings
import hdf5storage
import h5py
warnings.simplefilter('error', RuntimeWarning)

# Add details about imaging experiment (unit, expt, plane)
disk = 'G:'
animal = 'AF9'   # AF4, AF3
unit =   '004' #'004'  # '003'
expt_orisf = '003'  # '005' # '001'
expt_hue = '000'   # '008'   # '009'
plane = '001'  #'001'   #'000'

# Check which cell is selective
oriaucThres = 0.7
diraucThres = 0.7
hueaucThres = 0.7

# Selectivity defined by circular variance
# oricvThres = 0.7
# dircvThres = 0.7
# huecvThres = 0.7

# Selectivity defined by selective index
# osiThres = 0.33
# dsiThres = 0.33
# huesiThres = 0.33
huecpiThres = -1  # -1: plot all, 0: more prefer hue, 0.33: strongly tuning to hue

numberofColor = 12
hueList = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
colorSpace = 'HSL'

plot_huehist = False
plot_ori = True
plot_hue = True
# plot_ori = False
# plot_hue = False

white_BK = True    # True for white/gray background, False for aligned 2P image background
BK_alpha = 0.5     # background alpha 
font_size = 15
ylimt_hist = 400

main_path = os.path.join(disk, animal, '2P_analysis')
featureData_path = os.path.join(main_path, 'Summary', 'DataExport')  # has ROIs,Bkgs, and data organized in feature dimenstion
unit_path = os.path.join(main_path, 'U%s' % unit)
planeData_path = os.path.join(unit_path,'_Summary', 'DataExport') # has data organized in plane dimension

expt_path_hue = os.path.join(unit_path, '%s_%s_%s' % (unit, expt_hue, plane), 'DataExport') # only has hue data
expt_path_ori = os.path.join(unit_path, '%s_%s_%s' % (unit, expt_orisf, plane),'DataExport') # only has ori/dir data
result_path_plane0 = os.path.join(unit_path,'_Summary', 'plane_%s'%plane)

result_path_plane = os.path.join(result_path_plane0, '0. Original maps')
result_path_hue = os.path.join(unit_path, '%s_%s_%s' % (unit, expt_hue, plane), 'Plots') # only has hue data
result_path_ori = os.path.join(unit_path, '%s_%s_%s' % (unit, expt_orisf, plane),'Plots') # only has ori/dir data

file_name_ori = os.path.join('%s_%s_%s_%s_ori' % (animal, unit, expt_orisf, plane))
file_name_sf = os.path.join('%s_%s_%s_%s_sf' % (animal, unit, expt_orisf, plane))
file_name_dir = os.path.join('%s_%s_%s_%s_dir' % (animal, unit, expt_orisf, plane))

file_name_hueori = os.path.join('%s_%s_%s_%s_hueori' % (animal, unit, expt_hue, plane))
file_name_huedir = os.path.join('%s_%s_%s_%s_huedir' % (animal, unit, expt_hue, plane))

file_name_hue = os.path.join('%s_%s_%s_%s_hue' % (animal, unit, expt_hue, plane))

roi_name = os.path.join('%s_u%s_plane%s' % (animal, unit, plane))

# Make folders
if not os.path.exists(result_path_plane0):
    os.mkdir(result_path_plane0)

if not os.path.exists(result_path_plane):
    os.mkdir(result_path_plane)

if (not os.path.exists(result_path_hue)) & plot_hue:
    os.mkdir(result_path_hue)
resultCell_path_hue = os.path.join(result_path_hue, 'cells')
if (not os.path.exists(resultCell_path_hue)) & plot_hue:
    os.mkdir(resultCell_path_hue)

if not os.path.exists(result_path_ori):
    os.mkdir(result_path_ori)
resultCell_path_ori = os.path.join(result_path_ori, 'cells')
if not os.path.exists(resultCell_path_ori):
    os.mkdir(resultCell_path_ori)


## Define the color map for ploting
colormappath = 'C:\\Users\\lyr19\\.julia\\dev\\NeuroAnalysis\\src\\Visualization\\colors'
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


## load background images and ROIs
# I use images from orisf_achromatic data as background image. All maps are plotted on this image.
roibkg_file = os.path.join(planeData_path,'%s_%s_%s_roibkg.jld2' %(animal, unit, plane))
# roibkg_file = os.path.join(expt_path_hue,'%s_%s_%s_%s_roibkg.jld2' %(animal, unit, expt_hue, plane))
roibkg = h5py.File(roibkg_file, 'r+')
bkgimg = roibkg['bkg'][()]
whitBK = np.ones(np.transpose(bkgimg).shape)*0.8

# load roi segements/contours
rois = roibkg["roi"][()]
# obj = roibkg[objref[0]]
# rois = obj.value
# Transform roi contours to patches
allROIs = []   # patches of all ROIs
for i in range(np.size(rois)):
    allROIs.append(0)    # append 0 if there is ROI
    #Find the segment/contour for this cell
    roibkg[rois[i]][()] = roibkg[rois[i]][()].astype(int) 
    allROIs[i] = Polygon([tuple(l) for l in list(np.transpose(roibkg[rois[i]][()]-1))])

## Load Result files
data_file = os.path.join(planeData_path,'%s_%s_%s_sum.csv' %(animal, unit, plane))
# data_file = os.path.join(result_path_ori,'%s_%s_%s_%s_result.csv' %(animal, unit, expt_orisf, plane))
plda = pd.read_csv(data_file, index_col=0)
# exp_params_data = hdf5storage.loadmat(exp_params)
# plane_data = h5py.File(data_file, 'r+')
# objref = plane_data[unit+'_'+plane][()]
# obj = [plane_data[ref] for ref in objref]
# # obj = plane_data[objref[0]]
# plda = obj.value

# plane_data_sum = h5py.File(data_file, 'r+')[unit+'_'+plane+'_'+'sum']
# objref = plane_data_sum[()]
# obj = plane_data_sum[objref[0]]
# plda_sum = obj.value

# visual responsive
visResp =  plda.visResp
cellId = plda.cellId
numCell = np.size(cellId)
if plot_hue:
    # hue selective
    hue_ax_auc = plda.hueaxauc
    hue_di_auc = plda.huediauc

    hue_ax_cv = plda.hueaxcv
    hue_di_cv = plda.huedicv

    cv_hue_ax = plda.cvhueax  # from vector summation
    cv_hue_di = plda.cvhuedi
    fit_hue_ax = plda.fithueax  # from fitting
    fit_hue_di = plda.fithuedi
    fit_hue_ax_si = plda.hueaxsi1 # selective index
    fit_hue_di_si = plda.huedisi1
    max_hue = plda.maxhue  # no vector summation or fitting
    max_hue_resp = plda.maxhueresp

    # orientation/direction/sf selective from Hue grating
    hueori_auc = plda.hueoriauc
    huedir_auc = plda.huedirauc

    hueori_cv = plda.hueoricv
    huedir_cv = plda.huedircv

    cv_hueori = plda.cvhueori  # from vector summation
    cv_huedir = plda.cvhuedir
    fit_hueori = plda.fithueori  # from fitting
    fit_huedir = plda.fithuedir
    fit_hueosi = plda.hueosi1
    fit_huedsi = plda.huedsi1

    # fit_huesf = plda.huefitsf
    # xloc = plda.xloc
    # yloc = plda.yloc
# orientation/direction/sf selective from achromatic grating
ori_auc = plda.oriauc
dir_auc = plda.dirauc

ori_cv = plda.oricv  # cv of orientation  (1-mag)
dir_cv = plda.dircv  # cv of direction (1-mag)

cv_ori = plda.cvori  # preferred orientation from vector summation/cv
cv_dir = plda.cvdir # preferred direction from vector summation/cv
fit_ori = plda.fitori  # preferred orientation from cv
fit_dir = plda.fitdir
fit_osi = plda.osi1
fit_dsi = plda.dsi1

# fit_sf = plda.fitsf
# xloc = plda.xloc - 1
# yloc = plda.yloc - 1


## histogram of hue selective cells.
if plot_huehist:
    ## Calculate CPI
    hueResp = max_hue_resp
    hueResp = 1-np.min(pd.concat([hue_ax_cv, hue_di_cv], axis=1),axis=1)
    achroResp = 1-ori_cv
    cpi = (hueResp - achroResp) / (hueResp + achroResp)

    fig, ax = plt.subplots()

    hue_cellNum = np.sum(((hue_ax_auc > hueaucThres) | (hue_di_auc > hueaucThres)) & (visResp > 0))# & (cpi>huecpiThres))
    # print(hue_cellNum)
    hue_hist = max_hue[((hue_ax_auc > hueaucThres) | (hue_di_auc > hueaucThres)) & (visResp > 0)] #& (cpi>huecpiThres)]

    if numberofColor == 12:
        Hist = hue_hist
    elif numberofColor == 13:   # including visual responsive cells as gray color
        visHist = visResp.copy()
        for i in np.arange(numCell):
            if visResp[i] in hue_hist.index.values:
                visHist.remove(visResp[i])
                visual_values =  np.ones(np.size(visHist))*27  # use 27 represent visual responsive cells
                visual_hist =  pd.Series(visual_values, index=visHist)
        Hist = hue_hist.append(visual_hist)

    histPatch = Hist.plot(kind='hist', bins = 12, rwidth=0.8)
    # Set colors on hist patches
    for patches in range(0, 12):
        if colorSpace == 'DKL':
            histPatch.patches[patches].set_color(cmap_patch_lidkl(patches / 12))
        elif colorSpace == 'HSL':
            histPatch.patches[patches].set_color(color_hex_keys(patches))
        # histPatch.patches[patches].set_color(cmap_patch_hue1(patches/36))

    # ax.axes.get_xaxis().set_visible(False)
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off
    plt.ylabel('Cell Number', fontsize = font_size)
    plt.xlabel('Preferred Hue', fontsize = font_size)
    plt.ylim(0, ylimt_hist)  #380
    plt.grid(b=True, which='major', axis='y', alpha = 0.5, linewidth=1)
    ax.set_axisbelow(True)

    # plt.savefig('%s/%s_cellNum%s_cpiThres%s_hist.svg' % (result_path_hue, file_name_hue, str(hue_cellNum), huecpiThres), dpi=300, format='svg')
    plt.savefig('%s/%s_cellNum%s_aucThres%s_cpiThres%s_hist.png' % (result_path_plane, file_name_hue, str(hue_cellNum), hueaucThres, huecpiThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_cellNum%s_aucThres%s_cpiThres%s_hist.png' % (result_path_hue, file_name_hue, str(hue_cellNum), hueaucThres, huecpiThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.show()


### Plot hue and visual responsive maps
if plot_hue:

    ### Plot only hue responsive cells using hue direction gives max response
    fig, ax = plt.subplots()

    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
        # file_name_hue = os.path.join('%s_whtBK%s' % (file_name_hue, str(white_BK)))
    else:
        ax.imshow(bkgimg, cmap='gray')
    for vi in np.arange(numCell):
        if ((hue_ax_auc[vi] >hueaucThres) | (hue_di_auc[vi]>hueaucThres)) & (visResp[vi]>0):
            if colorSpace == 'DKL':
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=cmap_patch_lidkl(max_hue[vi]/360)))
            elif colorSpace == 'HSL':
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=color_hex_keys(hueList.index(max_hue[vi]))))

    # ax.set_title(file_name_hue, fontsize=font_size)
    ax.set_rasterized(True)
    plt.axis('off')
    # plt.savefig('%s/%s_maxhue_auc.svg' % (result_path_hue, file_name_hue), dpi=300, format='svg')
    plt.savefig('%s/%s_maxhue_auc%s.png' % (result_path_plane, file_name_hue, hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_maxhue_auc%s.png' % (result_path_hue, file_name_hue, hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')


    ### Plot hue selective cell using hue direction gives max response
    fig, ax = plt.subplots()
    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
    else:
        ax.imshow(bkgimg, cmap='gray')

    for vi in np.arange(numCell):
        if visResp[vi]>0:
           ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
        #    ax.text(xloc[vi], yloc[vi], str(vi+1), horizontalalignment='center', fontsize=3)    # add cellId if needed
        if ((hue_ax_auc[vi] >hueaucThres) | (hue_di_auc[vi]>hueaucThres)) & (visResp[vi]>0):
            if colorSpace == 'DKL':
                 ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=cmap_patch_lidkl(max_hue[vi]/360))) # np.log10(hue_p[vi])/np.min(np.log10(hue_p))
        #    ax.text(xloc[vi], yloc[vi], str(vi+1), horizontalalignment='center', fontsize=3)
            elif colorSpace == 'HSL':
                 ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=color_hex_keys(hueList.index(max_hue[vi]))))

    # ax.set_title( file_name_hue, fontsize=font_size)
    ax.set_rasterized(True)
    plt.axis('off')
    # plt.savefig('%s/%s_maxhue_auc.svg' % (result_path_hue, file_name_hue), dpi=300, format='svg')
    plt.savefig('%s/%s_maxhue_auc%s_all.png' % (result_path_plane, file_name_hue, hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_maxhue_auc%s_all.png' % (result_path_hue, file_name_hue, hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')


    ### Plot hue axis selective cell using fitting parameters
    # fig, ax = plt.subplots()
    # if white_BK:
    #     ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
    # else:
    #     ax.imshow(bkgimg, cmap='gray')

    # for vi in np.arange(numCell):
    #     if visResp[vi]>0:
    #        ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
    #     # if (hue_ax_auc[vi] >hueaucThres) & (visResp[vi]>0):
    #     if ((hue_ax_auc[vi] >hueaucThres) | (hue_di_auc[vi]>hueaucThres)) & (visResp[vi]>0):
    #         if np.isnan(fit_hue_ax[vi]):
    #             ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=gray))
    #         else:
    #             ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=cmap_patch_hue(fit_hue_ax[vi]/180))) # np.log10(hue_p[vi])/np.min(np.log10(hue_p))

    # # ax.set_title( file_name_hue, fontsize=font_size)
    # ax.set_rasterized(True)
    # plt.axis('off')
    # # plt.savefig('%s/%s_fithueax_auc.svg' % (result_path_hue, file_name_hue), dpi=300, format='svg')
    # plt.savefig('%s/%s_fithueax_auc%s.png' % (result_path_plane, file_name_hue, hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    # plt.savefig('%s/%s_fithueax_auc%s.png' % (result_path_hue, file_name_hue, hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')



    ### Plot hue di selective cell using fitting parameters
    fig, ax = plt.subplots()
    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
    else:
        ax.imshow(bkgimg, cmap='gray')

    for vi in np.arange(numCell):
        if visResp[vi]>0:
           ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
        # if (hue_di_auc[vi]>hueaucThres) & (visResp[vi]>0):
        if ((hue_ax_auc[vi] >hueaucThres) | (hue_di_auc[vi]>hueaucThres)) & (visResp[vi]>0):
            if np.isnan(fit_hue_di[vi]):
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=gray))
            else:
                if colorSpace == 'DKL':
                     ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=cmap_patch_lidkl(fit_hue_di[vi]/360))) # np.log10(hue_p[vi])/np.min(np.log10(hue_p))
                elif colorSpace == 'HSL':
                     ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=color_hex_keys(hueList.index(max_hue[vi]))))
    # ax.set_title( file_name_hue, fontsize=font_size)
    ax.set_rasterized(True)
    plt.axis('off')
    # plt.savefig('%s/%s_fithuedi_auc.svg' % (result_path_hue, file_name_hue), dpi=300, format='svg')
    plt.savefig('%s/%s_fithuedi_auc%s.png' % (result_path_plane, file_name_hue, hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_fithuedi_auc%s.png' % (result_path_hue, file_name_hue, hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')


    ### Plot hue selective cell using fitting parameters (both axis and dir)
    fig, ax = plt.subplots()
    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
    else:
        ax.imshow(bkgimg, cmap='gray')

    for vi in np.arange(numCell):
        if visResp[vi]>0:
           ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
        if (hue_ax_auc[vi] >hueaucThres) & (visResp[vi]>0):
            if np.isnan(fit_hue_ax[vi]):
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=gray))
            else:
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=cmap_patch_lidkl(fit_hue_ax[vi]/180))) # np.log10(hue_p[vi])/np.min(np.log10(hue_p))
        if (hue_di_auc[vi]>hueaucThres) & (visResp[vi]>0):
            if np.isnan(fit_hue_di[vi]):
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=gray))
            else:
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=cmap_patch_lidkl(fit_hue_di[vi]/360))) # np.log10(hue_p[vi])/np.min(np.log10(hue_p))
        
    # ax.set_title( file_name_hue, fontsize=font_size)
    ax.set_rasterized(True)
    plt.axis('off')
    # plt.savefig('%s/%s_fithuedi_auc.svg' % (result_path_hue, file_name_hue), dpi=300, format='svg')
    plt.savefig('%s/%s_fithue_diax_auc%s.png' % (result_path_plane, file_name_hue, hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_fithue_diax_auc%s.png' % (result_path_hue, file_name_hue, hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')



    ### Plot CPI map of all visual responsive cells 

    ## Calculate CPI
    # hueResp = max_hue_resp
    # hueResp = 1-np.min(pd.concat([hue_ax_cv, hue_di_cv], axis=1),axis=1)
    # achroResp = 1-ori_cv
    # cpi = (hueResp - achroResp) / (hueResp + achroResp)

    # fig, ax = plt.subplots()
    # if white_BK:
    #     ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
    # else:
    #     ax.imshow(bkgimg, cmap='gray')
    # for vi in np.arange(numCell):
    #     if (visResp[vi]>0) & (~np.isnan(cpi[vi])) :
    #     #     ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
    #     # if ((hue_ax_auc[vi] >hueaucThres) | (hue_di_auc[vi]>hueaucThres)) & (visResp[vi]>0):
    #         ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=cmap_patch_cpi((cpi[vi]/2+0.5))))
    #
    # ax.set_rasterized(True)
    # plt.axis('off')
    # # plt.savefig('%s/%s_CPI_auc%s.svg' % (result_path_hue, file_name_hue, hueaucThres), dpi=300, format='svg')
    # plt.savefig('%s/%s_CPI_auc%s.png' % (result_path_plane, file_name_hue, hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    # plt.savefig('%s/%s_CPI_auc%s.png' % (result_path_hue, file_name_hue, hueaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')

    # plot orientation map from hue stimuli
    fig, ax = plt.subplots()
    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
    else:
        ax.imshow(bkgimg, cmap='gray')
    
    for vi in np.arange(numCell):
        if visResp[vi]>0:
           ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
        if (hueori_auc[vi]>oriaucThres) & (visResp[vi]>0):
            if np.isnan(fit_hueori[vi]):
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=gray))
            else:
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1 , color= ori_lut(fit_hueori[vi]/180)))  #alpha=np.log10(ori_p[vi])/np.min(np.log10(ori_p)) alpha=1-(ori_p[vi]-np.min(ori_p[:]))/(np.max(ori_p[:])-np.min(ori_p[:])),
           # ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=ori_lut(np.mod(fit_dir[vi]-90,180) / 180)))
    ax.set_rasterized(True)
    plt.axis('off')
    # plt.savefig('%s/%s_auc%s_fit.svg' % (result_path_ori,  file_name_hueori), facecolor='#ffffff', dpi=300, format='svg')
    plt.savefig('%s/%s_auc%s_fit.png' % (result_path_plane, file_name_hueori,oriaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_auc%s_fit.png' % (result_path_ori, file_name_hueori,oriaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
 
    # plot direction map from hue stimuli
    fig, ax = plt.subplots()
    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
    else:
        ax.imshow(bkgimg, cmap='gray')
    
    for vi in np.arange(numCell):
        if visResp[vi]>0:
           ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
        if (huedir_auc[vi]>diraucThres) & (visResp[vi]>0):
            if np.isnan(fit_huedir[vi]):
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=gray))
            else:
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=ori_lut(fit_huedir[vi]/360))) # direction
                # ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=ori_lut(np.mod(fit_huedir[vi]-90,180)/180))) # orientation inferred from direction
    ax.set_rasterized(True)
    plt.axis('off')
    # plt.savefig('%s/%s_auc%s_fit.svg' % (result_path_ori,  file_name_huedir), facecolor='#ffffff', dpi=300, format='svg')
    plt.savefig('%s/%s_auc%s_fit.png' % (result_path_plane, file_name_huedir,oriaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_auc%s_fit.png' % (result_path_ori, file_name_huedir,oriaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')


######

## Plot orientation, direction and sf maps from achromatic stimuli
if plot_ori:
    
# plot orientation map using fitting parameters
    fig, ax = plt.subplots()
    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
    else:
        ax.imshow(bkgimg, cmap='gray')
    
    for vi in np.arange(numCell):
        if visResp[vi]>0:
           ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
        if (ori_auc[vi]>oriaucThres) & (visResp[vi]>0):
            if np.isnan(fit_ori[vi]):
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=gray))
            else:
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1 , color=ori_lut(fit_ori[vi]/180)))  #alpha=np.log10(ori_p[vi])/np.min(np.log10(ori_p)) alpha=1-(ori_p[vi]-np.min(ori_p[:]))/(np.max(ori_p[:])-np.min(ori_p[:])),
    ax.set_rasterized(True)
    plt.axis('off')
    # plt.savefig('%s/%s_auc%s_fit.svg' % (result_path_ori,  file_name_ori), facecolor='#ffffff', dpi=300, format='svg')
    plt.savefig('%s/%s_auc%s_fit.png' % (result_path_plane, file_name_ori,oriaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_auc%s_fit.png' % (result_path_ori, file_name_ori,oriaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')


# plot orientation map using cv parameters
    fig, ax = plt.subplots()
    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
    else:
        ax.imshow(bkgimg, cmap='gray')
    
    for vi in np.arange(numCell):
        if visResp[vi]>0:
            ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
            # ax.text(xloc[cellId[vi]-1],yloc[cellId[vi]-1], str(cellId[vi]),horizontalalignment='center', fontsize=3)   # add cellId
        if (ori_auc[vi]>oriaucThres) & (visResp[vi]>0):
            ax.add_patch(PolygonPatch(allROIs[vi], alpha=1 , color=ori_lut(cv_ori[vi]/180)))  #alpha=np.log10(ori_p[vi])/np.min(np.log10(ori_p)) alpha=1-(ori_p[vi]-np.min(ori_p[:]))/(np.max(ori_p[:])-np.min(ori_p[:])),
            # ax.text(xloc[cellId[vi]-1],yloc[cellId[vi]-1], str(cellId[vi]),horizontalalignment='center', fontsize=3)   # add cellId
    ax.set_rasterized(True)
    plt.axis('off')
    # plt.savefig('%s/%s_auc%s_cv.svg' % (result_path_ori,  file_name_ori), facecolor='#ffffff', dpi=300, format='svg')
    plt.savefig('%s/%s_auc%s_cv_id.png' % (result_path_plane, file_name_ori,oriaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_auc%s_cv_id.png' % (result_path_ori, file_name_ori,oriaucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')


   # plot OSI map, osi is calculated after fitting
    fig, ax = plt.subplots()
    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
    else:
        ax.imshow(bkgimg, cmap='gray')
    for vi in np.arange(numCell):
        if visResp[vi]>0:
           ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
        if (ori_auc[vi] > oriaucThres) & (visResp[vi]>0):
           if np.isnan(fit_osi[vi]):
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=gray))
           else:
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=cmap_patch_osi(fit_osi[vi])))
    # ax.set_title( file_name_ori, fontsize=font_size)
    ax.set_rasterized(True)
    plt.axis('off')
    # plt.savefig('%s/%s_osi_fit.svg' % (result_path_ori,  file_name_ori), facecolor='#ffffff', dpi=300, format='svg')
    plt.savefig('%s/%s_osi_fit.png' % (result_path_plane, file_name_ori), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_osi_fit.png' % (result_path_ori, file_name_ori), dpi=300, bbox_inches='tight', pad_inches=0, format='png')


    # plot Circular Variance map
    fig, ax = plt.subplots()
    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
    else:
        ax.imshow(bkgimg, cmap='gray')
    for vi in np.arange(numCell):
        if visResp[vi] > 0:
            ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=gray))
        if (ori_auc[vi] > oriaucThres) & (visResp[vi] > 0):
            if np.isnan(ori_cv[vi]):
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=gray))
            else:
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=cmap_patch_osi(ori_cv[vi])))
    # ax.set_title( file_name_ori, fontsize=font_size)
    ax.set_rasterized(True)
    plt.axis('off')
    # plt.savefig('%s/%s_cv_fit.svg' % (result_path_ori,  file_name_ori), facecolor='#ffffff', dpi=300, format='svg')
    plt.savefig('%s/%s_cv_fit.png' % (result_path_plane, file_name_ori), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_cv_fit.png' % (result_path_ori, file_name_ori), dpi=300, bbox_inches='tight', pad_inches=0, format='png')


    # Plot spatial frequency preference map using fitting parameters
    fig, ax = plt.subplots()
    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
    else:
        ax.imshow(bkgimg, cmap='gray')
    
    for vi in np.arange(numCell):
        if visResp[vi]>0:
            if np.isnan(fit_sf[vi]):
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
            else:
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=cmap_patch_sf(fit_sf[vi] / np.max(fit_sf))))
    ax.set_rasterized(True)
    # ax.set_title('Spatial Frequency Tuning Map', fontsize=font_size)
    plt.axis('off')
    # plt.savefig('%s/%s_sf.svg' % (result_path_ori, file_name_sf), dpi=300, format='svg')
    plt.savefig('%s/%s_fit.png' % (result_path_plane, file_name_sf), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_fit.png' % (result_path_ori, file_name_sf), dpi=300, bbox_inches='tight', pad_inches=0, format='png')


   # plot Direction map using fitting parameters
    fig, ax = plt.subplots()
    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
        # file_name_dir = os.path.join('%s_whtBK%s' % (file_name_dir, str(white_BK)))
    else:
        ax.imshow(bkgimg, cmap='gray')
    
    for vi in np.arange(numCell):
        if visResp[vi]>0:
           ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
        if (dir_auc[vi]>diraucThres) & (visResp[vi]>0):
            if np.isnan(fit_dir[vi]):
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=gray))
            else:
                ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=ori_lut(fit_dir[vi]/360)))  # direction
                # ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=ori_lut(np.mod(fit_dir[vi]-90,180)/180)))  # orientation inferred from direction
    ax.set_rasterized(True)
    plt.axis('off')
    # plt.savefig('%s/%s_auc%s_fit.svg' % (result_path_ori,  file_name_dir), facecolor='#ffffff', dpi=300, format='svg')
    plt.savefig('%s/%s_auc%s_fit.png' % (result_path_plane, file_name_dir,diraucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_auc%s_fit.png' % (result_path_ori, file_name_dir,diraucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')


   # plot Direction map using cv parameters
    fig, ax = plt.subplots()
    if white_BK:
        ax.imshow(whitBK, cmap='gray', alpha=BK_alpha, vmin=0, vmax=1)
        # file_name_dir = os.path.join('%s_whtBK%s' % (file_name_dir, str(white_BK)))
    else:
        ax.imshow(bkgimg, cmap='gray')
    
    for vi in np.arange(numCell):
        if visResp[vi]>0:
           ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color= gray))
        if (dir_auc[vi]>diraucThres) & (visResp[vi]>0):
            ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=ori_lut(cv_dir[vi]/360)))  # direction
            # ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=ori_lut(np.mod(cv_dir[vi]-90,180)/180)))  # orientation inferred from direction
    ax.set_rasterized(True)
    plt.axis('off')
    # plt.savefig('%s/%s_auc%s_cv.svg' % (result_path_ori,  file_name_dir), facecolor='#ffffff', dpi=300, format='svg')
    plt.savefig('%s/%s_auc%s_cv.png' % (result_path_plane, file_name_dir,diraucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')
    plt.savefig('%s/%s_auc%s_cv.png' % (result_path_ori, file_name_dir,diraucThres), dpi=300, bbox_inches='tight', pad_inches=0, format='png')