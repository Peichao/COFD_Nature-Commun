
# Peichao's Notes:
# 1. Code was written for 2P data (Hartley) from Scanbox. Will export results (dataframe and csv) for plotting.
# 2. If you have multiple planes, it works with splited & interpolated data. Note results are slightly different.
# 3. If you have single plane, need to set interpolateData to false.

using NeuroAnalysis,Statistics,StatsBase,StatsPlots,DataFrames,DataFramesMeta,Mmap,Images,Interact
using CSV,MAT,DataStructures,HypothesisTests,StatsFuns,Random,Plots, ePPR

# Expt info
disk = ":"
subject = ""  # Animal
recordSession = "" # Unit
testId = ""  # Stimulus test

interpolatedData = true   # If you have multiplanes. True: use interpolated data; false: use uniterpolated data. Results are slightly different.
hartelyBlkId = 5641
stanorm = nothing
stawhiten = nothing

delays = -0.066:0.033:0.4
print(collect(delays))
isplot = false

## Prepare data & result path
exptId = join(filter(!isempty,[recordSession, testId]),"_")
dataFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), exptId)
metaFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), "metaFiles")

## load expt, scanning parameters
# metaFile=matchfile(Regex("[A-Za-z0-9]*_[A-Za-z0-9]*_[A-Za-z0-9]*$testId*_meta.mat"),dir=metaFolder,join=true)[1]
metaFile=matchfile(Regex("[A-Za-z0-9]*_[A-Za-z0-9]*_$testId*_ot_meta.mat"),dir=metaFolder,join=true)[1]
dataset = prepare(metaFile)

ex = dataset["ex"]
envparam = ex["EnvParam"]
coneType = getparam(envparam,"colorspace")
sbx = dataset["sbx"]["info"]
sbxft = ex["frameTimeSer"]   # time series of sbx frame in whole recording
# Condition Tests
envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
condon = ex["CondTest"]["CondOn"]
condoff = ex["CondTest"]["CondOff"]
condidx = ex["CondTest"]["CondIndex"]
# condtable = DataFrame(ex["Cond"])
condtable =  DataFrame(ex["raw"]["log"]["randlog_T1"]["domains"]["Cond"])
rename!(condtable, [:oridom, :kx, :ky,:bwdom,:colordom])

# find out blanks and unique conditions
blkidx = condidx .>= hartelyBlkId  # blanks start from 5641
cidx = .!blkidx
condidx2 = condidx.*cidx + blkidx.* hartelyBlkId

## Load data
if interpolatedData
    segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.segment"),dir=dataFolder,join=true)[1]
    signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.signals"),dir=dataFolder,join=true)[1]
else
    segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*.segment"),dir=dataFolder,join=true)[1]
    signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9].signals"),dir=dataFolder,join=true)[1]
end
segment = prepare(segmentFile)
signal = prepare(signalFile)
# sig = transpose(signal["sig"])   # 1st dimention is cell roi, 2nd is fluorescence trace
spks = transpose(signal["spks"])  # 1st dimention is cell roi, 2nd is spike train

##
# Prepare Imageset
downsample = haskey(param,:downsample) ? param[:downsample] : 2
sigma = haskey(param,:sigma) ? param[:sigma] : 1.5
# bgRGB = [getparam(envparam,"backgroundR"),getparam(envparam,"backgroundG"),getparam(envparam,"backgroundB")]
bgcolor = RGBA([0.5,0.5,0.5,1]...)
coneType = string(getparam(envparam,"colorspace"))
masktype = getparam(envparam,"mask_type")
maskradius = getparam(envparam,"mask_radius")
masksigma = 1#getparam(envparam,"Sigma")
hartleyscale = haskey(param,:hartleyscale) ? param[:hartleyscale] : 1
hartleynorm = haskey(param, :hartleynorm) ? param[:hartleynorm] : false
xsize = getparam(envparam,"x_size")
ysize = getparam(envparam,"y_size")
stisize = xsize
ppd = haskey(param,:ppd) ? param[:ppd] : 52
ppd = ppd/downsample
imagesetname = "Hartley_stisize$(stisize)_hartleyscalescale$(hartleyscale)_ppd$(ppd)"
maskradius = maskradius /stisize + 0.03


if !haskey(param,imagesetname)
    imageset = Dict{Any,Any}(:image =>map(i->GrayA.(hartley(kx=i.kx,ky=i.ky,bw=i.bwdom,stisize=stisize, ppd=ppd,norm=hartleynorm,scale=hartleyscale)),eachrow(condtable)))
    imageset[:sizepx] = size(imageset[:image][1])
    param[imagesetname] = imageset
end

# Prepare Image Stimuli
imageset = param[imagesetname]
bgcolor = oftype(imageset[:image][1][1],bgcolor)
imagestimuliname = "bgcolor$(bgcolor)_masktype$(masktype)_maskradius$(maskradius)_masksigma$(masksigma)" # bgcolor and mask define a unique masking on an image set
if !haskey(imageset,imagestimuliname)
    imagestimuli = Dict{Any,Any}(:stimuli => map(i->alphablend.(alphamask(i,radius=maskradius,sigma=masksigma,masktype=masktype)[1],[bgcolor]),imageset[:image]))
    imagestimuli[:unmaskindex] = alphamask(imageset[:image][1],radius=maskradius,sigma=masksigma,masktype=masktype)[2]
    imageset[imagestimuliname] = imagestimuli
end
imagestimuli = imageset[imagestimuliname]

## Load data
planeNum = size(segment["mask"],3)  # how many planes
if interpolatedData
    planeStart = vcat(1, length.(segment["seg_ot"]["vert"]).+1)
end
## Use for loop process each plane seperately
for pn in 1:planeNum
    # pn=1  # for test
    display("plane: $pn")
    # Initialize DataFrame for saving results
    recordPlane = string("00",pn-1)  # plane/depth, this notation only works for expt has less than 10 planes
    siteId = join(filter(!isempty,[recordSession, testId, recordPlane]),"_")
    dataExportFolder = joinpath(disk,subject, "2P_analysis", join(["U",recordSession]), siteId, "DataExport")
    resultFolder = joinpath(disk,subject, "2P_analysis", join(["U",recordSession]), siteId, "Plots")
    isdir(dataExportFolder) || mkpath(dataExportFolder)
    isdir(resultFolder) || mkpath(resultFolder)

    if interpolatedData
        cellRoi = segment["seg_ot"]["vert"][pn]
    else
        cellRoi = segment["vert"]
    end
    cellNum = length(cellRoi)
    display("Cell Number: $cellNum")

    if interpolatedData
        # rawF = sig[planeStart[pn]:planeStart[pn]+cellNum-1,:]
        spike = spks[planeStart[pn]:planeStart[pn]+cellNum-1,:]
    else
        # rawF = sig
        spike = spks
    end

    imagesize = imageset[:sizepx]
    xi = imagestimuli[:unmaskindex]

    # estimate RF using STA
    if :STA in param[:model]
        uci = unique(condidx2)
        ucii = map(i->findall(condidx2.==i),deleteat!(uci,findall(isequal(hartelyBlkId),uci)))   # find the repeats of each unique condition
        ubii = map(i->findall(condidx2.==i), [hartelyBlkId])

        uy = Array{Float64}(undef,cellNum,length(delays),length(ucii))
        ucy = Array{Float64}(undef,cellNum,length(delays),length(ucii))
        uby = Array{Float64}(undef,cellNum,length(delays),length(ubii))
        uŷest = Array{Float64}(undef,cellNum,length(delays),length(ucii))
        ugof = Array{Any}(undef,cellNum,length(delays))

        usta = Array{Float64}(undef,cellNum,length(delays),length(xi))
        cx = Array{Float64}(undef,length(ucii),length(xi))
        foreach(i->cx[i,:]=gray.(imagestimuli[:stimuli][uci[i]][xi]),1:size(cx,1))
        for d in eachindex(delays)
            display("Processing delay: $d")
            y,num,wind,idx = epochspiketrain(sbxft,condon.+delays[d], condoff.+delays[d],isminzero=false,ismaxzero=false,shift=0,israte=false)
            spk=zeros(size(spike,1),length(idx))
            for i =1:length(idx)
                spkepo = @view spike[:,idx[i][1]:idx[i][end]]
                spk[:,i]= mean(spkepo, dims=2)
            end
            for cell in 1:cellNum
                # display(cell)
                cy = map(i->mean(spk[cell,:][i]),ucii)  # response to grating
                bly = map(i->mean(spk[cell,:][i]),ubii) # response to blank, baseline
                ry = cy.-bly  # remove baseline
                csta = sta(cx,ry,norm=stanorm,whiten=stawhiten)  # calculate sta
                ucy[cell,d,:]=cy
                uby[cell,d,:]=bly
                uy[cell,d,:]=ry
                usta[cell,d,:]=csta
                ŷest=cx*csta
                uŷest[cell,d,:]=ŷest
                ugof[cell,d] = goodnessoffit(ry,ŷest,k=length(xi))

                if isplot
                    r = [extrema(csta)...]
                    title = "Unit_$(cell)_STA_$(delays[d])"
                    p = plotsta(csta,sizepx=imagesize,sizedeg=stisize,ppd=ppd,index=xi,title=title,r=r)
                    foreach(i->save(joinpath(resultFolder,"$title$i"),p),[".png"])
                end
            end
        end
        save(joinpath(dataExportFolder,join([subject,"_",siteId,"_",coneType,"_sta.jld2"])),"imagesize",imagesize,"cx",cx,"xi",xi,"xcond",condtable[uci,:],"uy",uy,"ucy",ucy,"uŷest",uŷest,"usta",usta,"uby",uby,"ugof",ugof,"delays",delays,"maskradius",maskradius,"stisize",stisize,"color",coneType)
    end
end
