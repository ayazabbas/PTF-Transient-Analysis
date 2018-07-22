##install astroquery with cmd command: conda install -c astropy astroquery
import time
import astropy.units as u
import numpy as np
import wget
import os
import shutil
import math
from astroquery.irsa import Irsa
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
from astropy import wcs

fig = plt.figure()
matplotlib.rcParams.update({'font.size': 4})
Irsa.ROW_LIMIT = -1
Irsa.TIMEOUT = 120

def getVizierCatalogue(query):
    catalog_list = Vizier.find_catalogs(query)
    v = Vizier(columns=["RAJ2000","DEJ2000","W(Ha)","pmRA","pmDE"])
    v.ROW_LIMIT = -1
    catalogs = v.get_catalogs(catalog_list.keys())
    data = catalogs[0]
    return data

def properMotions(data):
    positions = np.column_stack([(data[:,]["RAJ2000"] + data[:,]["pmRA"]/360), data[:,]["DEJ2000"] + data[:,]["pmDE"]/360])
    return positions

def downloadLightCurve(ra, dec): 
    coords = SkyCoord(ra, dec, unit='deg')
    table = Irsa.query_region(coordinates=coords, catalog='ptf_lightcurves', radius=5*u.arcsec)
    table = table.filled(-99) ##fill masked vaules as -99, e.g. when limitmag = null, it will equal -99   
    return table

def getNearbyStars(ra, dec):    
    result = queryDatabase(ra, dec)
    ctlgName = downloadCtlg(result[result.shape[0] - 1]['afilename3']) ##downloads last catalog in query result
    catHDU = fits.open(ctlgName)
    catData = catHDU[1].data
    catHDU.close()
    
    stars = np.column_stack([catData['X_WORLD'], catData['Y_WORLD'], catData['MAG_APER'][:,0]]) #retrieve necessary data into an array
    stars = stars[stars[:,2].argsort()] #sort array by magnitude
    return stars, ctlgName  
  
def queryDatabase(ra, dec):
    ##Removes query.txt from local folder so no query(1).txt appears
    if "query.txt" in os.listdir(os.getcwd()):
        os.remove("query.txt")
    ##Uses the RA and Dec with level1 search on ptf
    url = "http://irsa.ipac.caltech.edu/ibe/search/ptf/images/level1?POS="+str(ra)+","+str(dec)
    ##Downloads the query results as a txt file
    filename = wget.download(url, "query.txt")
    ##Downloaded txt file is IPAC encoded so ascii.io.read with ipac tag
    t = Table.read(filename, format = "ipac")
    return t.as_array()

def downloadCtlg(URLExtension):
    ##Takes URL extension from query file for sextractor data
    url = "http://irsa.ipac.caltech.edu/ibe/data/ptf/images/level1/"+URLExtension
    ##Split the url extension to get proper filename from url
    filenames = URLExtension.split("/")
    ## Download the ctlg file
    filename = wget.download(url, filenames[len(filenames)-1])
    return filename   
 
def downloadPTF(URLExtension, filename):
    ##Takes URL extension from query file for sextractor data
    url = "http://irsa.ipac.caltech.edu/ibe/data/ptf/images/level1/"+URLExtension
    ## Download the file
    return wget.download(url, filename)
    
def removeCtlg(filename):
    ##Removes a file in the current folder with the given filename
    if filename in os.listdir(os.getcwd()):
        os.remove(filename)

def plotLightCurves(data, refData, filename):
    rtime = np.empty([0])
    rmag = np.empty([0])
    rfwhm = np.empty([0])
    rlimitmag = np.empty([0])
    gtime = np.empty([0])
    gmag = np.empty([0])
    gfwhm = np.empty([0])
    glimitmag = np.empty([0])
    
    imageDates = []
    increased = False
    gcount = 0
    for x in range(0, data.shape[0]):
        if int(data[x]['fid']) == 2:
            rtime = np.insert(rtime, 0, data[x]['obsmjd'])
            rmag = np.insert(rmag, 0, data[x]['mag_autocorr'])
            rfwhm = np.insert(rfwhm, 0, data[x]['fwhmsex'])
            rlimitmag = np.insert(rlimitmag, 0, data[x]['limitmag'])
        else:
            gcount += 1
            gtime = np.insert(gtime, 0, data[x]['obsmjd'])
            gmag = np.insert(gmag, 0, data[x]['mag_autocorr'])
            gfwhm = np.insert(gfwhm, 0, data[x]['fwhmsex'])
            glimitmag = np.insert(glimitmag, 0, data[x]['limitmag'])
            
    if gcount < 10:
        rmedian = np.median(rmag)
        targaxr = fig.add_subplot(4, 1, 1)
        refaxr = fig.add_subplot(4, 1, 2, sharex=targaxr)
        fwhmaxr = fig.add_subplot(4, 1, 3, sharex=targaxr)
        limitmagaxr = fig.add_subplot(4, 1, 4, sharex=targaxr)
        ##clean data further by removing magnitude increases of 5 or more
        badIndices = []
        for x in range(0, rmag.shape[0]):            
            if rmag[x] - rmedian > 5:
                badIndices.extend([x])
        rmag = np.delete(rmag, badIndices, axis=0) #remove bad rows
        rtime = np.delete(rtime, badIndices, axis=0)
        rfwhm = np.delete(rfwhm, badIndices, axis=0)
        rlimitmag = np.delete(rlimitmag, badIndices, axis=0)
        rmedian = np.median(rmag)
        ##check for outliers
        rzscores = np.column_stack([rtime, rmag, (0.6745*(rmag - rmedian)) / MAD(rmag)])
        for rz in rzscores:
            date = [rz[0], rz[1], 'R']
            if abs(rz[2]) > 3.5 and abs(rz[1] - rmedian) > 0.5 and date not in imageDates: ##test modified z score and difference from median to flag for image download
                imageDates.extend([date])    
                if (rz[1] - rmedian) < -0.5:
                    increased = True
        ##Plot the data
        targaxr.plot(rtime, rmag - rmedian, 'r.', markersize=1)
        fwhmaxr.plot(rtime, rfwhm,  'k.', markersize=1)
        limitmagaxr.plot(rtime, rlimitmag,  'k.', markersize=1)
    else: ##if plotting both R and G
        rmedian = np.median(rmag)
        gmedian = np.median(gmag)
        targaxr = fig.add_subplot(4, 2, 1)
        refaxr = fig.add_subplot(4, 2, 3, sharex=targaxr)
        fwhmaxr = fig.add_subplot(4, 2, 5, sharex=targaxr)
        limitmagaxr = fig.add_subplot(4, 2, 7, sharex=targaxr)
        targaxg = fig.add_subplot(4, 2, 2)
        refaxg = fig.add_subplot(4, 2, 4, sharex=targaxg)
        fwhmaxg = fig.add_subplot(4, 2, 6, sharex=targaxg)
        limitmagaxg = fig.add_subplot(4, 2, 8, sharex=targaxg)
        ##clean data further by removing magnitude increases of 5 or more
        badIndicesR = []
        for x in range(0, rmag.shape[0]):
            if rmag[x] - rmedian > 5:
                badIndicesR.extend([x])
        rmag = np.delete(rmag, badIndicesR, axis=0) #remove bad rows
        rtime = np.delete(rtime, badIndicesR, axis=0)
        rfwhm = np.delete(rfwhm, badIndicesR, axis=0)
        rlimitmag = np.delete(rlimitmag, badIndicesR, axis=0)
        rmedian = np.median(rmag)
        ##clean g
        badIndicesG = []
        for x in range(0, gmag.shape[0]):
            if (gmag[x] - gmedian) > 5:
                badIndicesG.extend([x])
        gmag = np.delete(gmag, badIndicesG, axis=0) #remove bad rows
        gtime = np.delete(gtime, badIndicesG, axis=0)
        gfwhm = np.delete(gfwhm, badIndicesG, axis=0)
        glimitmag = np.delete(glimitmag, badIndicesG, axis=0)
        gmedian = np.median(gmag)
        ##check for outliers in R
        rzscores = np.column_stack([rtime, rmag, (0.6745*(rmag - rmedian)) / MAD(rmag)])
        for rz in rzscores:
            date = [rz[0], rz[1], 'R']
            if abs(rz[2]) > 3.5 and abs(rz[1] - rmedian) > 0.5 and date not in imageDates: ##test modified z score and difference from median to flag for image download
                imageDates.extend([date])
                if (rz[1] - rmedian) < -0.5:
                    increased = True
        
        ##check for outliers in G
        gzscores = np.column_stack([gtime, gmag, (0.6745*(gmag - gmedian)) / MAD(gmag)])
        for gz in gzscores:
            date = [gz[0], gz[1], 'G']
            if abs(gz[2]) > 3.5 and abs(gz[1] - gmedian) > 0.5 and date not in imageDates: ##test modified z score and difference from median to flag for image download
                imageDates.extend([date])
                if (gz[1] - gmedian) < -0.5:
                    increased = True
    ##Plot the data
        targaxr.plot(rtime, rmag - rmedian, 'r.', markersize=1)
        fwhmaxr.plot(rtime, rfwhm,  'k.', markersize=1)
        limitmagaxr.plot(rtime, rlimitmag,  'k.', markersize=1)
        targaxg.plot(gtime, gmag - gmedian, 'r.', markersize=1)
        fwhmaxg.plot(gtime, gfwhm,  'k.', markersize=1)
        limitmagaxg.plot(gtime, glimitmag,  'k.', markersize=1)
        
    ##add x to plot where images are to be downloaded
    if len(imageDates) > 0:
        for i in imageDates: 
            if i[2] == 'R':
                targaxr.plot(i[0], i[1] - rmedian, 'kx', markersize=2)
            if i[2] == 'G':
                targaxg.plot(i[0], i[1] - gmedian, 'kx', markersize=2)
    
    colour = ''
    for x in range(0, len(refData)):
        rRefMag = np.empty([0])
        rRefTime = np.empty([0])
        gRefMag = np.empty([0])
        gRefTime = np.empty([0])
        if x == 0:
            colour = 'b'
        if x == 1:
            colour = 'g'
        if x == 2:
            colour = 'y'
        for i in range(0, refData[x].shape[0]):
            if int(refData[x][i]['fid']) == 2:
                rRefMag = np.insert(rRefMag, 0, refData[x][i]['mag_autocorr'])
                rRefTime = np.insert(rRefTime, 0, refData[x][i]['obsmjd'])
            else:
                if gcount > 10:
                    gRefMag = np.insert(gRefMag, 0, refData[x][i]['mag_autocorr'])
                    gRefTime = np.insert(gRefTime, 0, refData[x][i]['obsmjd'])
        refaxr.plot(rRefTime, rRefMag - np.median(rRefMag), colour + '.', markersize=1)
        if gRefMag.shape[0] > 0:
            refaxg.plot(gRefTime, gRefMag - np.median(gRefMag), colour + '.', markersize=1)
    for ax in fig.get_axes():
        ax.invert_yaxis()
    
    if increased == True:
        plt.savefig('lightcurves\\increased\\' + filename + '.png', dpi=2000)
    else:
        plt.savefig('lightcurves\\' + filename + '.png', dpi=2000)
    plt.cla()
    plt.clf()
    return imageDates, increased
    
def cleanData(data):
    badIndices = []
    for x in range(0, data.shape[0]):
        if data[x]['limitmag'] == -99 or data[x]['fwhmsex'] == -99: #or float(data[x]['obsmjd']) < 55200: ##points before ~this date are often dodgy
            badIndices.extend([x])
    data = np.delete(data, badIndices, axis=0) ##remove the bad rows
    return data
    
def extractData(filename):
    t = Table.read(filename, format='ipac')
    data = t.as_array() ##convert astropy Table to a numpy array
    return data
    
def MAD(data): ##median absolute deviations (supply 1d numpy array)
    return np.median(abs(data - np.median(data)))
    
def pixToWCS(header, x, y):
    w = wcs.WCS(header)
    pixcrd = np.array([[x, y]], np.float_)
    # Convert pixel coordinates to world coordinates
    # The second argument is "origin" -- in this case we're declaring we
    # have 1-based coordinates.
    world = w.wcs_pix2world(pixcrd, 1)
    return world
    
def WCSToPix(header, ra, dec):
    w = wcs.WCS(header)
    world = np.array([[float(ra), float(dec)]])
    pixcrd = w.wcs_world2pix(world, 1)
    return pixcrd
    
def plotImages(path, targRa, targDec, refRa, refDec, increased):
    quality = 2500
    images = os.listdir(path)
    imgTotal = len(images)
    imgNumber = 1
    offset = 5
    if imgTotal <= 4:
        rows = 2
        cols = 2
    elif imgTotal <= 9:
        rows = 3
        cols = 3
    elif imgTotal <= 16:
        quality = 3000
        rows = 4
        cols = 4
    elif imgTotal <= 25:
        offset = 2
        matplotlib.rcParams.update({'font.size': 2})
        quality = 3500
        rows = 5
        cols = 5
    elif imgTotal <= 36:
        offset = 1
        matplotlib.rcParams.update({'font.size': 2})
        quality = 4000
        rows = 6
        cols = 6
    for i in images:
        fitsHDU = fits.open(path + i)
        header = fitsHDU[0].header
        imgData = fitsHDU[0].data
        fitsHDU.close()
        ax = fig.add_subplot(rows, cols, imgNumber)
        circle1 = plt.Circle((WCSToPix(header, targRa, targDec)[0]), color='r', fill=False, radius=8, linewidth=0.2)
        ax.add_artist(circle1)
        for x in range(0, len(refRa)):
            if x == 0:
                colour = 'b'
            if x == 1:
                colour = 'g'
            if x == 2:
                colour = 'y'
            circle = plt.Circle((WCSToPix(header, refRa[x], refDec[x])[0]), color=colour, fill=False, radius=8, linewidth=0.2)
            ax.add_artist(circle)
        
        ##get brightest pixel in 20x20 region around target for contrast
        maxValue = 0
        for n in range(0, len(refRa)):
            targetPixel = WCSToPix(header, refRa[n], refDec[n])[0]
            if targetPixel[0] > 25 and targetPixel[1] > 35 and targetPixel[0] < (imgData.shape[0] - 25) and targetPixel[1] < (imgData.shape[1] - 25):
                for x in range(int(targetPixel[0]-10), int(targetPixel[0]+11)):
                    for y in range(int(targetPixel[1]-10), int(targetPixel[1]+11)):
                        pixelCount = imgData[x][y]
                        if pixelCount > maxValue:
                            maxValue = pixelCount
        if maxValue == 0:
            targetPixel = WCSToPix(header, targRa, targDec)[0]
            if targetPixel[0] > 25 and targetPixel[1] > 35 and targetPixel[0] < (imgData.shape[0] - 25) and targetPixel[1] < (imgData.shape[1] - 25):
                for x in range(int(targetPixel[0]-10), int(targetPixel[0]+11)):
                    for y in range(int(targetPixel[1]-10), int(targetPixel[1]+11)):
                        pixelCount = imgData[x][y]
                        if pixelCount > maxValue:
                            maxValue = pixelCount
        if maxValue == 0:
            maxValue = 5000
        ##get average count to use as minimum threshold
        totalCount = 0
        pixels = 0
        lowestPixel = imgData[5][5]
        for x in range(0, imgData.shape[0]-1):
            for y in range(0, imgData.shape[1]-1):
                pixels += 1
                totalCount = totalCount + imgData[x][y]
                if imgData[x][y] < lowestPixel:
                    lowestPixel = imgData[x][y]
        minValue = totalCount / pixels
        if minValue > maxValue:
            minValue = maxValue - 1000
        ax.imshow(imgData, cmap='gray_r', origin='lower', vmin=minValue, vmax=maxValue)
        limitmag = 'null'
        for k in header.keys():
            if k == 'LIMITMAG':
                limitmag = str(header[k])
                break
            if k == 'LMGAPCZP':
                limitmag = str(header[k])
                break
        if i.endswith('o.fits'):
            ax.set_xlabel(str(header['OBSMJD'])+', '+str(header['SEEING'])+', '+limitmag+', '+str(header['FILTER'])+', o', labelpad=offset)
        else:
            ax.set_xlabel(str(header['OBSMJD'])+', '+str(header['SEEING'])+', '+limitmag+', '+str(header['FILTER']), labelpad=offset)
        imgNumber += 1
    for ax in fig.get_axes():
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
    if increased == True:
        plt.savefig('lightcurves\\increased\\' + str(targRa) + ',' + str(targDec) + 'imgs.png', dpi=quality)
    else:
        plt.savefig('lightcurves\\' + str(targRa) + ',' + str(targDec) + 'imgs.png', dpi=quality)
    plt.cla()
    plt.clf()
    matplotlib.rcParams.update({'font.size': 4})
    
def downloadImages(targRa, targDec, dates): ##dates is list of tuples with obsmjd and filter ('R' or 'G')
    result = queryDatabase(targRa, targDec)
    imagesToGet = []
    outlierDates = []
    path = 'images\\'+str(targRa)+','+str(targDec)+'\\'
    for d in dates:
        timeDifferences = np.column_stack([result[:,]['pfilename'], result[:,]['obsmjd'] - d[0], result[:,]['filter']])
        timeDifferences = timeDifferences[timeDifferences[:,1].argsort()]
        if abs(float(timeDifferences[0][1])) < 1 and timeDifferences[0][2] == d[2]:
            imagesToGet.extend([[timeDifferences[0][0], float(timeDifferences[0][1]) + d[0]]]) ##extend the list with the filename and obsmjd
            outlierDates.extend([float(timeDifferences[0][1]) + d[0]])
    if len(imagesToGet) > 0:
        ##add images before and after the outlier:
        previousDate = 0
        addedDate = False
        result = result[result[:,]['obsmjd'].argsort()]
        for r in result:
            date = [r['pfilename'], float(r['obsmjd'])]
            if addedDate == True and date not in imagesToGet:
                imagesToGet.extend([date])
                addedDate = False
            if date in imagesToGet and previousDate not in imagesToGet and previousDate != 0:
                imagesToGet.extend([previousDate])
                addedDate = True
            previousDate = [r['pfilename'], float(r['obsmjd'])]
        if not os.path.exists(path):
            os.makedirs(path)
        for i in imagesToGet:
            if i[1] in outlierDates:
                downloadPTF(str(i[0])+'?center='+str(targRa)+','+str(targDec)+'&size=600arcsec&gzip=false', path+str(i[1])+'o.fits')
            else:
                downloadPTF(str(i[0])+'?center='+str(targRa)+','+str(targDec)+'&size=600arcsec&gzip=false', path+str(i[1])+'.fits')
    return path
    
def deleteImages(path):
    if os.path.exists(path):
        try:
            shutil.rmtree(path)
        except:
            print 'Could not delete images'

def imageProcessing(attempts, positions):
    path = 'data\\'
    lightcurves = os.listdir('lightcurves\\')
    try:
        for i in range(attempts, len(positions)):
            print attempts
            downloaded = False
            targra = positions[i][0]
            targdec = positions[i][1]
            for l in lightcurves: ##check to see if data is already there in which case don't waste time redownloading
                if l.startswith(str(targra)+','+str(targdec)):
                    downloaded = True
            if downloaded == False:
                targetTable = downloadLightCurve(targra, targdec)
                targetCurve = targetTable.as_array()
                if targetCurve.shape[0] > 10: ##continue if the table has more than 10 rows
                    refsdone = 0            
                    refStars, ctlgName = getNearbyStars(targra, targdec) 
                    removeCtlg(ctlgName)
                    index = 0
                    while refsdone < 3:
                        if index == refStars.shape[0]:
                            break
                        refra = refStars[index][0]
                        refdec = refStars[index][1]
                        distance = math.sqrt((refra - targra)**2 + (refdec - targdec)**2)
                        if distance > 0.002 and distance < 0.07:
                            refTable = downloadLightCurve(refra, refdec)
                            refCurve = refTable.as_array()
                            if refCurve.shape[0] > 10 and refCurve[0]['medianmag'] > 13:
                                refsdone += 1
                                refTable.write('data\\references\\' + str(targra) + ',' + str(targdec) + ',' + 'ref,' + str(refra) + ',' + str(refdec) + '.txt', format='ipac')
                        index += 1
                    print 'Saving '+str(targra)+','+str(targdec)
                    filename = str(targra) + ',' + str(targdec) + '.txt'
                    targetTable.write(path + filename, format='ipac') ##save the table in ipac format as a .txt file
                    
                    ## lightcurve and image plotting
                    refData = []
                    targetData = cleanData(targetTable.as_array())
                    references = os.listdir(path + 'references\\')
                    refRa = []
                    refDec = []
                    for r in references:
                        if r.startswith(filename[:-4]):
                            refLoc = r.split(',')
                            refRa.extend([refLoc[3]])
                            refDec.extend([refLoc[4][:-4]])
                            refData.extend([cleanData(extractData(path + 'references\\' + r))])
                    if len(refData) >= 1:
                        imageDates, increased = plotLightCurves(targetData, refData, filename[:-4])
                        if len(imageDates) > 36: ##too many images are a pain to display
                            print 'Too many images for star: '+str(targra)+','+str(targdec)
                        elif len(imageDates) > 0:
                            imagePath = downloadImages(targra, targdec, imageDates)
                            plotImages(imagePath, targra, targdec, refRa, refDec, increased)
                            deleteImages(imagePath)
                    else:
                        print 'No comparison stars for '+str(targra)+','+str(targdec)
            attempts += 1
    except:
        print 'Encountered error: waiting 60s before continuing'
        time.sleep(60) ##Wait 60 seconds before retrying downloads
        imageProcessing(attempts+1, positions) ##Proceed with the next file to be downloaded
        
def main():
    targets = getVizierCatalogue('J/AJ/145/102') ##'J/AJ/145/102' 1500 stars 'J/AJ/141/97' 70k stars
    targets.sort('W_Ha_') ##Sort in terms of Hydrogen alpha line width
    positions = properMotions(targets.as_array()) ##Apply proper motions
    imageProcessing(0, positions) ##Select where in the catalogue to start processing, 28525

main()