"""
Gru: Jue Wang
Minion 1: Sajeeth Wimalasuriyan
Minion 2: Yifan Jiang
"""
# -*- coding: utf-8 -*-
import arcpy, arcinfo
import os
import time
import logging
import math
import pytz
import numpy as np
from arcpy import sa
import sys
from arcpy import env
from arcpy import da
from datetime import datetime, timedelta
from collections import namedtuple
import random

"""
Code below extracts variables from GUI.
"""
spatial_ref = arcpy.GetParameter(0) #Grabs raw GPS data.
yourOwnTemplate = arcpy.GetParameter(1) #Grabs optional user inputted scrambled data.
kernalBand = int(arcpy.GetParameter(3)) #Grabs kernal bandwidth parameter.
hotspotLocationsName = arcpy.GetParameter(2) #Grabs third parameter.
timing = str(arcpy.GetParameter(4)) #Grabs buffer variable to allow user to make either time or location more accurate.
stockBuffer = str(arcpy.GetParameter(6)) #Grabs buffer variable for DALK to determine how large of a 
										 #area around a point is considered close enough to conceal said point.
userProvidedHomeTime = str(arcpy.GetParameter(7)) #Gets the time at which the OP considers is suitable for the user to be at home.
userProvidedHomeDuration = float(arcpy.GetParameter(8)) #Gets the duration at which the OP considers is suitable for the user to be at home.
identifiedHome = []#Global variable used to store unique location ID of the home location.
class computationalCore(object):
	"""
	This class is responsible for the the pre processing of input GPS data and for the 
	identification of potential activity locations.
	"""
	def __init__(self, kernalBand):
		#Constant numerical constraints
		self.earthLength =  6378137.0 * 2 * math.pi #earth radius from semimajor axis of WGS 84 spheroid
		self.cellSize = 10. #raster cellsize of kernel in meters | should be small enough in regard to bandwidth
		self.maxRstSize = 30000 #define maximum column / line of raster
		self.resampleFrequency = 0 #0 #specify resampling frequency of GPS points (Hz) | 0 = no reseamping | resampling helps discarding smaller sinks
		self.minDuration2KeepHS = 150 #150 #Minimum duration of longest visit to Hotspot to keep it (in seconds)
		self.minVisitDuration = 300 #300 #Minimum duration of visit to keep it (in seconds)
		self.minVisitTimespan = 300 #300 #Minimum timespan between two visits to keep them distinct (in seconds)
		self.minLoopDuration = 300 #300 #Minimum duration of loop to keep it (in seconds)
		self.simplifyTolerance = 100 #100 #in meters
		self.smoothTolerance = 50 #50 #in meters
		# Variables below is what is used for calculations and are recycled
		self.kernelBandwidth = kernalBand
		self.srWgs84 = arcpy.SpatialReference(4326)
		#self.localtz = pytz.timezone(localtz)
		self.maxTimeVisit = {} #dictionary to store max time for each location (loc = key)
		self.visitSequence = {} #dictionary to store a list of each visit at each location (loc = key and list is composed of tuples (start, timespan))
		self.filtSequence = {}  #new dict to store filtered visits
		self.coordXYHotspot = {} #store xy coordinates of locations (later used to change GPS point coords)
		#Working variables
		self.wkp = ''
		self.kernelRST  = ''
		self.aoiKernelRST = ''
		self.invertedKernelRST = ''
		self.filteredKernelRST = ''
		self.flowDirRST = ''
		self.sinkRST = ''
		self.aoiBasinRST = ''
		self.zoneRST = ''
		self.timingDic = {}
		#Below are formatting and legacy variables not activly needed for program function.
		self.datefmt = '%Y-%m-%d %H:%M:%S'
		self.fld_date = 'utc_date'
		self.fld_locid = arcpy.ValidateFieldName('LocationID', self.wkp)
		self.fld_utcStart = arcpy.ValidateFieldName("UTCStartTime",self.wkp)
		self.fld_utcStop = arcpy.ValidateFieldName("UTCStopTime",self.wkp)
		self.fld_duration = arcpy.ValidateFieldName("Duration",self.wkp)
		self.fld_start = arcpy.ValidateFieldName("StartTime",self.wkp)
		self.fld_stop = arcpy.ValidateFieldName("StopTime",self.wkp)
		self.fld_fromLoc = arcpy.ValidateFieldName("FromLocID",self.wkp)
		self.fld_toLoc = arcpy.ValidateFieldName("ToLocID",self.wkp)
		self.fld_pathVertices = arcpy.ValidateFieldName("NbVertices",self.wkp)
		self.fld_interpolFlag = arcpy.ValidateFieldName("Interpolated",self.wkp)
		self.fld_ticks = arcpy.ValidateFieldName('NbTicks', self.wkp)
		self.fld_addr = arcpy.ValidateFieldName("Google_Address",self.wkp)
		self.fld_kernel = arcpy.ValidateFieldName("kernel",self.wkp)
		self.fld_mKernel = arcpy.ValidateFieldName('modifiedKernel',self.wkp)
		self.fld_snapFlag = arcpy.ValidateFieldName('Snap',self.wkp)
		

	def identifyTime(self, gpsFC,dissolvePolyFC):
		whatTimeIsHome =  datetime(2000,6,5).date()
		'''
			After initial hotspot are identification the sample is than resampled to identify critical
			points around activity areas. After these areas are identified the time stamps are gathered,
			parsed and the time spent at the location is identified.
		'''
		x = 'C:\DALKAY.gdb\LocationsBuffered'#Sets the location of the buffer

		arcpy.Buffer_analysis(dissolvePolyFC, x, timing + " meters","FULL","ROUND")#Creates buffer increasing the buffer increase the accuracy of the time measurement at the cost of locational percision.

		pointInPolygon = arcpy.SelectLayerByLocation_management(gpsFC, 'INTERSECT', 
                                                          x,) #Finds all the points within the buffer x
		polygonCoveringPoints = arcpy.SelectLayerByLocation_management(x, 'INTERSECT', 
                                                          gpsFC) #Created a polygon which overlaps the buffer x
		y = r'C:\DALKAY.gdb\p85'#Sets the location of the disolved version of polygonCoveringPoints
		dissolveFields = []
		arcpy.Dissolve_management(polygonCoveringPoints,y,dissolveFields, "", 
                          "SINGLE_PART")#Dissolves polygonCoveringPoints
		arcpy.AddField_management(y, "POINTZONE", "LONG")#Adds new field(variable) to attribute table of polygons in y.
		fields = ['POINTZONE']#Determines what field the cursor will look at
		"""
		Code below assigns a individual identification to each polygon
		"""
		uniqueIdentifier = 0
		with arcpy.da.UpdateCursor(y, fields) as cursor:
			for row in cursor:
				row[0] = uniqueIdentifier
				cursor.updateRow(row) 
				uniqueIdentifier += 1
		z = r'C:\DALKAY.gdb\xo2'#Saves location of hotspot points with their unique id
		arcpy.SpatialJoin_analysis(pointInPolygon, y, z)#Combines unique id from polygon to points within.

		"""
		Code below takes points from all the activity locations and determines the stay time which is returned 
		as a dictionary.
		"""
		pointDic = {}
		with arcpy.da.SearchCursor(z, ['POINTZONE','UTC_DATE']) as cursor:
			for row in cursor:
				if row[0] in pointDic:
					pointDic[row[0]].append(row[1])
				else:
					pointDic[row[0]] = [row[1]]
		for i in pointDic:
			pointDic[i].sort()
		finalTimeValues = {}
		for i in pointDic:
			dtx = datetime.strptime(pointDic[i][0],self.datefmt)
			dty = datetime.strptime(pointDic[i][(len(pointDic[i]) - 1)],self.datefmt)
			homeTime = dtx.replace(hour=int(userProvidedHomeTime));
			amt = dty - dtx
			amtih = amt.total_seconds() / 3600
			if (dtx < homeTime < dty) and (amtih > userProvidedHomeDuration):#This code determines wheather or not a location is a home location based on time and duration.
				global identifiedHome
				identifiedHome = i
			finalTimeValues[i] = amtih
		self.timingDic = finalTimeValues
		arcpy.AddMessage(finalTimeValues)
	
	def resample(self, data_fc, sample_fc, sampleFrequency, filter='', 
				timestampField='utc_date', timestampFormat="%Y-%m-%d %H:%M:%S"):

		#In case data_fc1 is layer, replace it with path to FeatureClass
		desc_fc = arcpy.Describe(data_fc)
		data_fc = desc_fc.CatalogPath
		
		#retrieve fields of featureclass
		fieldnames = ['SHAPE@XY', timestampField] #build list of fieldnames (apart from OID and Shape fields)
				
		#get the whole FC as a Numpy array
		#return array is a structured numpy array, with column names corresponding to field names
		rawGpsArray = da.FeatureClassToNumPyArray(data_fc, fieldnames, where_clause=filter, skip_nulls=True)
		
		#sort the array according to utcdate
		rawGpsArray = np.sort(rawGpsArray, order=timestampField)
		
		#extract the utcdatetime column and convert it to 
		def _convertDT(utcdate): #define a datetime convesion function, return None in case of error
			try:
				return datetime.strptime(utcdate, timestampFormat)
			except:
				return None
			
		utcdateArray = rawGpsArray[timestampField]
		utcdateArray = list(map(_convertDT, utcdateArray))#Done due to conversion from python 2 to 3 (map no longer returns indexible array)
		
		#calculate difference between two samples and extract those that are over "sampleFreq"
		# deltaDateArray = np.diff(utcdateArray)
		#deltaDateArray = np.where(deltaDateArray == np.array(None), -1, deltaDateArray) #replace None value by -1
		#deltaDateArray = map(lambda t: np.asscalar(t).total_seconds(), deltaDateArray.flatten()) #transform timedelta object to seconds (float)
		keeperArray = [True]
		prevUtcdate = utcdateArray[0]
		for i in range(1, len(utcdateArray)-1):
			try:
				if (utcdateArray[i] - prevUtcdate).total_seconds() >= sampleFrequency:
					keeperArray.append(True)
					prevUtcdate = utcdateArray[i]
				else:
					keeperArray.append(False)
			except:
				keeperArray.append(False)
		
		#extract subsampled rows and copy to new Featureclass
		subsmplGpsArray = rawGpsArray[np.array(keeperArray)-1] #Done due to conversion from python 2 to 3
		if arcpy.Exists(sample_fc): arcpy.Delete_management(sample_fc) #overwrite does not work with da.NumPyToFC
		da.NumPyArrayToFeatureClass(subsmplGpsArray, sample_fc, 'SHAPE@XY', desc_fc.SpatialReference)
		#print subsmplGpsArray.dtype
		#print subsmplGpsArray

	def interpolate(self, data_fc, dummy_fc, maxDelay, freq, maxDropTime, maxDistance, verbose=True):
		'''
		Interpolate GPS points for droped periods
		Criteria:
			- if delay with previous point over maxDelay > start interpolation
			- if droped period is over maxDropTime > do not interpolate
			- except if distance from last valid point is below maxDistance
			
		Parameters:
			- data_fc: data feature class
			- dummy_fc: output featureclass with dummy points
			- maxDelay: max delay before patching points (seconds)
			- freq: frequency at which points are patched (seconds)
			- maxDropTime: max duration between two points to interpolate (hours)
			- maxDistance: max distance between two points to start interpolating without checking duration (in meters)
		'''
		#In case data_fc1 is layer, replace it with path to FeatureClass
		desc_fc = arcpy.Describe(data_fc)
		data_fc = desc_fc.CatalogPath
	
		#display data
		if verbose:
			arcpy.AddMessage( "   MAX DELAY: " + str(maxDelay) + "s")
			arcpy.AddMessage( "   FREQUENCY: " + str(freq) + "s")
			arcpy.AddMessage( "   MAX DISTANCE: " + str(maxDistance) + "m")
			arcpy.AddMessage( "   MAX DURATION: " + str(maxDropTime) + "h\n")
	
		
		#create a new featureclass to store dummy points
		arcpy.CreateFeatureclass_management(os.path.dirname(dummy_fc),os.path.basename(dummy_fc),
										 "POINT",data_fc,"SAME_AS_TEMPLATE","SAME_AS_TEMPLATE",desc_fc.SpatialReference)
		#check if field [dummy] exist, add it otherwise
		fieldList = arcpy.ListFields(dummy_fc)
		notFound = True
		for fld in fieldList:
			if fld.name.upper() == self.fld_interpolFlag.upper():
				notFound = False
		if notFound:
			arcpy.AddField_management(dummy_fc, self.fld_interpolFlag, "short")
			
		#check if standard DOP fields exist to fill later
		hdopFld = ''
		vdopFld = ''
		pdopFld = ''
		for fld in fieldList:
			if 'HDOP' == fld.name.upper():
				hdopFld = 'HDOP'
				continue
			if 'VDOP' == fld.name.upper():
				vdopFld = 'VDOP'
				continue
			if 'PDOP' == fld.name.upper():
				pdopFld = 'PDOP'
				continue
		#
		#Create Dummy Points for missing bouts
		#-------------------------------------
		
		#build progressor
		result = arcpy.GetCount_management(data_fc)
		nbRows = int(result.getOutput(0))
		p = int(math.log10(nbRows))
		if not p:
			p = 1
		increment = int(math.pow(10, p-1))
		arcpy.SetProgressor("step", "Creating dummy points...", 0, nbRows, increment)
	
		with da.SearchCursor(data_fc, ['OID@', 'SHAPE@XY', self.fld_date], 
							sql_clause=(None, 'ORDER BY {}'.format(self.fld_date))) as inCursor:
			#intitialize variables
			while True:
				try:
					inRow = inCursor.next()
					ptPrevXY = inRow[1] #retrieve XY coordinate as a tuple
					X,Y = ptPrevXY
					if (X<-180. or X>180.) or (Y<-90. or Y>90.): #skip invalid coordinates
						continue
					timePrevStr = inRow[2]
					timePrev = datetime.strptime(timePrevStr, self.datefmt)
					break
				except StopIteration:
					break
							
			#create an InsertCursor to push dummy point to output featureclass
			fields = ['SHAPE@XY', self.fld_date, self.fld_interpolFlag]
			if hdopFld: fields.append(hdopFld)
			if vdopFld: fields.append(vdopFld)
			if pdopFld: fields.append(pdopFld)
			outCursor = da.InsertCursor(dummy_fc, fields)
			
			#iterate over original points to detect gaps to be interpolated
			i = 0
			while True:
				if (i % increment) == 0:
					arcpy.SetProgressorPosition(i) #increment progressor
				try: #manage error arising from misformated date stamp
					inRow = inCursor.next()
					try:
						ptCurXY = inRow[1] #get coordinate tuple
						X,Y = ptCurXY
						if (X<-180. or X>180.) or (Y<-90. or Y>90.):
							i += 1
							continue
						timeCurr = datetime.strptime(inRow[2],self.datefmt)
					except:
						if verbose:
							message = "OID #{}:Unable to parse data [{}, {}] - skipping.".format(inRow[0], inRow[1], inRow[2])
							arcpy.AddWarning(message)
						i += 1
						continue
					
					#get delay between points...
					delayCurr = timeCurr - timePrev 
					#... and test interpolation condition
					if delayCurr > timedelta(seconds=int(maxDelay)):#Had to make instead of long due to python 3
						if verbose:
							message = 'Point #{} recorded at {} (UTC)'.format(inRow[0], inRow[2])
							message += '\n   Delay: {}\n   Distance: {:.1f}m'.format(delayCurr,
																					self._getGreatCircleDistance(ptPrevXY, ptCurXY))
							arcpy.AddMessage(message)
						if (delayCurr < timedelta(seconds=3600*long(maxDropTime))):
							if verbose:
								message = "   Delay below threshold ({:.1f}h), interpolating...".format(maxDropTime)
								arcpy.AddMessage(message)
							interpolatedTimeLap = long(freq) 
							while interpolatedTimeLap < delayCurr.total_seconds() :
								ptNewXY = self._createDummyPoint(ptPrevXY, ptCurXY, delayCurr.total_seconds(),interpolatedTimeLap)
								timestamp = timePrev + timedelta(seconds=interpolatedTimeLap)
								outRow = (ptNewXY, timestamp.strftime(self.datefmt), 1)
								#if any DOP field exists, dummy value of 1 is instered
								if hdopFld: outRow += (1,)
								if vdopFld: outRow += (1,)
								if pdopFld: outRow += (1,)
								outCursor.insertRow(outRow)
								interpolatedTimeLap = interpolatedTimeLap + long(freq)
						elif  (self._getGreatCircleDistance(ptPrevXY, ptCurXY) < float(maxDistance)):
							if verbose:
								message = "   Distance below threshold ({:.1f}m), interpolating...".format(maxDistance)
								arcpy.AddMessage(message)
							interpolatedTimeLap = long(freq) 
							while interpolatedTimeLap < delayCurr.total_seconds() :
								ptNewXY = self._createDummyPoint(ptPrevXY, ptCurXY, delayCurr.total_seconds(),interpolatedTimeLap)
								timestamp = timePrev + timedelta(seconds=interpolatedTimeLap)
								outRow = (ptNewXY, timestamp.strftime(self.datefmt), 1)
								#if any DOP field exists, dummy value of 1 is instered
								if hdopFld: outRow += (1,)
								if vdopFld: outRow += (1,)
								if pdopFld: outRow += (1,)
								outCursor.insertRow(outRow)
								interpolatedTimeLap = interpolatedTimeLap + long(freq)
						else:
							if verbose:
								message = "   Conditions not met, no interpolation!"
								arcpy.AddWarning(message)
					timePrev = timeCurr
					ptPrevXY = ptCurXY
					i += 1
				except StopIteration:
					break
			del outCursor	
			arcpy.ResetProgressor()

		#Append DummyPoint FeatureClass
		arcpy.Append_management(data_fc, dummy_fc,"NO_TEST")
	
	def gatherCriticalPoints(self, GpsFC, hotspotFC, scratchWS, verbose=True):
		'''
		Extract hotspots from GPS track by identifying the local peaks 
		on a kernel density surface built from the GPS fixes. 
		'''
		self.wkp = os.path.dirname(hotspotFC)
		self.cellSize = self.cellSize * 360. / self.earthLength
		self.kernelBandwidth = self.kernelBandwidth * 360. / self.earthLength
		ext = u''
		self.kernelRST = os.path.splitext(hotspotFC)[0] + u'_kernel' + ext
		result = sa.KernelDensity(GpsFC, 'NONE', self.cellSize, self.kernelBandwidth, 'SQUARE_MAP_UNITS')
		kdesc = arcpy.Describe(result)
		kext = kdesc.Extent

	
		if math.isinf(kext.height) or math.isinf(kext.width):
# 			#invalid raster due to clustered points below raster reolsution
# 			#in that case,  we now calculate the centroid of the points
			self._tinyExtentCase(GpsFC, hotspotFC)
			return -1
		else:
			result.save(self.kernelRST)

		arcpy.AddMessage(u'Redefine AOI')
		try:
			#v0.2.5: it may crash if kernel raster is invalid due to clustered points, in which
			# case we swirch to TinyExtentProc 
			self.aoiKernelRST = sa.SetNull(self.kernelRST,self.kernelRST,"value = 0")
			#		aoiKernelRST.save(arcpy.CreateScratchName('tmp', 'aoiKernel', 'RasterDataset', os.path.dirname(hotspotFC)))
		except:
			#v0.2.5 : change procedure-> we now calculate the centroid of the points
			self._tinyExtentCase(GpsFC, hotspotFC)
			return -1

		'''INVERT KERNEL SO THAT PEAKS BECOME PITS'''
		arcpy.AddMessage(u'Invert kernel')
		self.invertedKernelRST = sa.Negate(self.aoiKernelRST)

		'''FILTERS OUTLIERS'''
		arcpy.AddMessage(u'Filter low values')
		kMean = arcpy.GetRasterProperties_management(self.aoiKernelRST, "MEAN") #get mean of kernel, for all values > 0
		kMean = str(kMean.getOutput(0)) #str is needed to be sure that resultClass output is well formated
		self.filteredKernelRST = sa.SetNull(self.kernelRST, self.invertedKernelRST, "value < {}".format(kMean))
		
		'''FIND FLOW DIRECTION (WATERSHED FUNCTION)'''
		arcpy.AddMessage(u'Find flow direction')
		self.flowDirRST = sa.FlowDirection(self.filteredKernelRST, "NORMAL", "")

		'''DETERMINE SINKS (WATERSHED FUNCTION)'''
		arcpy.AddMessage(u'Find sinks')
		self.sinkRST = sa.Sink(self.flowDirRST)
	
		'''CONVERT RASTER SINKS TO POINT FEATURECLASS'''
		''' THESE ARE THE HOTSPOTS '''
		arcpy.AddMessage(u'Convert sinks to locations')

		#convert raster to polygon featureclass
		try:
			outPolyFC = arcpy.CreateScratchName("tmpg_","sinkpoly","FeatureClass",scratchWS)
			arcpy.RasterToPolygon_conversion(self.sinkRST, outPolyFC, "NO_SIMPLIFY", "VALUE")
		except:
			#an error happen if no polygons can be created
			#in that case,  we now calculate the centroid of the points,
			#using the TinyExtent as it is the most probable source of error
			self._tinyExtentCase(GpsFC, hotspotFC)
			return -1

		#dissolve polgygons that have same sink id
		dissolvePolyFC = arcpy.CreateScratchName("tmph_","dissolve","FeatureClass",scratchWS)
		try:
			arcpy.Dissolve_management(outPolyFC,dissolvePolyFC,"grid_code","#","MULTI_PART","DISSOLVE_LINES")
		except:
			arcpy.Dissolve_management(outPolyFC,dissolvePolyFC,"GRIDCODE","#","MULTI_PART","DISSOLVE_LINES")

		self.identifyTime(gpsFC,dissolvePolyFC)
		bufferedPolygon = r'C:\DALKAY.gdb\p85'#this is the location of the disolved polygon makes this edition of code more accurate than original.
		#calculate centroid of resulting polygons
		arcpy.FeatureToPoint_management(bufferedPolygon, hotspotFC, "CENTROID")
		arcpy.AddField_management(hotspotFC,self.fld_locid,"LONG")
		
		#clean up
		flds = arcpy.ListFields(hotspotFC)
		for fld in flds:
			if (fld.type not in ['OID','Geometry']) and (fld.name != self.fld_locid):
				arcpy.DeleteField_management(hotspotFC,fld.name)
		arcpy.Delete_management(outPolyFC)
		arcpy.Delete_management(dissolvePolyFC)


class SpatialK(object):
	"""
	This class is responsible for finding certain variables (ex) K) that relate to Spatial K.
	This class also utilizes all the variables it was given/found and performs the DalK calculation.
	"""

	def __init__(self):
		self.pointDicMask = {} #Saves the number of scrambled points within radius of original point.

	def calculation (self, locationOfSavedAL,maskedPoint):
		"""
		Responsible for collecting total false points ie the k value in 1/K for spatial K.
		"""
		
		actualPoint = locationOfSavedAL
		arcpy.AddField_management(actualPoint, "POINTZONE", "LONG")#Adds new field(variable) to attribute table of polygons in y.
		fields = ['POINTZONE'] #Unique identifiers for point within certain polygon.
		uniquePointIdentifier = 0
		with arcpy.da.UpdateCursor(actualPoint, fields) as cursor:
			for row in cursor:
				row[0] = uniquePointIdentifier
				cursor.updateRow(row) 
				uniquePointIdentifier += 1
		x = 'C:\DALKAY.gdb\goon3'#Sets the location of the buffer.

		arcpy.Buffer_analysis(actualPoint, x, str(stockBuffer)+" meters","FULL","ROUND")#Creates buffer around point at a distance that is preset to be close enough to mask original point.


		z = 'C:\DALKAY.gdb\outputOfSpatialJoin'# Stores Buffered polygons with unique id's.
		z2 = 'C:\DALKAY.gdb\outputOfPointsInRadi'# Stores points with unique ids that lay within buffered polygons.
		arcpy.SpatialJoin_analysis(x, actualPoint, z)
		arcpy.SpatialJoin_analysis(maskedPoint, z, z2)
		#Code below counts instances of dummy points close to real point (k).
		pointDic = {}#Saved dictionary of number of dummy points corresponding to a key which corresponds to the real point.
		with arcpy.da.SearchCursor(z2, ['POINTZONE']) as cursor:
			for row in cursor:
				if str(row[0]) != 'None':
					if row[0] in pointDic:
						pointDic[row[0]] = pointDic[row[0]] + 1
					else:
						pointDic[row[0]] = 1
		arcpy.AddMessage(pointDic)
		self.pointDicMask = pointDic 

	def dalKAnonymity(self,timingDic):
		"""
		This function calculates Dal K Anonymity based on self.pointDicMask which is the spatial k value and the timings 
		inputed through timingDic.
		"""
		arcpy.AddMessage("KARMA")
		arcpy.AddMessage(identifiedHome)
		arcpy.AddMessage(timingDic)
		arcpy.AddMessage(self.pointDicMask)

		saveHomeTime = [] #Holds time and position of a persons house.
		saveOtherTime = 0 #Saves part of Dalk formula that deals with calculating disclosure risk based on time of non home locations.
		arcpy.AddMessage(timingDic.keys())
		for i in timingDic.keys():
			#First if deals with identifying home.
			if identifiedHome == i:
				saveHomeTime = [i,float(timingDic[i])]
				arcpy.AddMessage(saveHomeTime)
			#else deals with non home areas
			else:
				saveOtherTime += (float(timingDic[i])/24) * (1 / self.pointDicMask[i])
				
		#finalCalc deals with computing part of Dalk that requires the home values and then it combines it with the non home calculations.
		finalCalc = saveOtherTime * (1 - (1/self.pointDicMask[saveHomeTime[0]])) + (1*(1/self.pointDicMask[saveHomeTime[0]]))
		arcpy.AddMessage('DAl K-Anonymity is:'+ str(finalCalc))

class Gaussian_Displacement(object):
	"""
	This class is used to perform a GPS point anonymizing algorithm called Gaussian Displacement which essentially
	randomizes input points by a provided radius.
	"""
	def Gaussian_Displacement(self, input_fc,workspace,radius,input_id,output_fc_name):
	    """
	    This function performs Gaussian Displacement a GPS point anonymizing algorithm.
	    """
	    """ 
	    Code below creates or stores working parameters.
	    """
	    env.workspace = workspace
	    arcpy.AddMessage(env.workspace)
	    env.overWriteOutput = True
	    # Create the array object and the point object
	    desc = arcpy.Describe(input_fc)
	    shapename = desc.shapeFieldName
	    old_cursors = arcpy.SearchCursor(input_fc)
	    arcpy.CreateFeatureclass_management(env.workspace, output_fc_name, "POINT", "", "", "", desc.spatialReference)
	    arcpy.AddField_management(output_fc_name, "UTC_DATE", "TEXT")
	   
	    xy = [] #Saves a list of randomized points.
	    # Get the coordinates of the original and masked points

	    """
	    for loop block looks through original data file and extracts the xy coordinates.
	    After the xy coordinates are extracted they are then randomized based on radius.
	    """
	    for old_cursor in old_cursors:
	        feat = old_cursor.getValue(shapename)
	        Date = old_cursor.getValue('UTC_DATE')
	        old_point = feat.getPart()
	        old_X = old_point.X
	        old_Y = old_point.Y
	        random_angle = random.uniform(0.0, math.pi*2)
	        random_number = random.gauss(0,1)
	        hypothenuse = random_number * radius
	        delta_X = (math.cos(random_angle)) * hypothenuse
	        delta_Y = (math.sin(random_angle)) * hypothenuse
	        new_X = old_X + delta_X
	        new_Y = old_Y + delta_Y
	        new_point = arcpy.Point(new_X,new_Y)
	        Distance = math.sqrt(math.pow(delta_X, 2) + math.pow(delta_Y, 2))
	        xy.append((str(Date), arcpy.Point(new_X, new_Y)))
	        
	    """
	    Block of code below creates a new .gdb file called Gaussian_Displacement then adds 
	    the randomized points created previous.
	    """
	    goon = r'C:\DALKAY.gdb\Gaussian_Displacement'
	    cursor = arcpy.da.InsertCursor(goon, ['UTC_DATE',"SHAPE@XY"])
	    for uv in xy:
	        cursor.insertRow(uv)
	    #Create the new feature class
	    del old_cursors


"""
Block of code below runs the program with user specified data from the GUI.
"""

if __name__ == '__main__':
	
	pwd_OutFolder = r'C:\DALKAY.gdb'
	if not os.path.exists(pwd_OutFolder): #Checks if file exists if not it creates new folder
		arcpy.CreateFileGDB_management("C:", "DALKAY.gdb", "9.2")

	if str(yourOwnTemplate) == "":
		radii = (float(arcpy.GetParameter(5)))/100000 # Convert KM to degrees for Gaussian Displacement.
		x = Gaussian_Displacement()
		x.Gaussian_Displacement(spatial_ref,pwd_OutFolder,radii,"ID",'Gaussian_Displacement')
		maskingTechniqueInput = r'C:\DALKAY.gdb\Gaussian_Displacement'
		gpsFC = maskingTechniqueInput
		filePathForHotspotSave = r'C:\DALKAY.gdb\maskingTechnique'  
		hotspotFC = filePathForHotspotSave
		scratchWS = r'C:\Users\sajee\Documents\ArcGIS\Default.gdb'
		kBW = kernalBand
		verbose = True
		cCore = computationalCore(kBW)
		cCore.gatherCriticalPoints(gpsFC,hotspotFC,scratchWS,verbose)
		arcpy.AddMessage(spatial_ref)
		yourOwnTemplate = r'C:\DALKAY.gdb\Gaussian_Displacement' 

	gpsFC = spatial_ref
	filePathForHotspotSave = 'C:\DALKAY.gdb\\' + hotspotLocationsName 
	hotspotFC = filePathForHotspotSave
	scratchWS = r'C:\DALKAY.gdb'
	kBW = kernalBand
	verbose = True
	activityLocation = computationalCore(kBW)
	activityLocation.gatherCriticalPoints(gpsFC,hotspotFC,scratchWS,verbose)
	arcpy.AddMessage(spatial_ref)
	if identifiedHome == []:
		arcpy.AddWarning("The program could not be completed due to the anonymized data not having any points in the DALK Buffer. Try adjusting (increasing) the DALK Buffer or reducing the level of anonymization/Stock Radius.")
	
	else:
		spatialK = SpatialK()
		spatialK.calculation(filePathForHotspotSave,yourOwnTemplate)
		spatialK.dalKAnonymity(activityLocation.timingDic)
