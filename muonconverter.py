#This script makes pytables, so that your plotter can understand your data
from icecube import icetray
from icecube import dataio,tableio,phys_services,hdfwriter,dataclasses,millipede,simclasses
from icecube.phys_services import converters
from icecube.phys_services.converters import I3EventInfoConverterFromRecoPulses 
from I3Tray import I3Tray
from icecube.hdfwriter import I3HDFTableService
from icecube.tableio import I3TableWriter
from icecube.dataclasses import converters
from icecube.dataclasses.converters import I3DoubleConverter
from glob import glob
import os
import sys

if len(sys.argv) < 3:
	print 'Usage: %s output.i3 input1.i3 [input2.i3] ...'
	sys.exit(1)

files = sys.argv[2:]
output = sys.argv[1]

tray = I3Tray()
#I3 file with stuff I want in a table
tray.AddModule('I3Reader','reader',FilenameList = files)
#Create the pytable
hdf5 = I3HDFTableService(output)

#Make Root and Pytable Writers with keys->Frames I want in the trees, SubEventStreams has to be set for Q frames for some reason... 
#tray.AddModule(I3TableWriter,'writer', tableservice = [hdf,roottree], SubEventStreams = ["nullsplit","in_ice"], keys=[things I want]

tray.AddModule(lambda frame: frame.Has('CylinderChargesInVeto'))

# Genie InteractionParticle Fix for Low Energy
def addMCTreeItem(frame):
	mostEnergetic = dataclasses.get_most_energetic_muon(frame["I3MCTree"])
	if  mostEnergetic:
		frame['MostEnergMuon'] = mostEnergetic
	else:
		mostEnergetic=dataclasses.get_most_energetic_primary(frame["I3MCTree"]) 
		frame['MostEnergMuon'] = mostEnergetic


tray.AddModule(addMCTreeItem,'addMCTree')

# All the things I want to analyze
tray.AddModule(I3TableWriter,'writer', tableservice = [hdf5], SubEventStreams = ["nullsplit","in_ice","InIceSplit"],keys=['ClosestDOMS','CylinderChargesInDC','CylinderChargesInVeto','CylinderChargesInDust','HitDistancesFromDCTrack','HitDistancesFromVetoTrack','HitDistancesInDust','DecayParticle','InIceNu','InteractionVertex','PrimaryNu','I3MCPrimary','InteractionParticle','MostEnergMuon', 'I3MCWeightDict', 'MMCTrackList', 'NclosestDOMS', 'HoerandelWeight', 'GaisserH4aWeight', 'MuonWeight','FilterMask'])

#Make everything nice
tray.AddModule('TrashCan','yeswecan')
#Run It 
tray.Execute()
#Done! 
tray.Finish()
