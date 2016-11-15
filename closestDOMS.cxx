#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <icetray/I3Tray.h>
#include "phys-services/I3Calculator.h"
#include <dataclasses/physics/I3EventHeader.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/geometry/I3Geometry.h>
#include <dataclasses/geometry/I3OMGeo.h>
#include <dataclasses/I3MapOMKeyMask.h>
#include <dataclasses/I3Vector.h>
#include <dataclasses/physics/I3MCTreeUtils.h>

static double trackCylinderRadius_;
static std::string spulses_;

//Check to see if the particle is a Muon
bool IsMuon(const I3Particle& p1){
	return p1.GetType()==I3Particle::MuPlus || p1.GetType()==I3Particle::MuMinus;
}

//Find the Most Energetic Muon, or if there are no muons
//return the most energetic particle
bool BestMuon(const I3Particle& p1, const I3Particle& p2){
	if(IsMuon(p1) && !IsMuon(p2)){
	return true;
	}
	if(IsMuon(p2) && !IsMuon(p1)){
	return false;
	}
	return p1.GetEnergy() > p2.GetEnergy();
}


//This is our main process. It collects and calculates all the values of interest
//from the particle interaction, and then returns them into the frame.
//The point is to find the DOMs (Digital Optical Modules or light detectors)
//within a cylinder of radius defined by the user and package thier information
//based on what part of the detector they are located in for further study. 
void CalculateCylinderValues(boost::shared_ptr<I3Frame> frame){
//Define values of interest
	I3VectorDouble hitDistancesFromVetoTrack;
	I3VectorDouble hitDistancesFromDCTrack;
	I3VectorDouble hitDistancesInDust;
	I3VectorDouble cylinderChargesInVeto;
	I3VectorDouble cylinderChargesInDC;
	I3VectorDouble cylinderChargesInDust;
	I3VectorDouble ClosestDOMS;
	I3VectorDouble NclosestDOMS;

//Ensure that the frame has the pulses specified for our calculations
	if(! frame->Has(spulses_))
		return;

	const std::string			stringtheorygeometry="I3Geometry";
	boost::shared_ptr<const I3Geometry>			geometry = frame->Get<boost::shared_ptr<const I3Geometry> >(stringtheorygeometry);
	boost::shared_ptr<const I3RecoPulseSeriesMap>	pulsesMap = frame->Get<boost::shared_ptr<const I3RecoPulseSeriesMap> >(spulses_);
	boost::shared_ptr<const I3MCTree>				tree = frame->Get<boost::shared_ptr<const I3MCTree> >("I3MCTree");
//The below section is for running on simulation files with Bad DOM lists that need to be specified
//	boost::shared_ptr<const BadDomsList>			baddoms = frame->Get<boost::shared_ptr<const I3Vector<OmKey>>("BadDomsList");
	
// When running on Corsika simulation this is needed to ensure the frame is not empty
//	if(frame->Has(CorsikaWeightMap) && frame['CorsikaWeightMap']['Multiplicity']==0)
//		return;

//Initializing variables
	double nDC = 0.;
	double nDust = 0.;
	double nVeto = 0.;
	double nqDC = 0.;
	double nqDust = 0.;
	double nqVeto = 0.;
	double smallestDC = 10000.;
	double smallestDust = 10000.;
	double smallestVeto = 10000.;
	double largest = 0.;
	double largestd = 0.;


//Find the most energetic muon, or if there are no muons, the most energetic
//particle in the event, save its energy and angle of incidence (zenith)
	I3Particle particle = *I3MCTreeUtils::GetBest(*tree, BestMuon);
	double energy = particle.GetEnergy();
	double zenith = particle.GetZenith();
	// Iterate over all DOMs.
	for(I3OMGeoMap::const_iterator geoMapCIter = geometry->omgeo.begin();
		geoMapCIter != geometry->omgeo.end(); ++geoMapCIter){
		const OMKey&             omKey  = geoMapCIter->first;
		const I3Position& omPos = geoMapCIter->second.position;
		double stringnom = omKey.GetString();
		double omnom = omKey.GetOM();

		// Check if the DOM is inside the track cylinder.
		double domCadToTrack = I3Calculator::ClosestApproachDistance(particle, omPos);
		if(std::isnan(domCadToTrack) || domCadToTrack > trackCylinderRadius_)
			continue;

		//Check if DOM is Bad for files with BadDomLists
		//bool isBad = (std::find(baddoms, end, omKey) != end);
		//if(isBad)
		//	continue;

		// Check if DOM is in DeepCore sub-detector, and if it is the dust layer where the ice has poor optical qualities.
		int deepCoreStrings [15] = {26,27,35,36,37,45,46,79,80,81,82,83,84,85,86};
		int* end = deepCoreStrings+15;
		bool isDC = (std::find(deepCoreStrings, end, omKey.GetString()) != end && omPos.GetZ() < -150.0);//omnom > 10);
		bool isDust = (omPos.GetZ() < 0. && omPos.GetZ() > -150.0);
	
		// Record the distance of this hit (light signal) from the track.
		// Keep track of the closest hit to the track in each detector region for the event.
		//double hitDistAlongTrack = I3Calculator::DistanceAlongTrack(particle, omPos);
		if(isDC){
			hitDistancesFromDCTrack.push_back(domCadToTrack);
			nDC++;
			if(domCadToTrack < smallestDC){
				smallestDC = domCadToTrack;
			}
		}
		else if(isDust){
			hitDistancesInDust.push_back(domCadToTrack);
			nDust++;
			if(domCadToTrack < smallestDust){
				smallestDust = domCadToTrack;
			}
		}
		else{
	        	hitDistancesFromVetoTrack.push_back(domCadToTrack);
			nVeto++;
			if(domCadToTrack < smallestVeto){
				smallestVeto = domCadToTrack;
			}
		}	

		I3RecoPulseSeriesMap::const_iterator pulsesMapCIter = pulsesMap->find(omKey);


		double qTotDom = 0.;
		double qTotDCDom=0.;

		// Sanity check: Check if there are pulses at all.
		if(pulsesMapCIter == pulsesMap->end())
		// Push values for empty event
		{
			if(isDC){
				cylinderChargesInDC.push_back(qTotDCDom);
				ClosestDOMS.push_back(1.0);
				ClosestDOMS.push_back(stringnom);
				ClosestDOMS.push_back(omnom);
				ClosestDOMS.push_back(domCadToTrack);
				ClosestDOMS.push_back(qTotDCDom);
				ClosestDOMS.push_back(energy);
				ClosestDOMS.push_back(zenith);
			}
			else if(isDust){
				cylinderChargesInDust.push_back(qTotDom);
				ClosestDOMS.push_back(3.0);
				ClosestDOMS.push_back(stringnom);
				ClosestDOMS.push_back(omnom);
				ClosestDOMS.push_back(domCadToTrack);
				ClosestDOMS.push_back(qTotDom);
				ClosestDOMS.push_back(energy);
				ClosestDOMS.push_back(zenith);
			}
			else{
				cylinderChargesInVeto.push_back(qTotDom);
				ClosestDOMS.push_back(0.0);
				ClosestDOMS.push_back(stringnom);
				ClosestDOMS.push_back(omnom);
				ClosestDOMS.push_back(domCadToTrack);
				ClosestDOMS.push_back(qTotDCDom);
				ClosestDOMS.push_back(energy);
				ClosestDOMS.push_back(zenith);
			}
			continue;
		}


		const I3RecoPulseSeries& pulses = pulsesMapCIter->second; //FIX ME

		// Iterate over the pulses of the DOM's pulse series.
		for(I3RecoPulseSeries::const_iterator pulsesCIter = pulses.begin();
		pulsesCIter != pulses.end();
		++pulsesCIter
		){
			const I3RecoPulse& pulse = *pulsesCIter;

			// Get the charge of the pulse and count it for the total charge of
			// the OM, if the charge has a physical value.
			double pulseCharge = pulse.GetCharge();
			if(! (std::isnan(pulseCharge) || std::isinf(pulseCharge))
			){
				if(isDC){
					qTotDCDom += pulseCharge;
				}
				else{
					qTotDom += pulseCharge;
				}
			}
		}//END FOR I3RecoPulseSeries
		//Record basic information for each registered hit.
		//Charge of the hit, what string the detector which registers the hit is on, the dom id that registers the hit,
		//The energy and zenith of the interaction that produced the hit, the distance from the DOM to the particle that
		//produced the light, and a numerical indicator or which zone the dom is in (1 for the Deep Core subarray, 
		//3 for the dust layer, and 0 for the rest of the detector, which acts as a veto region for this analysis).
		if(isDC){
			cylinderChargesInDC.push_back(qTotDCDom);
			ClosestDOMS.push_back(1.0);
			ClosestDOMS.push_back(stringnom);
			ClosestDOMS.push_back(omnom);
			ClosestDOMS.push_back(domCadToTrack);
			ClosestDOMS.push_back(qTotDCDom);
			ClosestDOMS.push_back(energy);
			ClosestDOMS.push_back(zenith);
			if(qTotDCDom > 0.){
				nqDC++;
				if(qTotDCDom > largest){
					largest = qTotDCDom;
					largestd = domCadToTrack;
				}
			}
		}
		else if(isDust){
			cylinderChargesInDust.push_back(qTotDom);
			ClosestDOMS.push_back(3.0);
			ClosestDOMS.push_back(stringnom);
			ClosestDOMS.push_back(omnom);
			ClosestDOMS.push_back(domCadToTrack);
			ClosestDOMS.push_back(qTotDom);
			ClosestDOMS.push_back(energy);
			ClosestDOMS.push_back(zenith);
			if(qTotDom > 0.){
				nqDust++;
				if(qTotDom > largest){
					largest = qTotDom;
					largestd = domCadToTrack;
				}
			}
		}
		else{
			cylinderChargesInVeto.push_back(qTotDom);
			ClosestDOMS.push_back(0.0);
			ClosestDOMS.push_back(stringnom);
			ClosestDOMS.push_back(omnom);
			ClosestDOMS.push_back(domCadToTrack);
			ClosestDOMS.push_back(qTotDom);
			ClosestDOMS.push_back(energy);
			ClosestDOMS.push_back(zenith);
			if(qTotDom > 0.){
				nqVeto++;
				if(qTotDom > largest){
					largest = qTotDom;
					largestd = domCadToTrack;
				}
			}
		}
	}//END FOR I3RecoPulseSeriesMap
	//Collect all the min/max/sum values for the entire event
	NclosestDOMS.push_back(nDC);
	NclosestDOMS.push_back(nDust);
	NclosestDOMS.push_back(nVeto);
	NclosestDOMS.push_back(nqDC);
	NclosestDOMS.push_back(nqDust);
	NclosestDOMS.push_back(nqVeto);
	NclosestDOMS.push_back(smallestDC);
	NclosestDOMS.push_back(smallestDust);
	NclosestDOMS.push_back(smallestVeto);
	NclosestDOMS.push_back(largest);
	NclosestDOMS.push_back(largestd);
	NclosestDOMS.push_back(energy);
	NclosestDOMS.push_back(zenith);
	//Output all the calculated values
	frame->Put("HitDistancesFromVetoTrack", boost::make_shared<const I3VectorDouble>(hitDistancesFromVetoTrack));
	frame->Put("HitDistancesFromDCTrack", boost::make_shared<const I3VectorDouble>(hitDistancesFromDCTrack));
	frame->Put("HitDistancesInDust", boost::make_shared<const I3VectorDouble>(hitDistancesInDust));
	frame->Put("CylinderChargesInVeto", boost::make_shared<const I3VectorDouble>(cylinderChargesInVeto));
	frame->Put("CylinderChargesInDC", boost::make_shared<const I3VectorDouble>(cylinderChargesInDC));
	frame->Put("CylinderChargesInDust", boost::make_shared<const I3VectorDouble>(cylinderChargesInDust));
	frame->Put("ClosestDOMS", boost::make_shared<const I3VectorDouble>(ClosestDOMS));
	frame->Put("NclosestDOMS", boost::make_shared<const I3VectorDouble>(NclosestDOMS));
}


//Input our user values, print error if there are not enough values, run over all the events/frames in the input file
// and produce an output file with the new calculated values.
//Includes exceptions for all the errors that were encountered during use.
int main(int argc, char* argv[]){
        if(argc<5){
                std::cerr << "Missing command line input: cylinder radius, input file, output file, pulse map" << std::endl;
                return(1);
        }
	trackCylinderRadius_= atof(argv[1]);
	spulses_= argv[4];
	
        I3::init_icetray_lib();
	try{
       		I3Tray tray;
        	tray.AddModule("I3Reader")("Filename", std::string(argv[2]));
        	tray.AddModule(CalculateCylinderValues);
		tray.AddModule("I3Writer")("Filename", std::string(argv[3]));

        	tray.Execute();
	}catch(boost::python::error_already_set& e){
		PyErr_Print();
	}
}
