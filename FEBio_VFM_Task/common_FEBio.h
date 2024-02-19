#pragma once

#include <FECore/FECoreTask.h>
#include <FECore/log.h>
#include <FECore/Callback.h>
#include <FECore/FEInitialCondition.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FESurface.h>
#include <FECore/FEFacetSet.h>
#include <FECore/DataStore.h>
#include <FECore/FEPlotData.h>
#include <FECore/FEPlotDataStore.h>
#include <FECore/DumpFile.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEModel.h>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include <FECore/FESolidDomain.h>
#include <FEBioMech/FEElasticMaterial.h>
#include <FEBioMech/FEIsotropicElastic.h>
#include <FEBioMech/FEPressureLoad.h>
#include <FECore/FELoadController.h>
#include <FECore/FELoadCurve.h>
#include <FEBioMech/FEViscoElasticMaterial.h>
#include <FEBioMech/FENeoHookean.h>
#include <FEBioMech/FEElasticSolidDomain.h>

#ifndef PI
#define PI 3.141592653589793
#endif

void write_to_log_2(FEModel* fem, const std::string& logs_string, std::ofstream& outFile);