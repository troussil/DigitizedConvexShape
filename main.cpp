#include <iostream>

//We use CLI11 to manage options
//see https://github.com/CLIUtils/CLI11
//and https://cliutils.gitlab.io/CLI11Tutorial/
#include "ext/CLI11.hpp"

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include <DGtal/io/viewers/PolyscopeViewer.h>

#include "managers/SurfaceManager.ih"
#include "managers/ConvexHullManager.ih"
#include "helpers/DigitizedCvxPolygon.ih"

using namespace std;
using namespace DGtal;
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;

//--------------------------------------------------------------------------------
//-- main 
//--------------------------------------------------------------------------------
int main (int argc, char* argv[]) {

  // parse command line using CLI
  CLI::App app;
  
  app.description("Demo: 'Geometry of Gauss digitized convex shapes'\n");

  double gridstep = 0.5;
  string polynomial = "ellipsoid";
  
  vector<string> polynomialList = { "ellipsoid", "rcube", "sphere9" }; 
  stringstream ss;
  ss << "A polynomial";
  if (polynomialList.size() > 0) {
    ss << " (or a name in the following list: "
       << *polynomialList.begin(); 
    for (auto it = ++polynomialList.begin(), itEnd = polynomialList.end();
	 it != itEnd; ++it) {
      ss << ", " << *it; 
    }
    ss << ")"; 
  }

  auto p_opt = app.add_option("-p", polynomial, ss.str(), true);
  auto g_opt = app.add_option("-g", gridstep, "Grid step", true);
 
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI

  // parameters
  auto params    = SH3::defaultParameters() | SHG3::defaultParameters();
  params( "polynomial", polynomial )( "gridstep", gridstep );
  
  //---------------------------------------------------------------------------- 
  trace.beginBlock("digital surface extraction"); 
  auto implicit_shape  = SH3::makeImplicitShape3D( params );
  auto digitized_shape = SH3::makeDigitizedImplicitShape3D( implicit_shape, params );
  auto K               = SH3::getKSpace( params );
  auto surfacePtr      = SH3::makeDigitalSurface( digitized_shape, K, params );

  // surface as mesh
  SurfaceManager<Z3i::KSpace> surface;
  surface.buildMesh(K, surfacePtr);
  trace.endBlock();

  //---------------------------------------------------------------------------- 
  trace.beginBlock("convex hull computation"); 
  ConvexHullManager<Z3i::KSpace> convexHull;
  convexHull.build(surface.myVertices);
  trace.endBlock();

  //---------------------------------------------------------------------------- 
  trace.beginBlock("Expected normals...");
  auto expectedNormals = SHG3::getNormalVectors( implicit_shape, K, SH3::getSurfelRange(surfacePtr), params );
  trace.endBlock();
  
  //---------------------------------------------------------------------------- 
  trace.beginBlock("normal estimation");
  set<Z3i::Cell> out, in; 
  vector<Z3i::RealVector> normals(surface.mySurfels.size(),
				  Z3i::Vector::diagonal(0.0));
  
  for (int k = 0; k < convexHull.myFaces.size(); ++k) {
    //trace.info() << "# Face #" << k << endl; 
    //face vertices
    vector<Z3i::Point> faceVertices;
    for (auto const& idx: convexHull.myFaces.at(k)) {
      faceVertices.push_back( convexHull.myVertices[idx] ); 
    }
    DigitizedCvxPolygon<Z3i::KSpace> cp(K, faceVertices);
    //face normal
    Z3i::RealVector n = cp.getUnitNormal(); 
    //face digitization
    set<Z3i::Cell> surfelSet;
    cp.digitize(inserter(surfelSet, surfelSet.end()));
    //trace.info() << "# " << surfelSet.size() << " cells" << endl; 

    //assign a normal to every surfel in the surface
    for (auto const& s: surfelSet) {
      auto it = surface.mySMap.find(s); 
      if (it != surface.mySMap.end()) { //if in surface
	in.insert(s); 
	Z3i::RealVector sum = normals.at(it->second) + n;
	normals.at(it->second) = sum.getNormalized(); 
      } else {
	out.insert(s); 
      }
    }
  }
  trace.info() << "# " << in.size() << " cells / "
	       << surface.mySurfels.size() << " recovered" << endl; 
  trace.endBlock();

  //---------------------------------------------------------------------------- 
  trace.beginBlock("Statistics..."); 
  ASSERT(normals.size() == expectedNormals.size());
  auto diff = SHG3::getVectorsAngleDeviation(expectedNormals, normals);
  auto stat = SHG3::getStatistic(diff);

  trace.info() << "#gridstep nbVtx min avg max stdVar" << endl; 
  trace.info() << gridstep << " "
	       << normals.size() << " "
	       << stat.min() << " "
	       << stat.mean() << " "
	       << stat.max() << " "
	       << sqrt((double) stat.variance()) << " "
	       << endl;
  trace.endBlock(); 

  //---------------------------------------------------------------------------- 
  // polyscope session
  PolyscopeViewer viewer;
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;

  //structures
  auto sPtr = polyscope::registerSurfaceMesh("Surface",
					     surface.myVertices,
					     surface.myFaces);
  sPtr->setEdgeWidth(1.0);

  auto cPtr = polyscope::registerSurfaceMesh("Convex hull",
					     convexHull.myVertices,
					     convexHull.myFaces);
  cPtr->setEdgeWidth(2.0);
  cPtr->setTransparency(0.75); 
  cPtr->setEnabled(false); 

  //quantities
  sPtr->addFaceVectorQuantity("Estimated Normals", normals); 
  sPtr->addFaceVectorQuantity("Expected Normals", expectedNormals);
  auto dPtr = sPtr->addFaceScalarQuantity("Deviation", diff);
  dPtr->setEnabled(true);

  viewer.show();
  
  return 0;
}
