#include "geometrycentral/pointcloud/point_cloud_heat_solver.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/transfer_functions.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::unique_ptr<pointcloud::PointCloud> cloud;
std::unique_ptr<pointcloud::PointPositionGeometry> pointGeom;
pointcloud::PointData<Vector3> pointPositions;

// Polyscope stuff
polyscope::SurfaceMesh* psMesh;
polyscope::PointCloud* psCloud;

enum MeshMode { Triangle = 0, Polygon, Points };
int MESH_MODE = MeshMode::Triangle;

// == Intrinsic triangulation stuff. All remeshing is performed on the manifold mesh.
std::unique_ptr<ManifoldSurfaceMesh> manifoldMesh;
std::unique_ptr<VertexPositionGeometry> manifoldGeom;
std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation> intTri;
std::vector<Curve> curvesOnManifold, curvesOnIntrinsic;
std::vector<Point> pointsOnManifold, pointsOnIntrinsic;
std::vector<std::vector<pointcloud::Point>> POINTS_CURVES;
std::vector<std::vector<Vertex>> POLYGON_CURVES;
float REFINE_AREA_THRESH = std::numeric_limits<float>::infinity();
float REFINE_ANGLE_THRESH = 25.;
int MAX_INSERTIONS = -1;

enum SolverMode { ExtrinsicMesh = 0, IntrinsicMesh };
int SOLVER_MODE = SolverMode::ExtrinsicMesh;
int CONSTRAINT_MODE = LevelSetConstraint::ZeroSet;

// Solvers & parameters
float TCOEF = 1.0;
float SOFT_WEIGHT = 1.0;
bool SOFT_CONSTRAINT = false;
bool SOLVE_AS_POINT_CLOUD = false;
bool EXPORT_RESULT = false;
bool VIZ = true;
bool VERBOSE, HEADLESS, IS_POLY, PIECEWISE;
std::unique_ptr<SignedHeatMethodSolver> SHM, intrinsicSolver;
std::unique_ptr<PointCloudHeatSolver> PCS;
VertexData<double> PHI_VERTICES;
CornerData<double> PHI_CORNERS;

// Program variables
std::string MESHNAME = "input mesh";
std::string MESH_FILEPATH, MESHROOT, OUTPUT_FILENAME;
std::string OUTPUT_DIR = "../export/";
bool USE_BOUNDS = false;
float LOWER_BOUND, UPPER_BOUND;

void callback() {

    if (ImGui::Button("Solve")) {

        if (MESH_MODE == MeshMode::Triangle) {
            if (SOLVER_MODE == SolverMode::ExtrinsicMesh) {
                SHM->setDiffusionTime(TCOEF);
                Vector<double> phi =
                    SHM->computeDistance(CURVES, POINTS, CONSTRAINT_MODE, SOFT_CONSTRAINT, SOFT_WEIGHT, PIECEWISE);

                if (phi.size() == mesh->nVertices()) {
                    PHI_VERTICES = VertexData<double>(*mesh, phi);
                    if (!HEADLESS) psMesh->addVertexSignedDistanceQuantity("GSD", PHI_VERTICES)->setEnabled(true);
                } else {
                    PHI_CORNERS = CornerData<double>(*mesh, phi);
                    if (!HEADLESS)
                        psMesh->addCornerScalarData("GSD", PHI_CORNERS)->setIsolinesEnabled(true)->setEnabled(true);
                }

            } else {
                // Create data structures, if they haven't been created already.
                if (intTri == nullptr && mesh->isManifold() && mesh->isOriented() &&
                    isSourceGeometryConstrained(CURVES, POINTS)) {
                    setIntrinsicSolver(*geometry, CURVES, POINTS, manifoldMesh, manifoldGeom, curvesOnManifold,
                                       pointsOnManifold, intTri, intrinsicSolver);
                }

                // Solve
                intrinsicSolver.reset(new SignedHeatMethodSolver(*intTri));
                intrinsicSolver->setDiffusionTime(TCOEF);
                Vector<double> phi = intrinsicSolver->solve(curvesOnIntrinsic, pointsOnIntrinsic, CONSTRAINT_MODE,
                                                            SOFT_CONSTRAINT, SOFT_WEIGHT, PIECEWISE);

                if (phi.size() == intTri->intrinsicMesh->nVertices()) {
                    PHI_VERTICES =
                        transferBtoA(*intTri, VertexData<double>(*(intTri->intrinsicMesh), phi), TransferMethod::L2);
                    if (!HEADLESS) psMesh->addVertexSignedDistanceQuantity("GSD", PHI_VERTICES)->setEnabled(true);
                } else {
                    // TODO: transfer & visualize CornerData
                }
            }

            if (EXPORT_RESULT) {
                exportCurves(geometry->vertexPositions, CURVES, POINTS, OUTPUT_DIR + "source.obj");
                if (phi.size() == mesh->nVertices()) {
                    exportSDF(*geometry, PHI_VERTICES, OUTPUT_FILENAME, USE_BOUNDS, LOWER_BOUND, UPPER_BOUND);
                } else {
                    exportSDF(*geometry, PHI_CORNERS, OUTPUT_FILENAME, USE_BOUNDS, LOWER_BOUND, UPPER_BOUND);
                }
            }

        } else if (MESH_MODE == MeshMode::POINTS) {
            // Reset point cloud solver in case parameters changed.
            PCS.reset(new PointCloudHeatSolver(*cloud, *pointGeom, TCOEF));
            pointcloud::PointData<double> phi = PCS->computeSignedDistance(POINTS_CURVES, CONSTRAINT_MODE);
            if (!HEADLESS) psCloud->addScalarQuantity("GSD", phi)->setIsolinesEnabled(true)->setEnabled(true);

            if (EXPORT_RESULT) {
                exportCurves(geometry->vertexPositions, CURVES, std::vector<Point>(), OUTPUT_DIR + "source.obj");
                exportSDF(pointGeom->positions, phi, OUTPUT_FILENAME);
            }
        } else if (MESH_MODE == MeshMode::Polygon) {
            // TODO: Set up polygon mesh solver
            // SHM->setDiffusionTime(TCOEF);
            // SHM->setSoftWeight(SOFT_WEIGHT);
            // Vector<double> phi = SHM->solve(CURVES, POINTS, CONSTRAINT_MODE, SOFT_CONSTRAINT);

            if (!HEADLESS) psMesh->addVertexSignedDistanceQuantity("GSD", phi)->setEnabled(true);
            if (EXPORT_RESULT) {
                exportCurves(geometry->vertexPositions, CURVES, POINTS, OUTPUT_DIR + "source.obj");
                exportSDF(*geometry, PHI_VERTICES, OUTPUT_FILENAME, USE_BOUNDS, LOWER_BOUND, UPPER_BOUND);
            }
        }
    }

    if (ImGui::TreeNode("Solve options")) {

        ImGui::InputFloat("tCoef", &TCOEF);

        ImGui::RadioButton("Constrain zero set", &CONSTRAINT_MODE, LevelSetConstraint::ZeroSet);
        ImGui::RadioButton("Constrain multiple levelsets", &CONSTRAINT_MODE, LevelSetConstraint::Multiple);
        ImGui::RadioButton("No levelset constraints", &CONSTRAINT_MODE, LevelSetConstraint::None);

        ImGui::Checkbox("Soft levelset constraint", &SOFT_CONSTRAINT);
        ImGui::InputFloat("soft weight", &SOFT_WEIGHT);

        ImGui::Checkbox("Export result", &EXPORT_RESULT);

        ImGui::TreePop();
    }

    if (MESH_MODE == MeshMode::TRIANGLE && mesh->isManifold() && intTri != nullptr) {
        ImGui::RadioButton("Solve on extrinsic mesh", &SOLVER_MODE, SolverMode::ExtrinsicMesh);
        ImGui::RadioButton("Solve on intrinsic mesh", &SOLVER_MODE, SolverMode::IntrinsicMesh);

        if (ImGui::TreeNode("Intrinsic mesh improvement")) {
            if (ImGui::Checkbox("Show intrinsic edges", &VIS_INTRINSIC_MESH)) {
                visualizeIntrinsicEdges(*intTri, *manifoldGeom, true);
            }

            if (ImGui::Button("Flip to Delaunay")) {
                intTri->flipToDelaunay();
                visualizeIntrinsicEdges(*intTri, *manifoldGeom, true);
            }

            ImGui::InputFloat("Angle thresh", &REFINE_ANGLE_THRESH);
            ImGui::InputFloat("Area thresh", &REFINE_AREA_THRESH);
            ImGui::InputInt("Max insert", &MAX_INSERTIONS);
            if (ImGui::Button("Delaunay refine")) {
                intTri->delaunayRefine(REFINE_ANGLE_THRESH, REFINE_AREA_THRESH, MAX_INSERTIONS);
                visualizeIntrinsicEdges(*intTri, *manifoldGeom, true);
            }

            ImGui::TreePop();
        }
    }
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("Solve for generalized signed distance.");
    args::HelpFlag help(parser, "help", "Display this help menu", {"help"});
    args::Positional<std::string> meshFilename(parser, "mesh", "A mesh file.");
    args::ValueFlag<std::string> inputFilename(parser, "input", "Input curve filepath", {"i", "input"});
    args::ValueFlag<std::string> outputFilename(parser, "output", "Output filename", {"o", "output"});

    args::Group group(parser);
    args::Flag points(group, "point", "Solve as point cloud.", {"pc", "points"});
    args::Flag piecewise(group, "piecewise", "Solve for a piecewise SDF.", {"L1", "piecewise"});
    args::Flag verbose(group, "verbose", "Verbose output", {"V", "verbose"});
    args::Flag headless(group, "headless", "Don't use the GUI.", {"h", "headless"});

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    if (!meshFilename) {
        std::cerr << "Please specify a mesh file as argument." << std::endl;
        return EXIT_FAILURE;
    }

    // Load mesh
    MESH_FILEPATH = args::get(meshFilename);
    DATA_DIR = getHomeDirectory(MESH_FILEPATH);
    OUTPUT_DIR = getHomeDirectory(OUTPUT_FILENAME);
    MESHROOT = polyscope::guessNiceNameFromPath(MESH_FILEPATH);

    OUTPUT_FILENAME = outputFilename ? args::get(outputFilename) : OUTPUT_DIR + "GSD.obj";
    HEADLESS = headless;
    VERBOSE = verbose;
    PIECEWISE = piecewise;

    if (points) {
        MESH_MODE = MeshMode::Points;
        readPointCloud(MESH_FILEPATH, pointPositions, *cloud, *pointGeom);
        PCS.reset(new PointCloudHeatSolver(*cloud, *pointGeom, TCOEF));
    } else {
        std::tie(mesh, geometry) = readSurfaceMesh(MESH_FILEPATH);
        if (!mesh->isTriangular()) MESH_MODE = MeshMode::Polygon;
    }

    // Load source geometry.
    if (inputFilename) {
        std::string filename = args::get(inputFilename);
        switch (MESH_MODE) {
            case (MeshMode::Triangle): {
                std::tie(CURVES, POINTS) = readInput(*mesh, filename);
                break;
            }
            case (MeshMode::Polygon): {
                POLYGON_CURVES = readCurveVertices(*mesh, filename);
                break;
            }
            case (MeshMode::Points): {
                POINTS_CURVES = readCurvePoints(*cloud, filename);
                break;
            }
        }
    }

    // Visualize data.
    if (!HEADLESS) {
        polyscope::init();
        polyscope::state::userCallback = callback;

        switch (MESH_MODE) {
            case (MeshMode::Triangle): {
                psMesh = polyscope::registerSurfaceMesh(MESHNAME, getGeom().inputVertexPositions,
                                                        getMesh().getFaceVertexList());
                psMesh->setAllPermutations(polyscopePermutations(getMesh()));
                displayInput(geometry->vertexPositions, CURVES, POINTS);
                break;
            }
            case (MeshMode::Polygon): {
                psMesh = polyscope::registerSurfaceMesh(MESHNAME, getGeom().inputVertexPositions,
                                                        getMesh().getFaceVertexList());
                displayInput(geometry->vertexPositions, POLYGON_CURVES);
                break;
            }
            case (MeshMode::Points): {
                psCloud = polyscope::registerPointCloud("point cloud", pointPositions);
                displayInput(pointGeom->positions, POINTS_CURVES);
                break;
            }
        }
        polyscope::show();
    } else {
        solve();
    }

    return EXIT_SUCCESS;
}