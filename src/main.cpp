#include "geometrycentral/pointcloud/point_cloud_heat_solver.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/transfer_functions.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

#include <chrono>
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

// == Geometry-central data
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::unique_ptr<pointcloud::PointCloud> cloud;
std::unique_ptr<pointcloud::PointPositionGeometry> pointGeom;
pointcloud::PointData<Vector3> pointPositions, pointNormals;

// Polyscope stuff
polyscope::SurfaceMesh* psMesh;
polyscope::PointCloud* psCloud;

enum MeshMode { Triangle = 0, Polygon, Points };
int MESH_MODE = MeshMode::Triangle;

// == Intrinsic triangulation stuff. All remeshing is performed on the manifold mesh.
VertexData<Vector3> csPositions;
std::unique_ptr<ManifoldSurfaceMesh> manifoldMesh;
std::unique_ptr<VertexPositionGeometry> manifoldGeom, csGeom;
std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation> intTri;
std::vector<Curve> CURVES, curvesOnManifold, curvesOnIntrinsic;
std::vector<SurfacePoint> POINTS, pointsOnManifold, pointsOnIntrinsic;
std::vector<std::vector<pointcloud::Point>> POINTS_CURVES;
std::vector<std::vector<Vertex>> POLYGON_CURVES;
float REFINE_AREA_THRESH = std::numeric_limits<float>::infinity();
float REFINE_ANGLE_THRESH = 25.;
int MAX_INSERTIONS = -1;
bool INTTRI_UPDATED = false;

enum SolverMode { ExtrinsicMesh = 0, IntrinsicMesh };
int SOLVER_MODE = SolverMode::ExtrinsicMesh;
int LAST_SOLVER_MODE;

// Solvers & parameters
float TCOEF = 1.;
bool TIME_UPDATED = false;
SignedHeatOptions SHM_OPTIONS;
int CONSTRAINT_MODE = static_cast<int>(LevelSetConstraint::ZeroSet);
bool SOLVE_AS_POINT_CLOUD = false;
bool EXPORT_RESULT = false;
bool VIZ = true;
bool VERBOSE, HEADLESS, IS_POLY;
std::unique_ptr<SignedHeatSolver> signedHeatSolver, intrinsicSolver;
std::unique_ptr<pointcloud::PointCloudHeatSolver> pointCloudSolver;
VertexData<double> PHI, PHI_CS;
pointcloud::PointData<double> PHI_POINTS;

// Program variables
std::string MESHNAME = "input mesh";
std::string MESH_FILEPATH, MESHROOT, OUTPUT_FILENAME, DATA_DIR;
std::string OUTPUT_DIR = "../export";
bool VIS_INTRINSIC_MESH = false;
bool COMMON_SUBDIVISION = true;
bool USE_BOUNDS = false;
float LOWER_BOUND, UPPER_BOUND;
bool CONSTRAINED_GEOM;
bool EXPORT_ON_CS = true;
std::chrono::time_point<high_resolution_clock> t1, t2;
std::chrono::duration<double, std::milli> ms_fp;


void ensureHaveIntrinsicSolver() {

    if (intrinsicSolver != nullptr) return;

    if (mesh->isManifold() && mesh->isOriented() & CONSTRAINED_GEOM) {
        if (VERBOSE) std::cerr << "Constructing intrinsic solver..." << std::endl;
        t1 = high_resolution_clock::now();
        setIntrinsicSolver(*geometry, CURVES, POINTS, manifoldMesh, manifoldGeom, curvesOnManifold, pointsOnManifold,
                           intTri, intrinsicSolver);
        t2 = high_resolution_clock::now();
        ms_fp = t2 - t1;
        if (VERBOSE) std::cerr << "Intrinsic solver construction time (s): " << ms_fp.count() / 1000. << std::endl;
    }
}

void solve() {

    if (MESH_MODE == MeshMode::Triangle) {
        if (SOLVER_MODE == SolverMode::ExtrinsicMesh) {

            SHM_OPTIONS.levelSetConstraint = static_cast<LevelSetConstraint>(CONSTRAINT_MODE);
            if (TIME_UPDATED) signedHeatSolver->setDiffusionTimeCoefficient(TCOEF);

            t1 = high_resolution_clock::now();
            PHI = signedHeatSolver->computeDistance(CURVES, POINTS, SHM_OPTIONS);
            t2 = high_resolution_clock::now();
            ms_fp = t2 - t1;
            if (VERBOSE) std::cerr << "Solve time (s): " << ms_fp.count() / 1000. << std::endl;

            if (!HEADLESS) {
                if (SHM_OPTIONS.levelSetConstraint != LevelSetConstraint::Multiple) {
                    psMesh->addVertexSignedDistanceQuantity("GSD", PHI)->setEnabled(true);
                } else {
                    // If there's multiple level sets, it's arbitrary which one should be "zero".
                    psMesh->addVertexScalarQuantity("GSD", PHI)->setIsolinesEnabled(true)->setEnabled(true);
                }
            }
            LAST_SOLVER_MODE = SolverMode::ExtrinsicMesh;

        } else {
            ensureHaveIntrinsicSolver();
            determineSourceGeometryOnIntrinsicTriangulation(*intTri, curvesOnManifold, pointsOnManifold,
                                                            curvesOnIntrinsic, pointsOnIntrinsic);
            if (INTTRI_UPDATED) {
                intrinsicSolver.reset(new SignedHeatSolver(*intTri));
                INTTRI_UPDATED = false;
            }
            SHM_OPTIONS.levelSetConstraint = static_cast<LevelSetConstraint>(CONSTRAINT_MODE);
            if (TIME_UPDATED) intrinsicSolver->setDiffusionTimeCoefficient(TCOEF);

            t1 = high_resolution_clock::now();
            PHI_CS = intrinsicSolver->computeDistance(curvesOnIntrinsic, pointsOnIntrinsic, SHM_OPTIONS);
            t2 = high_resolution_clock::now();
            ms_fp = t2 - t1;
            if (VERBOSE) std::cerr << "Solve time (s): " << ms_fp.count() / 1000. << std::endl;

            PHI = transferBtoA(*intTri, PHI_CS, TransferMethod::L2);
            if (!HEADLESS) {
                if (SHM_OPTIONS.levelSetConstraint != LevelSetConstraint::Multiple) {
                    psMesh->addVertexSignedDistanceQuantity("GSD", PHI)->setEnabled(true);
                    visualizeOnCommonSubdivision(*intTri, *manifoldGeom, csPositions, csGeom, PHI_CS, "GSD", true,
                                                 true);
                } else {
                    // If there's multiple level sets, it's arbitrary which one should be "zero".
                    psMesh->addVertexScalarQuantity("GSD", PHI)->setIsolinesEnabled(true)->setEnabled(true);
                    visualizeOnCommonSubdivision(*intTri, *manifoldGeom, csPositions, csGeom, PHI_CS, "GSD", true,
                                                 false);
                }
            }
            LAST_SOLVER_MODE = SolverMode::IntrinsicMesh;
        }
    } else if (MESH_MODE == MeshMode::Points) {
        // Reset point cloud solver in case parameters changed.
        if (TIME_UPDATED) pointCloudSolver.reset(new pointcloud::PointCloudHeatSolver(*cloud, *pointGeom, TCOEF));

        t1 = high_resolution_clock::now();
        PHI_POINTS = pointCloudSolver->computeSignedDistance(POINTS_CURVES, pointNormals, SHM_OPTIONS);
        t2 = high_resolution_clock::now();
        ms_fp = t2 - t1;
        if (VERBOSE) std::cerr << "Solve time (s): " << ms_fp.count() / 1000. << std::endl;

        if (!HEADLESS) psCloud->addScalarQuantity("GSD", PHI_POINTS)->setIsolinesEnabled(true)->setEnabled(true);

    } else if (MESH_MODE == MeshMode::Polygon) {
        throw std::logic_error("SHM on polygon meshes is not yet implemented - ETA mid-September 2024");
        // TODO: Set up polygon mesh solver
        // signedHeatSolver->setDiffusionTimeCoefficient(TCOEF);
        // Vector<double> phi = signedHeatSolver->solve(CURVES, POINTS, SHM_OPTIONS);

        // if (!HEADLESS) psMesh->addVertexSignedDistanceQuantity("GSD", phi)->setEnabled(true);
        // if (EXPORT_RESULT) {
        //     exportCurves(geometry->vertexPositions, CURVES, POINTS, OUTPUT_DIR + "/source.obj");
        //     exportSDF(*geometry, PHI, OUTPUT_FILENAME, USE_BOUNDS, LOWER_BOUND, UPPER_BOUND);
        // }
    }

    TIME_UPDATED = false;
}

void callback() {

    if (ImGui::Button("Solve")) {
        solve();
    }
    if (ImGui::Button("Export last solution")) {
        if (MESH_MODE == MeshMode::Triangle) {
            exportCurves(geometry->vertexPositions, CURVES, POINTS, OUTPUT_DIR);
            if (LAST_SOLVER_MODE == SolverMode::ExtrinsicMesh) {
                exportSDF(*geometry, PHI, OUTPUT_FILENAME, USE_BOUNDS, LOWER_BOUND, UPPER_BOUND);
            } else if (LAST_SOLVER_MODE == SolverMode::IntrinsicMesh) {
                ImGui::Checkbox("On common subdivision (vs. input mesh)", &EXPORT_ON_CS);
                if (!EXPORT_ON_CS) {
                    exportSDF(*geometry, PHI, OUTPUT_FILENAME, USE_BOUNDS, LOWER_BOUND, UPPER_BOUND);
                } else {
                    exportSDF(*intTri, *manifoldGeom, PHI_CS, OUTPUT_FILENAME, USE_BOUNDS, LOWER_BOUND, UPPER_BOUND);
                }
            }
        } else if (MESH_MODE == MeshMode::Points) {
            exportCurves(pointGeom->positions, POINTS_CURVES, OUTPUT_DIR + "/source.obj");
            exportSDF(pointGeom->positions, PHI_POINTS, OUTPUT_FILENAME);
        } else if (MESH_MODE == MeshMode::Polygon) {
            // TODO
        }
    }

    if (ImGui::TreeNode("Solve options")) {

        if (ImGui::InputFloat("tCoef", &TCOEF)) TIME_UPDATED = true;

        ImGui::Checkbox("Preserve source normals", &SHM_OPTIONS.preserveSourceNormals);
        ImGui::RadioButton("Constrain zero set", &CONSTRAINT_MODE, static_cast<int>(LevelSetConstraint::ZeroSet));
        ImGui::RadioButton("Constrain multiple levelsets", &CONSTRAINT_MODE,
                           static_cast<int>(LevelSetConstraint::Multiple));
        ImGui::RadioButton("No levelset constraints", &CONSTRAINT_MODE, static_cast<int>(LevelSetConstraint::None));

        ImGui::InputDouble("soft weight", &(SHM_OPTIONS.softLevelSetWeight));

        ImGui::Checkbox("Export result", &EXPORT_RESULT);
        ImGui::Checkbox("Specify upper/lower bounds for export", &USE_BOUNDS);
        if (ImGui::TreeNode("Bounds")) {
            ImGui::InputFloat("lower", &LOWER_BOUND);
            ImGui::InputFloat("upper", &UPPER_BOUND);
            ImGui::TreePop();
        }

        ImGui::TreePop();
    }

    if (MESH_MODE == MeshMode::Triangle && mesh->isManifold() && CONSTRAINED_GEOM) {
        ImGui::RadioButton("Solve on extrinsic mesh", &SOLVER_MODE, SolverMode::ExtrinsicMesh);
        ImGui::RadioButton("Solve on intrinsic mesh", &SOLVER_MODE, SolverMode::IntrinsicMesh);

        if (ImGui::TreeNode("Intrinsic mesh improvement")) {
            if (ImGui::Checkbox("Show intrinsic edges", &VIS_INTRINSIC_MESH)) {
                ensureHaveIntrinsicSolver();
                visualizeIntrinsicEdges(*intTri, *manifoldGeom, VIS_INTRINSIC_MESH);
                INTTRI_UPDATED = true;
            }

            if (ImGui::Button("Flip to Delaunay")) {
                ensureHaveIntrinsicSolver();
                intTri->flipToDelaunay();
                VIS_INTRINSIC_MESH = true;
                visualizeIntrinsicEdges(*intTri, *manifoldGeom, VIS_INTRINSIC_MESH);
                INTTRI_UPDATED = true;
            }

            ImGui::InputFloat("Angle thresh", &REFINE_ANGLE_THRESH);
            ImGui::InputFloat("Area thresh", &REFINE_AREA_THRESH);
            ImGui::InputInt("Max insert", &MAX_INSERTIONS);
            if (ImGui::Button("Delaunay refine")) {
                ensureHaveIntrinsicSolver();
                intTri->delaunayRefine(REFINE_ANGLE_THRESH, REFINE_AREA_THRESH, MAX_INSERTIONS);
                VIS_INTRINSIC_MESH = true;
                visualizeIntrinsicEdges(*intTri, *manifoldGeom, VIS_INTRINSIC_MESH);
                INTTRI_UPDATED = true;
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
    args::Flag points(group, "point", "Solve as point cloud.", {"p", "points"});
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
    MESHROOT = polyscope::guessNiceNameFromPath(MESH_FILEPATH);

    if (outputFilename) OUTPUT_DIR = getHomeDirectory(args::get(outputFilename));
    OUTPUT_FILENAME = OUTPUT_DIR + "/GSD.obj";
    HEADLESS = headless;
    VERBOSE = verbose;

    if (points) {
        MESH_MODE = MeshMode::Points;
        // Point cloud surface domain needs normals; it's up to the user to compute a consistent assignment of
        // normals. In the meantime just load in a mesh and inherit normals from the mesh. std::vector<Vector3>
        // positions = readPointCloud(MESH_FILEPATH);
        std::tie(mesh, geometry) = readSurfaceMesh(MESH_FILEPATH);
        size_t nPts = mesh->nVertices();
        cloud = std::unique_ptr<pointcloud::PointCloud>(new pointcloud::PointCloud(nPts));
        pointPositions = pointcloud::PointData<Vector3>(*cloud);
        pointNormals = pointcloud::PointData<Vector3>(*cloud);
        geometry->requireVertexNormals();
        for (size_t i = 0; i < nPts; i++) {
            pointPositions[i] = geometry->vertexPositions[i];
            pointNormals[i] = geometry->vertexNormals[i];
        }
        geometry->unrequireVertexNormals();
        pointGeom = std::unique_ptr<pointcloud::PointPositionGeometry>(
            new pointcloud::PointPositionGeometry(*cloud, pointPositions));
        pointCloudSolver.reset(new pointcloud::PointCloudHeatSolver(*cloud, *pointGeom, TCOEF));
    } else {
        std::tie(mesh, geometry) = readSurfaceMesh(MESH_FILEPATH);
        if (!mesh->isTriangular()) MESH_MODE = MeshMode::Polygon;
        signedHeatSolver = std::unique_ptr<SignedHeatSolver>(new SignedHeatSolver(*geometry));
    }

    // Load source geometry.
    if (inputFilename) {
        std::string filename = args::get(inputFilename);
        switch (MESH_MODE) {
            case (MeshMode::Triangle): {
                std::tie(CURVES, POINTS) = readInput(*mesh, filename);
                CONSTRAINED_GEOM = isSourceGeometryConstrained(CURVES, POINTS);
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
                psMesh = polyscope::registerSurfaceMesh(MESHNAME, geometry->vertexPositions, mesh->getFaceVertexList());
                psMesh->setAllPermutations(polyscopePermutations(*mesh));
                displayInput(geometry->vertexPositions, CURVES, POINTS);
                break;
            }
            case (MeshMode::Polygon): {
                psMesh = polyscope::registerSurfaceMesh(MESHNAME, geometry->vertexPositions, mesh->getFaceVertexList());
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