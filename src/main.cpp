#include "geometrycentral/pointcloud/point_cloud_heat_solver.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/polygon_mesh_heat_solver.h"
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
bool VIZ = true;
bool VERBOSE, HEADLESS;
std::unique_ptr<SignedHeatSolver> signedHeatSolver, intrinsicSolver;
std::unique_ptr<pointcloud::PointCloudHeatSolver> pointCloudSolver;
std::unique_ptr<PolygonMeshHeatSolver> polygonSolver;
VertexData<double> PHI, PHI_CS;
pointcloud::PointData<double> PHI_POINTS;

// Program variables
std::string MESHNAME = "input mesh";
std::string MESH_FILEPATH, OUTPUT_FILENAME, DATA_DIR;
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
            if (VERBOSE) {
                std::cerr << "min: " << PHI.toVector().minCoeff() << "\tmax: " << PHI.toVector().maxCoeff()
                          << std::endl;
                std::cerr << "Solve time (s): " << ms_fp.count() / 1000. << std::endl;
            }

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
            if (VERBOSE) {
                std::cerr << "min: " << PHI_CS.toVector().minCoeff() << "\tmax: " << PHI_CS.toVector().maxCoeff()
                          << std::endl;
                std::cerr << "Solve time (s): " << ms_fp.count() / 1000. << std::endl;
            }

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
        SHM_OPTIONS.levelSetConstraint = static_cast<LevelSetConstraint>(CONSTRAINT_MODE);
        // Reset point cloud solver in case parameters changed.
        if (TIME_UPDATED) pointCloudSolver.reset(new pointcloud::PointCloudHeatSolver(*cloud, *pointGeom, TCOEF));

        t1 = high_resolution_clock::now();
        PHI_POINTS = pointCloudSolver->computeSignedDistance(POINTS_CURVES, pointNormals, SHM_OPTIONS);
        t2 = high_resolution_clock::now();
        ms_fp = t2 - t1;
        if (VERBOSE) {
            std::cerr << "min: " << PHI_POINTS.toVector().minCoeff() << "\tmax: " << PHI_POINTS.toVector().maxCoeff()
                      << std::endl;
            std::cerr << "Solve time (s): " << ms_fp.count() / 1000. << std::endl;
        }

        if (!HEADLESS) psCloud->addScalarQuantity("GSD", PHI_POINTS)->setIsolinesEnabled(true)->setEnabled(true);

    } else if (MESH_MODE == MeshMode::Polygon) {

        if (TIME_UPDATED) polygonSolver.reset(new PolygonMeshHeatSolver(*geometry, TCOEF));

        t1 = high_resolution_clock::now();
        PHI = polygonSolver->computeSignedDistance(POLYGON_CURVES, static_cast<LevelSetConstraint>(CONSTRAINT_MODE));
        t2 = high_resolution_clock::now();
        ms_fp = t2 - t1;
        if (VERBOSE) {
            std::cerr << "min: " << PHI.toVector().minCoeff() << "\tmax: " << PHI.toVector().maxCoeff() << std::endl;
            std::cerr << "Solve time (s): " << ms_fp.count() / 1000. << std::endl;
        }
        if (!HEADLESS) psMesh->addVertexSignedDistanceQuantity("GSD", PHI)->setEnabled(true);
    }

    TIME_UPDATED = false;
}

Vector<double> unsignedHeatMethodStepsIIandIII(IntrinsicGeometryInterface& geom, const Vector<double>& u) {

    SurfaceMesh& mesh = geom.mesh;
    size_t V = mesh.nVertices();
    Vector<double> div = Vector<double>::Zero(V);
    for (Face f : mesh.faces()) {
        // taken from geometry-central implementation
        Vector2 gradUDir = Vector2::zero();
        for (Halfedge he : f.adjacentHalfedges()) {
            Vector2 ePerp = geom.halfedgeVectorsInFace[he.next()].rotate90();
            gradUDir += ePerp * u(geom.vertexIndices[he.vertex()]);
        }
        gradUDir = gradUDir.normalizeCutoff();

        for (Halfedge he : f.adjacentHalfedges()) {
            double val = geom.halfedgeCotanWeights[he] * dot(geom.halfedgeVectorsInFace[he], gradUDir);
            div[geom.vertexIndices[he.tailVertex()]] += val;
            div[geom.vertexIndices[he.tipVertex()]] += -val;
        }
    }
    SparseMatrix<double> C = geom.cotanLaplacian;
    Vector<double> phi = solvePositiveDefinite(C, div);
    double shift = phi.minCoeff();
    phi -= shift * Vector<double>::Ones(V);
    return phi;
}

VertexData<double> unsignedHeatMethod(IntrinsicGeometryInterface& geom, const std::vector<Curve>& curves,
                                      const std::vector<SurfacePoint>& points) {

    SurfaceMesh& mesh = geom.mesh;

    geom.requireVertexLumpedMassMatrix();
    geom.requireCotanLaplacian();
    geom.requireEdgeLengths();
    geom.requireVertexIndices();
    geom.requireHalfedgeCotanWeights();
    geom.requireHalfedgeVectorsInFace();

    size_t V = mesh.nVertices();

    // Compute mean node spacing.
    double h = 0;
    for (Edge e : mesh.edges()) {
        h += geom.edgeLengths[e];
    }
    h /= mesh.nEdges();
    double t = h * h;

    // Get operators
    SparseMatrix<double> C = geom.cotanLaplacian;
    SparseMatrix<double> M = geom.vertexLumpedMassMatrix;
    SparseMatrix<double> heatOp = M + t * C;

    PositiveDefiniteSolver<double> heatSolver(heatOp);

    Vector<double> RHS = Vector<double>::Zero(V);
    // Project curves, which are always assumed to be unsigned
    for (const Curve& curve : curves) {
        size_t nNodes = curve.nodes.size();
        for (size_t i = 0; i < nNodes - 1; i++) {
            const SurfacePoint& pA = curve.nodes[i];
            const SurfacePoint& pB = curve.nodes[i + 1];
            Edge commonEdge = sharedEdge(pA, pB);
            if (commonEdge != Edge()) {
                double length = geom.edgeLengths[commonEdge];
                Vertex vA = commonEdge.firstVertex();
                Vertex vB = commonEdge.secondVertex();
                size_t idxA = geom.vertexIndices[vA];
                size_t idxB = geom.vertexIndices[vB];
                RHS[idxA] += 0.5 * length;
                RHS[idxB] += 0.5 * length;
                continue;
            }
            Face commonFace = sharedFace(pA, pB);
            SurfacePoint pAInFace = pA.inFace(commonFace);
            SurfacePoint pBInFace = pB.inFace(commonFace);
            Vector3 mCoords = 0.5 * (pAInFace.faceCoords + pBInFace.faceCoords);
            BarycentricVector vec(pAInFace, pBInFace);
            double length = vec.norm(geom);
            size_t rIdx = 0;
            for (Vertex v : commonFace.adjacentVertices()) {
                size_t vIdx = geom.vertexIndices[v];
                RHS[vIdx] += length * mCoords[rIdx];
                rIdx++;
            }
        }
    }

    // Project points (assumed to be always at vertices, and always unsigned)
    for (const SurfacePoint& p : points) {
        if (p.type != SurfacePointType::Vertex) continue;
        size_t vIdx = geom.vertexIndices[p.vertex];
        RHS[vIdx] += 1;
        // TODO: handle barycentric points (will probably never need to compare against this though)
    }
    Vector<double> u_neumann = heatSolver.solve(RHS);

    Vector<bool> setAMembership = Vector<bool>::Ones(V); // true if interior
    for (Vertex v : mesh.vertices()) {
        if (v.isBoundary()) {
            setAMembership[geom.vertexIndices[v]] = false;
        }
    }
    int nB = V - setAMembership.cast<int>().sum(); // # of boundary vertices
    Vector<double> bcVals = Vector<double>::Zero(nB);
    BlockDecompositionResult<double> decomp = blockDecomposeSquare(heatOp, setAMembership, true);
    Vector<double> rhsValsA, rhsValsB;
    decomposeVector(decomp, RHS, rhsValsA, rhsValsB);
    Vector<double> combinedRHS = rhsValsA;
    shiftDiagonal(decomp.AA, 1e-8);
    Vector<double> Aresult = solvePositiveDefinite(decomp.AA, combinedRHS);
    Vector<double> u_dirichlet = reassembleVector(decomp, Aresult, bcVals);

    Vector<double> u = 0.5 * (u_neumann + u_dirichlet);
    Vector<double> phi = unsignedHeatMethodStepsIIandIII(geom, u);
    return VertexData<double>(mesh, phi);
}

void callback() {

    if (ImGui::Button("Solve (unsigned heat method)")) {
        if (SOLVER_MODE == SolverMode::ExtrinsicMesh) {
            PHI = unsignedHeatMethod(*geometry, CURVES, POINTS);
            psMesh->addVertexSignedDistanceQuantity("UDF", PHI)->setEnabled(true);
            LAST_SOLVER_MODE = SolverMode::ExtrinsicMesh;
        } else {
            ensureHaveIntrinsicSolver();
            determineSourceGeometryOnIntrinsicTriangulation(*intTri, curvesOnManifold, pointsOnManifold,
                                                            curvesOnIntrinsic, pointsOnIntrinsic);
            PHI_CS = unsignedHeatMethod(*intTri, curvesOnIntrinsic, pointsOnIntrinsic);
            PHI = transferBtoA(*intTri, PHI_CS, TransferMethod::L2);
            psMesh->addVertexSignedDistanceQuantity("UDF", PHI)->setEnabled(true);
            visualizeOnCommonSubdivision(*intTri, *manifoldGeom, csPositions, csGeom, PHI_CS, "UDF", true, true);
            LAST_SOLVER_MODE = SolverMode::IntrinsicMesh;
        }
    }

    if (ImGui::Button("Solve")) {
        solve();
    }
    if (ImGui::Button("Export last solution")) {
        if (MESH_MODE == MeshMode::Triangle) {
            exportCurves(geometry->vertexPositions, CURVES, POINTS, OUTPUT_DIR);
            if (LAST_SOLVER_MODE == SolverMode::ExtrinsicMesh) {
                exportSDF(*geometry, PHI, OUTPUT_FILENAME, USE_BOUNDS, LOWER_BOUND, UPPER_BOUND);
            } else if (LAST_SOLVER_MODE == SolverMode::IntrinsicMesh) {
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
            exportCurves(geometry->vertexPositions, POLYGON_CURVES, OUTPUT_DIR);
            exportSDF(*geometry, PHI, OUTPUT_FILENAME, USE_BOUNDS, LOWER_BOUND, UPPER_BOUND);
        }
    }
    if (LAST_SOLVER_MODE == SolverMode::IntrinsicMesh)
        ImGui::Checkbox("On common subdivision (vs. input mesh)", &EXPORT_ON_CS);

    ImGui::Text("Solve options");
    ImGui::Separator();

    if (ImGui::InputFloat("tCoef", &TCOEF)) TIME_UPDATED = true;

    ImGui::Checkbox("Preserve source normals", &SHM_OPTIONS.preserveSourceNormals);
    ImGui::RadioButton("Constrain zero set", &CONSTRAINT_MODE, static_cast<int>(LevelSetConstraint::ZeroSet));
    ImGui::RadioButton("Constrain multiple levelsets", &CONSTRAINT_MODE,
                       static_cast<int>(LevelSetConstraint::Multiple));
    ImGui::RadioButton("No levelset constraints", &CONSTRAINT_MODE, static_cast<int>(LevelSetConstraint::None));

    ImGui::InputDouble("soft weight", &(SHM_OPTIONS.softLevelSetWeight));

    ImGui::Checkbox("Specify upper/lower bounds for export", &USE_BOUNDS);
    if (ImGui::TreeNode("Bounds")) {
        ImGui::InputFloat("lower", &LOWER_BOUND);
        ImGui::InputFloat("upper", &UPPER_BOUND);
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
    std::string ext = MESH_FILEPATH.substr(MESH_FILEPATH.find_last_of(".") + 1);
    MESH_MODE = (ext == "pc") ? MeshMode::Points : MeshMode::Triangle;

    if (outputFilename) OUTPUT_DIR = getHomeDirectory(args::get(outputFilename));
    OUTPUT_FILENAME = OUTPUT_DIR + "/GSD.obj";
    HEADLESS = headless;
    VERBOSE = verbose;

    if (MESH_MODE == MeshMode::Points) {
        std::vector<Vector3> positions, normals;
        std::tie(positions, normals) = readPointCloud(MESH_FILEPATH);
        size_t nPts = positions.size();
        cloud = std::unique_ptr<pointcloud::PointCloud>(new pointcloud::PointCloud(nPts));
        pointPositions = pointcloud::PointData<Vector3>(*cloud);
        pointNormals = pointcloud::PointData<Vector3>(*cloud);
        for (size_t i = 0; i < nPts; i++) {
            pointPositions[i] = positions[i];
            pointNormals[i] = normals[i];
        }
        pointGeom = std::unique_ptr<pointcloud::PointPositionGeometry>(
            new pointcloud::PointPositionGeometry(*cloud, pointPositions));
        pointCloudSolver.reset(new pointcloud::PointCloudHeatSolver(*cloud, *pointGeom, TCOEF));
    } else {
        std::tie(mesh, geometry) = readSurfaceMesh(MESH_FILEPATH);

        // // Center and scale.
        // Vector3 bboxMin, bboxMax;
        // std::tie(bboxMin, bboxMax) = boundingBox(*geometry);
        // double diag = (bboxMin - bboxMax).norm();
        // Vector3 center = centroid(*geometry);
        // for (Vertex v : mesh->vertices()) {
        //     Vector3 p = geometry->vertexPositions[v];
        //     geometry->vertexPositions[v] = (p - center) / diag;
        // }
        // geometry->refreshQuantities();

        if (!mesh->isTriangular()) {
            MESH_MODE = MeshMode::Polygon;
            polygonSolver = std::unique_ptr<PolygonMeshHeatSolver>(new PolygonMeshHeatSolver(*geometry));
        } else {
            signedHeatSolver = std::unique_ptr<SignedHeatSolver>(new SignedHeatSolver(*geometry));
        }
    }

    // Load source geometry.
    if (inputFilename) {
        std::string filename = args::get(inputFilename);
        switch (MESH_MODE) {
            case (MeshMode::Triangle): {
                std::tie(CURVES, POINTS) = readInput(*mesh, filename);
                CURVES = preprocessCurves(CURVES);
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