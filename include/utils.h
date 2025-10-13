#pragma once

#include "geometrycentral/pointcloud/point_position_geometry.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/poisson_disk_sampler.h"
#include "geometrycentral/surface/signed_heat_method.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

const Vector3 COLOR_NEGATIVE = {0.0627451, 0.517647, 0.94902}; // pasadena blue
const Vector3 COLOR_POSITIVE = {0.992157, 0.431373, 0.337255}; // grapefruit

int roundToNearestInt(double x);

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// ================ GEOMETRIC

Vector3 centroid(VertexPositionGeometry& geometry);

double radius(VertexPositionGeometry& geometry, const Vector3& c);

std::tuple<Vector3, Vector3> boundingBox(VertexPositionGeometry& geometry);

// ================ TOPOLOGICAL

bool isSourceGeometryConstrained(const std::vector<Curve>& curves, const std::vector<SurfacePoint>& points);

Halfedge determineHalfedgeFromVertices(const Vertex& vA, const Vertex& vB);

SurfacePoint reinterpretTo(const SurfacePoint& p, SurfaceMesh& otherMesh);

// ================ I/O

std::string getHomeDirectory(const std::string& filepath);

std::tuple<std::vector<Vector3>, std::vector<Vector3>> readPointCloud(const std::string& filename);

std::tuple<std::vector<Curve>, std::vector<SurfacePoint>> readInput(SurfaceMesh& mesh, const std::string& filename);
std::vector<std::vector<Vertex>> readCurveVertices(SurfaceMesh& mesh, const std::string& filename);
std::vector<std::vector<pointcloud::Point>> readCurvePoints(pointcloud::PointCloud& cloud, const std::string& filename);
std::vector<Curve> preprocessCurves(const std::vector<Curve>& curves);

std::vector<std::vector<std::array<size_t, 2>>>
getCurveComponents(SurfaceMesh& mesh, const std::vector<SurfacePoint>& curveNodes,
                   const std::vector<std::array<size_t, 2>>& curveEdges);

/* `phi` is data on vertices */
std::vector<Curve> extractLevelsetAsCurves(IntrinsicGeometryInterface& geom, const Vector<double>& phi,
                                           double isoval = 0.);

/* Export curves as OBJ. */
void exportCurves(const VertexData<Vector3>& vertexPositions, const std::vector<Curve>& curves,
                  const std::vector<SurfacePoint>& points, const std::string& dir = "../output");

void exportCurves(const VertexData<Vector3>& vertexPositions, const std::vector<std::vector<Vertex>>& curves,
                  const std::string& dir = "../output");

void exportCurves(const pointcloud::PointData<Vector3>& positions,
                  const std::vector<std::vector<pointcloud::Point>>& curves, const std::string& dir = "../output");

/* Export SDF on an extrinsic surface mesh. */
void exportSDF(EmbeddedGeometryInterface& geom, const VertexData<double>& u, const std::string& filename,
               bool useBounds = false, double lowerBound = -1, double upperBound = -1);

void exportSDF(EmbeddedGeometryInterface& geom, const CornerData<double>& u, const std::string& filename,
               bool useBounds = false, double lowerBound = -1, double upperBound = -1);

/* Export SDF on common subdivision. */
void exportSDF(IntegerCoordinatesIntrinsicTriangulation& intTri, VertexPositionGeometry& manifoldGeom,
               const VertexData<double>& u, const std::string& filename, bool useBounds = false, double lowerBound = -1,
               double upperBound = -1);

void exportSDF(IntegerCoordinatesIntrinsicTriangulation& intTri, VertexPositionGeometry& manifoldGeom,
               const CornerData<double>& u, const std::string& filename, bool useBounds = false, double lowerBound = -1,
               double upperBound = -1, bool normalize = true);

/* Export SDF on a point cloud. */
void exportSDF(const pointcloud::PointData<Vector3>& pointPositions, const pointcloud::PointData<double>& u,
               const std::string& filename);

/* Normalize SDF data so it maps onto (custom) divergent colormap correctly. */
Vector<double> normalizeSDF(const Vector<double>& u, const std::string& dir = "", bool useBounds = false,
                            double lowerBound = -1, double upperBound = -1);

void exportVectors(VertexPositionGeometry& geometry, const VertexData<Vector2>& vectorField,
                   const std::string& filename, double sampleSparsity = -1);

void exportVectors(VertexPositionGeometry& geometry, const VertexData<Vector3>& vectorField,
                   const std::string& filename, double sampleSparsity = -1);

// ===================== MESH MUTATION

std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>
createIntrinsicTriangulation(VertexPositionGeometry& geometry, ManifoldSurfaceMesh& mesh,
                             const std::vector<Curve>& curves);

void setIntrinsicSolver(VertexPositionGeometry& geometry, const std::vector<Curve>& curves,
                        const std::vector<SurfacePoint>& points, std::unique_ptr<ManifoldSurfaceMesh>& manifoldMesh,
                        std::unique_ptr<VertexPositionGeometry>& manifoldGeom, std::vector<Curve>& curvesOnManifold,
                        std::vector<SurfacePoint>& pointsOnManifold,
                        std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>& intTri,
                        std::unique_ptr<SignedHeatSolver>& solver);

void determineSourceGeometryOnIntrinsicTriangulation(IntrinsicTriangulation& intTri,
                                                     const std::vector<Curve>& curvesOnManifold,
                                                     const std::vector<SurfacePoint>& pointsOnManifold,
                                                     std::vector<Curve>& curvesOnIntrinsic,
                                                     std::vector<SurfacePoint>& pointsOnIntrinsic);

// ================ VISUALIZATION

void displayInput(const VertexData<Vector3>& vertexPositions, const std::vector<Curve>& curves,
                  const std::vector<SurfacePoint>& pointSources, const std::string& name = "input",
                  bool display = true);

void displayInput(const VertexData<Vector3>& vertexPositions, const std::vector<std::vector<Vertex>>& curves,
                  const std::string& name = "input", bool display = true);

void displayInput(const pointcloud::PointData<Vector3>& positions,
                  const std::vector<std::vector<pointcloud::Point>>& curves, const std::string& name = "input",
                  bool display = true);

void visualizeIntrinsicEdges(IntegerCoordinatesIntrinsicTriangulation& intTri, VertexPositionGeometry& manifoldGeom,
                             bool display = true);

VertexData<double> visualizeOnCommonSubdivision(IntegerCoordinatesIntrinsicTriangulation& intTri,
                                                VertexPositionGeometry& manifoldGeom, VertexData<Vector3>& csPositions,
                                                std::unique_ptr<VertexPositionGeometry>& csGeom,
                                                const VertexData<double>& phi, const std::string& quantityName,
                                                bool rebuild = true, bool divergingColormap = true);