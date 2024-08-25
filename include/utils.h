#pragma once

#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/signed_heat_method.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// ================ GEOMETRIC

Vector3 centroid(VertexPositionGeometry& geometry);

double radius(VertexPositionGeometry& geometry, const Vector3& c);

std::tuple<Vector3, Vector3> boundingBox(VertexPositionGeometry& geometry);

// ================ TOPOLOGICAL

bool isSourceGeometryConstrained(const std::vector<Curve>& curves, const std::vector<Point>& points);

Halfedge determineHalfedgeFromVertices(const Vertex& vA, const Vertex& vB);

SurfacePoint reinterpretTo(const SurfacePoint& p, SurfaceMesh& otherMesh);

// ================ I/O

std::string getHomeDirectory(const std::string& filepath);

void readPointCloud(const std::string& filename, PointData<Vector3>& pointPositions, pointcloud::PointCloud& cloud,
                    PointPositionGeometry& pointGeom);

std::tuple<std::vector<Curve>, std::vector<Point>> readInput(SurfaceMesh& mesh, const std::string& filename);
std::vector<std::vector<Vertex>> readCurveVertices(SurfaceMesh& mesh, const std::string& filename);
std::vector<std::vector<Point>> readCurvePoints(pointcloud::PointCloud& cloud, const std::string& filename);

std::vector<std::vector<std::array<size_t, 2>>>
getCurveComponents(SurfaceMesh& mesh, const std::vector<SurfacePoint>& curveNodes,
                   const std::vector<std::array<size_t, 2>>& curveEdges);

/* `phi` is data on vertices */
std::vector<Curve> extractLevelsetAsCurves(IntrinsicGeometryInterface& geom, const Vector<double>& phi,
                                           double isoval = 0.);

/* Export curves as OBJ. */
void exportCurves(const VertexData<Vector3>& vertexPositions, const std::vector<Curve>& curves,
                  const std::vector<Point>& points, const std::string& filename = "../export/curves.obj");

/* Export SDF on an extrinsic surface mesh. */
void exportSDF(EmbeddedGeometryInterface& geom, const VertexData<double>& u, const std::string& filename,
               bool useBounds = false, double lowerBound = -1, double upperBound = -1);

/* Export SDF on a point cloud. */
void exportSDF(const PointData<Vector3>& pointPositions, const PointData<double>& u, const std::string& filename);

// ===================== MESH MUTATION

std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>
createIntrinsicTriangulation(VertexPositionGeometry& geometry, ManifoldSurfaceMesh& mesh,
                             const std::vector<Curve>& curves);

void setIntrinsicSolver(VertexPositionGeometry& geometry, const std::vector<Curve>& curves,
                        const std::vector<Point>& points, std::unique_ptr<ManifoldSurfaceMesh>& manifoldMesh,
                        std::unique_ptr<VertexPositionGeometry>& manifoldGeom, std::vector<Curve>& curvesOnManifold,
                        std::vector<Point>& pointsOnManifold,
                        std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>& intTri,
                        std::unique_ptr<SignedHeatMethodSolver>& solver);

void determineSourceGeometryOnIntrinsicTriangulation(IntrinsicTriangulation& intTri,
                                                     const std::vector<Curve>& curvesOnManifold,
                                                     const std::vector<Point>& pointsOnManifold,
                                                     std::vector<Curve>& curvesOnIntrinsic,
                                                     std::vector<Point>& pointsOnIntrinsic);

// ================ VISUALIZATION

void displayInput(const VertexData<Vector3>& vertexPositions, const std::vector<Curve>& curves,
                  const std::vector<Point>& pointSources, const std::string& name = "input", bool display = true);

void displayInput(const VertexData<Vector3>& vertexPositions, const std::vector<std::vector<Vertex>>& curves,
                  const std::string& name = "input", bool display = true);

void displayInput(const PointData<Vector3>& positions, const std::vector<std::vector<Point>>& curves,
                  const std::string& name = "input", bool display = true);