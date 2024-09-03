#include "utils.h"

int roundToNearestInt(double x) {

    double y = x + 0.5 - (x < 0);
    return (int)y;
}

// ================ GEOMETRIC

Vector3 centroid(VertexPositionGeometry& geometry) {

    Vector3 c = {0, 0, 0};
    SurfaceMesh& mesh = geometry.mesh;
    for (Vertex v : mesh.vertices()) {
        c += geometry.vertexPositions[v];
    }
    c /= mesh.nVertices();
    return c;
}

double radius(VertexPositionGeometry& geometry, const Vector3& c) {

    double r = 0;
    SurfaceMesh& mesh = geometry.mesh;
    for (Vertex v : mesh.vertices()) {
        r = std::max(r, (c - geometry.vertexPositions[v]).norm());
    }
    return r;
}

/*
 * Return one of the two pairs of furthest-away corners in the axis-aligned  bounding box of the mesh.
 *
 * I.e., the distance between bboxMin and bboxMax is the length of the bbox's  diagonal, and the length of the bbox in
 * the i-th dimension is |bboxMin[i] -  bboxMax[i]|.
 */
std::tuple<Vector3, Vector3> boundingBox(VertexPositionGeometry& geometry) {

    SurfaceMesh& mesh = geometry.mesh;
    const double inf = std::numeric_limits<double>::infinity();
    Vector3 bboxMin = {inf, inf, inf};
    Vector3 bboxMax = {-inf, -inf, -inf};
    for (Vertex v : mesh.vertices()) {
        Vector3& pos = geometry.vertexPositions[v];
        for (int i = 0; i < 3; i++) {
            if (pos[i] <= bboxMin[i]) bboxMin[i] = pos[i];
            if (pos[i] >= bboxMax[i]) bboxMax[i] = pos[i];
        }
    }
    return std::make_tuple(bboxMin, bboxMax);
}

// ================ TOPOLOGICAL

bool isSourceGeometryConstrained(const std::vector<Curve>& curves, const std::vector<SurfacePoint>& points) {

    for (const Curve& curve : curves) {
        for (const SurfacePoint& p : curve.nodes) {
            if (p.type != SurfacePointType::Vertex) return false;
        }
    }
    for (const SurfacePoint& p : points) {
        if (p.type != SurfacePointType::Vertex) return false;
    }
    return true;
}

Halfedge determineHalfedgeFromVertices(const Vertex& vA, const Vertex& vB) {

    for (Halfedge he : vA.outgoingHalfedges()) {
        if (he.tipVertex() == vB) return he;
    }
    return Halfedge();
}

/*
 * A SurfacePoint on the original input mesh, re-interpret to a new mesh with the same topology.
 */
SurfacePoint reinterpretTo(const SurfacePoint& p, SurfaceMesh& otherMesh) {

    switch (p.type) {
        case (SurfacePointType::Vertex):
            return SurfacePoint(otherMesh.vertex(p.vertex.getIndex()));
            break;
        case (SurfacePointType::Edge): {
            Edge e = p.edge;
            Vertex vA = otherMesh.vertex(e.firstVertex().getIndex());
            Vertex vB = otherMesh.vertex(e.secondVertex().getIndex());
            Edge otherEdge = determineHalfedgeFromVertices(vA, vB).edge();
            double tEdge = (e.firstVertex().getIndex() == otherEdge.firstVertex().getIndex()) ? p.tEdge : 1. - p.tEdge;
            return SurfacePoint(otherEdge, tEdge);
            break;
        }
        case (SurfacePointType::Face):
            return SurfacePoint(otherMesh.face(p.face.getIndex()), p.faceCoords);
            break;
    }
    throw std::logic_error("Bad switch");
}

// ================ I/O

std::string getHomeDirectory(const std::string& filepath) {

    std::string dir(filepath.begin(), filepath.begin() + filepath.find_last_of("/") + 1);
    return dir;
}

std::tuple<std::vector<Vector3>, std::vector<Vector3>> readPointCloud(const std::string& filepath) {

    std::ifstream curr_file(filepath.c_str());
    std::string line;
    std::string X;
    double x, y, z;
    std::vector<Vector3> positions, normals;
    if (curr_file.is_open()) {
        while (!curr_file.eof()) {
            getline(curr_file, line);
            // Ignore any newlines
            if (line == "") {
                continue;
            }
            std::istringstream iss(line);
            iss >> X;
            if (X == "v") {
                iss >> x >> y >> z;
                positions.push_back({x, y, z});
            } else if (X == "vn") {
                iss >> x >> y >> z;
                normals.push_back({x, y, z});
            }
        }
        curr_file.close();
    } else {
        std::cerr << "Could not open file <" << filepath << ">." << std::endl;
    }
    return std::make_tuple(positions, normals);
}

std::tuple<std::vector<Curve>, std::vector<SurfacePoint>> readInput(SurfaceMesh& mesh, const std::string& filename) {

    std::vector<Curve> curves;
    std::vector<SurfacePoint> pointSources;

    std::ifstream curr_file(filename.c_str());
    std::string line;
    std::string X;
    SurfacePoint pt;
    size_t idx;
    bool flip = false; // to flip orientation or not
    bool curveMode = true;
    double tEdge, a, b, c;
    std::vector<bool> curveFlips;

    if (curr_file.is_open()) {
        while (!curr_file.eof()) {
            getline(curr_file, line);
            // Ignore any newlines
            if (line == "") {
                continue;
            }
            std::istringstream iss(line);
            iss >> X;
            if (X == "v") {
                iss >> idx;
                pt = SurfacePoint(mesh.vertex(idx));
                if (curveMode) {
                    curves.back().nodes.push_back(pt);
                } else {
                    pointSources.push_back(pt);
                }
            } else if (X == "e") {
                iss >> idx >> tEdge;
                pt = SurfacePoint(mesh.edge(idx), tEdge);
                if (curveMode) {
                    curves.back().nodes.push_back(pt);
                } else {
                    pointSources.push_back(pt);
                }
            } else if (X == "f") {
                iss >> idx >> a >> b >> c;
                Vector3 faceCoords = {a, b, c};
                pt = SurfacePoint(mesh.face(idx), faceCoords);
                if (curveMode) {
                    curves.back().nodes.push_back(pt);
                } else {
                    pointSources.push_back(pt);
                }
            } else if (X == "l") {
                while (true) {
                    if (iss.eof()) break;
                    iss >> idx;
                    idx -= 1; // OBJ elements are 1-indexed, whereas geometry-central uses 0-indexing
                    curves.back().nodes.push_back(SurfacePoint(mesh.vertex(idx)));
                }
            } else if (X == "signed_curve") {
                iss >> flip;
                curveFlips.push_back(flip);
                curveMode = true;
                curves.emplace_back();
                curves.back().isSigned = true;
            } else if (X == "unsigned_curve") {
                iss >> flip;
                curveFlips.push_back(flip);
                curveMode = true;
                curves.emplace_back();
                curves.back().isSigned = false;
            } else if (X == "unsigned_point") {
                curveMode = false;
            } else if (X == "signed_point") {
                curveMode = false;
            }
        }
        curr_file.close();
        for (size_t i = 0; i < curves.size(); i++) {
            if (curveFlips[i]) {
                std::reverse(curves[i].nodes.begin(), curves[i].nodes.end());
            }
        }
    } else {
        std::cerr << "Could not open file <" << filename << ">." << std::endl;
    }
    return std::make_tuple(curves, pointSources);
}

std::vector<std::vector<Vertex>> readCurveVertices(SurfaceMesh& mesh, const std::string& filename) {

    std::vector<std::vector<Vertex>> curves;
    std::ifstream curr_file(filename.c_str());
    std::string line;
    std::string X;
    size_t idx;
    bool flip = false;
    bool read = false;
    std::vector<bool> curveFlips;

    if (curr_file.is_open()) {
        while (!curr_file.eof()) {
            getline(curr_file, line);
            // Ignore any newlines
            if (line == "") {
                continue;
            }
            std::istringstream iss(line);
            iss >> X;
            if (X == "signed_curve") {
                iss >> flip;
                curveFlips.push_back(flip);
                curves.emplace_back();
                read = true;
            } else if (X == "unsigned_curve") {
                read = false;
            } else if (X == "unsigned_point") {
                read = false;
            }
            if (read) {
                if (X == "v") {
                    iss >> idx;
                    curves.back().emplace_back(mesh.vertex(idx));
                } else if (X == "l") {
                    while (true) {
                        if (iss.eof()) break;
                        iss >> idx;
                        idx -= 1; // OBJ elements are 1-indexed, whereas geometry-central uses 0-indexing
                        curves.back().emplace_back(mesh.vertex(idx));
                    }
                }
            }
        }
        curr_file.close();

        for (size_t i = 0; i < curves.size(); i++) {
            if (curveFlips[i]) {
                std::reverse(curves[i].begin(), curves[i].end());
            }
        }

    } else {
        std::cerr << "Could not open file <" << filename << ">." << std::endl;
    }
    return curves;
}

std::vector<std::vector<pointcloud::Point>> readCurvePoints(pointcloud::PointCloud& cloud,
                                                            const std::string& filename) {

    std::vector<std::vector<pointcloud::Point>> curves;
    std::ifstream curr_file(filename.c_str());
    std::string line;
    std::string X;
    size_t idx;
    bool flip = false;
    bool read = false;
    std::vector<bool> curveFlips;

    if (curr_file.is_open()) {
        while (!curr_file.eof()) {
            getline(curr_file, line);
            // Ignore any newlines
            if (line == "") {
                continue;
            }
            std::istringstream iss(line);
            iss >> X;
            if (X == "signed_curve") {
                iss >> flip;
                curveFlips.push_back(flip);
                curves.emplace_back();
                read = true;
            } else if (X == "unsigned_curve") {
                read = false;
            } else if (X == "unsigned_point") {
                read = false;
            }
            if (read) {
                if (X == "v") {
                    iss >> idx;
                    curves.back().emplace_back(cloud.point(idx));
                } else if (X == "l") {
                    while (true) {
                        if (iss.eof()) break;
                        iss >> idx;
                        idx -= 1; // OBJ elements are 1-indexed, whereas geometry-central uses 0-indexing
                        curves.back().emplace_back(cloud.point(idx));
                    }
                }
            }
        }
        curr_file.close();

        for (size_t i = 0; i < curves.size(); i++) {
            if (curveFlips[i]) {
                std::reverse(curves[i].begin(), curves[i].end());
            }
        }

    } else {
        std::cerr << "Could not open file <" << filename << ">." << std::endl;
    }
    return curves;
}

std::vector<std::vector<std::array<size_t, 2>>>
getCurveComponents(SurfaceMesh& mesh, const std::vector<SurfacePoint>& curveNodes,
                   const std::vector<std::array<size_t, 2>>& curveEdges) {

    std::vector<std::array<size_t, 2>> edgesToAdd = curveEdges;
    std::vector<std::vector<std::array<size_t, 2>>> curves;
    size_t nSegs = curveEdges.size();
    while (edgesToAdd.size() > 0) {
        std::array<size_t, 2> startSeg = edgesToAdd.back();
        edgesToAdd.pop_back();
        curves.emplace_back();
        curves.back().push_back(startSeg);

        // Add segs to the front end until we can't.
        std::array<size_t, 2> currSeg = startSeg;
        while (true) {
            const SurfacePoint& front = curveNodes[currSeg[1]];
            bool didWeFindOne = false;
            for (size_t i = 0; i < edgesToAdd.size(); i++) {
                std::array<size_t, 2> otherSeg = edgesToAdd[i];
                if (curveNodes[otherSeg[0]] == front) {
                    currSeg = otherSeg;
                    curves.back().push_back(otherSeg);
                    edgesToAdd.erase(edgesToAdd.begin() + i);
                    didWeFindOne = true;
                    break;
                }
            }
            if (!didWeFindOne) break;
        }

        // Add segs to the back end until we can't.
        currSeg = startSeg;
        while (true) {
            const SurfacePoint& back = curveNodes[currSeg[0]];
            bool didWeFindOne = false;
            for (size_t i = 0; i < edgesToAdd.size(); i++) {
                std::array<size_t, 2> otherSeg = edgesToAdd[i];
                if (curveNodes[otherSeg[1]] == back) {
                    currSeg = otherSeg;
                    curves.back().insert(curves.back().begin(), otherSeg);
                    edgesToAdd.erase(edgesToAdd.begin() + i);
                    didWeFindOne = true;
                    break;
                }
            }
            if (!didWeFindOne) break;
        }
    }
    return curves;
}

std::vector<Curve> extractLevelsetAsCurves(IntrinsicGeometryInterface& geom, const Vector<double>& phi, double isoval) {

    SurfaceMesh& mesh = geom.mesh;
    geom.requireVertexIndices();
    std::vector<SurfacePoint> curveNodes;
    std::vector<std::array<size_t, 2>> curveEdges;

    // Find levelset.
    for (Face f : mesh.faces()) {
        std::vector<SurfacePoint> crossings;
        std::vector<std::array<double, 2>> vals;
        for (Halfedge he : f.adjacentHalfedges()) {
            size_t vA = geom.vertexIndices[he.tailVertex()];
            size_t vB = geom.vertexIndices[he.tipVertex()];
            double valA = phi[vA];
            double valB = phi[vB];
            if (!(valA <= isoval && isoval <= valB) && !(valB <= isoval && isoval <= valA)) continue;
            double t = (isoval - valA) / (valB - valA);
            vals.push_back({valA, valB});
            Edge e = he.edge();
            if (!he.orientation()) t = 1. - t;
            if (t >= 0. && t <= 1.) crossings.push_back(SurfacePoint(e, t));
        }
        if (crossings.size() != 2) continue;
        std::array<size_t, 2> seg;
        for (int i = 0; i < 2; i++) {
            auto iter = std::find(curveNodes.begin(), curveNodes.end(), crossings[i]);
            if (iter != curveNodes.end()) {
                seg[i] = iter - curveNodes.begin();
            } else {
                curveNodes.push_back(crossings[i]);
                seg[i] = curveNodes.size() - 1;
            }
        }
        // Switch order of segment, so that smaller values are on the "inside" of the curve (assuming CCW orientation.)
        if (vals[0][0] < vals[0][1]) {
            size_t temp = seg[0];
            seg[0] = seg[1];
            seg[1] = temp;
        }
        curveEdges.push_back(seg);
    }

    std::vector<std::vector<std::array<size_t, 2>>> curveComponents = getCurveComponents(mesh, curveNodes, curveEdges);
    std::vector<Curve> curves;
    for (size_t i = 0; i < curveComponents.size(); i++) {
        curves.emplace_back();
        Curve& curve = curves.back();
        curve.isSigned = true;
        size_t nSegs = curveComponents[i].size();
        for (const auto& seg : curveComponents[i]) {
            curve.nodes.push_back(curveNodes[seg[0]]);
        }
        curve.nodes.push_back(curveNodes[curveComponents[i][nSegs - 1][1]]);
    }
    geom.unrequireVertexIndices();

    return curves;
}

void exportCurves(const VertexData<Vector3>& vertexPositions, const std::vector<Curve>& curves,
                  const std::vector<SurfacePoint>& points, const std::string& dir) {

    std::string pFilename = dir + "/points.obj";
    std::string cFilename = dir + "/curves.obj";

    // Write points.
    std::fstream f;
    f.open(pFilename, std::ios::out | std::ios::trunc);
    if (f.is_open()) {
        for (const SurfacePoint& p : points) {
            Vector3 pos = p.interpolate(vertexPositions);
            f << "v " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
        }
        f.close();
        std::cerr << "File " << pFilename << " written succesfully." << std::endl;
    } else {
        std::cerr << "Could not save points '" << pFilename << "'!" << std::endl;
    }

    f.open(cFilename, std::ios::out | std::ios::trunc);
    if (f.is_open()) {
        // Assume curves have already been organized into components to be "as connected as possible".
        size_t offset = 0;
        for (const auto& curve : curves) {
            size_t nNodes = curve.nodes.size();
            bool isClosed = (curve.nodes[0] == curve.nodes[nNodes - 1]);
            size_t ub = isClosed ? nNodes - 1 : nNodes;
            for (size_t i = 0; i < ub; i++) {
                Vector3 pos = curve.nodes[i].interpolate(vertexPositions);
                f << "v " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
            }
            f << "l ";
            for (size_t i = 0; i < nNodes - 1; i++) {
                f << (offset + i) + 1 << " "; // OBJ files are 1-indexed
            }
            f << offset + (nNodes - 1) % ub + 1 << "\n";
            offset += ub;
        }
        f.close();
        std::cerr << "File " << cFilename << " written succesfully." << std::endl;
    } else {
        std::cerr << "Could not save curves '" << cFilename << "'!" << std::endl;
    }
}

void exportCurves(const pointcloud::PointData<Vector3>& positions,
                  const std::vector<std::vector<pointcloud::Point>>& curves, const std::string& dir) {

    std::string cFilename = dir + "/curves.obj";

    std::fstream f;
    f.open(cFilename, std::ios::out | std::ios::trunc);
    if (f.is_open()) {
        // Assume curves have already been organized into components to be "as connected as possible".
        size_t offset = 0;
        for (const auto& curve : curves) {
            size_t nNodes = curve.size();
            bool isClosed = (curve[0] == curve[nNodes - 1]);
            size_t ub = isClosed ? nNodes - 1 : nNodes;
            for (size_t i = 0; i < ub; i++) {
                Vector3 pos = positions[curve[i]];
                f << "v " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
            }
            f << "l ";
            for (size_t i = 0; i < nNodes - 1; i++) {
                f << (offset + i) + 1 << " "; // OBJ files are 1-indexed
            }
            f << offset + (nNodes - 1) % ub + 1 << "\n";
            offset += ub;
        }
        f.close();
        std::cerr << "File " << cFilename << " written succesfully." << std::endl;
    } else {
        std::cerr << "Could not save curves '" << cFilename << "'!" << std::endl;
    }
}

/* https://stackoverflow.com/a/141943 */
Vector3 redistributeRGB(const Vector3& color) {
    double r = color[0];
    double g = color[1];
    double b = color[2];
    double m = std::max(std::max(r, g), b);
    if (m <= 1.) return color;
    double total = r + g + b;
    if (total >= 3.) return {1., 1., 1.};
    double x = (3. - total) / (3. * m - total);
    double gray = 1. - x * m;
    return {gray + x * r, gray + x * g, gray + x * b};
}

/*
 * Create a custom diverging colormap to represent SDFs.
 * Saves a M x N PPM image that has sampled the colormap at N evenly-spaced locations (includes the extrema.)
 * If you have imagemagick, you can batch-convert ppm to png via: mogrify -format png *.ppm

 * Unlike typical (Matplotlib) diverging colormaps, negative/positive colors don't both fade to white near zero,
 * but rather fade in opposite direction to create a very clear demarcation where sign changes.
 *
 * Furthermore, normalize each positive/negative portion of the colormap separately, so they both "fade to white"
 * despite having different ranges. (Using uniform normalization doesn't look good)
 */
void SDFColormapToImageTexture(const Vector3& colorNegative, const Vector3& colorPositive, double split,
                               const std::string& dir, int M, int N) {

    std::vector<double> powScalings = {1.0, 0.8}; // go to white faster / slower
    int negSamples = roundToNearestInt(N * (split));
    int posSamples = N - negSamples;
    std::vector<Vector3> colorsN(negSamples);
    std::vector<Vector3> colorsP(posSamples);
    double tmin = 1.;
    double tmax = 1.5;
    for (int i = 0; i < negSamples; i++) {
        double t = i / (negSamples - 1.);
        t = (1. - t) * tmin + t * tmax;
        t = std::pow(t, powScalings[0]);
        colorsN[i] = redistributeRGB(t * colorNegative);
    }
    for (int i = 0; i < posSamples; i++) {
        double t = i / (posSamples - 1.);
        t = (1. - t) * tmin + t * tmax;
        t = std::pow(t, powScalings[1]);
        colorsP[i] = redistributeRGB(t * colorPositive);
    }

    // always have paler red adjacent to the curve
    std::reverse(colorsN.begin(), colorsN.end());
    std::reverse(colorsP.begin(), colorsP.end());

    std::fstream file;
    std::string filename = dir + "colormap.ppm";
    file.open(filename, std::ios::out | std::ios::trunc);
    if (file.is_open()) {
        file << "P3\n";
        file << N << " " << M << "\n"; // width x height of image in pixels
        file << "255\n";               // max value for each color
        for (int j = 0; j < M; j++) {
            for (int i = 0; i < N; i++) {
                Vector3 color = (i < negSamples) ? colorsN[i] : colorsP[i - negSamples];
                file << int(color[0] * 255.) << " " << int(color[1] * 255.) << " " << int(color[2] * 255.) << "\n";
            }
        }
        std::cerr << "File " << filename << " written succesfully." << std::endl;
        file.close();
    } else {
        std::cerr << "Could not save file '" << filename << "'!" << std::endl;
    }
}

Vector<double> normalizeSDF(const Vector<double>& u, const std::string& dir, bool useBounds, double lowerBound,
                            double upperBound) {

    const double inf = std::numeric_limits<double>::infinity();
    double maxNeg = 0.;
    double maxPos = 0.;
    double minNeg = inf;
    double minPos = inf;
    size_t N = u.size();
    for (size_t i = 0; i < N; i++) {
        if (std::isnan(u[i])) continue;
        double abs_u = std::abs(u[i]);
        if (u[i] <= 0.) {
            minNeg = std::min(minNeg, abs_u);
            maxNeg = std::max(maxNeg, abs_u);
        } else if (u[i] > 0.) {
            minPos = std::min(minPos, abs_u);
            maxPos = std::max(maxPos, abs_u);
        }
    }

    // Export custom colormap.
    double rangeNeg = maxNeg - minNeg;
    double rangePos = maxPos - minPos;
    if (std::isinf(rangeNeg)) rangeNeg = 0.;
    if (std::isinf(rangePos)) rangePos = 0.;

    double range, split;
    if (useBounds) {
        range = upperBound - lowerBound;
        split = (rangeNeg + (-maxNeg - lowerBound)) / range;
    } else {
        range = rangeNeg + rangePos;
        split = rangeNeg / range;
    }
    SDFColormapToImageTexture(COLOR_NEGATIVE, COLOR_POSITIVE, split, dir, 1, 1024);

    Vector<double> w(N);
    for (size_t i = 0; i < N; i++) {
        int s = sgn(u[i]);
        double abs_u = std::abs(u[i]);
        double minVal = (s < 0) ? minNeg : minPos;
        double val = s * ((abs_u - minVal) / range);
        val += split;
        w[i] = val;
    }

    return w;
}

CornerData<Vector2> toTexCoords(const VertexData<double>& u) {

    CornerData<Vector2> texCoords(*u.getMesh());
    for (Vertex v : u.getMesh()->vertices()) {
        for (Corner c : v.adjacentCorners()) {
            texCoords[c] = {u[v], 1.};
        }
    }
    return texCoords;
}

void exportFunction(EmbeddedGeometryInterface& geom, const VertexData<double>& u, const std::string& filename) {

    CornerData<Vector2> texCoords = toTexCoords(u);
    WavefrontOBJ::write(filename, geom, texCoords);
}

void exportFunction(IntegerCoordinatesIntrinsicTriangulation& intTri, VertexPositionGeometry& manifoldGeom,
                    const VertexData<double>& u, const std::string& filename) {

    CommonSubdivision& cs = intTri.getCommonSubdivision();
    cs.constructMesh(true, true);

    // Interpolate from intrinsic mesh to the common subdivision
    VertexData<double> cs_u = cs.interpolateAcrossB(u);
    CornerData<Vector2> texCoords = toTexCoords(cs_u);

    // Export mesh
    ManifoldSurfaceMesh& csMesh = *cs.mesh;
    VertexData<Vector3> csPositions = cs.interpolateAcrossA(manifoldGeom.vertexPositions);
    VertexPositionGeometry geom(csMesh, csPositions);
    writeSurfaceMesh(csMesh, geom, texCoords, filename);
}


void exportSDF(EmbeddedGeometryInterface& geom, const VertexData<double>& u, const std::string& filename,
               bool useBounds, double lowerBound, double upperBound) {

    std::string dir = getHomeDirectory(filename);
    Vector<double> w = normalizeSDF(u.toVector(), dir, useBounds, lowerBound, upperBound);
    exportFunction(geom, VertexData<double>(*u.getMesh(), w), filename);
}

void exportSDF(IntegerCoordinatesIntrinsicTriangulation& intTri, VertexPositionGeometry& manifoldGeom,
               const VertexData<double>& u, const std::string& filename, bool useBounds, double lowerBound,
               double upperBound) {

    std::string dir = getHomeDirectory(filename);
    Vector<double> w = normalizeSDF(u.toVector(), dir, useBounds, lowerBound, upperBound);
    exportFunction(intTri, manifoldGeom, VertexData<double>(*u.getMesh(), w), filename);
}

void exportSDF(const pointcloud::PointData<Vector3>& pointPositions, const pointcloud::PointData<double>& u,
               const std::string& filename) {

    std::string dir = getHomeDirectory(filename);
    Vector<double> w = normalizeSDF(u.toVector(), dir);
    size_t nPoints = pointPositions.size();

    // Write point positions.
    // Store distance values as the x-component in vertex normals.
    std::fstream file;
    file.open(filename, std::ios::out | std::ios::trunc);
    if (file.is_open()) {
        for (size_t i = 0; i < nPoints; i++) {
            Vector3 pos = pointPositions[i];
            file << "v " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
            file << "vn " << w[i] << " " << 0. << " " << 0. << "\n";
        }
        std::cerr << "File " << filename << " written succesfully." << std::endl;
        file.close();
    } else {
        std::cerr << "Could not save file '" << filename << "'!" << std::endl;
    }
}

// ===================== MESH MUTATION


std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>
createIntrinsicTriangulation(VertexPositionGeometry& geometry, ManifoldSurfaceMesh& mesh,
                             const std::vector<Curve>& curves) {

    std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation> intTri =
        std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>(
            new IntegerCoordinatesIntrinsicTriangulation(mesh, geometry));

    // If curve is given as a set of mesh edges, set these edges as "marked" in the intrinsic triangulation so they
    // never get flipped.
    if (curves.size() > 0) {

        EdgeData<bool> markedEdges(*intTri->intrinsicMesh, false);

        geometry.requireEdgeIndices();
        for (const Curve& curve : curves) {
            size_t nNodes = curve.nodes.size();
            for (size_t i = 0; i < nNodes - 1; i++) {
                const SurfacePoint& pA = curve.nodes[i];
                const SurfacePoint& pB = curve.nodes[i + 1];
                Halfedge he = determineHalfedgeFromVertices(pA.vertex, pB.vertex);
                size_t eIdx = geometry.edgeIndices[he.edge()];
                markedEdges[intTri->intrinsicMesh->edge(eIdx)] = true;
            }
        }
        geometry.unrequireEdgeIndices();

        intTri->setMarkedEdges(markedEdges);
    }
    return intTri;
}

void setIntrinsicSolver(VertexPositionGeometry& geometry, const std::vector<Curve>& curves,
                        const std::vector<SurfacePoint>& points, std::unique_ptr<ManifoldSurfaceMesh>& manifoldMesh,
                        std::unique_ptr<VertexPositionGeometry>& manifoldGeom, std::vector<Curve>& curvesOnManifold,
                        std::vector<SurfacePoint>& pointsOnManifold,
                        std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>& intTri,
                        std::unique_ptr<SignedHeatSolver>& solver) {

    // Set up corresponding manifold mesh & geometry (necessary for intrinsic triangulation operations.)
    SurfaceMesh& mesh = geometry.mesh;
    if (!mesh.isManifold() || !mesh.isOriented())
        throw std::logic_error(
            "Mesh must be manifold and orientable to convert to geometrycentral::ManifoldSurfaceMesh.");

    manifoldMesh = mesh.toManifoldMesh();
    manifoldMesh->compress();

    manifoldGeom = geometry.reinterpretTo(*manifoldMesh);
    manifoldGeom->refreshQuantities();

    intTri = std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>(
        new IntegerCoordinatesIntrinsicTriangulation(*manifoldMesh, *manifoldGeom));

    // If curve is given as a set of mesh edges, set these edges as "marked" in the intrinsic triangulation so they
    // never get flipped.
    curvesOnManifold.clear();
    pointsOnManifold.clear();
    if (curves.size() > 0) {

        EdgeData<bool> markedEdges(*intTri->intrinsicMesh, false);

        // Re-interpret curve on the newly-constructed manifold mesh.
        std::vector<Halfedge> curveHalfedgesOnManifold;
        geometry.requireVertexIndices();
        for (const Curve& curve : curves) {
            curvesOnManifold.emplace_back();
            curvesOnManifold.back().isSigned = curve.isSigned;
            std::vector<Halfedge> currHalfedges;
            size_t nNodes = curve.nodes.size();
            for (size_t i = 0; i < nNodes - 1; i++) {
                const SurfacePoint& pA = curve.nodes[i];
                const SurfacePoint& pB = curve.nodes[i + 1];
                if (pA.type != SurfacePointType::Vertex || pB.type != SurfacePointType::Vertex) {
                    throw std::logic_error("setIntrinsicSolver(): Curves must be constrained to mesh edges.");
                }
                Vertex vA = manifoldMesh->vertex(geometry.vertexIndices[pA.vertex]);
                Vertex vB = manifoldMesh->vertex(geometry.vertexIndices[pB.vertex]);
                Halfedge he = determineHalfedgeFromVertices(vA, vB);
                if (he == Halfedge()) {
                    std::cerr << "Could not find halfedge between vertices " << vA << " and " << vB << std::endl;
                    continue;
                }
                currHalfedges.push_back(he);
            }
            curveHalfedgesOnManifold.insert(curveHalfedgesOnManifold.end(), currHalfedges.begin(), currHalfedges.end());
            for (size_t i = 0; i < nNodes; i++) {
                curvesOnManifold.back().nodes.push_back(reinterpretTo(curve.nodes[i], *manifoldMesh));
            }
        }
        geometry.unrequireVertexIndices();

        manifoldGeom->requireEdgeIndices();
        for (Halfedge he : curveHalfedgesOnManifold) {
            size_t eIdx = manifoldGeom->edgeIndices[he.edge()];
            markedEdges[intTri->intrinsicMesh->edge(eIdx)] = true;
        }
        manifoldGeom->unrequireEdgeIndices();

        intTri->setMarkedEdges(markedEdges);
    }
    for (const SurfacePoint& p : points) {
        SurfacePoint pt = reinterpretTo(p, *manifoldMesh);
        pointsOnManifold.push_back(pt);
    }

    solver = std::unique_ptr<SignedHeatSolver>(new SignedHeatSolver(*intTri));
}

void determineSourceGeometryOnIntrinsicTriangulation(IntrinsicTriangulation& intTri,
                                                     const std::vector<Curve>& curvesOnManifold,
                                                     const std::vector<SurfacePoint>& pointsOnManifold,
                                                     std::vector<Curve>& curvesOnIntrinsic,
                                                     std::vector<SurfacePoint>& pointsOnIntrinsic) {

    // For some reason, need to call this first or else everything fails. Intrinsic mesh isn't built/indexed right.
    intTri.traceAllIntrinsicEdgesAlongInput();

    // Use traceInputHalfedgeAlongIntrinsic(), and verify that all SurfacePoints are vertex-type.
    curvesOnIntrinsic.clear();
    pointsOnIntrinsic.clear();
    for (const Curve& curve : curvesOnManifold) {
        curvesOnIntrinsic.emplace_back();
        curvesOnIntrinsic.back().isSigned = curve.isSigned;
        size_t nNodes = curve.nodes.size();
        for (size_t i = 0; i < nNodes - 1; i++) {
            const SurfacePoint& pA = curve.nodes[i];
            const SurfacePoint& pB = curve.nodes[i + 1];
            Halfedge he = determineHalfedgeFromVertices(pA.vertex, pB.vertex);
            std::vector<SurfacePoint> pts = intTri.traceInputHalfedgeAlongIntrinsic(he);
            for (SurfacePoint pt : pts) assert(pt.type == SurfacePointType::Vertex);
            size_t n = pts.size();
            for (size_t j = 0; j < n - 1; j++) {
                curvesOnIntrinsic.back().nodes.push_back(pts[j]);
            }
            if (i == nNodes - 2 && curve.nodes[0] != curve.nodes[nNodes - 1]) {
                curvesOnIntrinsic.back().nodes.push_back(pts[n - 1]);
            }
        }
        if (curvesOnIntrinsic.back().nodes.size() < 2) curvesOnIntrinsic.pop_back();
    }

    for (const SurfacePoint& p : pointsOnManifold) {
        SurfacePoint pt = reinterpretTo(p, *(intTri.intrinsicMesh));
        pointsOnIntrinsic.push_back(pt);
    }
}

// ================ VISUALIZATION

void displayInput(const VertexData<Vector3>& vertexPositions, const std::vector<Curve>& curves,
                  const std::vector<SurfacePoint>& pointSources, const std::string& name, bool display) {

    std::vector<Vector3> nodes;
    std::vector<std::array<size_t, 2>> edges;
    std::vector<Vector3> edgeColors;

    // Display signed geometry in orange, unsigned geometry in black.
    Vector3 signedColor = {1, 0.3, 0};
    Vector3 unsignedColor = {0, 0, 0};
    size_t offset = 0;
    for (const Curve& curve : curves) {
        size_t N = curve.nodes.size();
        for (size_t i = 0; i < N - 1; i++) {
            nodes.push_back(curve.nodes[i].interpolate(vertexPositions));
            edges.push_back({offset + i, offset + i + 1});
            edgeColors.push_back(curve.isSigned ? signedColor : unsignedColor);
        }
        nodes.push_back(curve.nodes[N - 1].interpolate(vertexPositions));
        offset += N;
    }
    for (const SurfacePoint& p : pointSources) {
        nodes.push_back(p.interpolate(vertexPositions));
        edges.push_back({offset, offset});
        edgeColors.push_back(unsignedColor);
        offset += 1;
    }

    auto psCurve = polyscope::registerCurveNetwork(name, nodes, edges);
    // psCurve->setRadius(0.001);
    psCurve->setEnabled(display);
    psCurve->addEdgeColorQuantity(name, edgeColors)->setEnabled(true);
}

void displayInput(const VertexData<Vector3>& vertexPositions, const std::vector<std::vector<Vertex>>& curves,
                  const std::string& name, bool display) {

    std::vector<Vector3> nodes;
    std::vector<std::array<size_t, 2>> edges;

    size_t offset = 0;
    for (const auto& curve : curves) {
        size_t N = curve.size();
        for (size_t i = 0; i < N - 1; i++) {
            nodes.push_back(vertexPositions[curve[i]]);
            edges.push_back({offset + i, offset + i + 1});
        }
        nodes.push_back(vertexPositions[curve[N - 1]]);
        offset += N;
    }

    auto psCurve = polyscope::registerCurveNetwork(name, nodes, edges);
    psCurve->setEnabled(display);
    psCurve->setColor({1, 0.3, 0});
}

void displayInput(const pointcloud::PointData<Vector3>& positions,
                  const std::vector<std::vector<pointcloud::Point>>& curves, const std::string& name, bool display) {

    std::vector<Vector3> nodes;
    std::vector<std::array<size_t, 2>> edges;

    size_t offset = 0;
    for (const auto& curve : curves) {
        size_t N = curve.size();
        for (size_t i = 0; i < N - 1; i++) {
            nodes.push_back(positions[curve[i]]);
            edges.push_back({offset + i, offset + i + 1});
        }
        nodes.push_back(positions[curve[N - 1]]);
        offset += N;
    }

    auto psCurve = polyscope::registerCurveNetwork(name, nodes, edges);
    psCurve->setEnabled(display);
    psCurve->setColor({1, 0.3, 0});
}

void visualizeIntrinsicEdges(IntegerCoordinatesIntrinsicTriangulation& intTri, VertexPositionGeometry& manifoldGeom,
                             bool display) {

    std::vector<Vector3> nodes;
    std::vector<std::array<size_t, 2>> edges;
    EdgeData<std::vector<SurfacePoint>> traces = intTri.traceAllIntrinsicEdgesAlongInput();
    for (Edge e : intTri.intrinsicMesh->edges()) {
        size_t N = nodes.size();
        size_t M = traces[e].size();
        for (const auto& pt : traces[e]) {
            nodes.push_back(pt.interpolate(manifoldGeom.vertexPositions));
        }
        for (size_t i = 0; i < M - 1; i++) edges.push_back({N + i, N + i + 1});
    }

    auto intEdgeQ = polyscope::registerCurveNetwork("intrinsic edges", nodes, edges);
    intEdgeQ->setEnabled(display);
    intEdgeQ->setColor(polyscope::render::RGB_ORANGE);
    intEdgeQ->setRadius(0.0005);
}

/*
 * Display VertexData quantity on the common subdivision of an intrinsic triangulation.
 */
VertexData<double> visualizeOnCommonSubdivision(IntegerCoordinatesIntrinsicTriangulation& intTri,
                                                VertexPositionGeometry& manifoldGeom, VertexData<Vector3>& csPositions,
                                                std::unique_ptr<VertexPositionGeometry>& csGeom,
                                                const VertexData<double>& w, const std::string& name, bool rebuild,
                                                bool divergingColormap) {

    CommonSubdivision& cs = intTri.getCommonSubdivision();
    cs.constructMesh();
    // Linearly interpolate data from intrinsic mesh to the common subdivision.
    VertexData<double> cs_w = cs.interpolateAcrossB(w);

    ManifoldSurfaceMesh& csMesh = *cs.mesh;
    csPositions = cs.interpolateAcrossA(manifoldGeom.vertexPositions);
    csGeom.reset(new VertexPositionGeometry(csMesh, csPositions));

    if (rebuild) {
        // Register and display with Polyscope
        polyscope::SurfaceMesh* psCsMesh =
            polyscope::registerSurfaceMesh("common subdivision", csPositions, csMesh.getFaceVertexList());
        psCsMesh->setAllPermutations(polyscopePermutations(csMesh));
        psCsMesh->setEnabled(true);
    }
    if (divergingColormap) {
        polyscope::getSurfaceMesh("common subdivision")->addVertexSignedDistanceQuantity(name, cs_w)->setEnabled(true);
    } else {
        polyscope::getSurfaceMesh("common subdivision")
            ->addVertexScalarQuantity(name, cs_w)
            ->setIsolinesEnabled(true)
            ->setEnabled(true);
    }
    return cs_w;
}