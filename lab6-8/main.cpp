#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "Vector.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <queue>
#include <random>
#include <omp.h>
#include <chrono>
#include "lbfgs.h"
#include <sstream>

using namespace std;

#define N_POINTS 50
#define FLUID 0.2

#define PRINT_VEC(v) (printf("%s: (%f %f %f)\n", #v, (v)[0], (v)[1], (v)[2]))

class Polygon
{
public:
    vector<Vector> vertices;
    Polygon()
    {
        this->vertices = {};
    }

    Polygon(vector<Vector> &vertices)
    {
        this->vertices = vertices;
    }

    void print()
    {
        for (auto v : vertices)
        {
            cout << "v " << v[0] << " " << v[1] << " ";
        }
        cout << endl;
    }

    double area()
    {
        double sum = 0;
        if (vertices.size() <= 2)
        {
            return 0;
        };
        for (int i = 0; i < vertices.size(); i++)
        {
            sum += (vertices[i][0] * vertices[i == vertices.size() - 1 ? 0 : i + 1][1]) - (vertices[i][1] * vertices[i == vertices.size() - 1 ? 0 : i + 1][0]);
        }
        return abs(sum) / 2.;
    }

    double integrate(Vector Pi)
    {
        double sum = 0;
        if (vertices.size() <= 2)
        {
            // cout << "small size" << endl;
            return 0;
        }
        int n = vertices.size();
        // for (int i = 1; i < vertices.size(); i++)
        // {
        //     double tmp_X = vertices[i == 0 ? n - 1 : i - 1][0] * vertices[i == 0 ? n - 1 : i - 1][0] + vertices[i == 0 ? n - 1 : i - 1][0] * vertices[i][0] + vertices[i][0] * vertices[i][0];
        //     double tmp_Y = vertices[i == 0 ? n - 1 : i - 1][1] * vertices[i == 0 ? n - 1 : i - 1][1] + vertices[i == 0 ? n - 1 : i - 1][1] * vertices[i][1] + vertices[i][1] * vertices[i][1];
        //     double tmp = (Pi[0] * (vertices[i == 0 ? n - 1 : i - 1][0] + vertices[i][0]) + Pi[1] * (vertices[i == 0 ? n - 1 : i - 1][1] + vertices[i][1]));
        //     sum += (vertices[i == 0 ? n - 1 : i - 1][0] * vertices[i][1] - vertices[i][1] * vertices[i == 0 ? n - 1 : i - 1][0]) * (tmp_X + tmp_Y - 4 * (tmp + 6 * Pi.norm2()));
        // }
        Vector c0 = vertices[0], c1, c2;
        for (int i = 1; i < vertices.size() - 1; i++)
        {
            c1 = vertices[i];
            c2 = vertices[i + 1];
            vector<Vector> cs = {c0, c1, c2};
            double S = 0.5 * cross(c1 - c0, c2 - c0).norm();
            for (int l = 0; l < 3; l++)
            {
                for (int k = 0; k <= l; k++)
                {
                    sum += dot(cs[k] - Pi, cs[l] - Pi) * S;
                }
            }
        }
        return sum / 6.;
    }

    Vector centroid()
    {
        Vector c;
        for (auto v : vertices)
        {
            c = c + v;
        }
        return c / (double)vertices.size();
    }
};

class Line
{
public:
    Vector u, N;
    Line(Vector u, Vector N)
    {
        N.normalize();
        this->u = u;
        this->N = N;
    }
};

void save_svg(const vector<Polygon> &polygons, string filename, string fillcol = "none")
{
    FILE *f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i = 0; i < polygons.size(); i++)
    {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].vertices.size(); j++)
        {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

void save_svg_w_points(const vector<Polygon> &polygons, const vector<Vector> &points, const string &filename, const string &fillcol = "none", const std::string &pointcol = "red")
{
    FILE *f = fopen(filename.c_str(), "w+");
    if (!f)
    {
        fprintf(stderr, "Error opening file %s for writing\n", filename.c_str());
        return;
    }

    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");

    // Draw polygons
    for (const auto &polygon : polygons)
    {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (const auto &vertex : polygon.vertices)
        {
            fprintf(f, "%3.3f, %3.3f ", (vertex[0] * 1000), (1000 - vertex[1] * 1000));
        }
        fprintf(f, "\" fill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }

    // Draw points
    for (const auto &point : points)
    {
        fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r = \"5\" fill = \"%s\" />\n", (point[0] * 1000), (1000 - point[1] * 1000), pointcol.c_str());
    }

    fprintf(f, "</svg>\n");
    fclose(f);
}
// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg_animated(const vector<Polygon> &polygons, string filename, int frameid, int nbframes)
{
    FILE *f;
    if (frameid == 0)
    {
        f = fopen(filename.c_str(), "w+");
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        fprintf(f, "<g>\n");
    }
    else
    {
        f = fopen(filename.c_str(), "a+");
    }
    fprintf(f, "<g>\n");
    for (int i = 0; i < polygons.size(); i++)
    {
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].vertices.size(); j++)
        {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
    }
    fprintf(f, "<animate\n");
    fprintf(f, "    id = \"frame%u\"\n", frameid);
    fprintf(f, "    attributeName = \"display\"\n");
    fprintf(f, "    values = \"");
    for (int j = 0; j < nbframes; j++)
    {
        if (frameid == j)
        {
            fprintf(f, "inline");
        }
        else
        {
            fprintf(f, "none");
        }
        fprintf(f, ";");
    }
    fprintf(f, "none\"\n    keyTimes = \"");
    for (int j = 0; j < nbframes; j++)
    {
        fprintf(f, "%2.3f", j / (double)(nbframes));
        fprintf(f, ";");
    }
    fprintf(f, "1\"\n   dur = \"5s\"\n");
    fprintf(f, "    begin = \"0s\"\n");
    fprintf(f, "    repeatCount = \"indefinite\"/>\n");
    fprintf(f, "</g>\n");
    if (frameid == nbframes - 1)
    {
        fprintf(f, "</g>\n");
        fprintf(f, "</svg>\n");
    }
    fclose(f);
}

/*
TODO:
- Calculate G
- Calculate lamda
- Calculate integral
- Calculate area of polygon
*/

const Vector C = Vector(0.5, 0.5, 0);

double lambda(Vector y)
{
    return exp(-(y - C).norm2() / 0.02);
}

double objectfunc(vector<double> W, vector<Vector> &points, vector<Polygon> &Pow)
{
    double sum = 0;
    for (int i = 0; i < W.size(); i++)
    {
        sum += Pow[i].integrate(points[i]) - W[i] * Pow[i].area() + lambda(points[i]) * W[i];
    }
    return sum;
}

Vector edgeLineIntersect(Vector &prevVertex, Vector &curVertex, Line clipEdge)
{
    double t = dot(clipEdge.u - prevVertex, clipEdge.N) / dot(curVertex - prevVertex, clipEdge.N);
    if (t < 0 || t > 1)
    {
        // cerr << "No intersection" << endl;
        return Vector(INFINITY, INFINITY, INFINITY);
    }
    return prevVertex + t * (curVertex - prevVertex);
}

Polygon polygonClip(Polygon subjectPolygon, Polygon clipPolygon)
{
    Polygon outPolygon;
    for (int i = 0; i < clipPolygon.vertices.size(); i++)
    {
        outPolygon = Polygon();
        Vector dir, normal;
        if (i < clipPolygon.vertices.size() - 1)
            dir = clipPolygon.vertices[i + 1] - clipPolygon.vertices[i];
        else
            dir = clipPolygon.vertices[0] - clipPolygon.vertices[i];
        normal = Vector();
        normal[0] = dir[1];
        normal[1] = -dir[0];
        normal.normalize();
        Line clipEdge = Line(clipPolygon.vertices[i], normal);

        for (int i = 0; i < subjectPolygon.vertices.size(); i++)
        {
            Vector curVertex = subjectPolygon.vertices[i];
            Vector prevVertex = subjectPolygon.vertices[i > 0 ? i - 1 : subjectPolygon.vertices.size() - 1];
            Vector intersection = edgeLineIntersect(prevVertex, curVertex, clipEdge);
            if (dot(curVertex - clipEdge.u, clipEdge.N) <= 0)
            {
                if (dot(prevVertex - clipEdge.u, clipEdge.N) > 0)
                {
                    outPolygon.vertices.push_back(intersection);
                }
                outPolygon.vertices.push_back(curVertex);
            }
            else if (dot(prevVertex - clipEdge.u, clipEdge.N) <= 0)
            {
                outPolygon.vertices.push_back(intersection);
            }
        }
        subjectPolygon = outPolygon;
    }
    return outPolygon;
}

Polygon constructDisk(Vector center, double R, int n_sample = 80)
{
    double d_angle = 2. * M_PI / (double)n_sample;
    double angle = M_PI;
    Polygon outPolygon;
    for (int i = 0; i < n_sample; i++)
    {
        double x, y;
        x = center[0] + R * cos(angle);
        y = center[1] + R * sin(angle);
        angle += d_angle;
        outPolygon.vertices.push_back(Vector(x, y, 0));
    }
    return outPolygon;
}

Polygon VoronoiPolygon(vector<Vector> &pointSet, vector<double> &weights, int index)
{
    vector<Vector> boundVertices = {Vector(0, 0, 0),
                                    Vector(0, 1, 0),
                                    Vector(1, 1, 0),
                                    Vector(1, 0, 0)};
    Polygon subjectPolygon = Polygon(boundVertices);
    Polygon outPolygon;
    int n = weights.size();

    for (int i = 0; i < pointSet.size(); i++)
    {
        if (i == index)
        {
            continue;
        }
        double maxDist = 0;

        outPolygon = Polygon();
        Vector M = (pointSet[i] + pointSet[index]) / 2.;
        M = M + (weights[index] - weights[i]) * (pointSet[i] - pointSet[index]) / (2. * (pointSet[index] - pointSet[i]).norm2());

        Line clipEdge = Line(M, pointSet[i] - pointSet[index]);

        for (int j = 0; j < subjectPolygon.vertices.size(); j++)
        {
            Vector curVertex = subjectPolygon.vertices[j];
            Vector prevVertex = subjectPolygon.vertices[j > 0 ? j - 1 : subjectPolygon.vertices.size() - 1];
            Vector intersection = edgeLineIntersect(prevVertex, curVertex, clipEdge);
            if ((curVertex - pointSet[index]).norm2() - weights[index] <= (curVertex - pointSet[i]).norm2() - weights[i])
            {
                if ((prevVertex - pointSet[index]).norm2() - weights[index] > (prevVertex - pointSet[i]).norm2() - weights[i])
                {
                    outPolygon.vertices.push_back(intersection);
                }
                outPolygon.vertices.push_back(curVertex);
            }
            else if ((prevVertex - pointSet[index]).norm2() - weights[index] <= (prevVertex - pointSet[i]).norm2() - weights[i])
            {
                outPolygon.vertices.push_back(intersection);
            }
        }
        subjectPolygon = outPolygon;
    }
    return outPolygon;
}

class objective_function
{
protected:
public:
    vector<Vector> points;
    vector<double> weights;
    lbfgsfloatval_t *m_x;

    objective_function(vector<Vector> &points)
    {
        this->points = points;
        m_x = NULL;
    }

    objective_function() : m_x(NULL)
    {
    }

    virtual ~objective_function()
    {
        if (m_x != NULL)
        {
            lbfgs_free(m_x);
            m_x = NULL;
        }
    }

    int run(int N)
    {
        lbfgsfloatval_t fx;
        lbfgsfloatval_t *m_x = lbfgs_malloc(N);

        if (m_x == NULL)
        {
            printf("ERROR: Failed to allocate a memory block for variables.\n");
            return 1;
        }

        /* Initialize the variables. */
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.8, 1.);
        weights.reserve(N);
        for (int i = 0; i < N; ++i)
        {
            weights[i] = dis(gen);
        }
        weights[N - 1] = 0.6;

        /*
            Start the L-BFGS optimization; this will invoke the callback functions
            evaluate() and progress() when necessary.
         */
        int ret = lbfgs(N, &weights[0], &fx, _evaluate, _progress, this, NULL);

        /* Report the result. */
        printf("L-BFGS optimization terminated with status code = %d\n", ret);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, weights[0], weights[1]);

        return ret;
    }

protected:
    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step)
    {
        return reinterpret_cast<objective_function *>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step)
    {
        lbfgsfloatval_t fx = 0.0;

        int n_fluid = n - 1;

        vector<double> lambdas(n_fluid);
        vector<Polygon> polygons(n_fluid);
        double sum_lambdas = 0;
        double estimated_vfluid = 0;
        for (int i = 0; i < n_fluid; i += 1)
        {
            lambdas[i] = lambda(points[i]);
            sum_lambdas += lambdas[i];
            Polygon Pow = VoronoiPolygon(points, weights, i);
            double R = sqrt(weights[i] - weights[n - 1]);
            polygons[i] = polygonClip(Pow, constructDisk(points[i], R, 200));
            estimated_vfluid += polygons[i].area();
        }

        for (int i = 0; i < n_fluid; i += 1)
        {
            // Pow.print()
            g[i] = polygons[i].area() - FLUID / n_fluid;
            // cout << polygons[i].integrate(points[i]) << " " << polygons[i].area() << " " << sum_lambdas << endl;
            fx += polygons[i].integrate(points[i]) - weights[i] * polygons[i].area() + (FLUID / n_fluid) * weights[i];
        }
        g[n - 1] = (1 - estimated_vfluid) - (1 - FLUID);
        fx += weights[n - 1] * ((1 - FLUID) - (1 - estimated_vfluid));
        return -fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls)
    {
        return reinterpret_cast<objective_function *>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls)
    {
        // printf("Iteration %d:\n", k);
        // printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        // printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        // printf("\n");
        return 0;
    }
};

int sgn(double d)
{
    if (d == 0)
        return 0;
    if (d < 0)
        return -1;
    return 1;
}

void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0)
{
    int W = 500, H = 500;
    std::vector<unsigned char> image(W * H * 3, 255);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++)
    {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++)
        {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W - 1., std::max(0., W * bminx));
        bminy = std::min(H - 1., std::max(0., H * bminy));
        bmaxx = std::max(W - 1., std::max(0., W * bmaxx));
        bmaxy = std::max(H - 1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++)
        {
            for (int x = bminx; x < bmaxx; x++)
            {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++)
                {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0) * (y1 - y0) - (y - y0) * (x1 - x0);
                    int sign = sgn(det);
                    if (prevSign == 0)
                        prevSign = sign;
                    else if (sign == 0)
                        sign = prevSign;
                    else if (sign != prevSign)
                    {
                        isInside = false;
                        break;
                    }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
                    double distEdge = std::abs(det) / edgeLen;
                    double dotp = (x - x0) * (x1 - x0) + (y - y0) * (y1 - y0);
                    if (dotp < 0 || dotp > edgeLen * edgeLen)
                        distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside)
                {
                    if (i < N_POINTS)
                    { // the N first particles may represent fluid, displayed in blue
                        image[((H - y - 1) * W + x) * 3] = 0;
                        image[((H - y - 1) * W + x) * 3 + 1] = 0;
                        image[((H - y - 1) * W + x) * 3 + 2] = 255;
                    }
                    if (mindistEdge <= 2)
                    {
                        image[((H - y - 1) * W + x) * 3] = 0;
                        image[((H - y - 1) * W + x) * 3 + 1] = 0;
                        image[((H - y - 1) * W + x) * 3 + 2] = 0;
                    }
                }
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

void GMScheme(vector<Vector> &points, vector<Vector> &v, double mass, int max_num_frames = 50)
{
    for (int i = 0; i < max_num_frames; i++)
    {
        objective_function obj(points);
        obj.run(N_POINTS + 1);
        double eps = 0.004;
        double dt = 0.008;
        Vector g = Vector(0, -9.8, 0);

        vector<Polygon> voronoi(points.size());
        double sum_area = 0;
        vector<Vector> new_points(N_POINTS);

        for (int i = 0; i < points.size(); i++)
        {
            voronoi[i] = VoronoiPolygon(points, obj.weights, i);
            double R = sqrt(obj.weights[i] - obj.weights[N_POINTS]);
            voronoi[i] = polygonClip(voronoi[i], constructDisk(points[i], R));

            Vector F_spring = (voronoi[i].centroid() - points[i]) / (eps * eps);
            // PRINT_VEC(F_spring);
            Vector F = F_spring + mass * g;
            v[i] = v[i] + F * dt / mass;
            new_points[i] = points[i] + v[i] * dt;
        }

        points = new_points;

        // printf("Frame %d done!\n", i);

        save_svg_animated(voronoi, "./fluid.svg", i, max_num_frames);
    }
}

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Vector> points;
    std::vector<Vector> vs;
    points.reserve(N_POINTS); // Reserve space for 500 points
    vs.reserve(N_POINTS);
    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Generate 500 random points
    for (int i = 0; i < N_POINTS; ++i)
    {
        double x = dis(gen);
        double y = dis(gen);
        points.emplace_back(x, y);
        vs.emplace_back(Vector(0, 0, 0));
    }

    GMScheme(points, vs, 200, 200);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
}