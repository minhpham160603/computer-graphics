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

using namespace std;

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

Polygon VoronoiPolygon(vector<Vector> &pointSet, vector<double> &weights, int index)
{
    vector<Vector> boundVertices = {Vector(0, 0, 0),
                                    Vector(0, 1, 0),
                                    Vector(1, 1, 0),
                                    Vector(1, 0, 0)};
    Polygon subjectPolygon = Polygon(boundVertices);
    Polygon outPolygon;

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
        std::uniform_real_distribution<> dis(0.0, 0.1);
        weights.reserve(N);
        for (int i = 0; i < N; ++i)
        {
            weights[i] = 1;
        }

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

        vector<double> lambdas(n);
        double sum_lambdas = 0;
        for (int i = 0; i < n; i += 1)
        {
            lambdas[i] = lambda(points[i]);
            sum_lambdas += lambdas[i];
        }

        for (int i = 0; i < n; i += 1)
        {
            Polygon Pow = VoronoiPolygon(points, weights, i);
            // Pow.print();
            g[i] = Pow.area() - (lambdas[i] / sum_lambdas);
            // cout << Pow.integrate(points[i]) << " " << Pow.area() << " " << sum_lambdas << endl;
            fx += Pow.integrate(points[i]) - weights[i] * Pow.area() + (lambdas[i] / sum_lambdas) * weights[i];
        }
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

#define N 500

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Vector> points;
    std::vector<double> weights;
    points.reserve(N); // Reserve space for 500 points
    weights.reserve(N);
    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Generate 500 random points
    for (int i = 0; i < N; ++i)
    {
        double x = dis(gen);
        double y = dis(gen);
        points.emplace_back(x, y);
        weights.emplace_back(dis(gen));
    }

    objective_function obj(points);
    obj.run(N);

    vector<Polygon> voronoi(points.size());
    double sum_area = 0;
    for (int i = 0; i < points.size(); i++)
    {
        voronoi[i] = (VoronoiPolygon(points, obj.weights, i));
        sum_area += voronoi[i].area();
    }

    cout << "area " << sum_area << endl;
    // save_svg({subject, clip}, "test.svg");

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";

    save_svg_w_points(voronoi, points, "power_diagram.svg");
}