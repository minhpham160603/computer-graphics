# CSE306 - Project 2

## Minh Pham

### Lecture 6: Geometry Processing - Clipping and Voronoi Diagrams

#### Key Implementations

1. **Polygon Clipping**
   - **Function:** `polygonClip`
   - **Description:** This function receives a subject polygon and a clip polygon. The output is a polygon clipped from the subject polygon, conforming to the shape of the clip polygon.

2. **Voronoi Polygon Construction**
   - **Function:** `VoronoiPolygon`
   - **Description:** This function takes a set of points and a specific index of a point, then constructs the Voronoi polygon around that point. The starting polygon is the bounding square of the whole surface (0,0) to (1,1). This is iteratively updated per iteration.
     - The outer loop iterates through every edge between the point index and every other point in the set.
     - The inner loop checks the intersection between this line and every edge of the subject polygon.
   - **Voronoi Diagram Construction:** The Voronoi diagram is constructed by iteratively performing the Voronoi Polygon function for all points in the set. OpenMP can be used at this stage as each Voronoi polygon of a point is independent of the others.

#### Results

- **Voronoi Diagrams Generation:**
  - Generated Voronoi diagrams for a set of 1000 random points.
  - The Sutherland-Hodgman algorithm clipped the polygons to construct the Voronoi cells.
  - **Performance:**
    - 1000 points Voronoi Diagram takes 0.023 seconds to complete with OpenMP and 0.15 seconds without OpenMP.

### Lecture 7: Geometry Processing - Semi-Discrete Optimal Transport in 2D

#### Key Implementations

1. **Power Diagram:**
   - Instead of taking the line going through the center of the edge between point index and other points, point M is updated according to additional weights (see lecture notes for the formula).

2. **Weight Update:**
   - The next step is to update the weights to match the requirements (the ratio of the polygon’s area). This step requires the integration of the L-BFGS library to optimize an objective function.

3. **Polygon Area Calculation:**
   - Implemented the `area()` method to calculate the area of the polygon and the `integrate(Vector P)` method to calculate the integral over the polygon surface of a given point P inside the polygon. Used the triangulation integral for accuracy.

4. **L-BFGS Optimization:**
   - Implemented the `objective_function` class following the C++ sample on their repository.
   - The `run()` method calls the `lbfgs` function that runs the optimization routine.
   - The `evaluate()` method computes the value of the objective function (fx) and its derivative (g) with respect to each variable, following the lecture notes formula. Both fx and g are negated as we aim to maximize the objective function.

#### Results

- **Power Diagram + L-BFGS Optimization:**
  - **Performance:**
    - 100 points take 0.61 seconds to render with OpenMP and O3 compiling flag.
    - 500 points take 19.16 seconds.

### Lab 8 Report: Optimizing Voronoi Diagram for Fluid Simulation

#### Key Implementations

1. **Fluid Simulation:**
   - The Power Diagram with L-BFGS already implemented simplifies fluid simulation.
   - Each polygon is clipped with a disk of radius $R = \sqrt{\text{weights[i]} - \text{weights[air]}}$.
   - Created a `constructDisk` method to change the number of sample points per disk later on.

2. **Objective Function Modification:**
   - Modified the `evaluate()` method in the `objective_function` class to accommodate fluid volume proportion (FLUID) in the whole surface.

3. **Physics Integration:**
   - Added the physics formula into a `GMScheme` function that iteratively modifies the position of the particle and recomputes the diagram.

#### Results

- **Fluid Simulation:**
  - 50 points with $p_{\text{fluid}} = 0.2$, mass = 200, 200 frames, rendered in 30 seconds.
  - **Example:**
<img src="/lab6-8/fluid.svg" alt="Test Polygon Clip Algorithm Input" style="margin-right: 10px; width=30px;height=auto">

