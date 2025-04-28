# Certified Curve Tracking
This is a Git repo for the paper "Certified algebraic curve projections by path tracking" by Burr, Byrd and Lee
Certified numerical tracking of solution curves for polynomial systems, with interval arithmetic validation and LaTeX/TikZ output.

## How to run codes

* Clone or download this repository to your local machine.
* Open your Julia REPL, and move to the directory where the `certified_curve_tracking.jl` file is located.
* Load the main tracking library:

```julia
include("certified_curve_tracking.jl")
```

* Prepare your polynomial system `F`, initial point `p`, and initial radius `r`, following the examples below.
* Run the tracking with the function `track_curve`.

## Running an Example

Example (Figure 3 style):

```julia
CCi = ComplexField()
R, (x,y,z) = CCi["x","y","z"]
F = hcat([-z-x^3+2.7x, y^2-2+z])
p = [CCi(2.3947), CCi(3.17), CCi(-8.04)]

x, r, A = track_curve(F, p, 0.1, 600, "~/Documents/GitHub/certified_curve_projections/small_example";
    show_tubular_neighborhood=true,
    box_thickness=0.2,
    line_thickness=0.1
)
```

This generates a TikZ `.tex` file under the specified path, visualizing the certified curve.

* The `.tex` file can be compiled using any LaTeX engine that supports TikZ.

## Input format

There is no separate input file required. You manually define the following inside Julia:

### Polynomials

* Define your polynomial system using variables over `ComplexField()`.

```julia
R, (x,y,z) = CCi["x","y","z"]
```

* The system `F` must be given as a Julia array using `hcat([...])`.

Example:

```julia
F = hcat([-z - x^3 + 2.7x, y^2 - 2 + z])
```

### Points

* The starting point `p` must be a Julia array of complex numbers.

Example:

```julia
p = [CCi(2.3947), CCi(3.17), CCi(-8.04)]
```

* Each coordinate corresponds to the order of variables in the polynomial system.

## Output file format

* A `.tex` file is generated automatically by `track_curve`.
* It contains:
  - Red lines connecting certified midpoints (solution path)
  - (Optional) Blue boxes representing the certified tubular neighborhoods.

Example snippet from output:

```latex
\draw[color=red,line width=.2mm] (x1,y1) -- (x2,y2);
\draw[color=blue,line width=.1mm] (box edges)
```

* You can compile this `.tex` file using `pdflatex` or Overleaf to get a high-quality figure.

## Examples included

* **Figure 3 Example**: Small simple 3D curve
* **KRTZ 2023 Figure 2**: Complicated parametric curve
* **MGGJ 2023 Section 4.1**: High-degree algebraic curve
* **Projected Curve Example**: Projection of higher-dimensional curves

Each example shows how to set up the system, point, and run `track_curve` to generate corresponding figures.

## Notes

* Tracking is designed for **codimension-1** curves (systems where the Jacobian drops rank by 1).
* Interval radius `r` is automatically adjusted during tracking.
* Certification of steps is performed by Krawczyk operator-based method.
* Works both in low and moderately high precision (default is standard `ComplexField`).

## Requirements

* Julia 1.10+
* Packages:
  - `Nemo.jl`
  - `AbstractAlgebra.jl`
  - `MultivariatePolynomials.jl`
  - `LinearAlgebra`

(Install missing packages via `import Pkg; Pkg.add("PackageName")`.)

## Contact

For questions or issues, feel free to open an issue on GitHub or contact Kisun Lee (kisunl@clemson.edu).

---

Please check the example outputs in the `certified_curve_projections` directory after running the examples.
