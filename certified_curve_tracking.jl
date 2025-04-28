# ========== Import Packages ==========

using Nemo
using AbstractAlgebra
using LinearAlgebra
using MultivariatePolynomials

# ========== System and Interval Utilities ==========
####### functions for systems and intervals
function jac_nsquare(system)
    ring = parent(system[1]);
    var = gens(ring);
    mat = zero_matrix(ring,length(system),length(var));
    for j in 1:length(system)
        for i in 1:length(var)
            d = derivative(system[j],var[i]);
            if d == 0
                mat[j,i] = zero(ring);
            else
                mat[j,i] = d;
            end
        end
    end
    Matrix(mat)
end

function jac(system)
    ring = parent(system[1]);
    var = gens(ring);
    mat = zero_matrix(ring,length(system),length(system));
    for j in 1:length(system)
        for i in 1:length(system)
            d = derivative(system[j],var[i]);
            if d == 0
                mat[j,i] = zero(ring);
            else
                mat[j,i] = d;
            end
        end
    end
    Matrix(mat)
end



function convert_to_double_int(interval)
    Base.convert(Float64,real(interval)) + Base.convert(Float64,imag(interval))*im
end

function convert_to_double_matrix(M)
    sz = size(M);
    nr = sz[1];
    nc = sz[2];
    result = zeros(Complex{Float64},nr,nc);
    for j in 1:nr
        for i in 1:nc
            result[j,i] = convert_to_double_int(M[j,i]);
        end
    end
    result
end

function convert_to_box_int(vec,ring)
    n = length(vec);
    result = [];
    for i in 1:n
        push!(result,ring(real(vec[i]),imag(vec[i])));
    end
    result
end

function max_int_norm(interval)
    r = real(interval);
    i = imag(interval);
    rmax = convert(Float64, abs(Nemo.midpoint(r)) + Nemo.radius(r));
    imax = convert(Float64, abs(Nemo.midpoint(r)) + Nemo.radius(r));
    maximum([rmax, imax])
end

function max_norm(intvec)
    sz = size(intvec);
    nr = sz[1];
    nc = sz[2];
    abs_list = [];
    for i in 1:nr
        for j in 1:nc
            push!(abs_list,max_int_norm(intvec[i,j]));
        end
    end 
    maximum(abs_list)
end


function evaluate_matrix(m, vec)
    ring = base_ring(m[1]);
    sz = size(m);
    nr = sz[1];
    nc = sz[2];
    mat = zero_matrix(ring,nr,nc);
#    println(m)
    for j in 1:nr
        for i in 1:nc
            mat[j,i] = AbstractAlgebra.evaluate(m[j,i],vec);
        end
    end
    mat#matrix(Matrix(mat))
end


function jacobian_inverse(system,vec)
    j = jac(system);
    eval_jac = evaluate_matrix(j, vec);
    Jinv, factor = pseudo_inv(eval_jac);
    1/factor * Jinv
end

function midpoint_complex_int(int)
    ring = parent(int);
    ring(midpoint(real(int)),midpoint(imag(int)));
end

function midpoint_complex_box(int)
    n = length(int);
    ring = parent(int[1]);
    result = zeros(ring, n);
    for i in 1:n
        result[i]=midpoint_complex_int(int[i]);
    end
    result
end


function convert_to_box_int(vec,ring)
    n = length(vec);
    result = [];
    for i in 1:n
        push!(result,ring(real(vec[i]),imag(vec[i])));
    end
    result
end

function convert_to_box_matrix(M,ring)
    sz = size(M);
    nr = sz[1];
    nc = sz[2];
    result = zero_matrix(ring,nr,nc);
    for j in 1:nr
        for i in 1:nc
            result[j,i] = ring(real(M[j,i]),imag(M[j,i]));
        end
    end
    Matrix(result)
end








####### functions for Krawczyk test

function krawczyk_operator(system, p, r, A)
    n = length(system);
    j = jac(system);
    CC = base_ring(system[1]);
    mat = zeros(CC,n,1);
    B = CC("+/- 1", "+/-1");
    for i in 1:n
        mat[i,1] = B;
    end
    id = identity_matrix(CC,n);
    eval_sys = evaluate_matrix(system, p);
    eval_jac = evaluate_matrix(j, vec(p+r*mat));
    K = zeros(CC,n,1);
    oper = (-1/r) * A * transpose(eval_sys) + (id - A* eval_jac)*matrix(mat);
    for i in 1:n
        K[i,1] = oper[i];
    end
    K
end



function krawczyk_test(system, point, r, A, rho)
    K = krawczyk_operator(system,point,r, A);
    max_norm(K) < rho
end

function refine_moore_box(f, x, r, A, rho)
    y = x;
    CR = parent(x[1]);
    n = size(A)[1];
    while krawczyk_test(f, y, r, A, rho) == false 
        d = A * transpose(evaluate_matrix(f, y));
        if max_norm(d) <= (1/64)*rho*r
            r = (1/2)*r;
        else
            y = midpoint_complex_box(y-d[:,1]);
        end
        A = jacobian_inverse(f, y);
    end
    while 2*r <= 1 && krawczyk_test(f, x, 2*r, A, rho)
        r = 2*r;
    end

    [y, r, A]
end




##### predictor parts
function speed_vector(H, x, A)
    ring = parent(H[1]);
    d_var = gens(ring)[end];
    n = length(H);

    result = zeros(ring, 1, n);
    for i = 1:n
        result[i]=derivative(H[i],d_var);
    end
    result = evaluate_matrix(result, x);

    midpoint_complex_box(-A*transpose(result))
end


function linear_predictor(H, v, x)

    eR = parent(H[1]);
    genseR = gens(eR);
    η =  genseR[end];

    n = size(v)[1];
    result = zeros(eR,n);
    for i = 1:n
        result[i] = v[i]*η;
    end
    x+ result
end

function taylor_model(H, lp, tval,A,r)

    n = length(H);
    eR = parent(H[1]);
    CC = coefficient_ring(eR);
    η = gens(eR)[end];
    lp = push!(lp,η);

    elim_var = gens(eR);
    elim_var[end] = tval+η;

    eH = zeros(eR,1,n);
    for i = 1:n
        eH[i] = AbstractAlgebra.evaluate(H[i], elim_var);
    end
    ejac = jac(eH);
    eHjac = zeros(eR, n,n);

    mat = zeros(CC,n+1);
    B = CC("+/- 1", "+/-1");
    for i in 1:n
        mat[i] = B;
    end



    for i = 1:n
        eH[i] = AbstractAlgebra.evaluate(eH[i],lp);
        for j = 1:n
            eHjac[i,j] = AbstractAlgebra.evaluate(ejac[i,j], lp+r*mat);
        end
    end


    id = identity_matrix(CC,n);
    (-1/r) * A * matrix(transpose(eH)) + (id - A* matrix(eHjac))*matrix(mat[1:n])
end


function proceeding_step(h, CCi, n, tm, K)

    while max_norm(K) > 7/8
        h = 1/2 * h;
        radii = h/2;
        if abs(h) < 10^(-10)
            error("h is too small!")
        end

        T = CCi("$radii +/- $radii");
        input = zeros(CCi, n+1);
        input[n+1] = T;
        K = evaluate_matrix(tm, input);
    end

    h
end

function refine_step(H, Ft, x, r, A, h)

    n = size(x)[1];
    trunc_x = x[1:n-1];
    trunc_x,r,A = refine_moore_box(Ft, trunc_x, r, A, 1/8);
    x = push!(trunc_x, x[end]);
    v = speed_vector(H,x,A);
    h = (5/4)*h;
    radii = h/2;

    [x, r, A, v, h, radii]
 
end


function hermite_predictor(H, x, xprev, v, vprev, hprev)

    eR = parent(H[1]);
    genseR = gens(eR);
    η =  genseR[end];

    n = size(v)[1];
    result = [];
    for i = 1:n
        result = push!(result, v[i]*η+ (3*v[i]/hprev-(v[i]-vprev[i])/hprev-3*(x[i]-xprev[i])/hprev^2)*(η^2)+
        (2*v[i]/hprev^2-(v[i]-vprev[i])/hprev^2-2*(x[i]-xprev[i])/hprev^3)*(η^3));
    end
    x+ result 
end

# ========== Transformation Utilities ==========
function interval_svd(M)
    ring = parent(M[1,1])
    M = convert_to_double_matrix(M)
    u, s, vt = svd(M; full=true)
    [convert_to_box_matrix(u, ring), convert_to_box_matrix(hcat(s), ring), convert_to_box_matrix(vt, ring)]
end

function origin_transformation(system, p)
    ring = parent(system[1])
    gensR = gens(ring)
    n = size(system)[2]

    transformed_system = zeros(ring, 1, n)
    for i in 1:n
        transformed_system[i] = AbstractAlgebra.evaluate(system[i], gensR + p)
    end
    transformed_system
end

function unitary_transformation(system, k_vec)
    ring = parent(system[1])
    gensR = gens(ring)
    n = size(system)[2]

    u, s, vt = interval_svd(transpose(kernel(k_vec)))
    u_transform = u * gensR

    transformed_system = zeros(ring, 1, n)
    for i in 1:n
        transformed_system[i] = AbstractAlgebra.evaluate(system[i], u_transform)
    end
    [u, transformed_system * vt]
end

function system_transform(F, x)
    ring = parent(F[1])
    CCi = base_ring(F[1])
    n = length(x)

    k_prev = kernel(transpose(evaluate_matrix(jac_nsquare(F), x)))
    t_F = origin_transformation(F, x)
    J = transpose(evaluate_matrix(jac_nsquare(t_F), zeros(CCi, n)))
    k = matrix(midpoint_complex_box(transpose(kernel(J))))

    u, u_F = unitary_transformation(t_F, k)

    elim_var = gens(ring)
    elim_var[end] = ring(0)
    eval_uF = zeros(ring, 1, n-1)
    for i in 1:n-1
        eval_uF[i] = AbstractAlgebra.evaluate(u_F[i], elim_var)
    end

    [u, u_F, eval_uF]
end

function system_transform_with_prev(F, x, xprev)
    ring = parent(F[1])
    CCi = base_ring(F[1])
    n = length(x)

    k_prev = kernel(transpose(evaluate_matrix(jac_nsquare(F), x)))
    t_F = origin_transformation(F, x)

    J = transpose(evaluate_matrix(jac_nsquare(t_F), zeros(CCi, n)))
    k = matrix(midpoint_complex_box(transpose(kernel(J))))

    u, u_F = unitary_transformation(t_F, k)

    elim_var = gens(ring)
    elim_var[end] = ring(0)
    eval_uF = zeros(ring, 1, n-1)
    for i in 1:n-1
        eval_uF[i] = AbstractAlgebra.evaluate(u_F[i], elim_var)
    end

    [u, u_F, eval_uF]
end



function safe_path(file_name::String)
    if startswith(file_name, "~")
        return expanduser(file_name)
    elseif isabspath(file_name)
        return file_name
    else
        return joinpath(pwd(), file_name)
    end
end

function pretty_print_status(x, iter, total_iter)
    output = "[Iteration $iter/$total_iter] Tracking point: x = ["
    coords = []
    for xi in x
        mid_real = round(Float64(midpoint(real(xi))), digits=6)
        err_real = round(Float64(radius(real(xi))), sigdigits=2)
        push!(coords, "$(mid_real) ± $(err_real)")
    end
    output *= join(coords, ", ") * "]"
    print("\r", output)
    flush(stdout)
end

# ========== Main Tracking Function ==========
function track_curve(F, x, r, max_iter, file_name; show_tubular_neighborhood=false, figure_scale=1, box_thickness=0.1, line_thickness=0.2)

    # Initialization
    ring = parent(F[1])
    CCi = base_ring(F[1])
    n = length(x)
    h = 1/2
    iter = 0
    init_x = x

    # Step 1: First transformation
    u, u_F, eval_uF = system_transform(F, x)
    A = jacobian_inverse(eval_uF, zeros(CCi, n))
    u_x, r, A, v, h, radii = refine_step(u_F, eval_uF, zeros(CCi, n), r, A, h)

    X = linear_predictor(u_F, v, u_x[1:n-1])
    t = u_x[end]
    tm = taylor_model(u_F, X, t, A, r)

    T = CCi("$(radii) +/- $(radii)")
    input = zeros(CCi, n+1)
    input[end] = T
    K = evaluate_matrix(tm, input)

    h = proceeding_step(h, CCi, n, tm, K)
    u_x[end] += h

    xprev = x
    x = transpose(transpose(u_x) * Matrix(inv(matrix(u)))) + x

    hprev = h
    vprev = v
    dx_prev = zeros(CCi, n)
    dx = x - xprev

    iter = 0
    rev_count = 0
    prev_count = 0

    path = safe_path("$file_name.tex")
    open(path, "w") do file
        write(file, "\\begin{tikzpicture}[scale=$figure_scale]\n")

        # Main tracking loop
        while iter < max_iter
            u, u_F, eval_uF = system_transform_with_prev(F, x, xprev)
            A = jacobian_inverse(eval_uF, zeros(CCi, n))
            u_x, r, A, v, h, radii = refine_step(u_F, eval_uF, zeros(CCi, n), r, A, h)

            xprev = transpose(transpose(x - xprev) * u)
            X = hermite_predictor(u_F, u_x[1:n-1], xprev[1:n-1], v, vprev, hprev)
            t = u_x[end]
            tm = taylor_model(u_F, X, t, A, r)

            T = CCi("$(radii) +/- $(radii)")
            input = zeros(CCi, n)
            input[end] = T
            K = evaluate_matrix(tm, input)

            b = vec(Matrix(evaluate_matrix(hcat(X), input))) * inv(matrix(u)) + x


            h = proceeding_step(h, CCi, n, tm, K)
            u_x[end] += h

            prev_count = rev_count
            if real(transpose(xprev) * u_x) < 0
                rev_count += 1
                u_x[end] -= 2*h
                radii = -radii
                T = CCi("$(radii) +/- $(radii)")
                input = zeros(CCi, n)
                input[end] = T
                b = vec(Matrix(evaluate_matrix(hcat(X), input))) * inv(matrix(u)) + x
            end

            input = zeros(CCi, n)
            input[n] = CCi(h)

            xprev = x
            x = midpoint_complex_box(transpose(transpose(u_x) * Matrix(inv(matrix(u))))) + x

            pretty_print_status(x, iter, max_iter)

            # Draw output
            if show_tubular_neighborhood && iter > 5
                if prev_count != rev_count
                    write(file,"\\draw[color=blue,line width=$(box_thickness)mm] ($(real(convert_to_double_int(b[1]))-convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))-convert(Float64,radius(real(b[2]))))) --($(real(convert_to_double_int(b[1]))+convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))-convert(Float64,radius(real(b[2])))));\n")
                    write(file,"\\draw[color=blue,line width=$(box_thickness)mm] ($(real(convert_to_double_int(b[1]))-convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))-convert(Float64,radius(real(b[2]))))) --($(real(convert_to_double_int(b[1]))-convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))+convert(Float64,radius(real(b[2])))));\n")
                    write(file,"\\draw[color=blue,line width=$(box_thickness)mm] ($(real(convert_to_double_int(b[1]))-convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))+convert(Float64,radius(real(b[2]))))) --($(real(convert_to_double_int(b[1]))+convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))+convert(Float64,radius(real(b[2])))));\n")
                    write(file,"\\draw[color=blue,line width=$(box_thickness)mm] ($(real(convert_to_double_int(b[1]))+convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))-convert(Float64,radius(real(b[2]))))) --($(real(convert_to_double_int(b[1]))+convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))+convert(Float64,radius(real(b[2])))));\n")
                else
                    write(file,"\\draw[color=blue,line width=$(box_thickness)mm] ($(real(convert_to_double_int(b[1]))-convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))-convert(Float64,radius(real(b[2]))))) --($(real(convert_to_double_int(b[1]))+convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))-convert(Float64,radius(real(b[2])))));\n")
                    write(file,"\\draw[color=blue,line width=$(box_thickness)mm] ($(real(convert_to_double_int(b[1]))-convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))-convert(Float64,radius(real(b[2]))))) --($(real(convert_to_double_int(b[1]))-convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))+convert(Float64,radius(real(b[2])))));\n")
                    write(file,"\\draw[color=blue,line width=$(box_thickness)mm] ($(real(convert_to_double_int(b[1]))-convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))+convert(Float64,radius(real(b[2]))))) --($(real(convert_to_double_int(b[1]))+convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))+convert(Float64,radius(real(b[2])))));\n")
                    write(file,"\\draw[color=blue,line width=$(box_thickness)mm] ($(real(convert_to_double_int(b[1]))+convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))-convert(Float64,radius(real(b[2]))))) --($(real(convert_to_double_int(b[1]))+convert(Float64,radius(real(b[1])))),$(real(convert_to_double_int(b[2]))+convert(Float64,radius(real(b[2])))));\n")
                end
            end

            write(file,"\\draw[color=red,line width=$(line_thickness)mm] ($(real(convert_to_double_int(x[1]))),$(real(convert_to_double_int(x[2])))) --($(real(convert_to_double_int(xprev[1]))),$(real(convert_to_double_int(xprev[2]))));\n")

            dx_prev = dx
            dx = x - xprev
            hprev = h
            vprev = v

            iter += 1
        end

        write(file,"\\end{tikzpicture}");
    end

    [x, r, A]
end
