include("certified_curve_tracking")

CCi = ComplexField()

R, (x,y,z) = CCi["x","y","z"]
F =[x^2+y^2-3-z z-x]#(x^3+y^2-3)*z+z^2+z-1]
p = [CCi(0),CCi(-1.7320508075688774),CCi(0)]#,CCi(.61803398874989479)]

x,r,A=track_curve(F,p,.1)


function interval_svd(M)
    ring = parent(M[1,1]);
    M = convert_to_double_matrix(M);
    u,s,vt = svd(M;full=true);
    [convert_to_box_matrix(u,ring),convert_to_box_matrix(hcat(s),ring),convert_to_box_matrix(vt,ring),convert_to_box_matrix(u',ring),convert_to_box_matrix(vt',ring)]
end

function unitary_transformation2(system, p)
    ring = parent(system[1]);
    gensR = gens(ring);
    n = size(system)[2];
#    println(k_vec);
    M = evaluate_matrix(jac_nsquare(system),p);
#    M = convert_to_double_matrix(M);

    u,s,vt,ut,v = interval_svd(M);
    v_transform = vt*gensR;

    transformed_system = zeros(ring, 1, n);
    for i = 1:n
        transformed_system[i] = AbstractAlgebra.evaluate(system[i], v_transform);
    end
    [vec(ut*transpose(transformed_system)),v]
end

tF, v= unitary_transformation2(F,p)
transpose(kernel(transpose(evaluate_matrix(jac_nsquare(tF),v*p))))

ring = CCi
M = evaluate_matrix(jac_nsquare(F),p)
M = convert_to_double_matrix(M);
u,s,vt = svd(M;full=true)
u,s,vt =     [convert_to_box_matrix(u,ring),convert_to_box_matrix(hcat(s),ring),convert_to_box_matrix(vt,ring)]
gensR = gens(R)
ngens = vt*gensR
n=3
system = F
for i = 1:2
    system[i] = AbstractAlgebra.evaluate(F[i], ngens);
end
evaluate_matrix(u*transpose(system), transpose(vt)*p)
usys = vec(Matrix(u*transpose(system)))
kernel(transpose(evaluate_matrix(jac_nsquare(usys),transpose(vt)*p)))


R, (x,y,z,w) = CCi["x","y","z","w"]
F = [x^2+y^2-3*x*z+w w*z^2+onei(CCi)*z*w-x*y+31 y^2*x^2-w^2*x^2]
p = [CCi(1), CCi(-3.44158,.210174), CCi(3.11957,-.412162), CCi(-3.44158,.210174)]
r = .1

F = [x^10+y^2-3*x*z+w w^3*z^2+z*w-x*y+31 y^3*x^2-w^2*x^2]
p = [CCi(1.3426,.263165), 1, CCi(-.475736,5.42291), 1]

x,r,A=track_curve2(F,p,.1)
evaluate_matrix(F,x)


CCi = ComplexField()

R, (x,y,t) = CCi["x","y","t"]
F =[x-t^8+8*t^6-20*t^4+16*t^2-2 y-t^7+7*t^5-14*t^3+7*t]
p = [CCi(3),CCi(-3), CCi(-2.3)]
p = [CCi(-1.772923292092279974),CCi(1.37361544115583327), CCi(-1.890422726995060203)]

x, r,xprev, v, vprev, hprev=track_curve2(F,p,.1)
x, r,A=track_curve2(F,p,.1)
evaluate_matrix(F,x)
x, r,A=track_curve2(F,x,.1)
x,r,A=track_curve(F,p,.1)
 


    CCi = ComplexField()

R, (x,y,t) = CCi["x","y","t"]
F =[(29*t^3+98*t^2-23*t+10)*x-37*t^3+23*t^2-87*t-44 (11*t^3-49*t^2-47*t+40)*y+61*t^3+8*t^2+29*t-95]
p = [CCi(1.2932),CCi(0), CCi(.986694)]


x,r,A=track_curve(F,p,.1)
evaluate_matrix(F,x)
x,r,A=track_curve(F,x,.1)



CCi = ComplexField()

R, (x,y) = CCi["x","y"]
F = hcat([y^8-x^7-8*y^6+7*x^5+20*y^4-14*x^3-16*y^2+7*x+2])
p = [CCi(2),CCi(0)]

evaluate_matrix(F,p)
x,r,A=track_curve(F,p,.1)
evaluate_matrix(F,x)


CCi = ComplexField()

R, (x,y) = CCi["x","y"]
e = .99
F = hcat([x^8-(1-e)*x^6+4*x^6*y^2-(3+15*e)*x^4*y^2+6*x^4*y^4-(3-15*e)*x^2*y^4+4*x^2*y^6-(1+e)*y^6+y^8])
p = [CCi(sqrt(1-e)),CCi(0)]

evaluate_matrix(F,p)
x,r,A=track_curve(F,p,.1)
evaluate_matrix(F,x)


CCi = ComplexField()

R, (x,y) = CCi["x","y"]
F =hcat([x^2+y^2-3])
p = [CCi(0),CCi(-1.7320508075688774)]



x,r,A=track_curve(F,p,.1)


CCi = ComplexField()

R, (x,y,z) = CCi["x","y","z"]
F =hcat([-z-x^3+2.7*x y^2-2+z])
p = [CCi(2.3947),CCi(3.17),CCi(-8.04)]
p = [CCi(-2),CCi(.1),CCi(sqrt(1.99))]



x,r,A=track_curve(F,p,.1)


CCi = ComplexField()

R, (x,y) = CCi["x","y"]
F =hcat([x^3-2.7*x-y^2+2])
p = [CCi(3),CCi(-4.57)]
p = [CCi(-2),CCi(.1),CCi(sqrt(1.99))]



x,r,A=track_curve(F,p,.1)


CCi = ComplexField()

R, (x,y,z) = CCi["x","y","z"]
F =hcat([x+z^5-1.3*z^3 y-z^3+z])
p = [CCi(1.81),CCi(-1.34),CCi(-1.4)]
x,r,A=track_curve(F,p,.1)


p = [CCi(-2),CCi(.1),CCi(sqrt(1.99))]
