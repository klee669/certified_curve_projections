include("certified_curve_tracking.jl")

# Figure 3 example
CCi = ComplexField()
R, (x,y,z) = CCi["x","y","z"]
F =[-z-x^3+2.7*x y^2-2+z]
p = [CCi(2.3947),CCi(3.17),CCi(-8.04)]

x,r,A=track_curve(F,p,.1,600,"~/Documents/GitHub/certified_curve_projections/small_example";
                    show_tubular_neighborhood=true, 
                    box_thickness=0.2,
                    line_thickness=0.1
                )



# Figure 6 example (Katsamaki, Rouillier, Tsigaridas, Zafeirakopoulos. 2023 Figure 2)
CCi = ComplexField()
R, (x,y,t) = CCi["x","y","t"]
F =[x-t^8+8*t^6-20*t^4+16*t^2-2 y-t^7+7*t^5-14*t^3+7*t]
p = [CCi(3),CCi(-3), CCi(-2.3)]

x,r,A=track_curve(F,p,.1,5700,"~/Documents/GitHub/certified_curve_projections/KRTZ";)
 


# example (Martin, Goldztejn, Granvilliers, Jermann. 2023 Section 4.1)
CCi = ComplexField()
R, (x,y) = CCi["x","y"]
e = .99
F = hcat([x^8-(1-e)*x^6+4*x^6*y^2-(3+15*e)*x^4*y^2+6*x^4*y^4-(3-15*e)*x^2*y^4+4*x^2*y^6-(1+e)*y^6+y^8])
p = [CCi(sqrt(1-e)),CCi(0)]

x,r,A=track_curve(F,p,.1,2500,"~/Documents/GitHub/certified_curve_projections/MGGJ"; figure_scale=3)


# Figure 5 example
CCi = ComplexField()
R, (x,y,z) = CCi["x","y","z"]
F = [x+z^5-1.3*z^3 y-z^3+z]
p = [CCi(1.81),CCi(-1.34),CCi(-1.4)]

x,r,A=track_curve(F,p,.1,550,"~/Documents/GitHub/certified_curve_projections/projected_curve"; figure_scale=3)