### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 6f222ff3-57b0-4048-bbec-6f1d5742654e
begin
	import Pkg
	Pkg.activate()
	using JuMP, AmplNLWriter, Couenne_jll, Plots, Clipper
end

# ╔═╡ 86a7842c-855b-46fb-b4e8-afccada2996c
include("Clipper_utilisation.jl")

# ╔═╡ aa8dbdcb-2f02-4733-a783-0a9da14b4090
include("fonction.jl")

include("polygon_subdivision.jl")
include("test_functions.jl")
include("basic_planar_and_polygonal_functions.jl")

# ╔═╡ 0db861fd-2d46-49c6-9ae4-3dfe9c0c89ae
#Creation structure Triangle
#Convention Premier élément du tableau coordonnées x puis y
struct Triangle
	coord_S1::Array{Float64}
	coord_S2::Array{Float64}
	coord_S3::Array{Float64}
end


# ╔═╡ 9fe67875-c0d5-42a1-a887-2e21f39df090
#=triangle =  Any[Main.workspace2.Triangle([1.3750000020262492, 2.000000003514397, -2.1093750084854035], [0.5, 0.5, 0.0], [0.5, 2.0, -3.75]), Main.workspace2.Triangle([2.2499999999523683, 3.5, -7.187500000214342], [1.3750000020262492, 2.000000003514397, -2.1093750084854035], [0.5, 2.0, -3.75]), Main.workspace2.Triangle([2.2499999999523683, 3.5, -7.187500000214342], [0.5, 2.0, -3.75], [0.5, 3.5, -12.0]), Main.workspace2.Triangle([3.999999999999999, 3.5, 3.749999999999993], [0.5, 0.5, 0.0], [1.3749999969217472, 1.9999999947638223, -2.4427083208538174]), Main.workspace2.Triangle([3.999999999999999, 3.5, 3.749999999999993], [1.3749999969217472, 1.9999999947638223, -2.4427083208538174], [2.2499999999523683, 3.5, -7.187500000214342]), Main.workspace2.Triangle([2.2499999878828376, 1.2499999948069305, 3.499999958455443], [0.5, 0.5, 0.0], [3.999999999999999, 3.5, 3.749999999999993]), Main.workspace2.Triangle([3.9999999935373927, 1.9999999972303113, 11.999999959377895], [2.2499999878828376, 1.2499999948069305, 3.499999958455443], [3.999999999999999, 3.5, 3.749999999999993]), Main.workspace2.Triangle([5.749999999999994, 3.5, 20.81249999999993], [3.9999999935373927, 1.9999999972303113, 11.999999959377895], [3.999999999999999, 3.5, 3.749999999999993]), Main.workspace2.Triangle([5.749999994866336, 2.749999997799858, 25.499999953063647], [3.9999999935373927, 1.9999999972303113, 11.999999959377895], [5.749999999999994, 3.5, 20.81249999999993]), Main.workspace2.Triangle([7.5, 3.5, 44.0], [5.749999994866336, 2.749999997799858, 25.499999953063647], [5.749999999999994, 3.5, 20.81249999999993]), Main.workspace2.Triangle([6.62499999566402, 1.9999999925668905, 39.89062497228071], [7.5, 3.5, 44.0], [7.5, 1.9999999999999991, 52.25]), Main.workspace2.Triangle([5.750000000000001, 0.5, 32.81250000000001], [6.62499999566402, 1.9999999925668905, 39.89062497228071], [7.5, 1.9999999999999991, 52.25]), Main.workspace2.Triangle([5.750000000000001, 0.5, 32.81250000000001], [7.5, 1.9999999999999991, 52.25], [7.5, 0.5, 56.0]), Main.workspace2.Triangle([4.0, 0.5, 15.75], [7.5, 3.5, 44.0], [6.625000008854893, 2.000000015179816, 39.55729172327473]), Main.workspace2.Triangle([4.0, 0.5, 15.75], [6.625000008854893, 2.000000015179816, 39.55729172327473], [5.750000000000001, 0.5, 32.81250000000001]), Main.workspace2.Triangle([5.750000013302918, 2.7500000057012506, 25.500000121626677], [7.5, 3.5, 44.0], [4.0, 0.5, 15.75]), Main.workspace2.Triangle([4.000000008870701, 2.0000000038017287, 12.00000005575869], [5.750000013302918, 2.7500000057012506, 25.500000121626677], [4.0, 0.5, 15.75]), Main.workspace2.Triangle([2.2499999999999933, 0.5, 4.81249999999997], [4.000000008870701, 2.0000000038017287, 12.00000005575869], [4.0, 0.5, 15.75]), Main.workspace2.Triangle([2.250000006607323, 1.2500000028317098, 3.5000000226536794], [4.000000008870701, 2.0000000038017287, 12.00000005575869], [2.2499999999999933, 0.5, 4.81249999999997]), Main.workspace2.Triangle([0.5, 0.5, 0.0], [2.250000006607323, 1.2500000028317098, 3.5000000226536794], [2.2499999999999933, 0.5, 4.81249999999997])]=#


# ╔═╡ cbc4c2ec-64c4-4012-8283-a237c55eb6b7
function passtriangle(triangle)
	poly = []
	for i in triangle
		push!(poly, [i.coord_S1, i.coord_S2, i.coord_S3])
	end
	return poly
end

# ╔═╡ 06d9059a-e9ff-4416-a070-be7c24f83f3d
function passPolytope(p)
	poly = []
	for i in p
		coor = []
		for j in 1:length(i)
			if length(i[j])==2
				push!(coor, [i[j][1],i[j][2],0.0])
			else
				push!(coor, i[j])
			end
		end
		push!(poly,coor)
	end
	return poly
end

# ╔═╡ e78c2699-4df0-4bdf-9b74-1ce09bdc3c2f
#test = passtriangle(triangle)

# ╔═╡ fa6ae86f-fbb7-43fa-af56-e0d7e2240d26
test = Any[[[0.5, 0.5, 0.3317653298563168], [0.8795712336832989, 0.5, 0.2917685103559552], [1.5877216839915005, 0.951581759990485, 0.06950073395567472], [0.5, 1.2544314124496994, 0.0850994726676645]], [[0.5, 2.0, 0.026504606781503653], [0.5, 1.2544314124496994, 0.0850994726676645], [1.5877216839915005, 0.951581759990485, 0.06950073395567472], [2.0, 1.2144882757852618, 0.033904956035722134], [2.0, 2.0, -0.027829074744194993]], [[2.0, 0.5, 0.031992189819741335], [2.0, 0.8367928883910293, -0.011324280200369269], [1.5877216839915005, 0.951581759990485, 0.06950073395567472], [0.8795712336832989, 0.5, 0.2917685103559552]], [[2.0, 0.8367928883910293, -0.011324280200369269], [2.0, 1.2144882757852618, 0.033904956035722134], [1.5877216839915005, 0.951581759990485, 0.06950073395567472]]]


function createDomain(listeDomain)
	return [[listeDomain[1],listeDomain[2]], [listeDomain[3],listeDomain[4]]]
end

# ╔═╡ 5d56e2eb-fbfb-471b-ada2-ae1f5c208a34
#On supprime un élément particulier du tableau
function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end


# ╔═╡ 4293c8df-e0c8-42e3-b61c-696a5055de12
#Calcul des coefficients directeurs des côtés du triangle
function slopeCoefficientCheck(poly)
	slopecoefficient = []
	for i in 1:length(poly)
		if i == length(poly)
			push!(slopecoefficient, (poly[1][2] - poly[i][2])/(poly[1][1] - poly[i][1]))
		else
			push!(slopecoefficient, (poly[i+1][2] - poly[i][2])/(poly[i+1][1] - poly[i][1]))
		end
	end
	return slopecoefficient
end

# ╔═╡ 354996bb-d6dd-41a8-bc4c-e023792ebdce
#Calcul des ordonnées à l'origine
function interceptCheck(poly, slopecoefficient)
	intercept = []
	for i in 1:length(poly)
		calcul = poly[i][2] - (slopecoefficient[i]*poly[i][1])
		push!(intercept, calcul)
	end
	return intercept
end

# ╔═╡ 618b0169-e701-4aca-945e-798494015bf0
#Configuration du modèle pour le solveur non linéaire
function configModelNL(solver, analysedfunction, domain, delta, poly, slopecoefficient, intercept, bigDelta, valx, valy)

	#Creation du modèle
	if isequal(solver, "couenne")
		model = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe))
	else
		model = Model(() -> AmplNLWriter.Optimizer(Ipopt.amplexe))
	end

	@variable(model, 0.0 <= lambda[1:(length(poly))] <= 1.0)
	@expression(model, aux, sum(lambda[i] * poly[i] for i in 1:(length(poly))))
	@variable(model, Pt[1:2])

	#Contraindre à rester au sein du triangle
	@constraint(model, Pt[1] == aux[1])
	@constraint(model, Pt[2] == aux[2])
	@constraint(model, sum(lambda[i] for i in 1:length(poly)) == 1.0)

	#Définition de plusieurs expressions de la fonction étudiée afin de bien bâtir le modèle et notre approximation linéaire l de la fonction
	@expression(model, fonctcoord[i = 1:(length(poly))], poly[i][3])

	#A chaque itération de la boucle, on ne sauvegarde pas le l dans le modèle
	unregister(model, :l)


	#L'optimisation est activée ou non
	@NLexpression(model, l, sum(lambda[i] * (fonctcoord[i]) for i in 1:(length(poly))))

	println("---------- avant la sélection de la fonction objectif ------------")

	if analysedfunction == f1
		@NLobjective(model, Max, abs(l- Pt[1]^2 + Pt[2]^2))
	elseif analysedfunction == f2
		@NLobjective(model, Max, abs(l- Pt[1]^2 - Pt[2]^2))
	elseif analysedfunction == f3
		@NLobjective(model, Max, abs(l- Pt[1] * Pt[2]))
	elseif analysedfunction == f4
		@NLobjective(model, Max, abs(l- Pt[1] * exp(-Pt[1]^2-Pt[2]^2)))
	elseif analysedfunction == f5
		@NLobjective(model, Max, abs(l- Pt[1] * sin(Pt[2])))
	elseif analysedfunction == f6
		@NLobjective(model, Max, abs(l- sin(Pt[1])/Pt[1] * Pt[2]^2))
	elseif analysedfunction == f7
		@NLobjective(model, Max, abs(l- Pt[1] * sin(Pt[1]) * sin(Pt[2])))
	elseif analysedfunction == f8
		@NLobjective(model, Max, abs(l- (Pt[1]^2 - Pt[2]^2)^2))
	elseif analysedfunction == f9
		@NLobjective(model, Max, abs(l- exp(-10 * (Pt[1]^2 - Pt[2]^2)^2)))
	elseif analysedfunction == f0
		@NLobjective(model, Max, abs(l- sqrt(Pt[1]^2 + Pt[2]^2)))
	else
		println("---------- sélection de la fonction objectif ratée, on est dans le else ------------")
	end

	println("---------- après la sélection de la fonction objectif ------------")

	#println(model)
	optimize!(model)
	#println(JuMP.optimize!(model))
	term_status = JuMP.termination_status(model)
	println("statut de terminaison de l'optimisation du modèle : $term_status")

	#On regarde si la solution est faisable où non pour retourner le résultat du solveur ou pas
	if(term_status != MOI.INFEASIBLE && term_status != MOI.OTHER_ERROR )
		bigDelta = objective_value(model)
		valx = JuMP.value.(Pt[1])
		valy = JuMP.value.(Pt[2])
	end

	return bigDelta, valx, valy
end

function newConfigModelNL(solver, analysedfunction, domain, delta, poly, slopecoefficient, intercept, bigDelta, valx, valy)
	# résolution d'un modèle plus petit en remplaçant l'utilisation de combinaisons convexes des sommets par une variable
	# donnant directement la position limitée par des contraintes linéaires

	# calculer [a,b,c] donnant la fonction linéaire z = ax+by+c
	a,b,c = get_planar_equation(poly)
	println("équation du plan poly : $poly\na,b,c = $a,$b,$c")

	#Creation du modèle
	if isequal(solver, "couenne")
		model = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe))
	else
		model = Model(() -> AmplNLWriter.Optimizer(Ipopt.amplexe))
	end

	@variable(model, Pt[1:2])

	#Contraindre Pt à rester au sein du polygone
	n = length(poly)
	for i in 1:n
		v1 = poly[i][1:2]
		v2 = poly[mod(i,n)+1][1:2]
		edge = [v2[1]-v1[1], v2[2]-v1[2]]
		normal_vec = [-edge[2], edge[1]]
		normal_vec = normal_vec ./ norm(normal_vec)
		@constraint(model, normal_vec[1]*Pt[1]+normal_vec[2]*Pt[2] >= normal_vec[1]*v1[1]+normal_vec[2]*v1[2])
	end

	#A chaque itération de la boucle, on ne sauvegarde pas le l dans le modèle
	unregister(model, :l)

	#L'optimisation est activée ou non
	@NLexpression(model, l, a*Pt[1] + b*Pt[2] + c)

	println("---------- avant la sélection de la fonction objectif ------------")

	if analysedfunction == f1
		@NLobjective(model, Max, abs(l- Pt[1]^2 + Pt[2]^2))
	elseif analysedfunction == f2
		@NLobjective(model, Max, abs(l- Pt[1]^2 - Pt[2]^2))
	elseif analysedfunction == f3
		@NLobjective(model, Max, abs(l- Pt[1] * Pt[2]))
	elseif analysedfunction == f4
		@NLobjective(model, Max, abs(l- Pt[1] * exp(-Pt[1]^2-Pt[2]^2)))
	elseif analysedfunction == f5
		@NLobjective(model, Max, abs(l- Pt[1] * sin(Pt[2])))
	elseif analysedfunction == f6
		@NLobjective(model, Max, abs(l- sin(Pt[1])/Pt[1] * Pt[2]^2))
	elseif analysedfunction == f7
		@NLobjective(model, Max, abs(l- Pt[1] * sin(Pt[1]) * sin(Pt[2])))
	elseif analysedfunction == f8
		@NLobjective(model, Max, abs(l- (Pt[1]^2 - Pt[2]^2)^2))
	elseif analysedfunction == f9
		@NLobjective(model, Max, abs(l- exp(-10 * (Pt[1]^2 - Pt[2]^2)^2)))
	elseif analysedfunction == f0
		@NLobjective(model, Max, abs(l- sqrt(Pt[1]^2 + Pt[2]^2)))
	else
		println("---------- sélection de la fonction objectif ratée, on est dans le else ------------")
	end

	println("---------- après la sélection de la fonction objectif ------------")

	#println(model)
	optimize!(model)
	#println(JuMP.optimize!(model))
	term_status = JuMP.termination_status(model)
	println("statut de terminaison de l'optimisation du modèle : $term_status")

	#On regarde si la solution est faisable où non pour retourner le résultat du solveur ou pas
	if(term_status != MOI.INFEASIBLE && term_status != MOI.OTHER_ERROR )
		bigDelta = objective_value(model)
		valx = JuMP.value.(Pt[1])
		valy = JuMP.value.(Pt[2])
	end

	return bigDelta, valx, valy
end

function betterConfigModelNL(solver, analysedfunction, domain, delta, poly, slopecoefficient, intercept, bigDelta, valx, valy)
	# résolution d'un modèle plus petit en remplaçant l'utilisation de combinaisons convexes des sommets par une variable
	# donnant directement la position limitée par des contraintes linéaires

	# calculer [a,b,c] donnant la fonction linéaire z = ax+by+c
	a,b,c = get_planar_equation(poly)
	println("équation du plan poly : $poly\na,b,c = $a,$b,$c")

	#Creation du modèle
	if isequal(solver, "couenne")
		model = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe))
	else
		model = Model(() -> AmplNLWriter.Optimizer(Ipopt.amplexe))
	end

	@variable(model, Pt[1:2])

	#Contraindre Pt à rester au sein du polygone
	n = length(poly)
	for i in 1:n
		v1 = poly[i][1:2]
		v2 = poly[mod(i,n)+1][1:2]
		edge = [v2[1]-v1[1], v2[2]-v1[2]]
		normal_vec = [-edge[2], edge[1]]
		normal_vec = normal_vec ./ norm(normal_vec)
		@constraint(model, normal_vec[1]*Pt[1]+normal_vec[2]*Pt[2] >= normal_vec[1]*v1[1]+normal_vec[2]*v1[2])
	end

	#A chaque itération de la boucle, on ne sauvegarde pas le l dans le modèle
	unregister(model, :l)

	#L'optimisation est activée ou non
	@NLexpression(model, l, a*Pt[1] + b*Pt[2] + c)

	println("---------- avant la sélection de la fonction objectif ------------")

	eps = 2e-6

	if analysedfunction == f1
		@NLconstraint(model, abs(l - Pt[1]^2 + Pt[2]^2) >= delta+eps)
		#@NLconstraint(model, Pt[1]^2 - Pt[2]^2 - l >= delta+eps)
		#@NLobjective(model, Max, abs(l- Pt[1]^2 + Pt[2]^2))
	elseif analysedfunction == f2
		@NLconstraint(model, abs(l - Pt[1]^2 - Pt[2]^2) >= delta+eps)
		#@NLconstraint(model, Pt[1]^2 + Pt[2]^2 - l >= delta+eps)
		#@NLobjective(model, Max, abs(l- Pt[1]^2 - Pt[2]^2))
	elseif analysedfunction == f3
		@NLconstraint(model, abs(l - Pt[1] * Pt[2]) >= delta+eps)
		#@NLconstraint(model, Pt[1] * Pt[2] - l >= delta+eps)
		#@NLobjective(model, Max, abs(l- Pt[1] * Pt[2]))
	elseif analysedfunction == f4
		@NLconstraint(model, abs(l - Pt[1] * exp(-Pt[1]^2-Pt[2]^2)) >= delta+eps)
		#@NLconstraint(model, Pt[1] * exp(-Pt[1]^2-Pt[2]^2) - l >= delta+eps)
		#@NLobjective(model, Max, abs(l- Pt[1] * exp(-Pt[1]^2-Pt[2]^2)))
	elseif analysedfunction == f5
		@NLconstraint(model, abs(l - Pt[1] * sin(Pt[2])) >= delta+eps)
		#@NLconstraint(model, Pt[1] * sin(Pt[2]) - l >= delta+eps)
		#@NLobjective(model, Max, abs(l- Pt[1] * sin(Pt[2])))
	elseif analysedfunction == f6
		@NLconstraint(model, abs(l - sin(Pt[1])/Pt[1] * Pt[2]^2) >= delta+eps)
		#@NLconstraint(model, sin(Pt[1])/Pt[1] * Pt[2]^2 - l >= delta+eps)
		#@NLobjective(model, Max, abs(l- sin(Pt[1])/Pt[1] * Pt[2]^2))
	elseif analysedfunction == f7
		@NLconstraint(model, abs(l - Pt[1] * sin(Pt[1]) * sin(Pt[2])) >= delta+eps)
		#@NLconstraint(model, Pt[1] * sin(Pt[1]) * sin(Pt[2]) - l >= delta+eps)
		#@NLobjective(model, Max, abs(l- Pt[1] * sin(Pt[1]) * sin(Pt[2])))
	elseif analysedfunction == f8
		@NLconstraint(model, abs(l - (Pt[1]^2 - Pt[2]^2)^2) >= delta+eps)
		#@NLconstraint(model, (Pt[1]^2 - Pt[2]^2)^2 - l >= delta+eps)
		#@NLobjective(model, Max, abs(l- (Pt[1]^2 - Pt[2]^2)^2))
	elseif analysedfunction == f9
		@NLconstraint(model, abs(l - exp(-10 * (Pt[1]^2 - Pt[2]^2)^2)) >= delta+eps)
		#@NLconstraint(model, exp(-10 * (Pt[1]^2 - Pt[2]^2)^2) - l >= delta+eps)
		#@NLobjective(model, Max, abs(l- exp(-10 * (Pt[1]^2 - Pt[2]^2)^2)))
	elseif analysedfunction == f0
		@NLconstraint(model, abs(l - sqrt(Pt[1]^2 + Pt[2]^2)) >= delta+eps)
		#@NLconstraint(model, sqrt(Pt[1]^2 + Pt[2]^2) - l >= delta+eps)
		#@NLobjective(model, Max, abs(l- sqrt(Pt[1]^2 + Pt[2]^2)))
	else
		println("---------- sélection de la fonction objectif ratée, on est dans le else ------------")
	end

	println("---------- après la sélection de la fonction objectif ------------")

	#println(model)
	optimize!(model)
	#println(JuMP.optimize!(model))
	term_status = JuMP.termination_status(model)
	println("statut de terminaison de l'optimisation du modèle : $term_status")

	#On regarde si la solution est faisable où non pour retourner le résultat du solveur ou pas
	if(term_status != MOI.INFEASIBLE && term_status != MOI.OTHER_ERROR )
		valx = JuMP.value.(Pt[1])
		valy = JuMP.value.(Pt[2])
		valpwl = a*valx + b*valy + c
		valf = analysedfunction(valx,valy)
		bigDelta = abs(valf-valpwl)
		file = open("error_checker.txt", "a")
		println(file, "erreur checker pour polygone $poly plan  z = $a x + $b y + $c ")
		println(file, "pour (x,y) = ($valx,$valy) et bigDelta $bigDelta car f(x,y) = $valf, valpwl = $valpwl")
		close(file)
	elseif term_status == MOI.INFEASIBLE
		bigDelta = 0
		valx = 2
		valy = 2
	end

	return bigDelta, valx, valy
end

# ╔═╡ 509f5ea5-bafa-47a9-9b77-9db708dcb122
function affichage(polytrue, polyfalse, domain, fast_exit=0)
	f = plot()

	Domain = [[[domain[1][1],domain[2][1]], [domain[1][2],domain[2][1]],[domain[1][2],domain[2][2]], [domain[1][1],domain[2][2]]]]

	#=polytrue = [[[[3.0,0.5],[4.0,0.5],[4.0,3.5],[3.0,3.5]]],
				#[[[0.5,0.5],[3.0,3.5],[0.5,3.5]]],
				[[[0.5,0.5],[3.0,0.5],[3.0,3.5]]],
				[[[4.0,0.5],[7.5,0.5],[4.0,3.5]]],
	]
	polyfalse = [[[[7.5,0.5],[7.5,3.5],[4.0,3.5]]]]=#
	if fast_exit == 0
		for i in polytrue
			x = []
			y = []
			#println("i : $i")
			for j in 1:length(i[1])
				push!(x, i[1][j][1])
				push!(y, i[1][j][2])
			end
			if i==first(polytrue)
				f = plot!(x,y, color=:green, seriestype=:shape, fillalpha=0.4, label="Polytope respecting the tolerance threshold")
			else
				f = plot!(x,y, color=:green, seriestype=:shape, fillalpha=0.4, label="")
			end
			#println(length(domain))
			#println("Domaine: $Domain")

			k=1
			newDomain = []
			while k <= length(Domain)
				#println("$(i[1])")
				# old version with problems : Domain = general_polygon_difference_clipper(Domain[k], i[1])
				new_polygons = general_polygon_difference_clipper(Domain[k], i[1])
				for j in 1:length(new_polygons)
					new_polygons[j] = check_polygon(new_polygons[j])
				end
				append!(newDomain, new_polygons)
				#println("Domainefin: $Domain")
				k+=1
			end
			Domain = newDomain
		end
	end

	for i in polyfalse
		x = []
		y = []

		for j in 1:length(i[1])
			push!(x, i[1][j][1])
			push!(y, i[1][j][2])
		end

		if i==first(polyfalse)
			f = plot!(x,y, color=:red, seriestype=:shape, fillalpha=0.4, label="Polytope not respecting the tolerance threshold")
		else
			f = plot!(x,y, color=:red, seriestype=:shape, fillalpha=0.4, label="")
		end
		newDomain = []
		k=1
		while k <= length(Domain)
			# old version with problems : Domain = general_polygon_difference_clipper(Domain[k], i[1])
			new_polygons = general_polygon_difference_clipper(Domain[k], i[1])
			for j in 1:length(new_polygons)
				new_polygons[j] = check_polygon(new_polygons[j])
			end
			append!(newDomain, new_polygons)
			k+=1
		end
		Domain = newDomain
	end

	for i in Domain
		x = []
		y = []
		#println("domainei : $i")
		for j in 1:length(i)
			push!(x, i[j][1])
			push!(y, i[j][2])
		end

		if i==first(Domain)
			f = plot!(x,y, color=:black, seriestype=:shape, fillalpha=1.0, label="Polytope not covered")
		else
			f = plot!(x,y, color=:black, seriestype=:shape, fillalpha=1.0, label="")
		end
	end

	p = plot!(ylims = (domain[2][1]-0.25, domain[2][2]+0.75))

	coveredDomain = false
	if length(Domain) == 0
		coveredDomain = true
	end

	return f, coveredDomain, Domain
end

# ╔═╡ 3f8ff1ac-e579-11eb-1f8e-279ea9199229
#Algorithme permettant de checker si un ensemble de polytopes respecte nos critères de tolérance
function checker(solver, analysedfunction, domain, delta, Polygones, fast_exit = 0, check_covering = true)

	# paramètre eps pour vérifier bigDelta < delta*(1+eps) pour éviter les erreurs 0.50001 > 0.5 :
	eps = 4e-4 # plus grosse erreur de précision remarquée jusqu'à présent : 3e-4

	#Polygone = passPolytope(Polygones)
	polyfalse = []
	polytrue = []
	juste = true

	if typeof(Polygones[1])==Triangle
		Polygones = passtriangle(Polygones)
	end

	for poly in Polygones
		#Récupération des coordonées du triangles
		#=coordpoly = []
		for j in 1:length(poly)
			push!(coordpoly, poly[j])
	end=#

		#Calcul coefficient directeur
		slopecoefficient = slopeCoefficientCheck(poly)
		#println(slopecoefficient)

		#Calcul des ordonnées à l'origine
		intercept = interceptCheck(poly, slopecoefficient)
		#println(intercept)

		#=#Ajustement du delta
		if modeDecoupage == 1
			delta = delta / 2
	end=#

		#Valeurs à retourner
		bigDelta = Inf
		valx = Inf
		valy = Inf
		lambda = [0,0,0,0,0,0,0,0]

		#bigDelta, valx, valy = configModelNL(solver, analysedfunction, domain, delta, poly, slopecoefficient, intercept, bigDelta, valx, valy)
		bigDelta, valx, valy = newConfigModelNL(solver, analysedfunction, domain, delta, poly, slopecoefficient, intercept, bigDelta, valx, valy)
		#bigDelta, valx, valy = betterConfigModelNL(solver, analysedfunction, domain, delta, poly, slopecoefficient, intercept, bigDelta, valx, valy)
		#=if bigDelta == -1
			bigDelta, valx, valy = newConfigModelNL(solver, analysedfunction, domain, delta, poly, slopecoefficient, intercept, bigDelta, valx, valy)
		end=#

		#Si bigDelta est supérieur à delta, alors on va relancer le solveur en mettant toutes les variables d'écart libres à 0
		if bigDelta > delta*(1+eps) && bigDelta != Inf
			println("Pavage ne satisfaisant pas nos critères d'erreurs")
			println("INSATISFAIT")
			println("//////////////////////////////////////////////////////////////")
			juste = false
			if fast_exit == 1
				return juste, "Pavage ne satisfaisant pas nos critères d'erreurs",bigDelta, poly
			end
			push!(polyfalse, [poly, bigDelta, valx, valy])
		else

			push!(polytrue, [poly, bigDelta, valx, valy])
			println("triangle qui satisfait nos critères d'erreurs")
			println("CORRECT")
			println("//////////////////////////////////////////////////////////////")
		end

	end

	if check_covering
		f, coveredDomain, uncoveredDomain = affichage(polytrue, polyfalse, domain, fast_exit)
	else
		f = plot()
		coveredDomain = "not computed"
		uncoveredDomain = "not computed"
	end

	recap = false

	if juste==true && coveredDomain== true
		recap = true
	elseif juste==true && coveredDomain=="not computed"
		recap = true
	end

	return recap, juste, polytrue, polyfalse, f, coveredDomain, uncoveredDomain
end

# ╔═╡ 13f46772-7495-4702-b719-d325f382b049
"""
	checker(solver, analysedfunction, domain, delta, Polygone, rapide)

	Function allowing to check if a polytopes set respects a tolerance threshold. To make verification easier, colors are used.
In green, these are the correct polytopes, respecting the tolerance threshold.
In red, these are the incorrect polytopes, not respecting the tolerance threshold.
In white, this is the domain not covered by the polytope set entered.

# Arguments
...

-'solver::String' : Choice of the solver.
-'analysedfunction::Function' : Function to be tested.
-'domain::Array{Array{Float64}(undef,2)}(undef,2)' : universe on which the function is tested.
-'delta::Float64' : tolerance threshold
-'Polygones::Array{Array{Array{Float64}(undef,2)}(undef,m)}(undef,n)' : set of polygones where n is the number of polygones and m is the number of vertices for each polygone.
-'fast_exit::Int8' : parameter allowing to exit the checker as soon as a polygone is incorrect. (0: no exit (default), 1:Interruption)

...

# Return
...
-'recap::Boolean' : indicates if all the domain is covered and all polygones are correct.
-'juste::Boolean' : indicates if all polygones are correct.
-'polytrue' : set of correct polygones with some information concerning the resolution of the model.
-'polyfalse' : set of incorrect polygones with some information concerning the resolution of the model.
Information return by the arrays polytrue and polyfalse are :
[[[xvertex1, yvertex1,zvertex1],[xvertex2, yvertex2,zvertex2],[xvertex3, yvertex3,zvertex3]], errormax, xerrormax, yerrormax],
. . . ]
where [[xvertex1, yvertex1,zvertex1],[xvertex2, yvertex2,zvertex2],[xvertex3, yvertex3,zvertex3]] is the triangle analysed,
errormax is the value given by the difference of the linear approximation at the vertices and the evaluation of the nonlinear function at the point (xerrormax,yerrormax)
-'f' : Plot that shows if the triangulation has a maximal error under the tolerance threshold or not thanks to colors: green if it is correct else red.
-'coveredDomain::Boolean' : indicates if all the domain is covered.
-'uncoveredDomain' : Domain not covered.
...
"""

# ╔═╡ ffa54005-76ba-4a0a-be4c-14a808190514

#TT = Any[Any[[2.0, 1.465965965966, 3.4837129998110985], [1.534034034034, 1.0, 1.4693941083222306], [2.0, 1.0, 8.509021784547306]], Any[[2.0, 1.465965965966, 3.129831880267256], [2.0, 2.0, -0.19997842568628954], [1.97742586735, 2.0, -0.3892459346485069], [1.583844458223, 1.049810424189, 2.2354876058897446]], Any[[1.583844458223, 1.049810424189, 1.635332444009415], [1.97742586735, 2.0, 0.28452411284621615], [1.871952099879, 2.0, -0.10093116354584275], [1.20377346196, 1.0, 0.3925699306285373], [1.534034034034, 1.0, 1.5995115535793283]], Any[[1.824572565665, 1.929091516062, 0.3861645687663413], [1.0, 1.488348297315, 1.1265529420678941], [1.0, 1.0, 0.3274880427171787], [1.20377346196, 1.0, -0.03370143256630298]], Any[[1.824572565664, 1.929091516062, -0.13918942464809625], [1.871952099879, 2.0, 0.25301693001854986], [1.239786827467, 2.0, 5.570750269399362], [1.399214046971, 1.701732651049, 0.9034296529047676]], Any[[1.399214046971, 1.701732651049, 0.6349198771797706], [1.239786827467, 2.0, 5.94236009096589], [1.147485308798, 2.0, 6.736713096362809], [1.0, 1.513806094585, 1.5910419889624592], [1.0, 1.488348297315, 1.2551467279639361]], Any[[1.147485308798, 2.0, 6.705944550959444], [1.0, 2.0, 8.520939901223038], [1.0, 1.548058835778, 1.6209585742303645], [1.00737426544, 1.538115789856, 1.3784040218755749]], Any[[1.00737426544, 1.538115789856, 1.3253472812495772], [1.0, 1.548058835778, 1.4506506751020014], [1.0, 1.513806094585, 1.168804131599173]]]
#recap,juste, polytrue, polyfalse, f, coveredDomain, uncoveredDomain = checker("couenne", f8, createDomain([1.0,2.0,1.0,2.0]), 0.5, TT)

#######################################################################recap,juste, polytrue, polyfalse, f, coveredDomain, uncoveredDomain = checker("couenne", f4, createDomain([0.5,2.0,0.5,2.0]), 0.03, test)
#println(polyfalse)
# ╔═╡ d6c735ba-cff7-49ae-8428-c4a16a0b0a1c


# ╔═╡ Cell order:
# ╠═6f222ff3-57b0-4048-bbec-6f1d5742654e
# ╠═86a7842c-855b-46fb-b4e8-afccada2996c
# ╠═aa8dbdcb-2f02-4733-a783-0a9da14b4090
# ╠═0db861fd-2d46-49c6-9ae4-3dfe9c0c89ae
# ╠═9fe67875-c0d5-42a1-a887-2e21f39df090
# ╠═cbc4c2ec-64c4-4012-8283-a237c55eb6b7
# ╠═06d9059a-e9ff-4416-a070-be7c24f83f3d
# ╠═e78c2699-4df0-4bdf-9b74-1ce09bdc3c2f
# ╠═fa6ae86f-fbb7-43fa-af56-e0d7e2240d26
# ╠═b6654944-4fa0-47d7-9631-81cfb882f99d
# ╠═5d56e2eb-fbfb-471b-ada2-ae1f5c208a34
# ╠═4293c8df-e0c8-42e3-b61c-696a5055de12
# ╠═354996bb-d6dd-41a8-bc4c-e023792ebdce
# ╠═618b0169-e701-4aca-945e-798494015bf0
# ╠═3f8ff1ac-e579-11eb-1f8e-279ea9199229
# ╠═509f5ea5-bafa-47a9-9b77-9db708dcb122
# ╠═13f46772-7495-4702-b719-d325f382b049
# ╠═ffa54005-76ba-4a0a-be4c-14a808190514
# ╠═d6c735ba-cff7-49ae-8428-c4a16a0b0a1c
