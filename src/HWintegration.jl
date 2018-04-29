
module HWintegration

	const A_SOL = 4  # enter your analytic solution. -18 is wrong.

	# imports
	using FastGaussQuadrature
	using Roots
	using Sobol
	using Plots
	using Distributions
	using NullableArrays
	using DataStructures  # OrderedDict


	# set random seed
	srand(12345)

	# demand function
	function demand(p)
		return 2 * p^-0.5
	end
	# gauss-legendre adjustment factors for map change
	function gausstransform(a,b,f,p)
		#a = 1
		#b = 4
		#p = 3
		x = (b - a)/2
		j = (b + a)/2
		newp = x * p + j
		newq = f(newp)
		return (newp, newq)
	end
	# eqm condition for question 2

	# makes a plot for questions 1
	function plot_q1()
		global plot1 = plot(demand, yaxis = ("Q", (0, 10)), title="Demand of Good", xaxis=("P", (0,5)), labels = "Demand function", color="magenta")
		hline!([2], label = "Quantity at P = 1", color="turquoise")
		hline!([1], label = "Quantity at P = 4", color="blue")
	end


	function question_1b(n)
		gl = gausslegendre(n)
		GL_SOL = 3/2*sum([gausstransform(1,4,demand,p)[2] for p in gl[1]] .* gl[2])
		error = GL_SOL - A_SOL
		glx = [gausstransform(1,4,demand,p)[1] for p in gl[1]]
		gly = [gausstransform(1,4,demand,p)[2] for p in gl[1]]
		global glplot = plot(demand, labels = "True demand function", yaxis = ("Q", (0, 10)))
		scatter!(glx, gly, labels = "Gauss-Legendre", title = "n = $n")
		push!(plotall,glplot)
		println("Using the Gauss-Legendre method, we get solution $GL_SOL with error $error")
	end

	function question_1c(n)
		mcx = rand(Uniform(1,4), n)
		mcy = demand.(mcx)
		MC_SOL = 3/n * sum([demand(p) for p in mcx])
		error = MC_SOL - A_SOL
		global mcplot = plot(demand, labels = "True demand function", yaxis = ("Q", (0, 10)))
		scatter!(mcx, mcy, labels = "Monte-Carlo", title = "n = $n")
		push!(plotall,mcplot)
		println("Using the Monte-Carlo method, we get solution $MC_SOL with error $error")
	end


	function question_1d(n)
		s = SobolSeq(1,1,4)
		sbx = [hcat([next(s) for i = 1:n])[x][1] for x = 1:n]#[x][1] for x = 1:n]+1
		sby = demand.(sbx)
		SB_SOL = 3/n * sum([demand(p) for p in sbx])
		error = SB_SOL - A_SOL
		global sbplot = plot(demand, labels = "True demand function", yaxis = ("Q", (0, 10)))
		scatter!(sbx, sby, labels = "Sobol-based Quasi MonteCarlo", title = "n = $n")
		push!(plotall,sbplot)
		println("Using the Quasi-Monte-Carlo method with Sobol sequences, we get solution $SB_SOL with error $error")
	end

	# question 2

	function price(θ1, θ2)
		d(p) = θ1 * p^-1 + θ2 * p^-0.5 - 2
		fzero(d, [0,1000])
	end

	function question_2a(n)
		μ = [0,0]
		Σ = [0.02 0.01;0.01  0.02]
		logθ = rand(MvNormal(μ, Σ), n)
		θs = exp.(logθ)
		θ1 = θs[1,:]
		θ2 = θs[2,:]
		gh = gausshermite(10)
		ghx = kron(gh[1], gh[1])
		Ep = (1/sqrt(π) * sum(gh[2] .* price.(θ1,θ2)))
		Varp = (1/sqrt(π) * sum(gh[2] .* ((price.(θ1,θ2)- Ep).^2)))
		println("Using the Gauss-Hermite method, we get E[p] = $Ep and Var[p] = $Varp")
	end

	function question_2b(n)
		n = 10
		μ = [0,0]
		Σ = [0.02 0.01;0.01  0.02]
		logθ = rand(MvNormal(μ, Σ), n)
		θs = exp.(logθ)
		θ1 = θs[1,:]
		θ2 = θs[2,:]
		Ep = 1/n * sum(price.(θ1,θ2))
		Varp = 1/n * sum((price.(θ1,θ2) - Ep).^2)
		println("Using the Monte-Carlo method, we get E[p] = $Ep and Var[p] = $Varp")
	end
	#
	# function question_2bonus(n)
	#
	#
	#
	#
	# end

	# function to run all questions
	function runall()
		info("Running all of HWintegration")
		info("question 1:")
		plot_q1()
		default(size=(600,600))
		global plotall = Any[]
		for n in (10,15,20)
			info("============================")
			info("now showing results for n=$n")
			info("question 1b:")
			question_1b(n)	# make sure your function prints some kind of result!
			info("question 1c:")
			question_1c(n)
			info("question 1d:")
			question_1d(n)
		end
		plot(plotall...)
		println("")


		info("question 2a:")
		q2 = question_2a(10)
		println(q2)
		info("question 2b:")
		q2b = question_2b(10)
		println(q2b)

		#info("bonus question: Quasi monte carlo:")
		#q2bo = question_2bonus(n)
		#println(q2bo)
		#end
	end
	info("end of HWintegration")

end
