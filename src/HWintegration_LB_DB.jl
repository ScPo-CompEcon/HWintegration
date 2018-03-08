module HWintegration

	const A_SOL = 4
	println("Analytically, the change of CR equals 4")
	# imports
	using FastGaussQuadrature
	using Roots
	using Sobol
	using Plots
	using Distributions
	using NullableArrays
	using DataFrames
	using DataStructures  # OrderedDict

	# set random seed
	srand(12345)

	# demand function
	function f(x)
		2*x^(-.5)
	end

	# makes a plot for question 1
	function plot_q1()
	global dem = plot(f, 0.1, 10, title="Demand", xaxis=("p"), yaxis=("quantitiy"), labels = "f(x)")
	hline!([2, 1], label = "")
	end

	function question_1a(np)
		for n in (10,15,20)
			np = gausslegendre(n)
			glx = 3/2*np[1] + 3/2
			gly = [f(3/2*x + 5/2) for x in np[1]]
			integ = 3/2*sum([f(3/2*x + 5/2) for x in np[1]] .* np[2])
			global gl = scatter(glx, gly, labels = "GL", title = "$n nodes")
			error = A_SOL-integ
			println("GL-M: integral = $integ and distance to optimal sol. = $error (with $n nodes)")
	end
end


	function question_1b(n)
		for n in (10,15,20)
			mcx = rand(Uniform(1,4), n)
			mcy = f.(mcx)
			integ = 3/n * sum([f(x) for x in mcx])
			global mc = scatter(mcx, mcy, labels = "MC", title = "$n nodes")
			error = A_SOL-integ
			println("MC: integral = $integ and distance to optimal sol. = $error (with $n nodes)")
	end
end

	function question_1c(n)
		for n in (10,15,20)
			s = SobolSeq(1)
			qmcx = 3*[ hcat([next(s) for i = 1:n])[x][1] for x = 1:n]+1
			qmcy = f.(qmcx)
			integ = 3/n * sum([f(x) for x in qmcx])
			global qmc = scatter(qmcx, qmcy, labels = "Quasi MC", title = "$n nodes")
			println("Quasi MC: integral = $integ and distance to optimal sol. = $error (with $n nodes)")
	end
end

	# question 2
# theta1*p^(-1) + theta2*p^(-0.5) = 2

	function question_2a()
		function p(t = [0,0])
			function eq(x)
				exp(t[1])* float(x)^(-1) + exp(t[2])*float(x)^(-0.5) - 2
			end
			fzero(eq, [0,1000])
		end
		immutable MvNormal{Cov<:AbstractPDMat,Mean<:Vector} <: AbstractMvNormal
		            μ::Mean
		            Σ::Cov
		end
		μ = [0,0]
		Σ = [0.02   0.01 ; 0.01  0.01]
		rand(MvNormal(μ, Σ), n)

		nodes = Any[]
		push!(nodes,repeat(rules["hermite"][1],inner=[1],outer=[10]))  # dim1
		push!(nodes,repeat(rules["hermite"][1],inner=[2],outer=[5]))  # dim2
		weights = kron(rules["hermite"][2],kron(rules["hermite"][2],rules["hermite"][2]))
		df = hcat(DataFrame(weights=weights),DataFrame(nodes,[:dim1,:dim2]))

	end

	function question_2b(n)
		for n in (10,15,20)
		μ = [0,0]
		Σ = [0.02   0.01 ; 0.01  0.01]
		draws = rand(MvNormal(μ, Σ),  n)
		draws = [[draws[1,i], draws[2,i]] for i in 1:n]
		function p(t = [0,0])
			function eq(x)
				exp(t[1])* float(x)^(-1) + exp(t[2])*float(x)^(-0.5) - 2
			end
			fzero(eq, [0,1000])
		end
		expec= 1/n * sum(p.(draws))
		
		var = 1/n * sum((p.(draws) - expec).^2)
		println("By Monte Carlo, the expectation is $expec and the variance is $var")
	end
end

	function question_2bonus()

	end

	# function to run all questions
	function runall()
		info("Running all of HWintegration")
		info("question 1:")
		plot_q1()
		for n in (10,15,20)
			info("============================")
			info("now showing results for n=$n")
			info("question 1b:")
			question_1a(n)
			info("question 1c:")
			question_1b(n)
			info("question 1d:")
			question_1c(n)
			println("")
			allplots = plot(dem, gl, mc, qmc, layout = (2,2))
			display(allplots)
			savefig("n=$n.png")
		end
			info("question 2a:")
			q2 = question_2a(n)
			println(q2)
			info("question 2b:")
			for n in (10,15,20)
			info("Now showing results for n = $n")
			question_2b(n)
			println(q2b)
		end
			info("bonus question: Quasi monte carlo:")
			q2bo = question_2bonus(n)
			println(q2bo)
		end
	info("end of HWintegration")
end
